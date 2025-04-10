"
Functions used during the simulation process. 
File does not need to be loaded separately as it is loaded when functions.R is. 
"


get_nac_tau = function(tauvec){
  # Based on 3 taus, return largest and build average above the other two
  # Assumption: 3rd is largest tau
  tau = c(mean(tauvec[1:2]), tauvec[3])
  return(tau)
}



nac_tau2theta = function(family_index, tau){
  return(
    c(
      HAC::tau2theta(tau = tau[1], type = family_index),
      HAC::tau2theta(tau = tau[2], type = family_index)
    )
  )
}
nac_theta2tau = function(family_index, theta){
  theta = unname(theta)
  return(
    c(
      HAC::theta2tau(theta = theta[1], type = family_index),
      HAC::theta2tau(theta = theta[2], type = family_index)
    )
  )
}

fit_ac = function(mat, cops, debug = FALSE){
  # HAC does not deal well with non-nested ACs, i.e. symmetric ACs. The fit seems to work, but then the logLikelihood cannot be evaluated
  # Solution: copula Package: Can fit ONLY non-nested ACs for more than 2 variables
  # Copula package uses bbmle to fit (see: https://cran.r-project.org/web/packages/bbmle/bbmle.pdf)
  # Due to the implementation in this package: 1) fit possible models 2) Use AIC to select the best one 
  if (debug) browser()
  ac_fits = lapply(cops, function(cop) copula::emle(u = mat, cop = copula::onacopula(family = cop, nacStructure = C(1, 1:3))))
  
  # Select the one with smallest AIC
  # NOTE HERE: Since p identical, we can just select smallest negative loglikelihood
  ac_lls = lapply(ac_fits, function(fit) - fit@min) # min gives the NEGATIVE loglikelihood
  ac_best_fit = which.max(ac_lls)
  # if (is.integer(ac_best_fit) && length(ac_best_fit) == 0) ac_best_fit = sample(1:3, 1) # TODO: TEMP SOL: In case of it breaking down, select one copula at random
  ac_mle = ac_fits[[ac_best_fit]]@coef[[1]]
  
  # Actual copula
  ac = copula::onacopulaL(family = cops[ac_best_fit], nacList = list(ac_fits[[ac_best_fit]]@coef, 1:3))
  attr(ac, "logLik") = ac_lls[ac_best_fit]
  
  return(ac)
}

fit_nac = function(mat, families){
  # Estimate copula using implemented selection method
  # HAC package does not offer way to select the copula family
  # Solution: Fit all considered families and use AIC to select accordingly
  
  # Allowed copula families
  types = unlist(unname(families))
  
  # Select one with smallest AIC. NOTE: p const --> use logLik
  nac_fits = lapply(types, function(type) HAC::estimate.copula(mat, type = type))
  # If any of the NACs suggests a non-nested structure, remove it. 
  # Two reasons: 
  # 1) I already fit a non-nested structure
  # 2) The HAC logLik function breaks down for these functions. And since I already fit one anyway, I use the simplest way to deal with this issue
  # Non-nested NACs have a tree structure of length 4. Use this to filter them out:
  nac_fits = nac_fits[!lapply(nac_fits, function(nac) length(nac$tree)) == 4]
  
  nac_lls = lapply(nac_fits, function(nac) get_nac_loglik(mat, nac))
  nac_best_fit = which.max(nac_lls)
  
  # Best fit
  nac = nac_fits[[nac_best_fit]]
  
  attr(nac, "logLik") = nac_lls[nac_best_fit]
  attr(nac, "theta") = get_nac_estimates(nac)
  attr(nac, "tau") = nac_theta2tau(nac$type, attr(nac, "theta"))
  
  return(nac)
}

get_nac_loglik = function(mat, nac_mdl){
  return(HAC::to.logLik(X = mat, hac = nac_mdl, eval = TRUE))
}

get_nac_estimates = function(nac_mdl){
  return(c(nac_mdl$tree[[3]], nac_mdl$tree[[1]][[3]]))
}

fit_vine = function(mat, families){
  # The VineCopula package is faster and applicable as long as we do not also estimate covariates
  # Families: 3, 4, 5 (see docu)
  vine = VineCopula::RVineStructureSelect(data = mat, rotations = FALSE, familyset = families)
  return(vine)
}

get_nac_name = function(idx, inverse_copula_families = list("1" = "Gumbel", "3" = "Clayton", "5" = "Frank")){
  return(inverse_copula_families[as.character(idx)][[1]])
}

get_ac_estimate = function(ac) unname(ac@copula@tau( ac@copula@theta ))

run_one_nac = function(
    seed,
    n,
    cop,
    taus,
    copula_families = list("Gumbel" = 1, "Clayton" = 3, "Frank" = 5) # Copula families contained in the package
    ){
  # Simulate from true model ------------------------------------------------
  # For reproducibility
  set.seed(seed)
  
  # Draw random tau vector and determine nac structure from this
  rtau = sample(taus, 1)
  river = names(rtau)
  tau = get_nac_tau(tauvec = rtau[[1]])
  theta = nac_tau2theta(family_index = copula_families[cop], tau = tau)
  
  # Create HAC-object
  mdl = HAC::hac.full(type = copula_families[cop], y = c("v3", "v2", "v1"), theta = theta)
  # Draw sample according to true_mdl (simulated sample)
  mat = HAC::rHAC(n, mdl) 
  
  # Fit models --------------------------------------------------------------
  ac = fit_ac(mat, names(copula_families))
  
  nac = fit_nac(mat, copula_families)
 
  vine = fit_vine(mat, families = c(3, 4, 5))
  
  # Save results ------------------------------------------------------------
  res = data.frame(
    list(
      # True values
      seed = seed,
      river = river,
      true_tau_outer = tau[1],
      true_theta_outer = theta[1],
      true_tau_inner = tau[2],
      true_theta_inner = theta[2],
      # Estimated values
      fit_cop = get_nac_name(nac),
      fit_theta_outer = attr(nac, "theta")[1],
      fit_theta_inner = attr(nac, "theta")[2],
      fit_tau_outer = attr(nac, "tau")[1],
      fit_tau_inner = attr(nac, "tau")[2],
      nac_aic = ll2aic(ll = attr(nac, "logLik")[[1]], p = 2),
      nac_kl = klMonteCarlo(mdl, nac),
      # Other copulas / Misspecifications
      # Results of AC fit
      ac_tau = get_ac_estimate(ac),
      ac_aic = ll2aic(ll = attr(ac, "logLik")[[1]], p = 1),
      ac_kl = klMonteCarlo(mdl, HAC::nacopula2hac(ac), est_mdl_AC = TRUE),
      # Results for Vine
      vine_aic = vine$AIC,
      vine_kl = klMonteCarlo(mdl, vine, est_mdl_vine = TRUE)
    )
  )
  return(res)
}




# For each sample size (n), each copula family, each dependence structure (AC, NAC, Vine), tau:
#   Repeat B times: 
#     Draw random sample  
#   For each sample, fit AC, NAC and Vine copula and evaluate
# run_one_n_ac = function(
#     seed, 
#     n,
#     cop,
#     dep,
#     copula_families = list("Gumbel" = 1, "Clayton" = 3, "Frank" = 5), # Copula families contained in the package
#     beta_a = 2.5, 
#     beta_b = 1.5
#   ){
#   # Simulate from true model ------------------------------------------------
#   # For reproducibility
#   set.seed(seed)
#   
#   # If symmetric, the values of the 'inner' and 'outer' tau bzw. copula parameter are simply identical
#   tau_inner = rbeta(n = 1, shape1 = beta_a, shape2 = beta_b)
#   tau_outer = tau_inner # base case is symmetric AC
#   if (dep == "nac") tau_outer = 1/2 * tau_inner # Simple solution to ensure outer parameter is smaller than inner
#   
#   # Define a selection of possible tau
#   # Derive a selection from the empirical observations? 
#   # AC: Take average of the observed tau
#   # NAC: Take most nested and average of other 2
#   # Vine: Take the observed taus
#   
#   # Package uses integer for copulas. Get this integer 
#   fam = copula_families[cop]
#   
#   # Create HAC-object
#   true_mdl = HAC::hac.full(
#     type = fam, 
#     y = c("v3", "v2", "v1"), 
#     theta = c(HAC::tau2theta(tau = tau_outer, type = fam), HAC::tau2theta(tau = tau_inner, type = fam))
#   )
#   
#   # Draw sample according to true_mdl (simulated sample)
#   mat = HAC::rHAC(n, true_mdl) 
#   attr(mat, "true_mdl") = true_mdl 
#   attr(mat, "tau") = c(outer = tau_outer, inner = tau_inner)
#   
#   # Fit models --------------------------------------------------------------
#   
#   # I) (symmetric) Archimedean copulas ----
#   # Funny enough, HAC does not deal well with non-nested ACs, i.e. symmetric ACs. The fit seems to work, but then the logLikelihood cannot be evaluated
#   # Lucky me, the copula package can fit ONLY non-nested ACs for more than 2 variables
#   # Copula packages are a mess, wtf...
#   # Also, due to the implementation in this package, I first fit all models, then use AIC to select the best one and finally create the best model
#   # Copula package uses bbmle to fit (see: https://cran.r-project.org/web/packages/bbmle/bbmle.pdf)
#   ac_fits = lapply(names(copula_families), function(name) copula::emle(u = mat, cop = copula::onacopula(family = name, nacStructure = C(1, 1:3))))
#   # Select the one with smallest AIC
#   # NOTE HERE: Since p identical, we can just select smallest negative loglikelihood
#   ac_lls = lapply(ac_fits, function(fit) - fit@min) # min gives the NEGATIVE loglikelihood
#   ac_best_fit = which.max(ac_lls)
#   ac_mle = ac_fits[[ac_best_fit]]@coef[[1]]
#   # Actual copula
#   ac_mdl = copula::onacopulaL(family = names(copula_families)[ac_best_fit], nacList = list(ac_fits[[ac_best_fit]]@coef, 1:3))
# 
#   
#   # II) Nested Archimedean copulas ----
#   # Estimate copula using implemented selection method
#   nac_mdl = HAC::estimate.copula(mat) # TODO There is no structure selection!! Again, solve this using AIC use type 1, 3, 5
#   # Calculate ll to evaluate AIC
#   nac_ll  = HAC::to.logLik(X = mat, hac = nac_mdl, eval = TRUE)
#   
#   estimates = c(outer = nac_mdl$tree[[3]], inner = nac_mdl$tree[[1]][[3]])
#   
#   # III) Vine copulas ----
#   # The VineCopula package is faster and applicable as long as we do not also estimate covariates
#   vine_mdl = VineCopula::RVineStructureSelect(
#     data = mat, 
#     rotations = FALSE, familyset = c(3, 4, 5) # Only allow for the considered copula families
#   )
#   vine_mdl
#   
#   # Save results in df ----
#   res = data.frame(
#     list(
#       seed = seed,
#       true_tau_inner = tau_inner,
#       true_tau_outer = tau_outer,
#       # Results of AC fit
#       ac_selected_cop = names(copula_families)[ac_best_fit],
#       ac_theta = ac_mle,
#       ac_tau = HAC::theta2tau(ac_mle, type = copula_families[names(copula_families)[ac_best_fit]]),
#       ac_aic = ll2aic(ll = ac_lls[[ac_best_fit]], p = 1),
#       ac_bic = ll2bic(ll = ac_lls[[ac_best_fit]], p = 1),
#       ac_kl = klMonteCarlo(true_mdl, HAC::nacopula2hac(ac_mdl), est_mdl_AC = TRUE),
#       # Results of NAC fit
#       nac_selected_cop = HAC::hac2nacopula(nac_mdl)@copula@name,
#       nac_theta_outer = estimates["outer"][[1]],
#       nac_theta_inner = estimates["inner"][[1]],
#       nac_tau_outer = HAC::theta2tau(theta = estimates["outer"][[1]], type = nac_mdl$type),
#       nac_tau_inner = HAC::theta2tau(theta = estimates["inner"][[1]], type = nac_mdl$type),
#       nac_aic = ll2aic(ll = nac_ll, p = 2),
#       nac_bic = ll2bic(ll = nac_ll, p = 2),
#       nac_kl = klMonteCarlo(true_mdl, nac_mdl),
#       # Results of Vine fit
#       vine_aic = vine_mdl$AIC,
#       vine_bic = vine_mdl$BIC,
#       vine_kl = klMonteCarlo(true_mdl = true_mdl, est_mdl = vine_mdl, est_mdl_vine = TRUE)
#     )
#   )
# 
#   return(res)
# }
run_one_vine = function(
    seed,
    n,
    taus,
    vine_copula_families = list("Gumbel" = 4, "Clayton" = 3, "Frank" = 5), # Copula families contained in the VineCopula package
    hac_copula_families = list("Gumbel" = 2, "Clayton" = 3, "Frank" = 5), # Copula families contained in HAC copula package
    vine_colnames = c("1-3", "2-3", "1-2")
    ){
  # Simulate from true model ------------------------------------------------
  # For reproducibility
  set.seed(seed)

  # Draw random tau vector and determine nac structure from this
  rtau = sample(taus, 1)
  river = names(rtau)
  tau = unlist(rtau)
  names(tau) = vine_colnames
  
  # Vine matrix defining the (un)conditional copulas in the model
  # 1 - 2 - 3 
  vine_matrix = matrix(
    c(
      1, 0, 0,
      3, 2, 0,
      2, 3, 3
    ),
    nrow = 3, ncol = 3,
    byrow = TRUE
  )

  # Randomly drawing the combination of copula families
  fams = unlist(sample(vine_copula_families, size = 3, replace = T))
  names(fams) = vine_colnames
  # Save families in family copula matrix
  family_matrix = matrix(0, nrow = 3, ncol = 3)
  family_matrix[2, 1] = fams["1-3"]
  family_matrix[3, 1:2] = fams[c("1-2", "2-3")]

  params = matrix(0, nrow = 3, ncol = 3)
  theta = c(
    "1-3" = VineCopula::BiCopTau2Par(fams["1-3"], tau["1-3"]),
    "1-2" = VineCopula::BiCopTau2Par(fams["1-2"], tau["1-2"]),
    "2-3" = VineCopula::BiCopTau2Par(fams["2-3"], tau["2-3"])
  )
  params[2, 1] = theta["1-3"]
  params[3, 1:2] = c(theta["1-2"], theta["2-3"])

  rvmat = VineCopula::RVineMatrix(Matrix = vine_matrix, family = family_matrix, par = params)
  mat = VineCopula::RVineSim(n, RVM = rvmat)

  attr(mat, "rvm") = rvmat
  colnames(mat) = c("v1", "v2", "v3")
  
  # Fit models --------------------------------------------------------------
  ac = fit_ac(mat, names(hac_copula_families))

  nac = fit_nac(mat, hac_copula_families)

  vine = fit_vine(mat, families = c(3, 4, 5))

  # Save results ------------------------------------------------------------
  res = data.frame(
    list(
      # True values
      seed = seed,
      river = river,
      true_tau_12 = tau["1-2"],
      true_fam_12 = fams["1-2"],
      true_theta_12 = theta["1-2"],
      true_tau_23 = tau["2-3"],
      true_fam_23 = fams["2-3"],
      true_theta_23 = theta["2-3"],
      true_tau_13 = tau["1-3"],
      true_fam_13 = fams["1-3"],
      true_theta_13 = theta["1-3"],
      # Estimated values
      fit_tau_12 = vine$tau[3, 1],
      fit_fam_12 = vine$family[3, 1],
      fit_theta_12 = vine$par[3, 1],
      fit_tau_23 = vine$tau[3, 2],
      fit_fam_23 = vine$family[3, 2],
      fit_theta_23 = vine$par[3, 2],
      fit_tau_13 = vine$tau[2, 1],
      fit_fam_13 = vine$family[2, 1],
      fit_theta_13 = vine$par[2, 1],
      vine_aic = vine$AIC,
      vine_kl = vine_klMonteCarlo(rvmat, vine, est = "vine"),
      # Other copulas / Misspecifications
      # NAC fit
      nac_tau_outer = attr(nac, "tau")[1],
      nac_tau_inner = attr(nac, "tau")[2], 
      nac_aic = ll2aic(ll = attr(nac, "logLik")[[1]], p = 2),
      nac_kl = vine_klMonteCarlo(rvmat, nac, est = "nac"),
      # AC fit
      ac_aic = ll2aic(ll = attr(ac, "logLik")[[1]], p = 1),
      ac_kl = vine_klMonteCarlo(rvmat, HAC::nacopula2hac(ac), est = "ac")
    )
  )
  return(res)
}

vine_klMonteCarlo = function(
    rvmat,
    est_mdl, 
    est = "vine", 
    values =  seq(from = 0.01, to = 0.99, length.out = 25) # Keep from (to) relatively high (low) simmplifying numerical stability. Downside: KL not considering full support..........lol.
  ){
  grid = expand.grid(v1 = values, v2 = values, v3 = values)
  
  true_d = VineCopula::RVinePDF(grid, rvmat)
  
  # ifelse cannot return matrices, aparently. Thus, keep it in two if-statements...
  if (est == "vine") est_d = VineCopula::RVinePDF(grid, est_mdl)
  if (est == "ac") for (i in 1:3) est_mdl$tree[[i]] = paste("v", i, sep = "")
  if (est != "vine") est_d = HAC::dHAC(as.matrix(grid), est_mdl)

  return(mean(log(true_d / est_d) ))
}
 
# TODO
# run_one_vine = function(
#     seed, 
#     n,
#     cop,
#     dep
#   ){
#   # Simulate from true model ------------------------------------------------
#   # Tau for the dependence structure
#   tau_13 = rbeta(n = 1, shape1 = beta_a, shape2 = beta_b)
#   tau_23 = 3/5 * tau_13
#   tau_12 = 1/3 * tau_13
#   
#   # Vine matrix defining the (un)conditional copulas in the model
#   # 1 - 2 - 3
#   vine_matrix = matrix(
#     c(
#       1, 0, 0, 
#       3, 2, 0,
#       2, 3, 3
#     ),
#     nrow = 3, ncol = 3,
#     byrow = TRUE
#   )
#   
#   # Randomly drawing the combination of copula families
#   fams = unlist(sample(copula_families, size = 3, replace = T))
#   # Save families in family copula matrix
#   family_matrix = matrix(0, nrow = 3, ncol = 3)
#   family_matrix[2, 1] = fams[1]
#   family_matrix[3, 1:2] = fams[2:3]
#   
#   params = matrix(0, nrow = 3, ncol = 3)
#   params[2, 1] = VineCopula::BiCopTau2Par(fams[1], tau_13)
#   params[3, 1:2] = c(
#     VineCopula::BiCopTau2Par(fams[2], tau_12),
#     VineCopula::BiCopTau2Par(fams[3], tau_23)
#   )
#   
#   rvmat = VineCopula::RVineMatrix(Matrix = vine_matrix, family = family_matrix, par = params)
#   mat = VineCopula::RVineSim(n, RVM = rvmat)
#   
#   attr(mat, "rvm") = rvmat
#   colnames(mat) = c("v1", "v2", "v3")
#   
#   # Fit models --------------------------------------------------------------
#   
#   
#   # I) (symmetric) Archimedean copulas ----
#   # Funny enough, HAC does not deal well with non-nested ACs, i.e. symmetric ACs. The fit seems to work, but then the logLikelihood cannot be evaluated
#   # Lucky me, the copula package can fit ONLY non-nested ACs for more than 2 variables
#   # Copula packages are a mess, wtf...
#   # Also, due to the implementation in this package, I first fit all models, then use AIC to select the best one and finally create the best model
#   # Copula package uses bbmle to fit (see: https://cran.r-project.org/web/packages/bbmle/bbmle.pdf)
#   ac_fits = lapply(names(copula_families), function(name) copula::emle(u = mat, cop = copula::onacopula(family = name, nacStructure = C(1, 1:3))))
#   # Select the one with smallest AIC
#   # NOTE HERE: Since p identical, we can just select smallest negative loglikelihood
#   ac_lls = lapply(ac_fits, function(fit) - fit@min) # min gives the NEGATIVE loglikelihood
#   ac_best_fit = which.max(ac_lls)
#   ac_mle = ac_fits[[ac_best_fit]]@coef[[1]]
#   # Actual copula
#   ac_mdl = copula::onacopulaL(family = names(copula_families)[ac_best_fit], nacList = list(ac_fits[[ac_best_fit]]@coef, 1:3))
#   
#   # II) Nested Archimedean copulas ----
#   # Estimate copula using implemented selection method
#   nac_mdl = HAC::estimate.copula(mat)
#   # Calculate ll to evaluate AIC
#   nac_ll = HAC::to.logLik(X = mat, hac = nac_mdl, eval = TRUE)
#   
#   # III) Vine copulas ----
#   # The VineCopula package is faster and applicable as long as we do not also estimate covariates
#   vine_mdl = VineCopula::RVineStructureSelect(
#     data = mat, 
#     rotations = FALSE, familyset = c(3, 4, 5) # Only allow for the considered copula families
#   )
#   vine_mdl
#   
#   # Save result df ----

# Simulation Analysis -----------------------------------------------------
klplots = function(
    df,
    cop_name,
    dens_alpha = 1,
    col_ac = "red",
    col_nac = "blue",
    col_vine = "green",
    scales = "fixed"
  ){
  p = df |> 
    dplyr::filter(cop == cop_name) |> 
    dplyr::select(n, contains("_kl")) |> 
    tidyr::pivot_longer(
      cols = contains("_kl"),
      names_to = "dep",
      values_to = "kld"
    ) |>
    dplyr::mutate(dep = as.factor(stringr::str_remove(dep, "_kl"))) |>
    ggplot() + 
    geom_density(aes(x = kld, color = dep), alpha = dens_alpha) +
    geom_vline(xintercept = 0, color = "black") + 
    facet_wrap(~ n, scale = scales) +
    labs(
      title = "Kullback Leibler by Sample Size and Fit",
      x = "Kullback Leibler Divergence",
      y = "Density"
    ) +
    theme(legend.position = "bottom") 
    
  return(p)
}

get_vine_famname = function(idx){
  vine_famlist = list("4" = "Gumbel", "3" = "Clayton", "5" = "Frank")
  return(unname(unlist(lapply(idx, function(i) vine_famlist[as.character(i)]))))
}


