#-------------------
#----DESCRIPTION----
#-------------------
#Distributions: Normal & T (truncated or not)
#Analysis: GLM & ANOVA (optional - Bayesian Anova)
#Design: Between-Subjects

#This set of functions can be used to simulate the datasets using Drift Diffusion Model.
#The data is being sampled from the multivariate normal distribution, to account for 
#correlations between the parameters
#The sampled parameters are then passed to rdiffusion to generate the data 
#The data is then analysed using GLM and ANOVA
#-----------------Sampling the parameter values----
#Function accepts the vector of means(mu) and covariance matrix(sigma) of the selcted parameter 
#values and returns the vector of parameter values sampled from a best-fit distribution
par_from_val <- function(mu,sigma) {
                 gen <- rtmvnorm(1, mean = mu, 
                  sigma = sigma,
                  lower=c(0,0,0,0,-10,0), 
                  upper=c(Inf, Inf, 10, Inf, 10, 1),
                  H = NULL)
  gen <- as.data.frame(gen)
  gen <- gen %>% 
    rename(
      a = V1,
      st0 = V2,
      sv = V3,
      t0 = V4,
      v = V5,
      z = V6)
}
#-----------------Simulating the data-------
simulate_dt <- function(mu, sigma, n, pp, group) {
  params <- function(mu, sigma) {
    values <- par_from_val(mu, sigma)
    trials <- rdiffusion(
      n = n,
      v = values$v,
      a = values$a,
      t0 = values$t0,
      sv = values$sv,
      st0 = values$st0,
      z = values$z,
      stop_on_error = FALSE
    ) #simulates n trials for 1 participant
  }
  repeat
  {result <- rerun(pp, params(mu, sigma)) %>%
    rbindlist(., idcol = TRUE) %>% #adds id column
    mutate(group = rep(group))  #adds group number
  result <- result %>% rename("id" = ".id") %>% 
    mutate(id = paste0(group, "_", id)) %>%
    group_by(id) %>%
    filter(mean(rt) != 0) %>% 
    ungroup() #removes participants for whom rdiffusion produced an error
  if (as.numeric(as.character(summarise(result, n_distinct(id)))) == pp) break
  } #if r diffuaion produced an error and there are fewer than 30 participants, 
  #reruns the simulation
  return(result)
}

#-----------------Analysing the simulated data----
#Analysing 1 dataset with 'pp' participants and 'n' number of trials 
data_analysis <- function(mu,sigma,n,pp) {
  s1 <- simulate_dt(
    mu = mu,
    sigma = sigma,
    n = n,
    pp = pp,
    group = 1  #adds column with a group number
  ) 
  pp_g1 <- as.numeric(as.character(summarise(s1, n_distinct(id)))) #calculates final number of pp's
  s2 <- simulate_dt(
    mu = mu, 
    sigma = sigma, 
    n = n,
    pp = pp,
    group = 2
  )
  pp_g2 <- as.numeric(as.character(summarise(s2, n_distinct(id))))
  ss <- rbind(s1, s2) %>%  #binds simulated datasets together
    mutate(response = ifelse(response == "upper", 1, 0))  #transforms response into numeric
  mean_prop <- ss %>%
    group_by(group) %>%
    summarise(diff_props = mean(response)) #proportion of upper responses
  g1_prop <- mean_prop$diff_props[1] #proportion of upper for group 1
  g2_prop <- mean_prop$diff_props[2] #proportion of upper for group 2
  mean_prop_real <- ss %>%
    summarise(diff_props = mean(response))
  mean_prop_real <- as.numeric(mean_prop_real) #proportion of upper responses
  #for both groups combined
  diff_props <- g1_prop - g2_prop #difference in proportion
  if (mean_prop_real == 1) {
    aov_p = 1}#in cases when all responses are the same, set p-value to 1
  #so that anova does not produce an error
  #particularly important when number of trials is low
  else {
    aov_ss <- aov_ez(
      id = "id",
      dv = "response",
      data = ss,
      between = "group", 
      fun_aggregate = mean
    ) #runs anova
    aov_p <- summary(aov_ss)[["Pr(>F)"]][[1]] #extracting p-value 
  }
  ss_for_glm <- ss %>% 
    group_by(id,group) %>% 
    summarise(resp_prop = mean(response), 
              n_trials = n()) %>% 
    ungroup %>% 
    mutate(group = factor(group))
  glm_ss <- glm(
    resp_prop ~ group,
    data = ss_for_glm,
    weights = n_trials,
    family = binomial
  ) #runs glm
  glm_anova <- car::Anova(glm_ss, type = 3)
  glm_p <- glm_anova$`Pr(>Chisq)`  #extracting p-value
  # bf_ss <- lmBF(response ~ group,
  #               data = ss) #runs bayesian anova
  #runs Bayesian Anova
  # bf <- extractBF(bf_ss)$bf #extracts bayes factor
  data <-
    tibble(
      aov_p = aov_p,
      glm_p = glm_p,
      #bf = bf,
      diff_props = diff_props,
      mean_prop_real = mean_prop_real,
      g1_prop = g1_prop,
      g2_prop = g2_prop,
      n_g1 = n,
      n_g2 = n,
      pp_g1 = pp_g1,
      pp_g2 = pp_g2,
      a_mean = mu[1],
      a_sd = sqrt(sigma[1,1]/mu[1]),
      st0_mean = mu[2],
      st_sd = sqrt(sigma[2,2]/mu[2]),
      sv_mean = mu[3],
      sv_sd = sqrt(sigma[3,3]/1),
      t0_mean = mu[4],
      t0_sd = sqrt(sigma[4,4]/mu[4]),
      v_mean = mu[5],
      v_sd = sqrt(sigma[5,5]/0.25),
      z_mean = mu[6],
      z_sd = sqrt(sigma[6,6]/mu[6])
    ) #produces list with all the values required for further analysis
  data
}

reruns <- function(mu, sigma,n,pp,runs) {
  result <- rerun(runs,
                  data_analysis(mu, sigma,n,pp)) %>% map_dfr(., as_tibble)
  return(as.data.frame(result))
}
