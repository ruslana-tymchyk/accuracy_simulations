#-------------------
#----DESCRIPTION----
#-------------------
#Distributions: Normal & T (truncated or not)
#Analysis: GLM & ANOVA (optional - Bayesian Anova)
#Design: Within-Subjects (Repeated Measures)

#This set of function can be used to simulate and analyse the datasets using Drift Diffusion Model.
#The data is being sampled from a distribution that best fits the given parameter 
#The sampled parameters are then passed to rdiffusion to generate the data 
#The data is then analysed using GLM and ANOVA

#-----------------Defining Distributions-------
#This code is required for the t-distribution, as the t-distribution specified in R does not permit to set location and scale parameters
dmyt <- function(x, location, scale, df) {
  1/scale * dt((x - location)/scale, df)
}
pmyt <- function(q, location, scale, df) {
  pt((q-location)/scale, df)
}

qmyt <- function(p, location, scale, df) {
  qt(p,df)*scale + location
}
rmyt <- function(n, location, scale, df) {
  (rt(n, df)*scale) + location
}

#Function that defines a truncated normal distribution
gen_norm <- function(par, a, b) {
  rtrunc(
    1,
    spec = "norm",
    a = a, #lower limit
    b = b, #upper limit
    mean = par[1],
    sd = par[2]
  )
}

#Function that defines a truncated T-distribution
gen_myt <- function(par, a, b) {
  rtrunc(
    1,
    spec = "myt",
    a = a, #lower limit
    b = b, #upper limit
    df = par[1],
    location = par[2],
    scale = par[3]
  )
}

#-----------------Sampling the parameter values-------
#Function accepts the summary statistics of each paramter as an argument and 
#returns the vector of parameter values sampled from a best-fit distribution
par_from_val <-  function(a, z, st0, v, t0, sv) {
  c(
    #defining limits of truncation and distribution for each parameter
    a = gen_norm(par = a, a = 0, b = Inf),
    z = gen_norm(par = z, a = 0, b = 1),
    st0 = unlist(gen_norm(
      par = st0, a = 0, b = Inf
    )),
    v = unlist(gen_myt(
      par = v,
      a = -10 + v[2],
      b = 10 - v[2]
      #addition and substraction of location parameter required 
      #to ensure than the upper and lower thresholds deviate by 10 
      #points from the central point of distribution rather than 
      #from zero 
    )),
    t0 = unlist(gen_myt(
      par = t0, a = 0, b = Inf
    )),
    sv = unlist(gen_myt(
      par = sv, a = 0, b = 10
    ))
  )
}

#-----------------Simulating the data-------
simulate_dt <- function(a, z, st0, v, t0, sv, n, pp) {
  params <- function(a, z, st0, v, t0, sv) {
    repeat{
      values <- par_from_val(a, z, st0, v, t0, sv)
      trials1 <- rdiffusion(
        n = n/2,
        v = values[["v"]],
        a = values[["a"]],
        t0 = values[["t0"]],
        sv = values[["sv"]],
        st0 = values[["st0"]],
        z = values[["z"]],
        stop_on_error = FALSE
      ) 
      trials2 <- rdiffusion(
        n = n/2,
        v = values[["v"]],
        a = values[["a"]],
        t0 = values[["t0"]],
        sv = values[["sv"]],
        st0 = values[["st0"]],
        z = values[["z"]],
        stop_on_error = FALSE
      ) 
      trials1 <- trials1 %>% mutate(frame = 1)
      trials2 <- trials2 %>% mutate(frame = 2)
      trials <- rbind(trials1, trials2)
      if (mean(trials$rt) != 0) {break}
      #when combination of sampled values produces error in rdiffusion, the 
      #values are resampled and the simulation is run again
    }
    return(trials)
  }
  result <- rerun(pp, params(a, z, st0, v, t0, sv)) %>%
    rbindlist(., idcol = TRUE)  %>% #adds id column 
    rename("id" = ".id") 
  return(result)
}

#-----------------Analysing the simulated data-------
#Analysing 1 dataset with 'pp' participants and 'n' number of trials 
data_analysis <- function(a,z,st0,v,t0,sv,
                          n,pp) {
  repeat{
    ss <- simulate_dt(
      a = a,
      z = z,
      st0 = st0,
      v = v,
      t0 = t0,
      sv = sv,
      n = n,
      pp = pp
    ) 
    ss <- ss %>% mutate(response = ifelse(response == "upper", 1, 0))  
    #transforms response into numeric
    mean_prop_frame <- ss %>%
      group_by(frame) %>%
      summarise(diff_props = mean(response)) #proportion of upper responses
    g1_prop <- mean_prop_frame$diff_props[1] #proportion of upper for group 1
    g2_prop <- mean_prop_frame$diff_props[2] #proportion of upper for group 2
    mean_prop_real <- ss %>%
      summarise(diff_props = mean(response))
    mean_prop_real <- as.numeric(mean_prop_real)
    diff_props <- g1_prop - g2_prop #difference in proportion
    if (mean_prop_real == 1 || mean_prop_real == 0) {
      aov_p = 1}
    #in cases when all responses are the same, set p-value to 1
    #so that anova does not produce an error
    #particularly important when number of trials is low
    else {
      aov_ss <- aov_ez(
        id = "id",
        dv = "response",
        data = ss,
        within = "frame", #the repeated measures factor is the within subjects factor
        fun_aggregate = mean
      ) #runs anova
      aov_p <- aov_ss$anova_table[["Pr(>F)"]][[1]] #extracting p-value 
    }
    ss_for_glm <- ss %>% 
      group_by(id,frame) %>% 
      summarise(resp_prop = mean(response), 
                n_trials = n()) %>% 
      ungroup %>% 
      mutate(frame = factor(frame))
    if (mean_prop_real == 1 || mean_prop_real == 0) {
      glm_p = 1
      sing = 0}
    else if(g1_prop == 1 & g2_prop == 0) {
      glm_p = 0
      sing = 1
    }
    else if(g2_prop == 1 & g1_prop == 0) {
      glm_p = 0
      sing = 1
    }
    #the two 'else if' statements ensure that there is no error due to complete separation in glm
    else {
    glm_ss <- glmer(
      resp_prop ~ frame + (frame|id),
      data = ss_for_glm,
      weights = n_trials, 
      family = binomial
    ) #runs glm
    sing <- ifelse(isSingular(glm_ss) == TRUE, 1, 0)
    glm_anova <- car::Anova(glm_ss, type = 3)
    glm_p <- glm_anova$`Pr(>Chisq)`[2] 
    } #extracting p-value
    # bf_ss <- lmBF(response ~ group,
    #               data = ss) #runs bayesian anova
    # bf <- extractBF(bf_ss)$bf #extracts bayes factor
    data <-
      tibble(
        aov_p = aov_p,
        glm_p = glm_p,
        sing = sing,
        #bf = bf,
        diff_props = diff_props,
        mean_prop_real = mean_prop_real,
        g1_prop = g1_prop,
        g2_prop = g2_prop,
        n_g1 = n,
        n_g2 = n,
        a_mean = a[1],
        a_sd = a[2],
        z_mean = z[1],
        z_sd = z[2],
        st0_mean = st0[1],
        st0_sd = st0[2],
        v_df = v[1],
        v_loc = v[2],
        v_scale = v[3],
        t0_df = t0[1],
        t0_loc = t0[2],
        t0_scale = t0[3],
        sv_df = sv[1],
        sv_loc = sv[2],
        sv_scale = sv[3]
      ) #produces list with all the values required for further analysis
    if (data$sing == 0) {break}
    #if sing != 0 glm reached singular fit and that analysis is discarded an re-run
  } 
  data
}

#Repeating analysis for 'runs' number of datasets
reruns <- function(a,z,st0,v,t0,sv,
                   n,pp,runs) {
  result <- rerun(runs,
                  data_analysis(a,z,st0,v,t0,sv,
                                n,pp)
  ) %>% map_dfr(., as_tibble)
  return(as.data.frame(result))
}

