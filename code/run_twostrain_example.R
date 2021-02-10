library(tidyverse)
library(pomp)

sapply(c("code/covid_variant_pomp.R", 
         'code/helper_fxns.R'), source)


init_parms <- get_init_two_strain_parms()
accum_vars <- get_twostrain_accum_vars()
data_list <- get_data_and_covariate_table(rt_shape = 'cyclical-med')
 
pomp(data = data_list$data,
     times="day",
     t0 = 0,
     covar = data_list$covars,
     rinit=rinit_twostrain,
     rprocess = pomp::euler(two_strain_rprocess, 1),
     statenames= two_strain_statenames,
     paramnames= two_strain_paramnames,
     params = init_parms,
     accumvars = accum_vars
     ) -> two_strain_model



sims <- simulate(two_strain_model,
         nsim = 1000) %>% 
  as('data.frame') %>% 
  as_tibble() 

sims %>%   
  mutate(newI1 = NI01 + NI21,
         newI2 = NI02 + N_mutant + NI12) %>% 
  select(.id, day, newI1, newI2) %>% 
  gather(key, value, newI1:newI2) %>% 
  ggplot(aes(day, value, group = .id)) + 
    geom_line(alpha = .3) +
    facet_wrap(~key)

