require(pomp)

Csnippet('
          double N, beta2;
          N = S + E01 + E02 + I01 + I02 + P1 + P2 + S12 + S21 + E12 + E21 + I12 + I21 + R;
          beta2 = beta1t*omega; 
          
          //---------------------------------
          // Calculate disease 1 pathway
          // initial infection pathway
          double new_E01 = rbinom(S, 1 - exp(-(beta1t/N*I01 + beta1t/N*I21)*dt));
          double leaving_E01 = rbinom(E01, 1 - exp(-eta*dt));
          double new_I01 = rbinom(leaving_E01, eta / (eta + mu));
          double new_I02_mutant = leaving_E01 - new_I01;
          
          // recovery and short-term protection
          double leaving_I01 = rbinom(I01, 1 - exp(-gamma*dt));
          double leaving_P1 = rbinom(P1, 1 - exp(-alpha*dt));
          
          // secondary infection pathway
          double new_E12 = rbinom(S12, 1 - exp(-(rho*beta2/N*I02 + rho*beta2/N*I12)*dt));
          double leaving_E12 = rbinom(E12, 1 - exp(-eta*dt));
          double leaving_I12 = rbinom(I12, 1 - exp(-gamma*dt));
          
          //---------------------------------
          // Calculate disease 2 pathway
          // initial infection pathway
          double new_E02 = rbinom(S, 1 - exp(-(beta2/N*I02 + beta2/N*I12)*dt));
          double leaving_E02 = rbinom(E02, 1 - exp(-eta*dt));
          
          // recovery and short-term protection
          double leaving_I02 = rbinom(I02, 1 - exp(-gamma*dt));
          double leaving_P2 = rbinom(P2, 1 - exp(-alpha*dt));
          
          // secondary infection pathway
          double new_E21 = rbinom(S21, 1 - exp(-(rho*beta1t/N*I01 + rho*beta1t/N*I21)*dt));
          double leaving_E21 = rbinom(E21, 1 - exp(-eta*dt));
          double leaving_I21 = rbinom(I21, 1 - exp(-gamma*dt));
          
          
          // Accounting all compartments
          S = S - new_E01 - new_E02;
          R = R + leaving_I12 + leaving_I21;
          
          // Disease 1 pathway
          E01 = E01 + new_E01 - leaving_E01;
          I01 = I01 + new_I01 - leaving_I01;
          P1 = P1 + leaving_I01 - leaving_P1;
          S12 = S12 + leaving_P1 - new_E12;
          I12 = I12 + leaving_E12 - leaving_I12;
          NI01 = new_I01;
          NI12 = leaving_E12;
          
          // Disease 2 pathway
          E02 = E02 + new_E02 - leaving_E02;
          I02 = I02 + new_I02_mutant + leaving_E02 - leaving_I02;
          P2 = P2 + leaving_I02 - leaving_P2;
          S21 = S21 + leaving_P2 - new_E21;
          I21 = I21 + leaving_E21 - leaving_I21;
          N_mutant = new_I02_mutant;
          NI02 = new_I02_mutant + leaving_E02;
          NI21 = leaving_E21;
        ') -> two_strain_rprocess

Csnippet("
  S = S_0;
  R = R_0;
  E01 = E01_0;
  I01 = I01_0;
  P1 = P1_0;
  S12 = S12_0;
  E12 = E12_0;
  I12 = I12_0;
  NI01 = NI01_0;
  NI12 = NI12_0;
  E02 = E02_0;
  I02 = I02_0;
  P2 = P2_0;
  S21 = S21_0;
  E21 = E21_0;
  I21 = I21_0;
  NI02 = NI02_0;
  NI21 = NI21_0;
  N_mutant = N_mutant_0;
  ") -> rinit_twostrain
c("eta", #rate to become infectious
  # "beta1", # transmission rate
  "mu", # rate which mutants form
  "gamma", # rate to recover from disease 
  "omega", # relative transmission of disease 2 compared to disease 1
  "alpha", # rate which individuals lose generalized protection
  "rho", # susceptibility to other disease following infection
  "S_0",
  "R_0",
  "E01_0",
  "I01_0",
  "P1_0",
  "S12_0",
  "E12_0",
  "I12_0",
  "NI01_0",
  "NI12_0",
  "E02_0",
  "I02_0",
  "P2_0",
  "S21_0",
  "E21_0",
  "I21_0",
  "NI02_0",
  "NI21_0",
  "N_mutant_0") -> two_strain_paramnames

c("S",
  "R",
  "E01",
  "I01",
  "P1",
  "S12",
  "E12",
  "I12",
  "NI01",
  "NI12",
  "E02",
  "I02",
  "P2",
  "S21",
  "E21",
  "I21",
  "NI02",
  "NI21",
  "N_mutant"
) -> two_strain_statenames

get_init_two_strain_parms <- function(parm_vals = NULL){
    c(eta = 1/2.9, #rate to become infectious
      # beta1 = 0.43, # transmission rate
      mu = 1/1000000, # rate which mutants form
      gamma = 1/7, # rate to recover from disease 
      omega = 1.5, # relative transmission of disease 2 compared to disease 1
      alpha = 1/42, # rate which individuals lose generalized protection
      rho = 0.9, # susceptibility to other disease following infection
      S_0 = 999950,
      R_0 = 0,
      E01_0 = 50,
      I01_0 = 0,
      P1_0 = 0,
      S12_0 = 0,
      E12_0 = 0,
      I12_0 = 0,
      NI01_0 = 0,
      NI12_0 = 0,
      E02_0 = 0,
      I02_0 = 0,
      P2_0 = 0,
      S21_0 = 0,
      E21_0 = 0,
      I21_0 = 0,
      NI02_0 = 0,
      NI21_0 = 0,
      N_mutant_0 = 0
    ) -> baseline_parms
  if(!is.null(parm_vals)){
    replace_parms(parm_vals, baseline_parms)
  } else {
    baseline_parms
  }
}

get_twostrain_accum_vars <- function(){
  c('NI01', 'NI12', 'NI02', 'NI21', 'N_mutant')
}

get_data_and_covariate_table <- function(rt_shape = 'constant', 
                                         gamma = 1/7,
                                         num_days = 400){
  ## rt_shape should be 'constant-high', 'constant-med', 'constant-low', 'cyclical-med', 'cyclical-high'
  day_vec <- 0:num_days
  
  if(rt_shape == 'constant-high'){
    betat <- rep(3*gamma, length(day_vec))
  } else if(rt_shape == 'constant-med'){
    betat <- rep(1.5*gamma, length(day_vec))
  } else if(rt_shape == 'constant-low'){
    betat <- rep(0.95*gamma, length(day_vec))
  } else if(rt_shape == 'cyclical-med'){
    betat <- 1*(1+.2*cos(2*pi*(day_vec)/90))*gamma
  } else if(rt_shape == 'cyclical-high'){
    betat <- 2*(1+.5*cos(2*pi*(day_vec)/90))*gamma
  } else{
    stop("Incorrect input for rt_shape")
  }
  
  list(data = tibble(day = day_vec),
       covars = covariate_table(day = day_vec,
                                beta1t = betat, 
                                times = 'day'))
}
