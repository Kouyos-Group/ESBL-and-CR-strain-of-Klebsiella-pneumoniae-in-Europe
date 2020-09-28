setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#Loading libraries
library(deSolve)
library(parallel)
library(Rcpp)
library(numDeriv)
#Loading the model written in C
sourceCpp("SCCH.cpp") 
#Main function, recieve initial parameters and give back number of infected individuals each year(at the end of the year)
sourceCpp("SCCII2.cpp")
#Loading the csv files with data
input_names<- c("Croatia","Denmark","Finland","France", "Greece", "Hungary", "Italy", "Netherlands", "Norway","Portugal", "Sweden")
for (c in 1:length(input_names)){
  name_curr <- paste(input_names[c], "3rd_4th_gen.csv", sep="_")
  name <- paste(input_names[c], "data_res_cons", sep="_")
  data_all1      <<- read.csv(name_curr, header =  TRUE)
  assign(name, data_all1)
}
#The function gives a log likelihood from the given parameters
one_country_likelihood <- function(input_one){
  #Extracting the parameters
  #ESBL fitness cost
  fce1  <- (0.5+atan(input_one$vec[1])/pi)*(0.15)
  #CR fitness cost
  fcc1  <- (0.5+atan(input_one$vec[2])/pi)*(0.15)
  #Time_of treatment with cephalosporins
  tdce1 <-  0.7
  #Time_of treatment with carbapenems
  tdcc1 <-  0.7
  #Characteristical time of reaction of microbiota to treatment
  ttc1  <-  0.7
  #Additional suseptability to colonization during antibiotic treatment
  mu1   <- (0.5+atan(input_one$vec[5])/pi)*100
  #Size of the import reservoir ESBL
  ime1  <- (0.5+atan(input_one$vec[3])/pi)*10000
  #Size of the import reservoir CR
  imc1  <- (0.5+atan(input_one$vec[4])/pi)*10000
  #Proportion of people who lose the resistance due to the replacement strain by the wild rype strain
  alpha <- (0.5+atan(input_one$vec[6])/pi) * 0.5
  itc <-  1
  #Loss of the resistance due to the decolonization
  tau = 365/(3*(0.5 + alpha))
  #proportion of people colonized with Klebsiella pneumoniae in the community
  base_c <- 0.2
  #Increased proportion of development of infection in hospital
  mult <- 25
  #Horizontal gene transfer
  HGT1 <- (0.5+atan(input_one$vec[7])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
  #Transmission rate in the community
  beta_tr1 <- 1/((1-base_c)*tau)
  {
    #Loading the name of the country
    name_curr <- paste(input_one$name, "data_res_cons", sep="_")
    #Loading data for country
    data_all1      <<- get(name_curr)
    #extractiong particular parameters
    hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
    hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
    consumption_community_cephalosporins <-  c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
    consumption_hospital_cephalosporins  <-  c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
    consumption_community_carbapenems    <-  c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
    consumption_hospital_carbapenems     <-  c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
    #Loading the initiol prevalences of ESBL and CR
    prev_ESBL   <- (0.5+atan(input_one$vec[9])/pi)*100
    prev_CRE    <- (0.5+atan(input_one$vec[10])/pi)*(0.5+atan(input_one$vec[9])/pi)*100
    #Transmission level in the hospital
    beta_tr2 <- (0.5+atan(input_one$vec[8])/pi) * 50 * beta_tr1
    #Helping function which calculated the prevalence of colonized induviduals in hospital
    a1 <- SCCH(99997,1,1,1,
               beta_tr1,
               beta_tr2, 
               1/tau,
               hospital_size1/(hospital_stay1*(100000-hospital_size1)),
               1/hospital_stay1)
    col_c <- a1[2]/(a1[2]+a1[1])
    col_h <- a1[4]/(a1[4]+a1[3])
    #Generation of initial parameters
    {inits <<- c(Community1_Susceptible_UT      = (1-col_c)*(a1[1]+a1[2] - consumption_community_cephalosporins[1]*100 - consumption_community_carbapenems[1]*100),
                 Community1_Susceptible_T_CEPH  = (1-col_c)*(consumption_community_cephalosporins[1]*100),
                 Community1_Susceptible_T_CAR   = (1-col_c)*(consumption_community_carbapenems[1]*100),
                 
                 Community1_Colonized_S_UT        = col_c * (100 - prev_ESBL)/100 * (a1[1]+a1[2] - consumption_community_cephalosporins[1]*100 - consumption_community_carbapenems[1]*100),
                 Community1_Colonized_S_T_CEPH    = col_c * (100 - prev_ESBL)/100 * (consumption_community_cephalosporins[1]*100),
                 Community1_Colonized_S_T_CAR     = col_c * (100 - prev_ESBL)/100 * (consumption_community_carbapenems[1]*100),
                 
                 Community1_Colonized_ESBL_UT     = col_c * (prev_ESBL - prev_CRE)/100 * (a1[1]+a1[2] - consumption_community_cephalosporins[1]*100 - consumption_community_carbapenems[1]*100),
                 Community1_Colonized_ESBL_T_CEPH = col_c * (prev_ESBL - prev_CRE)/100 * (consumption_community_cephalosporins[1]*100),
                 Community1_Colonized_ESBL_T_CAR  = col_c * (prev_ESBL - prev_CRE)/100 * (consumption_community_carbapenems[1]*100),
                 
                 Community1_Colonized_CRE_UT      = col_c * prev_CRE/100 * (a1[1]+a1[2] - consumption_community_cephalosporins[1]*100 - consumption_community_carbapenems[1]*100),
                 Community1_Colonized_CRE_T_CEPH  = col_c * prev_CRE/100 * (consumption_community_cephalosporins[1]*100),
                 Community1_Colonized_CRE_T_CAR   = col_c * prev_CRE/100 * (consumption_community_carbapenems[1]*100),
                 
                 
                 
                 Hospital1_Susceptible_UT     = (1-col_h) * (a1[3]+a1[4] - consumption_hospital_cephalosporins[1]*100 - consumption_hospital_carbapenems[1]*100),
                 Hospital1_Susceptible_T_CEPH = (1-col_h) * consumption_hospital_cephalosporins[1]*100,
                 Hospital1_Susceptible_T_CAR  = (1-col_h) * consumption_hospital_carbapenems[1]*100,
                 
                 Hospital1_Colonized_S_UT     = col_h * (100 - prev_ESBL)/100 * (a1[3]+a1[4] - consumption_hospital_cephalosporins[1]*100 - consumption_hospital_carbapenems[1]*100),
                 Hospital1_Colonized_S_T_CEPH = col_h * (100 - prev_ESBL)/100 * consumption_hospital_cephalosporins[1]*100,
                 Hospital1_Colonized_S_T_CAR  = col_h * (100 - prev_ESBL)/100 * consumption_hospital_carbapenems[1]/1000,
                 
                 Hospital1_Colonized_ESBL_UT     = col_h * (prev_ESBL - prev_CRE)/100 * (a1[3]+a1[4] - consumption_hospital_cephalosporins[1]*100 - consumption_hospital_carbapenems[1]*100),
                 Hospital1_Colonized_ESBL_T_CEPH = col_h * (prev_ESBL - prev_CRE)/100 * consumption_hospital_cephalosporins[1]*100,
                 Hospital1_Colonized_ESBL_T_CAR  = col_h * (prev_ESBL - prev_CRE)/100 * consumption_hospital_carbapenems[1]*100,
                 
                 Hospital1_Colonized_CRE_UT     = col_h * prev_CRE/100 * (a1[3]+a1[4] - consumption_hospital_cephalosporins[1]*100 - consumption_hospital_carbapenems[1]*100),
                 Hospital1_Colonized_CRE_T_CEPH = col_h * prev_CRE/100 * consumption_hospital_cephalosporins[1]*100,
                 Hospital1_Colonized_CRE_T_CAR  = col_h * prev_CRE/100 * consumption_hospital_carbapenems[1]*100,
                 
                 Hospital1_Infected_S_T    = (100 - prev_ESBL)/100,
                 Hospital1_Infected_ESBL_T = (prev_ESBL - prev_CRE)/100,
                 Hospital1_Infected_CRE_T  = prev_CRE/100)}
    #Generation of parameter vector
    parms <- c(beta_d = beta_tr1,
               mu = mu1,
               hospitalization_rate = hospital_size1/(hospital_stay1*(100000-hospital_size1)),
               discharge_rate = 1/hospital_stay1,
               time_disease_development_community = (mult*4000*365),
               time_disease_development_hospital = (4000*365),
               time_recovery = 10,
               fitness_cost_S     = 0.0,
               fitness_cost_ESBL  = fce1,
               fitness_cost_CRE   = fce1 + fcc1,
               R_H_C = beta_tr2/beta_tr1,
               treatment_time_clearance = 10 * ttc1, 
               treatment_duration_community_CEPH = 10 * tdce1,
               treatment_duration_community_CAR  = 10 * tdcc1,
               treatment_duration_hospital_CEPH  = 10 * tdce1,
               treatment_duration_hospital_CAR   = 10 * tdcc1,
               conversion_rate = 0.0,
               influence_infected = 1,
               t_con = 1,
               P_CEPH =  1,
               P_CAR  =  1,
               clearance_rate = 1/tau,
               Import_ESBL = ime1,
               Import_CRE  = imc1,
               HGT  = HGT1,
               HGT_D = (0.5-alpha)/(0.5+alpha),
               ITC = itc)
    #running the function
    sim1 <- SCCII(inits[1],inits[2],inits[3],
                  inits[4],inits[5],inits[6],
                  inits[7],inits[8],inits[9],
                  inits[10],inits[11],inits[12],
                  inits[13],inits[14],inits[15],
                  inits[16],inits[17],inits[18],
                  inits[19],inits[20],inits[21],
                  inits[22],inits[23],inits[24],
                  inits[25],inits[26],inits[27],
                  parms[1],
                  parms[2],
                  parms[3],
                  parms[4],
                  parms[5],
                  parms[6],
                  parms[7],
                  parms[8],
                  parms[9],
                  parms[10],
                  parms[11],
                  parms[12],
                  parms[13],
                  parms[14],
                  parms[15],
                  parms[16],
                  parms[17],
                  parms[18],
                  parms[19],
                  parms[20],
                  parms[21],
                  parms[22],
                  parms[23],
                  parms[24],
                  parms[25],
                  parms[26],
                  parms[27],
                  consumption_community_cephalosporins,consumption_hospital_cephalosporins,
                  consumption_community_carbapenems,consumption_hospital_carbapenems)
    #Extracting the results of the function
    Total_Infected = ( sim1[1,]+sim1[2,]+sim1[3,])
    Total_Infected_ESBL = (sim1[2,]+sim1[3,])
    Total_Infected_CRE =  (sim1[3,])
    #creating the vectors of modelled prevalences
    M_E_R <- c(Total_Infected_ESBL[1]/Total_Infected[1],sapply(2:11, function(x) (Total_Infected_ESBL[x])/Total_Infected[x]))
    M_C_R <- c(Total_Infected_CRE[1]/Total_Infected[1],sapply(2:11, function(x) (Total_Infected_CRE[x])/Total_Infected[x]))
    
    #Extraction of the collected data from the file
    ESBL_res_num   <- as.numeric(as.character(data_all1[2:12,2]))
    CRE_res_num    <- as.numeric(as.character(data_all1[2:12,4]))
    ESBL_res_prev  <- as.numeric(as.character(data_all1[2:12,3]))
    CRE_res_prev   <- as.numeric(as.character(data_all1[2:12,5]))
    
    ESBL_res_num   <- ESBL_res_num[1:11]
    CRE_res_num    <- CRE_res_num[1:11]
    ESBL_res_prev  <- ESBL_res_prev[1:11]
    CRE_res_prev   <- CRE_res_prev[1:11]
    M_E_R <- M_E_R[1:11]
    M_C_R <- M_C_R[1:11]
    #Generation of log-likelihood based on binomial distribution
    summ_e <- sum(sign(ESBL_res_num) * dbinom(ESBL_res_prev, ESBL_res_num, M_E_R, log=TRUE), na.rm = TRUE)
    summ_c <- sum(sign(CRE_res_num) * dbinom(CRE_res_prev, CRE_res_num, M_C_R, log=TRUE), na.rm = TRUE)
    summ  <- - summ_e - summ_c
  }
  return(summ)
}
#The function gives a log likelihood from the given parameters
curr_list <- list(name = "Country name",
                  vec = c(0,0,0,0,0,
                          0,0,
                          0,0,0))
input_fc <- rep(list(curr_list),11)
all_countries_run <- function(input_all){
  for(i in 1:11){
    vec <- c(input_all[1],input_all[2],input_all[3],input_all[4],input_all[5],
             input_all[6],input_all[40],
             input_all[7],input_all[17+i],input_all[28+i])
    curr_list <-  list(name = input_names[i],
                       vec = vec)
    input_fc[[i]] <- curr_list
  }
  ALL <- mclapply(X = input_fc,FUN = one_country_likelihood, mc.cores = 11)
  summ <- 0
  for (i in order){
    summ <- summ + as.numeric(ALL[[i]])
  }
  return(summ)
}

#Then we give an example of fitting the parameters for all 11 countries together and for the "one-out" counterfactual scenario
#Output values
out_lll <- rep(0,12)
out_vvv <- rep(3,12)
output_without_one <- matrix(0, 11,40)
#Countries which are taken into account in the current simulation
basic <- matrix(0,11,10)
basic[1,] <-c(2:11) 
for (i in 2:10){
  basic[i,] <- c(1:(i-1),(i+1):11)
}
basic[11,] <-c(1:10)
as.numeric(a[12,2:41])
out_lll <- rep(0,12)
out_vvv <- rep(3,12)
order <- c(1:11)
input_func <- rep(0,40)
for (t in 1:5){
  output_func <- optim(input_func, all_countries_run, method = "BFGS", control = list(maxit=5,trace=2, REPORT =2))
  input_func <- output_func$par
}
all_countries_run(input_func)
output_without_one[12,] <- input_func
out_lll[12] <- output_func$value
out_vvv[12] <- output_func$convergence
for (i in c(1:11)){
  order <- basic[i,]
  input_func <- output_without_one[12,]
  for (t in 1:2){
    output_func <- optim(input_func, all_countries_run, method = "BFGS", control = list(maxit=70,trace=10, REPORT =2))
    input_func <- output_func$par
  }
  output_without_one[i,] <- input_func
  out_lll[i] <- output_func$value
  out_vvv[i] <- output_func$convergence
}
