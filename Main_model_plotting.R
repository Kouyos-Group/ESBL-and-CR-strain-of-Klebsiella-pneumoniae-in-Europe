library(deSolve)
library(parallel)
library(Rcpp)
library(numDeriv)
library(grid)
library(vcd)
require(graphics)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sourceCpp("SCCH.cpp")   
sourceCpp("SCCII2.cpp") 
input_names<- c("Croatia","Denmark","Finland","France", "Greece", "Hungary", "Italy", "Netherlands", "Norway","Portugal", "Sweden")
#loading results of the fit (uniform rate)
uniform_results <- read.csv("Results_original_and_one_out_uniform.csv")
uniform_results <- uniform_results[c(1:12),c(2:41)]
#loading results of the fit (unique rate)
unique_results <- read.csv("Results_original_and_one_out_unique.csv")
unique_results <- unique_results[c(1:12),c(2:41)]
#extracting results from the table
fitness_cost1 <- as.numeric(unique_results[12,])
fitness_cost2 <- as.numeric(uniform_results[12,])
fitness_cost2[c(7:17)] <- fitness_cost2[7]
input <- list(names = input_names,
              fitness_cost1 = fitness_cost1,
              fitness_cost2 = fitness_cost2)
#loading the consumption data
for (c in 1:length(input_names)){
  name_curr <- paste(input_names[c], "3rd_4th_gen.csv", sep="_")
  name <- paste(input_names[c], "data_res_cons", sep="_")
  data_all1      <<- read.csv(name_curr, header =  TRUE)
  assign(name, data_all1)
}
#Plot
plot_vars <- function(input){layout(matrix(c(1,2,2,3,3,4,4,4,5,5,5,
                                             1,2,2,3,3,4,4,4,5,5,5,
                                             
                                             1,6,6,7,7,8,8,9,9,0,0,
                                             1,6,6,7,7,8,8,9,9,0,0,
                                             1,10,10,11,11,12,12,13,13,14,14,
                                             1,10,10,11,11,12,12,13,13,14,14,
                                             0,0,0,0,0,0,0,0,0,0,0), 7, 11, byrow = TRUE))
  par(mar=c(0,0,2.5,0.1))
  
  plot.new()
  g1 = grid_legend("left", labels=c("Prevalence of the resistant strains %"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
  grid.draw(grobTree(g1, vp=viewport(x=0.025, y=0.375,angle=90)))
  y = grid_legend("center", labels=c("Year"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
  grid.draw(grobTree(y, vp=viewport(x=0.5, y=0.025,angle=0)))
  #Increased proportion of development of infection in hospital
  mult <- 25
  #Time_of treatment with cephalosporins
  tdce1 <-0.7
  #Time_of treatment with carbapenems
  tdcc1 <-0.7
  #Characteristical time of reaction of microbiota to treatment
  ttc1 <-0.7
  t <- 1
  for (c in 1:length(input_names)){
    if ((input_names[c]=="Italy")||(input_names[c]=="Greece")){
      #Extracting the parameters
      #ESBL fitness cost
      fce1  <- (0.5+atan(input$fitness_cost1[1])/pi)*(0.15)
      #CR fitness cost
      fcc1  <-(0.5+atan(input$fitness_cost1[2])/pi)*(0.15)
      #Additional suseptability to colonization during antibiotic treatment
      mu1   <- (0.5+atan(input$fitness_cost1[5])/pi)*100
      #Size of the import reservoir ESBL (In the SCCII2 function it will be devided by 100)
      ime1  <- (0.5+atan(input$fitness_cost1[3])/pi)*10000
      #Size of the import reservoir CR (In the SCCII2 function it will be devided by 100)
      imc1  <- (0.5+atan(input$fitness_cost1[4])/pi)*10000
      #Proportion of people who lose the resistance due to the replacement strain by the wild rype strain
      alpha <- (0.5+atan(input$fitness_cost1[6])/pi) * 0.5
      #Loss of the resistance due to the decolonization
      tau = 365/(3*(0.5 + alpha))
      #proportion of people colonized with Klebsiella pneumoniae in the community
      base_c <- 0.2
      #Horizontal gene transfer
      HGT1 <- (0.5+atan(input$fitness_cost1[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
      #Transmission rate in the community
      beta_tr1 <- 1/((1-base_c)*tau)
      #Loading the name of the country
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      #Loading data for country
      data_all1      <<- get(name_curr)
      #extractiong particular parameters
      hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
      hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
      consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
      consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
      consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
      consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
      #Loading the initiol prevalences of ESBL and CR
      prev_ESBL   <- (0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
      prev_CRE    <- (0.5+atan(input$fitness_cost1[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
      #Transmission level in the hospital
      beta_tr2 <- (0.5+atan(input$fitness_cost1[6+c])/pi)*50*beta_tr1
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
                 time_disease_development_community = (mult*100*365),
                 time_disease_development_hospital = (100*365),
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
                 itc = 1)
      #running the function
      { sim1 <- SCCII(inits[1],inits[2],inits[3],
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
                      consumption_community_carbapenems,consumption_hospital_carbapenems)}
      #Extracting the results of the function
      Total_Infected1 = ( sim1[1,]+sim1[2,]+sim1[3,])
      
      Total_Infected_ESBL1 = (sim1[2,]+sim1[3,])
      Total_Infected_CRE1 =  (sim1[3,])
      #Same process for the uniform transmission rate
      fce1  <- (0.5+atan(input$fitness_cost2[1])/pi)*(0.15)
      fcc1  <-(0.5+atan(input$fitness_cost2[2])/pi)*(0.15)
      mu1   <- (0.5+atan(input$fitness_cost2[5])/pi)*100
      ime1  <- (0.5+atan(input$fitness_cost2[3])/pi)*10000
      imc1  <- (0.5+atan(input$fitness_cost2[4])/pi)*10000
      alpha <- (0.5+atan(input$fitness_cost2[6])/pi) * 0.5
      tau = 365/(3*(0.5 + alpha))
      base_c <- 0.2
      HGT1 <- (0.5+atan(input$fitness_cost2[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
      beta_tr1 <- 1/((1-base_c)*tau)
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      data_all1      <<- get(name_curr)
      hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
      hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
      consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
      consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
      consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
      consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
      prev_ESBL   <- (0.5+atan(input$fitness_cost2[6+11+c])/pi)*100
      prev_CRE    <- (0.5+atan(input$fitness_cost2[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost2[6+11+c])/pi)*100
      beta_tr2 <- (0.5+atan(input$fitness_cost2[6+c])/pi)*50*beta_tr1
      a1 <- SCCH(99997,1,1,1,
                 beta_tr1,
                 beta_tr2, 
                 1/tau,
                 hospital_size1/(hospital_stay1*(100000-hospital_size1)),
                 1/hospital_stay1)
      col_c <- a1[2]/(a1[2]+a1[1])
      col_h <- a1[4]/(a1[4]+a1[3])
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
      parms <- c(beta_d = beta_tr1,
                 mu = mu1,
                 hospitalization_rate = hospital_size1/(hospital_stay1*(100000-hospital_size1)),
                 discharge_rate = 1/hospital_stay1,
                 time_disease_development_community = (mult*100*365),
                 time_disease_development_hospital = (100*365),
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
                 itc = 1)
      { sim1 <- SCCII(inits[1],inits[2],inits[3],
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
                      consumption_community_carbapenems,consumption_hospital_carbapenems)}
      Total_Infected2 = ( sim1[1,]+sim1[2,]+sim1[3,])
      
      Total_Infected_ESBL2 = (sim1[2,]+sim1[3,])
      Total_Infected_CRE2 =  (sim1[3,])
      # Generating times for plot
      times <- seq(0,3650,365)
      #Extaction data
      N_ESBL <<- as.numeric(as.character(data_all1[2:12,2]))
      N_CRE  <<- as.numeric(as.character(data_all1[2:12,4]))
      ESBL_res <- as.numeric(as.character(data_all1[2:12,3]))
      CRE_res  <- as.numeric(as.character(data_all1[2:12,5]))
      years   <- c(2005:2015)
      E_R     <- 100*ESBL_res/N_ESBL
      C_R     <- 100*CRE_res/N_CRE
      if (t==1)
      {
        par(mar=c(0,0,2.5,0.2))
        t <- t+1
        plot((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", xlim = c(2005,2015), ylim = c(0,85), type="l", lwd=2,ylab='',xlab='', xaxt = 'n',cex.axis=1.5)
        lines((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, 100 * (Total_Infected_CRE1)/(Total_Infected1), lty=1, col = "orange", lwd=2)
        lines((times/365)+2005, 100 * (Total_Infected_ESBL2)/(Total_Infected2), lty=3, col = "blue", lwd=2)
        lines((times/365)+2005, 100 * (Total_Infected_CRE2)/(Total_Infected2), lty=3, col = "orange", lwd=2)
      }
      else
      {
        par(mar=c(0,0,2.5,0.2))
        t <- t+1
        plot((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", xlim = c(2005,2015), ylim = c(0,85), type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
        lines((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, 100 * (Total_Infected_CRE1)/(Total_Infected1), lty=1, col = "orange", lwd=2)
        lines((times/365)+2005, 100 * (Total_Infected_ESBL2)/(Total_Infected2), lty=3, col = "blue", lwd=2)
        lines((times/365)+2005, 100 * (Total_Infected_CRE2)/(Total_Infected2), lty=3, col = "orange", lwd=2)
      }
      title(main =paste(input_names[c]),cex.main=2)
      years <- c(2005:2015)
      lines(years, E_R, type ="o", lty=0, col = "blue", lwd=2)
      lines(years, C_R, type ="o", lty=0, col = "orange", lwd=2)  
    }}
  plot.new()
  legend("topleft",title = "ESBL strain", legend=c("Variable rate","Uniform rate","Reported data"),
         col=c("blue"), lty=c(1,3,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
         box.lty=0)
  plot.new()
  legend("topleft",title = "CRK strain", legend=c("Variable rate","Uniform rate","Reported data"),
         col=c("orange"), lty=c(1,3,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
         box.lty=0)
  t <- 1
  for (c in 1:length(input_names)){
    if ((input_names[c]=="Croatia")||(input_names[c]=="Hungary")||(input_names[c]=="France")||(input_names[c]=="Portugal"))
    {
    #Extracting the parameters
    #ESBL fitness cost
    fce1  <- (0.5+atan(input$fitness_cost1[1])/pi)*(0.15)
    #CR fitness cost
    fcc1  <-(0.5+atan(input$fitness_cost1[2])/pi)*(0.15)
    #Additional suseptability to colonization during antibiotic treatment
    mu1   <- (0.5+atan(input$fitness_cost1[5])/pi)*100
    #Size of the import reservoir ESBL (In the SCCII2 function it will be devided by 100)
    ime1  <- (0.5+atan(input$fitness_cost1[3])/pi)*10000
    #Size of the import reservoir CR (In the SCCII2 function it will be devided by 100)
    imc1  <- (0.5+atan(input$fitness_cost1[4])/pi)*10000
    #Proportion of people who lose the resistance due to the replacement strain by the wild rype strain
    alpha <- (0.5+atan(input$fitness_cost1[6])/pi) * 0.5
    #Loss of the resistance due to the decolonization
    tau = 365/(3*(0.5 + alpha))
    #proportion of people colonized with Klebsiella pneumoniae in the community
    base_c <- 0.2
    #Horizontal gene transfer
    HGT1 <- (0.5+atan(input$fitness_cost1[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
    #Transmission rate in the community
    beta_tr1 <- 1/((1-base_c)*tau)
    #Loading the name of the country
    name_curr <- paste(input$names[c], "data_res_cons", sep="_")
    #Loading data for country
    data_all1      <<- get(name_curr)
    #extractiong particular parameters
    hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
    hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
    consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
    consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
    consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
    consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
    #Loading the initiol prevalences of ESBL and CR
    prev_ESBL   <- (0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
    prev_CRE    <- (0.5+atan(input$fitness_cost1[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
    #Transmission level in the hospital
    beta_tr2 <- (0.5+atan(input$fitness_cost1[6+c])/pi)*50*beta_tr1
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
               time_disease_development_community = (mult*100*365),
               time_disease_development_hospital = (100*365),
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
               itc = 1)
    #running the function
    { sim1 <- SCCII(inits[1],inits[2],inits[3],
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
                    consumption_community_carbapenems,consumption_hospital_carbapenems)}
    #Extracting the results of the function
    Total_Infected1 = ( sim1[1,]+sim1[2,]+sim1[3,])
    
    Total_Infected_ESBL1 = (sim1[2,]+sim1[3,])
    Total_Infected_CRE1 =  (sim1[3,])
    
    summ <- rep(0,length(input$names))
    #Same process for the uniform transmission rate
    fce1  <- (0.5+atan(input$fitness_cost2[1])/pi)*(0.15)
    fcc1  <-(0.5+atan(input$fitness_cost2[2])/pi)*(0.15)
    mu1   <- (0.5+atan(input$fitness_cost2[5])/pi)*100
    ime1  <- (0.5+atan(input$fitness_cost2[3])/pi)*10000
    imc1  <- (0.5+atan(input$fitness_cost2[4])/pi)*10000
    alpha <- (0.5+atan(input$fitness_cost2[6])/pi) * 0.5
    tau = 365/(3*(0.5 + alpha))
    base_c <- 0.2
    HGT1 <- (0.5+atan(input$fitness_cost2[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
    beta_tr1 <- 1/((1-base_c)*tau)
    name_curr <- paste(input$names[c], "data_res_cons", sep="_")
    data_all1      <<- get(name_curr)
    hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
    hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
    consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
    consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
    consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
    consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
    prev_ESBL   <- (0.5+atan(input$fitness_cost2[6+11+c])/pi)*100
    prev_CRE    <- (0.5+atan(input$fitness_cost2[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost2[6+11+c])/pi)*100
    beta_tr2 <- (0.5+atan(input$fitness_cost2[6+c])/pi)*50*beta_tr1
    a1 <- SCCH(99997,1,1,1,
               beta_tr1,
               beta_tr2, 
               1/tau,
               hospital_size1/(hospital_stay1*(100000-hospital_size1)),
               1/hospital_stay1)
    col_c <- a1[2]/(a1[2]+a1[1])
    col_h <- a1[4]/(a1[4]+a1[3])
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
    parms <- c(beta_d = beta_tr1,
               mu = mu1,
               hospitalization_rate = hospital_size1/(hospital_stay1*(100000-hospital_size1)),
               discharge_rate = 1/hospital_stay1,
               time_disease_development_community = (mult*100*365),
               time_disease_development_hospital = (100*365),
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
               itc = 1)
    { sim1 <- SCCII(inits[1],inits[2],inits[3],
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
                    consumption_community_carbapenems,consumption_hospital_carbapenems)}
    Total_Infected2 = ( sim1[1,]+sim1[2,]+sim1[3,])
    
    Total_Infected_ESBL2 = (sim1[2,]+sim1[3,])
    Total_Infected_CRE2 =  (sim1[3,])
    # Generating times for plot
    times <- seq(0,3650,365)
    #Extaction data
    N_ESBL <<- as.numeric(as.character(data_all1[2:12,2]))
    N_CRE  <<- as.numeric(as.character(data_all1[2:12,4]))
    ESBL_res <- as.numeric(as.character(data_all1[2:12,3]))
    CRE_res  <- as.numeric(as.character(data_all1[2:12,5]))
    years   <- c(2005:2015)
    E_R     <- 100*ESBL_res/N_ESBL
    C_R     <- 100*CRE_res/N_CRE
    if (t==1){
      par(mar=c(0,0,2.5,0.2))
      t <- t+1
      plot((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", xlim = c(2005,2015), ylim = c(0,65), type="l", lwd=2,ylab='',xlab='',cex.axis=1.5, xaxt = 'n')
      lines((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE1)/(Total_Infected1), lty=1, col = "orange", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_ESBL2)/(Total_Infected2), lty=3, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE2)/(Total_Infected2), lty=3, col = "orange", lwd=2)
    }
    else{
      par(mar=c(0,0,2.5,0.2))
      t <- t+1
      plot((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", xlim = c(2005,2015), ylim = c(0,65), type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
      lines((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE1)/(Total_Infected1), lty=1, col = "orange", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_ESBL2)/(Total_Infected2), lty=3, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE2)/(Total_Infected2), lty=3, col = "orange", lwd=2)
    }
    title(main =paste(input_names[c]),cex.main=2)
    years <- c(2005:2015)
    lines(years, E_R, type ="o", lty=0, col = "blue", lwd=2)
    lines(years, C_R, type ="o", lty=0, col = "orange", lwd=2)  
    }
  }
  t <- 1
  for (c in 1:length(input_names)){
    if ((input_names[c]=="Netherlands")||(input_names[c]=="Norway")||(input_names[c]=="Denmark")||(input_names[c]=="Finland")||(input_names[c]=="Sweden"))
    {#Extracting the parameters
      #ESBL fitness cost
      fce1  <- (0.5+atan(input$fitness_cost1[1])/pi)*(0.15)
      #CR fitness cost
      fcc1  <-(0.5+atan(input$fitness_cost1[2])/pi)*(0.15)
      #Additional suseptability to colonization during antibiotic treatment
      mu1   <- (0.5+atan(input$fitness_cost1[5])/pi)*100
      #Size of the import reservoir ESBL (In the SCCII2 function it will be devided by 100)
      ime1  <- (0.5+atan(input$fitness_cost1[3])/pi)*10000
      #Size of the import reservoir CR (In the SCCII2 function it will be devided by 100)
      imc1  <- (0.5+atan(input$fitness_cost1[4])/pi)*10000
      #Proportion of people who lose the resistance due to the replacement strain by the wild rype strain
      alpha <- (0.5+atan(input$fitness_cost1[6])/pi) * 0.5
      #Loss of the resistance due to the decolonization
      tau = 365/(3*(0.5 + alpha))
      #proportion of people colonized with Klebsiella pneumoniae in the community
      base_c <- 0.2
      #Horizontal gene transfer
      HGT1 <- (0.5+atan(input$fitness_cost1[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
      #Transmission rate in the community
      beta_tr1 <- 1/((1-base_c)*tau)
      #Loading the name of the country
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      #Loading data for country
      data_all1      <<- get(name_curr)
      #extractiong particular parameters
      hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
      hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
      consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
      consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
      consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
      consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
      #Loading the initiol prevalences of ESBL and CR
      prev_ESBL   <- (0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
      prev_CRE    <- (0.5+atan(input$fitness_cost1[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
      #Transmission level in the hospital
      beta_tr2 <- (0.5+atan(input$fitness_cost1[6+c])/pi)*50*beta_tr1
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
                 time_disease_development_community = (mult*100*365),
                 time_disease_development_hospital = (100*365),
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
                 itc = 1)
      #running the function
      { sim1 <- SCCII(inits[1],inits[2],inits[3],
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
                      consumption_community_carbapenems,consumption_hospital_carbapenems)}
      #Extracting the results of the function
      Total_Infected1 = ( sim1[1,]+sim1[2,]+sim1[3,])
      
      Total_Infected_ESBL1 = (sim1[2,]+sim1[3,])
      Total_Infected_CRE1 =  (sim1[3,])
      
      summ <- rep(0,length(input$names))
      #Same process for the uniform transmission rate
      fce1  <- (0.5+atan(input$fitness_cost2[1])/pi)*(0.15)
      fcc1  <-(0.5+atan(input$fitness_cost2[2])/pi)*(0.15)
      mu1   <- (0.5+atan(input$fitness_cost2[5])/pi)*100
      ime1  <- (0.5+atan(input$fitness_cost2[3])/pi)*10000
      imc1  <- (0.5+atan(input$fitness_cost2[4])/pi)*10000
      alpha <- (0.5+atan(input$fitness_cost2[6])/pi) * 0.5
      tau = 365/(3*(0.5 + alpha))
      base_c <- 0.2
      HGT1 <- (0.5+atan(input$fitness_cost2[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
      beta_tr1 <- 1/((1-base_c)*tau)
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      data_all1      <<- get(name_curr)
      hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
      hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
      consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
      consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
      consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
      consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
      prev_ESBL   <- (0.5+atan(input$fitness_cost2[6+11+c])/pi)*100
      prev_CRE    <- (0.5+atan(input$fitness_cost2[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost2[6+11+c])/pi)*100
      beta_tr2 <- (0.5+atan(input$fitness_cost2[6+c])/pi)*50*beta_tr1
      a1 <- SCCH(99997,1,1,1,
                 beta_tr1,
                 beta_tr2, 
                 1/tau,
                 hospital_size1/(hospital_stay1*(100000-hospital_size1)),
                 1/hospital_stay1)
      col_c <- a1[2]/(a1[2]+a1[1])
      col_h <- a1[4]/(a1[4]+a1[3])
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
      parms <- c(beta_d = beta_tr1,
                 mu = mu1,
                 hospitalization_rate = hospital_size1/(hospital_stay1*(100000-hospital_size1)),
                 discharge_rate = 1/hospital_stay1,
                 time_disease_development_community = (mult*100*365),
                 time_disease_development_hospital = (100*365),
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
                 itc = 1)
      { sim1 <- SCCII(inits[1],inits[2],inits[3],
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
                      consumption_community_carbapenems,consumption_hospital_carbapenems)}
      Total_Infected2 = ( sim1[1,]+sim1[2,]+sim1[3,])
      
      Total_Infected_ESBL2 = (sim1[2,]+sim1[3,])
      Total_Infected_CRE2 =  (sim1[3,])
      # Generating times for plot
      times <- seq(0,3650,365)
      #Extaction data
      N_ESBL <<- as.numeric(as.character(data_all1[2:12,2]))
      N_CRE  <<- as.numeric(as.character(data_all1[2:12,4]))
      ESBL_res <- as.numeric(as.character(data_all1[2:12,3]))
      CRE_res  <- as.numeric(as.character(data_all1[2:12,5]))
      years   <- c(2005:2015)
      E_R     <- 100*ESBL_res/N_ESBL
      C_R     <- 100*CRE_res/N_CRE
    if (t==1){
      par(mar=c(0,0,2.5,0.2))
      t <- t+1
      plot((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='')
      lines((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE1)/(Total_Infected1), lty=1, col = "orange", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_ESBL2)/(Total_Infected2), lty=3, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE2)/(Total_Infected2), lty=3, col = "orange", lwd=2)
    }
    else{
      par(mar=c(0,0,2.5,0.2))
      t <- t+1
      plot((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', yaxt = 'n')
      lines((times/365)+2005, 100 * (Total_Infected_ESBL1)/(Total_Infected1), lty=1, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE1)/(Total_Infected1), lty=1, col = "orange", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_ESBL2)/(Total_Infected2), lty=3, col = "blue", lwd=2)
      lines((times/365)+2005, 100 * (Total_Infected_CRE2)/(Total_Infected2), lty=3, col = "orange", lwd=2)
    }
    title(main =paste(input_names[c]),cex.main=2)
    years <- c(2005:2015)
    lines(years, E_R, type ="o", lty=0, col = "blue", lwd=2)
    lines(years, C_R, type ="o", lty=0, col = "orange", lwd=2)  
    }
  }
}

plot_vars(input)

