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
#Loading the original fit
orig_ESBL <- matrix(nrow=11,ncol=11)
orig_CRK <- matrix(nrow=11,ncol=11)
mult <- 25
tdce1 <- 0.7
tdcc1 <- 0.7
ttc1  <- 0.7
t <- 1
fitness_cost1 <- as.numeric(unique_results[12,])
fitness_cost2 <- as.numeric(uniform_results[12,])
fitness_cost2[7:17] <- fitness_cost2[7]
input <- list(names = input_names,
              fitness_cost1 = fitness_cost2,
              fitness_cost2 = fitness_cost2)
#calculating the original fit
for (c in 1:11){
  {
    #See the description in Main_model.R
    fce1  <- (0.5+atan(input$fitness_cost1[1])/pi)*(0.15)
    fcc1  <-(0.5+atan(input$fitness_cost1[2])/pi)*(0.15)
    mu1   <- (0.5+atan(input$fitness_cost1[5])/pi)*100
    ime1  <- (0.5+atan(input$fitness_cost1[3])/pi)*10000
    imc1  <- (0.5+atan(input$fitness_cost1[4])/pi)*10000
    alpha <- (0.5+atan(input$fitness_cost1[6])/pi) * 0.5
    tau = 365/(3*(0.5 + alpha))
    base_c <- 0.2
    HGT1 <- (0.5+atan(input$fitness_cost1[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
    beta_tr1 <- 1/((1-base_c)*tau)
    
    name_curr <- paste(input$names[c], "data_res_cons", sep="_")
    data_all1      <<- get(name_curr)
    hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
    hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
    consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
    consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
    consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
    consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
    prev_ESBL   <- (0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
    prev_CRE    <- (0.5+atan(input$fitness_cost1[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
    beta_tr2 <- (0.5+atan(input$fitness_cost1[6+c])/pi)*50*beta_tr1
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
    Total_Infected1 = ( sim1[1,]+sim1[2,]+sim1[3,])
    
    Total_Infected_ESBL1 = (sim1[2,]+sim1[3,])
    Total_Infected_CRE1 =  (sim1[3,])
    M_E_R <- c(Total_Infected_ESBL1[1]/Total_Infected1[1],sapply(2:11, function(x) (Total_Infected_ESBL1[x])/Total_Infected1[x]))
    M_C_R <- c(Total_Infected_CRE1[1]/Total_Infected1[1],sapply(2:11, function(x) (Total_Infected_CRE1[x])/Total_Infected1[x]))
    orig_ESBL[c,] <- M_E_R
    orig_CRK[c,]  <- M_C_R
  }}
orig_ESBL <- 100*orig_ESBL
orig_CRK <- 100*orig_CRK
#setting boundaries of the possible trajectories to the initial values
min_val_ESBL <- orig_ESBL
min_val_CRK  <- orig_CRK
max_val_ESBL <- orig_ESBL
max_val_CRK  <- orig_CRK
#plotting original plot
{layout(matrix(c(1,2,2,3,3,4,4,4,5,5,5,
                 1,2,2,3,3,4,4,4,5,5,5,
                 
                 1,6,6,7,7,8,8,9,9,0,0,
                 1,6,6,7,7,8,8,9,9,0,0,
                 1,10,10,11,11,12,12,13,13,14,14,
                 1,10,10,11,11,12,12,13,13,14,14,
                 0,0,0,0,0,0,0,0,0,0,0), 7, 11, byrow = TRUE))
  par(mar=c(0,0,2.5,0))
  
  plot.new()
  g1 = grid_legend("left", labels=c("Prevalence of the resitant strains %"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
  grid.draw(grobTree(g1, vp=viewport(x=0.025, y=0.375,angle=90)))
  y = grid_legend("center", labels=c("Year"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
  grid.draw(grobTree(y, vp=viewport(x=0.5, y=0.025,angle=0)))
  t <- 1
  for (c in 1:length(input_names)){
    if ((input_names[c]=="Italy")||(input_names[c]=="Greece")){
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      data_all1      <<- get(name_curr)
      
      times <- seq(0,3650,365)
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
        plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,85),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n')
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
        lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
      }
      else
      {
        par(mar=c(0,0,2.5,0.2))
        t <- t+1
        plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,85),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
        lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
      }
      title(main =paste(input_names[c]),cex.main=2)
      years <- c(2005:2015)
      lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
      lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)   
    }}
  par(mar=c(0,0,0,0))
  plot.new()
  legend("topleft",title = "ESBL strain", legend=c("Original fit","Possible \n trajectories","Reported data"),
         col=c("blue"), fill=c("white", rgb(0, 0, 1,0.2), "white"),density=c(0, NA, 0),border=c("white", "blue", "white"), lty=c(1,0,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
         box.lty=0)
  plot.new()
  legend("topleft",title = "CRK strain", legend=c("Original fit","Possible \n trajectories","Reported data"),
         col=c("orange"), fill=c("white", rgb(255/256, 165/256, 0,0.2), "white"),density=c(0, NA, 0),border=c("white", "orange", "white"), lty=c(1,0,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
         box.lty=0)
  par(mar=c(0,0,2.5,0))
  t <- 1
  for (c in 1:length(input_names)){
    if ((input_names[c]=="Croatia")||(input_names[c]=="Hungary")||(input_names[c]=="France")||(input_names[c]=="Portugal"))
    {
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      data_all1      <<- get(name_curr)
      
      times <- seq(0,3650,365)
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
        plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,65),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n')
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
        lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
      }
      else
      {
        par(mar=c(0,0,2.5,0.2))
        t <- t+1
        plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,65),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
        lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
      }
      title(main =paste(input_names[c]),cex.main=2)
      years <- c(2005:2015)
      lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
      lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)    
    }
  }
  t <- 1
  for (c in 1:length(input_names)){
    if ((input_names[c]=="Netherlands")||(input_names[c]=="Norway")||(input_names[c]=="Denmark")||(input_names[c]=="Finland")||(input_names[c]=="Sweden"))
    {
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      data_all1      <<- get(name_curr)
      
      times <- seq(0,3650,365)
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
        plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='')
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
        lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
      }
      else
      {
        par(mar=c(0,0,2.5,0.2))
        t <- t+1
        plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', yaxt = 'n')
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
        polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
        lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
        lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
      }
      title(main =paste(input_names[c]),cex.main=2)
      years <- c(2005:2015)
      lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
      lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)    
    }
  }
}
#Narrow range
#Loading the generated latin hypercube
V <- read.csv("LH1.csv")
V <- V[1:20,2:6]
loadnprint <- function(name){
  min_val_ESBL1 <- orig_ESBL
  min_val_CRK1  <- orig_CRK
  max_val_ESBL1 <- orig_ESBL
  max_val_CRK1  <- orig_CRK
  #reading the file
  a <- read.csv(name)
  #creating a table of fitte parameters
  OUTPUT <- a[1:20,2:41]
  OUTPUT[7:17] <- OUTPUT[7]
  #Running the model
  for (i in c(1:20)){
    #loading a point from the hypercube
    V1 <- as.numeric(V[i,])
    #loading a vector fitted to this point
    fitness_cost1 <- as.numeric(OUTPUT[i,])
    input <- list(names = input_names,
                  fitness_cost1 = fitness_cost1,
                  fitness_cost2 = fitness_cost2)
    #Running for all 11 countries
    for (c in 1:11){
      #Running for all 11 countries
      fce1  <- (0.5+atan(input$fitness_cost1[1])/pi)*(0.15)
      fcc1  <-(0.5+atan(input$fitness_cost1[2])/pi)*(0.15)
      #adjusted parameters
      #time of treatment
      tdce1 <-  0.3+0.7*(V1[1])
      tdcc1 <-  0.3+0.7*(V1[1])
      #time of clearance
      ttc1  <-  0.1+0.9*(V1[2])
      mu1   <- (0.5+atan(input$fitness_cost1[5])/pi)*100
      ime1  <- (0.5+atan(input$fitness_cost1[3])/pi)*10000
      imc1  <- (0.5+atan(input$fitness_cost1[4])/pi)*10000
      alpha <- (0.5+atan(input$fitness_cost1[6])/pi) * 0.5
      tau = (365/(2*(0.5 + alpha)) - 365/(4*(0.5 + alpha))) * V1[5] + 365/(4*(0.5 + alpha))
      #Basic colonization rate
      base_c <- 0.15 + 0.1*(V1[3])
      #Multiplicator of develompnent of the infection in the hospital(the higher mult the higher the hospital onset)
      mult <- 1+199*(V1[4])
      HGT1 <- (0.5+atan(input$fitness_cost1[40])/pi) * (1-base_c)/base_c * (1 - 0.5 - alpha)/(0.5 + alpha)
      beta_tr1 <- 1/((1-base_c)*tau)
      #loading the consumption
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      data_all1      <<- get(name_curr)
      hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
      hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
      consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
      consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
      consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
      consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
      prev_ESBL   <- (0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
      prev_CRE    <- (0.5+atan(input$fitness_cost1[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
       beta_tr2 <- (0.5+atan(input$fitness_cost1[6+c])/pi)*50*beta_tr1
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
      #exctraction of the simulation
      Total_Infected1 = ( sim1[1,]+sim1[2,]+sim1[3,])
      
      Total_Infected_ESBL1 = (sim1[2,]+sim1[3,])
      Total_Infected_CRE1 =  (sim1[3,])
      M_E_R <- 100*c(Total_Infected_ESBL1[1]/Total_Infected1[1],sapply(2:11, function(x) (Total_Infected_ESBL1[x])/Total_Infected1[x]))
      M_C_R <- 100*c(Total_Infected_CRE1[1]/Total_Infected1[1],sapply(2:11, function(x) (Total_Infected_CRE1[x])/Total_Infected1[x]))
      #reevaluation of the possible trajectories
      for (j in 1:11){
        min_val_ESBL1[c,j] <- min(M_E_R[j],min_val_ESBL1[c,j])
        min_val_CRK1[c,j]  <- min(M_C_R[j],min_val_CRK1[c,j])
        max_val_ESBL1[c,j] <- max(M_E_R[j],max_val_ESBL1[c,j])
        max_val_CRK1[c,j]  <- max(M_C_R[j],max_val_CRK1[c,j])
      }
      }
  }
  min_val_ESBL <<- min_val_ESBL1
  min_val_CRK  <<- min_val_CRK1
  max_val_ESBL <<- max_val_ESBL1
  max_val_CRK  <<- max_val_CRK1
  #plotting
  {layout(matrix(c(1,2,2,3,3,4,4,4,5,5,5,
                   1,2,2,3,3,4,4,4,5,5,5,
                   
                   1,6,6,7,7,8,8,9,9,0,0,
                   1,6,6,7,7,8,8,9,9,0,0,
                   1,10,10,11,11,12,12,13,13,14,14,
                   1,10,10,11,11,12,12,13,13,14,14,
                   0,0,0,0,0,0,0,0,0,0,0), 7, 11, byrow = TRUE))
    par(mar=c(0,0,2.5,0))
    
    plot.new()
    g1 = grid_legend("left", labels=c("Prevalence of the resitant strains %"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
    grid.draw(grobTree(g1, vp=viewport(x=0.025, y=0.375,angle=90)))
    y = grid_legend("center", labels=c("Year"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
    grid.draw(grobTree(y, vp=viewport(x=0.5, y=0.025,angle=0)))
    t <- 1
    for (c in 1:length(input_names)){
      if ((input_names[c]=="Italy")||(input_names[c]=="Greece")){
        name_curr <- paste(input$names[c], "data_res_cons", sep="_")
        data_all1      <<- get(name_curr)
        
        times <- seq(0,3650,365)
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
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,85),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        else
        {
          par(mar=c(0,0,2.5,0.2))
          t <- t+1
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,85),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        title(main =paste(input_names[c]),cex.main=2)
        years <- c(2005:2015)
        lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
        lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)   
      }}
    par(mar=c(0,0,0,0))
    plot.new()
    legend("topleft",title = "ESBL strain", legend=c("Original fit","Possible \n trajectories","Reported data"),
           col=c("blue"), fill=c("white", rgb(0, 0, 1,0.2), "white"),density=c(0, NA, 0),border=c("white", "blue", "white"), lty=c(1,0,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
           box.lty=0)
    plot.new()
    legend("topleft",title = "CRK strain", legend=c("Original fit","Possible \n trajectories","Reported data"),
           col=c("orange"), fill=c("white", rgb(255/256, 165/256, 0,0.2), "white"),density=c(0, NA, 0),border=c("white", "orange", "white"), lty=c(1,0,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
           box.lty=0)
    par(mar=c(0,0,2.5,0))
    t <- 1
    for (c in 1:length(input_names)){
      if ((input_names[c]=="Croatia")||(input_names[c]=="Hungary")||(input_names[c]=="France")||(input_names[c]=="Portugal"))
      {
        name_curr <- paste(input$names[c], "data_res_cons", sep="_")
        data_all1      <<- get(name_curr)
        
        times <- seq(0,3650,365)
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
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,65),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        else
        {
          par(mar=c(0,0,2.5,0.2))
          t <- t+1
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,65),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        title(main =paste(input_names[c]),cex.main=2)
        years <- c(2005:2015)
        lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
        lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)    
      }
    }
    t <- 1
    for (c in 1:length(input_names)){
      if ((input_names[c]=="Netherlands")||(input_names[c]=="Norway")||(input_names[c]=="Denmark")||(input_names[c]=="Finland")||(input_names[c]=="Sweden"))
      {
        name_curr <- paste(input$names[c], "data_res_cons", sep="_")
        data_all1      <<- get(name_curr)
        
        times <- seq(0,3650,365)
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
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        else
        {
          par(mar=c(0,0,2.5,0.2))
          t <- t+1
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', yaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        title(main =paste(input_names[c]),cex.main=2)
        years <- c(2005:2015)
        lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
        lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)    
      }
    }
  }
}

start_time <- Sys.time()
name <- 'Parameters_narrow_range_uniform.csv'
loadnprint(name)
end_time <- Sys.time()
end_time - start_time

V <- read.csv("LH2.csv")
V <- V[1:20,2:6]
loadnprint2 <- function(name){
  a <- read.csv(name)
  OUTPUT <- a[1:20,2:41]
  OUTPUT[7:17] <- OUTPUT[7]
  for (i in 1:20){
    V1 <- as.numeric(V[i,])
    fitness_cost1 <- as.numeric(OUTPUT[i,])
    input <- list(names = input_names,
                  fitness_cost1 = fitness_cost1,
                  fitness_cost2 = fitness_cost2)
    for (c in 1:11){
      summ <- rep(0,length(input$names))
      fce1  <- (0.5+atan(input$fitness_cost1[1])/pi)*(0.15)
      fcc1  <-(0.5+atan(input$fitness_cost1[2])/pi)*(0.15)
      tdce1 <-  0.3+1.7*(V1[1])
      tdcc1 <-  0.3+1.7*(V1[1])
      ttc1  <-  0.1+1.9*(V1[2])
      mu1   <- (0.5+atan(input$fitness_cost1[5])/pi)*100
      ime1  <- (0.5+atan(input$fitness_cost1[3])/pi)*10000
      imc1  <- (0.5+atan(input$fitness_cost1[4])/pi)*10000
      alpha <- (0.5+atan(input$fitness_cost1[6])/pi) * 0.5
      tau = (365/(2*(0.5 + alpha)) - 365/(4*(0.5 + alpha))) * V1[5] + 365/(4*(0.5 + alpha))
      base_c <- 0.15 + 0.1*(V1[3])
      mult <- 1+199*(V1[4])
      HGT1 <- (0.5+atan(input$fitness_cost1[40])/pi) * (1-base_c)/base_c * (1 - 0.5 - alpha)/(0.5 + alpha)
      beta_tr1 <- 1/((1-base_c)*tau)
      
      name_curr <- paste(input$names[c], "data_res_cons", sep="_")
      data_all1      <<- get(name_curr)
      hospital_size1 <<- as.numeric(as.character(data_all1[15,2]))
      hospital_stay1 <<- as.numeric(as.character(data_all1[16,2]))
      consumption_community_cephalosporins <- c(as.numeric(as.character(data_all1[2:11,6])),as.numeric(as.character(data_all1[11,6])))
      consumption_hospital_cephalosporins  <- c(as.numeric(as.character(data_all1[2:11,7])),as.numeric(as.character(data_all1[11,7])))
      consumption_community_carbapenems    <- c(as.numeric(as.character(data_all1[2:11,8])),as.numeric(as.character(data_all1[11,8])))
      consumption_hospital_carbapenems     <- c(as.numeric(as.character(data_all1[2:11,9])),as.numeric(as.character(data_all1[11,9])))
      prev_ESBL   <- (0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
      prev_CRE    <- (0.5+atan(input$fitness_cost1[6+11+11+c])/pi)*(0.5+atan(input$fitness_cost1[6+11+c])/pi)*100
      beta_tr2 <- (0.5+atan(input$fitness_cost1[6+c])/pi)*50*beta_tr1
      
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
      Total_Infected1 = ( sim1[1,]+sim1[2,]+sim1[3,])
      
      Total_Infected_ESBL1 = (sim1[2,]+sim1[3,])
      Total_Infected_CRE1 =  (sim1[3,])
      M_E_R <- 100*c(Total_Infected_ESBL1[1]/Total_Infected1[1],sapply(2:11, function(x) (Total_Infected_ESBL1[x])/Total_Infected1[x]))
      M_C_R <- 100*c(Total_Infected_CRE1[1]/Total_Infected1[1],sapply(2:11, function(x) (Total_Infected_CRE1[x])/Total_Infected1[x]))
      for (j in 1:11){
        min_val_ESBL[c,j] <<- min(M_E_R[j],min_val_ESBL[c,j])
        min_val_CRK[c,j]  <<- min(M_C_R[j],min_val_CRK[c,j])
        max_val_ESBL[c,j] <<- max(M_E_R[j],max_val_ESBL[c,j])
        max_val_CRK[c,j]  <<- max(M_C_R[j],max_val_CRK[c,j])
      }}
  }
  #plotting
  {layout(matrix(c(1,2,2,3,3,4,4,4,5,5,5,
                   1,2,2,3,3,4,4,4,5,5,5,
                   
                   1,6,6,7,7,8,8,9,9,0,0,
                   1,6,6,7,7,8,8,9,9,0,0,
                   1,10,10,11,11,12,12,13,13,14,14,
                   1,10,10,11,11,12,12,13,13,14,14,
                   0,0,0,0,0,0,0,0,0,0,0), 7, 11, byrow = TRUE))
    par(mar=c(0,0,2.5,0))
    
    plot.new()
    g1 = grid_legend("left", labels=c("Prevalence of the resitant strains %"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
    grid.draw(grobTree(g1, vp=viewport(x=0.025, y=0.375,angle=90)))
    y = grid_legend("center", labels=c("Year"),gp_labels = gpar(fontsize=20), frame = FALSE, draw=FALSE)
    grid.draw(grobTree(y, vp=viewport(x=0.5, y=0.025,angle=0)))
    t <- 1
    for (c in 1:length(input_names)){
      if ((input_names[c]=="Italy")||(input_names[c]=="Greece")){
        name_curr <- paste(input$names[c], "data_res_cons", sep="_")
        data_all1      <<- get(name_curr)
        
        times <- seq(0,3650,365)
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
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,85),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        else
        {
          par(mar=c(0,0,2.5,0.2))
          t <- t+1
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,85),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        title(main =paste(input_names[c]),cex.main=2)
        years <- c(2005:2015)
        lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
        lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)   
      }}
    par(mar=c(0,0,0,0))
    plot.new()
    legend("topleft",title = "ESBL strain", legend=c("Original fit","Possible \n trajectories","Reported data"),
           col=c("blue"), fill=c("white", rgb(0, 0, 1,0.2), "white"),density=c(0, NA, 0),border=c("white", "blue", "white"), lty=c(1,0,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
           box.lty=0)
    plot.new()
    legend("topleft",title = "CRK strain", legend=c("Original fit","Possible \n trajectories","Reported data"),
           col=c("orange"), fill=c("white", rgb(255/256, 165/256, 0,0.2), "white"),density=c(0, NA, 0),border=c("white", "orange", "white"), lty=c(1,0,0), cex=2, lwd = c(2,2,2),pch = c(26,26,1),pt.cex=1,
           box.lty=0)
    par(mar=c(0,0,2.5,0))
    t <- 1
    for (c in 1:length(input_names)){
      if ((input_names[c]=="Croatia")||(input_names[c]=="Hungary")||(input_names[c]=="France")||(input_names[c]=="Portugal"))
      {
        name_curr <- paste(input$names[c], "data_res_cons", sep="_")
        data_all1      <<- get(name_curr)
        
        times <- seq(0,3650,365)
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
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,65),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        else
        {
          par(mar=c(0,0,2.5,0.2))
          t <- t+1
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,65),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', xaxt = 'n', yaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        title(main =paste(input_names[c]),cex.main=2)
        years <- c(2005:2015)
        lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
        lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)    
      }
    }
    t <- 1
    for (c in 1:length(input_names)){
      if ((input_names[c]=="Netherlands")||(input_names[c]=="Norway")||(input_names[c]=="Denmark")||(input_names[c]=="Finland")||(input_names[c]=="Sweden"))
      {
        name_curr <- paste(input$names[c], "data_res_cons", sep="_")
        data_all1      <<- get(name_curr)
        
        times <- seq(0,3650,365)
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
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        else
        {
          par(mar=c(0,0,2.5,0.2))
          t <- t+1
          plot((times/365)+2005, 100 * M_E_R, lty=1, col = "white", xlim = c(2005,2015), ylim = c(0,15),cex.axis=1.5, type="l", lwd=2,ylab='',xlab='', yaxt = 'n')
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_ESBL[c,] ,rev(min_val_ESBL[c,])),
                  border = rgb(0, 0, 1,0.1), col = rgb(0, 0, 1,0.2) )
          polygon(c((times/365)+2005, rev((times/365)+2005)), c(max_val_CRK[c,] ,rev(min_val_CRK[c,])),
                  border = rgb(255/256, 165/256, 0,0.2) , col = rgb(255/256, 165/256, 0,0.2) )
          lines((times/365)+2005, orig_ESBL[c,], lty=1, col = "blue", lwd=2)
          lines((times/365)+2005, orig_CRK[c,], lty=1, col = "orange", lwd=2)
        }
        title(main =paste(input_names[c]),cex.main=2)
        years <- c(2005:2015)
        lines(years, E_R, type ="o", lty=0, col = "blue", lwd=1)
        lines(years, C_R, type ="o", lty=0, col = "orange", lwd=1)    
      }
    }
  }
}
start_time <- Sys.time()
name <- 'Parameters_wide_range_uniform.csv'
loadnprint2(name)
end_time <- Sys.time()
end_time - start_time
