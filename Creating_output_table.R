setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
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
#the output matrix
output <- matrix(nrow = 42, ncol = 2)
#Recalculation of parameters
output[1,1] <- "Parameter"
output[1,2] <- "Value"
#Fitness cost ESBL
output[2,1] <- "Fitness cost ESBL"
output[2,2] <- signif((0.5+atan(fitness_cost1[1])/pi)*(0.15),4)*100
#Fitness cost CRK
output[3,1] <- "Fitness cost CRK"
output[3,2] <- signif((0.5+atan(fitness_cost1[1])/pi)*(0.15)+(0.5+atan(fitness_cost1[2])/pi)*(0.15),4)*100
#Import values
output[4,1] <- "Import of ESBL (reservoir size per 100000)"
output[5,1] <- "Import of ESBL (reservoir size per 100000)"
output[4,2] <- signif((0.5+atan(fitness_cost1[3])/pi)*1000,4)
output[5,2] <- signif((0.5+atan(fitness_cost1[4])/pi)*1000,4)
#helping calculation
alpha <-  (0.5+atan(fitness_cost1[6])/pi) * 0.5
tau = 365/(3*(0.5 + alpha))
beta_tr1 <- 1/((1-0.2)*tau)
#Colonization value
output[6,1] <- "Colonization rate"
output[6,2] <- signif(1/((1-0.2)*tau),4)
#Super-colonization coefficient
output[7,1] <- "Super-colonization coefficient (nu)"
output[7,2] <- (0.5+atan(fitness_cost1[40])/pi) * 4 * (1 - 0.5 - alpha)/(0.5 + alpha)
#Use this one if you want to calculate for the basic colonization, which differs from 20%
#base_c <- 0.2
#output[7,2] <- (0.5+atan(input$fitness_cost1[40])/pi) * (1-base_c)/base_c * (1 - 0.5 - alpha)/(0.5 + alpha)
#Increased susceptibility by treatment
output[8,1] <- "Increased susceptibility by treatment (mu)"
output[8,2] <- signif((0.5+atan(fitness_cost1[5])/pi),4)*100
output[9,1] <- "Displacement/loss of plasmid rate (relation to natural decolonization rate)"
output[9,2] <- (0.5-alpha)/(0.5+alpha)
#Reordering countries by classes
order_c <- c(4,7,8,5,1,6,2,9,10,3,11)
#Parameters, which depends on country
for (c in 1:11){
  name_curr <- paste(input_names[c], "3rd_4th_gen.csv", sep="_")
  beta_tr2 <- (0.5+atan(fitness_cost1[6+c])/pi)*50*beta_tr1
  output[9+order_c[c], 1] <-  paste("Hospital transmission rate in",input_names[c])
  output[9+order_c[c], 2] <-  signif((0.5+atan(fitness_cost1[6+c])/pi)*50,6)
  output[9+11+order_c[c], 1] <-  paste("Initial prevelence of ESBL strain in",input_names[c])
  output[9+11+order_c[c], 2] <-  signif((0.5+atan(fitness_cost1[6+11+c])/pi)*100,4)
  output[9+11+11+order_c[c], 1] <-  paste("Initial prevelence of CRK strain in",input_names[c])
  output[9+11+11+order_c[c], 2] <-  signif((0.5+atan(fitness_cost1[6+11+11+c])/pi)*(0.5+atan(fitness_cost1[6+11+c])/pi)*100,4)
}
