#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
// NumericMatrix
NumericMatrix SCCII(double Community1_Susceptible_UT, double Community1_Susceptible_T_CEPH, double Community1_Susceptible_T_CAR,
                    double Community1_Colonized_S_UT, double Community1_Colonized_S_T_CEPH, double Community1_Colonized_S_T_CAR,
                    double Community1_Colonized_ESBL_UT, double Community1_Colonized_ESBL_T_CEPH, double Community1_Colonized_ESBL_T_CAR,
                    double Community1_Colonized_CRE_UT, double Community1_Colonized_CRE_T_CEPH, double Community1_Colonized_CRE_T_CAR,
                    double Hospital1_Susceptible_UT, double Hospital1_Susceptible_T_CEPH, double Hospital1_Susceptible_T_CAR,
                    double Hospital1_Colonized_S_UT, double Hospital1_Colonized_S_T_CEPH, double Hospital1_Colonized_S_T_CAR,
                    double Hospital1_Colonized_ESBL_UT, double Hospital1_Colonized_ESBL_T_CEPH, double Hospital1_Colonized_ESBL_T_CAR,
                    double Hospital1_Colonized_CRE_UT, double Hospital1_Colonized_CRE_T_CEPH, double Hospital1_Colonized_CRE_T_CAR,
                    double Hospital1_Infected_S_T, double Hospital1_Infected_ESBL_T, double Hospital1_Infected_CRE_T,
                    double beta_d,
                    double mu,
                    double hospitalization_rate,
                    double discharge_rate,
                    double time_disease_development_community,
                    double time_disease_development_hospital,
                    double time_recovery,
                    double fitness_cost_S,
                    double fitness_cost_ESBL,
                    double fitness_cost_CRE,
                    double R_H_C,
                    double treatment_time_clearance, 
                    double treatment_duration_community_CEPH,
                    double treatment_duration_community_CAR,
                    double treatment_duration_hospital_CEPH,
                    double treatment_duration_hospital_CAR,
                    double conversion_rate,
                    double influence_infected,
                    double t_con,
                    double P_CEPH,
                    double P_CAR,
                    double clearance_rate,
                    double Import_ESBL,
                    double Import_CRE,
                    double HGT,
                    double HGT_D,
                    double ITC,
                    NumericVector consumption_community_cephalosporins,
                    NumericVector consumption_hospital_cephalosporins,
                    NumericVector consumption_community_carbapenems,
                    NumericVector consumption_hospital_carbapenems) {
  double Community1_Total = ( 
    Community1_Susceptible_UT    + Community1_Susceptible_T_CEPH + Community1_Susceptible_T_CAR
  
  + Community1_Colonized_S_UT    + Community1_Colonized_S_T_CEPH    + Community1_Colonized_S_T_CAR
  + Community1_Colonized_ESBL_UT + Community1_Colonized_ESBL_T_CEPH + Community1_Colonized_ESBL_T_CAR
  + Community1_Colonized_CRE_UT  + Community1_Colonized_CRE_T_CEPH  + Community1_Colonized_CRE_T_CAR
  );
  
  
  double  Hospital1_Total = (
    Hospital1_Susceptible_UT    + Hospital1_Susceptible_T_CEPH + Hospital1_Susceptible_T_CAR
    
    + Hospital1_Colonized_S_UT    + Hospital1_Colonized_S_T_CEPH    + Hospital1_Colonized_S_T_CAR
    + Hospital1_Colonized_ESBL_UT + Hospital1_Colonized_ESBL_T_CEPH + Hospital1_Colonized_ESBL_T_CAR
    + Hospital1_Colonized_CRE_UT  + Hospital1_Colonized_CRE_T_CEPH  + Hospital1_Colonized_CRE_T_CAR
    
    + Hospital1_Infected_S_T
    + Hospital1_Infected_ESBL_T
    + Hospital1_Infected_CRE_T
  );
  
  
  
  
  double  resistance_S_CEPH = 0;
  double  resistance_S_CAR  = 0;
  
  double  resistance_ESBL_CEPH = 1;
  double  resistance_ESBL_CAR  = 0;
  
  double  resistance_CRE_CEPH = 1;
  double  resistance_CRE_CAR  = 1;
  
  double  Hor_com  = 1;
  double  Dec_com  = 1;
  
  double  Hor_hosp = 1;
  double  Dec_hosp = 1;
  int mult = 50;
  NumericMatrix variables(27,3650*mult+1);
  variables(0,0) = Community1_Susceptible_UT;
  variables(1,0) = Community1_Susceptible_T_CEPH;
  variables(2,0) = Community1_Susceptible_T_CAR;
  variables(3,0) = Community1_Colonized_S_UT;
  variables(4,0) = Community1_Colonized_S_T_CEPH;
  variables(5,0) = Community1_Colonized_S_T_CAR;
  variables(6,0) = Community1_Colonized_ESBL_UT;
  variables(7,0) = Community1_Colonized_ESBL_T_CEPH;
  variables(8,0) = Community1_Colonized_ESBL_T_CAR;
  variables(9,0) = Community1_Colonized_CRE_UT;
  variables(10,0) = Community1_Colonized_CRE_T_CEPH;
  variables(11,0) = Community1_Colonized_CRE_T_CAR;
  
  variables(12,0) = Hospital1_Susceptible_UT;
  variables(13,0) = Hospital1_Susceptible_T_CEPH;
  variables(14,0) = Hospital1_Susceptible_T_CAR;
  variables(15,0) = Hospital1_Colonized_S_UT;
  variables(16,0) = Hospital1_Colonized_S_T_CEPH;
  variables(17,0) = Hospital1_Colonized_S_T_CAR;
  variables(18,0) = Hospital1_Colonized_ESBL_UT;
  variables(19,0) = Hospital1_Colonized_ESBL_T_CEPH;
  variables(20,0) = Hospital1_Colonized_ESBL_T_CAR;
  variables(21,0) = Hospital1_Colonized_CRE_UT;
  variables(22,0) = Hospital1_Colonized_CRE_T_CEPH;
  variables(23,0) = Hospital1_Colonized_CRE_T_CAR;
  
  variables(24,0) = Hospital1_Infected_S_T;
  variables(25,0) = Hospital1_Infected_ESBL_T;
  variables(26,0) = Hospital1_Infected_CRE_T;
  double t_step = 1.0/mult;
  double  beta_community1_S;
  double  beta_community1_ESBL;
  double  beta_community1_CRE;
  double  beta_hospital1_S;
  double  beta_hospital1_ESBL;
  double  beta_hospital1_CRE;
  double  treatment_rate_community1_CEPH;
  double  treatment_rate_community1_CAR;
  double  treatment_rate_hospital1_CEPH;
  double  treatment_rate_hospital1_CAR;
  
  NumericVector years(3650*mult+1);
  years(0) = 0;
  for (int i = 1; i <= 3650*mult; i += 1){
    beta_community1_S    = beta_d * (1 - fitness_cost_S   ) * ((variables(3,i-1) + variables(4,i-1) + variables(5,i-1))) / (Community1_Total);
    beta_community1_ESBL = beta_d * (1 - fitness_cost_ESBL) * ((variables(6,i-1) + variables(7,i-1) + variables(8,i-1) + 0.01*Import_ESBL)) / (Community1_Total); 
    beta_community1_CRE  = beta_d * (1 - fitness_cost_CRE ) * ((variables(9,i-1) + variables(10,i-1) + variables(11,i-1) + 0.01*Import_CRE)) / (Community1_Total);
    
    beta_hospital1_S     = R_H_C * beta_d * (1 - fitness_cost_S   ) * ((variables(15,i-1) + variables(16,i-1) + variables(17,i-1) + influence_infected * variables(24,i-1))) / (Hospital1_Total);
    beta_hospital1_ESBL  = R_H_C * beta_d * (1 - fitness_cost_ESBL) * ((variables(18,i-1) + variables(19,i-1) + variables(20,i-1) + influence_infected * variables(25,i-1))) / (Hospital1_Total); 
    beta_hospital1_CRE   = R_H_C * beta_d * (1 - fitness_cost_CRE ) * ((variables(21,i-1) + variables(22,i-1) + variables(23,i-1) + influence_infected * variables(26,i-1))) / (Hospital1_Total);
    
    
    int t = i/ mult;
    int curr_step = 0 + ((t) / 365);
    years(i) = curr_step;
    double C_U;
    double C_T_1;
    double C_T_2;
    double H_U;
    double H_T_1;
    double H_T_2;
    double discharge_rate_1;
    C_U = Community1_Total - (consumption_community_cephalosporins[curr_step] * (Community1_Total + Hospital1_Total)) / 1000 - (consumption_community_carbapenems[curr_step] * (Community1_Total + Hospital1_Total)) / 1000;
    C_T_1 = (consumption_community_cephalosporins[curr_step] * (Community1_Total + Hospital1_Total)) / 1000;
    C_T_2 = (consumption_community_carbapenems[curr_step] * (Community1_Total + Hospital1_Total)) / 1000;
    H_U = Hospital1_Total - (consumption_hospital_cephalosporins[curr_step] * (Community1_Total + Hospital1_Total)) / 1000 - (consumption_hospital_carbapenems[curr_step] * (Community1_Total + Hospital1_Total)) / 1000;
    H_T_1 = (consumption_hospital_cephalosporins[curr_step] * (Community1_Total + Hospital1_Total)) / 1000;
    H_T_2 = (consumption_hospital_carbapenems[curr_step] * (Community1_Total + Hospital1_Total)) / 1000;
    discharge_rate_1 = discharge_rate;
    discharge_rate = hospitalization_rate * (C_U + C_T_1 + C_T_2) / H_U;
    
    treatment_rate_community1_CEPH = (1 / (treatment_duration_community_CEPH * C_U)) * (C_T_1 * ( hospitalization_rate * treatment_duration_community_CEPH + 1));
    treatment_rate_community1_CAR  = (1 / (treatment_duration_community_CAR  * C_U)) * (C_T_2 * ( hospitalization_rate * treatment_duration_community_CAR  + 1));
    treatment_rate_hospital1_CEPH  = (1 / (treatment_duration_hospital_CEPH  * H_U)) * (H_T_1 - hospitalization_rate * treatment_duration_hospital_CEPH  * C_T_1);
    treatment_rate_hospital1_CAR   = (1 / (treatment_duration_hospital_CAR   * H_U)) * (H_T_2 - hospitalization_rate * treatment_duration_hospital_CAR   * C_T_2);
    
    //treatment_rate_community1_CEPH = 1 / treatment_duration_community_CEPH * consumption_community_cephalosporins[curr_step] * (Community1_Total + Hospital1_Total) / (1000 * Community1_Total-(Community1_Total + Hospital1_Total)*(consumption_community_cephalosporins[curr_step]+consumption_community_carbapenems[curr_step]));
    //treatment_rate_community1_CAR  = 1 / treatment_duration_community_CAR  * consumption_community_carbapenems[curr_step] * (Community1_Total + Hospital1_Total) / (1000 * Community1_Total-(Community1_Total + Hospital1_Total)*(consumption_community_cephalosporins[curr_step]+consumption_community_carbapenems[curr_step]));
    
    //treatment_rate_hospital1_CEPH = 1 / treatment_duration_hospital_CEPH * consumption_hospital_cephalosporins[curr_step] * (Community1_Total + Hospital1_Total) / (1000 * Hospital1_Total-(Community1_Total + Hospital1_Total)*(consumption_hospital_cephalosporins[curr_step]+consumption_hospital_carbapenems[curr_step]));
    //treatment_rate_hospital1_CAR  = 1 / treatment_duration_hospital_CAR  * consumption_hospital_carbapenems[curr_step]  * (Community1_Total + Hospital1_Total) / (1000 * Hospital1_Total-(Community1_Total + Hospital1_Total)*(consumption_hospital_cephalosporins[curr_step]+consumption_hospital_carbapenems[curr_step]));
    
    variables(0,i) = variables(0,i-1) + t_step * (
      - variables(0,i-1) * hospitalization_rate + variables(12,i-1) * discharge_rate
      - variables(0,i-1) * treatment_rate_community1_CEPH + variables(1,i-1) / treatment_duration_community_CEPH
      - variables(0,i-1) * treatment_rate_community1_CAR + variables(2,i-1) / treatment_duration_community_CAR
      
      - beta_community1_S    * variables(0,i-1)
      - beta_community1_ESBL * variables(0,i-1)
      - beta_community1_CRE  * variables(0,i-1)
      
      + variables(3,i-1)    * clearance_rate
      + variables(6,i-1) * clearance_rate
      + variables(9,i-1)  * clearance_rate
      
      + variables(24,i-1)    / time_recovery
      + variables(25,i-1) / time_recovery
      + variables(26,i-1)  / time_recovery
    );
    
    variables(1,i) = variables(1,i-1) + t_step * ( 
      - variables(1,i-1) * hospitalization_rate
      //+ variables(13,i-1) * discharge_rate
      + variables(0,i-1) * treatment_rate_community1_CEPH - variables(1,i-1) / treatment_duration_community_CEPH
      
      - ITC * beta_community1_ESBL * variables(1,i-1)
      - ITC * beta_community1_CRE  * variables(1,i-1)
      
      + variables(4,i-1)    * P_CEPH / treatment_time_clearance  * (1 - resistance_S_CEPH)
      + variables(7,i-1) * P_CEPH / treatment_time_clearance  * (1 - resistance_ESBL_CEPH)
      + variables(10,i-1)  * P_CEPH / treatment_time_clearance  * (1 - resistance_CRE_CEPH)
      
      + variables(4,i-1)    * clearance_rate
      + variables(7,i-1) * clearance_rate
      + variables(10,i-1)  * clearance_rate
      
    );
    
    variables(2,i) = variables(2,i-1) + t_step * (
      - variables(2,i-1) * hospitalization_rate 
      //+ variables(14,i-1) * discharge_rate
      + variables(0,i-1) * treatment_rate_community1_CAR - variables(2,i-1) / treatment_duration_community_CAR
      
      
      - ITC * beta_community1_CRE   * variables(2,i-1)
      
      + variables(5,i-1)    * P_CAR / treatment_time_clearance  * (1 - resistance_S_CAR)
      + variables(8,i-1) * P_CAR / treatment_time_clearance  * (1 - resistance_ESBL_CAR)
      + variables(11,i-1)  * P_CAR / treatment_time_clearance  * (1 - resistance_CRE_CAR)
      
      + variables(5,i-1)    * clearance_rate
      + variables(8,i-1) * clearance_rate
      + variables(11,i-1)  * clearance_rate
    );
    
    
    
    variables(3,i) = variables(3,i-1) + t_step * (
      - variables(3,i-1) * hospitalization_rate + variables(15,i-1) * discharge_rate
      - variables(3,i-1) * treatment_rate_community1_CEPH
      - variables(3,i-1) * treatment_rate_community1_CAR 
      
      + beta_community1_S    * variables(0,i-1)
      
      - Hor_com * HGT * beta_community1_ESBL * variables(3,i-1)
      - Hor_com * HGT * beta_community1_CRE  * variables(3,i-1)
      + Dec_com * HGT_D * clearance_rate * variables(6,i-1)
      + Dec_com * HGT_D * fitness_cost_CRE / fitness_cost_ESBL * clearance_rate * variables(9,i-1)
      
      + variables(4,i-1) / treatment_duration_community_CEPH
      + variables(5,i-1)  / treatment_duration_community_CAR
      
      - variables(3,i-1)         * clearance_rate
      - variables(3,i-1) / time_disease_development_community
    );
    
    variables(4,i) = variables(4,i-1) + t_step * (
      
      - variables(4,i-1) * hospitalization_rate 
      //+ variables(16,i-1) * discharge_rate
      + variables(3,i-1) * treatment_rate_community1_CEPH
      
      
      - mu * beta_community1_ESBL * variables(4,i-1)
      - mu * beta_community1_CRE  * variables(4,i-1)
      
      - variables(4,i-1) * resistance_ESBL_CEPH * conversion_rate
      - variables(4,i-1) * resistance_CRE_CEPH * conversion_rate
      
      - Hor_com * HGT * beta_community1_ESBL * variables(4,i-1) * resistance_ESBL_CEPH
      - Hor_com * HGT * beta_community1_CRE  * variables(4,i-1) * resistance_CRE_CEPH
      + Dec_com * HGT_D * clearance_rate * variables(7,i-1) * resistance_S_CEPH
      + Dec_com * HGT_D * fitness_cost_CRE / fitness_cost_ESBL * clearance_rate * variables(10,i-1) * resistance_S_CEPH
      
      - variables(4,i-1)  * P_CEPH / treatment_time_clearance * (1 - resistance_S_CEPH)
      - variables(4,i-1) / treatment_duration_community_CEPH
      
      - variables(4,i-1)     * clearance_rate
      - variables(4,i-1) / time_disease_development_community
    );
    
    variables(5,i) = variables(5,i-1) + t_step * (
      - variables(5,i-1) * hospitalization_rate 
      //+ variables(17,i-1) * discharge_rate
      + variables(3,i-1) * treatment_rate_community1_CAR 
      
      - mu * beta_community1_CRE  * variables(5,i-1)
      
      - variables(5,i-1) * resistance_ESBL_CAR * conversion_rate
      - variables(5,i-1) * resistance_CRE_CAR * conversion_rate
      
      - Hor_com * HGT * beta_community1_ESBL * variables(5,i-1) * resistance_ESBL_CAR
      - Hor_com * HGT * beta_community1_CRE  * variables(5,i-1) * resistance_CRE_CAR
      + Dec_com * HGT_D * clearance_rate * variables(8,i-1) * resistance_S_CAR
      + Dec_com * HGT_D * fitness_cost_CRE / fitness_cost_ESBL * clearance_rate * variables(11,i-1) * resistance_S_CAR
      
      - variables(5,i-1)   * P_CAR / treatment_time_clearance * (1 - resistance_S_CAR)
      - variables(5,i-1)  / treatment_duration_community_CAR
      
      - variables(5,i-1) * clearance_rate
      - variables(5,i-1) / time_disease_development_community
    );
    
    
    variables(6,i) = variables(6,i-1) + t_step * (
      - variables(6,i-1) * hospitalization_rate + variables(18,i-1) * discharge_rate
      - variables(6,i-1) * treatment_rate_community1_CEPH
      - variables(6,i-1) * treatment_rate_community1_CAR 
      
      + beta_community1_ESBL * variables(0,i-1)
      
      + Hor_com * HGT * beta_community1_ESBL * variables(3,i-1)
      - (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_com * HGT * beta_community1_CRE  * variables(6,i-1)
      - Dec_com * HGT_D * clearance_rate * variables(6,i-1)
      + Dec_com * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(9,i-1)
      
      + variables(7,i-1) / treatment_duration_community_CEPH
      + variables(8,i-1)  / treatment_duration_community_CAR
      - variables(6,i-1)      * clearance_rate
      - variables(6,i-1) / time_disease_development_community
    );
    
    variables(7,i) = variables(7,i-1) + t_step * (
      - variables(7,i-1) * hospitalization_rate 
      //+ variables(19,i-1) * discharge_rate
      + variables(6,i-1) * treatment_rate_community1_CEPH
      
      + ITC * beta_community1_ESBL * variables(1,i-1)
      
      + mu * beta_community1_ESBL * variables(4,i-1)
      
      + variables(4,i-1) * resistance_ESBL_CEPH * conversion_rate
      
      - (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_com * HGT * beta_community1_CRE  * variables(7,i-1) * resistance_CRE_CEPH
      + Hor_com * HGT * beta_community1_ESBL * variables(4,i-1) * resistance_ESBL_CEPH
      - Dec_com * HGT_D * clearance_rate * variables(7,i-1) * resistance_S_CEPH
      + Dec_com * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(10,i-1) * resistance_ESBL_CEPH
      
      - variables(7,i-1)  * P_CEPH / treatment_time_clearance * (1 - resistance_ESBL_CEPH)
      - variables(7,i-1) / treatment_duration_community_CEPH
      
      - variables(7,i-1)      * clearance_rate
      - variables(7,i-1) / time_disease_development_community
    );
    
    variables(8,i) = variables(8,i-1)+ t_step * ( 
      - variables(8,i-1) * hospitalization_rate 
      //+ variables(20,i-1) * discharge_rate
      + variables(6,i-1) * treatment_rate_community1_CAR 
      
      - mu *  beta_community1_CRE  * variables(8,i-1)
      
      + variables(5,i-1) * resistance_ESBL_CAR * conversion_rate
      - variables(8,i-1) * resistance_CRE_CEPH * conversion_rate
      
      - (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_com * HGT * beta_community1_CRE  * variables(8,i-1) * resistance_CRE_CAR
      + Hor_com * HGT * beta_community1_ESBL * variables(5,i-1) * resistance_ESBL_CAR
      - Dec_com * HGT_D * clearance_rate * variables(8,i-1) * resistance_S_CAR
      + Dec_com * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(11,i-1) * resistance_ESBL_CAR
      
      - variables(8,i-1)  * P_CAR / treatment_time_clearance * (1 - resistance_ESBL_CAR)
      - variables(8,i-1)  / treatment_duration_community_CAR
      
      - variables(8,i-1)      * clearance_rate
      - variables(8,i-1) / time_disease_development_community
    );
    
    
    
    variables(9,i) = variables(9,i-1) + t_step * ( 
      - variables(9,i-1) * hospitalization_rate + variables(21,i-1) * discharge_rate
      - variables(9,i-1) * treatment_rate_community1_CEPH
      - variables(9,i-1) * treatment_rate_community1_CAR 
      
      + beta_community1_CRE    * variables(0,i-1)
      
      + Hor_com * HGT * beta_community1_CRE  * variables(3,i-1)
      + (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) *  Hor_com * HGT * beta_community1_CRE  * variables(6,i-1)
      - Dec_com * HGT_D * fitness_cost_CRE / fitness_cost_ESBL * clearance_rate * variables(9,i-1)
      - Dec_com * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(9,i-1)
      
      + variables(10,i-1) / treatment_duration_community_CEPH
      + variables(11,i-1)  / treatment_duration_community_CAR
      
      - variables(9,i-1)      * clearance_rate
      - variables(9,i-1) / time_disease_development_community
    );
    
    variables(10,i) = variables(10,i-1) + t_step * ( 
      - variables(10,i-1) * hospitalization_rate 
      //+ variables(22,i-1) * discharge_rate
      + variables(9,i-1) * treatment_rate_community1_CEPH
      
      + ITC * beta_community1_CRE    * variables(1,i-1)   
      
      + mu * beta_community1_CRE  * variables(4,i-1)
      
      + variables(4,i-1) * resistance_CRE_CEPH * conversion_rate
      
      + Hor_com * HGT * beta_community1_CRE  * variables(4,i-1) * resistance_CRE_CEPH
      + (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_com * HGT * beta_community1_CRE  * variables(7,i-1) * resistance_CRE_CEPH
      - Dec_com * HGT_D * fitness_cost_CRE / fitness_cost_ESBL * clearance_rate * variables(10,i-1) * resistance_S_CEPH
      - Dec_com * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(10,i-1) * resistance_ESBL_CEPH
      
      - variables(10,i-1) / treatment_duration_community_CEPH
      - variables(10,i-1)  * P_CEPH / treatment_time_clearance     * (1 - resistance_CRE_CEPH)
      
      - variables(10,i-1)      * clearance_rate
      - variables(10,i-1)  / time_disease_development_community
    );
    
    variables(11,i) = variables(11,i-1) + t_step * (
      - variables(11,i-1) * hospitalization_rate 
      //+ variables(23,i-1) * discharge_rate
      + variables(9,i-1) * treatment_rate_community1_CAR 
      
      + ITC * beta_community1_CRE    * variables(2,i-1)
      
      + mu * beta_community1_CRE  * variables(5,i-1)
      + mu * beta_community1_CRE  * variables(8,i-1)
      
      + variables(5,i-1) * resistance_CRE_CAR * conversion_rate
      + variables(8,i-1) * resistance_CRE_CAR * conversion_rate
      
      + Hor_com * HGT * beta_community1_CRE  * variables(5,i-1) * resistance_CRE_CAR
      + (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_com * HGT * beta_community1_CRE  * variables(8,i-1) * resistance_CRE_CAR
      - Dec_com * HGT_D * fitness_cost_CRE / fitness_cost_ESBL * clearance_rate * variables(11,i-1) * resistance_S_CAR
      - Dec_com * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(11,i-1) * resistance_ESBL_CAR
      
      - variables(11,i-1)  / treatment_duration_community_CAR
      - variables(11,i-1)  * P_CAR / treatment_time_clearance      * (1 - resistance_CRE_CAR)
      
      - variables(11,i-1)      * clearance_rate
      - variables(11,i-1)  / time_disease_development_community
    );
    
    variables(12,i) = variables(12,i-1) + t_step * (
      + variables(0,i-1) * hospitalization_rate - variables(12,i-1) * discharge_rate
      - variables(12,i-1) * treatment_rate_hospital1_CEPH + variables(13,i-1) / treatment_duration_hospital_CEPH
      - variables(12,i-1) * treatment_rate_hospital1_CAR + variables(14,i-1) / treatment_duration_hospital_CAR
      
      - beta_hospital1_S    * variables(12,i-1)
      - beta_hospital1_ESBL * variables(12,i-1)
      - beta_hospital1_CRE  * variables(12,i-1)
      
      + variables(15,i-1)    * clearance_rate
      + variables(18,i-1) * clearance_rate
      + variables(21,i-1)  * clearance_rate
    );
    
    variables(13,i) = variables(13,i-1) + t_step * (
      + variables(1,i-1) * hospitalization_rate 
      //- variables(13,i-1) * discharge_rate
      + variables(12,i-1) * treatment_rate_hospital1_CEPH - variables(13,i-1) / treatment_duration_hospital_CEPH
      
      - ITC * beta_hospital1_ESBL * variables(13,i-1)
      - ITC * beta_hospital1_CRE  * variables(13,i-1)
      
      + variables(16,i-1)  * P_CEPH / treatment_time_clearance * (1 - resistance_S_CEPH)
      + variables(19,i-1)  * P_CEPH / treatment_time_clearance * (1 - resistance_ESBL_CEPH)
      + variables(22,i-1)  * P_CEPH / treatment_time_clearance * (1 - resistance_CRE_CEPH)
      
      + variables(16,i-1)    * clearance_rate
      + variables(19,i-1) * clearance_rate
      + variables(22,i-1)  * clearance_rate
    );
    
    variables(14,i) = variables(14,i-1) + t_step * (
      + variables(2,i-1) * hospitalization_rate 
      //- variables(14,i-1) * discharge_rate
      +  variables(12,i-1) * treatment_rate_hospital1_CAR - variables(14,i-1) / treatment_duration_hospital_CAR
      
      
      
      - ITC * beta_hospital1_CRE   * variables(14,i-1)
      
      + variables(17,i-1) * P_CAR / treatment_time_clearance    * (1 - resistance_S_CAR)
      + variables(20,i-1) * P_CAR  / treatment_time_clearance   * (1 - resistance_ESBL_CAR)
      + variables(23,i-1) * P_CAR / treatment_time_clearance    * (1 - resistance_CRE_CAR)
      
      + variables(17,i-1) * clearance_rate
      + variables(20,i-1) * clearance_rate
      + variables(23,i-1) * clearance_rate
    );
    
    variables(15,i) = variables(15,i-1) + t_step * ( 
      + variables(3,i-1)* hospitalization_rate - variables(15,i-1) * discharge_rate
      - variables(15,i-1) * treatment_rate_hospital1_CEPH
      - variables(15,i-1) * treatment_rate_hospital1_CAR 
      
      + beta_hospital1_S       * variables(12,i-1)
      
      - Hor_hosp * HGT * beta_hospital1_ESBL * variables(15,i-1)
      - Hor_hosp * HGT * beta_hospital1_CRE  * variables(15,i-1) 
      + Dec_hosp * HGT_D * clearance_rate * variables(18,i-1)
      + Dec_hosp * HGT_D * (fitness_cost_CRE) / fitness_cost_ESBL * clearance_rate * variables(21,i-1)
      
      + variables(16,i-1) / treatment_duration_hospital_CEPH
      + variables(17,i-1)  / treatment_duration_hospital_CAR 
      
      - variables(15,i-1)         * clearance_rate
      - variables(15,i-1) / time_disease_development_hospital
    );
    
    variables(16,i) = variables(16,i-1) + t_step * (
      + variables(4,i-1)* hospitalization_rate 
      //- variables(16,i-1) * discharge_rate
      + variables(15,i-1) * treatment_rate_hospital1_CEPH
      
      
      - mu * beta_hospital1_ESBL    * variables(16,i-1)
      - mu * beta_hospital1_CRE     * variables(16,i-1)
      
      - variables(16,i-1) * resistance_ESBL_CEPH * conversion_rate
      - variables(16,i-1) * resistance_CRE_CEPH * conversion_rate
      
      - Hor_hosp * HGT * beta_hospital1_ESBL * variables(16,i-1) * resistance_ESBL_CEPH
      - Hor_hosp * HGT * beta_hospital1_CRE  * variables(16,i-1) * resistance_CRE_CEPH
      + Dec_hosp * HGT_D * clearance_rate * variables(19,i-1) * resistance_S_CEPH
      + Dec_hosp * HGT_D * (fitness_cost_CRE) / fitness_cost_ESBL * clearance_rate * variables(22,i-1) * resistance_S_CEPH
      
      - variables(16,i-1) / treatment_duration_hospital_CEPH
      - variables(16,i-1) * P_CEPH / treatment_time_clearance     * (1 - resistance_S_CEPH)
      
      - variables(16,i-1)     * clearance_rate
      - variables(16,i-1) / time_disease_development_hospital
    );
    
    variables(17,i) = variables(17,i-1) + t_step * ( 
      + variables(5,i-1)* hospitalization_rate 
      //- variables(17,i-1) * discharge_rate
      + variables(15,i-1) * treatment_rate_hospital1_CAR 
      
      - mu * beta_hospital1_CRE     * variables(17,i-1)
      
      - variables(17,i-1) * resistance_ESBL_CAR * conversion_rate
      - variables(17,i-1) * resistance_CRE_CAR * conversion_rate
      
      - Hor_hosp * HGT * beta_hospital1_ESBL * variables(17,i-1) * resistance_ESBL_CAR
      - Hor_hosp * HGT * beta_hospital1_CRE  * variables(17,i-1) * resistance_CRE_CAR
      + Dec_hosp * HGT_D * clearance_rate * variables(20,i-1) * resistance_S_CAR
      + Dec_hosp * HGT_D * (fitness_cost_CRE) / fitness_cost_ESBL * clearance_rate * variables(23,i-1) * resistance_S_CAR
      
      - variables(17,i-1)  / treatment_duration_hospital_CAR
      - variables(17,i-1) * P_CAR  / treatment_time_clearance      * (1 - resistance_S_CAR)
      
      - variables(17,i-1)      * clearance_rate
      - variables(17,i-1)/ time_disease_development_hospital
    );
    
    variables(18,i) = variables(18,i-1) + t_step * ( 
      + variables(6,i-1) * hospitalization_rate - variables(18,i-1) * discharge_rate
      - variables(18,i-1) * treatment_rate_hospital1_CEPH
      - variables(18,i-1) * treatment_rate_hospital1_CAR 
      
      + beta_hospital1_ESBL    * variables(12,i-1)
      
      + Hor_hosp * HGT * beta_hospital1_ESBL  * variables(15,i-1)
      - (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_hosp * HGT * beta_hospital1_CRE   * variables(18,i-1) 
      + Dec_hosp * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(21,i-1)
      - Dec_hosp * HGT_D * clearance_rate * variables(18,i-1)
      
      + variables(19,i-1) / treatment_duration_hospital_CEPH
      + variables(20,i-1)  / treatment_duration_hospital_CAR
      
      - variables(18,i-1)      * clearance_rate
      - variables(18,i-1) / time_disease_development_hospital
    );
    
    variables(19,i) = variables(19,i-1) + t_step * ( 
      + variables(7,i-1) * hospitalization_rate 
      //- variables(19,i-1) * discharge_rate
      + variables(18,i-1) * treatment_rate_hospital1_CEPH
      
      + ITC * beta_hospital1_ESBL    * variables(13,i-1)
      
      + mu * beta_hospital1_ESBL * variables(16,i-1)
      
      + variables(16,i-1) * resistance_ESBL_CEPH * conversion_rate
      
      + Hor_hosp * HGT * beta_hospital1_ESBL * variables(16,i-1) * resistance_ESBL_CEPH
      - (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_hosp * HGT * beta_hospital1_CRE  * variables(19,i-1) * resistance_CRE_CEPH
      + Dec_hosp * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(22,i-1) * resistance_ESBL_CEPH
      - Dec_hosp * HGT_D * clearance_rate * variables(19,i-1) * resistance_S_CEPH
      
      - variables(19,i-1) * P_CEPH / treatment_duration_hospital_CEPH
      - variables(19,i-1) / treatment_time_clearance     * (1 - resistance_ESBL_CEPH)
      
      - variables(19,i-1)      * clearance_rate
      - variables(19,i-1) / time_disease_development_hospital
    );
    
    variables(20,i) = variables(20,i-1) + t_step * ( 
      + variables(8,i-1) * hospitalization_rate 
      //- variables(20,i-1) * discharge_rate
      + variables(18,i-1) * treatment_rate_hospital1_CAR 
      
      
      - mu * beta_hospital1_CRE * variables(20,i-1)
      
      + variables(17,i-1) * resistance_ESBL_CAR * conversion_rate
      - variables(20,i-1) * resistance_CRE_CAR * conversion_rate
      
      + Hor_hosp * HGT * beta_hospital1_ESBL  * variables(17,i-1) * resistance_ESBL_CAR
      - (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_hosp * HGT * beta_hospital1_CRE   * variables(20,i-1) * resistance_CRE_CAR
      + Dec_hosp * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(23,i-1) * resistance_ESBL_CAR
      - Dec_hosp * HGT_D * clearance_rate * variables(20,i-1) * resistance_S_CAR
      
      - variables(20,i-1)  / treatment_duration_hospital_CAR
      - variables(20,i-1)  * P_CAR / treatment_time_clearance      * (1 - resistance_ESBL_CAR)
      
      - variables(20,i-1)      * clearance_rate
      - variables(20,i-1) / time_disease_development_hospital
    );
    
    
    
    variables(21,i) = variables(21,i-1) + t_step * ( 
      + variables(9,i-1) * hospitalization_rate - variables(21,i-1) * discharge_rate
      - variables(21,i-1) * treatment_rate_hospital1_CEPH
      - variables(21,i-1) * treatment_rate_hospital1_CAR 
      
      + beta_hospital1_CRE    * variables(12,i-1)
      
      + Hor_hosp * HGT * beta_hospital1_CRE  * variables(15,i-1)
      + (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_hosp * HGT * beta_hospital1_CRE  * variables(18,i-1)
      - Dec_hosp * HGT_D * (fitness_cost_CRE) / fitness_cost_ESBL * clearance_rate * variables(21,i-1)
      - Dec_hosp * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(21,i-1) 
      
      + variables(22,i-1) / treatment_duration_hospital_CEPH
      + variables(23,i-1)  / treatment_duration_hospital_CAR
      
      - variables(21,i-1)      * clearance_rate
      - variables(21,i-1) / time_disease_development_hospital
    );
    
    variables(22,i) = variables(22,i-1) + t_step * ( 
      + variables(10,i-1) * hospitalization_rate 
      //- variables(22,i-1) * discharge_rate
      + variables(21,i-1) * treatment_rate_hospital1_CEPH
      
      + ITC * beta_hospital1_CRE    * variables(13,i-1)
      
      + mu * beta_hospital1_CRE * variables(16,i-1)
      
      + variables(16,i-1) * resistance_CRE_CEPH * conversion_rate
      
      + Hor_hosp * HGT * beta_hospital1_CRE  * variables(16,i-1) * resistance_CRE_CEPH
      + (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_hosp * HGT * beta_hospital1_CRE  * variables(19,i-1) * resistance_CRE_CEPH
      - Dec_hosp * HGT_D * (fitness_cost_CRE) / fitness_cost_ESBL * clearance_rate * variables(22,i-1) * resistance_S_CEPH
      - Dec_hosp * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(22,i-1) * resistance_ESBL_CEPH
      
      - variables(22,i-1) / treatment_duration_hospital_CEPH
      - variables(22,i-1)  * P_CEPH / treatment_time_clearance     * (1 - resistance_CRE_CEPH)
      
      - variables(22,i-1)      * clearance_rate
      - variables(22,i-1) / time_disease_development_hospital
    );
    
    variables(23,i) = variables(23,i-1) + t_step * ( 
      + variables(11,i-1) * hospitalization_rate 
      //- variables(23,i-1) * discharge_rate
      + variables(21,i-1) * treatment_rate_hospital1_CAR 
      
      + ITC * beta_hospital1_CRE    * variables(14,i-1)
      
      + mu * beta_hospital1_CRE    * variables(17,i-1)
      + mu * beta_hospital1_CRE    * variables(20,i-1)
      
      + variables(17,i-1) * resistance_ESBL_CAR * conversion_rate
      + variables(20,i-1) * resistance_CRE_CAR * conversion_rate
      
      + Hor_hosp * HGT * beta_hospital1_CRE  * variables(17,i-1) * resistance_CRE_CAR
      + (1 - fitness_cost_CRE+fitness_cost_ESBL) / (1 - fitness_cost_CRE) * Hor_hosp * HGT * beta_hospital1_CRE  * variables(20,i-1) * resistance_CRE_CAR
      - Dec_hosp * HGT_D * (fitness_cost_CRE) / fitness_cost_ESBL * clearance_rate * variables(23,i-1) * resistance_S_CAR
      - Dec_hosp * HGT_D * (fitness_cost_CRE-fitness_cost_ESBL) / fitness_cost_ESBL * clearance_rate * variables(23,i-1) * resistance_ESBL_CAR
      
      - variables(23,i-1)  / treatment_duration_hospital_CAR
      - variables(23,i-1)  * P_CAR / treatment_time_clearance      * (1 - resistance_CRE_CAR)
      - variables(23,i-1)  * clearance_rate
      - variables(23,i-1) / time_disease_development_hospital
    );
    
    variables(24,i)    = variables(24,i-1) + t_step * (
      + variables(3,i-1) / time_disease_development_community
      + variables(4,i-1) / time_disease_development_community
      + variables(5,i-1) / time_disease_development_community
      + variables(15,i-1) / time_disease_development_hospital
      + variables(16,i-1) / time_disease_development_hospital
      + variables(17,i-1) / time_disease_development_hospital
      - variables(24,i-1) / time_recovery
    );
    
    variables(25,i) = variables(25,i-1) + t_step * (
      + variables(6,i-1) / time_disease_development_community
      + variables(7,i-1) / time_disease_development_community
      + variables(8,i-1)  / time_disease_development_community
      + variables(18,i-1) / time_disease_development_hospital
      + variables(19,i-1)  / time_disease_development_hospital
      + variables(20,i-1)  / time_disease_development_hospital
      - variables(25,i-1) / time_recovery
    );
    variables(26,i)  = variables(26,i-1) + t_step * (
      + variables(9,i-1) / time_disease_development_community
      + variables(10,i-1)  / time_disease_development_community
      + variables(11,i-1)  / time_disease_development_community
      + variables(21,i-1) / time_disease_development_hospital
      + variables(22,i-1)  / time_disease_development_hospital
      + variables(23,i-1)  / time_disease_development_hospital
      - variables(26,i-1) / time_recovery
    );
  }
  NumericMatrix out(3,11);
  for(int i = 0; i <= (10); i += 1) {
    out(0,i) = variables(24,i*365*mult);
    out(1,i) = variables(25,i*365*mult);
    out(2,i) = variables(26,i*365*mult);
    //out(0,i) = variables(24,i*365*mult);
    //out(1,i) = variables(25,i*365*mult);
    //out(2,i) = variables(26,i*365*mult);
    
    //   for(int j=0; j<=26; j+= 1){
    //      out(3,i) = out(3,i) + variables(j,i*365*mult);
    //   }
    //out(3,i) = consumption_community_cephalosporins[i];
  }
  return(out);
}
