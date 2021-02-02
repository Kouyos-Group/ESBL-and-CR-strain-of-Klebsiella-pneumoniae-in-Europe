# Identifying the drivers of multidrug-resistant Klebsiella pneumoniae at a European level
Code for Kachalov VN, Nguyen H, Balakrishna S, Salazar-Vizcaya L, Sommerstein R, Kuster SP, et al. (2021) Identifying the drivers of multidrug-resistant Klebsiella pneumoniae at a European level. PLoS Comput Biol 17(1): e1008446. https://doi.org/10.1371/journal.pcbi.1008446.

If you want to use the model, see the paper and supporting information for the full model description.
#Files description
## Main code
Main_model.R consists of the main code.

Main_model_plotting.R is the R code for plotting the main plot from the paper.

SCCII2.cpp is the C++ code for the simulation.

SCCH.cpp is the C++ code for determining the number of colonized individuals in the community and hospital with the given transmission rates.


Results_original_and_one_out_uniform.csv and Results_original_and_one_out_unique.csv are the results of the original fit and the one-out sensitivity analysis.
Presented values should be converted from the (-Inf, Inf) to the (lower_boundary, upper_boundary) (see the code for details).

## Data
{Country_name}_3rd_4th_gen.csv are the files with the data for the R code.

{Country_name}_3rd_4th_gen_detailed.csv are the files with the data presented with more details. The distribution of consumption by the generation of cephalosporins and the parameters which were used to fill the missing data are given.

Hospital beds and length of stay.xlsx consist the data about mean length of stay in hospital and number of beds per 100000, according to WHO.

Presented data on antibiotic consumption and prevalence of the resistance from The European Surveillance System – TESSy, provided by (Croatia, Denmark, Finland, France, Greece, Hungary, Italy, Netherlands, Norway, Portugal, Sweden) and released by ECDC. As required by the ECDC, we confirm that “the views and opinions of the authors expressed herein do not necessarily state or reflect those of the ECDC. The accuracy of the authors' statistical analysis and the findings they report are not the responsibility of ECDC. ECDC is not responsible for conclusions or opinions drawn from the data provided. ECDC is not responsible for the correctness of the data and for data management, data merging and data collation after provision of the data. ECDC shall not be held liable for improper or incorrect use of the data”

## Multivariate sensitivity analysis
Printing_multivar_sens_unique.R and Printing_multivar_sens_uniform.R could be used to print the results of multivariate sensitivity analysis.

##Citation
Kachalov VN, Nguyen H, Balakrishna S, Salazar-Vizcaya L, Sommerstein R, Kuster SP, et al. (2021) Identifying the drivers of multidrug-resistant Klebsiella pneumoniae at a European level. PLoS Comput Biol 17(1): e1008446. https://doi.org/10.1371/journal.pcbi.1008446.

##Contact
Viacheslav Kachalov: viacheslav.kachalov@usz.ch or kachalov93@gmail.com
