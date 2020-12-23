library(readstata13)
library(tidyr)
library(MASS)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(purrr)
library(iterators)
library(doParallel)
## library(SparseM)

## Define Function ##

simulation_process = function(parameter_vector) {

year_start = years_with_shock[1]
shock_vector = shock_input_manual %>% filter(Year==year_start) %>% select(-Year) %>% slice(1) %>% as.numeric
shock_matrix = diag(shock_vector)
year_num = length(year_id[year_id > year_start])
single_result = matrix(0,nrow=1,ncol=country_num)
single_result_long_dataframe = as.data.frame(diag(rep(0,country_num)))
y_real = rep(0, country_num)
y_real_matrix = diag(y_real)
year_start_position = which(year_id == year_start)
years_with_effect = year_id[year_id >= year_start]
if (year_start == 1900){
   current_year = 18995
}else{
   current_year = year_start
}

for (j in (year_start_position+1):60) {
    current_year = current_year+20
    y_star = parameter_vector[3]*(shock_vector+y_real)
    y_star_matrix = parameter_vector[3]*(shock_matrix + y_real_matrix)
    ## start_index = (j-1)*country_num+1
    ## end_index = j*country_num
    W_N_current = W_N_shaped %>% filter(str_detect(cyr, as.character(current_year))) %>% select(-cyr) %>% as.matrix()
    #W_A_current = W_A_shaped %>% filter(str_detect(cyr, as.character(current_year))) %>% select(-cyr) %>% as.matrix()
    #LRSS_multiplier = solve(identity_m - parameter_vector[1]*W_N_current- parameter_vector[2]*W_A_current)
    LRSS_multiplier = solve(identity_m - parameter_vector[1]*W_N_current)
    ## LRSS_multiplier = solve(identity_m - draws[j,1]*W_N_SAm[start_index:end_index, start_index:end_index]- draws[j,2]*W_A_SAm[start_index:end_index, start_index:end_index])
    delta_multiplier = parameter_vector[1]*W_N_current
    ## delta_multiplier = draws[i,1]*W_N_SAm[start_index:end_index, start_index:end_index] + draws[i,2]*W_A_SAm[start_index:end_index, start_index:end_index]
    delta_temp = LRSS_multiplier %*% delta_multiplier %*% y_star
    delta_temp_matrix = LRSS_multiplier %*% delta_multiplier %*% y_star_matrix
    if (current_year %in% years_with_shock){
      shock_vector = (shock_input_manual%>%filter(Year==current_year)%>%select(-Year)%>%slice(1)%>%as.numeric())
    }else{
      shock_vector = rep(0, country_num)
    }

    shock_matrix = diag(shock_vector)

    y_real = delta_temp
    y_real_matrix = delta_temp_matrix
      ## single_result = rbind(single_result,(t(delta_temp)+shock_vector))
    single_result = rbind(single_result,(t(delta_temp)))
    single_result_long_dataframe = rbind(single_result_long_dataframe, delta_temp_matrix)
  }
  single_result = rbind(matrix(0, (year_start_position-1), country_num), single_result)
  single_result = cbind(year_id, single_result)
  colnames(single_result) = c('Year',country_id)

  single_result_long_dataframe = rbind(matrix(0, country_num*(year_start_position-1), country_num), single_result_long_dataframe)
  single_result_long_dataframe = cbind(rep(year_id,country_num) %>% sort, single_result_long_dataframe)
  colnames(single_result_long_dataframe) = c('Year',country_id)
  rownames(single_result_long_dataframe) =NULL

  ## return(list(single_result = single_result, single_result_long_dataframe = single_result_long_dataframe))
  return(single_result)
}

core_num = 24
c1<-makeCluster(core_num)
registerDoParallel(c1)



weight_matrix = readRDS("IDs_Weight_matrices.rds")
country_id = weight_matrix[[1]]
year_id = weight_matrix[[2]]
cyr = weight_matrix[[3]]
W_N_shaped = weight_matrix[[4]]
#W_A_shaped = weight_matrix[[5]]

## Preparation for simulation
vcvm = read.csv('Vcm_N.csv', stringsAsFactors = FALSE)[c(1,16),c(1,16)] %>%
  as.matrix()
rho_W_N = 0.0188566
phi = 0.9442287
num_draws <- 5000  #number of simulations to draw
coefficient_star = c(rho_W_N, phi)

library(dplyr)
shock_input_manual = read.csv('shock_input.csv', stringsAsFactors = FALSE) %>%
  filter_at(vars(!contains('Year')), any_vars(!is.na(.)))
shock_brief = shock_input_manual %>%
  pivot_longer(-Year, names_to="variable", values_to="value") %>%
  drop_na('value') %>%
  arrange(Year)
shock_input_manual = shock_input_manual %>%
                     mutate_all(~replace(., is.na(.), 0))


print('The briefing of input shock is: ')
print(shock_brief)

draws <- mvrnorm(n=num_draws,coefficient_star,vcvm,empirical=TRUE)
years_with_shock = shock_input_manual %>% pull(Year)


country_num = length(country_id)
identity_m <- diag(1,country_num,country_num)  #nxn identity matrix

#i is the draw number
#j is the country number, used to define the first and last row and
# column numbers in each submatrix below.

itx = iter(draws, by='row')

result = foreach(coefficient_vector = itx, .packages = c('tidyverse', 'doParallel', 'iterators')) %dopar% {

  simulation_process(coefficient_vector)

}

saveRDS(result, "simulated_effectsN.rds")
