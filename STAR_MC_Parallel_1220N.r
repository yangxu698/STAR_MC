library(readstata13)
library(tidyr)
library(MASS)
library(tidyverse)
library(reshape2)
library(iterators)
library(doParallel)
## library(SparseM)

## Define Function ##
rm(list=ls())
simulation_process = function(parameter_vector) {

year_start = years_with_shock[1]
shock_vector = shock_input_manual %>% filter(Year==year_start) %>% dplyr::select(-Year) %>% slice(1) %>% as.numeric
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
    y_star = parameter_vector[2]*(shock_vector+y_real)
    y_star_matrix = parameter_vector[2]*(shock_matrix + y_real_matrix)
    ## start_index = (j-1)*country_num+1
    ## end_index = j*country_num
    W_N_current = W_N_shaped %>% filter(str_detect(cyr, as.character(current_year))) %>% dplyr::select(-cyr) %>% as.matrix()
    pop_weights_current = pop_weights %>% filter(str_detect(cyr, as.character(current_year))) %>%
      arrange(cyr) %>% dplyr::pull(GMpoplog) %>% rep(181) %>% matrix(nrow=181,ncol=181, byrow=TRUE)
    W_N_current = W_N_current * pop_weights_current
    #W_A_current = W_A_shaped %>% filter(str_detect(cyr, as.character(current_year))) %>% select(-cyr) %>% as.matrix()
    #LRSS_multiplier = solve(identity_m - parameter_vector[1]*W_N_current- parameter_vector[2]*W_A_current)
    LRSS_multiplier = solve(identity_m - parameter_vector[1]*W_N_current)
    ## LRSS_multiplier = solve(identity_m - draws[j,1]*W_N_SAm[start_index:end_index, start_index:end_index]- draws[j,2]*W_A_SAm[start_index:end_index, start_index:end_index])
    delta_multiplier = parameter_vector[1]*W_N_current
    ## delta_multiplier = draws[i,1]*W_N_SAm[start_index:end_index, start_index:end_index] + draws[i,2]*W_A_SAm[start_index:end_index, start_index:end_index]
    delta_temp = LRSS_multiplier %*% delta_multiplier %*% y_star
    delta_temp_matrix = LRSS_multiplier %*% delta_multiplier %*% y_star_matrix
    if (current_year %in% years_with_shock){
      shock_vector = (shock_input_manual%>%filter(Year==current_year)%>%dplyr::select(-Year)%>%slice(1)%>%as.numeric())
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

pop_weights = readRDS("Pop_weights.rds")
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
num_draws <- 5000 #number of simulations to draw
coefficient_star = c(rho_W_N, phi)

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

core_num = 24
c1<-makeCluster(core_num)
registerDoParallel(c1)


itx = iter(draws, by='row')

result = foreach(coefficient_vector = itx, .packages = c('tidyverse', 'doParallel', 'iterators')) %dopar% {

  simulation_process(coefficient_vector)

}

arraymedian <- apply(simplify2array(result), 1:2, median)
arraylower <- apply(simplify2array(result), 1:2, quantile, probs = 0.025)
arrayupper <- apply(simplify2array(result), 1:2, quantile, probs = 0.975)
arraylower90 <- apply(simplify2array(result), 1:2, quantile, probs = 0.05)
arrayupper90 <- apply(simplify2array(result), 1:2, quantile, probs = 0.95)

medianlong <- as.data.frame(arraymedian) %>%
  rename(c("year" = "Year")) %>%
  pivot_longer(-c(year),
               names_to = "country_text_id",
               values_to = "median") %>%
  arrange(country_text_id, year)

lowerlong <- as.data.frame(arraylower) %>%
  rename(c("year" = "Year")) %>%
  pivot_longer(-c(year),
               names_to = "country_text_id",
               values_to = "lower") %>%
  arrange(country_text_id, year)

upperlong <- as.data.frame(arrayupper) %>%
  rename(c("year" = "Year")) %>%
  pivot_longer(-c(year),
               names_to = "country_text_id",
               values_to = "upper") %>%
  arrange(country_text_id, year)

lowerlong90 <- as.data.frame(arraylower90) %>%
  rename(c("year" = "Year")) %>%
  pivot_longer(-c(year),
               names_to = "country_text_id",
               values_to = "lower90") %>%
  arrange(country_text_id, year)

upperlong90 <- as.data.frame(arrayupper90) %>%
  rename(c("year" = "Year")) %>%
  pivot_longer(-c(year),
               names_to = "country_text_id",
               values_to = "upper90") %>%
  arrange(country_text_id, year)

MedLowUplong <- cbind(medianlong, lowerlong$lower, upperlong$upper, lowerlong90$lower90, upperlong90$upper90) %>%
  rename(c("lower95" = "lowerlong$lower", "upper95" = "upperlong$upper", "lower90" = "lowerlong90$lower90", "upper90" = "upperlong90$upper90"))

MedLowUplong$year10 <- as.numeric(MedLowUplong$year)/10
MedLowUplong <- MedLowUplong %>%
  mutate(year = ifelse(year10 == 190.0, 1900, year10))
MedLowUplong$year10 <- NULL

MedLowUplong$country_text_id <- factor(MedLowUplong$country_text_id,
                                       labels = c("Afghanistan", "Angola", "Albania", "United Arab Emirates", "Argentina", "Armenia", "Australia", "Austria", "Azerbaijan", "Burundi", "Belgium", "Benin", "Burkina Faso", "Bangladesh", "Bulgaria", "Bahrain", "Bosnia and Herzegovina", "Belarus", "Bolivia", "Brazil", "Barbados", "Bhutan", "Botswana", "Central African Republic", "Canada", "Switzerland", "Chile", "China", "Ivory Coast", "Cameroon", "Democratic Republic of the Congo", "Republic of the Congo", "Colombia", "Comoros", "Cape Verde", "Costa Rica", "Cuba", "Cyprus", "Czech Republic", "German Democratic Republic", "Germany", "Djibouti", "Denmark", "Dominican Republic", "Algeria", "Ecuador", "Egypt", "Eritrea", "Spain", "Estonia", "Ethiopia", "Finland", "Fiji", "France", "Gabon", "United Kingdom", "Georgia", "Ghana", "Guinea", "The Gambia", "Guinea-Bissau", "Equatorial Guinea", "Greece", "Guatemala", "Guyana", "Hong Kong", "Honduras", "Croatia", "Haiti", "Hungary", "Indonesia", "India", "Ireland", "Iran", "Iraq", "Iceland", "Israel", "Italy", "Jamaica", "Jordan", "Japan", "Kazakhstan", "Kenya", "Cambodia", "South Korea", "Kuwait", "Laos", "Lebanon", "Liberia", "Libya", "Sri Lanka", "Lesotho", "Lithuania", "Luxembourg", "Latvia", "Morocco", "Moldova", "Madagascar", "Maldives", "Mexico", "Macedonia", "Mali", "Burma/Myanmar", "Montenegro", "Mongolia", "Mozambique", "Mauritania", "Mauritius", "Malawi", "Malaysia", "Namibia", "Niger", "Nigeria", "Nicaragua", "Netherlands", "Norway", "Nepal", "New Zealand", "Oman", "Pakistan", "Panama", "Peru", "Philippines", "Papua New Guinea", "Poland", "North Korea", "Portugal", "Paraguay", "Palestine-British Mandate", "Palestine-West Bank", "Palestine-Gaza", "Qatar", "Romania", "Russia", "Rwanda", "Saudi Arabia", "Sudan", "Senegal", "Singapore", "Solomon Islands", "Sierra Leone", "El Salvador", "Somaliland", "Somalia", "Serbia", "South Sudan", "Sao Tome and Principe", "Suriname", "Slovakia", "Slovenia", "Sweden", "Swaziland", "Seychelles", "Syria", "Chad", "Togo", "Thailand", "Tajikistan", "Turkmenistan", "Timor Leste", "Trinidad and Tobago", "Tunisia", "Turkey", "Taiwan", "Tanzania", "Uganda", "Ukraine", "Uruguay", "United States", "Uzbekistan", "Vietnam, Democratic Republic of", "Venezuela", "Vietnam", "Vanuatu", "Kosovo", "Yemen", "Democratic Yemen", "South Africa", "Zambia", "Zimbabwe", "Zanzibar"))

saveRDS(MedLowUplong, "simulated_effects_upA.rds")
## saveRDS(result, "simulated_effectsN.rds")
