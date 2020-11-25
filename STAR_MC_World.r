library(readstata13)
library(tidyr)
library(MASS)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(SparseM)

## library(spdep)
## library(spatialreg)

## Choose either of the directory ##
setwd("~/Intake/Coppedge/ForMC")
setwd("C:/Users/mcoppedg/Dropbox/Diffusion/data/W/Levels")
####################################

## There are two ways to prepare the required data for simulation ##
## One: Read the un-row-standardized W_N and W_A matrix ##
W_A= read.dta13("W_A.dta") %>% select(-"v1", -'id')
W_N= read.dta13("W_N.dta") %>% select(-'id')

country_id = W_A %>%
    pull(country_text_id) %>%
    unique() %>% sort()

year_id = W_N %>%
  pull(seq) %>%
  unique() %>%
  sort()

W_N_shaped = data.frame()
W_A_shaped = data.frame()
for ( i in year_id){
   temp1 = W_N %>% filter(str_detect(seq, as.character(i))) %>%
      select(matches(paste0(as.character(i),'|cyr'))) %>%
      arrange(cyr)
   N_names = temp1 %>% colnames()
   N_names = gsub("[0-9]","", N_names)
   colnames(temp1) = N_names
   W_N_shaped = W_N_shaped  %>% bind_rows(temp1)
   temp2 = W_A %>% filter(str_detect(seq, as.character(i))) %>%
      select(matches(paste0(as.character(i),'|cyr'))) %>%
      arrange(cyr)
   A_names = temp1 %>% colnames()
   A_names = gsub("[0-9]","", A_names)
   colnames(temp2) = A_names
   W_A_shaped = W_A_shaped %>% bind_rows(temp2)
}

cyr = W_N_shaped %>% select(cyr)
W_N_shaped = W_N_shaped %>%
  select(-cyr) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate_all(~ ./row_sum) %>%
  select(-row_sum) %>%
  bind_cols(all_of(cyr)) %>%
  select(cyr, everything())

cyr = W_A_shaped %>% select(cyr)
W_A_shaped = W_A_shaped %>%
  select(-cyr) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate_all(~ ./row_sum) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(-row_sum) %>%
  bind_cols(all_of(cyr)) %>%
  select(cyr, everything())


saveRDS(list(country_id, year_id, cyr, W_N_shaped,W_A_shaped),"IDs_Weight_matrices.rds")
###################################

## Code to plot the sparse matrices ##
## Might not function for now due to the extreme large matrix ##
W_N = W_N %>% select(-cyr, -country_text_id, -seq) %>% select(-10861) %>% as.matrix()
W_A = W_A %>% select(-cyr, -seq) %>% as.matrix()
neighbor_list <- mat2listw(x=W_N,style="W") #Generate neighbor list, "W" = row standardized
W_N <- listw2mat(neighbor_list)
neighbor_list <- mat2listw(x=W_A,style="W") #Generate neighbor list, "W" = row standardized
W_A <- listw2mat(neighbor_list)
image(as.matrix.csr(W_N))
image(as.matrix.csr(W_A))
###############################


## Two: load RDS file that saved from the method one and use it directly ##

weight_matrix = readRDS("IDs_Weight_matrices.rds")
country_id = weight_matrix[[1]]
year_id = weight_matrix[[2]]
cyr = weight_matrix[[3]]
W_N_shaped = weight_matrix[[4]]
W_A_shaped = weight_matrix[[5]]

## Preparation for simulation
vcvm = read.dta13("VCVM.dta")[1:3,1:3] %>% as.matrix()

rho_W_N = 0.3577922
rho_W_A = 0.0015372
phi = 0.6500221
#theta_W_N = -0.2424203
#theta_W_A = 0.0252452
num_draws <- 1000  #number of simulations to draw
coefficient_star = c(rho_W_N, rho_W_A, phi)


shock_input = data.frame(matrix(ncol=length(country_id), nrow=length(year_id)))
shock_input = cbind(year_id, shock_input)
colnames(shock_input) = c('Year', country_id)
write.csv(shock_input,'shock_input.csv', row.names = FALSE)

shock_input_manual = read.csv('shock_input.csv', stringsAsFactors = FALSE) %>%
  filter_at(vars(!contains('Year')), any_vars(!is.na(.)))
shock_brief = shock_input_manual %>% melt(id='Year') %>%
  drop_na('value') %>%
  arrange(Year)

shock_input_manual = shock_input_manual %>%
                     mutate_all(~replace(., is.na(.), 0))

print('The briefing of input shock is: ')
print(shock_brief)

draws <- mvrnorm(n=num_draws,coefficient_star,vcvm,empirical=TRUE)

#Designate unit experiencing counterfactual shock
#All 1s:
#The initial shock

#Actual changes from 1900 to 1901.5:

years_with_shock = shock_input_manual %>% pull(Year)

## shock_in_matrix = diag(c(0, -0.0095, -0.0035, -0.0225, 0, -0.015, 0.09, 0, 0, 0, -0.021, 0),12,12)
#Change of 1 in Argentina only:
## shock_in_matrix = diag(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),12,12)

country_num = length(country_id)
identity_m <- diag(1,country_num,country_num)  #nxn identity matrix
simulated_effects_vector = vector("list", num_draws)
simulated_effects_matrix = vector("list", num_draws)

#i is the draw number
#j is the country number, used to define the first and last row and
# column numbers in each submatrix below.
year_start = years_with_shock[1]
shock_vector = shock_input_manual %>% filter(Year==year_start) %>% select(-Year) %>% slice(1) %>% as.numeric
shock_matrix = diag(shock_vector)
year_num = length(year_id[year_id > year_start])
##single_result = matrix(shock_vector,nrow=1,ncol=country_num)
single_result = matrix(0,nrow=1,ncol=country_num)
## single_result_long_dataframe = as.data.frame(diag(shock_vector))
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

for (i in 1:2)
{

  for (j in (year_start_position+1):60) {
      current_year = current_year+20
      y_star = draws[i,3]*(shock_vector+y_real)
      y_star_matrix = draws[i,3]*(shock_matrix + y_real_matrix)
      ## start_index = (j-1)*country_num+1
      ## end_index = j*country_num
      W_N_current = W_N_shaped %>% filter(str_detect(cyr, as.character(current_year))) %>% select(-cyr) %>% as.matrix()
      W_A_current = W_A_shaped %>% filter(str_detect(cyr, as.character(current_year))) %>% select(-cyr) %>% as.matrix()
      LRSS_multiplier = solve(identity_m - draws[j,1]*W_N_current- draws[j,2]*W_A_current)
      ## LRSS_multiplier = solve(identity_m - draws[j,1]*W_N_SAm[start_index:end_index, start_index:end_index]- draws[j,2]*W_A_SAm[start_index:end_index, start_index:end_index])
      delta_multiplier = draws[i,1]*W_N_current + draws[i,2]*W_A_current
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

    simulated_effects_vector[[i]] = single_result
    simulated_effects_matrix[[i]] = single_result_long_dataframe
    year_start = years_with_shock[1]
    shock_vector = shock_input_manual %>% filter(Year==year_start) %>% select(-Year) %>% slice(1) %>% as.numeric
    shock_matrix = diag(shock_vector)
    country_num = length(country_id)
    year_num = length(year_id[year_id > year_start])
    single_result = matrix(0,nrow=1,ncol=country_num)
    single_result_long_dataframe = as.data.frame(diag(rep(0,country_num)))
    y_real = rep(0, country_num)
    y_real_matrix = diag(y_real)
    year_start_position = which(year_id == year_start)
    if (year_start == 1900){
      current_year = 18995
    }else{
      current_year = year_start
    }
}


#This uses just the first draw from simulated_effects_vector.
dataframe = simulated_effects_vector[[1]] %>% as.data.frame()

## Due to the large number of countries, plotting the effect in one single image is slow and dense.
## Thus I set up two way to plot
## One: select certian coutries in the line below, in the parenthesis, Year should always be selected
##      any other countries selected will be plotted, this is good to manually pick countries and check the effects
dataframe %>% select(Year, afg:ben) %>% melt(id='Year') %>%
  mutate(country = as.character(variable)) %>%
  mutate(Year = ifelse(Year==1900, 19000, Year)) %>%
  ggplot(aes(x=as.numeric(Year)/10, y=value, color='Grey')) +
  geom_line(size=1) +
  facet_wrap(~ variable, ncol = 4, scales = "fixed")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color = "white"),
        panel.background = element_rect(fill="gray95")) +
  ## scale_x_continuous(limits=c(year_start/10,20175), breaks = seq(year_start/10,2017,2)) +
  ## scale_y_continuous(limits=c(0,.12), breaks = seq(0,.12,.02)) +
  labs(x="Year", y="Effect", title="Effect of Shocks in Multiple Country/Period") +
  theme(axis.text.x=element_text(angle=45, hjust=1), plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), legend.position='none')
ggsave("WorldSimulationInOne.png", width=49,height=49)

    #This just adds the country names and years to dataframe.
    ## colnames(dataframe) = c(country_id, 'Year')
    ## dataframe = dataframe %>% select(c('Year', all_of(country_id)))

    #This melts the df to year-country observations with the variables
    # year, country, and value; and makes a facet graph:


## Two: loop through all countries and plot every 12 countries in one image

for (i in 1:15){

dataframe %>% select(Year, (2+(i-1)*12):(12*i+1)) %>%
  melt(id='Year') %>%
  mutate(country = as.character(variable)) %>%
  mutate(Year = ifelse(Year==1900, 19000, Year)) %>%
  ## ggplot(aes(x=Year, y=value, color='Grey')) +
  ggplot(aes(x=as.numeric(Year)/10, y=value, color='Grey')) +
  geom_line(size=1) +
  facet_wrap(~ variable, ncol = 4, scales = "fixed")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color = "white"),
        panel.background = element_rect(fill="gray95")) +
  ## scale_x_continuous(limits=c(year_start/10,20175), breaks = seq(year_start/10,2017,2)) +
  ## scale_y_continuous(limits=c(0,.12), breaks = seq(0,.12,.02)) +
  labs(x="Year", y="Effect", title="Effect of Shocks in Multiple Country/Period") +
  theme(axis.text.x=element_text(angle=45, hjust=1), plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), legend.position='none') -> plot1
  print(plot1)
}

dataframe %>% select(Year, 182) %>%
  melt(id='Year') %>%
  mutate(country = as.character(variable)) %>%
  mutate(Year = ifelse(Year==1900, 19000, Year)) %>%
  ## ggplot(aes(x=Year, y=value, color='Grey')) +
  ggplot(aes(x=as.numeric(Year)/10, y=value, color='Grey')) +
  geom_line(size=1) +
  facet_wrap(~ variable, ncol = 4, scales = "fixed")+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_line(color = "white"),
        panel.background = element_rect(fill="gray95")) +
  ## scale_x_continuous(limits=c(year_start/10,20175), breaks = seq(year_start/10,2017,2)) +
  ## scale_y_continuous(limits=c(0,.12), breaks = seq(0,.12,.02)) +
  labs(x="Year", y="Effect", title="Effect of Shocks in Multiple Country/Period") +
  theme(axis.text.x=element_text(angle=45, hjust=1), plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), legend.position='none') -> plot1
  print(plot1)

ggsave("LatinAmerican12.png", width=49,height=49)

#This gathers the draws for one country, all years into a data frame.
#To change the country, change the column number in the
# simulated_effects_vector[[i]][,k])
#line.
#Key: 1 "arg" 2 "bol" 3 "bra" 4 "chl" 5 "col" 6 "ecu" 7 "guy"
# 8 "per" 9 "pry" 10 "sur" 11 "ury" 12 "ven"

#data_both_all = data.frame(year=character(), median.x=double(), lower.x=double(), upper.x=double(),
#median.y=double(), lower.y=double(), upper.y=double(), cnt=character(), stringsAsFactors = FALSE)
datalist <- list()

for (k in 1:country_num) {
  one_country_all_simulations = c()

for (i in 1:10)
  {
  temp = simulated_effects_vector[[i]][,c(1,(k+1))] %>%
    as.data.frame() %>% select(-Year) %>% t()
  one_country_all_simulations = rbind(one_country_all_simulations,temp)
  }

#This labels the years (in columns)
colnames(one_country_all_simulations) = year_id

#This reshapes the one-country draws & saves as df data_to_plot.
#Var1 is the draw number, Var2 is the year, Value is the MRP,
# and L1 is the line color for each estimate.
data_to_plot = one_country_all_simulations %>%
  melt() %>%
  mutate(Var2 = as.character(Var2))

#Note: Below this I have overwritten some of Yang's code.
#This calculates the median MRP across all draws for one country,
# reshapes it from wide to long, creates Var1 = 1001 for all years,
# and saves it all as new object median_line.
one_country_all_simulations %>% apply(2,median) %>%
  melt() %>% tibble::rownames_to_column("year") %>%
  rename(c("median" = "value")) -> median_line

one_country_all_simulations %>%
  apply(2, quantile, probs=c(0.025)) %>%
  melt() %>% tibble::rownames_to_column("year") %>%
  rename(c("lower" = "value")) -> lowerbound

one_country_all_simulations %>%
  apply(2, quantile, probs=c(0.975)) %>%
  melt() %>% tibble::rownames_to_column("year") %>%
  rename(c("upper" = "value")) -> upperbound

#This merges the three files
bounds <- merge(lowerbound, upperbound, by=c("year"))
data_to_plot2 <- merge(median_line, bounds, by=c("year"))

#Calculating the cumulative sum
#one_country_all_simulations %>%
  OCAS_t <-   as.data.frame(x = t(one_country_all_simulations), stringsAsFactors = FALSE)
  OCAS_t[1,] <- 0
  COCAS_t <- cumsum(OCAS_t)
  COCAS <-  as.data.frame(x = t(COCAS_t), stringsAsFactors = FALSE)

COCAS %>% apply(2,median) %>%
    melt() %>% tibble::rownames_to_column("year") %>%
    rename(c("median" = "value")) -> c_median_line

COCAS %>%
  apply(2, quantile, probs=c(0.025)) %>%
  melt() %>% tibble::rownames_to_column("year") %>%
  rename(c("lower" = "value")) -> c_lowerbound

COCAS %>%
  apply(2, quantile, probs=c(0.975)) %>%
  melt() %>% tibble::rownames_to_column("year") %>%
  rename(c("upper" = "value")) %>%
  mutate(cnt = k) -> c_upperbound

#This merges the three files
c_bounds <- merge(c_lowerbound, c_upperbound, by=c("year"))
c_data_to_plot <- merge(c_median_line, c_bounds, by=c("year"))

#This merges the yearly and cumulative effects for country k
data_both <- merge(data_to_plot2, c_data_to_plot, by=c("year"))

#This is supposed to append the new country's median and
#percentiles to the accumulated values for countries that
#preceed it alphabetically, but it keeps only Venezuela.
datalist[[k]] <- data_both
#data_both_all <- as.data.frame(rbind(data_both_all, data_both))
}

data_both_all = do.call(rbind, datalist)


data_both_all$cnt <- factor(data_both_all$cnt,
                        labels = c("Argentina", "Bolivia", "Brazil",
                        "Chile", "Colombia", "Ecuador", "Guyana",
                        "Peru", "Paraguay", "Suriname", "Uruguay", "Venezuela"))


########     PLOTS      ##################
#Plotting both together
ggplot(data = subset(data_both_all, year>1901 & year<2017)) +
  geom_ribbon(aes(x=as.numeric(year)/10,
                  ymin = lower.x, ymax = upper.x, group=1),
              fill = "mediumblue", alpha=.2) +
  geom_ribbon(aes(x=as.numeric(year)/10,
                  ymin = lower.y, ymax = upper.y, group=1),
              fill = "firebrick", alpha=.2) +
  geom_line(aes(x=as.numeric(year)/10, y = median.x, group=1),
            color="mediumblue", size=1) +
  geom_line(aes(x=as.numeric(year)/10, y = median.y, group=1),
            color="firebrick", size=1) +
  facet_wrap(~ reorder(cnt, -median.y), ncol = 3, scales = "fixed")+
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(color = "white"),
        panel.background = element_rect(fill="gray97"),
        strip.text = element_text(face = "bold"),
        plot.background = element_rect(fill = "cornsilk")) +
  scale_x_continuous(limits=c(1901,2017), breaks = seq(1902,2017,2)) +
  labs(x="", y="Effect", caption="95% of the estimates lie within the shadows.",
       title="LRSS effect of a unit shock from Argentina 1900. Red is cumulative, Blue is annual increment.")

####Maps#############
library(maps)
library(mapdata)
library(mapproj)
world_map <- map_data("world")
SAm <- map_data("world",
      region =c("Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
                "Ecuador", "Guyana", "Paraguay", "Peru", "Suriname",
                "Uruguay", "Venezuela"))

SAmmap <- merge(SAm, data_both_all, by.x = "region", by.y = "cnt")
SAmmap <- arrange(SAmmap, region)
View(SAmmap)

# Specify year for seq
specyear <- c(1900, 19015, 19035, 19055, 19075, 19095, 19115)


  ggplot(filter(SAmmap, year==1900), aes(x=long, y=lat,
                        group=group, fill = median.x)) +
  geom_polygon(color="black") +
  coord_map(projection = "mercator") +
  theme_void() +
    labs(title="1900", size=3) +
  scale_fill_gradient(low="white", high = "darkblue",
                      limits = c(0,.15),
                      name="incremental\neffect")
    ggsave(file="Map0.png")


############################
#Superceded code:
ggplot(data = subset(data_to_plot2, year>1901 & year<1913)) +
  geom_ribbon(aes(x=as.numeric(year)/10,
                  ymin = lower, ymax = upper, group=1),
              fill = "cornsilk") +
  geom_line(aes(x=as.numeric(year)/10, y = median, group=1),
            color="firebrick", size=1) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(color = "gray90"),
        panel.background = element_rect(fill="white")) +
  scale_x_continuous(limits=c(1901,1912), breaks = seq(1902,1912,2)) +
  labs(x="Uruguay", y="marginal response path")

ggplot(data = subset(c_data_to_plot, year>1901 & year<1913)) +
  geom_ribbon(aes(x=as.numeric(year)/10,
                  ymin = lower, ymax = upper, group=1),
              fill = "cornsilk") +
  geom_line(aes(x=as.numeric(year)/10, y = median, group=1),
            color="firebrick", size=1) +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_line(color = "gray90"),
        panel.background = element_rect(fill="white")) +
  scale_x_continuous(limits=c(1901,1912), breaks = seq(1902,1912,2)) +
  scale_y_continuous(limits=c(0,.16)) +
  labs(x="Uruguay", y="Cumulative effect on")
#Yang's code:
#This plots the many draws and the medians in the same line graph.
data_to_plot %>%
  ggplot(aes(Var2, value, group = factor(Var1),stat='density',
             color=factor(L1), size=factor(L1),alpha=factor(L1))) +
  geom_line() +
  scale_color_manual(values=c('grey','black')) +
  scale_size_manual(values=c(0.4,0.8)) +
  scale_alpha_manual(values=c(0.1,0.9)) +
  theme_classic() +
  labs(x="Year", y="Effect", title="ARG Effect in 1000 Simulations") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank(), legend.position="none")
ggsave("ARG1000Simulations.png", width=16,height=8)

dataframe %>% write_csv("12_countries_in_one_simulation.csv")

  ## X_no_constant = seqtrial_MI_SAm %>% filter(stringr::str_detect(cyr, toString(year)))
  ## X_no_constant = X_no_constant[,X_columns]
  ## y_star = draws[i,3]*shock + as.matrix(X_no_constant) %*% as.matrix(beta) + constant
  ## y_star = phi * rep(1,12)
  ## LRSS_multiplier = solve(identity_m - draws[i,1]*W_N_SAm[1:12,1:12] - draws[i,2]*W_A_SAm[1:12,1:12] )
  ## delta_multiplier = draws[i,1]*W_N_SAm[1:12,1:12] + draws[i,2]*W_A_SAm[1:12,1:12]
  ## delta_temp = LRSS_multiplier %*% delta_multiplier %*% y_star
