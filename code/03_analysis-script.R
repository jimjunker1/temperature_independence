# analysis script #
source("./code/02_data-import-clean.R")

## mean temperatures
ann_stream_temps
# A tibble: 6 x 3
# site_id tempkt tempC
# <fct>    <dbl> <dbl>
# hver      38.6  27.2
# st6       39.9  17.6
# st9       40.8  11.2
# st7       41.6   5.8
# oh2       41.6   5.5
# st14      41.7   5  

## correlations among stream temp-light: Table S1 ##
cor(overkt_to_C(tempkt_light.doy$Hver_tempC), tempkt_light.doy$Intensity)#0.67
cor(overkt_to_C(tempkt_light.doy$st9_tempC), tempkt_light.doy$Intensity)#0.48
cor(overkt_to_C(tempkt_light.doy$st6_tempC), tempkt_light.doy$Intensity)#-0.09
cor(overkt_to_C(tempkt_light.doy$st14_tempC), tempkt_light.doy$Intensity)#0.47
cor(overkt_to_C(tempkt_light.doy$st7_tempC), tempkt_light.doy$Intensity)#0.41
cor(overkt_to_C(tempkt_light.doy$oh2_tempC), tempkt_light.doy$Intensity)#0.57

## Mean estimates and bounds of annual P across streams ##
annual_summ %>% select(site_id, tempkt, prod_mg_m_y, prod_quant2.5, prod_quant97.5) %>% arrange(tempkt)
# A tibble: 6 x 5
# site_id tempkt prod_mg_m_y prod_quant2.5 prod_quant97.5
# <fct>    <dbl>       <dbl>         <dbl>          <dbl>
# hver      38.6      19931.        14872.         25187.
# st6       39.9       3966.         3225.          4692.
# st9       40.8       3037.         1836.          4596.
# oh2       41.6       3529.         3192.          3872.
# st7       41.6       7330.         6627.          8074.
# st14      41.7        451.          256.           703.

## Estimated temperature dependence of sec prod among streams ##
mean(prod_boot_df$tempkt_stand)#Epann = -0.67
quantile(prod_boot_df$tempkt_stand, c(0.025,0.5,0.975))
# Epann 0.57, 0.6746, 0.78

## Estimated change per degree C ##
exp(mean(prod_boot_C$tempC))#1.099
exp(quantile(prod_boot_C$tempC, c(0.025,0.975)))#1.084  1.115

## Estimated temperature dependence of size corrected B among streams ##
mean(bio_boot_df$tempkt_stand)#Epann = -0.91
quantile(bio_boot_df$tempkt_stand, c(0.025,0.5,0.975))
# Ebann -1.09 -0.90 -0.75

## Estimated temperature dependence of size corrected PB among streams ##
mean(pb_boot_df$tempkt_stand)#Epann = 1.52
quantile(pb_boot_df$tempkt_stand, c(0.025,0.5,0.975))
# Epbann 1.34 1.53 1.67

#### Within stream prod patterns ####
## within stream minimum daily prod ##
int_df %>% arrange(prod_d) %>% distinct(site_id, .keep_all = TRUE) %>%  arrange(mean_tempkt) %>% select(site_id,prod_d)
# site_id   prod_d
# hver      13.20799
#  st6       2.24333
#  st9       0.56808
#  oh2       0.54034
#  st7       6.44632
# st14       0.08285

## within stream maximum daily prod ##
int_df %>% arrange(-prod_d) %>% distinct(site_id, .keep_all = TRUE) %>% arrange(mean_tempkt) %>% select(site_id, prod_d)
# site_id    prod_d
#  hver      188.63583
#   st6       29.10403
#   st9       25.13090
#   st7       54.42227
#   oh2       35.35611
#  st14        4.32870

## the apparent temperature dependence of within stream P ##
prod_d_stream_coef_df %>% select(site_id, prod_slope_est, prod_quant_2.5, prod_quant_97, prod_r2_avg) 
# site_id prod_slope_est prod_quant_2.5 prod_quant_97 prod_r2_avg
#  hver      1.6993274      1.3266608      2.121315   0.7871004
#   st6      0.8207585      0.5231223      1.114395   0.0720197
#   st9      3.8565432      3.0940706      4.542700   0.8245517
#   st7      4.0865444      3.6683887      4.527877   0.7900215
#   oh2      2.6746914      2.5430067      2.804918   0.6089898
#  st14      1.7017471      0.9042982      2.349233   0.4974649

## the apparent temperature dependence of within stream corrected B##
bio_stream_coef_df %>% select(site_id, bio_slope_est, bio_quant_2.5, bio_quant_97, bio_r2_avg) 
# site_id bio_slope_est bio_quant_2.5 bio_quant_97 bio_r2_avg
#  hver     2.8998883    2.50126385     3.315717  0.7791279
#   st6     0.6009578   -0.14078399     1.443333  0.0530126
#   st9     2.2284265   -0.02957069     4.595299  0.2379404
#   st7     2.5920669    1.00451845     3.996570  0.2718937
#   oh2     0.9007695    0.57758775     1.240454  0.0920658
#  st14     1.5504234    0.04146219     2.923530  0.1332558

## the apparent temperature dependence of within stream corrected PB ##
pb_stream_coef_df %>% select(site_id, pb_slope_est, pb_quant_2.5, pb_quant_97, pb_r2_avg)
# site_id pb_slope_est pb_quant_2.5  pb_quant_97 pb_r2_avg
#  hver  -1.71488351   -2.2260548 -1.174894227 0.3666217
#   st6   0.08937757   -0.6673858  0.743688861 0.0350413
#   st9   1.00411955   -1.1034463  3.147449548 0.1232841
#   st7   0.31048068   -0.9920129  1.862352531 0.0406347
#   oh2   1.80380073    1.4860550  2.111168611 0.2426608
#  st14  -1.42681741   -2.7638153 -0.005966353 0.0863348

######  ACCOUNTING FOR RESOURCE VARIATION  ######
## Estimating the temperature dependence of chla among streams ##
### bootstrapped chla data ###
chla_temp_coefs_df$tempkt_est#0.53
chla_temp_coefs_df[,c("tempkt_quant2.5","tempkt_quant97.5")]
# tempkt_quant2.5 tempkt_quant97.5
#  0.3918432        0.6886159
chla_temp_coefs_df$adj_r2_mean#0.52

## chlorophyll a vs P ##
na.rm_mean(chla_prod_lmboots$chla)#1.19
quantile(chla_prod_lmboots$chla, c(0.025,0.975))
# 2.5%     97.5% 
# 0.8816   1.5770

## model selection on annual prod across streams ##
ann_mod_sel_df = ann_mod_sel_boot(annual_df, 'prod_y', nboot = 1000)
ann_mod_summary = ann_mod_sel_df %>% group_by(Model) %>% count() %>% ungroup() %>% arrange(desc(n))
ann_mod_summary
## A tibble: 4 x 2
# Model       n
# <chr>   <int>
# chla.lm   853
# full.lm   112
# add.lm     27
# temp.lm     8

#### 

## resource corrected annual P ##
na.rm_mean(chlatemp_prod_lmboots$tempkt_stand)#-0.05
quantile(chlatemp_prod_lmboots$tempkt_stand, c(0.025,0.975))
#     2.5%      97.5% 
# -0.3743697  0.3970664

## model selection on seasonal production within streams ##
mod_sel_summ = mod_sel_list[[1]] %>% group_by(Model) %>% count() %>% ungroup() %>% mutate(n = rev(sort(n)))
mod_sel_summ
# A tibble: 3 x 2
# Model             n
#   <chr>         <int>
# 1 prod_full_add   384
# 2 prod_full_int   321
# 3 temp_prod.int   295

## refit with REML is top_model_centerfit object
# summarise the goodness of fit
top_mod_summ = top_mod_centerfit[[1]] %>% select(contains("r2")) %>% mutate_all(as.numeric) %>%summarise_all(mean)
top_mod_summ
#      r2_m     r2_c
#  0.5701477 0.886565

#estimated apparent temperature dependence of sec prod within streams
na.rm_mean(prod_Ea_vec)#1.52
quantile(prod_Ea_vec, probs = c(0.025,0.5,0.975))
# 1.04 1.54 1.91
## estimated EaGPP from bootstrapped analysis
# full model in supplementary analysis
na.rm_mean(gpp_temp_vec)#1.37
quantile(gpp_temp_vec, probs = c(0.025,0.5, 0.975))
#    2.5%      50%    97.5% 
# 1.362931 1.372968 1.383397 

# Estimated within stream resource-adjusted temperature dependence 
na.rm_mean(gpp_adj_prod)#0.14
quantile(gpp_adj_prod, probs = c(0.025,0.5,0.975))
#    2.5%        50%      97.5% 
# -0.3340680  0.1638090  0.5322748

### Supplementary analysis
#Values for Table s2
gpp_chla_vec = sapply(gpp_boot_list, function(x) summary(x)[["coefficients"]]["log(chla_mean)","Estimate"])
na.rm_mean(gpp_chla_vec)#
quantile(gpp_chla_vec, probs = c(0.025,0.5,0.975))
#    2.5%        50%      97.5% 
# -0.1498082 -0.1475098 -0.1452362
gpp_temp_vec = sapply(gpp_boot_list, function(x) summary(x)[["coefficients"]]["tempkt_stand","Estimate"])
na.rm_mean(gpp_temp_vec)#
quantile(gpp_temp_vec, probs = c(0.025,0.5, 0.975))
#
gpp_light_vec = sapply(gpp_boot_list, function(x) summary(x)[["coefficients"]]["log(mean_light)","Estimate"])
na.rm_mean(gpp_light_vec)#
quantile(gpp_light_vec, probs = c(0.025,0.5,0.975))
# 
