source("./code/01_supporting-functions.R")
### derived data objects  ###
# data frame ordering variables.
stream_order = factor(c("hver", "st6","st9", "st7","oh2","st14"))#stream ordering
stream_order_list = stream_order %>% as.list() %>% setNames(.,stream_order) #stream ordering for lists
stream_temp_labels = c("27.2","17.6","11.2","5.8","5.5","5.0")#stream annual temperature labels
names(stream_temp_labels) = stream_order #setting named character vector of stream labels
#temperature data
#Temp light for all streams
tempkt_full = read.table(file = "./data/raw-data/tempkt_full.txt", header = TRUE, sep = "\t") # all streams temp file for study period
tempkt_light.doy = read.table(file = "./data/raw-data/doy_tempkt_light.txt", header = TRUE, sep = "\t")
tempkt_light.doy = na.omit(tempkt_light.doy)

light.doy_long = tempkt_light.doy %>%
  select(doy, Intensity) %>%
  group_by(doy) %>% 
  summarise(Intensity = na.rm_mean(Intensity, na.rm = TRUE)) %>%
  ungroup()
temp_light.doy = tempkt_light.doy %>% 
  mutate_at(vars(contains("tempC")), overkt_to_C) 

light.doy_long = tempkt_light.doy %>%
  select(doy, Intensity) %>%
  group_by(doy) %>% 
  summarise(Intensity = na.rm_mean(Intensity, na.rm = TRUE)) %>%
  ungroup()

temp.doy_long = temp_light.doy %>% select(-Intensity, -light) %>%
  gather(stream, temperature, contains("tempC")) %>%
  group_by(stream,doy) %>% summarise(temperature = na.rm_mean(temperature, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(stream = str_remove(tolower(stream), "_tempc")) %>%
  mutate(stream = factor(stream, levels = stream_order))

ann_stream_temps <- readRDS(file = "./data/derived-data/annual_tempkt.rds") %>% #annual mean boltzmann temperature
 pivot_longer(contains('temp'),names_to = 'site_id', values_to = 'tempkt') %>%
  mutate(site_id = tolower(gsub("(.*)_.*$", "\\1", site_id)),
         site_id = factor(site_id, levels = stream_order),
         tempC = round(overkt_to_C(tempkt),1),
         tempkt = round(tempkt, 1)) %>%
  arrange(tempkt) 

stream_sampling_meta = read.table("./data/raw-data/stream_sampling_meta.txt", header = TRUE, sep = "\t") %>% rename(DATE = 'date_id')#stream_temps_bound
stream_temps_list = stream_sampling_meta %>% junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list)) #interval temperature list

####
chla = read.table("./data/raw-data/chla_samples.txt", header = TRUE, sep = "\t") %>% #full chla file
  rename(DATE = 'date_id') %>%
  mutate(JULIAN = julian(as.Date(DATE), origin = as.Date("2010-01-01")),
         DATE = as.Date(DATE)) #convert dates to JULIAN

# stream_sampling_meta
# stream_temps = list()
# sites_list = unique((chla$site_id))
# debugonce(temp_subset)
# for(site in sites_list){#apply temp_subset function to all streams 
#   TEMP = tempkt_full %>% select(DATE,JULIAN, contains(as.character(site))) %>% mutate(DATE = as.Date(DATE))
#   bugs_df =  chla %>% filter(site_id == site)
#   stream_temps[[site]] <- temp_subset(bugs_df, TEMP, wrap = TRUE)
# }
# rm(TEMP)
# stream_temps = map(stream_temps, ~.x %>% mutate(days = c(diff(JULIAN),NA)))
# stream_temps = map2(stream_temps, sites_list, ~.x %>% mutate(site_id = rep(.y, nrow(.x))))

chla_summ = chla %>% group_by(site_id, DATE) %>%
  summarise(mean_chla_mg_m = na.rm_mean(chla_mg_m, na.rm = TRUE), doy = unique(doy)) %>% ungroup()%>%
  left_join(stream_sampling_meta %>% select(site_id, DATE, contains('temp')) %>% mutate(DATE = as.Date(DATE))) %>%
  mutate(site_id = factor(site_id, levels = stream_order)) %>% na.omit()

chla_ann_summ = chla_summ %>% group_by(site_id) %>%
  summarise(chla_mean_ann = na.rm_mean(mean_chla_mg_m, na.rm = TRUE),
            chla_sd_ann = sd(mean_chla_mg_m, na.rm = TRUE)) %>%
  left_join(ann_stream_temps)


set.seed(123)
chla_boot_list = chla %>% select(site_id, DATE, chla_mg_m) %>% group_by(site_id, DATE) %>%#keep @ 1000 to merge with prod,bio,pb
  sample_n(1000, replace = TRUE) %>% ungroup() %>% mutate(site_id = factor(site_id, levels = stream_order)) %>%
  group_split(site_id) %>% set_names(chla %>% group_keys(site_id) %>% pull(1)) %>%
  rlist::list.subset(names(stream_order_list)) %>%
  map(., ~.x %>% group_by(DATE) %>% mutate(id = 1:n()) %>% ungroup() %>%
        pivot_wider(names_from = DATE, values_from = chla_mg_m) %>%
        select(-id) %>% transmute(site_id,chla_ann_mean = rowMeans(select(.,-site_id), na.rm = TRUE))) 

set.seed(123)
chla_int_df = chla %>% select(site_id, DATE, chla_mg_m) %>% group_by(site_id, DATE) %>%
  sample_n(1000, replace = TRUE) %>% mutate(id = 1:n()) %>% 
  pivot_wider(id_cols = c(site_id,DATE), values_from = chla_mg_m, names_from = id ) %>% ungroup() %>% mutate(site_id = factor(site_id, levels = stream_order))

chla_boot_df= chla_boot_list %>%
  bind_rows() %>% group_by(site_id) %>%
  mutate(boot_id = paste0("X",1:n())) %>% ungroup %>%
  arrange(factor(site_id, rev(stream_order)))


chla_summ_ann = map(chla_boot_list, ~.x %>%
                      summarise(site_id = unique(site_id),
                                ann_chla_mean = na.rm_mean(chla_ann_mean, na.rm = TRUE),
                                ann_chla_quant2.5 = quantile(chla_ann_mean, 0.025, na.rm = TRUE),
                                ann_chla_quant25 = quantile(chla_ann_mean, 0.25, na.rm = TRUE),
                                ann_chla_quant75 = quantile(chla_ann_mean, 0.75, na.rm = TRUE),
                                ann_chla_quant97 = quantile(chla_ann_mean, 0.975, na.rm = TRUE),
                                ann_chla_min = min(chla_ann_mean, na.rm = TRUE),
                                ann_chla_max = max(chla_ann_mean, na.rm = TRUE))) %>%
  bind_rows()
rm(chla_boot_list)
#bootstrapped chla temp relationships
chla_temp_boots_df = readRDS("./data/derived-data/chla_temp_boots.rds")
chla_temp_coefs_df = chla_temp_boots_df %>% # summarised
  summarise(tempkt_est = na.rm_mean(tempkt_stand, na.rm = TRUE),
            tempkt_quant2.5 = quantile(tempkt_stand, 0.025, na.rm = TRUE),
            tempkt_quant25 = quantile(tempkt_stand, 0.25,na.rm = TRUE),
            tempkt_quant75 = quantile(tempkt_stand, 0.75,na.rm = TRUE),
            tempkt_quant97.5 = quantile(tempkt_stand,0.975, na.rm =TRUE),
            tempkt_min = min(tempkt_stand, na.rm = TRUE),
            temptk_max = max(tempkt_stand, na.rm = TRUE),
            adj_r2_mean = na.rm_mean(adj_r2, na.rm = TRUE),
            intercept_est = na.rm_mean(Intercept, na.rm = TRUE))
# rm(chla_temp_boots_df)

## 
#### annual evenness of population production objects ####
## bootstrap within stream temperature dependences from interval bootstraps ##
#load interval boot objects
stream_prod_list = read.table( "./data/derived-data/stream_prod_df.txt", header = TRUE, sep = "\t") %>%  #full interval production bootstrap list
  rename(DATE = "date_id") %>% pivot_wider( names_from = 'boot_id', values_from = 'prod_mg_m_int') %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list))
stream_prod_d_list = read.table( "./data/derived-data/stream_prod_d_df.txt", header = TRUE, sep = "\t") %>%  #full interval daily production bootstap list
  rename(DATE = "date_id") %>% pivot_wider( names_from = 'boot_id', values_from = 'prod_mg_m_d') %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list)) 

## stream-taxa-date mean body size list. Large file ~68.1 MB when loaded ##
## this was used to estimate the biomass-weighted body size corrections in B and PB ##
# stream_M_int = readRDS("./data/derived-data/stream_M_fin.rds") # full interval body size list
# stream_bio_list = read.table( "./data/derived-data/stream_bio_df.txt", header = TRUE, sep = "\t") %>%  # full interval raw biomass list
#   rename(DATE = "date_id") %>% pivot_wider( names_from = 'boot_id', values_from = 'bio_mg_m') %>%
#   junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list)) 
# stream_pb_list = read.table( "./data/derived-data/stream_pb_df.txt", header = TRUE, sep = "\t") %>%  # full interval raw P:B list
#   rename(DATE = "date_id") %>% pivot_wider( names_from = 'boot_id', values_from = 'pb_int') %>%
#   junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list)) 
# stream_N_list = read.table( "./data/derived-data/stream_N_df.txt", header = TRUE, sep = "\t") %>%  # full interval N list
#   rename(DATE = "date_id") %>% pivot_wider( names_from = 'boot_id', values_from = 'n_ind_m') %>%
#   junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list))
# 
# ## the above files are transition files to get to corrected files below ##

stream_biocorr_list = read.table( "./data/derived-data/stream_biocorr_df.txt", header = TRUE, sep = "\t") %>%  # full interval mass-corrected biomass list
  rename(DATE = "date_id") %>% pivot_wider( names_from = 'boot_id', values_from = 'corrbio_mg_m') %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list)) 
stream_pbcorr_list = read.table( "./data/derived-data/stream_pbcorr_df.txt", header = TRUE, sep = "\t") %>%  # full interval mass-corrected P:B list
  rename(DATE = "date_id") %>% pivot_wider( names_from = 'boot_id', values_from = 'corrpb_int') %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list)) 

## load the bootstrapped apparent temp dependences of B, PB, and Prod ##
bio_stream_coef_lists = read.table( "./data/derived-data/bio_stream_coefs_df.txt", header = TRUE, sep = "\t") %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list))
pb_stream_coef_lists = read.table( "./data/derived-data/pb_stream_coefs_df.txt", header = TRUE, sep = "\t") %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list))
prod_stream_coef_lists = read.table( "./data/derived-data/prod_stream_coefs_df.txt", header = TRUE, sep = "\t") %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list))
prod_d_stream_coef_lists= read.table( "./data/derived-data/prod_d_stream_coefs_df.txt", header = TRUE, sep = "\t") %>%
  junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list))
 

## summarise the estimated temperature dependences across streams ##
### Bootstrapped estimates of within-stream temperature dependence of P ###
prod_d_stream_coef_bounds =map(prod_d_stream_coef_lists, ~.x %>% summarise(prod_slope_est = na.rm_mean(tempkt_stand, na.rm = TRUE),
                                                                           prod_int_est = na.rm_mean(Intercept, na.rm = TRUE),
                                                                           prod_quant_25 = quantile(tempkt_stand, 0.25, na.rm = TRUE),
                                                                           prod_quant_75 = quantile(tempkt_stand, 0.75, na.rm = TRUE),
                                                                           prod_quant_50 = quantile(tempkt_stand, 0.5, na.rm = TRUE),
                                                                           prod_quant_2.5 = quantile(tempkt_stand, 0.025, na.rm = TRUE),
                                                                           prod_quant_97 = quantile(tempkt_stand, 0.975, na.rm = TRUE),
                                                                           prod_max_est = max(tempkt_stand, na.rm = TRUE),
                                                                           prod_min_est = min(tempkt_stand, na.rm = TRUE),
                                                                           prod_r2_avg = na.rm_mean(adj_r2, na.rm = TRUE),
                                                                           prod_r2_med = median(adj_r2, na.rm = TRUE)))
#create prod_d df
prod_d_stream_coef_df = list()
for(i in 1:length(prod_d_stream_coef_bounds)){
  site_id = names(prod_d_stream_coef_bounds)[i]
  prod_d_stream_coef_df[[i]] = data.frame(site_id = site_id, prod_d_stream_coef_bounds[[i]])
}

prod_d_stream_coef_df = prod_d_stream_coef_df %>% bind_rows() %>%
  mutate(site_id = factor(site_id, levels = rev(stream_order)))

### Bootstrapped estimates of within-stream temperature dependence of B ###

bio_stream_coef_bounds =map(bio_stream_coef_lists, ~.x %>% summarise(bio_slope_est = na.rm_mean(tempkt_stand, na.rm = TRUE),
                                                                     bio_int_est = na.rm_mean(Intercept, na.rm = TRUE),
                                                                     bio_quant_50 = quantile(tempkt_stand, 0.5, na.rm = TRUE),
                                                                     bio_quant_25 = quantile(tempkt_stand, 0.25, na.rm = TRUE),
                                                                     bio_quant_75 = quantile(tempkt_stand, 0.75, na.rm = TRUE),
                                                                     bio_quant_2.5 = quantile(tempkt_stand, 0.025, na.rm = TRUE),
                                                                     bio_quant_97 = quantile(tempkt_stand, 0.975, na.rm = TRUE),
                                                                     bio_max_est = max(tempkt_stand, na.rm = TRUE),
                                                                     bio_min_est = min(tempkt_stand, na.rm = TRUE),
                                                                     bio_r2_avg = na.rm_mean(adj_r2, na.rm = TRUE),
                                                                     bio_r2_med = median(adj_r2, na.rm = TRUE)))
#create bio df
bio_stream_coef_df = list()
for(i in 1:length(bio_stream_coef_bounds)){
  site_id = names(bio_stream_coef_bounds)[i]
  bio_stream_coef_df[[i]] = data.frame(site_id = site_id, bio_stream_coef_bounds[[i]])
}
bio_stream_coef_df = bio_stream_coef_df %>%
  bind_rows() %>% mutate(site_id = factor(site_id, levels = rev(stream_order)))

pb_stream_coef_bounds =map(pb_stream_coef_lists, ~.x %>% summarise(pb_slope_est = na.rm_mean(tempkt_stand, na.rm = TRUE),
                                                                   pb_int_est = na.rm_mean(Intercept, na.rm = TRUE),
                                                                   pb_quant_50 = quantile(tempkt_stand, 0.5, na.rm = TRUE),
                                                                   pb_quant_25 = quantile(tempkt_stand, 0.25, na.rm = TRUE),
                                                                   pb_quant_75 = quantile(tempkt_stand, 0.75, na.rm = TRUE),
                                                                   pb_quant_2.5 = quantile(tempkt_stand, 0.025, na.rm = TRUE),
                                                                   pb_quant_97 = quantile(tempkt_stand, 0.975, na.rm = TRUE),
                                                                   pb_max_est = max(tempkt_stand, na.rm = TRUE),
                                                                   pb_min_est = min(tempkt_stand, na.rm = TRUE),
                                                                   pb_r2_avg = na.rm_mean(adj_r2, na.rm = TRUE),
                                                                   pb_r2_med = median(adj_r2, na.rm = TRUE)))

#create pb df
pb_stream_coef_df = list()
for(i in 1:length(pb_stream_coef_bounds)){
  site_id = names(pb_stream_coef_bounds)[i]
  pb_stream_coef_df[[i]] = data.frame(site_id = site_id, pb_stream_coef_bounds[[i]])
}
pb_stream_coef_df = pb_stream_coef_df %>%
  bind_rows() %>% mutate(site_id = factor(site_id, levels = rev(stream_order)))

###create annual_df and annual_summ files based on bootstrapped files
prod_ann_df = map(stream_prod_list, ~.x %>% 
                        select(site_id, contains("X")) %>% 
                        group_by(site_id, .drop = FALSE) %>%
                        summarise_at(vars(contains("X")), sum)) %>%
  bind_rows() %>% pivot_longer(cols = contains("X"), names_to = "boot_id", values_to = "prod_y") %>%
  arrange(factor(site_id, rev(stream_order)))
corrbio_ann_df = map(stream_biocorr_list, ~.x %>%
                           select(site_id, contains("X")) %>%
                           group_by(site_id, .drop = FALSE) %>%
                           summarise_at(vars(contains("X")), na.rm_mean, na.rm = TRUE)) %>%
  bind_rows() %>% pivot_longer(cols = contains("X"), names_to = "boot_id", values_to = "corrbio_mg_m") %>%
  arrange(factor(site_id, rev(stream_order)))
corrpb_ann_df = map(stream_pbcorr_list, ~.x %>%
                          select(site_id, contains("X")) %>%
                          group_by(site_id, .drop = FALSE) %>%
                          summarise_at(vars(contains("X")), na.rm_mean, na.rm = TRUE)) %>%
  bind_rows() %>% pivot_longer(cols = contains("X"), names_to = "boot_id", values_to = "corrpb_y") %>%
  arrange(factor(site_id, rev(stream_order)))

annual_df = left_join(prod_ann_df, corrbio_ann_df) %>% left_join(corrpb_ann_df) %>% left_join(chla_boot_df) %>% select(-boot_id)
rm(list = ls()[ls() %in% c('prod_ann_df', 'corrbio_ann_df', 'corrpb_ann_df', 'chla_boot_df')])
annual_df = annual_df %>% left_join(ann_stream_temps) %>% mutate(tempkt_stand = C_to_overkt_stand15(tempC)) %>%
left_join(chla_ann_summ %>% select(site_id, contains('chla')))

annual_summ = annual_df %>% group_by(site_id) %>%
  summarise(prod_mg_m_y = na.rm_mean(prod_y, na.rm = TRUE),
            prod_quant2.5 = quantile(prod_y, 0.025, na.rm = TRUE),
            prod_quant25 = quantile(prod_y, 0.25, na.rm = TRUE),
            prod_quant75 = quantile(prod_y, 0.75, na.rm = TRUE),
            prod_quant97.5 = quantile(prod_y, 0.975, na.rm = TRUE),
            prod_max = max(prod_y, na.rm = TRUE),
            prod_min = min(prod_y, na.rm = TRUE),
            bio_mg_m = na.rm_mean(corrbio_mg_m, na.rm = TRUE),
            bio_quant2.5 = quantile(corrbio_mg_m, 0.025, na.rm = TRUE),
            bio_quant25 = quantile(corrbio_mg_m, 0.25, na.rm = TRUE),
            bio_quant75 = quantile(corrbio_mg_m, 0.75, na.rm = TRUE),
            bio_quant97.5 = quantile(corrbio_mg_m, 0.975, na.rm = TRUE),
            bio_max = max(corrbio_mg_m, na.rm = TRUE),
            bio_min = min(corrbio_mg_m, na.rm = TRUE),
            pb = na.rm_mean(corrpb_y, na.rm = TRUE),
            pb_quant2.5 = quantile(corrpb_y, 0.025, na.rm = TRUE),
            pb_quant25 = quantile(corrpb_y, 0.25, na.rm = TRUE), 
            pb_quant75 = quantile(corrpb_y, 0.75, na.rm = TRUE), 
            pb_quant97.5 = quantile(corrpb_y, 0.975, na.rm = TRUE),
            pb_max = max(corrpb_y, na.rm = TRUE),
            pb_min = min(corrpb_y, na.rm = TRUE),
            tempkt = unique(tempkt),
            tempC = unique(tempC),
            tempkt_stand = unique(tempkt_stand)) %>%
  ungroup() %>% 
  mutate(site_id = factor(site_id, levels = stream_order),
         prod_adj = prod_mg_m_y/365,
         prod_quant25_adj = prod_quant25/365,
         prod_quant75_adj = prod_quant75/365,
         prod_quant2.5_adj = prod_quant2.5/365,
         prod_quant97.5_adj = prod_quant97.5/365,
         prod_max_adj = prod_max/365,
         prod_min_adj = prod_min/365) %>%
  left_join(chla_summ_ann %>% select(site_id, contains('chla'))) %>%
  mutate(site_id = factor(site_id, levels = stream_order))


#### bootstrapping annual temperature dependence of secondary production, biomass, & turnover ####
#standardized temperature @ 15C
# set.seed(123);tic();prod_boot_df = tempkt_lmCI_boot(annual_df, 'prod_mg_m_y', nboot = 10000, standardize = TRUE);toc()#slope, intercept, and r2  = ~16secs
# set.seed(123);tic();prod_boot_C = tempkt_lmCI_boot(annual_df_update, 'prod_mg_m_y', nboot = 10000, standardize = FALSE);toc()# slope, intercept and r2 = ~17sec
# set.seed(123);bio_boot_df = tempkt_lmCI_boot(annual_df_update, 'corrbio_mg_m_new', nboot = 10000, standardize = TRUE)
# set.seed(123);bio_boot_C = tempkt_lmCI_boot(annual_df_update, 'corrbio_mg_m_new', nboot = 10000, standardize = FALSE)
# set.seed(123);pb_boot_df = tempkt_lmCI_boot(annual_df_update, 'corrpb_new', nboot = 10000, standardize = TRUE)
# set.seed(123);pb_boot_C = tempkt_lmCI_boot(annual_df_update, 'corrpb_new', nboot = 10000, standardize = FALSE)
# saveRDS(prod_boot_df, "./data/derived-data/prod_boot_df.rds")
# saveRDS(prod_boot_C, "./data/derived-data/prod_boot_C.rds")
# saveRDS(bio_boot_df, "./data/derived-data/bio_boot_df.rds")
# saveRDS(bio_boot_C, "./data/derived-data/bio_boot_C.rds")
# saveRDS(pb_boot_df, "./data/derived-data/pb_boot_df.rds")
# saveRDS(pb_boot_C, "./data/derived-data/pb_boot_C.rds")

prod_boot_df = readRDS(file = "./data/derived-data/prod_boot_df.rds") 
prod_boot_C = readRDS(file=  "./data/derived-data/prod_boot_C.rds")
bio_boot_df = readRDS(file = "./data/derived-data/bio_boot_df.rds")
pb_boot_df = readRDS(file = "./data/derived-data/pb_boot_df.rds")
bio_boot_C = readRDS(file=  "./data/derived-data/bio_boot_C.rds")
pb_boot_C = readRDS(file=  "./data/derived-data/pb_boot_C.rds")

## the bootstrapped relationship between chla and production
chla_prod_lmboots = read.table("./data/derived-data/chla_prod_lmboots.txt",sep = "\t", header = TRUE)

## Create predicted dataframe of mean relationship for plotting ##
chla_prod_predict = data.frame(chla = seq(min(annual_summ$ann_chla_mean), max(annual_summ$ann_chla_mean), length.out = 100)) %>%
  mutate(prod= exp(mean(chla_prod_lmboots[,'Intercept'])+mean(chla_prod_lmboots[,'chla'])*log(chla)),
         prod_adj = prod/365)

## bootstrapped resource adjusted temp dependence of prod across streams ##
chlatemp_prod_lmboots = read.table("./data/derived-data/chlatemp_prod_lmboots.txt",sep= "\t", header = TRUE)

chla_tempdep_coef = chlatemp_prod_lmboots %>%
  summarise(tempkt_est = na.rm_mean(tempkt_stand, na.rm = TRUE),
            tempkt_quant2.5 = quantile(tempkt_stand, 0.025, na.rm = TRUE),
            tempkt_quant25 = quantile(tempkt_stand, 0.25,na.rm = TRUE),
            tempkt_quant75 = quantile(tempkt_stand, 0.75,na.rm = TRUE),
            tempkt_quant97.5 = quantile(tempkt_stand,0.975, na.rm =TRUE),
            tempkt_min = min(tempkt_stand, na.rm = TRUE),
            temptk_max = max(tempkt_stand, na.rm = TRUE),
            adj_r2_mean = na.rm_mean(adj_r2, na.rm = TRUE),
            intercept_est = na.rm_mean(Intercept, na.rm = TRUE))

## files for working with seasonal within stream production
int_df = read.table("./data/derived-data/int_df.txt", sep = "\t", header = TRUE) %>%
  mutate(tempC= overkt_to_C(tempkt), 
         tempkt_stand = C_to_overkt_stand15(tempC))
int_df = int_df %>% select(site_id, DATE, JULIAN, doy, days, tempC, tempkt, tempkt_stand,int_adj)
# isolate mean light intensity by sampling dates across streams
int_df = bind_cols(int_df, light_subset(int_df, tempkt_light.doy) %>% data.frame %>% dplyr::slice(1:(n()-1)) %>%
  setNames(.,c("mean_light","sum_light"))) %>% mutate(DATE = as.Date(DATE), daily_light = sum_light/days)

## add in chla data for stream-intervals ##
## match some date offsets from sampling ##
chla_mod = chla %>% select(-JULIAN) %>% 
  mutate(doy = ifelse(site_id == "st14" & DATE == as.Date("2012-06-30"), doy-1, doy),
         doy = ifelse(site_id == "st9" & DATE == "2012-06-30", doy - 11,doy),
         doy = ifelse(site_id == "hver" & DATE == "2012-06-30", doy - 12,doy),
         doy = ifelse(site_id == "st6" & DATE == "2012-06-30", doy - 11,doy),
         DATE = if_else(site_id == "st14" & DATE == as.Date("2012-06-30"), as.Date("2012-06-29"), as.Date(DATE)),
         DATE = if_else(site_id == "st9" & DATE == "2012-06-30", as.Date("2012-06-19"), as.Date(DATE)),
         DATE = if_else(site_id == "hver" & DATE == "2011-07-01", as.Date("2012-06-19"), as.Date(DATE)),
         DATE = if_else(site_id == "st6" & DATE == "2012-06-30", as.Date("2012-06-19"), as.Date(DATE)))

int_df = int_df %>% left_join(chla_mod %>% 
  group_by(site_id, DATE) %>%
  summarise(chla_mean = na.rm_mean(chla_mg_m, na.rm = TRUE))) %>%
    ungroup %>% mutate(site_id = factor(site_id, levels = stream_order))


int_df = int_df %>% left_join(chla_ann_summ %>% select(site_id, chla_mean_ann))
int_df = int_df %>% 
  left_join(ann_stream_temps %>% select(site_id, mean_tempkt = 'tempkt')) %>%
  mutate(tempkt_center = mean_tempkt - tempkt)

#prediction to plot within stream estimates#
prod_within_pred = do.call(rbind, 
                           lapply(split(int_df, int_df$site_id),
                                  function(i) with(i, expand.grid(site_id = unique(site_id), 
                                                                  tempkt_stand = seq(max(tempkt_stand), min(tempkt_stand), 
                                                                                     length = 100)))))
## from bootstrapped estimates ##
prod_within_pred = prod_within_pred %>% left_join(prod_d_stream_coef_df) %>%
  left_join(bio_stream_coef_df) %>% left_join(pb_stream_coef_df) %>%
  mutate(prod_d_new = prod_int_est + prod_slope_est*tempkt_stand,
         bio_new = bio_int_est + bio_slope_est*tempkt_stand,
         pb_new = pb_int_est + pb_slope_est*tempkt_stand) %>%
  select(c(site_id:tempkt_stand, prod_d_new:pb_new)) %>% mutate(site_id = factor(site_id, levels = stream_order))


prod_d_df = map(stream_prod_d_list, ~.x %>% select(site_id, DATE, contains("X")) %>%
                      pivot_longer(cols = contains("X"), values_to = "prod_d", names_to = "boot") %>% select(-boot) %>% 
                      group_by(site_id, DATE) %>% summarise(prod_d = na.rm_mean(prod_d, na.rm = TRUE))) %>%
  bind_rows() %>% ungroup() %>% mutate(site_id = factor(site_id, levels = stream_order), DATE = as.Date(DATE))
bio_mg_m_df = map(stream_biocorr_list, ~.x %>% select(site_id, DATE, contains('X')) %>%
                    pivot_longer(cols = contains('X'), values_to = 'bio_mg_m', names_to = 'boot') %>% select(-boot) %>%
                    group_by(site_id, DATE) %>% summarise(corrbio_mg_m = na.rm_mean(bio_mg_m, na.rm = TRUE))) %>%
  bind_rows() %>% ungroup() %>% mutate(site_id = factor(site_id, levels = stream_order), DATE = as.Date(DATE))
pb_df = map(stream_pbcorr_list, ~.x %>% select(site_id, DATE, contains('X')) %>%
              pivot_longer(cols = contains('X'), values_to = 'pb', names_to = 'boot') %>% select(-boot) %>%
              group_by(site_id, DATE) %>% summarise(corrpb_d = na.rm_mean(pb, na.rm = TRUE))) %>%
  bind_rows() %>% ungroup() %>% mutate(site_id = factor(site_id, levels = stream_order), DATE = as.Date(DATE))
int_df = int_df %>% left_join(prod_d_df) %>% left_join(bio_mg_m_df) %>% left_join(pb_df) %>% mutate(site_id = factor(site_id, levels = stream_order))

#######
#### bootstrap model fitting for mixed models ####
#### full model fitting/selection takes ~10min to run ####
#### and produces ~ 31 MB object ####
# tic();set.seed(123);mod_sel_list = suppressMessages(suppressWarnings(mod_sel_boot(stream_prod_d_df, chla_int_df %>% mutate(DATE = as.Date(DATE)), int_df, nboot = 1000)));toc()
mod_sel_list <<- readRDS(file = "./data/derived-data/mod_sel_list.rds")
## top model is full additive w/ 1|SITE ##
#### fit models with REML on bootstrapped data ####
### bootstrap model fitting of mixed models ###

temp_var = "tempkt_stand"
nboot = 1000

## model refit with REML on standardized boltzmann temp ##
top_mod_centerfit <<- readRDS("./data/derived-data/top_mod_centerfit.rds")

##### primary production w/ variability #####
prim_prod <- read.csv("./data/raw-data/3fA_ModeledDailyGPP.csv", header = TRUE) %>%
  select(-1) %>% mutate(stream = tolower(stream))
prim_prod =  prim_prod %>% select(Pdt,DailyIntPAR_umolM2d,EstGPP_gCm2d_med, 
                                  EstGPP_gCm2d_U95per, EstGPP_gCm2d_L95per, stream) %>%
  mutate(Pdt = as.Date(Pdt)) %>%
  filter(Pdt < "2011-10-26" & Pdt >= "2010-10-26") %>% 
  left_join(int_df %>% filter(site_id %in% c("st7","oh2")) %>% 
              rename(Pdt = DATE, stream = site_id) %>%
              select(Pdt, stream, chla_mean) %>%
              mutate(Pdt = as.Date(Pdt))) %>%
  inner_join(tempkt_full %>% select(DATE, contains(c("st7", "oh2"))) %>%
               rename_at(vars(contains("tempkt")), ~ str_remove(.,"_tempkt")) %>%
               rename(Pdt = DATE) %>% mutate(Pdt = as.Date(Pdt)) %>% 
               pivot_longer(-Pdt, names_to = "stream", values_to = "tempkt")) %>%
  mutate(tempC = overkt_to_C(tempkt),
         chla_approx = chla_mean) 

prim_prod_list = prim_prod %>% junkR::named_group_split(stream) %>%
  map(.,function(x){
    if(is.na(x[dim(x)[1], "chla_approx"])){
      x[dim(x)[1],"chla_approx"] = x[1, "chla_mean"]
      x <- x %>% mutate(chla_approx = zoo::na.approx(chla_approx))
      return(x)} else {
        x <- x %>% mutate(chla_approx = zoo::na.approx(chla_approx))
        return(x)}
  })

prim_prod_boot <<- readRDS("./data/derived-data/prim_prod_boot.rds")
#create seasonal warming data frame
sites_list = int_df %>% filter(site_id %in% c("st7","oh2")) %>% junkR::named_group_split(site_id) %>% rlist::list.subset(names(stream_order_list))
sites_list = sites_list[names(sites_list) %in% c("st7","oh2")]
gpp_list = prim_prod_boot %>% map(~.x %>% rename(DATE = "Pdt") %>% mutate(DATE = as.Date(DATE)))
prod_df = stream_prod_d_list[names(stream_prod_d_list) %in% c("st7","oh2")] %>% bind_rows()


var = names(gpp_list[[1]])[grepl("EstGPP_",names(gpp_list[[1]]))]
warming_gpp = mapply(gpp_subset, sites_list, gpp_list, MoreArgs = list(var = var), SIMPLIFY = FALSE)
warming_gpp[[1]]$site_id = names(warming_gpp[1])
warming_gpp[[2]]$site_id = names(warming_gpp[2])

warming_gpp = bind_rows(warming_gpp)

warming_intdf = bind_rows(sites_list) %>% left_join(warming_gpp) %>%
  mutate(gpp_d = EstGPP_gCm2d_med/days)
warming_gpp_chla = lm(log(gpp_d)~log(mean_light)+tempkt_stand+log(chla_mean), data = warming_intdf);summary(warming_gpp_chla)

warming_df <<- warming_intdf %>%
  mutate(gpp_d_chla_warming = exp(predict(warming_gpp_chla, newdata = warming_intdf)),
         EstGPP_gCm2d_warming = gpp_d_chla_warming,
         EstGPP_gCm2_warming = gpp_d_chla_warming*days)

warming_intdf = warming_intdf %>% 
  select(-EstGPP_gCm2d_L95per, -EstGPP_gCm2d_U95per, -EstGPP_gCm2d_med, -gpp_d) %>%
  mutate_at(vars(contains('EstGPP_')), ~(./days))

warming_intdf <<- warming_intdf %>% 
  left_join(ann_stream_temps %>% rename(mean_tempkt = 'tempkt') %>% select(site_id, mean_tempkt)) %>%
  mutate(tempkt_center = mean_tempkt - tempkt) %>% select(-mean_tempkt)

#create matrix of estimated EaGPP
gpp_boot_cols = names(warming_intdf[grepl("EstGPP_",names(warming_intdf))])
gpp_boot_list = vector('list', nboot)
gpp_Ea_mat = matrix(NA, ncol = 3, nrow = length(gpp_boot_cols))
for(i in 1:nrow(gpp_Ea_mat)){
  gpp = warming_intdf %>% select(site_id,days, tempkt_center, tempkt_stand,chla_mean, mean_light,!!as.symbol(gpp_boot_cols[i])) %>% mutate(gpp_d = !!as.symbol(gpp_boot_cols[i]))
  gen_lm = suppressMessages(suppressWarnings(lme4::lmer(log(gpp_d)~log(chla_mean)+tempkt_stand+log(mean_light)+(1|site_id), data = gpp)))
 gpp_Ea_mat[i,] = c(summary(gen_lm)[['coefficients']]['tempkt_stand', 'Estimate'], MuMIn::r.squaredGLMM(gen_lm))
  gpp_boot_list[[i]] = gen_lm
}

gpp_temp_vec = sapply(gpp_boot_list, function(x) summary(x)[["coefficients"]][temp_var,"Estimate"])
prod_Ea_vec = sapply(top_mod_centerfit[[2]], function(x) summary(x)[["coefficients"]][temp_var,"Estimate"])
# prod_int_vec = sapply(top_mod_centerfit[[2]], function(x) summary(x)[["coefficients"]]["(Intercept)","Estimate"])
# mean_int = mean(prod_int_vec)
# adj_int = prod_int_vec-mean_int

gpp_adj_prod = prod_Ea_vec - gpp_Ea_mat[1]
# quantile(prod_Ea_vec, probs = c(0.025,0.5,0.975))
# na.rm_mean(prod_Ea_vec)
# quantile(gpp_Ea_mat[,1], probs = c(0.025, 0.5, 0.975))
# na.rm_mean(gpp_Ea_mat[,1])
# na.rm_mean(gpp_adj_prod)
# quantile(gpp_adj_prod, probs = c(0.025,0.5,0.975))
# na.rm_mean(boot_mat[,2])

seas_EA_df = data.frame(time = 'Within-streams', apparent = prod_Ea_vec, adjusted = gpp_adj_prod)
ann_Ea_df = data.frame(time = 'Among-streams', apparent = prod_boot_df$tempkt_stand, adjusted = chlatemp_prod_lmboots$tempkt_stand)

full_Ea_df = bind_rows(seas_EA_df, ann_Ea_df)

full_Ea_df_long = full_Ea_df %>%
  pivot_longer(apparent:adjusted, names_to = "type", values_to = "Ea") %>%
  mutate(type = factor(type, levels = c("apparent", "adjusted")),
         time = factor(time, levels = c("Among-streams", "Within-streams")))
##plotting attributes###
load("./data/ocecolors.rda")#color scheme object
oce_temp_disc = c("#E5FA6A","#CF696C","#544685","#072C44","#082A40","#0B222E")#color codes
oce_temp_pos = c(256,212,168,124,80,1)#color positions in 'temperature' list of ocecolors
breaks_in_C <- scales::extended_breaks()(overkt_to_C(c(42,41,40,39,38)))#> [1]  0 10 20 30 
breaks_in_C_plus = c(0,5,10,15,20,25,30,35, 40)
breaks_in_kt <- round(C_to_overkt(breaks_in_C_plus),1)
# breaks_in_kt
# [1] 42.5 41.7 41.0 40.3 39.6 38.9 38.3 37.7
breaks_in_kt_stand = round(C_to_overkt_stand15(breaks_in_C_plus),1)
#breaks_in_kt_stand
#[1] -2.2 -1.4 -0.7  0.0  0.7  1.4  2.0  2.6
breaks_in_logchla = c(log(0.1),log(0.5),log(1),log(5),log(10),log(50),log(100))
breaks_in_chla = exp(breaks_in_logchla)

breaks_in_logprod = c(log(0.001),log(0.005),log(0.01),log(0.05),log(1),log(5),log(10),log(50),log(100),log(500),log(1000),log(5000),log(10000),log(50000))
breaks_in_prod = exp(breaks_in_logprod)

breaks_in_log10prod <- c(log(0.001),log(0.01),log(0.1),log(1),log(10),log(100),log(1000),log(10000),log(100000))
#> [1]  -5.0 -2.5 0 2.5 5.0 7.5
breaks_in_prod10 <- exp(breaks_in_log10prod)

rm(list = ls()[ls() %in% c()])