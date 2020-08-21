load_packages = function(){
  if(!require("pacman")) install.packages("pacman")
  library(pacman)
  package.list <- c("twitteR","data.table", "RCurl","plyr","tidyverse","furrr","rriskDistributions",
                    "tictoc","chron","lubridate","httr","TTR", "lmerTest", "rlist",
                    "grid","gridExtra", "ggridges", "lme4", "nlme", "sjPlot", "VGAM",
                    "viridis", "broom","broom.mixed", "bbmle","ggthemes", "ggeffects")
  p_load(char = package.list, install = T)
  rm("package.list")
  ins_julian_path = getURL("https://raw.githubusercontent.com/jimjunker1/secprod_workflow/master/len_freq/ins_julian_function.txt",ssl.verifypeer = FALSE)
  ins_julian <<- eval(parse(text = ins_julian_path))
  pacman::p_load_gh(c("jimjunker1/junkR"))
  #  len_freq_path = getURL("https://raw.githubusercontent.com/jimjunker1/secprod_workflow/master/len_freq/len_freq_function.R",ssl.verifypeer = FALSE)
  #  len_freq <<- eval(parse(text = len_freq_path))
  cbbPalette <<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  
  C_to_overkt <- function(a){1/(8.61733*10^-5*(a+273.15))}#overkt function
  overkt_to_C <- function(a){1/(a*(8.61733*10^-5)) - 273.15}
  C_to_overkt_stand15 <- function(a){(1/(8.61733e-5*(15+273.15)) - (1/(8.61733e-5*(a+273.15))))}
  overkt_stand15_to_C <<- function(a){1/((C_to_overkt(15)-a)*(8.61733*10^-5)) - 273.15}
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  na.rm_mean <- function(...,na.rm=TRUE){mean(c(...),na.rm=na.rm)}
  options(stringsAsFactors = F)
  options(max.print = 1000000)
  
  mass_corrPB <<- function(a,b){a*(b^-0.25)}
  mass_corrBIO <<- function(a,b){a*(b^0.25)}
  
  temp_corrPB <<- function(a,b){a*(b^0.65)}
  temp_corrBIO <<- function(a,b){a*(b^-0.65)}
  
  '%ni%' <<- Negate('%in%')
  
  myspread <<- function(df, key, value) {
    # quote key
    keyq <- rlang::enquo(key)
    # break value vector into quotes
    valueq <- rlang::enquo(value)
    s <- rlang::quos(!!valueq)
    df %>% gather(variable, value, !!!s) %>%
      unite(temp, !!keyq, variable) %>%
      spread(temp, value)
  }
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  temp_subset <<- function(bugs_df, TEMP, wrap = TRUE){#generic function for 
    #subsetting temperatures between sampling dates and averaging
    na.rm_mean <- function(...,na.rm=TRUE){mean(c(...),na.rm=na.rm)}
    dates = data.frame(bugs_df %>% group_by(DATE) %>% summarise(JULIAN = unique(JULIAN)))
    dates <- dates %>% arrange(JULIAN)
    int_tempkt = matrix(NA, ncol = 1, nrow = nrow(dates)+1)
    if (wrap == TRUE){
      last.date = which(dates$JULIAN == max(dates$JULIAN))
      first.date = which(dates$JULIAN == min(dates$JULIAN))
      j.wrap = dates[last.date,'JULIAN'] + (365-(dates$JULIAN[last.date]-dates$JULIAN[first.date]))
      d.julian = unlist(month.day.year(j.wrap, origin. = c(month = 01, day = 01, year = 2010)))
      d.wrap = as.POSIXct(paste0(d.julian['year'],"-",d.julian['day'],"-",d.julian['month']), format = "%Y-%d-%m", tz = "UTC")
      wraps = data.frame(DATE = as.Date(d.wrap), JULIAN = j.wrap)
      dates = bind_rows(dates, wraps)
      dates <- dates %>% mutate(DATE = as.POSIXct(DATE, format = "%Y-%m-%d", tz = "UTC"))
      dates$tempkt = NA
      julians = c(dates$JULIAN ,j.wrap)
    }
    for(i in seq_along(julians)){
      temp_sub = subset(TEMP, JULIAN >= julians[i] & JULIAN <= julians[i+1])
      dates[i,'tempkt'] = temp_sub %>% summarise_at(vars(contains('temp')), funs(na.rm_mean))
    }
    dates = dates[-dim(dates)[1],]
    return(data.frame(dates))
  }
  tempkt_lmCI_boot<<- function(x, measure, nboot = 500, standardize = FALSE){
    x_sub = x %>% 
      select(site_id, measure) %>%
      group_by(site_id) %>% mutate(id = 1:n()) %>%
      ungroup() %>%
      spread(key = site_id, value = measure) %>%
      select(-id)
    
    boot_df = matrix(ncol = dim(x_sub)[2], nrow = nboot)
    for(col in 1:dim(x_sub)[2]){
      boot_df[,col] = sample(t(x_sub[,col]), nboot, replace = TRUE)
    }
    if(standardize == FALSE){
      temp_sub = x %>% group_by(site_id) %>%
        summarize(tempC = unique(tempC))
      boot_df =  boot_df %>% data.frame() %>%
        setNames(., names(x_sub)) %>%
        t() %>% data.frame() %>%
        rownames_to_column("site_id") %>%
        left_join(temp_sub)
      
      slope_df = matrix(ncol = 3, nrow = nboot)
      for(col in 2:(dim(boot_df)[2]-1)){
        x_ann = boot_df[,col]
        tempC = boot_df[,'tempC']
        temp.lm = lm(log(x_ann)~tempC)
        slope_df[(col-1),1] = coef(temp.lm)['(Intercept)']
        slope_df[(col-1),2] = coef(temp.lm)['tempC']
        slope_df[(col-1),3] = summary(temp.lm)[['r.squared']]
        
      }
      slope_df <- slope_df %>%
        data.frame() %>% setNames(.,c("Intercept", "tempC", "adj_r2"))
      return(slope_df)
    } else{
      temp_sub = x %>% group_by(site_id) %>%
        summarize(tempkt_stand = unique(tempkt_stand))
      
      boot_df =  boot_df %>% data.frame() %>%
        setNames(., names(x_sub)) %>%
        t() %>% data.frame() %>%
        rownames_to_column("site_id") %>%
        left_join(temp_sub)
      
      slope_df = matrix(ncol = 3, nrow = nboot)
      for(col in 2:(dim(boot_df)[2]-1)){
        x_ann = boot_df[,col]
        tempkt_stand = boot_df[,'tempkt_stand']
        temp.lm = lm(log(x_ann)~tempkt_stand)
        slope_df[(col-1),1] = coef(temp.lm)['(Intercept)']
        slope_df[(col-1),2] = coef(temp.lm)['tempkt_stand']
        slope_df[(col-1),3] = summary(temp.lm)[['r.squared']]
      }
      slope_df <- slope_df %>%
        data.frame() %>% setNames(.,c("Intercept", "tempkt_stand", "adj_r2"))
      return(slope_df)
    }
  }
  myunlist <<- function(inlist) {
    x <- as.data.frame(do.call(rbind, inlist))
    x[] <- lapply(x, unlist)
    x
  }
  int_tempkt_lmCI_boot <<- function(x, nboot = 500, standardize = FALSE){
    stream_name = unique(x$site_id)
    tempkt = x$tempkt
    int_df = matrix(NA, nrow = dim(x)[1], ncol = nboot)
    for(row in 1:dim(x)[1]){
      for(boot in 1:nboot){
        int_df[row,boot] = sample(unlist(x[row,6:dim(x)[2]]), 1, replace = TRUE)
      }
    }
    int_df = cbind(tempkt,int_df)
    int_coef_df = matrix(ncol = 3, nrow = nboot)
    if(standardize == FALSE){
      for(col in 2:(dim(int_df)[2]-1)){
        x_int = int_df[,col]
        tempkt = int_df[,'tempkt']
        int_temp.lm = lm(log(x_int)~tempkt)
        int_coef_df[(col-1),1] = coef(int_temp.lm)['(Intercept)']
        int_coef_df[(col-1),2] = coef(int_temp.lm)['tempkt']
        int_coef_df[(col-1),3] = summary(int_temp.lm)[['adj.r.squared']]
      }
      int_coef_df = int_coef_df %>% data.frame() %>% setNames(., c('Intercept', 'tempkt_stand', 'adj_r2'))
      # int_coef_df = list(int_coef_df);names(int_coef_df) = stream_name
      return(int_coef_df)
    } else if(standardize == TRUE){
      for(col in 2:(dim(int_df)[2]-1)){
        x_int = int_df[,col]
        tempkt_stand = C_to_overkt_stand15(overkt_to_C(int_df[,'tempkt']))
        int_temp.lm = lm(log(x_int)~tempkt_stand)
        int_coef_df[(col-1),1] = coef(int_temp.lm)['(Intercept)']
        int_coef_df[(col-1),2] = coef(int_temp.lm)['tempkt_stand']
        int_coef_df[(col-1),3] = summary(int_temp.lm)[['r.squared']]
      }
      int_coef_df = int_coef_df %>% data.frame() %>% setNames(., c('Intercept', 'tempkt_stand', 'adj_r2'))
      # int_coef_df = list(int_coef_df)
      # names(int_coef_df) = stream_name
      return(int_coef_df)
    }}
  int_tempkt_rlmCI_boot <<- function(x, nboot = 500, standardize = FALSE){
    stream_name = unique(x$site_id)
    tempkt = x$tempkt
    int_df = matrix(NA, nrow = dim(x)[1], ncol = nboot)
    for(row in 1:dim(x)[1]){
      for(boot in 1:nboot){
        int_df[row,boot] = sample(unlist(x[row,6:dim(x)[2]]), 1, replace = TRUE)
      }
    }
    int_df = cbind(tempkt,int_df)
    int_coef_df = matrix(ncol = 3, nrow = nboot)
    if(standardize == FALSE){
      for(col in 2:(dim(int_df)[2]-1)){
        x_int = int_df[,col]
        tempkt = int_df[,'tempkt']
        int_temp.lm = MASS::rlm(log(x_int)~tempkt)
        int_coef_df[(col-1),1] = coef(int_temp.lm)['(Intercept)']
        int_coef_df[(col-1),2] = coef(int_temp.lm)['tempkt']
        int_coef_df[(col-1),3] = summary(int_temp.lm)[['adj.r.squared']]
      }
      int_coef_df = int_coef_df %>% data.frame() %>% setNames(., c('Intercept', 'tempkt_stand', 'adj_r2'))
      # int_coef_df = list(int_coef_df);names(int_coef_df) = stream_name
      return(int_coef_df)
    } else if(standardize == TRUE){
      for(col in 2:(dim(int_df)[2]-1)){
        x_int = int_df[,col]
        tempkt_stand = C_to_overkt_stand15(overkt_to_C(int_df[,'tempkt']))
        int_temp.lm = MASS::rlm(log(x_int)~tempkt_stand)
        int_coef_df[(col-1),1] = coef(int_temp.lm)['(Intercept)']
        int_coef_df[(col-1),2] = coef(int_temp.lm)['tempkt_stand']
        int_coef_df[(col-1),3] = summary(int_temp.lm)[['r.squared']]
      }
      int_coef_df = int_coef_df %>% data.frame() %>% setNames(., c('Intercept', 'tempkt_stand', 'adj_r2'))
      # int_coef_df = list(int_coef_df)
      # names(int_coef_df) = stream_name
      return(int_coef_df)
    }}
  chla_prod_boot_function <<- function(x, nboot, ...){
    prod_boots = x %>% select(site_id, contains('prod')) %>%
      group_by(site_id) %>% sample_n(nboot, replace = TRUE) %>%
      mutate(id = 1:n()) %>%
      pivot_wider(id_cols = id, names_from = site_id, values_from = prod_y) %>%
      select(-id)
    chla_boots = x %>% select(site_id, chla_ann_mean) %>%
      group_by(site_id) %>% sample_n(nboot, replace = TRUE) %>%
      mutate(id = 1:n()) %>% 
      pivot_wider(id_cols = id, names_from = site_id, values_from = chla_ann_mean) %>%
      select(-id)
    chla_prod_coef_df = matrix(ncol = 3, nrow = nboot)
    for(row in 1:(dim(chla_prod_coef_df)[1])){
      prod = unlist(prod_boots[row,])
      chla = unlist(chla_boots[row,])
      chla_prod.lm = lm(log(prod)~log(chla))
      chla_prod_coef_df[row,1] = coef(chla_prod.lm)['(Intercept)']
      chla_prod_coef_df[row,2] = coef(chla_prod.lm)['log(chla)']
      chla_prod_coef_df[row,3] = summary(chla_prod.lm)[['adj.r.squared']]
    }
    chla_prod_coef_df = chla_prod_coef_df %>% data.frame() %>% setNames(., c('Intercept', 'chla', 'adj_r2'))
    return(chla_prod_coef_df)
  }
  ## function to bootstrap the estimated temperature-dependence of annual production ##
  ## after accounting for the effects of chlorophyll a biomass ##
  chlatemp_prod_boot_function <<- function(x, nboot, ...){
    prod_boots = x %>% select(site_id, prod_mg_m_y) %>%
      group_by(site_id) %>% sample_n(nboot, replace = TRUE) %>%
      mutate(id = 1:n()) %>%
      pivot_wider(id_cols = id, names_from = site_id, values_from = prod_mg_m_y) %>%
      select(-id)
    chla_boots = x %>% select(site_id, chla_ann_mean) %>%
      group_by(site_id) %>% sample_n(nboot, replace = TRUE) %>%
      mutate(id = 1:n()) %>% 
      pivot_wider(id_cols = id, names_from = site_id, values_from = chla_ann_mean) %>%
      select(-id)
    chla_prod_coef_df = matrix(ncol = 4, nrow = nboot)
    for(row in 1:(dim(chla_prod_coef_df)[1])){
      prod = unlist(prod_boots[row,])
      chla = unlist(chla_boots[row,])
      tempkt_stand = unique(x$tempkt)
      chla_prod.lm = lm(log(prod)~log(chla)+tempkt_stand)
      chla_prod_coef_df[row,1] = coef(chla_prod.lm)['(Intercept)']
      chla_prod_coef_df[row,2] = coef(chla_prod.lm)['log(chla)']
      chla_prod_coef_df[row,3] = coef(chla_prod.lm)['tempkt_stand']
      chla_prod_coef_df[row,4] = summary(chla_prod.lm)[['adj.r.squared']]
    }
    chla_prod_coef_df = chla_prod_coef_df %>% data.frame() %>% setNames(., c('Intercept', 'chla', 'tempkt_stand', 'adj_r2'))
    return(chla_prod_coef_df)
  }
  ##### Annual model selection of chla versus temperature versus null ####
  ann_mod_sel_boot <<- function(x, measure, nboot = 1000,...){
    if(nboot > (dim(x)[1])) {print('not enough observations')} else{
      df_measure = x %>% select(site_id,tempkt_stand, one_of(measure)) %>% 
        group_by(site_id) %>% mutate(id = 1:n()) %>% 
        pivot_wider(id_cols = c(site_id, tempkt_stand), names_from = id, values_from = eval(measure)) %>%
        ungroup()
      df_chla = x %>% select(site_id, chla_ann_mean) %>%
        group_by(site_id) %>% mutate(id = 1:n()) %>%
        pivot_wider(id_cols = site_id, values_from = chla_ann_mean, names_from = id) %>%
        ungroup()
      model_sel_df = matrix(ncol = 4, nrow = nboot)
      for(col in 1:nboot){
        tempkt_stand = unlist(df_measure %>% select(tempkt_stand))
        x_grab = unlist(df_measure[,(col+2)])
        chla_grab = unlist(df_chla[,(col+1)])
        null.lm = lm(log(x_grab)~1)
        chla.lm = lm(log(x_grab)~log(chla_grab))
        temp.lm = lm(log(x_grab)~tempkt_stand)
        full.lm = lm(log(x_grab)~log(chla_grab)*tempkt_stand)
        add.lm = lm(log(x_grab)~log(chla_grab)+tempkt_stand)
        mod_an = AIC(null.lm,
                     chla.lm,
                     temp.lm,
                     full.lm,
                     add.lm)
        mod_min <- mod_an %>% rownames_to_column('Model') %>% filter(AIC == min(AIC))
        model_sel_df[col,1:3] <- unlist(mod_min)
        mod_name = mod_min %>% pull(Model)
        model_sel_df[col,4] = summary(eval(parse(text = mod_name)))[['adj.r.squared']]
      }
      model_sel_df = model_sel_df %>% data.frame() %>%
        setNames(., c("Model", "df","AIC","adj_r2"))
      return(model_sel_df)
    }}
  ##### Model selection of interval level ####
  mod_sel_boot <<- function(prod,chla,df,nboot,...){#input wide stream_prod_df and wide chla_int_df
    if(nboot > (dim(prod)[2]-5)) {print('not enough observations')} else{
      df = df %>% select(site_id, DATE, mean_light, tempkt_stand)
      model_sel_df = matrix(ncol = 7, nrow = nboot)
      model_list = vector('list', nboot)
      for(col in 1:nboot){
        prod_d_grab = prod %>% select(site_id, DATE, prod_d = (col+5))
        chla_grab = chla %>% select(site_id, DATE, chla_mean = (col+2))
        df_grab = list(df, prod_d_grab, chla_grab) %>% reduce(left_join, by = c('site_id','DATE'))
        #models
        full.model = lmer(log(prod_d)~log(chla_mean)*tempkt_stand*scale(mean_light)+(1|site_id), data = na.omit(df_grab), REML = FALSE)
        prod_null = lmer(log(prod_d)~1+(1|site_id), data = na.omit(df_grab), REML = FALSE)
        prod_full_int= lmer(log(prod_d)~log(chla_mean)*tempkt_stand*scale(mean_light)+(1|site_id), data = na.omit(df_grab), REML = FALSE)
        prod_full_add= lmer(log(prod_d)~log(chla_mean)+tempkt_stand+scale(mean_light)+(1|site_id), data = na.omit(df_grab), REML = FALSE)
        prod_full_tempxlight = lmer(log(prod_d)~log(chla_mean)+tempkt_stand*scale(mean_light)+(1|site_id), data = na.omit(df_grab), REML = FALSE)
        chla_prod.int = lmer(log(prod_d)~log(chla_mean)+(1|site_id), data = na.omit(df_grab), REML = FALSE)
        temp_prod.int = lmer(log(prod_d)~tempkt_stand + (1|site_id), data = na.omit(df_grab), REML = FALSE)
        light_prod.int  = lmer(log(prod_d)~scale(mean_light)+(1|site_id), data = na.omit(df_grab), REML = FALSE)
        #likelihood test for best model
        mod_an = anova(prod_null,
                       prod_full_int,
                       prod_full_add,
                       prod_full_tempxlight,
                       chla_prod.int,
                       temp_prod.int,
                       light_prod.int)
        mod_min <- mod_an %>% rownames_to_column('Model') %>% filter(AIC == min(AIC)) %>% select(Model:logLik)
        model_sel_df[col,1:5] <- unlist(mod_min)
        mod_name = mod_min %>% pull(Model)
        model_sel_df[col,6:7] = MuMIn::r.squaredGLMM(eval(parse(text = mod_name)))
        model_list[[col]] <- eval(parse(text = mod_name))
      }
      model_sel_df = model_sel_df %>% data.frame() %>%
        setNames(., c("Model", "df","AIC","BIC","logLik","r2_m","r2_c"))
      return(list(model_sel_df,model_list))
    }
  }
  top_mod_boot <<- function(prod,chla,df,nboot,temp_var = "tempkt_stand",...){#input wide stream_prod_df and wide chla_int_df
    if(nboot > (dim(prod)[2]-5)) {print('not enough observations')} else{
      df = df %>% select(site_id, DATE, mean_light, !!as.symbol(temp_var))
      model_sel_df = matrix(ncol = 6, nrow = nboot)
      model_list = vector('list', nboot)
      for(col in 1:nboot){
        prod_d_grab = prod %>% select(site_id, DATE, prod_d = (col+5))
        chla_grab = chla %>% select(site_id, DATE, chla_mean = (col+2))
        df_grab = list(df, prod_d_grab, chla_grab) %>% reduce(left_join, by = c('site_id','DATE'))
        #models
        formula = as.formula(paste("log(prod_d)~log(chla_mean)+",temp_var,"+scale(mean_light)+(1|site_id)"))
        # prod_full_add= lmer(log(prod_d)~log(chla_mean)+as.symbol(temp_var)+scale(mean_light)+(1|site_id), data = na.omit(df_grab), REML = TRUE)
        prod_full_add= lmer(formula, data = na.omit(df_grab), REML = TRUE)
        #likelihood test for best model
        model_sel_df[col,1] <- "prod_full_add"
        model_sel_df[col,2] <- AIC(prod_full_add)
        model_sel_df[col,3] <- BIC(prod_full_add)
        model_sel_df[col,4] <- logLik(prod_full_add)[1]
        model_sel_df[col,5:6] = MuMIn::r.squaredGLMM(prod_full_add)
        model_list[[col]] <- eval(parse(text = 'prod_full_add'))
      }
      model_sel_df = model_sel_df %>% data.frame() %>%
        setNames(., c("Model", "AIC","BIC","logLik","r2_m","r2_c"))
      return(list(model_sel_df,model_list))
    }
  }
  back_mod_sel_boot <<- function(prod,chla,df,nboot,...){#input wide stream_prod_df and wide chla_int_df
    if(nboot > (dim(prod)[2]-5)) {print('not enough observations')} else{
      df = df %>% select(site_id, DATE, mean_light, tempkt_stand)
      model_sel_df = matrix(ncol = 6, nrow = nboot)
      model_list = vector('list', nboot)
      for(col in 1:nboot){
        prod_d_grab = prod %>% select(site_id, DATE, prod_d = (col+5))
        chla_grab = chla %>% select(site_id, DATE, chla_mean = (col+2))
        df_grab = na.omit(list(df, prod_d_grab, chla_grab) %>% reduce(left_join, by = c('site_id','DATE')))
        #models
        model_str = as.formula("log(prod_d)~log(chla_mean)*tempkt_stand*scale(mean_light)+(tempkt_stand|site_id)+(log(chla_mean)|site_id)+(scale(mean_light)|site_id)")
        full.model = do.call("lmer",args = list(formula = model_str, data = df_grab, REML = FALSE))
        step_model = lmerTest::step(full.model)
        model_obj = get_model(step_model)         
        mod_min = model_obj@call[["formula"]]
        model_sel_df[col,1] = as.character(mod_min[3])
        model_sel_df[col,2] = AIC(logLik(model_obj))
        model_sel_df[col,3] = BIC(logLik(model_obj))
        model_sel_df[col,4] = logLik(model_obj)
        model_sel_df[col,5:6] = MuMIn::r.squaredGLMM(model_obj)
        model_list[[col]] <- model_obj
      }
      model_sel_df = model_sel_df %>% data.frame() %>%
        setNames(., c("Model","AIC","BIC","logLik","r2_m","r2_c"))
      return(list(model_sel_df,model_list))
    }
  }
  ## subsetting light availability for each sampling interval ##
  light_subset <<- function(int_df, tempkt_light.doy){
    sum_light = matrix(NA, ncol = 1, nrow = nrow(int_df)+1)
    mean_light = matrix(NA, ncol = 1, nrow = nrow(int_df)+1)
    for(i in 1:nrow(int_df)){
      if(as.character(int_df$site_id[i]) != as.character(int_df$site_id[i+1]) | is.na(int_df$site_id[i+1])){
        light_sub = subset(tempkt_light.doy, doy >= int_df$doy[i] & doy <= (int_df$doy[i]+int_df$days))
      } else if(int_df$doy[i] > int_df$doy[i+1]){
        light_sub = subset(tempkt_light.doy, doy >= int_df$doy[i] | doy <= int_df$doy[i+1])
      } else {
        light_sub = subset(tempkt_light.doy, doy >= int_df$doy[i] & doy <= int_df$doy[i+1])
      }
      mean_light[i] = mean(light_sub[,'Intensity'], na.rm = TRUE)
      sum_light[i] = sum(light_sub[,'Intensity'], na.rm = TRUE)
    }
    return(cbind(mean_light, sum_light))
  }
  ### subsetting the summing the GPP for st7 and oh2 ###
  gpp_subset <<- function(int_df, gpp_df, var, wrap = TRUE){
    dates <- int_df[,c("DATE","JULIAN")]
    dates <- dates %>% arrange(JULIAN)
    int_gpp = matrix(NA, ncol = 1, nrow = nrow(dates)+1)
    if (wrap == TRUE){
      last.date = which(dates$JULIAN == max(dates$JULIAN))
      first.date = which(dates$JULIAN == min(dates$JULIAN))
      j.wrap = as.numeric(dates[last.date,'JULIAN'] + (365-(dates$JULIAN[last.date]-dates$JULIAN[first.date])))
      d.julian = unlist(month.day.year(j.wrap, origin = c(month = 01, day = 01, year = 2010)))
      d.wrap = as.POSIXct(paste0(d.julian['year'],"-",d.julian['day'],"-",d.julian['month']), format = "%Y-%d-%m", tz = "UTC")
      wraps = data.frame(DATE = as.Date(d.wrap), JULIAN = j.wrap)
      dates = bind_rows(dates, wraps)
      dates <- dates %>% mutate(DATE = as.Date(DATE))
      dates[,var] = NA_real_
      julians = c(dates$JULIAN ,j.wrap)
    }
    for(i in seq_along(julians)){
      temp_sub = subset(gpp_df, julian(DATE,origin = as.Date("2010-01-01")) >= julians[i] & julian(DATE, origin = as.Date("2010-01-01")) <= julians[i+1])
      dates[i,var] = temp_sub %>% summarise_at(vars(all_of(var)), funs(sum), na.rm = TRUE)
    }
    dates = dates[-dim(dates)[1],]
    return(data.frame(dates))
  }
  #top_mod_boot_center
  top_mod_boot_center <<- function(prod,chla,df,nboot,temp_var = "tempkt",...){#input wide stream_prod_df and wide chla_int_df
    if(nboot > (dim(prod)[2]-5)) {print('not enough observations')} else{
      df = df %>% select(site_id, DATE, mean_light, !!as.symbol(temp_var))
      model_sel_df = matrix(ncol = 6, nrow = nboot)
      model_list = vector('list', nboot)
      for(col in 1:nboot){
        prod_d_grab = prod %>% select(site_id, DATE, prod_d = (col+5))
        chla_grab = chla %>% select(site_id, DATE, chla_mean = (col+2))
        df_grab = list(df, prod_d_grab, chla_grab) %>% reduce(left_join, by = c('site_id','DATE'))
        #models
        # prod_full_add= lmer(log(prod_d)~log(chla_mean)+tempkt_center+log(mean_light)+(1|site_id), data = na.omit(df_grab), REML = TRUE)
        formula = as.formula(paste("log(prod_d) ~ log(chla_mean) +",temp_var,"+ log(mean_light)+(1|site_id)"))
        # prod_full_add= lmer(log(prod_d)~temp_var+log(mean_light)+(1|site_id), data = na.omit(df_grab), REML = TRUE)
        prod_full_add= lmer(formula, data = na.omit(df_grab), REML = TRUE)
        #likelihood test for best model
        model_sel_df[col,1] <- "prod_full_add"
        model_sel_df[col,2] <- AIC(prod_full_add)
        model_sel_df[col,3] <- BIC(prod_full_add)
        model_sel_df[col,4] <- logLik(prod_full_add)[1]
        model_sel_df[col,5:6] = MuMIn::r.squaredGLMM(prod_full_add)
        model_list[[col]] <- eval(parse(text = 'prod_full_add'))
      }
      model_sel_df = model_sel_df %>% data.frame() %>%
        setNames(., c("Model", "AIC","BIC","logLik","r2_m","r2_c"))
      return(list(model_sel_df,model_list))
    }
  }
  ##model selection for warming 
  warming_mod_sel_boot <<- function(gpp, prod, temp_var = "tempkt_stand", nboot = 1000,...){
    prod_cols = names(prod[grepl('V', names(prod))])
    gpp_cols = names(gpp[grepl("EstGPP_",names(gpp))])
    gpp_df = matrix(ncol = 3, nrow = nboot)
    resid_df = matrix(ncol = 2, nrow = nboot)
    for(i in 1:nboot){
      gpp = warming_intdf %>% select(site_id,DATE, all_of(temp_var),chla_mean, mean_light,!!as.symbol(gpp_cols[i])) %>% mutate(gpp_d = !!as.symbol(gpp_cols[i]))
      prod = prod_df %>% select(site_id, DATE,!!as.symbol(prod_cols[i])) %>% mutate(prod_d = !!as.symbol(prod_cols[i]))
      gen_df = suppressMessages(left_join(gpp,prod))
      gpp.lm = lm(log(prod_d)~log(gpp_d), data = gen_df)
      x = augment(gpp.lm) %>% bind_cols(gen_df)
      temp.lm = lm(`.resid`~tempkt_stand, data = x)
      gpp_df[i,1] = summary(gpp.lm)[['coefficients']]['(Intercept)','Estimate']
      gpp_df[i,2] = summary(gpp.lm)[['coefficients']]['log(gpp_d)','Estimate']
      gpp_df[i,3] = summary(gpp.lm)[['adj.r.squared']]
      resid_df[i,1] = summary(temp.lm)[['coefficients']]['tempkt_stand','Estimate']
      resid_df[i,2] = summary(temp.lm)[['adj.r.squared']]
      
    }
    
    gpp_df = gpp_df %>% data.frame() %>%
      setNames(., c("int", "tempkt_stand","adj_r2"))
    resid_df = resid_df %>% data.frame %>%
      setNames(.,c("adjusted","r2"))
    return(list(gpp_df, resid_df))
  }
  ##### publication specific details #####
  ## Ecology Letters MS details ###
  # Figure guidelines: 82mm (1/2 pg width), 110mm (2/3 pg width), 173mm (full pg).
  # Manuscript guidelines: max 5000 words (not including Abstract, Acknowledgements, References, etc); 
  # 6 figures and tables, text boxes, etc. Abstract limit: 150 words
  half_page_fig <<- 82/25.4
  two_thirds_fig <<- 110/25.4
  full_fig <<- 173/25.4
  ###
}
load_packages()