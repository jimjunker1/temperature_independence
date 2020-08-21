source("./code/02_data-import-clean.R")
## Figure 1 code ##
# ####### annual scaling of production figure ####
### annual plots ###
#Eann = -0.67
#p0[15] = 8.5
## create all boot estimates df ##
ann_prod_lmvar = prod_boot_C %>%
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(prod_boot_C), temp_seq = seq(min(annual_summ$tempC),
                                                                 max(annual_summ$tempC), length.out = 4))) %>%
  mutate(lm_pred = tempC*temp_seq + Intercept)

## create full plot
annual_prod_plot = ggplotGrob(ggplot(annual_summ, aes(x = tempC, y = log(prod_mg_m_y))) +
                                geom_line(data = ann_prod_lmvar, aes(x = temp_seq, y = lm_pred, group = id), color = 'darkgrey', alpha = 0.02)+
                                geom_smooth(method = "lm", se = FALSE, colour = 'black')+
                                geom_errorbar(aes(ymin = log(prod_min), ymax = log(prod_max), colour = site_id), width = 0, size = 1) +
                                geom_errorbar(aes(ymin = log(prod_quant2.5), ymax = log(prod_quant97.5), colour = site_id), width = 0, size = 3.2)+
                                geom_point(aes(group = site_id, fill = site_id), size = 2.5, colour = c(rep('black',5),'white'), shape = 21)+
                                scale_x_continuous(name = expression("Temperature ("*degree*"C)"), limits = c(0,30),
                                                   sec.axis = dup_axis(
                                                     breaks = breaks_in_C_plus,
                                                     labels = breaks_in_kt_stand,
                                                     name = expression("Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])")
                                                   ), expand = c(0.03,0)) +
                                scale_y_continuous(name = expression("Secondary Production (mg AFDM "~m^-2~y^-1*")"),limits = c(4,11),
                                                   breaks = breaks_in_logprod, labels = breaks_in_prod,expand = c(0.03,0)) +
                                scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                                scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels, guide = FALSE)+
                                theme_tufte(ticks = TRUE) +
                                geom_rangeframe(sides = "tlb") +
                                annotate('text', x =0, y = log(49020.8), family = 'serif', label = 'A', hjust = 0)+
                                annotate('text', x = 0, y = log(30000), family = 'serif', hjust = 0, size = 3, label = expr(~italic(E)*p[among]~"="))+  
                                annotate('text', x = 0, y = log(21000), family = 'serif', 
                                         label = expr(!!round(mean(prod_boot_df$tempkt),2)~"("*!!round(quantile(prod_boot_df$tempkt, 0.025), 2)~"–"~!!round(quantile(prod_boot_df$tempkt, 0.975),2)*")"),
                                         hjust = 0, size = 3)+
                                theme(legend.title = element_blank(), legend.position = c(0.8,0.20),
                                      legend.spacing.y = unit(0.0003,'cm'),
                                      legend.key.size = unit(0.4,'cm')))

####
#Ebioann =
## create boot B and PB estimates for plotting ##
ann_bio_lmvar = bio_boot_C %>%
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(bio_boot_C), temp_seq = seq(min(annual_summ$tempC),
                                                                max(annual_summ$tempC), length.out = 4))) %>%
  mutate(lm_pred = tempC*temp_seq + Intercept)

ann_pb_lmvar = pb_boot_C %>%
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(pb_boot_C), temp_seq = seq(min(annual_summ$tempC),
                                                               max(annual_summ$tempC), length.out = 4))) %>%
  mutate(lm_pred = tempC*temp_seq + Intercept)

## create full plot ##
annual_biopb_plot = ggplotGrob(ggplot(annual_summ, aes(x = tempC, y = log(bio_mg_m))) +
                                 geom_line(data = ann_bio_lmvar, aes(x = temp_seq, y = lm_pred, group = id), color = 'darkgrey', alpha = 0.02)+
                                 geom_smooth(method = 'lm', se = FALSE, colour = 'black') +
                                 geom_line(data = ann_pb_lmvar, aes(x = temp_seq, y = lm_pred, group = id), color = 'darkgrey', alpha = 0.02)+
                                 geom_smooth(data = annual_summ, aes(x = tempC, y = log(pb)), method = 'lm', se = FALSE, colour = 'black')+
                                 geom_errorbar(aes(ymin = log(bio_min), ymax = log(bio_max), colour = site_id),width = 0, size = 1)+
                                 geom_errorbar(aes(ymin = log(bio_quant2.5), ymax = log(bio_quant97.5), colour = site_id),width = 0, size = 3.2)+
                                 geom_errorbar(aes(ymin = log(pb_min), ymax = log(pb_max), colour = site_id),width = 0, size = 1)+
                                 geom_errorbar(aes(ymin = log(pb_quant2.5), ymax = log(pb_quant97.5), colour = site_id),width = 0, size = 3.2)+
                                 geom_point(aes(group= site_id, fill = site_id), shape = 21, colour = c(rep('black',5),'white'), size = 2.5)+
                                 geom_point(data = annual_summ, aes(x= tempC, y = log(pb), group = site_id, fill = site_id), shape = 22, colour = c(rep('black',5), 'white'), size = 2.5)+
                                 scale_x_continuous(name = expression("Temperature ("*degree*"C)"), limits = c(0,30),
                                                    sec.axis = dup_axis(
                                                      breaks = breaks_in_C_plus,
                                                      labels = breaks_in_kt_stand,
                                                      name = expression("Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])")
                                                    ), expand = c(0.03,0)) +
                                 scale_y_continuous(name = expression("B | P:B "[mass-corrected]), limits = c(-3.5,11), position = "right",
                                                    breaks = breaks_in_logprod, labels = breaks_in_prod) +
                                 scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos]) +
                                 scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels, guide = FALSE)+
                                 theme_tufte(ticks = TRUE) +
                                 geom_rangeframe(sides = "trb") + geom_rangeframe(sides = "r",aes(x = tempC, y = log(pb))) +
                                 annotate('text', x =29, y = log(49020.8), family = 'serif', label = 'B',hjust = 1) +
                                 annotate('text', x = 17, y = log(21500), family = 'serif', label = expr(~italic(E)*b[among]~"="), hjust = 0, size = 3)+  
                                 annotate('text', x = 17, y = log(10000), family = 'serif',
                                          label = expr(!!round(mean(bio_boot_df$tempkt_stand),2)~"("*!!round(quantile(bio_boot_df$tempkt_stand,0.025),2)~"–"~!!round(quantile(bio_boot_df$tempkt_stand,0.975),2)*")"),
                                          hjust = 0, size = 3)+
                                 annotate('text', x = 17, y = log(0.22), family = 'serif', label = expr(~italic(E)*pb[among]~"="), hjust = 0, size = 3)+
                                 annotate('text', x = 17, y = log(0.10), family = 'serif', 
                                          label = expr(!!round(mean(pb_boot_df$tempkt_stand),2)~"("*!!round(quantile(pb_boot_df$tempkt_stand,0.025),2)~"–"~!!round(quantile(pb_boot_df$tempkt_stand,0.975),2)*")"),
                                          hjust = 0, size = 3)+
                                 theme(legend.title = element_blank(), legend.position = "none"))

maxwidth = grid::unit.pmax(annual_prod_plot$widths[2:5], annual_biopb_plot$widths[2:5])
annual_prod_plot$widths[2:5] = as.list(maxwidth)
annual_biopb_plot$widths[2:5] = as.list(maxwidth)

# tiff("./output/annual_prod_pub.tiff", res = 600, height = 4, width = full_fig, units = "in", compression = "lzw")
grid.draw(gridExtra::gtable_cbind(annual_prod_plot,annual_biopb_plot, size = "max"))
# dev.off()

## Figure 2 ##
### plotting within stream estimates ###
## annual prediction from linear regression on mean for plotting ##
temp_prod.lm = lm(log(prod_mg_m_y)~tempkt_stand, data = annual_summ);summary(temp_prod.lm)
annual_prod_predict = data.frame(tempkt_stand = seq(min(int_df$tempkt_stand), max(int_df$tempkt_stand), length.out = 100))
annual_prod_predict = annual_prod_predict %>%
  mutate(prod_d= exp(predict(temp_prod.lm, newdata = annual_prod_predict))/365,
         tempC= overkt_stand15_to_C(tempkt_stand))
ann_prod_Eavar = prod_boot_df %>% 
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(prod_boot_C), temp_seq = seq(min(annual_summ$tempkt_stand),
                                                                 max(annual_summ$tempkt_stand), length.out = 4))) %>%
  mutate(lm_pred = tempkt_stand*temp_seq + Intercept)
## create plot
int_tempprod_plot = ggplotGrob(ggplot(int_df, aes(x = overkt_stand15_to_C(tempkt_stand), y = log(prod_d), linetype = site_id)) + 
                                 geom_line(data = annual_prod_predict, aes(x = overkt_stand15_to_C(tempkt_stand), y = log(prod_d)), colour = "darkgrey", linetype = 'dashed')+
                                 geom_path(data = prod_within_pred, aes(x = overkt_stand15_to_C(tempkt_stand), y = prod_d_new, colour = site_id, group = site_id), size=0.9)+
                                 geom_point(data = int_df,aes(group = site_id, fill = site_id),size = 2.2, colour = "black", shape = 21) +
                                 scale_x_continuous(name = expression("Temperature ("*degree*"C)"), 
                                                    sec.axis = dup_axis(
                                                      breaks = breaks_in_C_plus,
                                                      labels = breaks_in_kt_stand,
                                                      name = expression("Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])")),
                                                    expand = c(0.03,0)) +
                                 scale_y_continuous(name = expression("Secondary Production (mg"~m^-2~d^-1*")"), #position = "right",
                                                    breaks = breaks_in_log10prod, labels = breaks_in_prod10) +
                                 coord_cartesian(xlim = c(0,35), ylim =  c(-3.5,6))+
                                 scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                                 scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels)+
                                 scale_linetype_manual(values = c('solid','solid','solid','solid','solid','solid')) +
                                 guides(linetype = FALSE)+
                                 theme_tufte(ticks = TRUE) +
                                 geom_rangeframe(sides = "tlb") +
                                 theme(legend.title = element_blank(), legend.position = c(0.1, 0.8),
                                       legend.spacing.y = unit(0.0003,'cm'),
                                       legend.key.size = unit(0.4,'cm')))

### create inset of within stream apparent temp dependence bounds  ###
limits_wide = aes(ymin = (prod_max_est), ymax = (prod_min_est), colour = site_id)
limits_narrow = aes(ymin = (prod_quant_2.5), ymax = (prod_quant_97), colour = site_id)
ann_prod_quants = quantile(ann_prod_Eavar[,'tempkt_stand'], probs = c(0,1))

prod_temp_dep_plot = ggplotGrob(ggplot(prod_d_stream_coef_df, aes(x = site_id, y = (prod_slope_est))) + 
                                  geom_rect(xmin =0,xmax = 7,
                                            ymin = min(ann_prod_quants), ymax = max(ann_prod_quants), fill = 'grey', alpha = 0.1)+
                                  geom_errorbar(limits_narrow, width = 0, size = 3)+
                                  geom_hline(aes(yintercept = (coef(temp_prod.lm)['tempkt_stand'])), size = 1.2, linetype = 'dashed', colour = "darkgrey")+
                                  geom_errorbar(limits_wide, width = 0, size = 1)+
                                  geom_point(aes(fill = site_id), colour = c(rep('black',5),'white'),size = 1.8, shape = 21) +
                                  coord_flip(ylim = c(0,5.5)) + 
                                  scale_y_continuous( expand = c(0.01,0))+
                                  scale_colour_manual(values = ocecolors[['temperature']][rev(oce_temp_pos)])+
                                  scale_fill_manual(values = ocecolors[['temperature']][rev(oce_temp_pos)])+
                                  theme_tufte(ticks = TRUE) +
                                  # geom_rangeframe(sides = "lb") +
                                  annotate("text", x = 5.7, y = 5.5, hjust = 1, vjust = 0, label = expr(italic(E)[p]), family = 'serif')+  
                                  theme(legend.position = 'none',strip.background = element_rect(fill = 'white'), 
                                        axis.title = element_blank(), axis.text.y = element_blank(), 
                                        text = element_text(family = 'serif'), plot.title = element_text(hjust = 0.5, size = 11),
                                        panel.background = element_rect(fill = 'transparent', colour = NA),
                                        plot.background = element_rect(fill = 'transparent', colour = NA)))

prod_gm = int_tempprod_plot
prod_gs = prod_temp_dep_plot
prod_id = 1
prod_panel = prod_gm$layout[prod_gm$layout$name == 'panel',][prod_id,]

prod_inset = grobTree(prod_gs, vp = viewport(width = 0.42, height = 0.50, x = 0.8, y = 0.25))
prod_gm = gtable::gtable_add_grob(prod_gm, prod_inset, l = prod_panel$l, t = prod_panel$t)

# tiff("./output/int_prod_pub.tiff", res = 600, height = two_thirds_fig, width = two_thirds_fig, units = "in", compression = "lzw")
grid.draw(prod_gm)
# dev.off()

## Figure 3 ##
### interval corrected biomass plot ###
## annual prediction from linear regression on mean for plotting ##
temp_bio_corr.lm = lm(log(bio_mg_m)~tempkt_stand, data = annual_summ)
annual_bio_predict = data.frame(tempkt_stand = seq(min(int_df$tempkt_stand), max(int_df$tempkt_stand), length.out = 100))
annual_bio_predict = annual_bio_predict %>%
  mutate(corrbio_mg_m= exp(predict(temp_bio_corr.lm, newdata = annual_bio_predict)),
         tempC= overkt_stand15_to_C(tempkt_stand))

ann_bio_Eavar = bio_boot_df %>% 
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(prod_boot_C), temp_seq = seq(min(annual_summ$tempkt_stand),
                                                                 max(annual_summ$tempkt_stand), length.out = 4))) %>%
  mutate(lm_pred = tempkt_stand*temp_seq + Intercept)

## create plot
int_bio_plot = ggplotGrob(ggplot(int_df, aes(x = overkt_stand15_to_C(tempkt_stand), y = log(corrbio_mg_m))) + 
                            geom_line(data = annual_bio_predict, colour = "darkgrey", linetype = "dashed")+
                            geom_path(data = prod_within_pred, 
                                      aes(x = overkt_stand15_to_C(tempkt_stand), y = bio_new, colour = site_id, group = site_id, linetype = site_id), size=0.9)+
                            geom_point(aes(group = site_id, fill = site_id),size = 2.2, colour = "black", shape = 21) +
                            scale_x_continuous(name = expression("Temperature ("*degree*"C)"), limits = c(0,35),
                                               sec.axis = dup_axis(
                                                 breaks = breaks_in_C_plus,
                                                 labels = breaks_in_kt_stand,
                                                 name = expression("Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])")
                                               ), expand = c(0.03,0)) +
                            scale_y_continuous(name = expression("Biomass (mg AFDM "~m^-2*")"[mass-corrected]), limits = c(1.5,12), 
                                               breaks = breaks_in_log10prod, labels = breaks_in_prod10) +
                            scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos]) +
                            scale_colour_manual(values  = ocecolors[['temperature']][oce_temp_pos]) +
                            scale_linetype_manual(values = c('solid','solid','solid','solid','solid','solid')) +
                            guides(linetype = FALSE) +
                            theme_tufte(ticks = TRUE) +
                            geom_rangeframe(sides = "tlb") +
                            annotate('text', x = 1, y = 12, family = 'serif', label = 'A') +
                            theme(legend.title = element_blank(), legend.position = "none"))

### biomass temp dependence inset ###
limits_wide = aes(ymin = (bio_max_est), ymax = (bio_min_est), colour = site_id)
limits_narrow = aes(ymin = (bio_quant_2.5), ymax = (bio_quant_97), colour = site_id)
ann_bio_quants = quantile(ann_bio_Eavar[,'tempkt_stand'], probs = c(0,1))

bio_temp_dep_plot = ggplotGrob(ggplot(bio_stream_coef_df, aes(x = site_id, y = (bio_slope_est))) + 
                                 geom_rect(xmin =0,xmax = 7,
                                           ymin = min(ann_bio_quants), ymax = max(ann_bio_quants), fill = 'grey', alpha = 0.1)+
                                 geom_errorbar(limits_narrow, width = 0, size = 3)+
                                 # geom_hline(aes(yintercept= z), size = 1.2, colour = "darkgrey") +# fixed effects estimate (mean of within-stream estimates)
                                 geom_hline(aes(yintercept = (coef(temp_bio_corr.lm)['tempkt_stand'])), size = 1.2, linetype = 'dashed', colour = "darkgrey")+
                                 geom_errorbar(limits_wide, width = 0, size = 1)+
                                 geom_point(aes(fill = site_id), colour = c(rep('black',5),'white'),size = 1.4, shape = 21) +
                                 coord_flip() + #limits = c(-1.5, 4.5)
                                 scale_y_continuous(expand = c(0.01,0))+
                                 scale_colour_manual(values = ocecolors[['temperature']][rev(oce_temp_pos)])+
                                 scale_fill_manual(values = ocecolors[['temperature']][rev(oce_temp_pos)])+
                                 annotate("text", x = 5.6, y = 6, hjust = 1, vjust = 0, label = expr(italic(E)[b]), family = 'serif')+  
                                 theme(legend.position = 'none',strip.background = element_rect(fill = 'white'), 
                                       axis.title = element_blank(), axis.text.y = element_blank(), 
                                       text = element_text(family = 'serif'), plot.title = element_text(hjust = 0.5, size = 11),
                                       panel.background = element_rect(fill = 'transparent', colour = NA),
                                       plot.background = element_rect(fill = 'transparent', colour = NA)))

bio_gm = int_bio_plot
bio_gs = bio_temp_dep_plot
bio_id = 1
bio_panel = bio_gm$layout[bio_gm$layout$name == 'panel',][bio_id,]

bio_inset = grobTree(bio_gs, vp = viewport(width = 0.42, height = 0.45, x = 0.78, y = 0.79))
bio_gm = gtable::gtable_add_grob(bio_gm, bio_inset, l = bio_panel$l, t = bio_panel$t)

### interval corrected pb plot ###
## annual prediction from linear regression on mean for plotting ##
temp_pb_corr.lm = lm(log(pb)~tempkt_stand, data = annual_summ)
annual_pb_predict = data.frame(tempkt_stand = seq(min(int_df$tempkt_stand), max(int_df$tempkt_stand), length.out = 100))
annual_pb_predict = annual_pb_predict %>%
  mutate(corrpb_d= exp(predict(temp_pb_corr.lm, newdata = annual_pb_predict)),
         tempC= overkt_stand15_to_C(tempkt_stand))

ann_pb_Eavar = pb_boot_df %>% 
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(prod_boot_C), temp_seq = seq(min(annual_summ$tempkt_stand),
                                                                 max(annual_summ$tempkt_stand), length.out = 4))) %>%
  mutate(lm_pred = tempkt_stand*temp_seq + Intercept)

## create plot
int_pb_plot = ggplotGrob(ggplot(int_df, aes(x = overkt_stand15_to_C(tempkt_stand), y = log(corrpb_d))) + 
                           geom_line(data = annual_pb_predict,colour = "darkgrey", linetype = "dashed")+
                           geom_path(data = prod_within_pred, 
                                     aes(x = overkt_stand15_to_C(tempkt_stand), y = pb_new, colour = site_id, group = site_id, linetype = site_id), size=0.9)+
                           geom_point(aes(group = site_id, fill = site_id),size = 2.2, colour = "black", shape = 21) +
                           scale_x_continuous(name = expression("Temperature ("*degree*"C)"), limits = c(0,35),
                                              sec.axis = dup_axis(
                                                breaks = breaks_in_C_plus,
                                                labels = breaks_in_kt_stand,
                                                name = expression("Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])")
                                              ), expand = c(0.03,0)) +
                           scale_y_continuous(name = expression("P:B ("*d^-1*")"[mass-corrected]), position = "right", breaks = breaks_in_log10prod, labels = breaks_in_prod10) +
                           scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                           scale_colour_manual(values  = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                           scale_linetype_manual(values = c('solid','solid','solid','solid','solid','solid')) +
                           guides(linetype = FALSE)+
                           theme_tufte(ticks = TRUE) +
                           geom_rangeframe(sides = "trb") +
                           annotate('text', x = 35, y = 5.5, family = 'serif', label = 'B') +
                           theme(legend.title = element_blank(), legend.position = c(0.85, 0.23),
                                 legend.spacing.y = unit(0.0003,'cm'),
                                 legend.key.size = unit(0.4,'cm')))

### pb temp dependence inset ###
limits_wide = aes(ymin = (pb_max_est), ymax = (pb_min_est), colour = site_id)
limits_narrow = aes(ymin = (pb_quant_2.5), ymax = (pb_quant_97), colour = site_id)
ann_pb_quants = quantile(ann_pb_Eavar[,'tempkt_stand'], probs = c(0,1))

pb_temp_dep_plot = ggplotGrob(ggplot(pb_stream_coef_df, aes(x = site_id, y = (pb_slope_est))) + 
                                geom_rect(xmin =0,xmax = 7,
                                          ymin = min(ann_pb_quants), ymax = max(ann_pb_quants), fill = 'grey', alpha = 0.1)+
                                geom_hline(aes(yintercept = (coef(temp_pb_corr.lm)['tempkt_stand'])), size = 1.2, linetype = 'dashed', colour = "darkgrey")+
                                geom_errorbar(limits_narrow, width = 0, size = 3)+
                                geom_errorbar(limits_wide, width = 0, size = 1)+
                                geom_point(aes(fill = site_id), colour = c(rep('black',5),'white'),size = 1.4, shape = 21) +
                                coord_flip() + #limits = c(-1, 2),
                                scale_y_continuous( expand = c(0.01,0))+
                                scale_x_discrete(position = 'top')+
                                scale_colour_manual(values = ocecolors[['temperature']][rev(oce_temp_pos)])+
                                scale_fill_manual(values = ocecolors[['temperature']][rev(oce_temp_pos)])+
                                annotate("text", x = 5.6, y = 5.8, hjust = 1, vjust = 0, label = expr(italic(E)[pb]), family = 'serif')+  
                                # facet_wrap(~facet, scales = 'free') + 
                                # ggtitle('Temperature dependence')+
                                theme(legend.position = 'none',strip.background = element_rect(fill = 'white'), 
                                      axis.title = element_blank(), axis.text.y = element_blank(), 
                                      text = element_text(family = 'serif'), plot.title = element_text(hjust = 0.5, size = 11),
                                      panel.background = element_rect(fill = 'transparent', colour = NA),
                                      plot.background = element_rect(fill = 'transparent', colour = NA)))#;grid.draw(pb_temp_dep_plot)

pb_gm = int_pb_plot
pb_gs = pb_temp_dep_plot
pb_id = 1
pb_panel = pb_gm$layout[pb_gm$layout$name == 'panel',][pb_id,]

pb_inset = grobTree(pb_gs, vp = viewport(width = 0.42, height = 0.45, x = 0.15, y = 0.79))
pb_gm = gtable::gtable_add_grob(pb_gm, pb_inset, l = pb_panel$l, t = pb_panel$t)

maxwidth = grid::unit.pmax(int_bio_plot$widths[2:5], int_pb_plot$widths[2:5])
int_bio_plot$widths[2:5] = as.list(maxwidth)
int_pb_plot$widths[2:5] = as.list(maxwidth)

# tiff("./output/int_biopb_pub.tiff", res = 600, height = 4, width = full_fig, units = "in", compression = "lzw")
grid.draw(gridExtra::gtable_cbind(bio_gm,pb_gm, size = "max"))
# dev.off()

## Figure 4 ##
### log chlorophyll and temperature ###

## create individual boot df for plotting ##
ann_chla_lmvar = chla_temp_boots_df %>%
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(chla_temp_boots_df),
                        temp_seq = seq(min(C_to_overkt_stand15(overkt_to_C(chla_ann_summ$tempkt))),
                                       max(C_to_overkt_stand15(overkt_to_C(chla_ann_summ$tempkt))),length.out = 4))) %>%
  mutate(lm_pred = tempkt_stand*temp_seq+Intercept)

## creat full plot ##
chla_temp = ggplotGrob(ggplot(chla_summ, aes(x = overkt_to_C(tempkt), y = log(mean_chla_mg_m))) +
                         geom_line(data = ann_chla_lmvar, aes(x = overkt_stand15_to_C(temp_seq), y = lm_pred, group = id), color = 'darkgrey', alpha = 0.01 )+
                         geom_smooth(data = chla_ann_summ, aes(x = overkt_to_C(tempkt), y = log(chla_mean_ann)), method= "lm", se = FALSE, colour = "black", size = 1)+
                         geom_point(aes(fill = site_id), size = 2, shape= 22, alpha = 0.85, colour = 'darkgrey') +
                         geom_errorbar(data = annual_summ, aes(x = overkt_to_C(tempkt), ymin = log(ann_chla_quant2.5),
                                                               ymax = log(ann_chla_quant97), colour = site_id), width = 0, size = 3.2, inherit.aes = FALSE) +
                         geom_errorbar(data = annual_summ, aes(x = overkt_to_C(tempkt), ymin = log(ann_chla_min),
                                                               ymax = log(ann_chla_max), colour = site_id), width = 0, size = 1, inherit.aes = FALSE) +
                         geom_point(data = chla_ann_summ, aes(x = overkt_to_C(tempkt), y = log(chla_mean_ann), group = site_id, fill = site_id),
                                    colour = c(rep("black",5),"white"), shape = 21, size = 1.4) +
                         scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                         scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels, guide = FALSE) +
                         scale_x_continuous(name = expression("Temperature ("*degree*"C)"), expand = c(0.03,0),
                                            limits = c(0,35),sec.axis = dup_axis(
                                              breaks = breaks_in_C_plus,
                                              labels = breaks_in_kt_stand,
                                              name = expression(" Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])") )) +
                         scale_y_continuous(name = expression("Chlorophyll "~italic(a)~"(mg"~m^-2*")"),limits = c(-3,6),
                                            breaks = breaks_in_logchla, labels = breaks_in_chla,expand = c(0.03,0))+
                         theme_tufte(ticks = TRUE) +
                         geom_rangeframe(sides = "tbl", colour = "black") +
                         annotate('text', x = 0, y = log(300), family = 'serif', label = 'A', hjust = 0)+
                         annotate('text', x = 0, y = log(180), family = 'serif', size = 3, label = expr(~italic(E)*b[chla]~"="~!!round(chla_temp_coefs_df$tempkt_est,2)), hjust = 0)+
                         annotate('text', x = 0, y = log(130), family = 'serif', size = 3, label = expr(~r^2~"="~!!round(chla_temp_coefs_df$adj_r2_mean,2)), hjust = 0)+
                         theme(legend.position = c(0.85, 0.23), legend.title = element_blank(),
                               legend.spacing.y = unit(0.0003,'cm'),
                               legend.key.size = unit(0.4,'cm')))

### log chlorophyll a and secondary production ###
limits_wide = aes(ymin = log(prod_min_adj), ymax = log(prod_max_adj), colour = site_id)
limits_narrow = aes(ymin = log(prod_quant2.5_adj), ymax = log(prod_quant97.5_adj), colour = site_id)
limits_h_wide = aes(xmin = log(ann_chla_quant2.5), xmax = log(ann_chla_quant97), colour = site_id)
limits_h_narrow = aes(xmin = log(ann_chla_min), xmax = log(ann_chla_max), colour = site_id)
mean(chla_prod_lmboots[,'adj_r2'])
mean(chla_prod_lmboots[,'chla'])
quantile(chla_prod_lmboots[,'chla'],c(0.025,0.975))

## create individual boot df for plotting ##
ann_chlaprod_lmvar = chla_prod_lmboots %>%
  mutate(id = 1:n()) %>%
  left_join(expand.grid(id = 1:nrow(chla_prod_lmboots), chla_seq = seq(log(min(annual_summ$ann_chla_mean)),
                                                                       log(max(annual_summ$ann_chla_mean)), length.out = 4))) %>%
  mutate(lm_pred = exp(chla*chla_seq + Intercept)/365)

## create full plot ##
log_chla_plot = ggplotGrob(ggplot(annual_summ, aes(x = log(ann_chla_mean), y = log(prod_adj))) +
                             geom_line(data = ann_chlaprod_lmvar ,aes(x = chla_seq, y = log(lm_pred), group = id),color = 'darkgrey', alpha = 0.02) +
                             geom_line(data = chla_prod_predict, aes(x = log(chla), y = log(prod_adj)), colour = 'black', size = 1, se = FALSE)+
                             geom_errorbar(limits_wide, width = 0, size = 1)+
                             geom_errorbar(limits_narrow, width = 0, size = 3.2)+
                             geom_errorbarh(limits_h_wide, height = 0, size = 3.2) +
                             geom_errorbarh(limits_h_narrow, height = 0, size = 1) +
                             geom_point(data = na.omit(int_df), aes(x = log(chla_mean),y=log(prod_d), fill = site_id),
                                        shape = 22, colour = 'darkgrey', size = 1.7, alpha = 0.9) + 
                             geom_point(aes(group = site_id, fill = site_id), shape = 21,colour = c(rep("black",5),"white"),size = 1.4) +
                             scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos]) +
                             scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], guide = FALSE) +
                             scale_y_continuous(name = expression("Secondary Production (mg AFDM  "~m^-2~d^-1~")"), expand = c(0.03,0),
                                                breaks = breaks_in_log10prod, labels = breaks_in_prod10, position = 'right') +
                             scale_x_continuous(name = expression("Chlorophyll "~italic(a)~"(mg"~m^-2*")"), expand = c(0.03,0),
                                                limits = c(-2,4),breaks = breaks_in_logchla, labels = breaks_in_chla ) +
                             theme_tufte(ticks = TRUE) +
                             geom_rangeframe(sides = 'br', na.rm = TRUE, colour = 'black') +
                             geom_rangeframe(data = na.omit(int_df), aes(x = log(chla_mean), y = log(prod_d)), sides = 'br', colour = 'black')+
                             annotate('text', x = log(50), y = log(450), family = 'serif', label = "B", hjust = 0)+
                             annotate('text', x = log(0.3), y = log(600), family = 'serif', size = 3,
                                      label = expr("log-log slope ="~!!round(mean(chla_prod_lmboots[,'chla']),1)~"("*!!round(quantile(chla_prod_lmboots[,'chla'],0.025),1)~"–"~!!round(quantile(chla_prod_lmboots[,'chla'],0.975),1)*")"), hjust = 0)+
                             annotate('text', x = log(0.3), y = log(400), family = 'serif', size = 3, label = expr(~r^2~"="~!!round(mean(chla_prod_lmboots[,'adj_r2']),2)), hjust = 0)+
                             theme(legend.position = "none"))

maxwidth = grid::unit.pmax(log_chla_plot$widths[2:5],
                           chla_temp$widths[2:5])
maxheight = grid::unit.pmax(log_chla_plot$heights[2:5],
                            chla_temp$widths)

# tiff(file = "./output/tempchla_prod_pub.tiff",res = 600, height = 4, width = full_fig, units = "in", compression = "lzw")
grid.draw(gridExtra::gtable_cbind(chla_temp, log_chla_plot,size = "max" ))
# dev.off()

## Figure 5 ##

Ea_dist_plot = ggplotGrob(ggplot(full_Ea_df_long, aes(x = time, y = Ea, group = interaction(time, type)))+
                            geom_hline(yintercept = 0, linetype = "dashed", color = 'darkgrey', size = 1.2) +
                            geom_boxplot(aes(fill = type), width = 0.3, outlier.alpha = 0.2, size = 0.6, outlier.shape = 21)+
                            scale_y_continuous(name = expression(italic(E)~"p (eV)"))+
                            scale_fill_manual(values = c("#c2a5cf","#008837"), labels = c("apparent", "resource-corrected"))+
                            scale_color_manual(values = c("#c2a5cf","#008837"), labels = c("apparent", "resource-corrected"))+
                            geom_rangeframe(sides = "lb") +
                            theme_tufte(ticks = TRUE) +
                            theme(legend.title = element_blank(), axis.title.x = element_blank(),
                                  legend.position = c(0,1),legend.justification = c(0,1), 
                                  legend.background = element_rect(fill = "transparent", colour = NA),
                                  legend.spacing.y = unit(0.0003,'cm'),
                                  legend.key.size = unit(0.4,'cm'),
                                  legend.text = element_text(size = 14, family = 'serif'),
                                  axis.text = element_text(color = 'black',family = 'serif', size = 14),
                                  axis.title.y = element_text(color = 'black', family = 'serif', size = 16)))

# tiff("./output/Ea_raw-corr_dist.tiff", res = 600, height = two_thirds_fig, width = two_thirds_fig, units = "in", compression = "lzw")
grid.draw(Ea_dist_plot)
# dev.off()


#################### SUPPLEMENTAL FIGURES ######################
## Figure S1 ##
t.plot_leg = ggplot(temp.doy_long, aes(x = doy, y = temperature)) +
  geom_line(aes(group = stream, colour = stream), size = 0.8)+
  scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels)+
  theme_tufte()+
  theme(legend.title = element_blank(), legend.spacing.y = unit(0.0003,'cm'),
        legend.key.size = unit(0.4,'cm'))
t.plot_leg
temp_leg = get_legend(t.plot_leg)

L.plot = ggplotGrob(ggplot(light.doy_long, aes(x = doy, y = Intensity/1000)) + 
                      geom_line(size = 0.6, color = 'black') + 
                      scale_y_continuous(name = "Light intensity (Lux/1000)") + 
                      scale_x_continuous(expand = c(0.03,0.1)) +
                      geom_rangeframe(sides = "l") +
                      theme_tufte(ticks = TRUE) +
                      annotate('text', x = 0, y = 65, family = 'serif', label = "A")+
                      annotation_custom(grob = temp_leg, ymin = 40, ymax = 65, xmin = 300, xmax = 365)+
                      theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(),
                            axis.ticks.x = element_blank()))

t.plot = ggplotGrob(ggplot(temp.doy_long, aes(x = doy, y = temperature)) + 
                      geom_line(aes(group = stream,colour = stream), size = 0.8) +
                      scale_y_continuous(name = expression("Temperature ("~degree*"C)"),
                                         sec.axis = dup_axis(
                                           breaks = breaks_in_C_plus,
                                           labels = breaks_in_kt_stand,
                                           name = "Standardized temperature (1/kT)"
                                         ), expand = c(0.03,0)) +
                      scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos]) +
                      scale_x_continuous(name = "Day of Year", expand = c(0.03,0)) +
                      theme_tufte(ticks = TRUE) + 
                      annotate('text', x = 0, y = 34, family = 'serif', label = 'B') +
                      geom_rangeframe(sides = "blr", na.rm = TRUE)+
                      # theme_black() +
                      theme(legend.position = "none"))

maxwidth = grid::unit.pmax(L.plot$widths[2:5], t.plot$widths[2:5])
L.plot$widths[2:5] = as.list(maxwidth)
t.plot$widths[2:5] = as.list(maxwidth)

# tiff("./output/temp_light_doy.tiff", res = 600, height = 5, width = 4.330709, units = "in", compression = "lzw")
grid.draw(gridExtra::gtable_rbind(L.plot,t.plot))
# dev.off()

## Figure S2 ##
## predicate multivariate responses
cohort_prod_summ = readRDS(file = "./data/derived-data/cohort_prod_summ.rds") %>% rename(site_id = 'SITE') %>% mutate(tempkt_stand = C_to_overkt_stand15(tempC))
IGR_lm = lm(log(IGR)~log(mean_mass)+tempkt_stand, data = cohort_prod_summ)#;summary(IGR_lm)
predicted_df <- data.frame(IGR = predict(IGR_lm, cohort_prod_summ), mean_mass=na.omit(cohort_prod_summ$mean_mass), tempkt = na.omit(cohort_prod_summ$tempkt))

## create plots
mass_corrIGR = ggplotGrob(ggplot(cohort_prod_summ, aes(x = log(mean_mass), y = log(IGR_tempcorr))) +
                            geom_smooth(method = 'lm', se = TRUE, alpha = 0.5, colour = 'black') +
                            geom_point(aes(fill = tempC),size = 2.2, colour = "black", shape = 21) + 
                            scale_x_continuous(name = "log(body mass) [mg AFDM]", limits = c(-8,4), breaks = c(-8,-6,-4,-2,0,2,4))+
                            scale_y_continuous(limits = c(-10,0), name = expression("log(Growth rate "~italic(e)^italic(E/kT)*") ["~d^-1*"]")) +
                            scale_fill_gradientn(name = expression("Temperature ("~degree*"C)"),colors = ocecolors[['temperature']])+
                            theme_tufte(ticks = TRUE) + geom_rangeframe(sides = "bl", show.legend = FALSE) +
                            annotate('text', x = -6.1, y = -7.8, hjust = 0.5, size = 3, label = expression("Temperature ("~degree*"C)"), family = 'serif')+
                            annotate('text', x = -8, y = 0, hjust = 0, size = 4, label = "A", family = 'serif')+
                            theme(legend.position = c(0.18, 0.08), legend.background = element_rect(colour = NA,fill = NA),
                                  legend.text = element_text(size = 8),legend.title = element_blank(),
                                  legend.direction = "horizontal"))


cohort_prod_summ = cohort_prod_summ %>% arrange(log(mean_mass))
temp_corrIGR = ggplotGrob(ggplot(cohort_prod_summ, aes(x = tempC, y = log(IGR_masscorr), fill = log(mean_mass))) + 
                            geom_smooth(method = "lm", colour = 'black', se = TRUE, alpha = 0.5) +    
                            geom_point(size = 2.2, colour = "black", shape = 21) + 
                            scale_x_continuous(limits = c(0,35), name = expression("Temperature ("~degree*"C)"), expand = c(0.03,0),
                                               sec.axis = dup_axis(
                                                 breaks = breaks_in_C_plus,
                                                 labels = breaks_in_kt,
                                                 name = expression("Standardized temperature (1/"*italic(kT)*")"))) +
                            scale_y_continuous(limits = c(-10,0), name = expression("log(Growth rate "~italic(M)^-0.25~")["*d^-1*"]"), position = 'right') +
                            scale_fill_viridis_c(name = "log(mass)\n[mg AFDM]", direction = 1) +
                            theme_tufte(ticks = TRUE) + geom_rangeframe(sides = "trb") +
                            annotate('text', x = 28.5, y = -7.9, hjust = 0.5, family = 'serif', size = 3, label ="log(mass) [mg AFDM]" )+
                            annotate('text', x = 34, y = 0, hjust = 0, family = 'serif', size = 4, label = "B") +
                            theme(legend.position = c(0.78,0.08), legend.background = element_rect(colour = NA,fill = NA),
                                  legend.text = element_text(size = 8),
                                  legend.title = element_blank(), legend.direction = "horizontal"))

maxwidth = grid::unit.pmax(mass_corrIGR$widths[2:5], temp_corrIGR$widths[2:5])
mass_corrIGR$widths[2:5] = as.list(maxwidth)
temp_corrIGR$widths[2:5] = as.list(maxwidth)
# tiff("./output/temp_mass_corrIGR.tiff", res = 600, height = 3.5, width = 10, units = "in", compression = "lzw")
grid.draw(gridExtra::gtable_cbind(mass_corrIGR,temp_corrIGR))
# dev.off()

## Figure 3 ##
spp_annual_df = read.table("./data/derived-data/spp_annual_df.txt", sep = "\t", header = TRUE) %>%
  mutate(site_id = factor(site_id, levels = stream_order))
#### annual evenness of population production figure ####
ann_even = ggplotGrob(ggplot(na.omit(spp_annual_df), aes(x = rank, y = log(prod_mg_m_y), group = site_id)) + 
                        geom_path(aes(colour = site_id), size = 1.5) +
                        geom_point(aes(fill = site_id),size = 2.5, colour = "black", shape = 21) + 
                        scale_y_continuous(name = expression("log(Secondary Production) [mg "~m^-2~yr^-1~"]"), 
                                           breaks = breaks_in_logprod, labels = breaks_in_prod) +
                        scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                        scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels)+
                        theme_tufte(ticks = TRUE) + geom_rangeframe(sides = "lb") +
                        theme(legend.title = element_blank(), legend.position = c(0.9,0.80), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
                              axis.title.x = element_blank(), legend.background = element_rect(fill = NA, colour = NA),
                              legend.spacing.y = unit(0.0003,'cm'),
                              legend.key.size = unit(0.4,'cm')))
# tiff("./output/prod_evenness.tiff", res = 600, height = 4, width = 5, units = "in", compression = "lzw")
grid.draw(ann_even)
# dev.off()

## Figure S4 ##
int_doyprod_plot = ggplotGrob(ggplot(int_df, aes(x = yday(DATE), y = prod_d, group = site_id)) + 
                                geom_smooth(aes(group = site_id, color = site_id), size = 1.1, se = FALSE, span = 0.6)+
                                geom_point(aes(group = site_id, fill = site_id),size = 2.2, colour = "black", shape = 21) +
                                scale_x_continuous(name = "Day of Year",
                                                   expand = c(0.03,0)) +
                                scale_y_continuous(name = expression("Secondary Production (mg"~m^-2~d^-1*")")) +
                                coord_cartesian()+
                                scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                                scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels)+
                                guides(linetype = FALSE)+
                                theme_tufte(ticks = TRUE) +
                                geom_rangeframe(sides = "tlb") +
                                theme(legend.title = element_blank(), legend.position = c(0.1, 0.8),
                                      legend.spacing.y = unit(0.0003,'cm'),
                                      legend.key.size = unit(0.4,'cm')))

# tiff("./output/doy_prod_pub.tiff", res = 600, height = two_thirds_fig, width = two_thirds_fig, units = "in", compression = "lzw")
grid.draw(int_doyprod_plot)
# dev.off()

## Figure S5 ##
chla_temp_int = ggplotGrob(ggplot(chla_summ, aes(x = overkt_to_C(tempkt), y = log(mean_chla_mg_m))) +
                             geom_smooth(aes(group=site_id, colour = site_id), method = 'lm', se = FALSE)+
                             geom_point(aes(fill = site_id), size = 2, shape= 22, colour = 'darkgrey') +
                             scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels) +
                             scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], labels = stream_temp_labels, guide = FALSE) +
                             scale_x_continuous(name = expression("Temperature ("*degree*"C)"), expand = c(0.03,0),
                                                limits = c(0,35),sec.axis = dup_axis(
                                                  breaks = breaks_in_C_plus,
                                                  labels = breaks_in_kt_stand,
                                                  name = expression("Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])") )) +
                             scale_y_continuous(name = expression("Chlorophyll "~italic(a)~"(mg"~m^-2*")"),limits = c(-3,6),
                                                breaks = breaks_in_logchla, labels = breaks_in_chla,expand = c(0.03,0))+
                             theme_tufte(ticks = TRUE) +
                             geom_rangeframe(sides = "tbl", colour = "black") +
                             annotate('text', x = 0, y = log(300), family = 'serif', label = 'A', hjust = 0)+
                             theme(legend.position = c(0.85, 0.23), legend.title = element_blank(),
                                   legend.spacing.y = unit(0.0003,'cm'),
                                   legend.key.size = unit(0.4,'cm')))

chla_prod_int = ggplotGrob(ggplot(na.omit(int_df), aes(x = log(chla_mean),y=log(prod_d))) +
                             geom_smooth(aes(group= site_id, colour = site_id), method = 'lm', se = FALSE)+
                             geom_point(aes(fill = site_id),shape = 22, colour = 'darkgrey', size = 1.7) + 
                             scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos]) +
                             scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos], guide = FALSE) +
                             scale_y_continuous(name = expression("Secondary Production (mg AFDM  "~m^-2~d^-1~")"), expand = c(0.03,0),
                                                breaks = breaks_in_log10prod, labels = breaks_in_prod10, position = 'right') +
                             scale_x_continuous(name = expression("Chlorophyll "~italic(a)~"(mg"~m^-2*")"), expand = c(0.03,0),
                                                limits = c(-2,4),breaks = breaks_in_logchla, labels = breaks_in_chla ) +
                             theme_tufte(ticks = TRUE) +
                             geom_rangeframe(sides = 'br', na.rm = TRUE, colour = 'black') +
                             geom_rangeframe(data = na.omit(int_df), aes(x = log(chla_mean), y = log(prod_d)), sides = 'br', colour = 'black')+
                             annotate('text', x = log(50), y = log(600), family = 'serif', label = "B", hjust = 0)+
                             theme(legend.position = "none"))

maxwidth = grid::unit.pmax(chla_temp_int$widths[2:5],
                           chla_prod_int$widths[2:5])
maxheight = grid::unit.pmax(chla_temp_int$heights[2:5],
                            chla_prod_int$widths)

# tiff(file = "./output/tempchla_prod_supp.tiff",res = 600, height = 4, width = full_fig, units = "in", compression = "lzw")
grid.draw(gridExtra::gtable_cbind(chla_temp_int, chla_prod_int))
# dev.off()

## Figure S6 ##
warming_time = ggplotGrob(ggplot(prim_prod %>% rename(DATE = 'Pdt', site_id = 'stream') %>% filter(site_id %in% c("st7","oh2")), aes(x = DATE, y = EstGPP_gCm2d_med, colour = site_id)) +
                            geom_line(size = 0.8) +
                            geom_ribbon(aes(x = DATE, ymin = EstGPP_gCm2d_L95per, ymax = EstGPP_gCm2d_U95per, colour = site_id, fill = site_id), alpha = 0.5)+
                            geom_col(data = warming_intdf %>% filter(site_id %in% c("st7","oh2")), aes(x = DATE, y = chla_mean/10, fill = site_id), width = 5, alpha = 0.6) +
                            scale_y_continuous(name = expression("GPP [g C"~m^-2~d^-1*"] | Chlorophyll"~italic(a)~"(mg"~m^-2*"/10)"), expand = c(0.003,0), limits = c(NA, 6)) +
                            scale_colour_manual(values = ocecolors[['temperature']][oce_temp_pos][4:5])+
                            scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos][4:5])+
                            theme_tufte(ticks = TRUE) +
                            geom_rangeframe(sides = "tblr", colour = "black") +
                            facet_grid(~site_id)+
                            theme(legend.position = 'none'))

# tiff("./output/warming_timeseries.tiff", res = 600, height = 5, width = full_fig, units = "in", compression = "lzw")
grid.draw(warming_time)
# dev.off()

## Figure S7
### GPP predictive model to estimate the temp dependence of within stream resource supply ###

light_gpp = ggplotGrob(ggplot(warming_df, aes(x = log(mean_light), y = log(gpp_d))) + 
                         geom_smooth(method = 'lm', se = FALSE, color = 'darkgrey')+
                         geom_point(aes(fill = site_id), shape = 21, size = 2)+
                         scale_y_continuous(name = expression(log[e]*"Gross Primary Production (g C"~m^-2~d^-1*")"), expand = c(0.03,0))+
                         scale_x_continuous(name = expression('Mean Interval Light (Lux'~d^-1*")"))+
                         scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos][4:5])+
                         theme_tufte(ticks = TRUE) +
                         geom_rangeframe(sides = "bl", colour = "black") +
                         annotate('text', x = -Inf, y = Inf, family = 'serif', label = 'A', size = 5, hjust = 0, vjust = 1)+
                         theme(legend.title = element_blank(), legend.position = c(1,0),
                               legend.justification = c(1,0),axis.title.y = element_blank(),
                               legend.background = element_rect(fill = 'transparent', colour = NA)));grid.draw(light_gpp)

temp_gpp = ggplotGrob(ggplot(warming_df, aes(x = tempkt_stand, y = log(gpp_d))) + 
                        geom_smooth(method = 'lm', se = FALSE, color = 'darkgrey')+
                        geom_point(aes(fill = site_id), shape = 21, size = 2)+
                        scale_y_continuous(name = expression(log[e]*"Gross Primary Production (g C"~m^-2~d^-1*")"), expand = c(0.03,0))+
                        scale_x_continuous(name = expression("Standardized temperature (1/["*italic(kT)[15]~"-"~italic(kT)*"])"))+
                        scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos][4:5])+
                        theme_tufte(ticks = TRUE) +
                        geom_rangeframe(sides = "bl", colour = "black") +
                        annotate('text', x = -Inf, y = Inf, family = 'serif', size = 5, label = 'B', hjust = 0, vjust = 1)+
                        theme(legend.title = element_blank(), legend.position = 'none',
                              axis.title.y = element_blank()))

pred_gpp = ggplotGrob(ggplot(na.omit(warming_df), aes(x = log(EstGPP_gCm2d_warming), y = log(gpp_d))) + 
                        geom_point(aes(fill = site_id), shape = 21, size = 2)+
                        geom_abline()+
                        scale_y_continuous(name = expression("Gross Primary Production (g C"~m^-2~d^-1*")"), expand = c(0.03,0))+
                        scale_x_continuous(name = expression(log[e]*"Predicted Gross Primary Production (g C"~m^-2~d^-1*")"), expand = c(0.03,0))+
                        scale_fill_manual(values = ocecolors[['temperature']][oce_temp_pos][4:5])+
                        theme_tufte(ticks = TRUE) +
                        geom_rangeframe(sides = "bl", colour = "black") +
                        annotate('text', x = -Inf, y = Inf, family = 'serif', label = 'D', hjust= 0, vjust = 1, size= 5)+
                        theme(legend.title = element_blank(), legend.position = 'none',
                              axis.title.y = element_blank()))


# tiff("./output/gpp_mod_pub.tiff", res = 600, height = full_fig , width = full_fig, units = "in", compression = "lzw")
grid.arrange(light_gpp, temp_gpp, pred_gpp, ncol = 2, left = textGrob(expression(log[e]*"-Observed Gross Primary Production (g C"~m^-2~d^-1*")"), rot = 90))
# dev.off()
