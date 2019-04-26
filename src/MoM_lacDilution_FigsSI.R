
#  xxx add comments + flow control
# remove lac_iptg
mytables[['expts_list']] %>% 
  filter(! condition %in% c("switch_ramp15min", "switch_lactose_priming"),
         ! str_detect(condition, "_stdIllum")) %>% 
  ungroup() %>% select(-condition) %>% rename(condition=label) %>% 
  mutate(condition=str_replace(condition, '>', ' to ')) %>% 
  (function(.df)
    knitr::kable(.df, "latex", booktabs=TRUE, #longtable = TRUE, 
                 col.names=c('condition', 'date', '# growth channels', '# full cell cycles', '# observations',
                             '# cells at switch', '# estimated lags', '# arrested cells at switch'),
                 caption='List of experiments used in this study with summary statistics. Experiments that have been discarded from further analysis are greyed out.') %>% 
     # kableExtra::kable_styling(full_width=TRUE) %>% 
     kableExtra::kable_styling(latex_options = c("striped", "scale_down")) %>%
     kableExtra::column_spec(3:6, width = "3.5em") %>% 
     kableExtra::column_spec(7:8, width = "5em") %>% 
     kableExtra::row_spec(which(.df$date %in% discarded_dates), italic=T, color="gray") %>% 
     identity()
  ) %>% 
  str_replace(fixed("{tab:}"), "{tab:expts-list}") %>% 
  str_replace(fixed("\\resizebox{\\linewidth}{!}"), "\\resizebox*{!}{0.9\\textheight}") %>% 
  write(here('plots', 'SI_figs', 'expts-list.tex'))



myplots[['gr_length_medians']](.article_ds=TRUE) %>% 
  save_plot(here("plots", "SI_figs", "gr-length-medians.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.35)

# plot_grid(
#   myplots[['lags_gr_before']](.size_n=FALSE) +
#     theme(legend.position=c(1,1), legend.justification=c(1,1),
#           axis.title.x=element_blank()) +
#     NULL,
#   myplots[['gr_hist_before']] +
#     theme(legend.position='none') +
#     NULL,
#   ncol=1, rel_heights=c(1.5, 1), align='v'
# )
# # gr-before-filtering


myplots[['naive_lags_per_expt']]() %>% 
  save_plot(here("plots", "SI_figs", "naive-lags-per-expt.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.2)


myplots[['naive_lags_per_pos']] %>% 
  save_plot(here("plots", "SI_figs", "naive-lags-per-pos.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.6)

myplots[['lags_hist_ramp']] %>% 
  save_plot(here("plots", "SI_figs", "lags-hist-ramp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)

myplots[['naive_lags_correl']] %>% 
save_plot(here("plots", "SI_figs", "naive-lags-correl.pdf"), .,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 2)

# (myplots[['high-illum-sensitivity']] <-
# bind_rows(
#   read.csv("../data/20180403_glu_lac_hiIllum_1_pos2_t215_MG1655_hist.csv") %>%
#     with(rep.int(value, count)) %>% tibble(strain='MG1655', int=.),
#   read.csv("../data/20180403_glu_lac_hiIllum_1_pos2_t215_ASC662_hist.csv") %>%
#     with(rep.int(value, count)) %>% tibble(strain='ASC662', int=.)
#   ) %>%
#   ggplot(aes(int, ..density.., col=strain)) +
#   geom_freqpoly(binwidth=50) +
#   scale_colour_discrete(breaks=c('MG1655', 'ASC662'), labels=c('MG1655 (70 cells)', 'ASC662 (93 cells)')) +
#   labs(x='pixel intensity (AU)') +
#   NULL)

plot_grid(
  plot_grid(
    myplots[['lags_gfp_scatter']](.fit=FALSE) +
      theme(legend.position = 'none',
            axis.title.x = element_blank()),
    myplots[['lags_gfp_diff_cdf']] +
      scale_y_continuous(breaks=c(0, 0.5, 1)) +
      theme(legend.position = 'none'),
    ncol=1, rel_heights = c(1.2, 0.8), align='v'),
  get_legend(myplots[['lags_gfp_diff_cdf']]) ,
  nrow=1, rel_widths = c(0.85, 0.15)
) %>% 
  save_plot(here("plots", "SI_figs", "lag-gfp-ini.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.35)


# myplots[['basal_perturb_with']]
# # basal-perturb-with


myplots[['basal_perturb_gfp']] %>% 
  save_plot(here("plots", "SI_figs", "basal-perturb-gfp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)

# myplots[['lags_hist_glyc_glcLac']]
# # lags-hist-glyc-glcLac


# myplots[['memory_cdfs_facets']] %>% 
#   save_plot(here("plots", "SI_figs", "memory-cdfs-facets.pdf"), .,
#             base_height=NULL, base_width=4.75 * 14/8, # 2 cols
#             base_aspect_ratio = 1.35)


myplots[['lag_memory']] %>% 
  save_plot(here("plots", "SI_figs", "lag-memory.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.2)


myplots[['memory_elapsed_divs']] %>% 
  save_plot(here("plots", "SI_figs", "memory-elapsed-divs.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)


myplots[['lags_inherited_gfp']] %>% 
  save_plot(here("plots", "SI_figs", "lags-inherited-gfp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2.3)


(myplots[['TMG_induction_gly04']] <- (function() {
  load('data/20180703_ASC662_M9gly04pc_TMG.RData')
  mydata %>% 
    filter(Ch==2) %>% 
    ggplot(aes(tmg, gfp)) +
    geom_point(alpha=.2, stroke=0, position=position_jitter(width=.1)) +
    stat_function(fun=function(c, ...) log(hill_fn(c, ...)) / log(2), args=as.list(summary(myfit_all)$coefficients[,1]), col="red") +
    stat_function(fun=function(c, ...) log(hill_fn(c, ...)) / log(2), args=as.list(summary(myfit_all)$coefficients[,1]), col="red", 
                  geom="point", n=1, xlim=log2(c(20, 20)), size=3) +
    scale_x_continuous(trans='log2', limits=c(5, 200), breaks=c(6, 12, 25, 50, 100, 200)) +
    scale_y_continuous(trans='log2', breaks=c(16, 128, 1024)) +
    labs(x="TMG concentration (µM)", y="LacZ-GFP concentration per cell (AU)") +
    NULL
  })()) %>% 
  save_plot(here("plots", "SI_figs", "induction-tmg-gly.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)


myplots[['TMG_switch_gr_hist']] %>% 
  save_plot(here("plots", "SI_figs", "growth-arrest-tmg.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)


plot_grid(
  myplots[['lacl_gfp_facets']] + theme(legend.position = 'right'),
  myplots[['lacl_gr_facets']] + theme(legend.position = 'right'),
  ncol=1, rel_heights = c(5, 4), align='v') %>% 
  save_plot(here("plots", "SI_figs", "lacl-gfp-gr-hists.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1)


myplots[['lacl_gr_lacz']] %>% 
  save_plot(here("plots", "SI_figs", "lacl-gr-lacz.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2)



