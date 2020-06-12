
#  xxx add comments + flow control
# remove lac_iptg
mytables[['expts_list']] %>% 
  filter_article_ds() %>% 
  ungroup() %>%
  mutate(
    condition=fct_relevel(factor(condition), 'mg1655', 'glucose', 'lactose'),
    condition=fct_relevel(factor(condition), 'switch_late', 'switch_ramp40min', 'switch_preIPTG5uM',
                          'switch_lacIoe', 'switch_lacIoe_preIPTG10uM', 
                          'switch_gly_lac', 'switch_glcLac_lac', after=200L),
    label=str_replace(label, '>', ' to '),
  ) %>%
  arrange(condition) %>% 
  select(-condition) %>% rename(condition=label) %>% 
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


myplots[['naive_lags_per_expt_facet']]() %>% 
  save_plot(here("plots", "SI_figs", "naive-lags-per-expt.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.2)


myplots[['naive_lags_per_pos']]() %>% 
  save_plot(here("plots", "SI_figs", "naive-lags-per-pos.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.6)

myplots[['naive_arrest_hist']] %>% 
  save_plot(here("plots", "SI_figs", "growth-lags-histo.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2.1)

plot_grid(
  plot_grid(
    myplots[['naive_lags_per_expt']](),
    myplots[['naive_arrest_hist']] +
      coord_cartesian(xlim=c(-10, 240), expand=FALSE) +
      NULL,
    ncol=1, labels = c('A', 'C'), align = 'v'),
  myplots[['naive_lags_per_pos']](c(1, 1)),
  nrow = 1, labels=c('', 'B'), align='hv', rel_widths = c(0.45, 0.6)) %>% 
  save_plot(here("plots", "SI_figs", "sc-lags-ctrls.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 1.6)


myplots[['lags_hist_ramp']] %>% 
  save_plot(here("plots", "SI_figs", "lags-hist-ramp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2.1)

(myplots[['glyc_mix_violin']] +
  theme(legend.position = 'none')) %>% 
  save_plot(here("plots", "SI_figs", "glyc-mix-violin.pdf"), .,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 2)

(myplots[['basal_perturb_gfp']] +
  labs(x='', y='LacZ-GFP at the switch\n(molecules / cell)') )%>% 
  save_plot(here("plots", "SI_figs", "basal-perturb-gfp.pdf"), .,
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


(myplots[['lags_inherited_gfp']] +
  scale_x_log10(labels=function(.x) formatC(.x, format="fg")) ) %>% 
  save_plot(here("plots", "SI_figs", "lags-inherited-gfp.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2.3)


plot_grid(
  plot_grid(
    myplots[['flim_gfp_criteria_cpm']],
    
    ggdraw() + 
      draw_text("Criteria 1:\nnorm CPM > mean + 5 s.d. of background", size=10, x=-0.5, y=.8, hjust=0, vjust=1),
          
    myplots[['flim_gfp_criteria_nb']] +
      theme(axis.title.x = element_blank()),
    
    ggdraw() + 
      draw_text("Criteria 2:\nfitted number of molecules < 20", size=10, x=-0.5, y=.8, hjust=0, vjust=1),
    
    myplots[['flim_gfp_criteria_diff']] +
      theme(axis.title.x = element_blank()),
    
    ggdraw() + 
      draw_text("Criteria 3:\n3ms < fitted diffusion time < 30ms", size=10, x=-0.5, y=.8, hjust=0, vjust=1),
    
    ncol=2, nrow=3, rel_widths = c(2,1), rel_heights = c(1.3, 1, 1), align='v'
  ),
  plot_grid(
    NULL, 
    get_legend(
      myplots[['flim_gfp_criteria_nb']] +
        labs(col='LacZ-GFP status\nbefore the switch') +
        theme(legend.position = 'bottom') +
        NULL),
    nrow=1, rel_widths = c(.125, 1)),
  ncol=1, rel_heights=c(10, 1) ) %>% 
# draw_grob(, 2/3, 2/3, 1/3, 0.5) # cf https://stackoverflow.com/a/41570754/576684
  save_plot(here("plots", "SI_figs", "flim-gfp-criteria.pdf"), .,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.4)


(myplots[['diauxie_gcs_all']] +
    # scale_color_manual(values = c('0µM'=ggCustomTJ::qual_cols[2], '200µM'=ggCustomTJ::qual_cols[1])) +
    theme_half_open() + # this is needed otherwise setting strip.text raises an error ?!?
    theme(legend.position = 'top', panel.border = element_rect(colour='gray50')) +
    # guides(colour=guide_legend(title.position="top")) +
    NULL) %>% 
  save_plot(here("plots", "SI_figs", "diauxie-gcs-all.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = .75)

(myplots[['diauxie_wash_ctrl']] +
  ggCustomTJ::scale_colour_discrete() +
  theme_half_open() + # this is needed otherwise setting strip.text raises an error ?!?
  theme(legend.position = 'top',
        strip.background = element_blank(), strip.text = element_blank()) +
  guides(colour=guide_legend(title.position="top")) +
  NULL) %>% 
  save_plot(here("plots", "SI_figs", "diauxie-wash-ctrl.pdf"), .,
            base_height=NULL, base_width=2.25 * 14/8, # 1 col
            base_aspect_ratio = 1.15)

(myplots[['2cs_qms']] +
    # scale_x_continuous(trans='log10', breaks=c(1, 1e2, 1e4), labels=c('1', '1e2', '1e4')) +
    scale_x_log10(breaks = c(1, 1e2, 1e4), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
    theme(legend.position = 'right') +
    NULL ) %>% 
  save_plot(here("plots", "SI_figs", "2cs-qms.pdf"), .,
            base_height=NULL, base_width=4.75 * 14/8, # 2 cols
            base_aspect_ratio = 2.3)

