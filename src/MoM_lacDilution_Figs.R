
#####
myfigs[[1]] <- 
  plot_grid(
    # plot_grid(
    #   NULL, 
    #   NULL,
    #   labels=c('A', 'B'), nrow=1, rel_widths=c(2, 1)
    # ),
    myplots[['raw_traces']](),
    plot_grid(
      plot_grid(
        myplots[['naive_lags_hist']],
        myplots[['memory_cdfs']],
        labels=c('B', 'C'), ncol=1, rel_heights=c(0.7, 1), align='v'
      ),
      myplots[['memory_frac_short']],
      labels=c('', 'D'), nrow=1, rel_widths=c(1.2, 1)
    ),
    labels=c('A', ''), ncol=1, rel_heights=c(1, 1.1)
  )

save_plot(here("plots", "MoM_lacDilution_fig1.pdf"), myfigs[[1]],
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.1
)

#####
(myfigs[[2]] <- (function(.relh=c(1.4, 1.1), .relw=c(4.5, 3), .prop_breaks=c(0, 0.5, 1))
  plot_grid(
    # first row
    plot_grid(
      get_legend(myplots[['basal_perturb_pre_violin']] +
                   theme(legend.justification = c(0.5, 0))),
      NULL,
      labels="AUTO", nrow=1, rel_widths=c(.relw)),
    # second row
    plot_grid(
      plot_grid(
        myplots[['basal_perturb_pre_violin']] +
          theme(axis.text.x=element_blank(),
                legend.position = 'none',
                # legend.position = c(1, 1), legend.justification = c(1, 1),
          ),
        myplots[['basal_perturb_pre_frac_short']] +
          scale_y_continuous(breaks=.prop_breaks) +
          labs(y='prop. of\nshort lags') +
          guides(fill='none') +
          NULL,
        ncol=1, rel_heights = .relh, align='v'),
      plot_grid(
        # NULL,
        myplots[['glyc_mix_violin']] +
          theme(axis.text.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position = 'none'
          ),
        myplots[['glyc_mix_frac_short']] +
          scale_y_continuous(breaks=.prop_breaks) +
          guides(fill='none') +
          theme(axis.title.y=element_blank()) +
          NULL,
        ncol=1, rel_heights = .relh, align='v'),
    nrow=1, rel_widths=.relw, align='h'),
    ncol=1, rel_heights=c(0.15, sum(.relh))
  ))())

save_plot(here("plots", "MoM_lacDilution_fig2.pdf"), myfigs[[2]],
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.6
)

#####
myfigs[[3]] <- 
  plot_grid(
    myplots[['naive_arrest_frac']],
    myplots[['lacl_arrest_induction_fracs']],
    NULL,
    labels="AUTO", ncol=1, rel_heights=c(0.45, 1, 1)
  )

save_plot(here("plots", "MoM_lacDilution_fig3.pdf"), myfigs[[3]],
          base_height=NULL, base_width=2.25 * 14/8, # 1 col
          base_aspect_ratio = 1/2.2
)

############
(pl1 <- myframes_switching %>% 
  ungroup() %>% 
  # filter(between(time_sec/3600, 6, 15.99), !discard_top) %>% 
  filter(condition %in% c('switch_glycerol_TMG20', 'switch_lactulose_TMG20_lowIllum')) %>%
  filter(switch_idx==1) %>% 
  filter(time_sec>=t_lac_switch, time_sec<=t_lac_switch+3*360+10) %>% 
  group_by(condition, date, pos, gl, id) %>%
  # summarise(n=n()) %>% qplot(n, data=.)
  filter(n()==4) %>% 
  do(mod_ll_t=lm(log(length_um)~time_sec, data=.)) %>% 
  mutate(islope_ll_time=mod_ll_t$coefficients[2], 
         islope_ll_time_sd=summary(mod_ll_t)$coefficients[2,2],
         # ir2_ll_time=summary(mod_ll_t)$r.squared,
  ) %>% 
   mutate(condition=fct_relevel(condition, 'switch_glycerol_TMG20')) %>% 
   # ggplot(aes(rate_sd, rate)) +
   # geom_point(stroke=0, alpha=.1)
   ggplot(aes(islope_ll_time/log(2)*3600, y=..density.., col=condition)) +
   facet_grid(condition~., scales='free_y') +
   geom_freqpoly(aes(lty='after 7h'), binwidth=.08,
                 data= myrates_lactulose %>% ungroup() %>% 
                   filter(between(time_sec/3600, 14, 15.99)) %>%
                   # filter(time_sec/3600==8+7) %>%
                   # filter(irate_npoints>3, ir2_ll_time>0.8) %>% 
                   left_join(myframes %>% select(condition, date, pos, gl, id, time_sec, length_um, gfp_nb)) %>% 
                   mutate(condition=fct_relevel(condition, 'switch_glycerol_TMG20')) %>% 
                   filter((condition=='switch_lactulose_TMG20_lowIllum' & gfp_nb/length_um > 250) | 
                            condition=='switch_glycerol_TMG20') ) +
   geom_freqpoly(aes(lty='at the switch'), binwidth=.08) +
   coord_cartesian(xlim=c(-0.2, 1.1)) +
   # scale_colour_discrete(labels=rename_conds, guide='none') +
   labs(x='inst. growth rate (dbl/h)', lty='sample') +
   guides(col='none') +
   theme(legend.position='top') +
   NULL
)

(pl2 <- myframes %>% select(condition, date, pos, gl, id, time_sec, length_um, gfp_nb) %>% 
    ungroup() %>% 
    filter(condition %in% c('switch_glycerol_TMG20', 'switch_lactulose_TMG20_lowIllum')) %>%
    filter(between(time_sec/3600, 8+7, 8+8)) %>%
    mutate(condition=fct_relevel(condition, 'switch_glycerol_TMG20')) %>% 
    ggplot(aes(gfp_nb/length_um, ..density..)) +
    facet_grid(condition~., scales='free_y') +
    geom_freqpoly(aes(col=condition), binwidth=50) +
    geom_vline(aes(xintercept=225), data.frame(condition='switch_lactulose_TMG20_lowIllum'), lty='dashed') +
    # scale_y_continuous(labels=scales::scientific) +
    # coord_cartesian(xlim=c(-0.2, 1)) +
    labs(x='[LacZ-GFP] (molecules/µm)') +
    guides(col='none') +
    theme(legend.position='top') +
    NULL
)

(tmp <- plot_grid(
  plot_grid(NULL, get_legend(pl1 + theme(legend.position = 'top', legend.justification = c(0.5, 0))),
            nrow=1, rel_widths=c(0.45, 2)), 
  plot_grid(  
    NULL, 
    pl1 + 
      theme(legend.position = 'none',
            strip.background = element_blank(),
            strip.text.y = element_blank()), 
    pl2 +
      theme(legend.position = 'none',
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_blank()), 
    nrow=1, rel_widths=c(0.45, 1, 1)),
  ncol=1, rel_heights=c(0.1, 1))
)
save_plot(here("plots", "MoM_lacDilution_fig3tmp.pdf"), tmp,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 2
)
