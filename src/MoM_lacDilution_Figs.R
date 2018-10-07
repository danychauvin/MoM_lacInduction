
theme_set(theme_cowplot())

#####
fig1 <- 
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

save_plot(here("plots", "MoM_lacDilution_fig1.pdf"), fig1,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.1
)

#####
(fig2 <- (function(.relh=c(1.4, 1.1), .relw=c(4.5, 3), .prop_breaks=c(0, 0.5, 1))
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

save_plot(here("plots", "MoM_lacDilution_fig2.pdf"), fig2,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.6
)

#####
fig3 <- 
  plot_grid(
    myplots[['naive_arrest_frac']],
    myplots[['lacl_arrest_induction_fracs']],
    NULL,
    labels="AUTO", ncol=1, rel_heights=c(0.45, 1, 1)
  )

save_plot(here("plots", "MoM_lacDilution_fig3.pdf"), fig3,
          base_height=NULL, base_width=2.25 * 14/8, # 1 col
          base_aspect_ratio = 1/2.2
)


