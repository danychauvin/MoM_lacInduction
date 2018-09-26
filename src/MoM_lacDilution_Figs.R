
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
fig2 <- 
  plot_grid(
    myplots[['basal_perturb_pre']],
    myplots[['glyc_frac_short']] ,
    labels="AUTO", nrow=1, rel_widths=c(2, 1), align='h'
  )

save_plot(here("plots", "MoM_lacDilution_fig2.pdf"), fig2,
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 2.2
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


