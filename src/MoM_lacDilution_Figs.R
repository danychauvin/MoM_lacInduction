
#####
(myfigs[[1]] <- 
  plot_grid(
    myplots[['raw_traces']](1.2),
    # second row (caption)
    plot_grid(
      # caption row
      NULL,
      get_legend(myplots[['basal_perturb_pre_violin']] +
                   theme(legend.title = element_text(size=rel(0.6)),
                         legend.text =  element_text(size=rel(0.6)),
                         legend.key.size = unit(0.6, "lines"),
                         legend.justification = c(0.9, 0),
                         legend.margin = margin(t=20))),
      NULL,
      nrow=1, rel_widths=c(1.5, 1.5, 1)
    ),
    # plots row
    plot_grid(
      myplots[['naive_lags_hist']] +
        theme(axis.title.x = element_text(margin=margin(t=-110))),
      myplots[['basal_perturb_pre_violin']] +
        labs(y=expression(paste(italic("lac"), " induction lag (min)    "))) +
        theme(legend.position = 'none',
              plot.margin = margin(t=20),
              # axis.title.y = element_text(margin=margin(b=-100)),
              axis.title.x = element_blank() ),
      myplots[['glyc_mix_violin']] +
        theme(legend.position = 'none',
              plot.margin = margin(t=20),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()),
      nrow=1, labels=c("B", "C", "D"), rel_widths=c(1.5, 1.5, 1), align='h'
    ),
    labels=c('A', ''), ncol=1, rel_heights=c(1, 0, 0.8)
  ))

save_plot(here("plots", "MoM_lacDilution_fig1.pdf"), myfigs[[1]],
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.1
)

#####
(myfigs[[2]] <- plot_grid(
  myplots[['memory_cdfs']] +
    theme(axis.title.y = element_text(margin=margin(l=-2, r=4))),
  myplots[['memory_frac_short']],
  labels="AUTO", ncol=1, rel_heights=c(1, 1), align='v'
  ))

save_plot(here("plots", "MoM_lacDilution_fig2.pdf"), myfigs[[2]],
          base_height=NULL, base_width=2.25 * 14/8, # 2 cols
          base_aspect_ratio = 1/1.5
)

#####
# myfigs[[3]] <- 
#   plot_grid(
#     myplots[['naive_arrest_frac']],
#     myplots[['lacl_arrest_induction_fracs']],
#     NULL,
#     labels="AUTO", ncol=1, rel_heights=c(0.45, 1, 1)
#   )
# 
# save_plot(here("plots", "MoM_lacDilution_fig3.pdf"), myfigs[[3]],
#           base_height=NULL, base_width=2.25 * 14/8, # 1 col
#           base_aspect_ratio = 1/2.2
# )

############
(myfigs[[3]] <- function() { # local envt
  pdftools_installed <- require(pdftools)
  plot_grid(
    plot_grid(myplots[['naive_arrest_cdf']] + labs(y='fraction of growth\narrested cells'), NULL,
                                                 nrow=1, labels=c('A', 'B'), rel_widths=c(1.5, 2) ),
    plot_grid(NULL, get_legend(myplots[['lacl_gr_hist']] + guides(col='legend') +
                                 scale_colour_discrete(name='nutrient', breaks=c('switch_glycerol_TMG20', 'switch_lactulose_TMG20_lowIllum'),
                                                       labels=c('glycerol', 'lactulose')) +
                                 theme_cowplot_legend_inset() +
                                 theme(legend.position = 'top', legend.justification = c(0, 0), 
                                       legend.box = 'vertical', legend.box.just = 'left', 
                                       legend.spacing = unit(0, 'mm'), legend.box.spacing = unit(0, 'mm'),
                                       legend.margin=margin(), legend.box.margin = margin(t=30),
                                 )),
              nrow=1, rel_widths=c(1, 1.42)), 
    plot_grid(
      if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "MoM_lacDilution_fig3_cartoon.pdf")), scale=1) + 
        theme(plot.margin = margin(t=4)),
      if (!pdftools_installed) NULL else ggdraw() + draw_image(here("material", "montage_TMG_glyc_lacl.jpg"), y=0.033, scale=0.77),
    myplots[['lacl_gr_hist']] + 
      theme(plot.margin = margin(t=30),
            legend.position = 'none',
            strip.background = element_blank(),
            strip.text.y = element_blank()), 
    myplots[['lacl_lacGFP_hist']] +
      theme(plot.margin = margin(t=30, l=10, r=10),
            legend.position = 'none',
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            strip.text.y = element_blank()), 
    nrow=1, labels=c('C', 'D', "E", ""), rel_widths=c(0.8, 0.25, 0.8, 1.1)),
  ncol=1, rel_heights=c(0.8, 0, 1.3))
}) ()

save_plot(here("plots", "MoM_lacDilution_fig3.pdf"), myfigs[[3]](),
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.4
)

