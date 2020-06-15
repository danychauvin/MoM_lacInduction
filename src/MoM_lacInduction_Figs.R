
#####
(myfigs[[1]] <- function() {
  pdftools_installed <- try(require(pdftools))
  plot_grid(
    # plots row
    plot_grid(
      if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "lacOperon_lactose.ai.pdf")), 
                                                               x=.0, y=.07, scale=.94) + # x=-.15, y=.12, scale=1.4) + 
        theme(),
      myplots[['naive_lags_hist']] +
        # theme(plot.margin = margin(28, 7, 26, 7)) +
        NULL,
      # myplots[['glyc_mix_violin']] +
      #   theme(legend.position = 'none',
      #         plot.margin = margin(t=24),
      #         axis.title.x = element_blank(),),
      myplots[['lags_types_correl']](2.5) +
        theme_cowplot_legend_inset(0.8) +
        guides(size=guide_legend(direction="horizontal", title.position = "top")) +
        theme(
          legend.title = element_text(size=rel(0.8)),
          legend.text =  element_text(size=rel(0.7)),
          legend.key.size = unit(0.7, "lines")
        ),
      nrow=1, labels=c("A", "C", "D"), rel_widths=c(0.95, 1.05, 1), align='h'
    ),
    plot_grid(
      myplots[['raw_traces']](1.22),
      myplots[['basal_perturb_pre_violin']](-35, 3) +
        coord_cartesian(ylim=c(-50, 240)) +
        guides(col= guide_legend(direction="horizontal", title.position = "top"),
               fill=guide_legend(direction="horizontal", title.position = "top")) +
        # labs(y=expression(paste(italic("lac"), " induction\nlag (min)"))) +
        # labs(y="lac induction\nlag (min)") +
        theme(legend.box.margin = margin(0, 0, 0, -10),
              plot.margin = margin(7, 7, -10, 7),
              legend.title = element_text(size=rel(0.8)),
              legend.text =  element_text(size=rel(0.7)),
              legend.key.size = unit(0.7, "lines"),
              axis.title.x = element_blank(),
        ),
      labels=c('B', 'E'), nrow=1, rel_widths=c(2, 1)
    ),
    ncol=1, rel_heights=c(0.7, 1), align='v'
  )
})()

save_plot(here("plots", "figs", "MoM_lacInduction_fig1.pdf"), myfigs[[1]](),
          base_height=NULL, base_width=7.5 * 12/8, # full width
          base_asp = 2
)

#####
(myfigs[[2]] <- plot_grid(
  myplots[['lags_gfp_diff_cdf']] +
    coord_cartesian(xlim=c(-40, 60)) +
    scale_x_continuous(breaks=c(-40, 0, 40, 80)) +
    # scale_y_continuous(breaks=c(0, 0.5, 1)) +
    theme_cowplot_legend_inset(0.8) +
    labs(x='normalized fluo. at the switch\n(equiv. GFP molecules)',
    # y='reverse cumulative probability') +
    y='rev. cumul. probability') +
    theme(legend.position = c(1,1),
          legend.justification = c(1,1)) +
    NULL,
  
  myplots[['memory_cdfs']] +
    labs(y='rev. cumul. probability') +
    # labs(y='reverse cumulative probability') +
    theme(axis.title.y = element_text(margin=margin(l=-2, r=4))),
  
  myplots[['memory_frac_short']],
  
  # labels="AUTO", ncol=1, rel_heights=c(1.05, 1, 1), axis='l', align='v'
  labels="AUTO", nrow=1, rel_widths=c(.9, 1.2, 0.9), axis='b', align='h'
))

save_plot(here("plots", "figs", "MoM_lacInduction_fig2.pdf"), myfigs[[2]],
          base_height=NULL, base_width=7.5 * 12/8, # full width
          base_aspect_ratio = 4
)

############
(myfigs[[3]] <- plot_grid(
  myplots[['FLIM_decay']] +
    theme_cowplot_legend_inset() +
    theme(axis.title.x = element_text(margin = margin(-12, 0, 0, 0, "pt"))) +
    NULL,
  
  myplots[['FLIM_snapshot_hist']] +
    scale_y_continuous(breaks=c(0, 0.15, .3, .45)) +
    NULL,
  
  myplots[['FLIM_lag_hist']] +
    scale_y_continuous(breaks=scales::pretty_breaks(n=2) ) +
    theme_cowplot_legend_inset() +
    theme(legend.position = c(.99, .99), legend.justification = c(1, 1)) +
    NULL,

    labels='AUTO', nrow=1, rel_widths=c(1, 0.75, 1.25), align='h'
))
  
save_plot(here("plots", "figs", "MoM_lacInduction_fig3.pdf"), myfigs[[3]],
          base_height=NULL, base_width=7.5 * 12/8, # 2 cols
          base_aspect_ratio = 4
)


############
(myfigs[[4]] <- function() { # local envt
  # browser()
  # pdftools_installed <- require(pdftools)
    plot_grid(
      plot_grid(
      NULL,
      
      myplots[['glyc_mix_violin']] +
        theme(axis.title.x = element_blank()),
     
       nrow=2, labels=c("A", "B"), rel_heights = c(0.9, 1)),

      plot_grid(
        myplots[['simul_gcs']] +
          # scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
          scale_y_continuous(trans='log10', breaks=c(2e8, 4e8, 8e8)) +
          theme(legend.position = 'none') +
          guides(col=guide_legend(ncol = 1)) +
          # labs(y='pop. size') +
          NULL,
        
        myplots[['simul_lags']] +
          scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
          scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
          labs(col='condition') +
          # labs(x="fraction of short\ngrowth lags") +
          scale_fill_manual(values=c(
            'gluc > lac (naive)'=ggCustomTJ::qual_cols[2], 'gluc > lac (full memory)'=ggCustomTJ::qual_cols[4], 
            'gluc + lac > lac'=ggCustomTJ::qual_cols[1],  'glyc > lac'=ggCustomTJ::qual_cols[3],  
            'gluc > lac (short only)'=ggCustomTJ::qual_cols[5], 'gluc > lac (long only)'=ggCustomTJ::qual_cols[7]),
            labels=c('gluc \u2794 lac (naive)', 'gluc \u2794 lac (full memory)', 
                     'gluc + lac \u2794 lac',  'glyc \u2794 lac',  'gluc \u2794 lac (short only)', 
                     'gluc \u2794 lac (long only)')) +
          guides(col=guide_legend(ncol = 2)) +
          theme(
            # legend.position = 'right',
            legend.position = 'bottom', legend.title=element_blank(),
            # legend.box.spacing = unit(0.25, "lines"), #plot.margin=margin(14, 7, 7, 7, "pt"),
            legend.text = element_text(size=rel(0.8)), legend.key.size = unit(0.6, "lines"),
          ) +
          
          expand_limits(y=-.1) +
          NULL, 
        
        nrow=2, labels=c("C", "D"), rel_heights = c(0.9, 1), align='v'),
      
      ncol=2)
}) ()

save_plot(here("plots", "figs", "MoM_lacInduction_fig4.pdf"), myfigs[[4]](),
          device=grDevices::cairo_pdf, base_height=NULL, base_width=5.2 * 12/8, # 1 col
          base_aspect_ratio = 1.5
)



(myfigs[[5]] <- plot_grid(
  myplots[['diauxie_gcs']] +
    xlim(NA, 8) +
    guides(alpha='none') +
    # theme(legend.position=c(0.02, 0.98), legend.justification=c(0, 1)) +
    theme(legend.position = "top") +
    annotation_custom(ggplotGrob(
      myplots[['diauxie_delay']] +
        theme_cowplot(font_size = 11) +
        annotate("segment", x=.08, xend=.18, y=-3, yend=-3, lty='dotted') +
        annotate("segment", x=.08, xend=.18, y=49, yend=49, lty='dotted') +
        annotate("segment", x=.2, xend=.2, y=-3, yend=49,
                 arrow=arrow(ends="both", type='closed', length = unit(0.02, "inches"))) +
        annotate('text', .12, 25, label='lag', hjust=.5, vjust=0, angle=90) +
        scale_x_log10(name='optical density', limits=c(8e-3, .2)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
        NULL), 
      xmin=4, xmax=8.5, ymin=log10(.0035), ymax=log10(0.05)) +
    NULL,
  myplots[['diauxie_lags']] +
    labs(y='population growth lag (min)') +
    theme(
      #legend.position = c(.01, .99), legend.justification = c(0,1),
      legend.position = 'top', 
      axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    NULL,
  ncol=2, labels=c("A", "B"), rel_widths = c(1.5,1), align="h")
)

save_plot(here("plots", "figs", "MoM_lacInduction_fig5.pdf"), myfigs[[5]],
          base_height=NULL, base_width=5.2 * 12/8, # 1 col
          base_aspect_ratio = 2.8
)


(myfigs[[6]] <- myplots[['2cs_qms_kinases']] +
    labs(y='rev. cumul. probability\n(all conditions)') +
  coord_cartesian(xlim=c(0.08, 500), ylim=c(0, 1.8)) +
  scale_x_continuous(labels = identity, trans = 'log10') +
  NULL)

save_plot(here("plots", "figs", "MoM_lacInduction_fig6.pdf"), myfigs[[6]],
          base_height=NULL, base_width=5.2 * 12/8, # 1 col
          base_aspect_ratio = 3
)
