
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
      myplots[['lags_types_correl']] +
        theme_cowplot_legend_inset(0.8) +
        guides(size=guide_legend(direction="horizontal", title.position = "top")) +
        theme(
          legend.title = element_text(size=rel(0.6)),
          legend.text =  element_text(size=rel(0.6)),
          legend.key.size = unit(0.6, "lines")
          ),
      nrow=1, labels=c("A", "C", "D"), rel_widths=c(1.2, 1.3, 1.1), align='h'
    ),
   plot_grid(
     myplots[['raw_traces']](1.25),
     myplots[['basal_perturb_pre_violin']](-35) +
       coord_cartesian(ylim=c(-50, 240)) +
       guides(col= guide_legend(direction="horizontal", title.position = "top"),
              fill=guide_legend(direction="horizontal", title.position = "top")) +
     # labs(y=expression(paste(italic("lac"), " induction\nlag (min)"))) +
       # labs(y="lac induction\nlag (min)") +
       theme(legend.box.margin = margin(0, 0, 0, -10),
             plot.margin = margin(7, 7, -10, 7),
             legend.title = element_text(size=rel(0.6)),
             legend.text =  element_text(size=rel(0.6)),
             legend.key.size = unit(0.6, "lines"),
             axis.title.x = element_blank(),
       ),
     labels=c('B', 'E'), nrow=1, rel_widths=c(2.5, 1.1)
   ),
   ncol=1, rel_heights=c(0.7, 1), align='v'
  )
})()

save_plot(here("plots", "figs", "MoM_lacInduction_fig1.pdf"), myfigs[[1]](),
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_asp = 1.5
)

#####
(myfigs[[2]] <- plot_grid(
  myplots[['lags_gfp_diff_cdf']] +
    coord_cartesian(xlim=c(-40, 60)) +
    scale_x_continuous(breaks=c(-40, 0, 40, 80)) +
    # scale_y_continuous(breaks=c(0, 0.5, 1)) +
    theme_cowplot_legend_inset(0.8) +
    labs(x='normalized fluo. at the switch\n(equiv. GFP molecules)',
         y='rev. cumul.\nprobability') +
    theme(legend.position = c(1,1),
          legend.justification = c(1,1)) +
    NULL,
  myplots[['memory_cdfs']] +
    labs(y='rev. cumul.\nprobability') +
    theme(axis.title.y = element_text(margin=margin(l=-2, r=4))),
  myplots[['memory_frac_short']],
  labels="AUTO", ncol=1, rel_heights=c(1.05, 1, 1), axis='l', align='v'
))

save_plot(here("plots", "figs", "MoM_lacInduction_fig2.pdf"), myfigs[[2]],
          base_height=NULL, base_width=2.25 * 14/8, # 1 col
          base_aspect_ratio = 1/1.7
)

############
(myfigs[[3]] <- function() { # local envt
  # browser()
  # pdftools_installed <- require(pdftools)
    plot_grid(
      NULL,
      
      myplots[['glyc_mix_violin']] +
        theme(axis.title.x = element_blank()),
      
      myplots[['simul_lags']] +
        scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
        labs(col='condition') +
        # labs(x="fraction of short\ngrowth lags") +
        scale_fill_manual(values=c(
          'gluc > lac (naive)'=ggCustomTJ::qual_cols[2], 'gluc > lac (full memory)'=ggCustomTJ::qual_cols[4], 
          'gluc + lac > lac'=ggCustomTJ::qual_cols[1],  'glyc > lac'=ggCustomTJ::qual_cols[3],  
          'gluc > lac (short only)'=ggCustomTJ::qual_cols[5], 'gluc > lac (long only)'=ggCustomTJ::qual_cols[7])) +
        guides(col=guide_legend(ncol = 2)) +
        theme(
          # legend.position = 'right',
          legend.position = 'bottom', legend.title=element_blank(),
          # legend.box.spacing = unit(0.25, "lines"), #plot.margin=margin(14, 7, 7, 7, "pt"),
          legend.text = element_text(size=rel(0.8)), legend.key.size = unit(0.6, "lines"),
        ) +
        
        expand_limits(y=-.1) +
        annotation_custom(ggplotGrob(
          myplots[['simul_gcs']] +
            theme_cowplot(font_size = 11) +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
            scale_y_continuous(trans='log10', breaks=c(2e8, 4e8, 8e8)) +
            geom_polygon(aes(fill=type), x=0, y=0) + # hack to change the legend appearance
            theme(legend.position = 'none') +
            guides(col=guide_legend(ncol = 1)) +
            labs(y='pop. size') +
            NULL), 
          xmin=-.04, xmax=0.5, ymin=-.18, ymax=0.8) +
        NULL
      , 
      
      nrow=3, labels="AUTO", rel_heights = c(0.9, 1, 1)) 
}) ()

save_plot(here("plots", "figs", "MoM_lacInduction_fig3.pdf"), myfigs[[3]](),
          base_height=NULL, base_width=2.25 * 14/8, # 1 col
          base_aspect_ratio = 1/2.1
)



(myfigs[[4]] <- function() { # local envt
  diauxie_env <- new.env()
  load('material/SC1ss_diauxieGC_plots.RData', envir=diauxie_env)
  # browser()
  
  plot_grid(
    plot_grid(
      diauxie_env$myplots[['diauxie_gcs']] +
        xlim(NA, 10) +
        guides(alpha='none') +
        # theme(legend.position=c(0.02, 0.98), legend.justification=c(0, 1)) +
        theme(legend.position = "top") +
        annotation_custom(ggplotGrob(
          diauxie_env$myplots[['diauxie_delay']] +
            theme_cowplot(font_size = 11) +
            annotate("segment", x=.08, xend=.18, y=-3, yend=-3, lty='dotted') +
            annotate("segment", x=.08, xend=.18, y=49, yend=49, lty='dotted') +
            annotate("segment", x=.2, xend=.2, y=-3, yend=49,
                     arrow=arrow(ends="both", type='closed', length = unit(0.02, "inches"))) +
            annotate('text', .12, 25, label='lag', hjust=.5, vjust=0, angle=90) +
            scale_x_log10(limits=c(8e-3, .2)) +
            scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
            NULL), 
          xmin=3, xmax=10.5, ymin=log10(.0055), ymax=log10(0.04)) +
        NULL,
      diauxie_env$myplots[['diauxie_lags']] +
        labs(y='population growth lag (min)') +
        theme(
          #legend.position = c(.01, .99), legend.justification = c(0,1),
          legend.position = 'top', 
          axis.text.x=element_text(angle=45, hjust=1, vjust=1)),
      ncol=2, labels=c("A", "B"), rel_widths = c(1,1), align="h"),
    
    myplots[['2cs_qms_kinases']] +
      coord_cartesian(xlim=c(0.08, 500), ylim=c(0, 1.8)) +
      scale_x_continuous(labels = identity, trans = 'log10') +
      NULL,
    nrow=2, labels=c("", "C"), rel_heights = c(1, 0.8))
}) ()

save_plot(here("plots", "figs", "MoM_lacInduction_fig4.pdf"), myfigs[[4]](),
          base_height=NULL, base_width=4 * 14/8, # 1 col
          base_aspect_ratio = 1.25
)

