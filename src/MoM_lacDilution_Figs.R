
#####
(myfigs[[1]] <- 
   plot_grid(
     myplots[['raw_traces']](1.3),
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
       nrow=1, rel_widths=c(2, 2)
     ),
     # plots row
     plot_grid(
       myplots[['naive_lags_hist']] +
         theme(axis.title.x = element_text(margin=margin(t=-110))),
       myplots[['glyc_mix_violin']] +
         theme(legend.position = 'none',
               plot.margin = margin(t=24),
               axis.title.x = element_blank(),),
       myplots[['basal_perturb_pre_violin']] +
         labs(y=expression(paste(italic("lac"), " induction lag (min)"))) +
         theme(legend.position = 'none',
               plot.margin = margin(t=24),
               # axis.title.y = element_text(margin=margin(b=-100)),
               axis.title.x = element_blank(),
               axis.title.y = element_blank()),
       nrow=1, labels=c("B", "C", "D"), rel_widths=c(1.5, 1.2, 1.3), align='h'
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

(myplots[['lac_model_induction']] <-
   tibble(path=list.files(here("material"), "lac_induc.*", full.names = TRUE)) %>% 
   mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
   unnest() %>% 
   extract(path, c('curve', 'branch'), ".*/lac_induc(\\d)_(.+)\\.txt") %>% 
   ggplot() +
   geom_line(aes(exp(X1), exp(X2), lty=branch, col=curve)) +
   scale_x_log10(breaks=c(1e-4, 1e-3, 1e-2, .1), expand=c(0, 0),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
   scale_y_log10(breaks=c(10, 1e3, 1e5),
                 # breaks = trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x))) +
   expand_limits(y=c(2.5, 1e5)) +
   labs(x=expression(paste('inducer level ', italic('b'))), y='LacY expression\n(molecules)') +
   scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dashed')) +
   scale_color_manual(values=c('1'='black', '2'='grey60')) +
   theme(legend.position = 'none') +
   NULL
)

(myplots[['lac_model_phase']] <-
    tibble(path=list.files(here("material"), "lac_phase.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest() %>% 
    extract(path, c('curve'), ".*/lac_phasediagram_(.+)\\.txt") %>% 
    (function(.df)
      ggplot(.df) +
       geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                    data=bind_rows(.df, tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper'))))) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_vline(xintercept = log(2)/0.5, lty='dashed') +
    geom_vline(xintercept = log(2)/0.5*10, lty='dashed', col='grey60') +
    annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-6, 1e-3, 1), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.3, y=5e-7) +
    labs(x='doubling time (h)', y=expression(paste('inducer level ', italic('b')))) +
    # scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dashed')) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL)
)


(myplots[['TMG-Miller-induction']] <- (function() {
  load('data/TMG_sensitivity_Miller.RData')
  mylacz %>% 
    filter(r2>0.8) %>%
    left_join(myplates, by=c("date", "row", "col")) %>% 
    filter(!is.na(replicate)) %>%
    filter(!(media=='M9+0.2gly' & date==20190321),) %>% 
    left_join(mygcs %>% filter(!is.na(replicate)) %>% filter(elapsed>-2*3600) %>% 
                group_by(date, media, replicate, row) %>% filter(sum(od>0) > 2) %>% nest() %>% 
                mutate(mod_gc=map(data, ~lm(log(od)~elapsed, data=.)),
                       gr=map_dbl(mod_gc, ~coef(.)[2]), gr_se=map_dbl(mod_gc, ~summary(.)$coefficients[2,2]),
                       od_log=map_dbl(mod_gc, ~predict(., data.frame(elapsed=5*60)) ), od=exp(od_log),
                ), by=c("date", "row", "media", "replicate") ) %>% 
    mutate(activity=1000*slope / (od*3) / 2 * gr*3600/log(2),
           media=fct_relevel(media, 'M9+0.2ara+CA'),
           media=fct_recode(media, '0.2% arabinose+CA'='M9+0.2ara+CA', '0.2% arabinose'='M9+0.2ara', '0.2% glycerol'='M9+0.2gly', '0.2% mannose'='M9+0.2man', '0.2% pyruvate'='M9+0.2pyr'),
           media=fct_rev(media)) %>%
    (function(.df)
    {
      .fits <- .df %>% 
        filter(!(media=='M9+0.2ara' & tmg==1000),
               !(media=='M9+0.2gly' & tmg==1000),
               !(media=='M9+0.2man' & tmg>200),
               !(media=='M9+0.2pyr' & tmg>200)) %>% 
        mutate(tmg=ifelse(tmg==0, .1, tmg)) %>%
        group_by(media) %>% nest() %>%
        # filter(!media %in% c('M9+0.2ara', 'M9+0.2man')) %>%
        mutate(
          #k_ini = c('M9+0.2ara+CA'=250, 'M9+0.2ara'=100, 'M9+0.2gly'=10, 'M9+0.2man'=5, 'M9+0.2pyr'=5)[media],
          mod=map(data, function(.d) {#browser();
            nls(activity~hill_fn(tmg, f, k, m, b), data=.d,
                algorithm='port', control=nls.control(maxiter = 1e3, minFactor = 1/2^15), 
                start=list(f=200, k=100, m=2.5, b=.002), 
                lower=list(f=20, k=3, m=2, b=1.8e-3),
                upper=list(f=2000, k=800, m=4, b=2.2e-3),
            ) }),
        )
      # # tweak from https://stackoverflow.com/a/31688907
      # se <- function(x) sd(x)/sqrt(length(x))
      # maxf <- function(dati) {mean(dati) + se(dati)}
      # minf <- function(dati) {max(10e-8, mean(dati) - se(dati))}
      
      # browser()
      .df %>% 
        mutate(tmg=ifelse(tmg==0, 1e-1, tmg)) %>% 
        ggplot(aes(tmg, activity, col=media)) +
        # geom_point(aes(shape=replicate)) +
        
        # stat_function(aes(col='M9+0.2ara+CA'), fun=hill_fn_log, args=filter(.fits, media=='M9+0.2ara+CA') %>% with(as.list(coef(mod[[1]]))) ) +
        # stat_function(aes(col='M9+0.2ara'), fun=hill_fn_log, args=filter(.fits, media=='M9+0.2ara') %>% with(as.list(coef(mod[[1]]))) ) +
        # stat_function(aes(col='M9+0.2gly'), fun=hill_fn_log, args=filter(.fits, media=='M9+0.2gly') %>% with(as.list(coef(mod[[1]]))) ) +
        # stat_function(aes(col='M9+0.2man'), fun=hill_fn_log, args=filter(.fits, media=='M9+0.2man') %>% with(as.list(coef(mod[[1]]))) ) +
        # stat_function(aes(col='M9+0.2pyr'), fun=hill_fn_log, args=filter(.fits, media=='M9+0.2pyr') %>% with(as.list(coef(mod[[1]]))) ) +
        
        geom_line(data=.fits %>% group_by(media) %>% mutate(data=list(data.frame(tmg=10^seq(-1, 3, length.out = 101))),
                                                            data=map2(mod, data, function(.m, .d) mutate(.d, activity=predict(.m, .d) )) ) %>% unnest(data)
        ) +
        stat_summary(fun.args=list(na.rm=TRUE), size=.4) +
        scale_x_log10(breaks=c(1, 10, 100)) +
        # dirty hack!
        coord_cartesian(xlim=c(.098, 1500), expand = 0) +
        # expand_limits(x=0.02) +
        scale_y_log10(breaks=c(0.01, 0.1, 1)) +
        labs(x='TMG concentration (µM)', y=expression(paste('p'['lac'],' activity (AU)')) ) +
        NULL
    } )
})()
) +
  guides(col=guide_legend(nrow=2, byrow=TRUE)) +
  theme(legend.position = 'top')


(myplots[['TMG-Miller-conc-dt']] <- (function() {
  load('data/TMG_sensitivity_Miller.RData')
  mylacz %>% 
    filter(r2>0.8) %>%
    left_join(myplates, by=c("date", "row", "col")) %>% 
    filter(!is.na(replicate)) %>%
    filter(!(media=='M9+0.2gly' & date==20190321),) %>% 
    left_join(mygcs %>% filter(!is.na(replicate)) %>% filter(elapsed>-2*3600) %>% 
                group_by(date, media, replicate, row) %>% filter(sum(od>0) > 2) %>% nest() %>% 
                mutate(mod_gc=map(data, ~lm(log(od)~elapsed, data=.)),
                       gr=map_dbl(mod_gc, ~coef(.)[2]), gr_se=map_dbl(mod_gc, ~summary(.)$coefficients[2,2]),
                       od_log=map_dbl(mod_gc, ~predict(., data.frame(elapsed=5*60)) ), od=exp(od_log),
                ), by=c("date", "row", "media", "replicate") ) %>% 
    mutate(activity=1000*slope / (od*3) / 2 * gr*3600/log(2)) %>%
    filter(!(media=='M9+0.2ara' & tmg==1000),
           !(media=='M9+0.2gly' & tmg==1000),
           !(media=='M9+0.2man' & tmg>200),
           !(media=='M9+0.2pyr' & tmg>200)) %>% 
    mutate(tmg=ifelse(tmg==0, .1, tmg)) %>%
    group_by(media) %>% nest() %>%
    # filter(!media %in% c('M9+0.2ara', 'M9+0.2man')) %>%
    mutate(
      #k_ini = c('M9+0.2ara+CA'=250, 'M9+0.2ara'=100, 'M9+0.2gly'=10, 'M9+0.2man'=5, 'M9+0.2pyr'=5)[media],
      mod=map(data, function(.d) {#browser();
        nls(activity~hill_fn(tmg, f, k, m, b), data=.d,
            algorithm='port', control=nls.control(maxiter = 1e3, minFactor = 1/2^15), 
            start=list(f=200, k=100, m=2.5, b=.002), 
            lower=list(f=20, k=3, m=2, b=1.8e-3),
            upper=list(f=2000, k=800, m=4, b=2.2e-3),
        ) }),
    ) %>% 
    mutate(k = map_dbl(mod, ~coef(.)['k'])) %>% 
    left_join(mygcs %>% filter(!is.na(replicate)) %>% filter(elapsed>-2*3600, row=='H') %>% 
                group_by(date, media, replicate, row) %>% filter(sum(od>0) > 2) %>% nest() %>% 
                mutate(mod_gc=map(data, ~lm(log(od)~elapsed, data=.)),
                       gr=map_dbl(mod_gc, ~coef(.)[2]), gr_se=map_dbl(mod_gc, ~summary(.)$coefficients[2,2]),
                ) %>% group_by(media) %>% summarise(gr=mean(gr) *3600/log(2)) ) %>% 
    mutate(media=fct_relevel(media, 'M9+0.2ara+CA'),
           media=fct_recode(media, '0.2% arabinose+CA'='M9+0.2ara+CA', '0.2% arabinose'='M9+0.2ara', '0.2% glycerol'='M9+0.2gly', '0.2% mannose'='M9+0.2man', '0.2% pyruvate'='M9+0.2pyr'),
           media=fct_rev(media)) %>%
    ggplot(aes(1/gr, k, col=media)) +
    stat_smooth(aes(group=1), col='black', alpha=.1, method='lm') +
    geom_point(size=3) +
    # scale_x_continuous(trans='log10', breaks=c(.4, .7, 1)) +
    scale_x_continuous(trans='log10', breaks=c(1, 2, 3)) +
    scale_y_continuous(trans='log10', breaks=c(10, 50, 100)) +
    labs(x='doubling time (h)', y='critical TMG \nconcentration (µM)') +
    # theme(legend.position = c(1, 0), legend.justification = c(1.05, -0.05)) +
    theme(legend.position = 'none') +
    NULL
})()
)

############
(myfigs[[3]] <- function() { # local envt
  pdftools_installed <- require(pdftools)
  plot_grid(
    plot_grid(myplots[['naive_arrest_cdf']] + 
                coord_cartesian(xlim=c(0, 210), ylim=c(0, 1.05)) + 
                labs(x='time after the switch (min)', y='fraction of \narrested cells') +
                theme(axis.title.x = element_text(hjust=0.9)), 
              myplots[['lac_model_induction']] ,
              myplots[['lac_model_phase']] + theme(axis.title.y = element_text(hjust=15)),
              nrow=1, labels=c('A', 'B', 'C'), rel_widths=c(0.85, 1, 1), align='h' ),
    plot_grid(myplots[['TMG-Miller-induction']] +
                theme_cowplot_legend_inset() +
                theme(plot.margin = margin(t=10, b=10, l=6, r=10),
                      axis.title.y = element_text(hjust=0))+
                NULL,
              myplots[['TMG-Miller-conc-dt']] +
                theme(plot.margin = margin(t=10, b=10, l=6, r=10)) +
                NULL,
              nrow=1, labels=c('D', 'E'), rel_widths=c(2, 1), align='h'),
    plot_grid(NULL, get_legend(myplots[['lacl_gr_hist']] + guides(col='legend') +
                                 scale_colour_discrete(name='nutrient', breaks=c('switch_glycerol_TMG20', 'switch_lactulose_TMG20'),
                                                       labels=c('glycerol', 'lactulose')) +
                                 theme_cowplot_legend_inset() +
                                 theme(legend.position = 'top', legend.justification = c(0, 0), 
                                       legend.box = 'vertical', legend.box.just = 'left', 
                                       legend.spacing = unit(0, 'mm'), legend.box.spacing = unit(0, 'mm'),
                                       legend.margin=margin(), legend.box.margin = margin(t=30),
                                 )),
              nrow=1, rel_widths=c(1, 1.8)), 
    plot_grid(
      if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "MoM_lacDilution_fig3_cartoon.pdf")), scale=1.0) + 
        theme(plot.margin = margin(t=4)),
      if (!pdftools_installed) NULL else ggdraw() + draw_image(here("material", "montage_TMG_glyc_lacl.jpg"), y=-0.0, scale=0.8),
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
      nrow=1, labels=c('F', 'G', "H", ""), rel_widths=c(0.7, 0.2, 0.85, 1.2)),
    ncol=1, rel_heights=c(0.7, 0.8, 0, 1.0))
}) ()

save_plot(here("plots", "MoM_lacDilution_fig3.pdf"), myfigs[[3]](),
          base_height=NULL, base_width=4.75 * 14/8, # 2 cols
          base_aspect_ratio = 1.2
)

####
(myplots[['autoactiv_model_phase']] <-
    tibble(path=list.files(here("material"), "autoactivating_phasediagram_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest() %>% 
    extract(path, c('curve'), ".*/autoactivating_phasediagram_(.+)\\.txt") %>% 
    (function(.df)
      ggplot(.df) +
       geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                    data=bind_rows(.df, tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                   tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) ))) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    geom_hline(yintercept = 0.014269, lty='dashed') +
    # annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 10, 100), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-4, 1e-2, 1), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.35, y=c(2e-5, 7)) +
    labs(x='doubling time (h)', y='promoter activity') +
    # scale_linetype_manual(values=c('high'='solid', 'low'='solid', 'unstable'='dashed')) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)


(myplots[['twocmps_model_phase']] <-
    tibble(path=list.files(here("material"), "twocomp_phasediagram_.*", full.names = TRUE)) %>% 
    mutate(data=map(path, ~read_delim(., delim='\t', col_names = FALSE))) %>% 
    unnest() %>% 
    extract(path, c('curve'), ".*/twocomp_phasediagram_(.+)\\.txt") %>% 
    # with(range(X1))
    (function(.df)
      ggplot(.df) +
       geom_polygon(aes(exp(X1), exp(X2), fill=curve), 
                    data=bind_rows(.df, tibble(X1=Inf, X2=Inf, curve=c('lower', 'upper')),
                                   tibble(X1=-Inf, X2=Inf, curve=c('lower', 'upper')) ))) +
    geom_line(aes(exp(X1), exp(X2), col=curve)) +
    # geom_hline(yintercept = 0.014269, lty='dashed') +
    # annotate("text", x=0, y=0, label='uninduced', hjust=-0.1, vjust=-1.1) +
    annotate("text", x=Inf, y=0, label='uninduced', hjust=1.1, vjust=-1.1) +
    annotate("text", x=Inf, y=Inf, label='induced', hjust=1.1, vjust=1.5) +
    scale_x_log10(breaks=c(1, 24, 1000), expand=c(0, 0)) +
    scale_y_log10(breaks=c(1e-3, 1, 1e3), expand=c(0, 0),
                  # breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    expand_limits(x=.15, y=c(2e-6, 7)) +
    # labs(x='doubling time (h)', y=expression(paste('signal strength ', italic('s'), '/', italic('s[0]')))) +
    labs(x='doubling time (h)', y=expression(paste('signal strength ', s / s[0] ))) +
    scale_fill_manual(values=qual_cols %>% hex_lighten(1.2) %>% hex_desaturate(.3)) +
    theme(legend.position = 'none') +
    NULL
)

# (myfigs[[4]] <- function(){
#   pdftools_installed <- require(pdftools)
#   plot_grid(
#     if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=1), scale=1.0),
#     myplots[['autoactiv_model_phase']] +
#       # theme(axis.title.y = element_text(margin=margin(l=-2, r=4))) +
#       NULL,
#     if (!pdftools_installed) NULL else ggdraw() + draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=3), scale=1.0),
#     myplots[['twocmps_model_phase']],
#     labels="AUTO", ncol=1, rel_heights=c(1, 1, 1), align='v'
#   )}
# )()
# 
# save_plot(here("plots", "MoM_lacDilution_fig4.pdf"), myfigs[[4]](),
#           base_height=NULL, base_width=2.25 * 14/8, # 2 cols
#           base_aspect_ratio = 1/2
# )

# (myfigs[[4]] <- function(){
#   pdftools_installed <- require(pdftools)
#   plot_grid(
#     ggdraw() +
#       ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=1), 
#                                                        x=-.105, y=-.1, scale=.4) ) +
#       draw_plot(myplots[['autoactiv_model_phase']]) +
#       NULL,
#     ggdraw() +
#       ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=3), 
#                                                       x=-.105, y=-.1, scale=.4) ) +
#       draw_plot(myplots[['twocmps_model_phase']]) +
#       NULL,
#     labels="AUTO", ncol=1, rel_heights=c(1, 1), align='v'
#   )}
# )()

(myfigs[[4]] <- function(){
  pdftools_installed <- require(pdftools)
  plot_grid(
    NULL, 
    myplots[['autoactiv_model_phase']] +
      # theme(axis.title.y = element_text(margin=margin(l=-2, r=4))) +
      ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=1),
                                                      x=-.22, y=-4, scale=2.5) ) +
      NULL,
    NULL,
    myplots[['twocmps_model_phase']] +
      ( if (!pdftools_installed) NULL else draw_image(magick::image_read_pdf(here("material", "autoactivation.ai.pdf"), pages=3),
                                                      x=-.2, y=-3.7, scale=5) ) +
      NULL,
    labels=c("A  auto-activating operon", "", "B  two-components system", ""), ncol=1, rel_heights=c(.1, 1.2, .1, 1.2), align='v', hjust=-.05
  )}
)()

save_plot(here("plots", "MoM_lacDilution_fig4.pdf"), myfigs[[4]](),
          base_height=NULL, base_width=2.25 * 14/8, # 2 cols
          base_aspect_ratio = 1/1.5
)

