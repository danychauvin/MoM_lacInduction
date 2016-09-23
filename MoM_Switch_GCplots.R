# # SAVE ENVIRONMENT (but `pls`)
# myvars <- ls(all.names = TRUE)
# save(list=myvars[which(myvars!='pls')], file=".RData", envir=.GlobalEnv)


# myframes %>% 
#   filter((date=='20150703' & pos==0 & gl==9)) %>% 
#   filter(!discard_top, vertical_top>vertical_cutoff) %>% 
#   mutate(time_sec=time_sec-2*3600) %>% 
#   plot_faceted_var_tracks(.var_col='length_um', .show_cellid=TRUE, .log=TRUE) +
#   # mask early frames (requires a dummy df!)
#   geom_rect(aes(xmin=-Inf, xmax=0, ymin=0, ymax=Inf, group=1), fill='white', alpha=.6, col='transparent', data=data.frame(a=1)) +
#   scale_x_hours(4) +
#   scale_y_continuous(breaks=2:4, trans='log2') +
#   labs(y='length (µm)')

# Plot overall experiment
pls <- myframes %>% 
  # filter((date=='20150703' & pos==0 & gl==9)) %>%
  #   filter(condition=='switch4h') %>% 
  group_by(condition, date, pos, gl) %>%
  do(pll=(function(.df){
    # browser()
    if (dim(filter(.df, !discard_top))[1] == 0) return(list())
    .cond <- unique(.df$condition)
    .fill <- RColorBrewer::brewer.pal(3, 'Set1')
    .hmin <- max(1.2,  min(filter(.df, !discard_top)$length_um) )
    .hrange <- (max(filter(.df, !discard_top)$length_um)-.hmin) / 2
    .fmin <- min(filter(.df, !discard_top)$gfp_nb)
    .frange <- (max(filter(.df, !discard_top)$gfp_nb)-.fmin) / 2
    
    custom_labels <- function (.str) {
      .labels <- paste('rank:', .str)
      .labels[.labels=='rank: -1'] <- 'all'
      .labels[.labels=='rank: 6'] <- 'rank: >= 6'
      return(.labels)
    }
    
    list(
      l = filter(.df, !discard_top, vertical_top>vertical_cutoff) %>% 
        mutate(time_sec=time_sec-2*3600) %>% 
        plot_faceted_var_tracks(.var_col='length_um', .show_all=TRUE, .show_cellid=TRUE, .log=TRUE, .facet_labeller=custom_labels) +
        # show medium bar
        geom_rect(aes(xmin=t_start - 2*3600, xmax=t_end - 2*3600, ymin=hmin, ymax=hmax, fill=medium, group=1), size=0.2, data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, hmin=.hmin-.25, hmax=.hmin-.1)) +
        geom_vline(aes(xintercept=t_start - 2*3600, group=1), alpha=.2, size=.5, data=filter(condition_ts, condition==.cond)) +
        geom_text(aes(x=t_start+(t_end-t_start)/2 - 2*3600, y=h, label=medium, group=1), col='white', size=2, hjust=0.5, vjust=0,
                  data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, h=.hmin-.2)) +
        scale_fill_manual(values=c('glucose'=.fill[1], 'lactose'=.fill[2]), guide='none') +
        # mask early frames (requires a dummy df!)
        geom_rect(aes(xmin=-Inf, xmax=0, ymin=0, ymax=Inf, group=1), fill='white', alpha=.6, col='transparent', data=data.frame(a=1)) +
        scale_x_hours(4) +
        scale_y_continuous(breaks=2:4, trans='log2') +
        labs(y='length (µm)'), 
      
      ft = filter(.df, !discard_top, vertical_top>vertical_cutoff) %>% 
        mutate(time_sec=time_sec-2*3600) %>% 
        plot_faceted_var_tracks(.var_col='gfp_nb', .show_all=TRUE, .show_cellid=TRUE, .facet_labeller=custom_labels) +
        # show medium bar
        geom_rect(aes(xmin=t_start - 2*3600, xmax=t_end - 2*3600, ymin=fmin, ymax=fmax, fill=medium), size=0.2, data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, fmin=.fmin-.frange/2.5, fmax=.fmin-.frange/5)) +
        geom_vline(aes(xintercept=t_start - 2*3600), alpha=.2, size=.5, data=filter(condition_ts, condition==.cond, t_start>0)) +
        geom_text(aes(x=(t_start+(t_end-t_start)/2 - 2*3600), y=f, label=medium), col='white', size=2, hjust=0.5, vjust=0,
                  data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, f=.fmin-.frange/3)) +
        scale_fill_manual(values=c('glucose'=.fill[1], 'lactose'=.fill[2]), guide='none') +
        # mask early frames (requires a dummy df!)
        geom_rect(aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf, group=1), fill='white', alpha=.6, col='transparent', data=data.frame(a=1)) +
        scale_x_hours(4) +
        # scale_y_continuous(breaks=2:4, trans='log2') +
        labs(y='total GFP per cell (molecules)'),
      
      fc = filter(.df, !discard_top, vertical_top>vertical_cutoff) %>% 
        mutate(time_sec=time_sec-2*3600, gfp_conc=gfp_nb/length_um) %>% 
        plot_faceted_var_tracks(.var_col='gfp_conc', .show_all=TRUE, .show_cellid=TRUE, .facet_labeller=custom_labels) +
        # show medium bar
        geom_rect(aes(xmin=t_start - 2*3600, xmax=t_end - 2*3600, ymin=fmin/2, ymax=fmax/2, fill=medium), size=0.2, data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, fmin=.fmin-.frange/2.5, fmax=.fmin-.frange/5)) +
        geom_vline(aes(xintercept=t_start - 2*3600), alpha=.2, size=.5, data=filter(condition_ts, condition==.cond, t_start>0)) +
        geom_text(aes(x=(t_start+(t_end-t_start)/2-120), y=f/2, label=medium), col='white', size=2, hjust=0.5, vjust=0, 
                  data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, f=.fmin-.frange/3)) +
        scale_fill_manual(values=c('glucose'=.fill[1], 'lactose'=.fill[2]), guide='none') +
        # mask early frames (requires a dummy df!)
        geom_rect(aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf, group=1), fill='white', alpha=.6, col='transparent', data=data.frame(a=1)) +
        scale_x_hours(4) +
        # scale_y_continuous(breaks=2:4, trans='log2') +
        labs(y='total GFP per cell (molecules)') 
    )
  })(.))


pdf('plots/switch_path_length_all.pdf', width=12, height=10)
for (i in 1:dim(pls)[1]) {
  if ('l' %in% names(pls[[i, 'pll']]))
    plot(pls[[i, 'pll']] [['l']] + 
           labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
}
dev.off()

pdf('plots/switch_path_fluotot_all.pdf', width=12, height=10)
for (i in 1:dim(pls)[1]) {
  if ('ft' %in% names(pls[[i, 'pll']]))
    plot(pls[[i, 'pll']] [['ft']] + 
           labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
}
dev.off()

pdf('plots/switch_path_fluoconc_all.pdf', width=12, height=10)
for (i in 1:dim(pls)[1]) {
  if ('fc' %in% names(pls[[i, 'pll']]))
    plot(pls[[i, 'pll']] [['fc']] + 
           labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
}
dev.off()


# Kymograph plots

# ggplot(data=swi_frames, aes(group=interaction(date, pos, gl, id))) + 
#   #   geom_ribbon(aes(frame, ymin=-hmin, ymax=-hmax, fill=factor(id)), alpha=.15) +
#   geom_ribbon(aes(frame, ymin=-(hcenter-height/2), ymax=-(hcenter+height/2), fill=factor(id)), alpha=.3) +
#   geom_path(aes(frame, -hcenter, col=factor(id)), alpha=.5) +
#   scale_colour_periodic_brewer() +
#   scale_fill_periodic_brewer() +
#   guides(col="none", fill="none")

plk <- group_by(swi_frames, date, pos, gl) %>%
  do(pl=(function(.df){
    #     browser()
    .df_div <- filter(swi_frames_divs, date==unique(.df$date), pos==unique(.df$pos), gl==unique(.df$gl))
    ggplot(data=.df, aes(group=interaction(date, pos, gl, id))) + 
      geom_rect(aes(xmin=t_start, xmax=t_end, ymin=-Inf, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), data=filter(condition_ts, condition=='switch4h', medium=='lactose')) +
      geom_rect(aes(xmin=dt*(frame-.5), xmax=dt*(frame+.5), ymin=-dl*(hcenter-height/2), ymax=-dl*(hcenter+height/2), fill=fluo_signal/height)) +
      geom_path(aes(dt*frame, -dl*hcenter)) +
      geom_path(aes(dt*frame, -dl*hcenter), data=.df_div, linetype="dotted") +
      geom_rect(aes(xmin=t_start, xmax=t_end, ymin=-Inf, ymax=Inf, group=1), col='red', fill="transparent", data=filter(condition_ts, condition=='switch4h', medium=='lactose')) +
      scale_fill_gradient2(low="gray50", high="green", midpoint=200) + 
      labs(x="time (min)", y="position (µm)", "GFP concentration")
  })(.))

pdf('plots/switch_kymo_fluoconc.pdf', width=12, height=6)
for (i in 1:dim(pls)[1])
  plot(plk[[i, 'pl']] + 
         labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
dev.off()


