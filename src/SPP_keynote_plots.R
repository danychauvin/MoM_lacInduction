(function(.ylabel=1.3)
  myframes %>% ungroup() %>% 
  filter(date==20150708, pos==5) %>% #pull(gl) %>% table()
  # group_by(gl) %>% nest() %>% slice(8L) %>% unnest() %>% 
  filter(gl==26, id!=144, length_um < 4.5) %>% 
  filter(!discard_start, !discard_top,
         time_sec>2*3600) %>%
  filter(total_cell_in_lane-cell_num_in_lane<5,
         end_type %in% c('div', 'lost')) %>% 
  # ungroup() %>% mutate(gfp_nb=gfp_nb - 2*min(gfp_nb)) %>% 
  ungroup() %>% mutate(gfp_nb=gfp_nb + 135) %>% 
  gather(variable, value, length_um, gfp_nb) %>% 
  ggplot() +
  # facet_grid(pos+gl~.) +
  facet_wrap(~variable, scales='free_y', ncol=1, strip.position='left',
             labeller=as_labeller(function(.str) {
               .str <- str_replace(.str, "gfp_nb", "LacZ-GFP\n(molecules)")
               .str <- str_replace(.str, "length_um", "cell length (µm)")
               return(.str)
             }) ) +
  geom_line(aes(time_sec-2*3600, value, col=ugen)) +
  scale_colour_periodic(.n=3, guide='none') +
  # show medium bar
  geom_vline(aes(xintercept=t_start - 2*3600, group=1), alpha=.2, size=.5, 
             data=filter(condition_ts, condition=='switch_04h')) +
  geom_rect(aes(xmin=t_start - 2*3600, xmax=t_end - 2*3600, ymin=.ylabel, ymax=1.5, fill=medium, group=1), show.legend=FALSE, #size=0.2,
            data=filter(condition_ts, condition=='switch_04h') %>% mutate(variable='length_um')) +
  geom_text(aes(x=t_start+(t_end-t_start)/2 - 2*3600, y=.ylabel, label=medium, group=1), size=4, col='white', hjust=0.5, vjust=-0.5, 
            data=filter(condition_ts, condition=='switch_04h') %>% mutate(variable='length_um', t_start=ifelse(t_start<2*3600, 2*3600, t_start))) +
  # scale_fill_manual(values=c('glucose'=.fill[1], 'lactose'=.fill[2]), guide='none') +
  coord_cartesian(xlim=c(-3600/2, 23.9*3600), expand=FALSE) +
  scale_x_hours(4) +
  scale_y_continuous(trans='log10', breaks=c(2, 3, 4, 100, 500, 2500)) +
  theme(strip.placement='outside', strip.background.y=element_blank(),
        axis.title.y=element_blank(), strip.text.y=element_text(size=rel(1.2))) +
  NULL
)()

save_plot(here("plots","keynote_raw_traces_log.pdf"), last_plot(), base_width=8, base_height=4)


myframes %>% ungroup() %>% 
  filter(date==20150708, pos==5) %>% #pull(gl) %>% table()
  # group_by(gl) %>% nest() %>% slice(8L) %>% unnest() %>% 
  filter(gl==26, id!=144, length_um < 4.5) %>% 
  filter(!discard_start, !discard_top,
         time_sec>2*3600) %>%
  filter(total_cell_in_lane-cell_num_in_lane<1,
         end_type %in% c('div', 'lost')) %>% 
  (function(.df, .ylabel=1.25) {
    p_length <- ggplot(.df) +
      geom_line(aes(time_sec-2*3600, length_um, group=ugen)) +
      scale_y_continuous(name="cell length (µm)", trans='log10', breaks=2:4) +
      scale_colour_periodic(.n=3, guide='none') +
      # show medium bar
      geom_vline(aes(xintercept=t_start - 2*3600, group=1), alpha=.2, size=.5, 
                 data=filter(condition_ts, condition=='switch_04h')) +
      geom_rect(aes(xmin=t_start - 2*3600, xmax=t_end - 2*3600, ymin=.ylabel, ymax=1.5, fill=medium, group=1), show.legend=FALSE, #size=0.2,
                data=filter(condition_ts, condition=='switch_04h')) +
      geom_text(aes(x=t_start+(t_end-t_start)/2 - 2*3600, y=.ylabel, label=medium, group=1), size=4, col='white', hjust=0.5, vjust=-0.5, 
                data=filter(condition_ts, condition=='switch_04h') %>% mutate(t_start=ifelse(t_start<2*3600, 2*3600, t_start))) +
      coord_cartesian(xlim=c(-3600/2, 23.9*3600), expand=FALSE) +
      scale_x_hours(4) +
      theme(plot.margin=margin(0, 14/2, 14/2, 14/2)) +
      NULL
    
    p_gfp <- ggplot(.df) +
      geom_line(aes(time_sec-2*3600, gfp_nb, group=ugen), col=qual_cols[3]) +
      scale_colour_periodic(.n=3, guide='none') +
      # show medium bar
      geom_vline(aes(xintercept=t_start - 2*3600, group=1), alpha=.2, size=.5, 
                 data=filter(condition_ts, condition=='switch_04h')) +
      coord_cartesian(xlim=c(-3600/2, 23.9*3600), expand=FALSE) +
      scale_x_hours(4) +
      labs(y="LacZ-GFP\n(molecules)") +
      # theme(plot.margin=margin(14/2, 14/2, 0, 14/2),
            # axis.title.x = element_blank(), axis.text.x = element_blank(),
            # axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
      NULL
    
    print(p_gfp)
    print(p_length)
    plot_grid(p_gfp    + coord_cartesian(xlim=c(0*3600, 8.2*3600)) + scale_y_continuous(position = 'right', limits=c(-1500, 4500), breaks=c(0, 2e3, 4e3)), 
              p_length + coord_cartesian(xlim=c(0*3600, 8.2*3600)),
              ncol=1, align = 'v', rel_heights = c(1, 1))
  })
save_plot(here("plots","keynote_trace_zoom.pdf"), last_plot(), base_width=5, base_height=4)

myframes %>% ungroup() %>% 
  filter(date==20150708, pos==5) %>% #pull(gl) %>% table()
  # group_by(gl) %>% nest() %>% slice(8L) %>% unnest() %>% 
  filter(gl==26, id!=144, length_um < 4.5) %>% 
  filter(!discard_start, !discard_top,
         time_sec>2*3600) %>%
  filter(total_cell_in_lane-cell_num_in_lane<5,
         end_type %in% c('div', 'lost')) %>% 
  (function(.df, .ylabel=1.28) {
    p_length <- ggplot(.df) +
      geom_line(aes(time_sec-2*3600, length_um, col=ugen)) +
      scale_y_continuous(name="cell length (µm)", trans='log10', breaks=2:4) +
      scale_colour_periodic(.n=3, guide='none') +
      # show medium bar
      geom_vline(aes(xintercept=t_start - 2*3600, group=1), alpha=.2, size=.5, 
                 data=filter(condition_ts, condition=='switch_04h')) +
      geom_rect(aes(xmin=t_start - 2*3600, xmax=t_end - 2*3600, ymin=.ylabel, ymax=1.5, fill=medium, group=1), show.legend=FALSE, #size=0.2,
                data=filter(condition_ts, condition=='switch_04h')) +
      geom_text(aes(x=t_start+(t_end-t_start)/2 - 2*3600, y=.ylabel, label=medium, group=1), size=4, col='white', hjust=0.5, vjust=-0.5, 
                data=filter(condition_ts, condition=='switch_04h') %>% mutate(t_start=ifelse(t_start<2*3600, 2*3600, t_start))) +
      coord_cartesian(xlim=c(-3600/2, 23.9*3600), expand=FALSE) +
      scale_x_hours(4) +
      theme(plot.margin=margin(0, 14/2, 14/2, 14/2)) +
      NULL
    
    p_gfp <- ggplot(.df) +
      geom_line(aes(time_sec-2*3600, gfp_nb, col=ugen)) +
      scale_colour_periodic(.n=3, guide='none') +
      # show medium bar
      geom_vline(aes(xintercept=t_start - 2*3600, group=1), alpha=.2, size=.5, 
                 data=filter(condition_ts, condition=='switch_04h')) +
      coord_cartesian(xlim=c(-3600/2, 23.9*3600), expand=FALSE) +
      scale_x_hours(4) +
      labs(y="LacZ-GFP\n(molecules)") +
      theme(plot.margin=margin(14/2, 14/2, 0, 14/2),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.line.x = element_blank(), axis.ticks.x = element_blank()) +
      NULL
    
    print(p_gfp)
    print(p_length)
    plot_grid(p_gfp, p_length,
              ncol=1, align = 'v', rel_heights = c(1, 1.25))
  })
save_plot(here("plots","keynote_raw_traces.pdf"), last_plot(), base_width=9, base_height=4)

myplots[['naive_lags_hist']] 
save_plot(here("plots","naive_lags_hist.pdf"), last_plot(), base_width=4, base_height=1.8)

myplots[['basal_perturb_pre_violin']] +
  labs(y=expression(paste(italic("lac"), " induction lag (min)    "))) +
  theme(legend.position = "none")
save_plot(here("plots","basal_perturb_pre_violin.pdf"), last_plot(), base_width=7, base_height=3)




mycells_switching %>% ungroup %>% 
    filter(!date %in% discarded_dates) %>%
    filter(str_detect(condition, '^switch_[0-9]+h$')) %>%
    filter(switch_idx==1, !is.na(lag_200)) %>% 
    mutate(type=ifelse(lag_200/60<50, 'short', 'long')) %>% 
    group_by(date) %>%
    mutate(gfp_ini=(gfp_ini-mean(gfp_ini))) %>%
    ggplot(aes(gfp_ini, lag_200/60)) + 
    geom_hline(yintercept = 50, lty='dashed') +
    geom_point(aes(col=type), size=1, alpha=.2, stroke=0) +
    stat_smooth(method=MASS::rlm, col='black', fullrange = TRUE) +
    labs(x='fluorescence at the switch\n(centered per experiment; equivalent GFP molecules)', y=lac_lags_label) +
    guides(col=guide_legend(override.aes = list(size=2, alpha=1))) +
    coord_cartesian(xlim=c(-60, 120)) +
    NULL

mycells_switching %>% ungroup %>% 
  filter(!date %in% discarded_dates) %>%
  filter(str_detect(condition, '^switch_[0-9]+h$')) %>%
  filter(switch_idx==1, !is.na(lag_200)) %>% 
  mutate(type=ifelse(lag_200/60<50, 'short', 'long')) %>%
  left_join(myframes_switching %>% ungroup %>%  #filter(ugen %in% unique(ugen)[1:5]) %>% 
              filter(time_sec<t_lac_switch) %>% 
              # ggplot() + geom_line(aes(time_sec-t_lac_switch, length_um, lty=type, col=ugen)) +
              # scale_colour_periodic_brewer(guide='none')
              group_by(ugen) %>% 
              filter(time_sec==min(time_sec), time_sec>0) %>% 
              select(ugen, length_um_birth=length_um) ) %>% 
  mutate(cc=(t_lac_switch-time_birth) / 
           # average growth rate in glucose
           ( mycells_constant %>% filter(condition=='glucose') %>% 
           mutate(dt=time_div-time_birth) %>% pull(dt) %>% mean ),
         dl=length_ini - length_um_birth,
         growth_rate=logl_time_slope_before/log(2)*3600) %>% 
  select(date, lag_200, type, gfp_ini, growth_rate, cc, dl) %>% 
  # group_by(date) %>% 
  # mutate(gfp_ini=(gfp_ini-mean(gfp_ini))) %>%
  filter(cc < 1.5) %>% 
  gather(variable, value, -date, -lag_200, -type) %>%
  group_by(variable) %>% 
  filter(between(value, quantile(value, 0.01, na.rm=TRUE), quantile(value, 0.99, na.rm=TRUE))) %>% 
  (function(.df)
    ggplot(.df, aes(value, lag_200/60)) + 
     facet_wrap(~variable, nrow=1, scales='free_x', strip.position='bottom',
                labeller=as_labeller(function(.str) {
                  .str <- str_replace(.str, "cc", "cell cycle")
                  .str <- str_replace(.str, "dl", "added length\n(µm)")
                  .str <- str_replace(.str, "gfp_ini", "LacZ-GFP")
                  .str <- str_replace(.str, "growth_rate", "growth rate\n(dbl/h)")
                  return(.str)
                })) +
     geom_hline(yintercept = 50, lty='dashed') +
     geom_point(aes(col=type), size=1, alpha=.2, stroke=0) +
     # stat_smooth(method=MASS::rlm, col='black', fullrange = TRUE) + 
     geom_text(aes(-Inf, 0, label=r2str), size=5, hjust=-0.1, vjust=0.8, parse=TRUE,
               data=.df %>% filter(is.finite(lag_200)) %>% group_by(variable) %>% summarise(r2=cor(value, lag_200)^2) %>% mutate(r2str=sprintf("r^2 == '%.4f'", round(r2, digits=4)))) +
     scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
     expand_limits(y=0) +
     coord_cartesian(ylim=c(-25, 240)) +
     labs(x='fluorescence at the switch\n(centered per experiment; equivalent GFP molecules)', y=expression(paste(italic("lac"), " induction lag (min)"))) +
     theme(legend.position = 'none',
           strip.placement='outside', strip.background.x=element_blank(),
           axis.title.x=element_blank(), strip.text.x=element_text(size=rel(1.2), vjust=1)) +
     NULL )
save_plot(here("plots","naive_lags_correl.pdf"), last_plot(), base_width=6, base_height=2.6)
