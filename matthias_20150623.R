setwd("~/Desktop/MM_Analysis_Manuscript/switching/20150703_pos0")

dt <- 3             # frame intervall (min)
dh <- 0.065         # pixel size (µm)
hmin_cutoff <- 50   # after it touched this coordinate a cell is discarded
switch_ts <- data.frame(t_start=c(0, 360, 600, 840, 1080, 1320), 
                        t_end=c(360, 600, 840, 1080, 1320, 1560), 
                        cond=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'))


swi_files <- list.files(".", ".*\\d+\\.csv", recursive=TRUE, full.names=TRUE)


# load perl scripts output to dataframes
# NB: creating a list of dataframes and rbind them at last is faster than using rbind in a loop
source("/Users/matthiaskaiser/Desktop/Julou_Test/MoM_toolbox.R")
l_data <- list(cells=list(), frames=list())
for (file in myfiles) {
  .dat <- load_timm_data(file, .verbose=TRUE)#, .force=TRUE)
  if (!is.null(.dat))
    l_data <- mapply(function(.x1, .x2) c(.x1, list(.x2)), 
                     l_data, .dat, SIMPLIFY = FALSE)
}

swi_data <- lapply(l_data, function(.var_df) do.call(rbind, .var_df) )
swi_cells <- swi_data$cells

swi_frames <- group_by(swi_data$frames, date, pos, gl, id) %>%
  mutate(discard_top=(function(.h) {
    .t <- which(.h<hmin_cutoff)
    if (length(.t) == 0) {
      return(rep(FALSE, length(.h)))
    } else {
      return(c(rep(FALSE, .t[1]),
               rep(TRUE, length(.h)-.t[1])))
    }
  })(hmin)) %>%
#   mutate(discard_top=FALSE) %>%
  mutate(b_rank=mean(tot_rank-rank)) %>%
  ungroup %>%
  filter(discard_top==FALSE) %>%
  mutate(hcenter=(hmin+hmax)/2,
         cell=interaction(date, pos, gl, id, drop=TRUE))

# compute connections to parent
swi_frames_divs <- group_by(swi_frames, date, pos, gl) %>%
  do( (function(.df1){
    # first loop on all growth lanes
    .df1 <- mutate(.df1, cell=interaction(date, pos, gl, id, drop=TRUE))
    group_by(.df1, id) %>%
      do( (function(.df1, .df2){
        # then loop on all cells
        if (unique(.df2$parent_id) < 0) return(data.frame())
        .dfp <- filter(.df1, 
                       cell==factor(interaction(unique(.df2$date), unique(.df2$pos), unique(.df2$gl), unique(.df2$parent_id)), 
                                    levels=levels(.df1$cell)))
        .out <- rbind(filter(.df2, frame==min(frame)),
                      filter(.dfp, frame==max(frame)))
        .out <- mutate(.out, id_ini=id, id=unique(.df2$id))
        return(.out)
      })(.df1, .) )
  })(.) )


# Plot overall experiment
pls <- group_by(swi_frames, date, pos, gl) %>%
  do(pll=(function(.df){
    .df_div <- filter(swi_frames_divs, date==unique(.df$date), pos==unique(.df$pos), gl==unique(.df$gl))
    list(
      hf = ggplot(data=mutate(.df, b_rank=-round(b_rank)) %>% filter(b_rank>-8)) + 
        geom_rect(aes(xmin=t_start, xmax=t_end, ymin=0, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), col='red', data=filter(switch_ts, cond=='lactose')) +
        geom_path(aes(dt*frame, dh*height, col=factor(id), group=interaction(date, pos, gl, id)), alpha=.85) +
        facet_grid(b_rank~.) +
        scale_y_log10(breaks=c(2, 4, 8)) +
        scale_colour_periodic_brewer() +
        labs(x="time (min)", y="height (µm)") +
        guides(col="none"),
      
      h = ggplot(data=mutate(.df, b_rank=-round(b_rank)) %>% filter(b_rank>-8)) + 
        geom_rect(aes(xmin=t_start, xmax=t_end, ymin=0, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), col='red', data=filter(switch_ts, cond=='lactose')) +
        geom_path(aes(dt*frame, dh*height, col=-b_rank, group=interaction(date, pos, gl, id)), alpha=.85) +
        geom_path(aes(dt*frame, dh*height, group=interaction(date, pos, gl, id)), data=.df_div, linetype="dotted") +
        scale_y_log10(breaks=c(2, 4, 8)) +
        labs(x="time (min)", y="height (µm)", col="position rank"),
      
      ftf = ggplot(data=mutate(.df, b_rank=-round(b_rank)) %>% filter(b_rank>-8)) + 
        geom_rect(aes(xmin=t_start, xmax=t_end, ymin=-Inf, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), col='red', data=filter(switch_ts, cond=='lactose')) +
        geom_path(aes(dt*frame, fluo_signal, col=factor(id), group=interaction(date, pos, gl, id)), alpha=.85) +
        facet_grid(b_rank~.) +
        scale_colour_periodic_brewer() +
        labs(x="time (min)", y="total GFP") +
        guides(col="none"),

      ft = ggplot(data=mutate(.df, b_rank=-round(b_rank)) %>% filter(b_rank>-8)) + 
        geom_rect(aes(xmin=t_start, xmax=t_end, ymin=0, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), col='red', data=filter(switch_ts, cond=='lactose')) +
        geom_path(aes(dt*frame, fluo_signal, col=-b_rank, group=interaction(date, pos, gl, id)), alpha=.85) +
        geom_path(aes(dt*frame, fluo_signal, group=interaction(date, pos, gl, id)), data=.df_div, linetype="dotted") +
        labs(x="time (min)", y="total GFP", col="position rank"),
      
      fcf = ggplot(data=mutate(.df, b_rank=-round(b_rank)) %>% filter(b_rank>-8)) + 
        geom_rect(aes(xmin=t_start, xmax=t_end, ymin=-Inf, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), col='red', data=filter(switch_ts, cond=='lactose')) +
        geom_path(aes(dt*frame, fluo_signal/height, col=factor(id), group=interaction(date, pos, gl, id)), alpha=.85) +
        facet_grid(b_rank~.) +
        scale_colour_periodic_brewer() +
        labs(x="time (min)", y="GFP concentration") +
        guides(col="none"),
      
      fc = ggplot(data=mutate(.df, b_rank=-round(b_rank)) %>% filter(b_rank>-8)) + 
        geom_rect(aes(xmin=t_start, xmax=t_end, ymin=0, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), col='red', data=filter(switch_ts, cond=='lactose')) +
        geom_path(aes(dt*frame, fluo_signal/height, col=-b_rank, group=interaction(date, pos, gl, id)), alpha=.85) +
        geom_path(aes(dt*frame, fluo_signal/height, group=interaction(date, pos, gl, id)), data=.df_div, linetype="dotted") +
        labs(x="time (min)", y="GFP concentration", col="position rank")
      
    )
  })(.))

pdf('switch_path_height.pdf', width=12, height=6)
for (i in 1:dim(pls)[1])
  plot(pls[[i, 'pll']] [['hf']] + 
         labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
dev.off()

pdf('switch_path_fluotot.pdf', width=12, height=6)
for (i in 1:dim(pls)[1])
  plot(pls[[i, 'pll']] [['ftf']] + 
         labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
dev.off()


# Chimograph plots
ggplot(data=swi_frames, aes(group=interaction(date, pos, gl, id))) + 
#   geom_ribbon(aes(frame, ymin=-hmin, ymax=-hmax, fill=factor(id)), alpha=.15) +
  geom_ribbon(aes(frame, ymin=-(hcenter-height/2), ymax=-(hcenter+height/2), fill=factor(id)), alpha=.3) +
  geom_path(aes(frame, -hcenter, col=factor(id)), alpha=.5) +
  scale_colour_periodic_brewer() +
  scale_fill_periodic_brewer() +
  guides(col="none", fill="none")

pls <- group_by(swi_frames, date, pos, gl) %>%
  do(pl=(function(.df){
#     browser()
    .df_div <- filter(swi_frames_divs, date==unique(.df$date), pos==unique(.df$pos), gl==unique(.df$gl))
    ggplot(data=.df, aes(group=interaction(date, pos, gl, id))) + 
      geom_rect(aes(xmin=t_start, xmax=t_end, ymin=-Inf, ymax=Inf, group=1), fill=rgb(1, 0, 0, .1), data=filter(switch_ts, cond=='lactose')) +
      geom_rect(aes(xmin=dt*(frame-.5), xmax=dt*(frame+.5), ymin=-dh*(hcenter-height/2), ymax=-dh*(hcenter+height/2), fill=fluo_signal/height)) +
      geom_path(aes(dt*frame, -dh*hcenter)) +
      geom_path(aes(dt*frame, -dh*hcenter), data=.df_div, linetype="dotted") +
      geom_rect(aes(xmin=t_start, xmax=t_end, ymin=-Inf, ymax=Inf, group=1), col='red', fill="transparent", data=filter(switch_ts, cond=='lactose')) +
      scale_fill_gradient2(low="gray50", high="green", midpoint=200) + 
      labs(x="time (min)", y="position (µm)", "GFP concentration")
    })(.))

pdf('switch_chimo_fluoconc.pdf', width=12, height=6)
for (i in 1:dim(pls)[1])
  plot(pls[[i, 'pl']] + 
         labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
dev.off()

