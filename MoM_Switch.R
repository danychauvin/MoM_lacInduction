
# SET VARIABLES & ENVIRONMENT ####

# SET VARIABLES
dt <- 3             # frame intervall (min)
dl <- 0.065         # pixel size (µm)
vertical_cutoff <- 4 / dl   # after it touched this coordinate a cell is discarded
use_eriks_params <- TRUE

proj_path <- "~/Documents/Biozentrum/Projects/MoM_Switch"
r_scripts_path <- c("~/Documents/Biozentrum/Projects/vngWetLabR/mother_machine",
                  "~/Documents/Biozentrum/Projects/vngWetLabR/ggplot")
perl_scripts_path <- "~/Documents/Biozentrum/Projects/vngWetLabR/mother_machine"
data2preproc <- function(.d) sub('/(data_matthias|data_thomas)/*', '/preproc/', .d) # store cache file in preproc subdir

# # local paths
# proj_path <- "~/Downloads/january/MoM_Switch"
# r_scripts_path <- c("~/Downloads/january/vngWetLabR/mother_machine",
#                   "~/Downloads/january/vngWetLabR/ggplot")
# perl_scripts_path <- "~/Downloads/january/vngWetLabR/mother_machine"

date_cond <- c("20150616"="glucose", "20150617"="glucose", 
               "20150624"="lactose", "20150630"="lactose", 
               "20160318"="lactose_lowillum", 
               "20150703"="switch_4h", "20150708"="switch_4h",
               "20151204"="switch_6h",
               "20151207"="switch_iptg",
               "20151218"="switch_8h",
               "20151221"="switch_m9" )

condition_ts <- rbind(data.frame(duration=1560, medium='glucose', condition='glucose'),
                      data.frame(duration=1560, medium='lactose', condition='lactose'),
                      data.frame(duration=1560, medium='lactose', condition='lactose_lowillum'),
                      data.frame(duration=c(360, 240, 240, 240, 240, 240),
                                 medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
                                 condition='switch_4h'),
                      data.frame(duration=c(360, 240, 360, 240, 360, 240), 
                                 medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
                                 condition='switch_6h'),
                      data.frame(duration=c(360, 240, 480, 240, 480, 240), 
                                 medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
                                 condition='switch_8h'),
                      data.frame(duration=c(360, 240, 240, 240, 240, 240),
                                 medium=c('glucose', 'lactose+IPTG', 'glucose', 'lactose+IPTG', 'glucose', 'lactose+IPTG'),
                                 condition='switch_iptg'),
                      data.frame(duration=c(360, 120, 360, 720, 360, 240), 
                                 medium=c('glucose', 'M9', 'glucose', 'M9', 'glucose', 'M9'),
                                 condition='switch_m9') 
)


# SET ENVIRONMENT
invisible(sapply(
  list.files(r_scripts_path, pattern="\\.[Rr]$", full.names=TRUE, ignore.case=TRUE), 
  source, .GlobalEnv))
setwd(proj_path)

library(parallel)
numCores <- min(30, detectCores()-1) # do not use more than 30 cores

library(multidplyr)
mycluster <- create_cluster(numCores)
for (.l in mylibs) { mycluster %>% cluster_library(.l) }
# set_default_cluster(mycluster)

# SAVE ENVIRONMENT (but `pls`)
# myvars <- ls(all.names = TRUE)
# save(list=myvars[which(myvars!='pls')], file=".RData", envir=.GlobalEnv)

condition_ts <- condition_ts %>% group_by(condition) %>% 
  mutate(t_start=cumsum(c(0, duration[-(length(duration))])) * 60,
         t_end=cumsum(duration) * 60 - 1e-5,
         m_cycle=value_occurence_index(medium))


# LOAD MG1655 DATA ####
mg_files <- list.files("./data_matthias/MG1655_glu_lac", ".*\\d+\\.csv", recursive=TRUE, full.names=TRUE)
nc <- min(length(mg_files), length(mycluster)) # this is a dirty hack because multiplyr crashes with less shards than cores

# load perl scripts output to dataframes (using multidplyr)
mg_frames <- dplyr::data_frame(path=mg_files) %>%
  group_by(path) %>%
  partition(path, cluster=mycluster[1:nc] %>%
              cluster_assign_func(data2preproc, load_timm_data, parse_frames_stats, compute_genealogy, which_touch_exit, which_to_progeny) %>%
              cluster_assign_obj(perl_scripts_path, dt, dl, vertical_cutoff) ) %>%
  do((function(.df){
    # browser()
    .l <- try( load_timm_data(.df$path, perl_scripts_path, .verbose=TRUE, .data2preproc=data2preproc, .force=FALSE))
    if ('frames' %in% names(.l)) {
      return( .l$frames %>%
                mutate(time_sec=frame*dt*60, length_um=length_pixel*dl,
                       discard_start=(time_sec < 2*3600)) %>%
                # remove frames after touching the exit
                group_by(id) %>%
                mutate(discard_top=which_touch_exit(vertical_top, vertical_cutoff)) %>%
                mutate(discard_top=ifelse(discard_start, FALSE, discard_top)) %>% # not in the preexpt step (2h)
                mutate(end_type=ifelse(any(discard_top), 'exit', end_type)) %>% # update end_type
                # remove daughters of cells that touched the exit
                ungroup %>%
                mutate(discard_top=which_to_progeny(discard_top, cid)) )
    } else {
      return (data.frame())
    }
  })(.)) %>%
  collect() %>% 
  group_by(date, pos, gl, id) %>%
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         b_rank=round(mean(total_cell_in_lane - cell_num_in_lane)),
         cell=paste(date, pos, gl, id, sep='.'))


# LOAD AS662 DATA ####
asc_files <- c("./data_matthias/glucose", "./data_matthias/lactose", 
               "./data_matthias/glu_lac_switch", "./data_thomas/") %>% 
  list.files(".*\\d+\\.csv", recursive=TRUE, full.names=TRUE)
nc <- min(length(asc_files), length(mycluster)) # this is a dirty hack because multiplyr crashes with less shards than cores

# load perl scripts output to dataframes (using multidplyr)
myframes <- dplyr::data_frame(path=asc_files) %>%
  group_by(path) %>%
  partition(path, cluster=mycluster[1:nc] %>%
              cluster_assign_func(data2preproc, load_timm_data, parse_frames_stats, compute_genealogy, which_touch_exit, which_to_progeny) %>%
              cluster_assign_obj(perl_scripts_path, dt, dl, vertical_cutoff) ) %>%
  do((function(.df){
    # browser()
    .l <- try( load_timm_data(.df$path, perl_scripts_path, .verbose=TRUE, .data2preproc=data2preproc, .force=FALSE))
    if ('frames' %in% names(.l)) {
      return( .l$frames %>%
                mutate(time_sec=frame*dt*60, length_um=length_pixel*dl,
                       discard_start=(time_sec < 2*3600)) %>%
                # remove frames after touching the exit
                group_by(id) %>%
                mutate(discard_top=which_touch_exit(vertical_top, vertical_cutoff)) %>%
                mutate(discard_top=ifelse(discard_start, FALSE, discard_top)) %>% # not in the preexpt step (2h)
                mutate(end_type=ifelse(any(discard_top), 'exit', end_type)) %>% # update end_type
                # remove daughters of cells that touched the exit
                ungroup %>%
                mutate(discard_top=which_to_progeny(discard_top, cid)) )
    } else {
      return (data.frame())
    }
  })(.)) %>%
  collect() %>% 
  arrange(date, pos, gl, id, frame) %>%  # sort data after `partition()`
  group_by(date, pos, gl, id) %>%
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         b_rank=round(mean(total_cell_in_lane - cell_num_in_lane))) %>% 
  #   filter(discard_top==FALSE) %>%
  ungroup %>% 
  mutate(condition=date_cond[as.character(date)],
         vertical_center=(vertical_bottom + vertical_top)/2,
         gl_id=paste(date, pos, gl, sep='.') %>% 
         cell=paste(gl_id, id, sep='.')) %>% 
  # propagate medium info
  group_by(condition) %>% 
  do((function(.df){
    .ts <- filter(condition_ts, condition==unique(.df$condition))
    .idx <- find_unique_interval(.df$time_sec, .ts$t_start, .ts$t_end)
    mutate(.df, medium=.ts$medium[.idx], m_start=.ts$t_start[.idx], m_cycle=.ts$m_cycle[.idx])
  })(.))  



# CURATION STATS ####
timm_files <- list.files("./data", "^\\d+_pos.*\\.timm$", recursive=TRUE, full.names=TRUE)

l_data <- lapply(timm_files, function(.f) try(parse_timm_curation(.f) %>%
                                                data.frame(path=.f, .)) )
timm_data <- lapply(l_data, function(.df) {
  if (class(.df) == 'try-error') return(data.frame())
  # keep only the first line of a frame for each SSCxAction events
  # combine it (using rbind) with all other events
  rbind(.df %>% 
          filter(type=='SSC') %>%
          group_by(frame, action) %>%
          summarise_each(funs(first)) %>%
          ungroup %>% 
          select(path, type, frame, action),
        filter(.df, type!='SSC') %>%
          select(path, type, frame, action) ) %>%
    arrange(type, frame)
  }) %>%
  do.call(rbind, .) %>%
  extract(path, c('date', 'pos', 'gl'), ".*/(\\d+)_pos(\\d+)_[^/]*GL0*(\\d+)\\.timm") %>%
  # to do: keep only files for which the output is used
  mutate(date=as.numeric(date), pos=as.numeric(pos), gl=as.numeric(gl))


mygl <- myframes %>%
  # start from the total number of observations
  group_by(condition, date, pos, gl) %>%
  filter(!discard_top, !discard_start) %>% 
  summarise(n_obs=n()) %>%
  # add the number of dividing cells
  left_join(myframes %>%
              group_by(date, pos, gl, id) %>% 
              filter(!discard_top, !any(discard_start), end_type=='div', row_number()==1) %>% 
              group_by(date, pos, gl) %>% 
              summarise(n_div_cells=n()) ) %>%
  # add the number of frames
  left_join(myframes %>%
              group_by(date, pos, gl) %>% 
              summarise(nframes=max(frame)) ) %>%
  # add curation times
  left_join(read.csv("./data_matthias/curation_times.csv", comment.char="#") %>%
              na.omit %>%
              mutate(pos=as.numeric(gsub('pos', '', pos)),
                     gl=as.numeric(gsub('GL', '', gl))) ) %>%
  # add the number of curated frames
  left_join(rbind(
    timm_data %>% 
      filter(type != 'none') %>%
      group_by(date, pos, gl) %>%
      summarise(n_cur_frames=length(unique(frame))),
    timm_data %>% 
      filter(type == 'none') %>%
      group_by(date, pos, gl) %>%
      summarise(n_cur_frames=0)) )


# kable(mygl %>%
#   group_by(condition) %>%
#   summarise(n_lanes=length(unique(interaction(date, pos, gl))),
#             n_div_cells=sum(n_div_cells),
#             n_obs=sum(n_obs),
#             time_avg=mean(time, na.rm=TRUE),
#             time_sd=sd(time, na.rm=TRUE) ))
#
# qplot(time/nframes*100, data=mygl, xlab='curation time (min; per 100 frames)', col=I('darkblue'), fill=I('darkblue'), alpha=I(.4))
# ggsave('plots/curation_times_hist.pdf', width=4, height=3)
# 
# ggplot(data=mygl, aes(n_cur_frames/nframes*100, time/nframes*100)) +
#   geom_smooth(method='lm', se=FALSE) +
#   geom_point(position=position_jitter(width=.1)) +
#   # geom_rug(sides='r', position=position_jitter(), size=5, alpha=.2) +
#   # ylim(0, 20) +
#   expand_limits(y=0) +
#   labs(x='fraction of frames curated (%)', y='curation time (min; per 100 frames)')
# ggsave('plots/curation_times_frames.pdf', width=4, height=3)
# 
# ggplot(data=filter(mygl, n_div_cells>40), aes(n_div_cells, time)) +
#   geom_point(position=position_jitter(width=.1)) +
#   geom_point(data=filter(mygl, n_div_cells>40) %>% ungroup %>% select(n_div_cells, time) %>% mutate_each(funs(as.numeric)) %>% summarise_each(funs(median(., na.rm=T))), col='red', pch='+', size=10) +
#   # geom_rug(sides='r', position=position_jitter(), size=5, alpha=.2) +
#   ylim(0, 20) +
#   labs(x='number of entire cell cycles', y='curation time (min)')
# ggsave('plots/curation_times_divcells.pdf', width=4, height=3)


# RENDER ANALYSIS FILES ####
# calling `render()` or `render_site()` from the command line allows to execute the function 
# in the global env() (hence inheriting existing variables and keeping newly created ones)...

knitr::opts_chunk$set(echo=FALSE, message=FALSE, warning=FALSE)
# rmarkdown::clean_site()

# rmarkdown::render_site('index.Rmd')
rmarkdown::render_site('MoM_Switch_GFP_Estimation.Rmd')
rmarkdown::render_site('MoM_Switch_Constant_Envts.Rmd')
rmarkdown::render_site('MoM_Switch_Switching_Envts.Rmd')
rmarkdown::render_site('MoM_Switch_FCM_Comparison.Rmd')



# CONTROL PLOTS ####

# Plot overall experiment
pls <- myframes %>% 
  # filter((date=='20150703' & pos==0 & gl==7))%>% 
#   filter(condition=='switch4h') %>% 
  group_by(condition, date, pos, gl) %>%
  do(pll=(function(.df){
    # browser()
    if (dim(filter(.df, !discard_top))[1] == 0) return(list())
    .cond <- unique(.df$condition)
    .fill <- brewer.pal(3, 'Set1')
    .hmin <- max(1.2,  min(filter(.df, !discard_top)$length_um) )
    .hrange <- (max(filter(.df, !discard_top)$length_um)-.hmin) / 2
    .fmin <- min(filter(.df, !discard_top)$gfp_nb)
    .frange <- (max(filter(.df, !discard_top)$gfp_nb)-.fmin) / 2
    
    # compute connections to parent
    .df <- mutate(.df, b_rank=ifelse(b_rank>6, 6, b_rank)) 
    .df_div <- filter(.df, !discard_top) %>% 
      group_by(id) %>% 
      do((function(.dfgl, .dfc) { # compute_div_between_facets
        # .dfgl: dataframe of the entire GL
        # .dfc: dataframe of the current cell
        if (unique(.dfc$parent_id) < 0) return(data.frame())
        .dfp <- filter(.dfgl, id==unique(.dfc$parent_id))
        
        if (unique(.dfc$b_rank) == unique(.dfp$b_rank)) {
          .length_ump <- filter(.dfp, time_sec==max(time_sec))$length_um
          .gfp_nbp <- filter(.dfp, time_sec==max(time_sec))$gfp_nb
          .out <- filter(.dfc, time_sec==min(time_sec)) %>%
            select(id, b_rank, time_sec, length_um, gfp_nb) %>%
            mutate(time_secp=time_sec-dt*60, length_ump=.length_ump, gfp_nbp=.gfp_nbp)
        }
        if (unique(.dfc$b_rank) > unique(.dfp$b_rank)) {
          .out <- bind_rows(
            filter(ungroup(.dfc), time_sec==min(time_sec)) %>%
              select(id, b_rank, time_sec, length_um, gfp_nb) %>%
              mutate(time_secp=time_sec-dt*60/2, length_ump=0, gfp_nbp=-Inf),
            filter(ungroup(.dfp), time_sec==max(time_sec)) %>%
              select(id, b_rank, time_sec, length_um, gfp_nb) %>%
              mutate(id=unique(.dfc$id), time_secp=time_sec+dt*60/2, length_ump=Inf, gfp_nbp=Inf) )
          if(diff(.out$b_rank) < -1) {
            .out <- rbind(.out,
                          data.frame(id=unique(.dfc$id), b_rank=(min(.out$b_rank)+1):(max(.out$b_rank)-1),
                                     time_sec=min(.dfc$time_sec)-dt*60/2, time_secp=min(.dfc$time_sec)-dt*60/2,
                                     length_um=0, length_ump=Inf, gfp_nb=-Inf, gfp_nbp=Inf))
          }
        }
        return(.out)
      })(.df, .))
    
    custom_labels <- function (.str) {
      .labels <- paste('rank:', .str)
      .labels[.labels=='rank: -1'] <- 'all'
      .labels[.labels=='rank: 6'] <- 'rank: >= 6'
      return(.labels)
    }
    
    list(
      l = ggplot() +
        facet_grid(b_rank~., scales='free_y', as.table=FALSE, labeller=as_labeller(custom_labels)) +
        # facets alternated background
        geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf), alpha=.05, data=data.frame(b_rank=seq(0, 6, 2))) +
        # show medium bar
        geom_rect(aes(xmin=t_start*60 - 2*3600, xmax=t_end*60 - 2*3600, ymin=hmin, ymax=hmax, fill=medium, group=1), size=0.2, data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, hmin=.hmin-.25, hmax=.hmin-.1)) +
        geom_vline(aes(xintercept=t_start*60 - 2*3600, group=1), alpha=.2, size=.5, data=filter(condition_ts, condition==.cond)) +
        geom_text(aes(x=60*(t_start+(t_end-t_start)/2-120), y=h, label=medium, group=1), col='white', data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, h=.hmin-.2), size=2, hjust=0.5, vjust=0) +
        # show divisions (lty='11' means densely dotted line)
        geom_segment(aes(x=time_sec - 2*3600, xend=time_secp - 2*3600, y=length_um, yend=length_ump, col=factor(id)), alpha=.3, lty='11', data=.df_div) + 
        # show cell traces
        geom_path(aes(time_sec - 2*3600, length_um, col=factor(id)), data=filter(.df, !discard_top, vertical_top>vertical_cutoff)) +
        geom_text(aes(time_sec - 2*3600, length_um, col=factor(id), label=id), size=2, hjust=0, vjust=1,
                  data=filter(.df, !discard_top) %>% group_by(id) %>% filter(row_number()==1) ) +
        # show all traces in one panel
        geom_path(aes(time_sec - 2*3600, length_um, col=factor(id)), data=filter(.df, !discard_top, vertical_top>vertical_cutoff) %>% mutate(b_rank=-1)) +
        # mask early frames (requires a dummy df!)
        geom_rect(aes(xmin=-Inf, xmax=0, ymin=0, ymax=Inf, group=1), fill='white', alpha=.6, col='transparent', data=data.frame(a=1)) +
        scale_colour_periodic_brewer(guide='none') + 
        scale_fill_manual(values=c('glucose'=.fill[1], 'lactose'=.fill[2]), guide='none') +
        scale_x_hours(4) +
        scale_y_continuous(trans='log2', breaks=2:4) +
        labs(y='length (µm)') +
        theme(panel.margin = unit(0, "lines"), panel.border=element_blank()),
      
      ft = ggplot() +
        facet_grid(b_rank~., scales='free_y', as.table=FALSE, labeller=as_labeller(custom_labels)) +
        # facets alternated background
        geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), alpha=.05, data=data.frame(b_rank=seq(0, 6, 2))) +
        # show medium bar
        geom_rect(aes(xmin=t_start*60 - 2*3600, xmax=t_end*60 - 2*3600, ymin=fmin, ymax=fmax, fill=medium), size=0.2, data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, fmin=.fmin-.frange/2.5, fmax=.fmin-.frange/5)) +
        geom_vline(aes(xintercept=t_start*60 - 2*3600), alpha=.2, size=.5, data=filter(condition_ts, condition==.cond, t_start>0)) +
        geom_text(aes(x=60*(t_start+(t_end-t_start)/2-120), y=f, label=medium), col='white', data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, f=.fmin-.frange/3), size=2, hjust=0.5, vjust=0) +
        # show divisions (lty='11' means densely dotted line using hex notation)
        geom_segment(aes(x=time_sec - 2*3600, xend=time_secp - 2*3600, y=gfp_nb, yend=gfp_nbp, col=factor(id)), alpha=.3, lty='11', data=.df_div) + 
        # show cell traces
        geom_path(aes(time_sec - 2*3600, gfp_nb, col=factor(id)), data=filter(.df, !discard_top, vertical_top>vertical_cutoff)) +
        geom_text(aes(time_sec - 2*3600, gfp_nb, col=factor(id), label=id), size=2, hjust=0, vjust=1,
                  data=filter(.df, !discard_top) %>% group_by(id) %>% filter(row_number()==1) ) +
        # show all traces in one panel
        geom_path(aes(time_sec - 2*3600, gfp_nb, col=factor(id)), data=filter(.df, !discard_top, vertical_top>vertical_cutoff) %>% mutate(b_rank=-1)) +
        # mask early frames (requires a dummy df!)
        geom_rect(aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf, group=1), fill='white', alpha=.6, data=data.frame(a=1)) +
        scale_colour_periodic_brewer(guide='none') + 
        scale_fill_manual(values=c('glucose'=.fill[1], 'lactose'=.fill[2]), guide='none') +
        scale_x_hours(4) +
        labs(y='total GFP per cell (molecules)') +
        theme(panel.margin = unit(0, "lines"), panel.border=element_blank()),

      fc = ggplot() +
        facet_grid(b_rank~., scales='free_y', as.table=FALSE, labeller=as_labeller(custom_labels)) +
        # facets alternated background
        geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), alpha=.05, data=data.frame(b_rank=seq(0, 6, 2))) +
        # show medium bar
        geom_rect(aes(xmin=t_start*60 - 2*3600, xmax=t_end*60 - 2*3600, ymin=fmin/2, ymax=fmax/2, fill=medium), size=0.2, data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, fmin=.fmin-.frange/2.5, fmax=.fmin-.frange/5)) +
        geom_vline(aes(xintercept=t_start*60 - 2*3600), alpha=.2, size=.5, data=filter(condition_ts, condition==.cond, t_start>0)) +
        geom_text(aes(x=60*(t_start+(t_end-t_start)/2-120), y=f/2, label=medium), col='white', data=filter(condition_ts, condition==.cond) %>% mutate(b_rank=-1, f=.fmin-.frange/3), size=2, hjust=0.5, vjust=0) +
        # show divisions (lty='11' means densely dotted line using hex notation)
        geom_segment(aes(x=time_sec - 2*3600, xend=time_secp - 2*3600, y=gfp_nb/length_um, yend=gfp_nbp/length_um, col=factor(id)), alpha=.3, lty='11', data=.df_div) + 
        # show cell traces
        geom_path(aes(time_sec - 2*3600, gfp_nb/length_um, col=factor(id)), data=filter(.df, !discard_top, vertical_top>vertical_cutoff)) +
        geom_text(aes(time_sec - 2*3600, gfp_nb/length_um, col=factor(id), label=id), size=2, hjust=0, vjust=1,
                  data=filter(.df, !discard_top) %>% group_by(id) %>% filter(row_number()==1) ) +
        # show all traces in one panel
        geom_path(aes(time_sec - 2*3600, gfp_nb/length_um, col=factor(id)), data=filter(.df, !discard_top, vertical_top>vertical_cutoff) %>% mutate(b_rank=-1)) +
        # mask early frames (requires a dummy df!)
        geom_rect(aes(xmin=-Inf, xmax=0, ymin=-Inf, ymax=Inf, group=1), fill='white', alpha=.6, data=data.frame(a=1)) +
        scale_colour_periodic_brewer(guide='none') + 
        scale_fill_manual(values=c('glucose'=.fill[1], 'lactose'=.fill[2]), guide='none') +
        scale_x_hours(4) +
        labs(y='GFP concentration (molecules/µm)') +
        theme(panel.margin = unit(0, "lines"), panel.border=element_blank())
    )
  })(.))


pdf('plots/switch_path_length.pdf', width=12, height=10)
for (i in 1:dim(pls)[1]) {
  if ('l' %in% names(pls[[i, 'pll']]))
    plot(pls[[i, 'pll']] [['l']] + 
           labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
}
dev.off()

pdf('plots/switch_path_fluotot.pdf', width=12, height=10)
for (i in 1:dim(pls)[1]) {
  if ('ft' %in% names(pls[[i, 'pll']]))
  plot(pls[[i, 'pll']] [['ft']] + 
         labs(title=sprintf("%s  pos:%02d  GL:%02d", pls[[i, "date"]], pls[[i, "pos"]], pls[[i, "gl"]])))
}
dev.off()

pdf('plots/switch_path_fluoconc.pdf', width=12, height=10)
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



