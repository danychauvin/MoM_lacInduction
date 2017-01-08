
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
data2preproc_dir <- function(.d)
  str_match(.d, '20\\d{6}') %>% na.omit %>% as.character %>% 
  file.path('.', 'preproc', .)
data2preproc_file <- function(.f)
  basename(.f) %>% sub("ExportedCellStats_", "", .) %>% 
  file_path_sans_ext %>% paste0("_frames.txt")
data2preproc <- function(.f)
  file.path(data2preproc_dir(.f), data2preproc_file(.f))
  
# # local paths
# proj_path <- "~/Downloads/thursday/MoM_Switch/"
# r_scripts_path <- c("/Users/julou/Documents/Biozentrum/Projects/vngWetLabR/mother_machine",
#                     "/Users/julou/Documents/Biozentrum/Projects/vngWetLabR/ggplot")


# EXPERIMENTAL CONDITIONS AND DATA PATHS
myconditions <- list(
  list(condition='mg1655', duration=c(720, 720), medium=c('glucose', 'lactose'),
       paths=c("./data_MoM_ms/MG1655_glu_lac")),
  list(condition='glucose', duration=1560, medium='glucose',
       paths=c("./data_MoM_ms/glucose")),
  list(condition='lactose', duration=1560, medium='lactose',
       paths=c("./data_MoM_ms/lactose")),
  # list(condition='lactose_lowillum', duration=1560, medium='lactose',
  #      paths=c("./data_copy/20160318", "./data_copy/20160427")),
  list(condition='switch_long_lac',
       duration=c(360, 1800),
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20161104/20161104_curated")),
  list(condition='switch_lac_withIPTG1uM',
       duration=c(360, 240),
       medium=c('glucose+IPTG', 'lactose+IPTG'),
       paths=c("./data_thomas/20161130/20161130_switch_IPTG1uM_curated")),
  list(condition='switch_lac_withIPTG5uM',
       duration=c(360, 240),
       medium=c('glucose+IPTG', 'lactose+IPTG'),
       paths=c("./data_thomas/20161207/20161207_curated")),
  list(condition='switch_04h',
       duration=c(360, 240, 240, 240, 240, 240),
       medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_MoM_ms/glu_lac_switch")),
  list(condition='switch_06h',
       duration=c(360, 240, 360, 240, 360, 240), 
       medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20151204/20151204_switch6h_curated")),
  list(condition='switch_08h', 
       duration=c(360, 240, 480, 240, 480, 240), 
       medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20151218/20151218_switch8h_curated")),
  list(condition='switch_12h', duration=c(240, 240, 720, 360), 
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20160526/20160526_curated")),
  list(condition='switch_16h', duration=c(360, 240, 960, 360), 
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20160912/20160912_curated")),
  list(condition='switch_24h', duration=c(360, 240, 1440, 360), 
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20161007/20161007_curated")),
  list(condition='switch_iptg', duration=c(360, 240, 240, 240, 240, 240),
       medium=c('glucose', 'lactose+IPTG', 'glucose', 'lactose+IPTG', 'glucose', 'lactose+IPTG'),
       paths=c("./data_thomas/20151207/20151207_switch_iptg_curated")) #,
  # list(condition='switch_m9',
  #      duration=c(360, 120, 360, 720, 360, 240), 
  #      medium=c('glucose', 'M9', 'glucose', 'M9', 'glucose', 'M9'),
  #      paths=c("./data_thomas/20151221")) 
)

# SET ENVIRONMENT
invisible(sapply(
  list.files(r_scripts_path, pattern="\\.[Rr]$", full.names=TRUE, ignore.case=TRUE), 
  source, .GlobalEnv))

setwd(proj_path)
dir.create("qlogs", showWarnings=FALSE) # create a directory to store logs from the queue

# set a parallel environment to run multidplyr
library(multidplyr)
numCores <- min(30, parallel::detectCores()-1) # do not use more than 30 cores
mycluster <- create_cluster(numCores)
for (.l in mylibs) { mycluster %>% cluster_library(.l) } # load all libraries on each core


# LOAD MoMA DATA ####
# find raw data files from myconditions and store them in a dataframe
myfiles <- myconditions %>% 
  # convert the relevant list items to a dataframe
  lapply(function(.l) .l[ - which(names(.l) %in% c("duration", "medium"))] %>% 
           as.data.frame(stringsAsFactors=FALSE) ) %>% 
  do.call(rbind, .) %>% 
  rename(data_path=paths) %>%
  # for each path, find all files matched by the pattern .*\\d+\\.csv (e.g. *20151023.csv)
  group_by(condition, data_path) %>% 
  do((function(.df)
    # list.files(.df$data_path, ".*\\d+\\.csv", recursive=TRUE, full.names=TRUE) %>% 
    find.files(.df$data_path, "ExportedCellStats_*.csv") %>% 
      data.frame(path=., stringsAsFactors=FALSE) )(.))  

# create condition_ts (describing each temporal change of each condition) from myconditions 
condition_ts <- myconditions %>% 
  # convert the relevant list items to a dataframe
  lapply(function(.l) .l[ - which(names(.l) == "paths")] %>% as.data.frame ) %>% 
  do.call(rbind, .) %>% 
  # add useful variables for each condition step
  group_by(condition) %>% 
  mutate(t_start=cumsum(c(0, duration[-(length(duration))])) * 60,
         t_end=cumsum(duration) * 60 - 1e-5,
         m_cycle=value_occurence_index(medium))

# load perl scripts output to dataframes (using parallel dplyr)
nc <- min(nrow(myfiles), length(mycluster)) # this is a dirty hack because multidplyr crashes with less shards than cores
myframes <- myfiles %>% 
  # process exported files on the cluster if required (otherwise return the list of paths)
  ungroup %>% 
  mutate(ppath=process_moma_data(path, .data2preproc=data2preproc, 
                                 .scripts_path=perl_scripts_path, .qsub_name="MMex_pl", .force=FALSE) ) %>% 
  # load perl scripts output to dataframes (in parallel, using multidplyr)
  partition(condition, path, cluster=mycluster[1:nc] %>%
              cluster_assign_func(parse_frames_stats, compute_genealogy, which_touch_exit, which_to_progeny) %>%
              cluster_assign_obj(dt, dl, vertical_cutoff)) %>%
  # group_by(condition, path) %>% # non-parallel alternative
  do((function(.df){
    parse_frames_stats(.df$ppath) %>% 
      mutate(time_sec=frame*dt*60, length_um=length_pixel*dl,
             discard_start=(time_sec < 2*3600)) %>%
      # remove frames after touching the exit
      group_by(id) %>%
      mutate(discard_top=which_touch_exit(vertical_top, vertical_cutoff)) %>%
      mutate(discard_top=ifelse(discard_start, FALSE, discard_top)) %>% # not in the preexpt step (2h)
      mutate(end_type=ifelse(any(discard_top), 'exit', end_type)) %>% # update end_type to exit for cells which have touched the vertical cutoff
      # remove daughters of cells that touched the exit
      ungroup %>%
      mutate(discard_top=which_to_progeny(discard_top, cid))
  })(.)) %>%
  collect() %>% 
  # append useful variables (per row)
  ungroup %>% 
  extract(path, c("date", "pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE, convert=TRUE) %>%
  # mutate_at(vars(date, pos, gl), factor) %>% # convert to factors (ad libidum)
  mutate(cell=paste(date, pos, gl, id, sep='.'),
         gl_id=paste(date, pos, gl, sep='.'),
         vertical_center=(vertical_bottom + vertical_top)/2) %>% 
  # propagate medium info
  group_by(condition) %>% 
  do((function(.df){
    .ts <- filter(condition_ts, condition==unique(.df$condition))
    .idx <- find_unique_interval(.df$time_sec, .ts$t_start, .ts$t_end)
    mutate(.df, medium=.ts$medium[.idx], m_start=.ts$t_start[.idx], m_end=.ts$t_end[.idx], m_cycle=.ts$m_cycle[.idx])
  })(.)) %>% 
  # sort and append useful variables (per cell)
  arrange(date, pos, gl, id, frame) %>%  # sort data after `partition()`
  group_by(date, pos, gl, id) %>%
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         b_rank=round(mean(total_cell_in_lane - cell_num_in_lane)))



# RENDER ANALYSIS FILES ####
# calling `render()` or `render_site()` from the command line allows to execute the function 
# in the global envt (hence inheriting existing variables and keeping the newly created ones)...

# render control plots of each GC
# source('MoM_Switch_GCplots.R')

knitr::opts_chunk$set(echo=FALSE, message=FALSE, warning=FALSE)
# rmarkdown::clean_site()

# rmarkdown::render_site('index.Rmd')
rmarkdown::render_site('MoM_Switch_GFP_Estimation.Rmd')
# rmarkdown::render_site('MoM_Switch_Constant_Envts.Rmd')
rmarkdown::render_site('MoM_Switch_Switching_Envts.Rmd')
# rmarkdown::render_site('MoM_Switch_FCM_Comparison.Rmd')



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
