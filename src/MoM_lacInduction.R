
# SET VARIABLES & ENVIRONMENT ####

# SET VARIABLES
dl <- 0.065                # pixel size (µm)
vertical_cutoff <- 4 / dl  # after it touched this coordinate a cell is discarded
min_growth_rate <- 4e-5    # growth rate threshold to discard non growing cells before the switch (sec-1)
use_eriks_params <- TRUE

data2preproc_dir <- function(.d)
  str_match(.d, '20\\d{6}') %>% na.omit %>% as.character %>% 
  file.path('.', 'preproc', .)
data2preproc_file <- function(.f)
  basename(.f) %>% sub("ExportedCellStats_", "", .) %>% 
  file_path_sans_ext %>% paste0("_frames.txt")
data2preproc <- function(.f)
  file.path(data2preproc_dir(.f), data2preproc_file(.f))

# scale_colour_discrete <- scale_colour_viridis_d
# scale_colour_continuous <- scale_colour_viridis_c

# EXPERIMENTAL CONDITIONS AND DATA PATHS
myconditions <- list(
  # CONSTANT CONDITIONS
  list(condition='mg1655', duration=c(720, 720), dt=180, medium=c('glucose', 'lactose'),
       paths=c("./data_MoM_ms/MG1655_glu_lac")),
  list(condition='glucose', duration=1560, dt=180, medium='glucose',
       paths=c("./data_MoM_ms/glucose")),
  list(condition='lactose', duration=1560, dt=180, medium='lactose',
       paths=c("./data_MoM_ms/lactose")),
  
  # GLC / LAC SWITCH (MEMORY)
  list(condition='switch_04h',
       duration=c(360, 240, 240, 240, 240, 240), dt=180,
       medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_MoM_ms/glu_lac_switch")),
  list(condition='switch_06h',
       duration=c(360, 240, 360, 240, 360, 240), dt=180, 
       medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20151204/20151204_switch6h_curated")),
  list(condition='switch_08h', 
       duration=c(360, 240, 480, 240, 480, 240), dt=180, 
       medium=c('glucose', 'lactose', 'glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20151218/20151218_switch8h_curated", "./data_thomas/20180206/20180206_glu_lac_switch8h_curated/")),
  list(condition='switch_12h', duration=c(360, 240, 720, 360), dt=180, 
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20180207/20180207_glu_lac_switch12h_curated/", "./data_thomas/20180216/20180216_glu_lac_switch12h_curated/")),
  list(condition='switch_16h', duration=c(360, 240, 960, 360), dt=180, 
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20160912/20160912_curated")),
  list(condition='switch_20h', duration=c(360, 240, 1200, 360), dt=180, 
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20161014/20161014_curated")),
  list(condition='switch_24h', duration=c(360, 240, 1440, 360), dt=180, 
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20161007/20161007_curated", "./data_thomas/20180313/20180313_glu_lac_switch24h_curated/")),
  # paths=c("./data_thomas/20161007/20161007_curated")),
  list(condition='switch_12h_old', duration=c(240, 240, 720, 360), dt=180,
       medium=c('glucose', 'lactose', 'glucose', 'lactose'),
       paths=c("./data_thomas/20160526/20160526_curated")),
  
  # GLC / LAC SWITCH (CONTROLS)
  list(condition='switch_lactose_priming', duration=c(240, 240), dt=180,
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20161212/20161212_curated")), 
  list(condition='switch_late', duration=c(1560, 480), dt=180, 
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20161021/20161021_curated", "./data_thomas/20170108/20170108_glu_lac_ctrl16h_curated/",
               "./data_thomas/20180516/20180516_glu_lac_ctrl26h_curated/", "./data_thomas/20180615/20180615_glu_lac_ctrl26h_curated/")), 
  # list(condition='switch_long_lac',
  #      duration=c(360, 1800), dt=180,
  #      medium=c('glucose', 'lactose'),
  #      paths=c("./data_thomas/20161104/20161104_curated")),
  # list(condition='switch_long_lac_hiExpr',
  #      duration=c(360, 1800, 240), dt=180,
  #      medium=c('glucose', 'lactose', 'glucose'),
  #      paths=c("./data_thomas/20170926/20170926_glu_lac_30h_hiExpr_curated/")),
  
  # LAC PROTEINS THRESHOLD
  list(condition='switch_pre_lac001_6h', duration=c(360, 1080), dt=c(180, 360), 
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20190605/20190605_glu_lowLac_curated/")),
  list(condition='switch_pre_lac001_10h', duration=c(600, 1080), dt=c(180, 360), 
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20190614/20190614_glu_lowLac_curated/")),
  
  # CATABOLITE REPRESSION
  list(condition='switch_gly_lac', duration=c(480, 360), dt=180, 
       medium=c('glycerol', 'lactose'),
       paths=c("./data_thomas/20170919/20170919_glyc_lac_curated/", "./data_thomas/20170920/20170920_glyc_lac_curated/")),
  
  # PERTURBATION OF REPRESSED EXPRESSION
  list(condition='switch_withIPTG1uM',
       duration=c(360, 240), dt=180,
       medium=c('glucose+IPTG', 'lactose+IPTG'),
       paths=c("./data_thomas/20161130/20161130_switch_IPTG1uM_curated")),
  list(condition='switch_withIPTG5uM',
       duration=c(360, 240), dt=180,
       medium=c('glucose+IPTG', 'lactose+IPTG'),
       paths=c("./data_thomas/20161207/20161207_curated", "./data_thomas/20180122/20180122_glu_lac_IPTG5uM_curated/")),
  list(condition='switch_lacIoe',
       duration=c(360, 720), dt=180,
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20180116/20180116_glu_lac_lacIoe_curated", "./data_thomas/20180123/20180123_lacIoe_curated/",
               "./data_thomas/20180214/20180214_lacIoe_curated/")),
  list(condition='switch_lacIoe_withIPTG10uM',
       duration=c(360, 480), dt=180,
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20180119/20180119_lacIoe_IPTG10uM_curated/")),
  list(condition='switch_lacIPTG5uM',
       duration=c(360, 360), dt=180,
       medium=c('glucose', 'lactose+IPTG'),
       paths=c("./data_thomas/20180316/20180316_glu_lacIPTG5uM_curated/")),
  list(condition='switch_preIPTG5uM',
       duration=c(360, 240), dt=180,
       medium=c('glucose+IPTG', 'lactose'),
       paths=c("./data_thomas/20180514/20180514_gluIPTG5uM_lac_curated/", "./data_thomas/20180531/20180531_gluIPTG5uM_lac_curated/")),
  list(condition='switch_lacIoe_preIPTG10uM',
       duration=c(360, 240), dt=180,
       medium=c('glucose+IPTG', 'lactose'),
       paths=c("./data_thomas/20180604/20180604_gluIPTG10uM_lac_lacIoe_curated/")),
  
  list(condition='switch_glcLac_lac', duration=c(720, 240, 960, 240), dt=180,
       medium=c('glucose+lac', 'lactose', 'glucose+lac', 'lactose'),
       paths=c("./data_thomas/20171114/20171114_glcLac_lac_switch_curated/", "./data_thomas/20180108/20180108_gluLac_lac_curated/", # only 12+4h
               "./data_thomas/20180606/20180606_gluLac_lac_switch16h_curated/")), # stopped before last switch
  list(condition='switch_ramp15min',
       duration=c(585, 255, 240), dt=180,
       medium=c('glucose', 'lactose', 'glucose'),
       paths=c("./data_thomas/20170901/20170901_switch_glu_lac_ramp15min_curated/")),
  list(condition='switch_ramp40min',
       duration=c(414, 246), dt=180,
       medium=c('glucose', 'lactose'),
       paths=c("./data_thomas/20171121/20171121_glu_lac_ramp40min_curated/", "./data_thomas/20180319/20180319_glu_lac_ramp40min_curated/"))
  
)

# SET ENVIRONMENT
# install.packages(c('tools', 'devtools', 'here', 'tidyverse', 'RcppArmadillo', 'svglite'))
# devtools::install_github(c('julou/ggCustomTJ', 'hadley/multidplyr'))
# devtools::install_github('vanNimwegenLab/vngMoM', auth_token='xxx')
mylibs <- c('here', 'tidyverse', 'cowplot', 'tools', 'RcppArmadillo', 'vngMoM', 'ggCustomTJ')
invisible( suppressPackageStartupMessages( # don't use %>% before loading dplyr
  lapply(mylibs, library, character.only=TRUE) ))
library(svglite)

dir.create(here("slogs"), showWarnings=FALSE) # create a directory to store logs from the queue

# set a parallel environment to run multidplyr
library(multidplyr)
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  create_cluster() %>% cluster_library(mylibs)       # load libraries on each core

# LOAD MoMA DATA ####
# find raw data files from myconditions and store them in a dataframe
myfiles <- myconditions %>% 
  # convert the relevant list items to a dataframe
  lapply(function(.l) .l[ which(names(.l) %in% c("condition", "paths"))] %>% 
           as.data.frame(stringsAsFactors=FALSE) ) %>% 
  do.call(rbind, .) %>% 
  rename(data_path=paths) %>%
  # for each path, find all files matched by the pattern .*\\d+\\.csv (e.g. *20151023.csv)
  group_by(condition, data_path) %>% 
  do((function(.df)
    # list.files(.df$data_path, ".*\\d+\\.csv", recursive=TRUE, full.names=TRUE) %>% 
    find.files(.df$data_path, "ExportedCellStats_*.csv") %>% 
      data.frame(path=., stringsAsFactors=FALSE) )(.))  

#  create condition_acq_times (describing acquisition times and temporal change of each condition) from myconditions
condition_acq_times <- myconditions %>% 
  # convert the relevant list items to a dataframe
  lapply(function(.l) .l[ - which(names(.l) == "paths")] %>% as_tibble ) %>% 
  bind_rows() %>% 
  group_by(condition) %>% 
  mutate(m_start=cumsum(c(0, duration[-(length(duration))])) * 60,
         m_end=cumsum(duration) * 60 - 1e-5, duration=NULL,
         m_cycle=value_occurence_index(medium), 
  ) %>% 
  group_by(condition, m_start) %>% 
  mutate(data=list(data.frame(time_sec=seq(m_start, m_end, dt)) )) %>% unnest() %>% 
  group_by(condition) %>% 
  mutate(frame=as.numeric(order(time_sec)-1))

# load perl scripts output to dataframes (using parallel dplyr)
nc <- min(nrow(myfiles), length(mycluster)) # this is a dirty hack because multidplyr crashes with less shards than cores
myframes <- myfiles %>% 
  # process exported files on the cluster if required (otherwise return the list of paths)
  ungroup %>% 
  mutate(ppath=process_moma_data(path, .data2preproc=data2preproc, .frames_pl_script="get_size_and_fluo_multich.pl", #.skip=TRUE,
                                 .qsub_name="MMex_pl", .force=FALSE, .skip=TRUE) ) %>% 
  filter(!is.na(ppath)) %>% 
  # load perl scripts output to dataframes (in parallel, using multidplyr)
  partition(condition, path, cluster=mycluster[1:nc] ) %>%
  # group_by(condition, path) %>% # non-parallel alternative
  do((function(.df){
    # browser()
    # print(.df$ppath)
    parse_frames_stats(.df$ppath)
  })(.)) %>%
  collect() %>%
  ungroup %>%
  extract(path, c("date", "pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE, convert=TRUE) %>%
  # propagate time and medium info
  left_join(condition_acq_times) %>%
  group_by(path) %>% 
  mutate(m_end=ifelse(m_end > max(time_sec), max(time_sec), m_end)) %>% 
  # fix end_type for pruned cells
  mutate(discard_start=(time_sec < 2*3600), length_um=length_pixel*dl) %>%
  group_by(path) %>%
  mutate(ndgt=compute_daughters_numbers(cid)) %>%
  mutate(end_type_moma=end_type,
         end_type=ifelse(ndgt==0, "lost", "weird"),
         end_type=ifelse(ndgt==2, "div", end_type)) %>%
  # remove frames after touching the exit
  group_by(id, add=TRUE) %>%
  mutate(discard_top=which_touch_exit(vertical_top, vertical_cutoff)) %>%
  mutate(discard_top=ifelse(discard_start, FALSE, discard_top)) %>% # not in the preexpt step (2h)
  mutate(end_type=ifelse(any(discard_top), 'touchtop', end_type)) %>% # update end_type to exit for cells which have touched the vertical cutoff
  # remove daughters of cells that touched the exit
  group_by(path) %>%
  mutate(discard_top=which_to_progeny(discard_top, cid)) %>%
  # append useful variables (per cell)
  mutate(
    cell=paste(date, pos, gl, id, sep='.'),
    ugen=paste(date, pos, gl, genealogy, sep='.'),
    gl_id=paste(date, pos, gl, sep='.'),
    vertical_center=(vertical_bottom + vertical_top)/2,
    mstep=paste(medium, m_cycle, sep='.'),
    strain='ASC662', strain=ifelse(condition=='mg1655', 'MG1655', strain),
    strain=ifelse(condition=='switch_long_lac_hiExpr', 'MG1655_pHi-GFP', strain),
    strain=ifelse(condition=='switch_∆lacA', 'AB460', strain),
  ) %>%
  group_by(date, pos, gl, id) %>%
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         b_rank=round(mean(total_cell_in_lane - cell_num_in_lane)),
         length_raw=(vertical_bottom-vertical_top)*dl,
         length_erik=length_um, length_um=length_raw)

# CONVERT FLUO UNITS ####
if (use_eriks_params)
  autofluo_predict <- function(.h) .h * 422.8
myframes <- myframes %>% ungroup() %>% mutate(
  fluo_amplitude=fluo_amplitude * ifelse(date==20190605 & time_sec/3600 >= 6, 5, 1),
  fluo_amplitude=fluo_amplitude * ifelse(date==20190614 & time_sec/3600 >= 10, 5, 1),
  
  fluogfp_amplitude = fluo_amplitude - autofluo_predict(length_um)
)


fp_per_oligomers <- 4 # lacZ is tetrameric
if (use_eriks_params)
  fp_per_dn <- 0.0361 * fp_per_oligomers
myframes <- myframes %>% ungroup() %>% 
  # convert to gfp units (after subtracting autofluorescence)
  mutate(gfp_nb = fluogfp_amplitude * fp_per_dn ) %>% 
  group_by(date, pos, gl, id)

# experiments to be discarded as identified by controls
discarded_dates <- c(
  # slow growth in initial condition
  20151218, # switch_08h
  20161021, # switch_late
  20160526, # switch_12h_old (long gr tail + bad r2)
  20170108, # switch_late
  20170901, # switch_ramp 15'
  20180313, # switch_24h only 10GLs analysed
  # other reasons
  20180123, # switch_lacIoe (weird lag distrib, but all longer than without overexpressing LacI)
  20180615,  # weird late switch control (all switch fast!)
  20180606  # weird late switch control (all switch fast!)
  
  # TODO: check whether the heritability at first switch can be used as a criteria
)


# RENDER ANALYSIS FILES ####
# calling `render()` or `render_site()` from the command line allows to execute the function 
# in the global envt (hence inheriting existing variables and keeping the newly created ones)...

rename_conds <- function (.str) {
# TODO: check NAs for all conditions
  .labels <- .str
  # .labels[.labels=='1'] <- 'switch:1'
  .labels <- str_replace(.labels, "switch_glcLac_lac", "glc+lac > lac")
  .labels <- str_replace(.labels, "switch_gly_lac", "glyc > lac")
  .labels <- str_replace(.labels, "switch_lactulose", "> lactulose")
  .labels <- str_replace(.labels, "switch_glycerol_", "> glyc ")
  # .labels <- str_replace(.labels, "glycerol", "glyc")
  .labels <- str_replace(.labels, "switch_0?(\\d+)h", "memory \\1h")
  .labels <- str_replace(.labels, "_old", " (old)")
  .labels <- str_replace(.labels, "^switch_", "")
  .labels <- str_replace(.labels, "_", " ")
  return(.labels)
}
lac_lags_label <- expression(paste(italic('lac'), ' induction lag (min)'))


myscales <- list()
myplots <- list()
mytables <- list()
theme_set(theme_cowplot() + theme(strip.text.x=element_text(margin=margin(t=1, b=2)),
                                  strip.text.y=element_text(margin=margin(l=2, r=1))) )
theme_cowplot_legend_inset <- function(.rel=0.7) theme(legend.title=element_text(size=rel(.rel), face='bold'),
                                                       legend.text=element_text(size=rel(.rel)))
knitr::opts_chunk$set(
  echo=FALSE, message=FALSE, warning=FALSE,
  dev="svglite"
)
# rmarkdown::clean_site()

# render control plots of each GC
# source('./src/MoM_lacInduction_GCplots.R')

# rmarkdown::render_site('./src/MoM_lacInduction_GFP_Estimation.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_Constant_Envts.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_Lags_Estimation.Rmd')
# rmarkdown::render_site('./src/MoM_lacInduction_Heritability_Estimation.Rmd')

# DISCARD SOME DATASETS
rmarkdown::render_site('./src/MoM_lacInduction_Controls.Rmd')

rmarkdown::render_site('./src/index.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_Native.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_PerturbRepressed.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_Sensitivity.Rmd')


# RENDER ARTICLE FILES ####

knitr::opts_chunk$set(
  echo=FALSE, message=FALSE, warning=FALSE
  )

myfigs <- list()
source('./src/MoM_lacInduction_Figs.R')
source('./src/MoM_lacInduction_FigsSI.R')


# rmarkdown::render("./manuscript/MoM_lacInduction_ms.Rmd", 
#                   bookdown::pdf_book(base_format=bookdown::pdf_book, fig_width=4.75, fig_height=2.25*14/9), 
#                   # bookdown::pdf_book(base_format=rticles::plos_article, fig_width=4.75, fig_height=2.25*14/9), 
#                   clean=FALSE)
# 
# rmarkdown::render("./manuscript/MoM_lacInduction_SM.Rmd", 
#                   bookdown::pdf_book(base_format=rticles::plos_article, fig_width=4.75, fig_height=2.25*14/9), 
#                   clean=FALSE)
# 
# # resolve references
# (function(.dir, .files) {
#   setwd("manuscript")
#   tinytex::pdflatex("MoM_lacInduction_ms.tex", clean=FALSE)
#   tinytex::pdflatex("MoM_lacInduction_SM.tex", clean=FALSE)
#   
#   tinytex::pdflatex("MoM_lacInduction_ms.tex", clean=TRUE)
#   tinytex::pdflatex("MoM_lacInduction_SM.tex", clean=TRUE)
#   setwd('..')
# })()

