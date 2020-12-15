# EVALUATION PARAMETERS ####
load_from_Rdata <- TRUE


# SET ENVIRONMENT ####
options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/centos7/latest"))
#install.packages(c('remotes', 'here', 'tidyverse', 'RcppArmadillo', 'svglite'))
# remotes::install_github(c('julou/ggCustomTJ', 'hadley/multidplyr'))
# remotes::install_github('vanNimwegenLab/vngMoM', auth_token='xxx')
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(cowplot)
  # library(tools)
  library(RcppArmadillo)
  library(vngMoM)
  library(ggCustomTJ)
})

dir.create(here("slogs"), showWarnings=FALSE) # create a directory to store logs from the queue

theme_set(theme_cowplot() + theme(title = element_text(size = rel(1/1.14)),
                                  strip.text.x=element_text(margin=margin(t=1, b=2)),
                                  strip.text.y=element_text(margin=margin(l=2, r=1))) )
theme_cowplot_legend_inset <- function(.rel=0.7) theme(legend.title=element_text(size=rel(.rel), face='bold'),
                                                       legend.text=element_text(size=rel(.rel)))

# set a parallel environment to run multidplyr
library(multidplyr)
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  create_cluster() %>% cluster_library( # load currently loaded packages on each core
    names(sessionInfo()$otherPkgs)
  )


if (!load_from_Rdata) {
    source("./src/MoM_lacInduction_dataLoad.R")
  } else {
    # select the RData file in the dialog window
    load(file.choose())
  }


# RENDER ANALYSIS FILES ####
# calling `render()` or `render_site()` from the command line allows to execute the function 
# in the global envt (hence inheriting existing variables and keeping the newly created ones)...
filter_article_ds <- function(.df, .filter=TRUE) {
  if (!.filter) return(.df)
  filter(.df,
         condition %in% c("mg1655", "glucose", "lactose", "switch_glcLac_lac", "switch_gly_lac", 
                          "switch_late", "switch_ramp40min", "switch_lacIoe") |
           # "switch_glycerol_",
           # "glycerol"
           str_detect(condition, "switch_lac001") | 
           str_detect(condition, "switch_0?(\\d+)h") | 
           str_detect(condition, "preIPTG")
  ) }


rename_conds <- function (.str) {
# TODO: check NAs for all conditions
  .labels <- .str
  # .labels[.labels=='1'] <- 'switch:1'
  # .labels <- str_replace(.labels, "switch_glcLac_lac", "glc+lac > lac")
  .labels <- str_replace(.labels, "switch_gly_lac", "glyc > lac")
  .labels <- str_replace(.labels, "switch_glycerol_", "> glyc ")
  # .labels <- str_replace(.labels, "glycerol", "glyc")
  .labels <- str_replace(.labels, "switch_0?(\\d+)h", "memory \\1h")
  .labels <- str_replace(.labels, "switch_lac001.*", "low [lactose]")
  .labels <- str_replace(.labels, "_old", " (old)")
  .labels <- str_replace(.labels, "^switch_", "")
  .labels <- str_replace(.labels, "_", " ")
  return(.labels)
}
lac_lags_label <- expression(paste(italic('lac'), ' induction lag (min)'))


myplots <- list()
mytables <- list()

library(svglite)
knitr::opts_chunk$set(
  echo=TRUE, message=FALSE, warning=FALSE,
  dev="svglite"
  # dev.args=list(),
)
# rmarkdown::clean_site()

# render control plots of each GC
# source('./src/MoM_lacInduction_GCplots.Rd')
rmarkdown::render_site('./src/MoM_lacInduction_Constant_Envts.Rmd') # to produce comparison values
rmarkdown::render_site('./src/MoM_lacInduction_Lags_Estimation.Rmd')

# DISCARD SOME DATASETS
rmarkdown::render_site('./src/MoM_lacInduction_Controls.Rmd')

rmarkdown::render_site('./src/index.Rmd')

rmarkdown::render_site('./src/MoM_lacInduction_Native.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_PerturbRepressed.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_lowLactose.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_PopLagSimul.Rmd')

rmarkdown::render_site('./src/MoM_lacInduction_FLIM.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_diauxieGrowthCurves.Rmd')
rmarkdown::render_site('./src/MoM_lacInduction_QMS.Rmd')


# RENDER ARTICLE FIGURES ####

knitr::opts_chunk$set(
  echo=FALSE, message=FALSE, warning=FALSE
  )

myfigs <- list()
source('./src/MoM_lacInduction_Figs.R')
source('./src/MoM_lacInduction_FigsSI.R')

# save.image(".xz.RData", compress="xz")
