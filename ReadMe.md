

This repository contains an Rstudio project. 
The R environment used for this project is managed using `renv` (running `renv::init()` should restore all necessary package, but you need to make sure that you use an appropriate R version — at best the same as described in the file `renv.lock`). 
Learn more about collaborating with `renv` at https://rstudio.github.io/renv/articles/collaborating.html#collaborating-with-renv.

Run `MoM_lacInduction.R` to load the data and render the analysis files to html.
Note that calling `render()` or `render_site()` from the command line allows to execute the function in the global env() (hence inheriting existing variables and keeping newly created ones).

These scripts rely heavily on `multidplyr`...

## Rmardown rendering
Designed as a Rmarkdown "site". hence requires rmarkdown ≥ 1.0

`index.Rmd` file is required for site_render() to execute.

in _site.yml, `exclude: ["*"]` is required to prevent all subdirectories to be copies (all the more so as symlinks are followed!)

Here is an example of the minimal YAML header to put in each Rmarkdown file.
NB: date syntax from http://stackoverflow.com/questions/23449319

```
---
title: "My relevant title"
author: Thomas Julou
date: "`r format(Sys.time(), '%d %B, %Y')`"
---
```

