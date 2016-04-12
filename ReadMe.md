

use `MoM_Switch.R` to load the data and render the analysis files to html.
calling `render()` or `render_site()` from the command line allows to execute the function in the global env() (hence inheriting existing variables and keeping newly created ones).

these scripts rely heavily on `multidplyr`...

## Rmardown rendering
designed as a Rmarkdown "site". hence requires rmarkdown dev version (as of April 2016).

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

