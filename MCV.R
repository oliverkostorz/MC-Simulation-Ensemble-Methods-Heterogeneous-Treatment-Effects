rm(list = ls(all.names = TRUE))
set.seed(0815)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pacman::p_load()

#Import custom functions
source('functions.R')


#Theme for graphs
own_theme <- theme_bw() +
  theme(legend.position = 'none', axis.title.x = element_blank(),
        axis.title.y = element_blank(), text = element_text(size = 20))