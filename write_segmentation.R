library(tidyverse)
library(reticulate)
library(magick)
reticulate::use_condaenv('C:/Users/jeffr/AppData/Local/r-miniconda/envs/cellprofiler')

reticulate::conda_create(getwd())
reticulate::evn
py_bin <- reticulate::conda_list()$python[6]
Sys.setenv(RETICULATE_PYTHON = py_bin)

writeLines(
  paste0("RETICULATE_PYTHON = ", getwd()), 
  paste0(here::here(), ".Rprofile")
)

reticulate::use_condaenv(reticulate::conda_list()$python[2])

reticulate::virtualenv_list()

reticulate::conda_install(
  packages = c("numpy", "matplotlib"),
  envname = getwd()
)


id <- str_extract(list.files('Output/', pattern = '*.tif', full.names = T),
                  '[:alnum:]+\\_Phase\\_[:upper:]\\d{1,2}\\_\\d\\_\\d{2}d\\d{2}h\\d{2}m\\_\\d')

header <- "---
title: 'Bright Field Segmentation'
output: powerpoint_presentation
date: '2023-09-19'
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = F)
library(tidyverse)
library(magick)

setwd('C:/Users/jeffr/Documents/CellProfiler-4.2.6')

object_files <- list.files('Output/', pattern = 'object', full.names = T)
overlay_files <- list.files('Output/', pattern = 'overlay', full.names = T)
image_files <- list.files('Images/', pattern = '*.tif', full.names = T, recursive = T)

measurement_files <- list.files('Output', pattern = '.csv', full.names = T)

image_data <- read_csv(measurement_files[2]) %>% 
  group_by(ImageNumber, Count_RunCellpose, URL_DNA) %>%
  summarise(Confluence = 100*AreaOccupied_AreaOccupied_RunCellpose/AreaOccupied_TotalArea_RunCellpose)

object_data <- read_csv(measurement_files[3]) %>%
  group_by(ImageNumber) %>%
  summarise(Mean_Area = mean(AreaShape_Area),
            Mean_FormFactor = mean(AreaShape_FormFactor),
            Mean_Eccentricity = mean(AreaShape_Eccentricity))

cellpose_results <- left_join(image_data, object_data) 

```
"


image_loop <- unlist(map(id, ~paste0(
  "
## ", .x,"

```{r}


image_temp <- image_read(str_subset(image_files, '", .x, "'))

overlay_temp <- image_read(str_subset(overlay_files, '", .x, "'))
object_temp <- image_transparent(image_read(str_subset(object_files, '", .x, "')), 'black')

results_temp <- cellpose_results %>%
  mutate(URL_DNA = str_extract(URL_DNA, '", .x, "')) %>%
  dplyr::filter(URL_DNA == '", .x, "') %>%
  tibble::column_to_rownames('URL_DNA') %>%
  t() %>%
  data.frame() %>%
  mutate_if(is.numeric, function(r) round(r, digits = 2))


mont_temp <- image_montage(c(image_temp, object_temp, overlay_temp), tile = '2x2',  geometry = 'x400+1+1')

str_table <- str_c(paste(rownames(results_temp), ':', results_temp[, 1]), collapse = '\n')

image_annotate(mont_temp, str_table, size = 30, location = '+600+450')

```")))

#writeLines(str_c(header, str_c(image_loop, collapse = '\n'), collapse = '\n\n'), 'segmentation.rmd')
#
#allImagesFiles <- list.files('Images', pattern = '.*tif', recursive = T, full.names = T)
#
#image_ids <- str_extract(allImagesFiles, '[:alnum:]+\\_Phase\\_[:upper:]\\d{1,2}\\_\\d\\_\\d{2}d\\d{2}h\\d{2}m') 
#
#
#for(i in image_ids){
#  
#  file_temp <- str_subset(allImagesFiles, i)
#  
#  image_group <- image_read(file_temp[c(1, 3, 2, 4)])
#  
#  mont_temp <- image_montage(image_group, tile = '2x2',  geometry = 'x400+0+0')
#  
#  image_write(mont_temp, path = paste0('Images/livecell_complete/', i, '.tif'), format = 'tif')
#  
#  
#  
#}


