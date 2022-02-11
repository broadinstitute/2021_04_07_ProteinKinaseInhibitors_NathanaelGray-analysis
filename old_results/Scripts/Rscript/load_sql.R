#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = T)

extends <- methods::extends


library(tidyverse)
library(cytominer)
library(stringr)
library(readr)
library(reshape2)


print("script started")

top_dir <-'/home/ubuntu/bucket/projects/'
proj_dir <- '2018_11_20_GeneCpdFollowup'
batch.name <- '2018_11_20_Batch1'


output <- "/home/ubuntu/bucket/projects/2018_11_20_GeneCpdFollowup/workspace/backend/2018_11_20_Batch1/"


plate.list <- c("BR00100032", "BR00100037")


sql.path <- NULL
for (pl in seq_along(plate.list)) {
  sql.path[pl] <- as.vector(paste0("/home/ubuntu/bucket/projects/2018_11_20_GeneCpdFollowup/workspace/backend/", batch.name, "/", plate.list[pl], "/", plate.list[pl], ".sqlite"))
}

print("Step1")

# reading sqlite
read_sql<- function(sql.path) {
  db <- DBI::dbConnect(RSQLite::SQLite(), sql.path)
  RSQLite::initExtension(db)
  
  image <- RSQLite::dbReadTable(conn = db, "Image")
  cells <- RSQLite::dbReadTable(conn = db, "Cells")
  nuclei <- RSQLite::dbReadTable(conn = db, "Nuclei")
  cytoplasm <- RSQLite::dbReadTable(conn = db, "Cytoplasm")
  
  dt <- cells %>%
    left_join(cytoplasm, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
    left_join(nuclei, by = c("TableNumber", "ImageNumber", "ObjectNumber")) %>%
    left_join(image, by = c("TableNumber", "ImageNumber"))
  
  return(dt)
  
}

print("Step2")

for (i in 1:length(sql.path)) {
  
  #Extracting metadata
  
  # reading sqlite file
  sql_data <- as.data.frame(lapply(sql.path[i], read_sql))
  
  selected_wells <- c('L02','N19','H10','O19','N17','A20','N22','D09','O05','O10','C05','A21',
                      'P08','O23','B06','N10','E03','E04','E05','E07','E08',
                      'E09','E16','E17','E18','E20','E21','E22','F03','F07','F16','F20','K05','K08',
                      'K18','K21','L03','L05','L06','L07','L08','L09','L16','L18','L19','L20','L21','L22')
  
  print("Step4")
  
  sql_data <- sql_data %>% dplyr::filter(Metadata_Well %in% selected_wells)
  
  
  
  #dmso <- merge(sql_data, meta, by.x="Image_Metadata_Well", by.y="Metadata_Well", all.x = TRUE)
  
  
  print("Step5") 
  
  
  
  readr::write_csv(sql_data, paste0(output,plate.list[i], "/", plate.list[i], "_selected_PIN_DCLK_wells", ".csv"))
  print("Successfully executed") 
}
