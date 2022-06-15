# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
IN_FILES = list("/hps/nobackup/birney/users/ian/MIKK_F2_tracking/with_metrics/novel_object/20211117_1326_R/20211117_1326_38-2_21-2_R_q2/0.2.csv",
                "/hps/nobackup/birney/users/ian/MIKK_F2_tracking/with_metrics/novel_object/20211117_1326_R/20211117_1326_38-2_21-2_R_q4/0.2.csv")

## True
IN_FILES = snakemake@input
OUT_FILE = snakemake@output[[1]]

# Get metadata from IN_FILES
names(IN_FILES) = purrr::map(IN_FILES, function(IN_FILE){
  elements = IN_FILE %>% 
    stringr::str_split(pattern = "/", simplify = T)
  n_el = length(elements)
  
  out = paste(elements[n_el - 3], elements[n_el - 1], sep = "_")
  
  return(out)
}) 

# Read in file and process

out = purrr::map(IN_FILES, function(IN_FILE){
  readr::read_csv(IN_FILE)
}) %>% 
  # bind into single DF
  dplyr::bind_rows(.id = "video") %>% 
  # separate metadata
  tidyr::separate(col = video,
                  into = c("assay1", "assay2", "date", "time", "pat_line", "mat_line", "tank_side", "quadrant"),
                  sep = "_") %>% 
  # unite assay columns
  tidyr::unite(col = "assay",
               assay1, assay2,
               sep = "_",
               remove = T) %>%
  # add `test_fish`
  tidyr::unite(col = "test_fish",
               pat_line, mat_line,
               sep = "x",
               remove = F) %>% 
  # add `ref_fish`
  dplyr::mutate(ref_fish = "iCab") %>%
  # adjust `pat_line` and `mat_line` for iCab reference fishes
  dplyr::mutate(pat_line = dplyr::if_else(fish == "ref",
                                          "iCab",
                                          pat_line),
                mat_line = dplyr::if_else(fish == "ref",
                                          "iCab",
                                          mat_line)) %>%
  # remove unnecessary columns
  dplyr::select(-c(x_lag1, y_lag1, x_lag2, y_lag2, distance_b)) %>% 
  # order columns
  dplyr::select(assay, date, time, ref_fish, test_fish, tank_side, quadrant, frame, seconds, fish, x, y, distance, angle, pat_line, mat_line) %>% 
  # drop NAs (they don't work with the HMM)
  tidyr::drop_na()

# Write to file

readr::write_csv(out, OUT_FILE)
