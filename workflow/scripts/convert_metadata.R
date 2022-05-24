# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#IN = here::here("config/20220509_New_Fish_Movies.xlsx")

## True
IN = snakemake@input[[1]]
OUT = snakemake@output[[1]]

# Read in file and process

df = readxl::read_xlsx(IN)

# Process

out = df %>% 
  # remove empty columns
  dplyr::select(-c(`...13`, `Remarks`, `...15`)) %>% 
  # fill empty values
  tidyr::fill(`Tank Side`,
              `Icab Age (reference fish)`,
              `Icabs from Tank`,
              `Date of filming`,
              `Time`,
              `movie name`,
              `parental Generation`,
              .direction = "down") %>% 
  # recode `quadrant` from `Tank`
  dplyr::mutate(quadrant = dplyr::recode(`Tank`,
                                         "A (OL)" = "q2",
                                         "B (OR)" = "q1",
                                         "C (UR)" = "q4",
                                         "D (UL)" = "q3")) %>% 
  # change tank side first rows
  dplyr::mutate(`Tank Side` = dplyr::recode(`Tank Side`, "R =right" = "R")) %>% 
  # split `movie name` into `date`, `time`, and `tank_side`
  tidyr::separate(`movie name`,
                  into = c("date", "time", "tank_side"),
                  sep = "_",
                  remove = F) %>% 
  # remove `.avi` from `tank_side`
  dplyr::mutate(tank_side = stringr::str_remove(tank_side, ".avi")) %>% 
  # split `parental Generation`
  tidyr::separate(`parental Generation`,
                  into = c("f0_1", "f0_2"),
                  sep = "x",
                  remove = F) %>% 
  # get `mat` and `pat` lines
  dplyr::mutate(pat = dplyr::case_when(stringr::str_detect(f0_1, " male") ~ f0_1,
                                       stringr::str_detect(f0_2, " male") ~ f0_2,
                                       TRUE ~ NA_character_),
                mat = dplyr::case_when(stringr::str_detect(f0_1, "female") ~ f0_1,
                                       stringr::str_detect(f0_2, "female") ~ f0_2,
                                       TRUE ~ NA_character_)) %>% 
  # remove leading and trailing white spaces
  dplyr::mutate(dplyr::across(c("pat", "mat"),
                              ~stringr::str_remove(.x, "^ ")),
                dplyr::across(c("pat", "mat"),
                              ~stringr::str_remove(.x, " $"))) %>% 
  # extract pat and mat line and generation
  tidyr::separate(pat,
                  into = c("pat_line", NA, "pat_gen"),
                  sep = " ",
                  remove = F) %>% 
  tidyr::separate(mat,
                  into = c("mat_line", NA, "mat_gen"),
                  sep = " ",
                  remove = F) %>% 
  # add sample column
  tidyr::unite(col = "sample",
               date, time, pat_line, mat_line, tank_side, quadrant, 
               sep = "_",
               remove = F)



# checks

## `tank_side`
all(out$`Tank Side` == out$tank_side)

## `pat` and `mat`
unique(out$pat)
unique(out$mat)

# rename and reorder columns

final = out %>% 
  dplyr::select(movie_name = `movie name`,
                date,
                time,
                quadrant,
                tank_side,
                sample,
                test_line = `Line`,
                alter = Alter,
                ref_dob = `Icab Age (reference fish)`,
                pat_line,
                pat_gen,
                mat_line,
                mat_gen) %>% 
  # order by date, time, and quadrant
  dplyr::arrange(date, time, quadrant)

# save to file

readr::write_csv(final, OUT)
