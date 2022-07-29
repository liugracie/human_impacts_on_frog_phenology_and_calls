# Code to calculate breeding phenology estimates (onset, offset and median)
# for each species and grid cell.

# Warning: Running this is very time-consuming! 
# The output has been provided in the folder "Data/phenesse_output/".

# The complete raw dataset used to obtain these phenology estimates 
# (i.e. "Data/frogid_dat.rds") has not been provided due to sensitivities 
# in relation to locations of rare or threatened species and 
# citizen scientist information. However, the data, with sensitive species'
# localities removed or buffered, are made available annually 
# (Rowley and Callaghan, 2020; data available at GBIF: 
# https://doi.org/10.15468/wazqft and 
# FrogID: https://www.frogid.net.au/explore).

library(tidyverse)
library(phenesse)

# The grid size to be used, in decimal degrees: 0.05, 0.1, or 0.25
grid_size_deg <- 0.1 

## Clean and load the data

get_clean_dat <- function(grid_size_deg) {
  
  dat_initial <- 
    read_rds("Data/frogid_dat.rds") %>% 
    select(id:capture_time_in_zone, user_id, 
           species, date, 
           paste0("grid_id_",grid_size_deg,"deg"), 
           paste0("GHM_",grid_size_deg,"deg")) %>% 
    rename(grid_id = paste0("grid_id_",grid_size_deg,"deg"), 
           GHM = paste0("GHM_",grid_size_deg,"deg")) %>% 
    mutate(nth_day = as.numeric(format(date, "%j"))) %>%
    distinct() %>% 
    # remove duplicate records, i.e. obs with the 
    # same species-grid-date-location (lat, lng) combo
    mutate(latlng = paste(lat, lng, sep=",")) %>% 
    group_by(species, grid_id, date, latlng) %>% 
    slice(1L) %>% # select the first
    ungroup() %>% 
    select(-latlng)
  
  # species-grid combo with enough data
  spp_grid_to_keep <- 
    dat_initial %>% 
    group_by(species, grid_id, nth_day) %>% 
    tally() %>% # n obs per day per grid for each species
    group_by(species, grid_id) %>% 
    summarise(n_days = n(), 
              n_obs = sum(n, na.rm = T),
              avg_obs_per_day = mean(n, na.rm = T), 
              sd_obs_days = sd(n, na.rm = T)) %>% 
    ungroup() %>% 
    # include only grids with 10 or more obs over 3 or more days
    filter(n_obs >= 10, n_days >= 3) %>% 
    group_by(species) %>% 
    mutate(n_grids = n()) %>% # n grids with that fulfil the min obs criteria
    ungroup() %>% 
    # include only spp with enough obs in 5 or more grids
    filter(n_grids >= 5)
  
  # only keep species-grid-obs that fulfil min obs and grid criteria above
  dat_clean <- 
    spp_grid_to_keep %>% 
    left_join(dat_initial, by = c("species", "grid_id"))
  
  # df containing non-breeding season for each species
  offseason_df <- 
    read_csv("Data/season_start.csv") %>% 
    mutate(
      # define the start of the season (julian day)
      season_start_day = as.Date(paste(2020, base_month, 1, sep = "-")) %>% 
        format(., "%j") %>% 
        as.numeric()) %>% 
    select(species, season_start_day) %>% 
    distinct() 
  
  # adjust obs nth day according to the season start (defined above)
  dat_clean <- 
    dat_clean %>% 
    left_join(offseason_df, by = "species") %>% 
    mutate(
      nth_day_adj = ifelse(nth_day < season_start_day, 
                           nth_day + (366 - season_start_day + 1), 
                           nth_day - (season_start_day - 1))
    )
  
  return(dat_clean)
  
}

dat_clean <- get_clean_dat(grid_size_deg)

## Get the estimates for each species and grid

# function to get estimates
get_estimate <- function(x, n_iter, perct){
  tibble(estimate = weib_percentile(observations = x$nth_day_adj, 
                                    iterations = n_iter, 
                                    percentile = perct))
}

lapply_with_error <- function(X,FUN,...){    
  lapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                      error=function(e) NULL))
}

# function to apply
do_analysis <- function(species_name){
  
  message(paste0("Analysing ", species_name))
  
  tmp <- dat_clean %>%
    dplyr::filter(species == species_name)
  
  # split by grid
  grid_function <- function(grid_id_to_analyse){
    
    message(paste0("Analysing ", grid_id_to_analyse))
    
    tmp2 <- tmp %>%
      dplyr::filter(grid_id == grid_id_to_analyse)
    
    onset <- get_estimate(tmp2, n_iter = 250, perct = 0.05) %>% 
      mutate(species = species_name) %>%
      mutate(grid_id = grid_id_to_analyse) %>% 
      rename_at(., vars(grid_id), ~paste0("grid_id_",grid_size_deg,"deg"))
    
    saveRDS(onset, paste0("Data/phenesse_output/grids_",grid_size_deg,"/onset/", gsub(" ", "_", species_name), "_", grid_id_to_analyse, ".RDS"))
    
    offset <- get_estimate(tmp2, n_iter = 250, perct = 0.95) %>%
      mutate(species = species_name) %>%
      mutate(grid_id = grid_id_to_analyse) %>% 
      rename_at(., vars(grid_id), ~paste0("grid_id_",grid_size_deg,"deg"))
    
    saveRDS(offset, paste0("Data/phenesse_output/grids_",grid_size_deg,"/offset/", gsub(" ", "_", species_name), "_", grid_id_to_analyse, ".RDS"))
    
    median <- get_estimate(tmp2, n_iter = 250, perct = 0.50) %>%
      mutate(species = species_name) %>%
      mutate(grid_id = grid_id_to_analyse) %>% 
      rename_at(., vars(grid_id), ~paste0("grid_id_",grid_size_deg,"deg"))
    
    saveRDS(median, paste0("Data/phenesse_output/grids_",grid_size_deg,"/median/", gsub(" ", "_", species_name), "_", grid_id_to_analyse, ".RDS"))
  }
  
  lapply_with_error(unique(tmp$grid_id), grid_function)
  
}

lapply(unique(dat_clean$species), do_analysis)
