# Functions to combine breeding phenology estimates 
# (onset, offset, median, duration) for each species and grid cell
# into a dataframe, then filter and clean the output.

# The 'load_estimates' function requires the raw dataset 
# (i.e. "Data/frogid_dat.rds"), which has not been provided due to 
# sensitivities in relation to locations of rare or threatened species 
# and citizen scientist information. However, the data, with sensitive 
# species' localities removed or buffered, are made available annually 
# (Rowley and Callaghan, 2020; data available at GBIF: 
# https://doi.org/10.15468/wazqft and 
# FrogID: https://www.frogid.net.au/explore).
# The  output (file containing filtered and clean phenology estimates) 
# has been provided: "Data/phenological_estimates_0.1grid.csv".

library(tidyverse)

# Load and combine onset, offset and median estimates for each species and grid cell. 
#' @param grid_size_deg The grid cell size to be used, in decimal degrees: 0.05, 0.1, or 0.25
#' @param min_obs_in_grid The cutoff of the number of observations in a grid cell to be used. 
#' Any grid cells that have less observations than this cutoff will be excluded.
#' @param n_grids_cutoff The cutoff of the number of grid cells with estimates to be used. 
#' Any species that have less grids than this cutoff will be excluded.
#' @param filter_estimates Filter out phenoestimates that are not between days 0 and 366?
#'
load_estimates <- function(grid_size_deg = 0.1, 
                           min_obs_in_grid = 20, 
                           n_grids_cutoff = 10,
                           filter_estimates = FALSE){
  
  grid_cell_estimates <- 
    # load the onset estimates
    list.files(paste0("Data/phenesse_output/grids_",grid_size_deg,"/onset/"), 
               full.names = T) %>% 
    map(read_rds) %>% 
    bind_rows() %>% 
    rename(onset = estimate) %>% 
    # add the offset estimates
    left_join(
      list.files(paste0("Data/phenesse_output/grids_",grid_size_deg,"/offset/"), 
                 full.names = T) %>% 
        map(read_rds) %>% 
        bind_rows() %>% 
        rename(offset = estimate)
    ) %>% 
    # add the median estimates
    left_join(
      list.files(paste0("Data/phenesse_output/grids_",grid_size_deg,"/median/"),
                 full.names = T) %>% 
        map(read_rds) %>% 
        bind_rows() %>% 
        rename(median = estimate)
    )
  
  # whether to filter out phenoestimates that are not between days 0 and 366 
  grid_cell_estimates <- 
    if(filter_estimates){
      grid_cell_estimates %>% 
        # remove estimates that are less than 0 or more than 366
        filter(onset >= 0 & onset <= 366) %>% 
        filter(offset >= 0 & offset <= 366) %>% 
        # calculate season duration
        mutate(duration = offset - onset) %>% 
        # reorder vars
        select(species, 
               paste0("grid_id_",grid_size_deg,"deg"), 
               median, onset, offset, duration)
    } else{
      grid_cell_estimates %>% 
        # calculate season duration
        mutate(duration = offset - onset) %>% 
        # reorder vars
        select(species, 
               paste0("grid_id_",grid_size_deg,"deg"), 
               median, onset, offset, duration)
    }
  
  # function to convert onset/offset days into actual days of the year 
  convert_to_actual_doy <- function(day, season_start_day){
    doy <- day + season_start_day
    actual_doy <- ifelse(doy > 366, 
                         doy - 366, 
                         doy)
    return(actual_doy)
  }
  
  output <- 
    grid_cell_estimates %>% 
    # add grid cell modification scores
    left_join(
      read_csv(paste0("Data/Aus_grid_files/aus_grids_",
                      format(signif(grid_size_deg, digits = 2), nsmall = 2),
                      "deg_modscores.csv")) %>% 
        select(paste0("grid_id_",grid_size_deg,"deg"), 
               paste0("GHM_",grid_size_deg,"deg")) %>% 
        distinct()
    ) %>% 
    # add the season start doy for each species
    left_join(
      read_csv("Data/season_start.csv") %>% 
        mutate(
          # define the start of the season (julian day)
          season_start_day = 
            as.Date(paste(2020, base_month, 1, sep = "-")) %>% 
            format(., "%j") %>% 
            as.numeric()
          ) %>% 
        select(species, season_start_day) %>% 
        distinct(),
      by = "species"
    )  %>% 
    # convert adjusted onset, offset and median days 
    # back to actual days of the year
    mutate(
      onset_doy = convert_to_actual_doy(onset, season_start_day),
      offset_doy = convert_to_actual_doy(offset, season_start_day),
      median_doy = convert_to_actual_doy(median, season_start_day)
    ) %>% 
    # rename original (adjusted) onset and offset days 
    rename(
      onset_adj = onset,
      offset_adj = offset,
      median_adj = median
    ) 
  
  # Now we want to remove any species-grid cell combos that 
  # do not meet the specified min_obs_in_grid and n_grids thresholds.
  
  # load a summary of the species-grid cell combos 
  # that we calculated phenoestimates for
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
    select(-latlng) %>% 
    group_by(species, grid_id, nth_day) %>% 
    tally() %>% # n obs per day per grid for each species
    group_by(species, grid_id) %>% 
    summarise(n_days = n(), 
              n_obs = sum(n, na.rm = T),
              avg_obs_per_day = mean(n, na.rm = T), 
              sd_obs_days = sd(n, na.rm = T)) %>% 
    ungroup() %>% 
    # include only grids with 10 or more obs over 3 or more days
    # note: this was our initial filtering criteria
    filter(n_obs >= 10, n_days >= 3) %>% 
    group_by(species) %>% 
    mutate(n_grids = n()) %>% # n grids that fulfil the min obs criteria
    ungroup() %>% 
    # include only spp with enough obs in 5 or more grids
    # note: this was also our initial filtering criteria
    filter(n_grids >= 5) 
  
  # update the output, i.e. filter out species-grid cell combos 
  # that do NOT meet the thresholds
  output_filtered <- 
    output %>% 
    # add a column with the number of obs in each species-grid cell combo
    left_join(dat_initial %>% 
                select(species, grid_id, n_obs) %>% 
                rename_at(vars(grid_id), ~paste0("grid_id_",grid_size_deg,"deg"))
    ) %>% 
    # filter species-grid cell combos based on min_obs_in_grid threshold
    filter(n_obs >= min_obs_in_grid) %>% 
    # count number of grids for each species
    group_by(species) %>% 
    mutate(n_grids = n()) %>% 
    ungroup() %>% 
    # include only species that have 
    # more grids with estimates than the specified cutoff
    filter(n_grids >= n_grids_cutoff)  
  
  return(output_filtered)
  
}


# Clean the output data where offset or duration > 366.
#' @param data The data file resulting from the load_estimates function above.
#' 
clean_output <- function(data){
  data <- 
    data %>% 
    mutate(
      # set any species-grid cells with duration >366 
      # to 366 to represent year-round breeding
      duration = if_else(duration >= 366, 366, duration),
      # set onset/offset/median to 0 for 
      # any species-grid cells with duration >366 
      onset_adj = if_else(duration == 366, 0, onset_adj),
      offset_adj = if_else(duration == 366, 0, offset_adj),
      median_adj = if_else(duration == 366, 0, median_adj),
      # set offset to 0 for any species-grid cells with offset >366
      offset_adj = if_else(offset_adj > 366, 0, offset_adj),
      # also set corresponding onset/offset DOYs to 0 
      onset_doy = if_else(onset_adj == 0, 0, onset_doy),
      offset_doy = if_else(offset_adj == 0, 0, offset_doy),
      median_doy = if_else(median_adj == 0, 0, median_doy)) %>% 
    # set any values that are 0 to NA
    na_if(0)
}
