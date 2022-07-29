# Code for main analyses and figures.

# Check if the required packages are installed.
# If not, install and load them.
{list.of.packages <- 
    c("tidyverse", "lme4", "lmerTest", "ggeffects", "mgcv", 
      "sf", "cowplot", "grid", "gridExtra")
  new.packages <- 
    list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages, require, character.only = TRUE)
  rm(list.of.packages, new.packages)}

# Calling phenology models -----

# load data (phenological estimates)
if(file.exists("Data/phenological_estimates_0.1grid.csv")){
  dat_pheno <- 
    read_csv("Data/phenological_estimates_0.1grid.csv")
} else {
  source("R/combine_phenesse_estimates.R")
  dat_pheno <- 
    load_estimates(grid_size_deg = 0.1, min_obs_in_grid = 20) %>% 
    clean_output()
}

# function to run the linear mixed effect models
run_lme <- function(response_var, grid_size_deg) {
  
  dat_pheno <- 
    dat_pheno %>% 
    rename(grid_id = paste0("grid_id_",grid_size_deg,"deg"), 
           GHM = paste0("GHM_",grid_size_deg,"deg")) 
  
  formula <- as.formula(paste0(
    response_var, "~ GHM + (1 | grid_id) + (1 | species) + (0 + GHM | species)"
  ))
  
  lmer(formula, 
       data = dat_pheno, REML = T,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  
}

# function to organise the results in tables
get_results <- function(grid_size) {
  
  mod_dur <- run_lme("duration", grid_size)
  
  mod_onset <- run_lme("onset_adj", grid_size)
  
  mod_offset <- run_lme("offset_adj", grid_size)
  
  mod_median <- run_lme("median_adj", grid_size)
  
  models = list(mod_dur = mod_dur, mod_onset = mod_onset, 
                mod_offset = mod_offset, mod_median = mod_median)
  
  list2env(models, envir = .GlobalEnv)
  
  # Join results in a table
  get_supp_table <- function(mylme){
    
    resp_var <- 
      deparse(substitute(mylme)) %>% 
      stringr::str_split("_") %>% .[[1]] %>% .[2]
    
    fixed_df <- 
      summary(mylme)[["coefficients"]] %>% 
      as.data.frame() %>% 
      rownames_to_column("Term") %>% 
      mutate(Effect = "Fixed", Response = resp_var) %>% 
      rename(SE = `Std. Error`, t.value = `t value`, p.value = `Pr(>|t|)`)
    
    ran_df <-
      summary(mylme)[["varcor"]] %>% 
      as.data.frame() %>% 
      select(grp, var1, sdcor) %>% 
      rename(Term = grp, Estimate = sdcor) %>% 
      mutate(Term = case_when(
        var1 == "(Intercept)" ~ paste0("(1 | ", Term, ")"),
        var1 == paste0("GHM") ~ paste0("(", var1, " | ", Term, ")"),
        is.na(var1) ~ Term
      )) %>% 
      select(-var1) %>% 
      mutate(Term = paste("SD of", Term),
             Effect = "Random",
             Response = resp_var)
    
    bind_rows(fixed_df, ran_df) %>% select(8,7,1:6) 
  }
  
  model_results <- 
    get_supp_table(mod_dur) %>% 
    bind_rows(get_supp_table(mod_onset)) %>% 
    bind_rows(get_supp_table(mod_median)) %>% 
    bind_rows(get_supp_table(mod_offset)) %>% 
    mutate_at(vars(Estimate, SE, p.value), function(x) round(x, digits = 3)) %>% 
    mutate_at(vars(df,t.value), function(x) round(x, digits = 2)) 
  
  return(model_results)
}

# summarise model results
Table_2_expanded <- get_results(0.1) #full results

Table_2 <- 
  Table_2_expanded %>% filter(Effect == "Fixed") %>% select(-Effect)

## Figure 1 -----

map_grids_ghm_main_insets <- function(grid_size_deg){
  
  # grid cells to plot
  grids_to_plot <- read_csv("Data/sampled_grid_cells_0.1.csv")
  
  # define grid cell corners and convert to a shapefile
  
  # function to make a polygon 
  # (define the lat/lon of the grid cell corners) 
  # based on the grid cell centroid
  define_polygon <- function(lng, lat, deg) {
    
    st_sfc(
      st_polygon(
        list(
          cbind(c(lng-deg/2,lng+deg/2,lng+deg/2,lng-deg/2,lng-deg/2),
                c(lat-deg/2,lat-deg/2,lat+deg/2,lat+deg/2,lat-deg/2)))
      )
    )
    
  } 
  
  # apply the above function to all the centroids of the grids
  grid_geometry <- vector("list", nrow(grids_to_plot))
  
  for (i in 1:nrow(grids_to_plot)) {
    
    grid_geometry[[i]] <- 
      st_sf(geometry = define_polygon(grids_to_plot$lng[i], 
                                      grids_to_plot$lat[i], 
                                      grid_size_deg))
    
  }
  
  grid_geometry <- bind_rows(grid_geometry)
  
  # now you have a column with all the polygon/grids, 
  # add this onto the grid cell data
  grids_to_plot <- bind_cols(grids_to_plot, grid_geometry) %>% st_sf()
  
  # read in the shapefile
  # for the Australian coast
  aus_coastline <- st_read("Data/Aus_shp_file/aust_cd66states.shp")
  
  # box to draw around areas 
  box_bris_syd_melb <- 
    st_sfc(st_polygon(list(cbind(c(143,143,154,154,143),
                                 c(-38.9,-24.5,-24.5,-38.9,-38.9)))))
  box_cairns <- 
    st_sfc(st_polygon(list(cbind(c(144,144,147.8,147.8,144),
                                 c(-19.6,-16,-16,-19.6,-19.6)))))
  box_darwin <- 
    st_sfc(st_polygon(list(cbind(c(130.1,130.1,133.3,133.3,130.1),
                                 c(-14.7,-12,-12,-14.7,-14.7)))))
  box_perth <- 
    st_sfc(st_polygon(list(cbind(c(114.7,114.7,118.3,118.3,114.7),
                                 c(-35.3,-31.2,-31.2,-35.3,-35.3)))))
  box_adel <- 
    st_sfc(st_polygon(list(cbind(c(137.1,137.1,139.9,139.9,137.1),
                                 c(-36.2,-33.5,-33.5,-36.2,-36.2)))))
  box_hob <- 
    st_sfc(st_polygon(list(cbind(c(145.5,145.5,148,148,145.5),
                                 c(-43.2,-40.6,-40.6,-43.2,-43.2)))))
  
  # set crs
  crs <- "+proj=longlat +ellps=GRS80 +no_defs"
  st_crs(aus_coastline) <- crs 
  st_crs(grids_to_plot) <- crs
  st_crs(box_bris_syd_melb) <- crs
  st_crs(box_cairns) <- crs
  st_crs(box_darwin) <- crs
  st_crs(box_perth) <- crs
  st_crs(box_adel) <- crs
  st_crs(box_hob) <- crs
  
  # now plot it
  overall_map <- 
    ggplot()+
    geom_sf(data = aus_coastline, color="black", fill="grey90") + 
    geom_sf(data = grids_to_plot, 
            aes(fill = GHM, colour = GHM)) + 
    geom_sf(data = box_bris_syd_melb, colour = "red", fill=NA) +
    geom_sf(data = box_cairns, colour = "red", fill=NA) +
    geom_sf(data = box_darwin, colour = "red", fill=NA) +
    geom_sf(data = box_perth, colour = "red", fill=NA) +
    geom_sf(data = box_adel, colour = "red", fill=NA) +
    geom_sf(data = box_hob, colour = "red", fill=NA) +
    scale_fill_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_color_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_x_continuous(limits = c(112.5, 154.2), expand = c(0,0)) + 
    scale_y_continuous(limits = c(-44.5,-8.5), expand = c(0,0)) +
    theme_void() +
    theme(legend.position = c(0.27,0.08), legend.direction = "horizontal") +
    annotate("text", y = -12.3, x = 129, label = "a", size=5) +
    annotate("text", y = -32, x = 113.54, label = "b", size=5) +
    annotate("text", y = -36.55, x = 136, label = "c", size=5) +
    annotate("text", y = -42, x = 149.52, label = "d", size=5) +
    annotate("text", y = -25.5, x = 144.2, label = "e", size=5) +
    annotate("text", y = -17, x = 149, label = "f", size=5) 
  
  # zoomed maps
  zoomed_eastcoast <-
    ggplot() +
    geom_sf(data = aus_coastline, color="black", fill="grey90") + 
    geom_sf(data = grids_to_plot, 
            aes(fill = GHM, colour = GHM)) + 
    geom_sf(data = box_bris_syd_melb, colour = "black", fill=NA) +
    scale_fill_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_color_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_x_continuous(limits = c(143,154), expand = c(0,0)) +
    scale_y_continuous(limits = c(-38.9,-24.49), expand = c(0,0)) + 
    theme_void() +
    theme(legend.position = "none") + 
    annotate("text", y = -25.1, x = 143.8, label = "(e)", size=5)
  
  zoomed_cairns <-
    ggplot()+
    geom_sf(data = aus_coastline, color="black", fill="grey90") + 
    geom_sf(data = grids_to_plot, 
            aes(fill = GHM, colour = GHM)) + 
    geom_sf(data = box_cairns, colour = "black", fill=NA) +
    scale_fill_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_color_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_x_continuous(limits = c(144,147.8), expand = c(0,0)) +
    scale_y_continuous(limits = c(-19.6,-16), expand = c(0,0)) + 
    theme_void() +
    theme(legend.position = "none") + 
    annotate("text", y = -16.2, x = 144.22, label = "(f)", size=5)
  
  zoomed_darwin <-
    ggplot()+
    geom_sf(data = aus_coastline, color="black", fill="grey90") + 
    geom_sf(data = grids_to_plot, 
            aes(fill = GHM, colour = GHM)) + 
    geom_sf(data = box_darwin, colour = "black", fill=NA) +
    scale_fill_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_color_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_x_continuous(limits = c(130.1,133.3), expand = c(0,0)) +
    scale_y_continuous(limits = c(-14.7,-12), expand = c(0,0)) + 
    theme_void() +
    theme(legend.position = "none") + 
    annotate("text", y = -12.15, x = 130.3, label = "(a)", size=5)
  
  zoomed_perth <-
    ggplot()+
    geom_sf(data = aus_coastline, color="black", fill="grey90") + 
    geom_sf(data = grids_to_plot, 
            aes(fill = GHM, colour = GHM)) + 
    geom_sf(data = box_perth, colour = "black", fill=NA) +
    scale_fill_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_color_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_x_continuous(limits = c(114.7,118.3), expand = c(0,0)) +
    scale_y_continuous(limits = c(-35.3,-31.2), expand = c(0,0)) + 
    theme_void() +
    theme(legend.position = "none") + 
    annotate("text", y = -31.43, x = 114.95, label = "(b)", size=5)
  
  zoomed_adel <-
    ggplot()+
    geom_sf(data = aus_coastline, color="black", fill="grey90") + 
    geom_sf(data = grids_to_plot, 
            aes(fill = GHM, colour = GHM)) + 
    geom_sf(data = box_adel, colour = "black", fill=NA) +
    scale_fill_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_color_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_x_continuous(limits = c(137.1,139.9), expand = c(0,0)) +
    scale_y_continuous(limits = c(-36.2,-33.5), expand = c(0,0)) + 
    theme_void() +
    theme(legend.position = "none") + 
    annotate("text", y = -33.65, x = 137.33, label = "(c)", size=5)
  
  zoomed_hobart <-
    ggplot()+
    geom_sf(data = aus_coastline, color="black", fill="grey90") + 
    geom_sf(data = grids_to_plot, 
            aes(fill = GHM, colour = GHM)) + 
    geom_sf(data = box_hob, colour = "black", fill=NA) +
    scale_fill_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_color_viridis_c(name = "GHM", limits = seq(0, 1)) + 
    scale_x_continuous(limits = c(145.5,148), expand = c(0,0)) +
    scale_y_continuous(limits = c(-43.2,-40.6), expand = c(0,0)) + 
    theme_void() +
    theme(legend.position = "none") + 
    annotate("text", y = -40.73, x = 145.67, label="(d)", size=5)
  
  # combine panels
  plot_grid(
    plot_grid(zoomed_darwin + 
                theme(plot.margin = margin(0,0,8,0,"pt")), 
              zoomed_perth + 
                theme(plot.margin = margin(8,0,0,0,"pt")), 
              align = "h", ncol = 1), 
    plot_grid(overall_map + 
                theme(plot.margin = margin(0,0,2,0,"pt")), 
              plot_grid(zoomed_adel + 
                          theme(plot.margin = margin(5,5,5,5,"pt")), 
                        zoomed_hobart + 
                          theme(plot.margin = margin(5,5,5,5,"pt")), 
                        align = "h", ncol = 2), 
              ncol = 1, rel_heights = c(1.5,1)),
    plot_grid(zoomed_cairns + 
                theme(plot.margin = margin(0,0,8,0,"pt")), 
              zoomed_eastcoast + 
                theme(plot.margin = margin(8,0,0,0,"pt")), 
              ncol = 1, align="h"),
    nrow = 1, 
    ncol = 3, 
    rel_widths = c(1,1.5,1), 
    align = "h", 
    axis = "b")
}

map_grids_ghm_main_insets(0.1)

## Figure 2 -----

ggeffects_plot <- function(response_var, grid_size_deg, plotcolour = "blue") {
  
  dat_pheno <- 
    dat_pheno %>% 
    rename(grid_id = paste0("grid_id_",grid_size_deg,"deg"), 
           GHM = paste0("GHM_",grid_size_deg,"deg")) 
  
  formula <- as.formula(paste0(
    response_var, "~ GHM + (1 | grid_id) + (1 | species) + (0 + GHM | species)"
  ))
  
  fit <- 
    lmer(formula, 
         data = dat_pheno, REML = T,
         control = lmerControl(optimizer = "bobyqa", 
                               optCtrl = list(maxfun = 2e5)))
  
  # plot
  df <- ggpredict(fit, terms = "GHM")
  
  ggplot(df, aes(x, predicted)) +
    geom_line(colour = plotcolour) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
                alpha = .2, 
                fill = plotcolour) + 
    labs(y = paste(str_to_title(sub("_.*", "", response_var))), 
         x ="GHM") +
    theme_classic() + 
    theme(panel.background = element_rect(fill = NULL, colour = "black"))
  
}

dur <- ggeffects_plot("duration", 0.1) + labs(y = "Duration (days)")
on <- ggeffects_plot("onset_adj", 0.1) + labs(y = "Onset (DOY)")
med <- ggeffects_plot("median_adj", 0.1) + labs(y = "Median (DOY)")
off <- ggeffects_plot("offset_adj", 0.1) + labs(y = "Offset (DOY)")

x.grob <- grid::textGrob("         Landscape-level global human modification index (anthropogenic modification gradient)", gp = grid::gpar(fontsize = 12))

Figure_2 <- 
  grid.arrange(arrangeGrob(
    plot_grid(dur + theme(axis.title.x = element_blank()), 
              on + theme(axis.title.x = element_blank()), 
              med + theme(axis.title.x = element_blank()), 
              off + theme(axis.title.x = element_blank()), 
              labels = c("(a)", "(b)", "(c)", "(d)"), 
              label_size = 12, hjust=-0.1), 
    bottom = x.grob))

## Figure 3 -----

species_slopes <- function(mylme){
  
  resp_var <- deparse(substitute(mylme)) %>% stringr::str_split("_") %>% .[[1]] %>% .[2]
  
  mylme %>% coef() %>% .$species %>% 
    rownames_to_column("species") %>% 
    rename(intercept = "(Intercept)", slope = GHM) %>% 
    mutate(response_var = resp_var)
  
}

species_results <-
  species_slopes(mod_dur) %>% 
  bind_rows(species_slopes(mod_onset)) %>% 
  bind_rows(species_slopes(mod_offset)) %>% 
  bind_rows(species_slopes(mod_median))

my_spp_order <- # specify the order to plot species based on duration trend
  species_results %>% 
  filter(response_var == "dur") %>% 
  arrange(slope) %>% 
  .$species

f_labels <- # labels to add to the panels
  data.frame(panel = c("dur", "onset_offset"), 
             label = c("(a)", "(b)"))

panel.labels <- # new panel labels
  c("Duration of breeding season", 
    "Breeding season onset and offset") 

names(panel.labels) <- c("dur", "onset_offset")

annot <- # annotations to add to graph
  data.frame(panel = c("dur", 
                       "onset_offset", 
                       "onset_offset"), 
             label = c("Longer breeding season", 
                       "Earlier \nonset / offset", 
                       "Later \nonset / offset"))

species_results %>% 
  filter(response_var != "median") %>% # leave median values out of the figure
  mutate(panel = ifelse(response_var == "dur", "dur", "onset_offset")) %>%
  mutate(species = factor(species, my_spp_order)) %>% 
  ggplot(aes(y = slope, 
             x = species, 
             colour = response_var)) + 
  geom_hline(yintercept = 0, 
             linetype = "dotted") +  
  geom_point() + 
  scale_color_manual(name = NULL, 
                     values = c("black", 
                                "orangered2", 
                                "blue")) +
  theme(axis.text.y = element_text(angle = 0, 
                                   hjust = 1, 
                                   vjust = 0.5, 
                                   face = "italic"),
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major.y = element_line(colour="grey90"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", 
                                    fill = "NA"), 
        legend.position = c(0.9,0.9), 
        legend.key = element_rect(fill = "transparent", 
                                  colour = "transparent"), 
        legend.background = element_rect(fill = "transparent"), 
        strip.text = element_text(colour = "white", 
                                  size = 10),
        strip.background = element_rect(fill = "grey15", 
                                        colour = "black")) + 
  labs(y = "Difference in days (phenological trend) along \nglobal human modification index (anthropogenic modification gradient)") + 
  coord_flip() + 
  facet_wrap(~panel, 
             scales = "free_x", 
             labeller = labeller(panel = panel.labels)) +
  # add a and b labels to panel (x and y axes are reversed)
  geom_text(y = c(1,-59), x = 41.5, 
            aes(label = label, fontface = "bold"), data = f_labels, 
            inherit.aes = F) +
  # add increasing/earlier/delayed annotations to panels 
  geom_text(y = c(18, -18, 28.5), x = 41.2, 
            aes(label = label), data = annot, size = 3.3,
            inherit.aes = F) +
  # add arrows to panels 
  geom_segment(aes(y = c(11, -8, 18.5), yend = c(25, -28, 38.5), 
                   x = c(40.1, 39.5, 39.5), xend = c(40.1, 39.5, 39.5)), 
               arrow = arrow(length = unit(0.02, "npc"), type = "closed"), 
               data = annot, 
               inherit.aes = F) +
  theme(legend.position = "none") # turn off legend


# Call characteristics models -----

# load data
dat_calls <- read_csv("Data/call_characteristics.csv")

Lit.per <- dat_calls %>% filter(species=="Litoria peronii")
Lim.per <- dat_calls %>% filter(species=="Limnodynastes peronii")
L.caer <- dat_calls %>% filter(species=="Litoria caerulea")

## Table 3 -----
## summary of call parameters

Table_3 <-
  dat_calls %>% 
  select(species, mean_call_dur.s) %>%
  filter(complete.cases(.)) %>% 
  group_by(species) %>% 
  summarise(
    Call_parameter = "Call duration (s)",
    Mean = mean(mean_call_dur.s, na.rm = T),
    SD = sd(mean_call_dur.s, na.rm = T),
    Minimum = min(mean_call_dur.s, na.rm = T),
    Maximum = max(mean_call_dur.s, na.rm = T),
    N = n()
  ) %>% 
  bind_rows(
    dat_calls %>% 
      select(species, call_RR.calls_per_min) %>%
      filter(complete.cases(.)) %>% 
      group_by(species) %>% 
      summarise(
        Call_parameter = "Call repetition rate (calls/min)",
        Mean = mean(call_RR.calls_per_min, na.rm = T),
        SD = sd(call_RR.calls_per_min, na.rm = T),
        Minimum = min(call_RR.calls_per_min, na.rm = T),
        Maximum = max(call_RR.calls_per_min, na.rm = T),
        N = n())
  ) %>% 
  bind_rows(
    dat_calls %>% 
      select(species, silent_interval.s) %>%
      filter(complete.cases(.)) %>% 
      group_by(species) %>% 
      summarise(Call_parameter = "Intercall interval (s)",
                Mean = mean(silent_interval.s, na.rm = T),
                SD = sd(silent_interval.s, na.rm = T),
                Minimum = min(silent_interval.s, na.rm = T),
                Maximum = max(silent_interval.s, na.rm = T),
                N = n())
  ) %>% 
  bind_rows(
    dat_calls %>% 
      select(species, mean_calling_effort.call_dur_div_call_period) %>%
      filter(complete.cases(.)) %>% 
      group_by(species) %>% 
      summarise(
        Call_parameter = "Calling effort (call duration/call period)",
        Mean = mean(mean_calling_effort.call_dur_div_call_period, na.rm = T),
        SD = sd(mean_calling_effort.call_dur_div_call_period, na.rm = T),
        Minimum = min(mean_calling_effort.call_dur_div_call_period, na.rm = T),
        Maximum = max(mean_calling_effort.call_dur_div_call_period, na.rm = T),
        N = n())
  ) %>% 
  bind_rows(
    dat_calls %>% 
      select(species, mean_dfreq.Hz) %>%
      filter(complete.cases(.)) %>% 
      group_by(species) %>% 
      summarise(Call_parameter = "Dominant frequency (Hz)",
                Mean = mean(mean_dfreq.Hz, na.rm = T),
                SD = sd(mean_dfreq.Hz, na.rm = T),
                Minimum = min(mean_dfreq.Hz, na.rm = T),
                Maximum = max(mean_dfreq.Hz, na.rm = T),
                N = n())
  ) %>% 
  bind_rows(
    dat_calls %>% 
      select(species, mean_notes_per_call) %>%
      filter(species == "Litoria peronii", complete.cases(.)) %>% 
      group_by(species) %>% 
      summarise(Call_parameter = "Notes per call",
                Mean = mean(mean_notes_per_call, na.rm = T),
                SD = sd(mean_notes_per_call, na.rm = T),
                Minimum = min(mean_notes_per_call, na.rm = T),
                Maximum = max(mean_notes_per_call, na.rm = T),
                N = n())
  ) %>% 
  bind_rows(
    dat_calls %>% 
      select(species, mean_note_RR.notes_per_s) %>%
      filter(species == "Litoria peronii", complete.cases(.)) %>% 
      group_by(species) %>% 
      summarise(Call_parameter = "Note repetition rate (notes/s)",
                Mean = mean(mean_note_RR.notes_per_s, na.rm = T),
                SD = sd(mean_note_RR.notes_per_s, na.rm = T),
                Minimum = min(mean_note_RR.notes_per_s, na.rm = T),
                Maximum = max(mean_note_RR.notes_per_s, na.rm = T),
                N = n())
  ) %>% 
  mutate_if(is.numeric, function(x) round(x, digits = 3))

## GAMs -----

# list of response variables
# for Limnodynastes peronii and Litoria caerulea 
# (all call characteristics that do not relate to notes)
response_vars_LimP_LC <- c("mean_call_dur.s",
                           "call_RR.calls_per_min",
                           "log(silent_interval.s)",
                           "log(mean_calling_effort.call_dur_div_call_period)",
                           "mean_dfreq.Hz")

# for Litoria peronii (all call characteristics)
response_vars_LitP <- c("mean_call_dur.s",
                        "log(call_RR.calls_per_min)",
                        "log(silent_interval.s)",
                        "log(mean_calling_effort.call_dur_div_call_period)",
                        "mean_dfreq.Hz",
                        "mean_notes_per_call",
                        "mean_note_RR.notes_per_s")

# function to run GAM with any response variable
run_gam <- function(response_var, dat){
  
  formula <- as.formula(paste(
    response_var, "~ GHM + s(nth.day, bs='cc') + s(lat,lng) + s(hourly_temp)"
  ))
  
  gam(formula, data = dat)
  
}

# function to get p-values from GAMs 
collate_p_values <- function(response_var, dat){
  
  formula <- as.formula(paste(
    response_var, "~ GHM + s(nth.day, bs='cc') + s(lat,lng) + s(hourly_temp)"
  ))
  
  gam_to_eval <- gam(formula, data = dat)
  
  p_linear <- 
    summary(gam_to_eval)$p.table %>% 
    as.data.frame() %>% 
    select(`Pr(>|t|)`) %>% rename(p.value = `Pr(>|t|)`) %>% 
    rownames_to_column(var = "term") %>% 
    mutate(p.value = round(p.value, digits=3))
  
  p_smooth <- 
    summary(gam_to_eval)$s.table %>% 
    as.data.frame() %>% 
    select("p-value") %>% rename(p.value = "p-value") %>% 
    rownames_to_column(var = "term") %>% 
    mutate(p.value = round(p.value, digits=3))
  
  results <- 
    bind_rows(p_linear, p_smooth) %>% 
    mutate(call_char = response_var, 
           species = dat$species %>% unique(),
           significance = ifelse(p.value < 0.05, 
                                 "significant", 
                                 "not-significant"),
           N = nobs(gam_to_eval))
  
  return(results)
}

# df of p-values from all gams
gam_pvalues <- 
  # p-values from Limnodynastes peronii gams
  response_vars_LimP_LC %>% 
  map_df(collate_p_values, dat = Lim.per) %>% 
  # add p-values from Litoria caerulea gams
  bind_rows(., response_vars_LimP_LC %>% 
              map_df(collate_p_values, dat = L.caer)) %>% 
  # add p-values from Litoria peronii gams
  bind_rows(., response_vars_LitP %>% 
              map_df(collate_p_values, dat = Lit.per)) %>% 
  select(call_char, species, term, p.value, N) %>% 
  pivot_wider(names_from = term, values_from = p.value) %>% 
  select(-"(Intercept)")

# arrange the table by specified order
ref_order <- 
  c("mean_call_dur.s",  
    "call_RR.calls_per_min",
    "log(call_RR.calls_per_min)",
    "log(silent_interval.s)",
    "log(mean_calling_effort.call_dur_div_call_period)",
    "mean_dfreq.Hz", 
    "mean_notes_per_call",
    "mean_note_RR.notes_per_s")

Table_4 <- 
  gam_pvalues[order(factor(gam_pvalues$call_char, levels = ref_order)), ] %>% 
  # rename call variables
  mutate(call_char = gsub("mean_call_dur.s", 
                          "Call duration", call_char),
         call_char = gsub("call_RR.calls_per_min", 
                          "Call repetition rate", call_char),
         call_char = gsub("silent_interval.s", 
                          "Intercall interval", call_char),
         call_char = gsub("mean_calling_effort.call_dur_div_call_period", 
                          "Calling effort", call_char),
         call_char = gsub("mean_dfreq.Hz", 
                          "Dominant frequency", call_char),
         call_char = gsub("mean_notes_per_call", 
                          "Notes per call", call_char),
         call_char = gsub("mean_note_RR.notes_per_s", 
                          "Note repetition rate", call_char)) 

# Get full result table of all gams
get_full_result_table <- function(response_var, dat){
  
  mygam <- run_gam(response_var, dat) %>% summary()
  
  p_linear <- 
    mygam %>% 
    .$p.table %>% 
    as.data.frame() %>% 
    rename(SE = `Std. Error`,
           p.value = `Pr(>|t|)`) %>% 
    rownames_to_column(var = "term") %>% 
    mutate_if(is.numeric, function(x) round(x, digits = 3)) %>% 
    mutate(call_char = response_var) %>% 
    pivot_wider(names_from = term, 
                values_from = `Estimate`:p.value, 
                names_vary = "slowest")
  
  p_smooth <- 
    mygam %>% 
    .$s.table %>% 
    as.data.frame() %>% 
    rename(p.value = `p-value`) %>% 
    rownames_to_column(var = "term") %>% 
    mutate_if(is.numeric, function(x) round(x, digits = 3)) %>% 
    select(-Ref.df) %>% 
    mutate(call_char = response_var) %>% 
    pivot_wider(names_from = term, 
                values_from = edf:p.value, 
                names_vary = "slowest")
  
  results <- 
    left_join(p_linear, p_smooth, by = "call_char") %>% 
    mutate(call_char = response_var, 
           species = dat$species %>% unique(),
           N = nobs(run_gam(response_var, dat)))
  
  return(results)
}

Table_S8 <- 
  # p-values from Limnodynastes peronii gams
  response_vars_LimP_LC %>% 
  map_df(get_full_result_table, dat = Lim.per) %>% 
  # add p-values from Litoria caerulea gams
  bind_rows(., response_vars_LimP_LC %>% 
              map_df(get_full_result_table, dat = L.caer)) %>% 
  # add p-values from Litoria peronii gams
  bind_rows(., response_vars_LitP %>% 
              map_df(get_full_result_table, dat = Lit.per)) %>% 
  select(call_char, species, `Estimate_(Intercept)`:N) %>% 
  .[order(factor(.$call_char, levels = ref_order)), ] %>% 
  mutate(call_char = gsub("mean_call_dur.s", 
                          "Call duration", call_char),
         call_char = gsub("call_RR.calls_per_min", 
                          "Call repetition rate", call_char),
         call_char = gsub("silent_interval.s", 
                          "Intercall interval", call_char),
         call_char = gsub("mean_calling_effort.call_dur_div_call_period", 
                          "Calling effort", call_char),
         call_char = gsub("mean_dfreq.Hz", 
                          "Dominant frequency", call_char),
         call_char = gsub("mean_notes_per_call", 
                          "Notes per call", call_char),
         call_char = gsub("mean_note_RR.notes_per_s", 
                          "Note repetition rate", call_char)) 
