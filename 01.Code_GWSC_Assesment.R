# Load Library 
library(raster)
library(tidyverse)
library(sf)
library(stringr)
library(ggplot2)
library(lubridate)
library(readxl)
library(rnaturalearth)
library(tidyterra)
library(terra)
library(fixest)
library(kableExtra)
library(patchwork)
library(viridis)
library(dplyr)
library(kableExtra)
library(corrplot)
library(car)

# Figure 1: Density of political violence and protests in Kenya, 
#including events related to environmental and agricultural issues.

# Che the density of conflit to determine is there any growth in conflicts at two years
df_conflict <- read_excel("GWSC Interview Data/Conflict Events/Africa_1997-2023_Aug04-2.xlsx")%>%
  filter(COUNTRY == "Kenya") # only consider Kenya

# Get Kenya boundary
kenya <- ne_countries(country = "Kenya", returnclass = "sf")%>%
  dplyr::select(name)

# Set grid cell size (in meters) ~ This is roughy resolution of Temp and Prcp Anamalities data -so to be consistent 
grid_size <- 25000 # 25 km x 25 km grid

# Projection
kenya_proj <- st_transform(kenya, 3857) 

# Create a grid covering Kenya with given cell size
kenya_grid <- st_make_grid(
  kenya_proj,
  cellsize = c(grid_size, grid_size),
  what = "polygons"
)

# Clip - to keep pixels within Kneya 
kenya_grid <- st_intersection(st_sf(geometry = kenya_grid), kenya_proj)

# Use WGS84 for mapping
kenya_grid <- st_transform(kenya_grid, 4326)

# Use year 1997 and 2023 (starting and end year) - the two extreme ponts of conflicts data 
df_conflict_filter <- df_conflict %>% filter(YEAR %in% c(1997, 2023))%>%
  select(EVENT_ID_CNTY,YEAR,LONGITUDE,LATITUDE)


# create sf object
df_conflict_sf <- st_as_sf(df_conflict_filter, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)


# chcek crs
df_conflict_sf <- st_transform(df_conflict_sf, st_crs(kenya_grid))

# Id for easy identification
kenya_grid <- kenya_grid %>% mutate(grid_id = row_number())

# Spatial join
joined <- st_join(kenya_grid,df_conflict_sf)%>%
  drop_na()

# conflicts per grid cell - by grid id and year
conflict_counts <- joined %>%
  st_drop_geometry() %>%
  group_by(grid_id,YEAR) %>%
  summarise(conflict_count = n(), .groups = "drop")

# counts back to the grid
kenya_grid_counts <- kenya_grid %>%
  left_join(conflict_counts, by = "grid_id") %>%
  mutate(conflict_count = replace_na(conflict_count, 0))


# define conflics bins
grid_binned <- kenya_grid_counts %>%
  mutate(conflict_bin = case_when(
    conflict_count == 0 ~ "0",
    conflict_count > 0 & conflict_count <= 5 ~ "1-5",
    conflict_count > 5 & conflict_count <= 10 ~ "6-10",
    conflict_count > 10 & conflict_count <= 20 ~ "11-20",
    conflict_count > 20 & conflict_count <= 50 ~ "21-50",
    conflict_count > 50 ~ ">50"
  ))

grid_binned$conflict_bin <- factor(
  grid_binned$conflict_bin,
  levels = c("0", "1-5", "6-10", "11-20", "21-50", ">50")
)


# data for 1997 and 2023
grid_binned_1997 <- grid_binned %>% filter(YEAR == 1997)
grid_binned_2023 <- grid_binned %>% filter(YEAR == 2023)

# legend
grid_binned_1997$conflict_bin <- factor(grid_binned_1997$conflict_bin,
                                        levels = c("0", "1-5", "6-10", "11-20", "21-50", ">50"))
grid_binned_2023$conflict_bin <- factor(grid_binned_2023$conflict_bin,
                                        levels = c("0","1-5", "6-10", "11-20", "21-50", ">50"))

# plot
plot_conflicts <- function(data, year) {
  ggplot(data) +
    geom_sf(aes(fill = conflict_bin), color = NA, show.legend = TRUE) +
    geom_sf(data = kenya_proj, fill = NA, color = "black", size = 0.5) +
    geom_sf(data = kenya_grid, fill = NA, color = "grey", size = 0.00001) +
    scale_fill_manual(
      values = c(
        "0" = "white",
        "1-5" = "#fde725",
        "6-10" = "#F97316",
        "11-20" = "#5ec962",
        "21-50" = "#21918c",
        ">50" = "#440154"
      ),
      drop = FALSE,
      name = "Number of Conflicts",
      guide = guide_legend(override.aes = list(color = "grey")) # outline legend boxes
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",  # hide legend in each subplot
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white")
    ) +
    labs(title = paste("Conflict Density (", year, ")", sep = ""))
}

# combine plot
p1 <- plot_conflicts(grid_binned_1997, 1997)
p2 <- plot_conflicts(grid_binned_2023, 2023)

conflict_density_image <- (p1 + p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")


# Save to PNG
ggsave(
  filename = "Figures/conflict_density_image.png",
  plot = conflict_density_image,  # Explicitly specify plot object
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

conflict_density_image


# Figure 2: Precipitation anomalies across Kenya (a) temporal patterns and (b) spatial distribution
# Precipitation Anomalies data
files <- list.files("GWSC Interview Data/Precip and Temp Monthly Anomalies/Prcp/",
                    pattern = "\\.tif$", full.names = TRUE)
# kenya
kenya_vect <- vect(kenya)

# process each raster file there are files from 2000 - 2020 at monthly frequent
process_raster <- function(file) {
  year <- str_extract(file, "(?<=time)\\d{4}")
  month <- str_extract(file, "(?<=time\\d{4}-)\\d{2}")
  r <- rast(file)
  r <- mask(crop(r, kenya_vect), kenya_vect)
  polys <- as.polygons(r, values = TRUE, dissolve = FALSE, na.rm = TRUE)
  polys_sf <- st_as_sf(polys)
  names(polys_sf)[1] <- "prcp_anomaly"
  polys_sf$year <- as.integer(year)
  polys_sf$month <- as.integer(month)
  polys_sf
}

# combine them into one file 
prcp_polygons_all <- lapply(files, process_raster) %>% do.call(rbind, .)

# time series graph - with annual mean across multiple collection sites at 0.25 degree resolution 
prcp_graph <- prcp_polygons_all %>%
  st_drop_geometry() %>%
  group_by(year, month) %>%
  summarize(
    ave_prcp_anomaly = mean(prcp_anomaly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    date = make_date(year, month, 1),
    month_name = month.abb[month]
  )

# Time series plot with - precipitation anomalies - montly average across Kenya 
prcp_anamoly_timeseries_image <- ggplot(prcp_graph, aes(x = date, y = ave_prcp_anomaly)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "",
    subtitle = "",
    x = "Year",
    y = "Precipitation Anomaly (mm)",
    caption = ""
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank()
  )

# Save to PNG
ggsave(
  filename = "Figures/prcp_anamoly_timeseries_image.png",
  plot = prcp_anamoly_timeseries_image, 
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

# Define years (2010-2020) - long period may hidden some recent variation so keep it for 3 year
selected_years <- 2018:2020


pixel_avg <- prcp_polygons_all %>%
  filter(year %in% selected_years) %>%
  group_by(geometry) %>%
  summarize(
    mean_prcp_anomaly = mean(prcp_anomaly, na.rm = TRUE),
    years_included = paste(unique(year), collapse = ", "),
    n_years = n_distinct(year), 
    .groups = "drop"
  ) %>%
  st_as_sf()

# clip to Kenya boundary 
pixel_avg <- st_intersection(pixel_avg, kenya)

# plot spatial variation in precepitation anamoly - from 2018-2023
prcp_anamoly_spatial_image <- ggplot(pixel_avg) +
  geom_sf(aes(fill = mean_prcp_anomaly), color = NA, alpha = 0.8) +
  geom_sf(data = kenya, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "maroon", mid = "white", high = "darkgreen",
    midpoint = 0,
    name = "Precip. Anomaly (mm)",
    limits = c(-0.2, 45),
    oob = scales::squish
  ) +
  labs(title = "") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
  )+
  coord_sf(crs = 4326) +
  theme_void()

# Save to PNG
ggsave(
  filename = "Figures/prcp_anamoly_spatial_image.png",
  plot = prcp_anamoly_spatial_image,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

prcp_anamoly_timeseries_image
prcp_anamoly_spatial_image



# Figure 3: Temperature anomalies across across Kenya (a) temporal patterns and (b) spatial distribution

# Temperature Anomalies data
files <- list.files("GWSC Interview Data/Precip and Temp Monthly Anomalies/Temp/",
                    pattern = "\\.tif$", full.names = TRUE)


process_raster <- function(file) {
  year <- str_extract(file, "(?<=time)\\d{4}")
  month <- str_extract(file, "(?<=time\\d{4}-)\\d{2}")
  r <- rast(file)
  r <- mask(crop(r, kenya_vect), kenya_vect)
  polys <- as.polygons(r, values = TRUE, dissolve = FALSE, na.rm = TRUE)
  polys_sf <- st_as_sf(polys)
  names(polys_sf)[1] <- "temp_anomaly"
  polys_sf$year <- as.integer(year)
  polys_sf$month <- as.integer(month)
  polys_sf
}

# combine into one big sf dataframe
temp_polygons_all <- lapply(files, process_raster) %>% do.call(rbind, .)

# # time series graph - with annual mean across multiple collection sites at 0.25 degree resolution 
temp_graph <- temp_polygons_all %>%
  st_drop_geometry() %>%
  group_by(year, month) %>%
  summarize(
    ave_temp_anomaly = mean(temp_anomaly, na.rm = TRUE),  # Correct column name
    .groups = "drop"
  ) %>%
  mutate(
    date = make_date(year, month, 1),
    month_name = month.abb[month]
  )

# Time series plot for temp. anomalies
temp_anamoly_timeseries_image <- ggplot(temp_graph, aes(x = date, y = ave_temp_anomaly)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "",
    subtitle = "",
    x = "Year",
    y = "Temperature Anomaly (°C)",  # Changed unit
    caption = ""
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background = element_blank()
  )

# Save to PNG
ggsave(
  filename = "Figures/temp_anamoly_timeseries_image.png",
  plot = temp_anamoly_timeseries_image,  # Explicitly specify plot object
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)


pixel_avg <- temp_polygons_all %>%
  filter(year %in% selected_years) %>%  
  group_by(geometry) %>%             
  summarize(
    mean_temp_anomaly = mean(temp_anomaly, na.rm = TRUE),
    years_included = paste(unique(year), collapse = ", "),  
    n_years = n_distinct(year),        
    .groups = "drop"
  ) %>%
  st_as_sf()

# clip to Kenya boundary
pixel_avg <- st_intersection(pixel_avg, kenya)

# plot spatial variation in temp. anamoly - from 2018-2023
temp_anamoly_spatial_image <- ggplot(pixel_avg) +
  geom_sf(aes(fill = mean_temp_anomaly), color = NA, alpha = 0.8) +
  geom_sf(data = kenya, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "maroon", mid = "white", high = "darkgreen",
    midpoint = 0,
    name = "Temp. Anomaly (mm)",
    limits = c(-0.6, 0.6),  # Adjust based on your data
    oob = scales::squish  # Handle out-of-bounds values
  ) +
  labs(title = "") +
  coord_sf(crs = 4326) +
  theme_void()

# Save to PNG
ggsave(
  filename = "Figures/temp_anamoly_spatial_image.png",
  plot = temp_anamoly_spatial_image,  # Explicitly specify plot object
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

temp_anamoly_timeseries_image
temp_anamoly_spatial_image

# Figure 4: Areas with changes in land classification between 2000 and 2015.

# land cover data for 2000 and 2015 - two years at extreme where to deterct any land conversion is going on 
lc_2000 <- rast("./GWSC Interview Data/Land Cover Classification 2000 through 2020/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7cds.nc")$lccs_class %>%
  crop(kenya) %>% mask(kenya)

lc_2018 <- rast("./GWSC Interview Data/Land Cover Classification 2000 through 2020/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7cds.nc")$lccs_class %>% 
  crop(kenya) %>% mask(kenya)

# create change rasters directly
change_raster <- lc_2000 

# calculate all changes
values(change_raster) <- case_when(
  # Grassland and Crop Land conversions
  lc_2000[] == 10 & lc_2018[] == 120 ~ 1,  
  lc_2000[] == 10 & lc_2018[] == 130 ~ 1,  
  lc_2000[] == 120 & lc_2018[] == 10 ~ 1,  
  lc_2000[] == 130 & lc_2018[] == 10 ~ 1, 
  lc_2000[] != lc_2018[] ~ 2,              
  TRUE ~ 0                                 
)

# categorical raster for plotting
change_cats <- change_raster
levels(change_cats) <- data.frame(
  id = 0:2,
  change_type = c("No Change", 
                  "Grassland and Crop Land conversion",
                  "Other Change")
)

# mask the raster by kenya
change_cats <- change_cats %>% 
  crop(kenya) %>%    
  mask(kenya) 

#  CRS (EPSG:4326) - map purpose
kenya <- st_transform(kenya, crs = 4326)
change_cats <- project(change_cats, "EPSG:4326")

# land class change plot - 2000 - 2015
land_classification <- ggplot() +
  geom_spatraster(data = change_cats) +
  geom_sf(data = kenya, fill = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(
    name = " ",
    values = c(
      "No Change" = "white",
      "Grassland and Crop Land conversion" = "green",
      "Other Change" = "maroon"
    ),
    na.value = NA,  # Makes NA values invisible
    drop = TRUE,    # Drops unused levels
    breaks = c("No Change", "Grassland and Crop Land conversion", "Other Change")  # Explicitly specify legend items
  ) +
  labs(title = "Areas with changes in land classification during 2000–2018") +
  theme_void()+
  theme(
    legend.position = "bottom",  # Move legend to bottom
    plot.title = element_text(hjust = 0.5, face = "bold")  # Center title
  )


# Save to PNG
ggsave(
  filename = "Figures/land_classification.png",
  plot = land_classification,  # Explicitly specify plot object
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

land_classification

#Figure 5: Rainfed Crop Production Area (All Crops, 2000)

# Load data Global data set of Monthly Irrigated and Rainfed Crop Areas around the year 2000 (MIRCA2000)
mirca_kenya <- raster("./GWSC Interview Data/Rainfed area/annual_area_harvested_rainfed_allcrops_ha_30mn.asc") %>%
  {
    crs(.) <- "+proj=longlat +datum=WGS84"  # Set CRS if missing
    .
  } %>%
  crop(
    ne_countries(country = "Kenya", returnclass = "sf") %>%
      st_transform(4326)  # Kenya in WGS84
  ) %>%
  mask(kenya)

# Resolution ~ 55km X 55Km
res(mirca_kenya )  # Returns c(x_resolution, y_resolution)


# column name
mirca_df <- as.data.frame(mirca_kenya, xy = TRUE, na.rm = TRUE) %>%
  rename(rainfed_area = 3)

# plot
mirca_2000 <- ggplot() +
  geom_raster(aes(x, y, fill = rainfed_area), data = mirca_df) +
  geom_sf(data = kenya, fill = NA, color = "red", linewidth = 0.7) +
  scale_fill_viridis(
    name = "Rainfed Crop Area (ha)",
    option = "plasma",
    limits = range(mirca_df$rainfed_area, na.rm = TRUE)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm")
  )


# Save to PNG
ggsave(
  filename = "Figures/mirca_2000.png",
  plot = mirca_2000,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

mirca_2000

# create analyzin dataframe
# load conflict data
df_conflict <- read_excel("GWSC Interview Data/Conflict Events/Africa_1997-2023_Aug04-2.xlsx") %>%
  filter(COUNTRY == "Kenya") %>%
  select(EVENT_ID_CNTY, YEAR,EVENT_DATE, LONGITUDE, LATITUDE)%>%
  mutate(
    EVENT_DATE = as.Date(EVENT_DATE),  # ensure it's a Date
    month = month(EVENT_DATE)        # numeric month (1–12)
  )
# convert as sf points
df_conflict_sf <- st_as_sf(df_conflict, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)


# use above created grid  0.25° cells - a constant grid is comparable across maps 
prcp_grid <- prcp_polygons_all %>%
  filter(year == 2000, month == 1) %>%
  st_geometry() %>%
  st_sf() %>%
  mutate(grid_id = row_number())

# Spatial join - assign each conflict point to a grid cell
df_conflict_sf <- st_transform(df_conflict_sf, st_crs(prcp_grid))
joined <- st_join(prcp_grid, df_conflict_sf) %>% drop_na()

# conflicts per grid cell per year
conflict_counts <- joined %>%
  st_drop_geometry() %>%
  group_by(grid_id, YEAR, month) %>%
  summarise(conflict_count = n(), .groups = "drop")

# counts back to grid
grid_counts <- prcp_grid %>%
  left_join(conflict_counts, by = "grid_id") %>%
  mutate(conflict_count = replace_na(conflict_count, 0))%>%
  filter(between(YEAR, 2000, 2020))

# centroids for each polygon
df_conflict_count <- grid_counts %>%
  mutate(centroid = st_centroid(geometry)) %>%           
  mutate(
    longitude = st_coordinates(centroid)[, 1],          
    latitude = st_coordinates(centroid)[, 2]            
  ) %>%
  st_drop_geometry() %>%                                 
  select(longitude, latitude, everything())%>%
  select(longitude,latitude,YEAR,month,conflict_count)%>%
  rename_with(tolower)


df_prcp_anomal <- prcp_polygons_all%>%
  mutate(centroid = st_centroid(geometry)) %>%           
  mutate(
    longitude = st_coordinates(centroid)[, 1],          
    latitude = st_coordinates(centroid)[, 2]            
  ) %>%
  st_drop_geometry() %>%                                 
  select(longitude, latitude, everything())%>%
  select(longitude,latitude,year,month,prcp_anomaly)


df_temp_anomal <- temp_polygons_all%>%
  mutate(centroid = st_centroid(geometry)) %>%           
  mutate(
    longitude = st_coordinates(centroid)[, 1],          
    latitude = st_coordinates(centroid)[, 2]            
  ) %>%
  st_drop_geometry() %>%                                 
  select(longitude, latitude, everything())%>%
  select(longitude,latitude,year,month,temp_anomaly)



df_final <- df_prcp_anomal%>%
  left_join(df_temp_anomal,by = c("longitude","latitude","year","month"))%>%
  left_join(df_conflict_count, by = c("longitude","latitude","year","month"))%>%
  mutate(conflict_count = replace_na(conflict_count, 0))%>%    # Replace NA with 0
  mutate(conflict_binary = ifelse(conflict_count > 0, 1,0))%>%
  group_by(longitude,latitude)%>%
  mutate(pixel_id = cur_group_id())

# Table 1 :Estimation results for the primary model specifications using a mixed logit model.

# Prepare data with lags - because coflict may higgly corelated with paste climate evens 
df_lags <- df_final %>%
  arrange(pixel_id, year, month) %>%
  group_by(pixel_id) %>%
  mutate(
    temp_anomaly_lag1 = dplyr::lag(temp_anomaly, n = 1),
    prcp_anomaly_lag1 = dplyr::lag(prcp_anomaly, n = 1),
    temp_anomaly_lag2 = dplyr::lag(temp_anomaly, n = 2),
    prcp_anomaly_lag2 = dplyr::lag(prcp_anomaly, n = 2)
  ) %>%
  ungroup()

# logit model
model_lags <- feglm(
  conflict_binary ~ temp_anomaly + temp_anomaly_lag1 + temp_anomaly_lag2 +
    prcp_anomaly + prcp_anomaly_lag1 + prcp_anomaly_lag2 |
    pixel_id + year + month,
  data = df_lags,
  family = binomial(link = "logit")
)

# table
model_summary <- tibble(
  Variable = c("Temperature Anomaly (t)",
               "Temperature Anomaly (t-1)",
               "Temperature Anomaly (t-2)",
               "Precipitation Anomaly (t)",
               "Precipitation Anomaly (t-1)",
               "Precipitation Anomaly (t-2)"),
  Coefficient = sprintf("%.3f", coef(model_lags)),
  `Std. Error` = sprintf("%.3f", se(model_lags)),
  `Odds Ratio` = sprintf("%.3f", exp(coef(model_lags))),
  Significance = ifelse(
    summary(model_lags)$coeftable[,4] < 0.01, "***",
    ifelse(summary(model_lags)$coeftable[,4] < 0.05, "**",
           ifelse(summary(model_lags)$coeftable[,4] < 0.1, "*", "")))
)


# Save as LaTeX table


latex_table <- model_summary %>%
  kbl(
    caption = "Fixed Effects Logit Model of Conflict Probability",
    align = c("l", "r", "r", "r", "r"),
    col.names = c("Variable", "Coeff.", "S.E.", "Odds Ratio", ""),
    format = "latex",
    booktabs = TRUE,  # Professional quality tables
    linesep = "",     # Remove extra space between rows
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    position = "center",
    font_size = 10
  ) %>%
  column_spec(1, width = "4cm") %>%  # Wider first column
  column_spec(2:4, width = "1.5cm") %>%
  column_spec(5, width = "0.8cm") %>%
  footnote(
    general = c(
      "Note: Coefficients represent log-odds. Odds ratios obtained by exponentiating coefficients.",
      "Fixed effects for pixel, year, and month included in all models.",
      "*** p < 0.01, ** p < 0.05, * p < 0.1"
    ),
    threeparttable = TRUE,  # Proper footnote formatting
    escape = FALSE
  )

# Save to .tex file
writeLines(as.character(latex_table), "Tables/model_summary_table.tex")

# For direct PDF output (requires tinytex/LaTeX installation)
# save_kable(latex_table, "model_summary_table.pdf")
model_summary

# Figure 6: Probability of conflict under temperature and precipitation anomalies.

# conflict odd ration during changes in future anomalies. since our model asess the lage efefct I use pase events
# to predict changes in 2021 
# Real observed data from 2018 onward
# Artificial scenarios for 2021 with extreme heat (temp_anomaly) and normal precipitation
# first Temp anomaly 
pred_data <- df_lags %>%
  filter(year >= 2018) %>%
  select(pixel_id, year, month, temp_anomaly, prcp_anomaly) %>%
  bind_rows(
    expand.grid(
      pixel_id = unique(df_lags$pixel_id),
      year = 2021,
      month = 1:12,
      temp_anomaly = 1,
      prcp_anomaly = 10
    )
  ) %>%
  group_by(pixel_id) %>%
  arrange(year, month) %>%
  mutate(
    temp_anomaly_lag1 = lag(temp_anomaly, 1),
    temp_anomaly_lag2 = lag(temp_anomaly, 2),
    prcp_anomaly_lag1 = lag(prcp_anomaly, 1),
    prcp_anomaly_lag2 = lag(prcp_anomaly, 2)
  ) %>%
  ungroup() %>%
  filter(year == 2021)

# account for fixed effect - we use fixed effct to control all othe unobservalbe factors
last_year <- max(df_lags$year)
last_month <- df_lags %>% 
  filter(year == last_year) %>% 
  pull(month) %>% max()

pred_data <- pred_data %>%
  left_join(
    df_lags %>%
      filter(year == last_year, month == last_month) %>%
      distinct(pixel_id, .keep_all = TRUE) %>%
      select(pixel_id, year_fe = year, month_fe = month),
    by = "pixel_id"
  )

# logit model 
model_lags <- feglm(
  conflict_binary ~ temp_anomaly + temp_anomaly_lag1 + temp_anomaly_lag2 +
    prcp_anomaly + prcp_anomaly_lag1 + prcp_anomaly_lag2 |
    pixel_id + year + month,
  data = df_lags,
  family = binomial(link = "logit")
)

replacement_year <- max(df_lags$year)

pred_data <- pred_data %>%
  mutate(year_for_fe = replacement_year, .after = year)

pred_data$prob <- predict(
  model_lags,
  newdata = pred_data %>% mutate(year = year_for_fe),
  type = "response"
)

# prediction data
pred_data <- pred_data %>%
  select(-year_for_fe) %>%
  mutate(
    prediction_year = year,
    fe_reference_year = replacement_year
  )


pred_data <- pred_data%>%
  select(pixel_id,prob)

pixel <- df_lags%>%
  filter(year==2020 & month ==1)%>%
  select(pixel_id,longitude,latitude)%>%
  left_join(pred_data)%>%
  st_as_sf(
    coords = c("longitude", "latitude"),  
    crs = 4326,                          
    remove = FALSE                        
  )

grid <- prcp_polygons_all%>%
  filter(year==2020 & month ==1)

# same CRS
pixel <- st_transform(pixel, st_crs(grid))

# spatial join 
joined_data <- st_join(grid, pixel, join = st_contains)

joined_data <- joined_data%>%
  group_by(pixel_id)%>%
  mutate(ave_prob = mean(prob))

# same CRS
st_crs(joined_data) <- st_crs(kenya)  # Match coordinate systems

# clip pixels to kenya's boundaries
joined_data <- st_intersection(joined_data, kenya)

joined_data <- joined_data %>%
  drop_na() 

conflit_prob <- ggplot(joined_data,) +
  geom_sf(aes(fill = prob), color = NA, alpha = 0.8) +
  geom_sf(data = kenya, fill = NA, color = "black", linewidth = 0.3) +
  scale_fill_distiller(
    name = "Conflict Probability",
    palette = "Spectral",  # Or "RdYlBu", "PuBuGn", etc.
    direction = -1,
    limits = c(0, 1),
    na.value = "white"
  )+
  labs(
    title = "",
    subtitle = "",
    caption = ""
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )


# Save to PNG
ggsave(
  filename = "Figures/conflit_prob.png",
  plot = conflit_prob,  # Explicitly specify plot object
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

conflit_prob

# Quality Check : correlation matrix

#var-cor matrix
numeric_vars <- df_lags %>% 
  select(temp_anomaly, temp_anomaly_lag1, temp_anomaly_lag2,
         prcp_anomaly, prcp_anomaly_lag1, prcp_anomaly_lag2)

# Compute correlation matrix
cor_matrix <- cor(numeric_vars, use = "complete.obs")  # Handles missing data
#print(cor_matrix)


corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black")

# vif model

model <- lm(conflict_binary ~ temp_anomaly + temp_anomaly_lag1 + temp_anomaly_lag2 +
              prcp_anomaly + prcp_anomaly_lag1 + prcp_anomaly_lag2, 
            data = df_lags)
vif(model)
