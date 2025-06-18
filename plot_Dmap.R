library(terra)
library(geodata)
library(sf)
#mapping ranges2
#download Little's range map shapefiles from https://github.com/wpetry/USTreeAtlas

qubi<-vect("<path_to_download>/querbico.shp")
quly<-vect("<path_to_download>/querlyra.shp")

quly_ig<-rast("<path_to_QULY_binary_SDMmodelProjection>/QULY_IG_binary.tif")
qubi_ig<-rast("<path_to_QUBI_binary_SDMmodelProjection>/QUBI_IG_binary.tif")
ig_overlap <- qubi_ig * 1 + quly_ig * 2  
us_states <- gadm("USA", level=1, path=tempdir())  # Level 1 = states
us_sf <- st_as_sf(us_states)
can<-gadm("Canada", level=0, path=tempdir())
can_sf <- st_as_sf(can)
raster_extent <- st_as_sfc(st_bbox(c(xmin=-100, xmax=-69, ymin=24, ymax=49), crs=st_crs(us_sf)))
# Clip US states to raster extent
can_sf_clipped <- st_crop(can_sf, raster_extent)
us_sf_clipped <- st_crop(us_sf, raster_extent)
single_sf <- dplyr::bind_rows(list(us_sf_clipped,can_sf_clipped))
dissolve_sf <- st_union(single_sf)
ig_overlap<-crop(ig_overlap, dissolve_sf)
df <- as.data.frame(ig_overlap, xy=TRUE)
colnames(df) <- c("x", "y", "presence")

dna<-read.csv("~/data/Dsuite_df.csv")
library(dplyr)
dna<-dplyr::select(dna, c("names","D","lon","lat"))
dnav<-vect(dna,crs="EPSG:4326" , c('lon','lat'))
# Convert SpatVector to sf
dnav_sf <- st_as_sf(dnav)

library(sf)
library(dplyr)
library(ggplot2)
library(ggnewscale)
# Adjusted bar size and more diagonal line layout
# Prepare data with more spacing for diagonal lines
dnav_df <- dnav_sf %>%
  mutate(X = st_coordinates(.)[,1],
         Y = st_coordinates(.)[,2]) %>%
  arrange(desc(Y)) %>%
  mutate(bar_id = factor(names, levels = names),
         bar_y = max(Y),
         bar_x = max(X) + 15,   # Push bars farther to the right for better diagonal lines
         bar_width = D*25)
dnav_df$bar_y<-c(seq(48.5,40.5,-1), seq(36,25,-1))

# Add species and source info as data frames for legend mapping
quly_sf<-st_as_sf(quly)
qubi_sf<-st_as_sf(qubi)
quly_sf$species <- "Q. lyrata"
qubi_sf$species <- "Q. bicolor  "
dnav_df$source <- "Primary Samples"
library(ggnewscale)  # Make sure this is loaded

# Add a new column for separate fill legend
dnav_df$source_fill <- "Primary Samples"

map <- ggplot() +
  geom_sf(data = quly_sf, aes(fill = species), color = 'black', alpha = 0.75) +
  geom_sf(data = qubi_sf, aes(fill = species), color = 'black', alpha = 0.75) +
  scale_fill_manual(name = "", values = c("Q. lyrata" = "#440154FF","Q. bicolor  " = "#FDE725FF"), labels=c(expression(italic("Q. bicolor  ")),expression(italic("Q. lyrata  "))))+
  ggnewscale::new_scale_fill() +  # Reset fill scale for points
  # Base map
  geom_sf(data = us_sf_clipped, fill = NA, color = "grey50", size = 0.25) +
  
  # Points (Primary Samples)
  geom_point(data = dnav_df, aes(x = X, y = Y, shape = source, fill = source_fill),
             size = 3, colour = 'black') +
  
  # Diagonal lines
  geom_segment(data = dnav_df,
               aes(x = X, y = Y, xend = bar_x - bar_width, yend = bar_y),
               color = "black", linetype = "dashed") +
  
  # Horizontal bars
  geom_rect(data = dnav_df,
            aes(xmin = bar_x - bar_width,
                xmax = bar_x,
                ymin = bar_y - 0.6,
                ymax = bar_y + 0.6),
            fill = "#21908CFF", linewidth = 0.5, colour = 'black') +
  
  # Bar labels
  geom_text(data = dnav_df,
            aes(x = bar_x + 0.5, y = bar_y, label = round(D, 2)),
            hjust = 0, size = 5, family = 'mono') +
  
  # Legends
  scale_fill_manual(name = "", values = c("Primary Samples" = "#21908CFF")) +
  scale_shape_manual(name = "", values = c("Primary Samples" = 24)) +
  guides(
    fill = guide_legend(override.aes = list(shape = 24, fill = "#21908CFF", colour = "black")),
    shape = "none"
  ) +
  annotate("text",
             x = max(dnav_df$bar_x) + 5,    # a bit to the right of the bars
             y = mean(range(dnav_df$bar_y)), # vertically centered
             label = "Patterson's D Statistics by Sampling Site",
             angle = 270,
             size = 6,
             family = "mono",
             hjust = 0.5,
             vjust = 1.2)+
  coord_sf(clip = "off") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.y = element_blank(),
    text = element_text(family = "mono"),
    legend.text = element_text(family = "mono", size=18),
    legend.title = element_text(family = "mono")
  )
map

# Save the final plot
ggsave("fig5.pdf", plot = map, width = 8, height = 8, dpi = 800)
ggsave("fig5.png", plot = map, width = 8, height = 8, dpi = 800, bg="white")

