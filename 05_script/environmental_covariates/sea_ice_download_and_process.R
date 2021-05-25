#==============================================================================#
#                                                                              #
#                  Sea ice data batch download & processing                    #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(RCurl)
library(raster)
library(ncdf4)
Sys.setenv(LANG = "en")



# A. Download sea ice data =====================================================

# + 1. Create a directory for the sea ice files of each year -------------------

years <- seq(from = 1991, to = 2020, by = 1)

for (k in 1:length(years)) {
  dir.create(path = paste0("04_raw_data/sea_ice/hamburg/", years[k], "/"))
}

# + 2. Create a file within each folder for the URL list -----------------------

for (k in 1:length(years)) {
  dir.create(path = paste0("04_raw_data/sea_ice/hamburg/", years[k], "/",
                           "list_filenames", "/"))
}

# + 3. Create the list of URLs (one per year) ----------------------------------

url_base <- "ftp://ftp-icdc.cen.uni-hamburg.de/asi_ssmi_iceconc/arc/"

url_list <- c()
for (k in 1:length(years)) {
  url_list <- c(url_list, 
                paste0(url_base, years[k], "/"))
}


# + 4. Retrieve the online file name and save it into yearly folder ------------

for (k in 1:length(url_list)) {   # = length(years)
  filenames_k <- getURL(url_list[k], 
                         ftp.use.epsv = FALSE, 
                         dirlistonly = TRUE) %>%
    strsplit(split = "\r\n") %>%
    unlist() %>%
    sort()
  # Save the list of URLs
  saveRDS(filenames_k, file = paste0("04_raw_data/sea_ice/hamburg/", years[k], "/",
                                      "list_filenames/",                    # name of the folder
                                      "list_filenames_per_day_", years[k], ".rds")) # Name of the file

}



# + 5. Download and save the .nc files -----------------------------------------

# In the lab (with the eduroam wifi), it takes ~55min to download sea ice .nc files 
# for a year
# At home, it takes ~30min

years <- seq(from = 1991, to = 2020, by = 1)
url_base <- "ftp://ftp-icdc.cen.uni-hamburg.de/asi_ssmi_iceconc/arc/"

start <- Sys.time()
for (i in 1:1) { #length(years)) {
  filenames.year.i <- readRDS(file = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/",
                                            "list_filenames/",                    # name of the folder
                                            "list_filenames_per_day_", years[i], ".rds"))
  
  for (j in 1:length(filenames.year.i)) {  
    download.file(url = paste0(url_base, years[i], "/", filenames.year.i[j]),
                  destfile = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/", 
                                   filenames.year.i[j]),
                  mode = "wb")
  }
}

end <- Sys.time()




# + 6. Download and save sea ice data for 1990-1991 ----------------------------

filenames_90_91_north_south <- read_csv("04_raw_data/sea_ice/NSIDC/sea_ice_1990-1991.txt") %>%
  filter(nchar(url) == 90) %>% # Delete the filenames that don't seam to correspond to daily sea ice concentration
  arrange(url) 

index_south_pole <- grep(pattern = "v1.1_s.bin", x = filenames_90_91_north_south$url) 

filenames_90_91 <- filenames_90_91_north_south[-index_south_pole, ] # Delete filenames corresponding to the south pole

k = 1
download.file(url = filenames_90_91$url[k],
              destfile = paste0("04_raw_data/sea_ice/NSIDC/1990/test.bin"),
              mode = "wb")






# B. process the sea ice data ==================================================

# Download, crop, stack, and save the files ------------------------------------

# This code includes the code to create a matrix containing the raster data to 
# accelerate the calculations of sea ice metrics later on

library(ff)
`%notin%` = Negate(`%in%`)
years <- seq(from = 1991, to = 2020, by = 1)
crop.square <- extent(-1000000, 2500000, -2000000, 1000000)

# 50 seconds per year 
start <- Sys.time()
for (i in 1:length(years)) {
  
  sea_ice_files <- sort(list.files(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/"))) # Get the list of all the .nc files for the year i
  sea_ice_files <- sea_ice_files[sea_ice_files %notin% c("list_filenames", "raster_stack")] # Remove the name of folders from the list of .nc files
  
  sea_ice_raster <- list()
  for (j in 1:(length(sea_ice_files))) { 
    print(j)
    ras_sub <- raster(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/", sea_ice_files[j]), 
                      ncdf = TRUE, 
                      varname = "sea_ice_area_fraction")
    
    # print("Loading raster")
    
    ss <- as.matrix(ras_sub)
    sss <- raster(ss, xmn = -3850000, xmx = 3700000, ymn = -5350000, ymx = 5840000)
    # sss <- raster(ss, xmn = -3850, xmx = 3750, ymn = -5350, ymx = 5850)
    projection(sss) <- CRS("+init=EPSG:3413")
    ssss <- raster::crop(sss, crop.square)
    
    sea_ice_raster[[j]] <- ssss
  }
  
  # Stack the layers
  sea_ice_raster_stack <- raster::stack(sea_ice_raster)
  # Assign the name of each layer of the stack
  names(sea_ice_raster_stack) <- paste0("sea_ice_", substring(sea_ice_files, 1, 4), "_",
                                        substring(sea_ice_files, 5, 6), "_",
                                        substring(sea_ice_files, 7, 8))
  
  # Create a directory to save the raster stack
  dir.create(path = paste0("04_raw_data/sea_ice/hamburg/",               
                           years[i], "/raster_stack/"))
  
  # Save the raster stack
  writeRaster(sea_ice_raster_stack, filename = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
                                                      "sea_ice_raster_stack_", years[i], ".tif"),
              overwrite = TRUE) # File name
  save(sea_ice_raster_stack, file = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
                                           "sea_ice_raster_stack_", years[i], ".RData")) # File name
  
  # Create a matrix stored on disk to accelerate the calculations later on
  mat <- ff(vmode = "double",
            dim = c(ncell(sea_ice_raster_stack),
                    nlayers(sea_ice_raster_stack)),
            filename = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/",
                              "sea_ice_matrix_stack_", years[i], ".ffdata"))
  for(j in 1:nlayers(sea_ice_raster_stack)) {
    mat[, j] <- sea_ice_raster_stack[[j]][]
  }
  save(mat, file = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
                          "sea_ice_matrix_stack_", years[i], ".RData"))
  save(mat, file = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
                          "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  rm(sea_ice_raster_stack, ssss, sss, ss, ras_sub, mat)
  
}
end <- Sys.time()















# C. Process 1990-1991 sea ice data ============================================

test <- readBin(con = "04_raw_data/sea_ice/NSIDC/1990-1991/130228830/nt_19891231_f08_v1.1_n.bin",
                what = "integer",
                n = 600000)
?readBin
column.names <- readBin(read.filename, character(),  n = 3)
04_raw_data\sea_ice\NSIDC\1990-1991\130228830


ncell(ras_sub)





















# OLD --------------------------------------------


`%notin%` = Negate(`%in%`)
years <- seq(from = 2000, to = 2019, by = 1)
crop.square <- extent(-1000000, 2500000, -2000000, 1000000)


# 40 seconds per year 
start <- Sys.time()
for (i in 1:1) { #length(years)) {
  
  sea_ice_files <- sort(list.files(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/"))) # Get the list of all the .nc files for the year i
  sea_ice_files <- sea_ice_files[sea_ice_files %notin% c("list_filenames", "raster_stack")] # Remove the name of folders from the list of .nc files

  sea_ice_raster <- list()
  for (j in 1:(length(sea_ice_files))) {
    print(j)
    ras_sub <- raster(paste0("04_raw_data/sea_ice/hamburg/", years[i], "/", sea_ice_files[j]), 
                      ncdf = TRUE, 
                      varname = "sea_ice_area_fraction")
    
    print("Loading raster")
    
    ss <- as.matrix(ras_sub)
    sss <- raster(ss, xmn = -3850000, xmx = 3700000, ymn = -5350000, ymx = 5840000)
    # sss <- raster(ss, xmn = -3850, xmx = 3750, ymn = -5350, ymx = 5850)
    projection(sss) <- CRS("+init=EPSG:3413")
    ssss <- raster::crop(sss, crop.square)
    
    sea_ice_raster[[j]] <- ssss
  }
  
  # Stack the layers
  sea_ice_raster_stack <- raster::stack(sea_ice_raster)
  # Assign the name of each layer of the stack
  names(sea_ice_raster_stack) <- paste0("sea_ice_", substring(sea_ice_files, 1, 4), "_",
                                        substring(sea_ice_files, 5, 6), "_",
                                        substring(sea_ice_files, 7, 8))
  
  # Create a directory to save the raster stack
  dir.create(path = paste0("04_raw_data/sea_ice/hamburg/",               
                           years[i], "/raster_stack/"))
  
  # Save the raster stack
  writeRaster(sea_ice_raster_stack, filename = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
                                                      "sea_ice_raster_stack_", years[i], ".tif")) # File name
  save(sea_ice_raster_stack, file = paste0("04_raw_data/sea_ice/hamburg/", years[i], "/raster_stack/", # Folder name
                                           "sea_ice_raster_stack_", years[i], ".RData")) # File name
  rm(sea_ice_raster_stack, ssss, sss, ss, ras_sub)
}
end <- Sys.time()




# This is the code I used to try and determine the right limits to attribute to 
# the sea ice rasters, so that it is well ajusted compared to the bathymetry data 
# and the subpopulation boundary data.

# Load a sea ice raster to do the trials
years <- seq(from = 1993, to = 2019, by = 1)
sea_ice <- sort(list.files(paste0("04_raw_data/sea_ice/hamburg/", years[1], "/")))

ras_sub <- raster(paste0("04_raw_data/sea_ice/hamburg/", years[1], "/", sea_ice[1]), 
                  ncdf = TRUE, 
                  varname = "sea_ice_area_fraction")

# Try different limits to the raster
# these are the limits attributed originally in Marie-Anne's code covnerted from
# km to m
xmin <- -3850000  
xmax <- 3850000
ymin <- -5350000
ymax <- 5850000

# These limits seem to work best
xmin <- -3850000
xmax <- 3700000
ymin <- -5350000
ymax <- 5840000

# Test the limits 
bb <- extent(xmin, xmax, ymin, ymax)
extent(ras_sub) <- bb

crop.square <- extent(-1000000, 2500000, -2000000, 1000000)
ras_sub <- raster::crop(ras_sub, crop.square)

plot(ras_sub)
plot(barents.sea.reproj,
     add = TRUE)

# Plot the bathymetry data with the polygon to compare (we know the limits of 
# the bathymetry data are right, so the way the bathemetry map is ajusted compared
# to the subpopulation polygon can serve as a reference)
bathymetry <- raster("06_processed_data/bathymetry_data/ICBAO_bathymetry_polar_cropped.tif")
plot(bathymetry)
plot(barents.sea.reproj,
     add = TRUE)





# barents.sea.reproj.df <- tidy(barents.sea.reproj)
# df_sea_ice_sep <- as.data.frame(ras_sub, 
#                                 xy = TRUE)
# 
# ggplot() +
#   geom_raster(data = df_sea_ice_sep, aes(x = x, y = y, fill = sea.ice.concentration)) +
#   geom_polygon(data = barents.sea.reproj.df, aes(x = long, y = lat, group = group), 
#                color = "red", alpha = 0.25) +
#   theme_classic() +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#         axis.text = element_text(size = 7),
#         legend.title = element_text(size = 7),
#         legend.position = "bottom",
#         legend.text = element_text(size = 7)) +
#   
#   labs(x = "",
#        y = "",
#        fill = "Sea ice \ncover") +
#   scale_fill_continuous(low = "#132b43", high = "white", 
#                         guide = "colorbar", na.value = "#cccccc",
#                         breaks = c(0, 50, 100),
#                         labels = c("0%", "50%", "100%"))
# 
# 
# range(df_sea_ice_sep$x)
# range(df_sea_ice_sep$y)
# ggsave("banquise_polygone_still_wrong_1.png",
#        width = 10, height = 15.9, unit = "cm")
# 
# 
# # bathymetrly + polygone
# barents.sea.reproj.df <- tidy(barents.sea.reproj)
# 
# 
# bathymetry <- raster("06_processed_data/bathymetry_data/ICBAO_bathymetry_polar_cropped.tif")
# 
# bathymetry <- aggregate(bathymetry, fact = 8)
# df_bathymetry <- as.data.frame(bathymetry, 
#                                xy = TRUE)
# 
# ggplot() +
#   geom_raster(data = df_bathymetry, aes(x = x, y = y, fill = ICBAO_bathymetry_polar_cropped)) +
#   geom_polygon(data = barents.sea.reproj.df, aes(x = long, y = lat, group = group), 
#                color = "red", alpha = 0.25) +
#   theme_classic() +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#         axis.text = element_text(size = 7),
#         legend.title = element_text(size = 7),
#         legend.position = "bottom",
#         legend.text = element_text(size = 7)) +
#   
#   labs(x = "",
#        y = "",
#        fill = "Depth") +
#   scale_fill_continuous(low = "black", high = "white", #"#132b43"
#                         guide = "colorbar", na.value = "#cccccc")
# 
# 
# x <- (range(df_bathymetry$x)[2] - range(df_bathymetry$x)[1])
# y <- (range(df_bathymetry$y)[2] - range(df_bathymetry$y)[1])
# max(c(x, y))/min(c(x, y))
# x>y
# range(df_bathymetry$y)
# ggsave("bathymetrie_polygone_1.png",
#        width = 10, height = 14.3, unit = "cm")









