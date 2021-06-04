#==============================================================================#
#                                                                              #
#                  Sea ice data batch download & processing                    #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(RCurl)
library(raster)
library(ncdf4)
library(ff)
Sys.setenv(LANG = "en")



# A. Download sea ice data =====================================================

# + 1. Create a directory for the sea ice files of each year -------------------

years <- seq(from = 1991, to = 2020, by = 1)

for (k in 1:length(years)) {
  dir.create(path = paste0("04_raw_data/sea_ice/", years[k], "/"))
}

# + 2. Create a file within each folder for the URL list -----------------------

for (k in 1:length(years)) {
  dir.create(path = paste0("04_raw_data/sea_ice/", years[k], "/",
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
  saveRDS(filenames_k, file = paste0("04_raw_data/sea_ice/", years[k], "/",
                                      "list_filenames/",                    # name of the folder
                                      "list_filenames_per_day_", years[k], ".rds")) # Name of the file

}



# + 5. Download and save the .nc files -----------------------------------------

# In the lab (with the eduroam wifi), it takes ~55min to download sea ice .nc files 
# for a year
# At home, it takes ~30min

years <- seq(from = 1990, to = 2021, by = 1)
url_base <- "ftp://ftp-icdc.cen.uni-hamburg.de/asi_ssmi_iceconc/arc/"

start <- Sys.time()
for (i in 32:32) { #length(years)) {
  filenames.year.i <- readRDS(file = paste0("04_raw_data/sea_ice/", years[i], "/",
                                            "list_filenames/",                    # name of the folder
                                            "list_filenames_per_day_", years[i], ".rds"))
  
  for (j in 1:length(filenames.year.i)) {  
    download.file(url = paste0(url_base, years[i], "/", filenames.year.i[j]),
                  destfile = paste0("04_raw_data/sea_ice/", years[i], "/", 
                                   filenames.year.i[j]),
                  mode = "wb")
  }
}

end <- Sys.time()






# B. process the sea ice data ==================================================

# Download, crop, stack, and save the files 

# This code includes the code to create a matrix containing the raster data to 
# accelerate the calculations of sea ice metrics later on

`%notin%` = Negate(`%in%`)
years <- seq(from = 1992, to = 2021, by = 1)
crop.square <- extent(-1000000, 2500000, -2000000, 1000000)

# 50 seconds per year 
start <- Sys.time()
for (i in 30:length(years)) {
  
  sea_ice_files <- sort(list.files(paste0("04_raw_data/sea_ice/", # Get the list of all the .nc files for the year i
                                          years[i], "/"))) 
  sea_ice_files <- sea_ice_files[sea_ice_files %notin% c("list_filenames", 
                                                         "raster_stack")] # Remove the name of folders or files other than .nc files
                                                                          # from the list of .nc files
  
  sea_ice_raster <- list()
  for (j in 1:(length(sea_ice_files))) { 
    print(j)
    ras_sub <- raster(paste0("04_raw_data/sea_ice/", years[i], "/", sea_ice_files[j]), 
                      ncdf = TRUE, 
                      varname = "sea_ice_area_fraction")
    
    # print("Loading raster")
    
    ss <- as.matrix(ras_sub)
    sss <- raster(ss, xmn = -3850000, xmx = 3700000, ymn = -5350000, ymx = 5840000)
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
  dir.create(path = paste0("04_raw_data/sea_ice/",               
                           years[i], "/raster_stack/"))
  
  # Save the raster stack
  writeRaster(sea_ice_raster_stack, filename = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                                                      "sea_ice_raster_stack_", years[i], ".tif"),
              overwrite = TRUE) # File name
  save(sea_ice_raster_stack, file = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                                           "sea_ice_raster_stack_", years[i], ".RData")) # File name
  
  # Create a matrix stored on disk to accelerate the calculations later on
  mat <- ff(vmode = "double",
            dim = c(ncell(sea_ice_raster_stack),
                    nlayers(sea_ice_raster_stack)),
            filename = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/",
                              "sea_ice_matrix_stack_", years[i], ".ffdata"))
  for(j in 1:nlayers(sea_ice_raster_stack)) {
    mat[, j] <- sea_ice_raster_stack[[j]][]
  }
  save(mat, file = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                          "sea_ice_matrix_stack_", years[i], ".RData"))
  save(mat, file = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                          "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  rm(sea_ice_raster_stack, ssss, sss, ss, ras_sub, mat)
  
}
end <- Sys.time()















# C. Process 1990-1991 sea ice data ============================================

# These files were uploaded from the NSIDC. I didn't use an R code to download
# them contrary to the sea ice files from Hamburg.

# + 1. Get list of all the folders ----------------------------------------------
folder_list <- data.frame(directory = list.dirs(path = "04_raw_data/sea_ice/NSIDC/1990-1991/", 
                                                full.names = TRUE, recursive = TRUE)) %>%
  filter(nchar(directory) == 46)


# + 2. Get the name of the sea ice file within each of these folders -----------
list_file_names <- c()
list_file_directories <- c()
for (i in 1:nrow(folder_list)) {
  files_i <- list.files(path = folder_list$directory[i])
  bin_SI_i <- files_i[nchar(files_i) == 26] # Keep only the .bin files (there are other files)
                                            # There are also folders that contain .bin files with 
                                            # an incomplete date (e.g. nt_199008_f08_v1.1_n.bin)
                                            # The names of these files are shorter than 26 characters
                                            # and will thus not be added to the list of files to
                                            # convert to rasters.
  if (length(bin_SI_i) == 1) {
    list_file_directories <- c(list_file_directories, 
                               paste0(folder_list$directory[i], "/",  bin_SI_i))
    list_file_names <- c(list_file_names, bin_SI_i)
    
    
  }
}


list_file_directories <- data.frame(directory = list_file_directories,
                                    file_name = list_file_names) %>%
  filter(substr(directory, 51, 54) != "1989") %>% # Remove file corresponding to 31/12/1989
  arrange(file_name)

x <- raster("04_raw_data/sea_ice/NSIDC/1990-1991//130233175/nt_19900101_f08_v1.1_n.bin")


# + 3. Convert to raster, crop, and stack the rasters by year ------------------

# I downloaded this function from https://github.com/cran/raster/blob/master/R/nsidcICE.R
get_raster_from_NSIDC_file <- function(x) {
  ## check name structure
  ## "nt_19781119_f07_v01_s.bin"
  
  bx <- basename(x)
  ## test that we can get a date from this
  ## (as POSIXct so that Z-comparisons are more natural)
  dts <- as.POSIXct(basename(x), format = "nt_%Y%m%d", tz = "UTC")
  ## test that we see _f and _v
  fyes <- tolower(substr(bx, 13L, 13L)) %in% c("f", "n")
  vyes <- tolower(substr(bx, 17L, 17L)) %in% c("v", "n")
  
  ## finally, it's north or south
  hemi <- tolower(substr(bx, 22L, 22L))
  hyes <- hemi %in% c("s", "n")
  if(!(!is.na(dts) & fyes & vyes & hyes)) return(NULL)
  
  ## NSIDC projection and grid size
  ## https://nsidc.org/data/polar_stereo/ps_grids.html
  ## http://spatialreference.org/ref/?search=nsidc
  ## Hughes 1980 ellipsoid, True Scale Lat is +/-70
  
  if (hemi == "s") {
    prj <-  "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"
    
    dims <- c(316L, 332L)
    ext <- c(-3950000, 3950000, -3950000, 4350000)
  } else {
    ## northern hemisphere
    prj <- "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs"
    dims <- c(304, 448)
    # ext <- c(-3837500, 3762500, -5362500, 5837500)
    ext <- c(-3850000, 3700000, -5350000, 5840000)
  }
  on.exit(close(con))
  con <- file(x, open = "rb")
  
  ## chuck the header
  try1 <- try(trash <- readBin(con, "integer", size = 1, n = 300))
  ## TODO: warnings that we thought it was NSIDC, but it did not work?
  if (inherits(try1, "try-error")) return(NULL)
  dat <- try(readBin(con, "integer", size = 1, n = prod(dims), endian = "little", signed = FALSE))
  if (inherits(dat, "try-error")) return(NULL)
  
  r100 <- dat > 250
  r0 <- dat < 1
  ##      if (rescale) {
  dat <- dat/2.5  ## rescale back to 100
  ##      }
  ##      if (setNA) {
  dat[r100] <- NA
  ## dat[r0] <- NA
  ##      }
  r <- raster(t(matrix(dat, dims[1])), 
              xmn = ext[1], xmx = ext[2], 
              ymn = ext[3], ymx = ext[4], 
              crs = CRS(prj))
  
  setZ(r, dts, name = "time")
  
}


# Run everything
years <- c(1990, 1991)
crop.square <- extent(-1000000, 2500000, -2000000, 1000000)

start <- Sys.time()
for (i in 1:length(years)) {
  list_files_year_i <- list_file_directories %>%
    filter(as.numeric(substr(file_name, 4, 7)) == years[i])
  
  # Convert to a raster a store all the rasters in a list
  sea_ice_raster <- list()
  for (j in 1:nrow(list_files_year_i)) {
    raster_j <- get_raster_from_NSIDC_file(list_files_year_i$directory[j])
    raster_j_cropped <- raster::crop(raster_j, crop.square)
    sea_ice_raster[[j]] <- raster_j_cropped
    print(j)
  }
  
  # Stack the layers
  sea_ice_raster_stack <- raster::stack(sea_ice_raster)
  # Assign the name of each layer of the stack
  names(sea_ice_raster_stack) <- paste0("sea_ice_", 
                                        substring(list_files_year_i$file_name, 4, 7), "_",
                                        substring(list_files_year_i$file_name, 8, 9), "_",
                                        substring(list_files_year_i$file_name, 10, 11))
  # Create a directory to save the raster stack
  dir.create(path = paste0("04_raw_data/sea_ice/",               
                           years[i], "/raster_stack/"))
  
  # Save the raster stack
  writeRaster(sea_ice_raster_stack, filename = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                                                      "sea_ice_raster_stack_", years[i], ".tif"),
              overwrite = TRUE) # File name
  save(sea_ice_raster_stack, file = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                                           "sea_ice_raster_stack_", years[i], ".RData")) # File name
  
  # Create a matrix stored on disk to accelerate the calculations later on
  mat <- ff(vmode = "double",
            dim = c(ncell(sea_ice_raster_stack),
                    nlayers(sea_ice_raster_stack)),
            filename = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/",
                              "sea_ice_matrix_stack_", years[i], ".ffdata"))
  
  for(j in 1:nlayers(sea_ice_raster_stack)) {
    mat[, j] <- sea_ice_raster_stack[[j]][]
  }
  save(mat, file = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                          "sea_ice_matrix_stack_", years[i], ".RData"))
  save(mat, file = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                          "sea_ice_matrix_stack_2_", years[i], ".ffdata"))
  
  rm(sea_ice_raster_stack, raster_j, raster_j_cropped, mat)
  
}
end <- Sys.time()


# raster_1990 <- raster::stack("04_raw_data/sea_ice/1990/raster_stack/sea_ice_raster_stack_1990.tif")
# raster_1990_01_01 <- raster_1990@layers[[1]]
# 
# raster_1992 <- raster::stack("04_raw_data/sea_ice/1992/raster_stack/sea_ice_raster_stack_1992.tif")
# raster_1992_01_01 <- raster_1992@layers[[1]]
# 
# plot(raster_1990_01_01)
# plot(raster_1992_01_01)












# OLD --------------------------------------------


`%notin%` = Negate(`%in%`)
years <- seq(from = 2000, to = 2019, by = 1)
crop.square <- extent(-1000000, 2500000, -2000000, 1000000)


# 40 seconds per year 
start <- Sys.time()
for (i in 1:1) { #length(years)) {
  
  sea_ice_files <- sort(list.files(paste0("04_raw_data/sea_ice/", years[i], "/"))) # Get the list of all the .nc files for the year i
  sea_ice_files <- sea_ice_files[sea_ice_files %notin% c("list_filenames", "raster_stack")] # Remove the name of folders from the list of .nc files

  sea_ice_raster <- list()
  for (j in 1:(length(sea_ice_files))) {
    print(j)
    ras_sub <- raster(paste0("04_raw_data/sea_ice/", years[i], "/", sea_ice_files[j]), 
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
  dir.create(path = paste0("04_raw_data/sea_ice/",               
                           years[i], "/raster_stack/"))
  
  # Save the raster stack
  writeRaster(sea_ice_raster_stack, filename = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                                                      "sea_ice_raster_stack_", years[i], ".tif")) # File name
  save(sea_ice_raster_stack, file = paste0("04_raw_data/sea_ice/", years[i], "/raster_stack/", # Folder name
                                           "sea_ice_raster_stack_", years[i], ".RData")) # File name
  rm(sea_ice_raster_stack, ssss, sss, ss, ras_sub)
}
end <- Sys.time()




# This is the code I used to try and determine the right limits to attribute to 
# the sea ice rasters, so that it is well ajusted compared to the bathymetry data 
# and the subpopulation boundary data.

# Load a sea ice raster to do the trials
years <- seq(from = 1993, to = 2019, by = 1)
sea_ice <- sort(list.files(paste0("04_raw_data/sea_ice/", years[1], "/")))

ras_sub <- raster(paste0("04_raw_data/sea_ice/", years[1], "/", sea_ice[1]), 
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









