#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                              #
#                     Data processing + visualization                          #
#                                                                              #
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

library(tidyverse)
library(lubridate)
Sys.setenv(LANG = "en")


# A. Clean data ================================================================

# Select relevant columns
capture_data_raw <- read_delim("04_raw_data/polar_bears/capturedata2020.csv", 
                              ";", escape_double = FALSE, trim_ws = TRUE)
colnames(capture_data_raw)


# Process columns
capture_data_all <- capture_data_raw[, c(2:4, 6:7, 11:12, 28:29, 39:41, 45, 72:73)] %>%
  
  mutate(date = as.Date(DATO, format = "%d / %m / %Y"),
         day_number = yday(date),
         lat = str_replace(LAT, ",", "."),
         lon = str_replace(LON, ",", ".")) %>%
  
  dplyr::select(date, 
                year = YEAR,
                day_number,
                ID_NR = `ID NR`,
                sex = SEX,
                merking = MERKING,
                age_class = kohort,
                age_for_analyses = `Age for analyses`,
                lat, long = LON,
                z_length = `Z LENGTH`,
                s_length = `S LENGTH`,
                girth = GIRTH,
                condition = CONDITION,
                cub_status = `Cub status`,
                cub_number = `Cub number`)


# Process rows

# 1. Only females with at least 1 cub
capture_data_females_w_cubs <- capture_data_all %>%
  filter(cub_status == "c",
         year >= 1992,
         day_number < 150,
         merking %in% c("n", "r"),
         !is.na(z_length)) %>%
  mutate(year = as.numeric(year),
         age_for_analyses = as.numeric(age_for_analyses),
         cub_number = as.numeric(cub_number),
         cub_number_2 = ifelse(cub_number %in% c(2, 3), 4, as.numeric(cub_number))) %>%
  arrange(date)

write_csv(capture_data_females_w_cubs, "06_processed_data/CR_data/CR_f_with_cubs_clean.csv")


# 2. Females with cubs and no cubs (but not females with yearlings and 2-year-olds)
capture_data_females <- capture_data_all %>%
  filter(sex == "f",
         cub_status %in% c("c", "n"), # Keep
         year >= 1992,
         day_number < 150,
         merking %in% c("n", "r"),
         !is.na(z_length)) %>%
  mutate(year = as.numeric(year),
         age_for_analyses = as.numeric(age_for_analyses),
         cub_number = ifelse(is.na(cub_number), 0, as.numeric(cub_number)),
         cub_number_2 = ifelse(cub_number %in% c(2, 3), 4, as.numeric(cub_number))) %>%
  arrange(date)

write_csv(capture_data_females, "06_processed_data/CR_data/CR_f_clean.csv")

ggplot(capture_data_females, aes(x = as.factor(cub_number))) +
  geom_bar()


# B. Visualization =============================================================

capture_data_females_w_cubs <- read_csv("06_processed_data/CR_data/CR_f_with_cubs_clean.csv")

# 1. Litter size accross years
(
plot_1 <- capture_data_females_w_cubs %>%
  mutate(cub_number = cub_number) %>%
  group_by(year) %>%
  summarise(average = mean(cub_number)) %>%
  ggplot(aes(x = year, y = average)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_classic() +
  labs(x = "Year",
       y = "Mean litter size")
)

ggsave("7_results/1_interim_results/polar_bear_data_exploration/03-04 mean litter size VS years.png",
       plot = plot_1, width = 6, height = 6, units = "cm")

 
(plot_1_2 <- capture_data_females_w_cubs %>%
    mutate(CUB_NUMBER = as.numeric(CUB_NUMBER)) %>%
    ggplot(aes(x = as.numeric(YEAR), y = CUB_NUMBER)) +
    geom_jitter(width = 0.25,
                height = 0.15) +
    geom_smooth(method = lm) +
    theme_bw() +
    labs(x = "Year",
         y = "Litter size")
)

ggsave("7_results/1_interim_results/polar_bear_data_exploration/03-04 jitter litter size VS years.png",
       plot = plot_1_2, width = 6, height = 6, units = "cm")


(
plot_1_3 <- capture_data_females_w_cubs %>%
  ggplot(aes(x = YEAR, fill = as.factor(CUB_NUMBER))) +
  geom_bar() +
  scale_fill_grey(start = 0.8, end = 0.2) +
  theme_classic() +
  labs(x = "Year",
       y = "Count",
       fill = "Litter size")
)

ggsave("7_results/1_interim_results/polar_bear_data_exploration/03-04 barplot litter size VS years.png",
       plot = plot_1_3, width = 8, height = 6, units = "cm")




data.temp <- capture_data_females_w_cubs %>%
  mutate(data_point = 1) %>%
  group_by(YEAR) %>%
  summarise(nbr_captures = sum(data_point))

data.temp.2 <- capture_data_females_w_cubs %>%
  count(CUB_NUMBER, YEAR)

percentage <- rep(NA, times = nrow(data.temps.2))
for (k in 1:nrow(data.temps.2)) {
  percentage[k] <- data.temps.2$n[k]/data.temp$nbr_captures[data.temp$YEAR == data.temps.2$YEAR[k]]
}
data.temp.3 <- data.frame(data.temp.2,
                           proportion = percentage)




(plot_1_4 <- ggplot() +
  geom_col(data = data.temp.3, aes(x = YEAR, y = proportion, fill = as.factor(CUB_NUMBER))) +
  geom_text(data = data.temp, 
            aes(x = YEAR, y = 1.05, label = nbr_captures),
            size = 3) +
    scale_fill_grey(start = 0.8, end = 0.2) +
    theme_classic() +
    labs(x = "Year",
         y = "proportion",
         fill = "Litter size")
)


ggsave("7_results/1_interim_results/polar_bear_data_exploration/03-04 barplot proportion litter size VS years.png",
       plot = plot_1_4, width = 8, height = 6, units = "cm")

rm(data.temp,
   data.temp.2,
   data.temp.3,
   percentage)



# Date of capture
(plot_2 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = DAY_NUMBER)) +
    geom_histogram(bins = 24) +
    theme_classic() +
    labs(x = "Day"))


(plot_2_1 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = YEAR, y = DAY_NUMBER)) +
    geom_point() +
    geom_smooth(method = lm) +
    theme_classic() +
    labs(x = "Year",
         y = "Day capture")
  )



# Size 
(plot_3 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = S_LENGTH)) +
    geom_histogram(bins = 30) +
    theme_classic() +
    labs(x = "Length"))


(plot_3_1 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = DAY_NUMBER, y = S_LENGTH)) +
    geom_point() +
    geom_smooth(method = lm) +
    theme_classic() +
    labs(x = "Day",
         y = "Length")
)

(plot_3_2 <- capture_data_females_w_cubs %>%
    group_by(YEAR) %>%
    summarise(AVERAGE = mean(S_LENGTH)) %>%
    ggplot(aes(x = YEAR, y = AVERAGE)) +
    geom_point() +
    geom_smooth(method = lm) +
    theme_classic() +
    labs(x = "Year",
         y = "Mean length")
)


(plot_3_3 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = YEAR, y = S_LENGTH)) +
    geom_violin() +
    geom_jitter(width = 0.25) +
    geom_smooth(method = lm) +
    scale_x_discrete(breaks = c("1992", "2000", "2010", "2019")) +
    theme_classic() +
    labs(x = "Year",
         y = "Length")
)


# Age
(plot_4 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = AGE_FOR_ANALYSES)) +
    geom_histogram(bins = 20) +
    theme_classic() +
    labs(x = "Age"))


(plot_4_1 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = as.numeric(YEAR), y = AGE_FOR_ANALYSES)) +
    geom_point() +
    geom_smooth(method = lm) +
    theme_classic() +
    labs(x = "Year",
         y = "Age")
)



# Litter size VS on day capture
(plot_5 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = as.factor(CUB_NUMBER), y = DAY_NUMBER)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_jitter(width = 0.25) +
    theme_classic() +
    labs(x = "Litter size",
         y = "Day")
)

ggsave("7_results/1_interim_results/polar_bear_data_exploration/03-04 litter size VS day capture.png",
       plot = plot_5, width = 6, height = 6, units = "cm")


# Litter size VS mother's age
(plot_6 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = as.factor(CUB_NUMBER), y = AGE_FOR_ANALYSES)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_jitter(width = 0.25) +
    theme_classic() +
    labs(x = "Litter size",
         y = "Mother's age")
)

ggsave("7_results/1_interim_results/polar_bear_data_exploration/03-04 litter size VS mother's age.png",
       plot = plot_6, width = 6, height = 6, units = "cm")


# Litter size VS mother's size
(plot_7 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = as.factor(CUB_NUMBER), y = S_LENGTH)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_jitter(width = 0.25) +
    theme_classic() +
    labs(x = "Litter size",
         y = "Mother's length")
)

ggsave("7_results/1_interim_results/polar_bear_data_exploration/03-04 litter size VS mother's size.png",
       plot = plot_7, width = 6, height = 6, units = "cm")


# Litter size VS date of capture
(plot_8 <- capture_data_females_w_cubs %>%
    ggplot(aes(x = as.factor(CUB_NUMBER), y = DAY_NUMBER)) +
    geom_violin() +
    geom_boxplot(width=0.1) +
    geom_jitter(width = 0.25) +
    theme_classic() +
    labs(x = "Litter size",
         y = "Day of capture")
)




# Plot the position of the carcasses +++++++++++++++++++++
library(raster)
library(tmap)

# Create a shapefile with spatial points from the coordinates in the csv file
shapefile_captures <- SpatialPointsDataFrame(capture_data_females_w_cubs[, c(10, 9)], #the coordinates
                                       capture_data_females_w_cubs,    #the R object to convert
                                       proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Load map of Svalbard
svalbard <- rgdal::readOGR(dsn = "06_processed_data/svalbard_map",
                           layer = "svalbard_map")
crs(svalbard)


tmap_mode("view")
(map_captures <- tm_shape(svalbard) +
    tm_polygons(col  =  "grey90", alpha = 1, boundary.col = "grey50") +
    tm_shape(shapefile_captures) +
    tm_symbols(size = 0.15, col = "#E69F00" , border.lwd = 0.3, border.col = "black") +
    tm_grid(lines = FALSE,
            labels.size = 1)
  
)

tmap_save(tm = map_captures,
          filename = "07_results/01_interim_results/polar_bear_data_exploration/capture_map.png")


# With ggplot
library(broom)
svalabrd.df <- tidy(svalbard)


ggplot() +
  geom_polygon(data = svalabrd.df, aes(x = long, y = lat, group = group), 
               fill="grey90", color = "grey50") +
  geom_point(data = capture_data_females_w_cubs, 
             aes(x = long, y = lat, color = year),
             size = 2) +
  theme_bw() +
  labs(x = "Longitude",
       y = "Latidude",
       color = "Year")

ggsave("10_meetings/2021-04-09 Meeting with Sarah/map_captures.png",
       width = 12, height = 8)



# Jon's explanations in his email from 19-02-2021 ------------------------------

-Data from 1993 only (first year NPI started to capture families including small cubs)

-data from spring only (March, April, May)

-sort out females with cubs of the year (and also yearlings?),

Do not remember what we agreed, but seems smart to make a data set with both, even if you start analysing only cubs of the year



Important columns are:
  
  YEAR

SEX (f=female)

MERKING (n=first capture, r=recapture, think that is the two you would want)

ID NR this is the name of the bear

LAT and LON = positions of capture

kohort : ad = adult, coy = cub of the year, yrlg = one year old, 2yr=2 year old, subad = 3-4 year old

Age for analyses : this is the age you want to use (not the earlier called AGE FIELD)

Z length : body length

S length : body length measured as a straight line above the bear

Condition : how fat a bear is 1-5 1 being very lean, 5 being very fat

Girth : measure of circumference around waist



ABEAR CLASS1 SEX1 and so on : data on other individuals that were with a bear (e.g. the cubs)

Cub status : for adult females, c = with cubs of the year, y=with yearlings,  2y = with two year old, n = alone (or with other adult, but not with any cubs)

Cub number : the number of cubs she had





As early females have no weight measure, we usually estimate it for all,

And we also use to calculate an index for body condition (the one above is based on a subjective field score)

See attached excel sheet for calculations of weight and body condition, and citations to papers they are based on

(based on length and girth)