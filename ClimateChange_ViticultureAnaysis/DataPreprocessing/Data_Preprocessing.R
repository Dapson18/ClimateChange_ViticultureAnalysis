
library(dismo)
library(tidyverse)
library(sp)
library(caret)
library(pROC)


data <- read.table('C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\0018035-240506114902167\\occurrence.txt', sep="\t", header=TRUE, fill = TRUE)

dim(data)
colnames(data)

data <- data %>% filter(basisOfRecord == 'OCCURRENCE' | basisOfRecord == 'HUMAN_OBSERVATION')
table(data$basisOfRecord)

###################################

###  Do we need to georeference??
### Check for rows where locality has values, but lat and lon are empty strings or NA
# Filter data based on the conditions and select the necessary columns
filtered_data <- data %>%
  filter(locality != "" & (decimalLatitude == "" | is.na(decimalLongitude))) %>%
  select(locality, decimalLatitude, decimalLongitude)

#### Check the row in our dataframe that satisfy this condition and confirm via 
## inspection to confirm the result
filtered_row <-  which(data$locality != "" & (data$decimalLatitude == "" | is.na
                                              (data$decimalLongitude)))

### Check for the row number where locality has values, but lat and lon are empty strings
filtered_indices <- which(data$locality != "" & data$decimalLatitude == "" & is.na(data$decimalLongitude))

filt_lat <- which(data$locality != "" & data$decimalLatitude == "")

filt_long <- which(data$locality != "" & is.na(data$decimalLongitude))


#### Create a subset of the dataframe to contain only rows with lon and lat
#gp <- subset(data, !is.na(decimalLongitude) & !is.na(decimalLatitude))
gp <- subset(data, !is.na(decimalLongitude) & decimalLatitude != "")
### Show the rows where lat and lon are not present
lp <- which(is.na(data$decimalLongitude) & data$decimalLatitude == "")

dim(gp)
head(gp)


library(maptools)

library(maps)

### Create a map for visualization

data(wrld_simpl)
plot(wrld_simpl, xlim=c(-124,-69), ylim=c(18,49), axes=TRUE, col="light yellow")
# restore the box around the map
box()
# add the points
points(gp$decimalLongitude, gp$decimalLatitude, col='orange', pch=20, cex=0.75)
# plot points again to add a border, for better visibility
points(gp$decimalLongitude, gp$decimalLatitude, col='red', cex=0.75)






#################################
### Removing duplicates

gp_clean <- gp %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)


####################################################################################
#########################################


library(mice)

# Perform multiple imputation
imputed_data <- mice(gp_clean[, c("decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters")], m=5, method='pmm', maxit=50)

# Create a complete dataset from the imputed data
complete_data <- complete(imputed_data, 1)

# Now 'complete_data' can be used for further SDM analysis


############ Check that the imputation done did not introduce bias to the analysis 

library(ggplot2)
library(gridExtra)

# Plot the distribution of original values
p1 <- ggplot(gp_clean, aes(x = coordinateUncertaintyInMeters)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
  labs(title = "Original Data Distribution", xlab = 'Coordinate Uncertainty (metres)')

# Plot the distribution of imputed values
p2 <- ggplot(complete_data, aes(x = coordinateUncertaintyInMeters)) +
  geom_histogram(bins = 30, fill = "red", alpha = 0.5) +
  labs(title = "Imputed Data Distribution", xlab = 'Coordinate Uncertainty (metres)')


# Combine the two plots
grid.arrange(p1, p2, ncol = 2)

#### Summary statistics
# Summary of original data
summary(gp_clean$coordinateUncertaintyInMeters)

# Summary of imputed data
summary(complete_data$coordinateUncertaintyInMeters)

###############################################
#### Duplicating our dataset
new_gp <- gp_clean

### Replacing coordinate uncertainity in the original dataset
# with that acquired by multiple imputation
new_gp$coordinateUncertaintyInMeters <- complete_data$coordinateUncertaintyInMeters

### Filtering out unreliable data
gp_cleaner <- new_gp %>% filter(coordinateUncertaintyInMeters < 1000)
gp_cleaner$coordinateUncertaintyInMeters
#########################################
########## Ensuring Location is correctly represented


### best options for reverse geocoding is the tidygeocoder package
library(tidygeocoder)

# Perform reverse geocoding using the reverse_geocode function
results <- complete_data %>% reverse_geocode(lat = decimalLatitude,
                                             long = decimalLongitude,
                                             method = 'osm')


print(results)

t.test(gp_clean$coordinateUncertaintyInMeters, complete_data$coordinateUncertaintyInMeters, na.rm = TRUE)

#### The result from the above is better than the locality in the original dataset
########################################
### Add the result from the address into our dataset
newer_gp <- new_gp %>% mutate(address = results$address)

## Assign the correct country to missing data
newer_gp$level0Gid[327] <- 'USA'

# Remove rows in which country or level0Gid cannot be identified by geocoding
gp_countryCorrected<- newer_gp %>% filter(!row_number() %in% 338)



#############################################


data(wrld_simpl)
plot(wrld_simpl, xlim=c(-123,-64), ylim=c(18,49), axes=TRUE, col="light yellow")
# restore the box around the map
box()
# add the points
points(gp_countryCorrected$decimalLongitude, gp_countryCorrected$decimalLatitude, col='orange', pch=20, cex=0.75)
# plot points again to add a border, for better visibility
points(gp_countryCorrected$decimalLongitude, gp_countryCorrected$decimalLatitude, col='red', cex=0.75)



############ ENVIRONMENTAL DATA   ###########################


# Provide the full path to the folder containing the .tif files
biovars <- stack(list.files(path = "C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\wc2.1_2.5m_bio", pattern = "*.tif$", full.names = TRUE))

# Check the names of the layers
names(biovars)

# Convert character columns to numeric
gp_countryCorrected$decimalLongitude <- as.numeric(gp_countryCorrected$decimalLongitude)
gp_countryCorrected$decimalLatitude <- as.numeric(gp_countryCorrected$decimalLatitude)

# Now create the extent object
map_extent <- extent(min(gp_countryCorrected$decimalLongitude), 
                     max(gp_countryCorrected$decimalLongitude),
                     min(gp_countryCorrected$decimalLatitude), 
                     max(gp_countryCorrected$decimalLatitude))

# Crop the global raster to the North America extent
cropped_raster <- crop(biovars, map_extent)

### Show each of the variable on a plot
plot(cropped_raster)

### plot of a single layer in a RasterStack
data("wrld_simpl")
# first layer of the RasterStack
plot(cropped_raster, 1)
# note the "add=TRUE" argument with plot
plot(wrld_simpl, add=TRUE)
# with the points function, "add" is implicit
points(gp_countryCorrected$decimalLongitude, gp_countryCorrected$decimalLatitude, col='blue')




#### PRESENCE DATA ###########

#### Create a dataframe for lon and lat alone
presence_data <- gp_countryCorrected %>% select(decimalLongitude, decimalLatitude)

# Make a copy of the presence_data
coords_presence <- presence_data

# Convert the matrix to SpatialPoints
sp_presence_points <- SpatialPoints(coords_presence, proj4string=CRS("+proj=longlat +datum=WGS84"))

# Extract environmental data for presence points
presence_env <- raster::extract(cropped_raster, sp_presence_points)

# Combine presence points with environmental variables
presence <- cbind(presence_data, presence_env)

### add binary indicator for  presence
presence <- presence %>% mutate(species = 1)



###### To create Background Data ########

# Plot the first layer to visually inspect
plot(cropped_raster[[1]])

# Check for NA in the first bioclim variable as see if it equates to water
na_raster <- is.na(cropped_raster[[1]])
## Count for the number of NA
na_count <- sum(values(na_raster), na.rm = TRUE)

# Plot the NA raster
plot(na_raster, main = "NA Values in Bioclimatic Variable")

# Invert the na_raster to create a land mask
land_mask <- calc(na_raster, fun = function(x) { ifelse(!x, 1, NA) })

# Ensure that the longitude values are in the correct order
min_longitude <- min(presence_data$decimalLongitude, na.rm = TRUE)
max_longitude <- max(presence_data$decimalLongitude, na.rm = TRUE)
min_latitude <- min(presence_data$decimalLatitude, na.rm = TRUE)
max_latitude <- max(presence_data$decimalLatitude, na.rm = TRUE)

# Now create the extent with the correct values
study_extent <- extent(min_longitude, max_longitude, min_latitude, max_latitude)


set.seed(123)

# Generate background points using the land mask
background_points <- randomPoints(land_mask, n = nrow(presence_data), ext = study_extent)

# Since 'background_points' is a matrix with longitude in the first column (x) and latitude in the second column (y)
coords <- background_points
colnames(coords) <- c("decimalLongitude", "decimalLatitude")

# Convert the matrix to SpatialPoints
sp_background_points <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84"))

## Check if there is an overlap between presence point and abscence point
overlaps <- over(sp_background_points, sp_presence_points)

# Now you can use the extract function
background_env <- raster::extract(cropped_raster, sp_background_points)

### convert background_points to a dataframe 
background_points_df <- as.data.frame(background_points)
colnames(background_points_df) <- c("decimalLongitude", "decimalLatitude")

# Combine background points with environmental variables
background_data <- cbind(background_points_df, background_env)

### add binary indicator for abscence
background_data <- background_data %>% mutate(species = 0)