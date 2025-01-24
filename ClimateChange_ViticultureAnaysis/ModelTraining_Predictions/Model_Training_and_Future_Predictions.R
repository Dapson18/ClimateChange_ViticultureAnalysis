# Fit the model using glm

# Combine presence and background data
combined_data <- rbind(presence, background_data)

combined_data <- na.omit(combined_data)

combined_data$species <- factor(combined_data$species)


set.seed(11)
# Fit a GLM with a binomial family and logit link function
# Exclude the longitude and latitude columns from the predictors
glm_model <- glm(species ~ ., data=combined_data[, -c(1,2)], family=binomial(link="logit"))

# Summary of the model to see the effects of each variable
summary(glm_model)

## Using the bioclim model which typically require only presence points defines the range of climate conditions under which a species is known to occur based on the climate data from locations where the species is present. 

bclim <- bioclim(presence[,c('wc2.1_2.5m_bio_10', 'wc2.1_2.5m_bio_13', 'wc2.1_2.5m_bio_14', 
                             'wc2.1_2.5m_bio_15', 'wc2.1_2.5m_bio_18', 'wc2.1_2.5m_bio_19',
                             'wc2.1_2.5m_bio_2', 'wc2.1_2.5m_bio_3', 'wc2.1_2.5m_bio_5', 
                             'wc2.1_2.5m_bio_9')])
pairs(bclim)

#########################################

# Perform stepwise model selection based on AIC
stepwise_model <- step(glm_model, direction="both")
summary(stepwise_model)



# Assuming 'combined_data' is your dataset with species presence/absence as the response variable
# and the WorldClim bioclimatic variables as predictors

# Refitting the GLM model with the selected variables from the step wise model
final_glm_model <- glm(species ~ wc2.1_2.5m_bio_10 + wc2.1_2.5m_bio_12 + 
                         wc2.1_2.5m_bio_13 + wc2.1_2.5m_bio_14 + wc2.1_2.5m_bio_15 + 
                         wc2.1_2.5m_bio_17 + wc2.1_2.5m_bio_19 + wc2.1_2.5m_bio_2 + 
                         wc2.1_2.5m_bio_6 + wc2.1_2.5m_bio_7 + wc2.1_2.5m_bio_9, data=combined_data, 
                       family=binomial(link="logit"))

summary(final_glm_model)

BIC(final_glm_model,glm_model)


### Create a reduced version of the final glm model that removes collinear variables
mod <- glm(species ~ wc2.1_2.5m_bio_10 +
             wc2.1_2.5m_bio_13 + wc2.1_2.5m_bio_14 + wc2.1_2.5m_bio_15 
           + wc2.1_2.5m_bio_19 + wc2.1_2.5m_bio_2 + 
             wc2.1_2.5m_bio_6 + wc2.1_2.5m_bio_7 + wc2.1_2.5m_bio_9, data=combined_data, 
           family=binomial(link="logit"))

summary(mod)
BIC(mod, final_glm_model)
AIC(mod, final_glm_model)

###########################################################################
############# Remove NA in cropped_raster ##############################

# Now you can use the cropped raster for predictions and visualizations
predictors_glm <- predict(cropped_raster, final_glm_model, type = 'response')
plot(predictors_glm)

predictors_mod <- predict(cropped_raster, mod, type = 'response')

plot(predictors_mod)

############### COMPARING PREDICTIVE POWER  ######################
set.seed(123) # Ensure reproducibility
splitIndex <- createDataPartition(combined_data$species, p = 0.8, list = FALSE)
glm_train_data <- combined_data[splitIndex, ]
glm_test_data <- combined_data[-splitIndex, ]

# Full model with all bioclimatic variables
full_model <- glm(species ~ ., family=binomial(link="logit"), data=glm_train_data[, -c(1,2)])

# Stepwise reduced model
step_model <- glm(species ~ wc2.1_2.5m_bio_10 + wc2.1_2.5m_bio_12 + 
                    wc2.1_2.5m_bio_13 + wc2.1_2.5m_bio_14 + wc2.1_2.5m_bio_15 + 
                    wc2.1_2.5m_bio_17 + wc2.1_2.5m_bio_19 + wc2.1_2.5m_bio_2 + 
                    wc2.1_2.5m_bio_6 + wc2.1_2.5m_bio_7 + wc2.1_2.5m_bio_9, 
                  family=binomial(link="logit"), data=glm_train_data)

mod <- glm(species ~ wc2.1_2.5m_bio_10 +
             wc2.1_2.5m_bio_13 + wc2.1_2.5m_bio_14 + wc2.1_2.5m_bio_15 
           + wc2.1_2.5m_bio_19 + wc2.1_2.5m_bio_2 + 
             wc2.1_2.5m_bio_6 + wc2.1_2.5m_bio_7 + wc2.1_2.5m_bio_9, data=glm_train_data, 
           family=binomial(link="logit"))

# Predicted probabilities for the full model
predicted_probs_full <- predict(full_model, newdata=glm_test_data, type="response")

# Predicted probabilities for the stepwise reduced model
predicted_probs_step <- predict(step_model, newdata=glm_test_data, type="response")

# Predicted probabilities for the stepwise reduced model
predicted_probs_mod <- predict(mod, newdata=glm_test_data, type="response")

# Confusion matrix for the full model
cm_full <- confusionMatrix(as.factor(round(predicted_probs_full)), glm_test_data$species, positive = '1')

# Confusion matrix for the stepwise reduced model
cm_step <- confusionMatrix(as.factor(round(predicted_probs_step)), glm_test_data$species, positive = '1')

# Confusion matrix for the mod reduced model
cm_mod <- confusionMatrix(as.factor(round(predicted_probs_mod)), glm_test_data$species, positive = '1')

# Calculate AUC for both models
roc_full <- roc(glm_test_data$species, predicted_probs_full)
roc_step <- roc(glm_test_data$species, predicted_probs_step)
roc_mod <- roc(glm_test_data$species, predicted_probs_mod)
auc_full <- auc(roc_full)
auc_step <- auc(roc_step)
auc_mod <- auc(roc_mod)

# Print AUC values
cat("AUC for Full Model:", auc_full, "\n")
cat("AUC for Stepwise Reduced Model:", auc_step, "\n")
cat("AUC for Mod Model:", auc_mod, "\n")
##### The result indicates that the mod model is best but because of ecologicial importance and uncertainty about future climate pattern

# Extract precision and recall for each model
precision_step <- cm_step$byClass['Precision']
recall_step <- cm_step$byClass['Recall']
f1_score_step <- 2 * (precision_step * recall_step) / (precision_step + recall_step)

precision_mod <- cm_mod$byClass['Precision']
recall_mod <- cm_mod$byClass['Recall']
f1_score_mod <- 2 * (precision_mod * recall_mod) / (precision_mod + recall_mod)

# Print the F1 scores and precision for each model
cat("F1 Score reduced model:", f1_score_step, "Precision:", precision_step, "\n")
cat("F1 Score for downsized model:", f1_score_mod, "Precision:", precision_mod, "\n")


AIC(step_model, mod, reduced_gam)
BIC(step_model, mod, reduced_gam)
###############################
summary(reduced_gam)

#### Using the vip package to create variable importance plots
library(vip)
vip(final_glm_model)

##### Comparing the predicted probabilities to the observed presence/absence.
combined <- combined_data
combined$predicted <- predict(gams_reduced, type="response")
ggplot(combined, aes(x=predicted, fill=factor(species))) +
  geom_histogram(position="identity", alpha=0.5, bins=30) +
  labs(x="Predicted Probability", y="Count", fill="Species Presence")


#### Generating a Receiver Operating Characteristic (ROC) curve to evaluate model performance
roc_curve <- roc(combined$species, combined$predicted)
plot(roc_curve, main="ROC Curve", col="#1c61b6")
abline(a=0, b=1, col="red")




############################################################
#############################################################


############## USING GAM ################################

###   Load the required library
library(mgcv)

set.seed(11)

# Fitting a GAM with all variables as smooth terms (excluding longitude and latitude)
gam_model <- gam(species ~ s(wc2.1_2.5m_bio_1) + s(wc2.1_2.5m_bio_2) + s(wc2.1_2.5m_bio_3) + 
                   s(wc2.1_2.5m_bio_4) + s(wc2.1_2.5m_bio_5) + s(wc2.1_2.5m_bio_6) + 
                   s(wc2.1_2.5m_bio_7) + s(wc2.1_2.5m_bio_8) + s(wc2.1_2.5m_bio_9) + 
                   s(wc2.1_2.5m_bio_10) + s(wc2.1_2.5m_bio_11) + s(wc2.1_2.5m_bio_12) + 
                   s(wc2.1_2.5m_bio_13) + s(wc2.1_2.5m_bio_14) + s(wc2.1_2.5m_bio_15) + 
                   s(wc2.1_2.5m_bio_16) + s(wc2.1_2.5m_bio_17) + s(wc2.1_2.5m_bio_18) + 
                   s(wc2.1_2.5m_bio_19), data=combined_data, 
                 family=binomial(link="logit"), method="REML", select=TRUE)


summary(gam_model)



#### Model containing the selected variables from the above 
gams_reduced <- gam(species ~ s(wc2.1_2.5m_bio_1) + s(wc2.1_2.5m_bio_2) + s(wc2.1_2.5m_bio_4) + 
                      s(wc2.1_2.5m_bio_5) + s(wc2.1_2.5m_bio_10) + s(wc2.1_2.5m_bio_11) + 
                      s(wc2.1_2.5m_bio_12) + s(wc2.1_2.5m_bio_13) + s(wc2.1_2.5m_bio_15) + 
                      s(wc2.1_2.5m_bio_19), 
                    data=combined_data[, -c(1,2)], family=binomial(link="logit"))

summary(gams_reduced)


##### After removing collinear variables
downsized_gam_model <- gam(species ~ s(wc2.1_2.5m_bio_1) + s(wc2.1_2.5m_bio_2) + s(wc2.1_2.5m_bio_5) 
                           + s(wc2.1_2.5m_bio_12) + s(wc2.1_2.5m_bio_15) + 
                             s(wc2.1_2.5m_bio_19), 
                           family=binomial(link="logit"), data=combined_data)
summary(downsized_gam_model)


### Check AIC and BIC
AIC(reduced_gam, downsized_gam_model, gam_model)
BIC(reduced_gam, downsized_gam_model, gam_model)

# Split the combined_data dataset into training and testing sets
set.seed(123) # for reproducibility
training_index <- createDataPartition(combined_data$species, p = 0.8, list = FALSE)
train_data <- combined_data[training_index, ]
test_data <- combined_data[-training_index, ]

# Fit your models on the training data
gam_full <- gam(species ~ s(wc2.1_2.5m_bio_1) + s(wc2.1_2.5m_bio_2) + s(wc2.1_2.5m_bio_3) + 
                  s(wc2.1_2.5m_bio_4) + s(wc2.1_2.5m_bio_5) + s(wc2.1_2.5m_bio_6) + 
                  s(wc2.1_2.5m_bio_7) + s(wc2.1_2.5m_bio_8) + s(wc2.1_2.5m_bio_9) + 
                  s(wc2.1_2.5m_bio_10) + s(wc2.1_2.5m_bio_11) + s(wc2.1_2.5m_bio_12) + 
                  s(wc2.1_2.5m_bio_13) + s(wc2.1_2.5m_bio_14) + s(wc2.1_2.5m_bio_15) + 
                  s(wc2.1_2.5m_bio_16) + s(wc2.1_2.5m_bio_17) + s(wc2.1_2.5m_bio_18) + 
                  s(wc2.1_2.5m_bio_19), data=train_data, 
                family=binomial(link="logit"))


reduced_gam <- gam(species ~ s(wc2.1_2.5m_bio_1) + s(wc2.1_2.5m_bio_2) + s(wc2.1_2.5m_bio_4) + 
                     s(wc2.1_2.5m_bio_5) + s(wc2.1_2.5m_bio_10) + s(wc2.1_2.5m_bio_11) + 
                     s(wc2.1_2.5m_bio_12) + s(wc2.1_2.5m_bio_13) + s(wc2.1_2.5m_bio_15) + 
                     s(wc2.1_2.5m_bio_19), 
                   family=binomial(link="logit"), data=train_data)

downsized_gam_model <- gam(species ~ s(wc2.1_2.5m_bio_1) + s(wc2.1_2.5m_bio_2) + 
                             s(wc2.1_2.5m_bio_5) +  
                             s(wc2.1_2.5m_bio_12) + s(wc2.1_2.5m_bio_15) + 
                             s(wc2.1_2.5m_bio_19), 
                           family=binomial(link="logit"), data=train_data)

# Calculate predicted probabilities for the test set
predicted_probs_full <- predict(gam_full, newdata=test_data, type="response")
predicted_probs_reduced <- predict(reduced_gam, newdata=test_data, type="response")
predicted_probs_downsized <- predict(downsized_gam_model, newdata=test_data, type="response")

# Calculate ROC curves
roc_curve_full <- roc(test_data$species, predicted_probs_full)
roc_curve_reduced <- roc(test_data$species, predicted_probs_reduced)
roc_curve_downsized <- roc(test_data$species, predicted_probs_downsized)

# Calculate AUC for each model
auc_full <- auc(roc_curve_full)
auc_reduced <- auc(roc_curve_reduced)
auc_downsized <- auc(roc_curve_downsized)

# Print AUC values
cat("AUC full model", auc_full, "\n")
cat("AUC reduced model", auc_reduced, "\n")
cat("AUC for downsized model:", auc_downsized, "\n")


# Calculate confusion matrices for each model
# Calculate the confusion matrix, specifying the positive class
cm_full <- confusionMatrix(as.factor(round(predicted_probs_full)), test_data$species, positive = "1")
cm_reduced <- confusionMatrix(as.factor(round(predicted_probs_reduced)), test_data$species, positive = "1")
cm_downsized <- confusionMatrix(as.factor(round(predicted_probs_downsized)), test_data$species, positive = "1")

# Now the 'Positive' class should be correctly considered as '1'


# Extract precision and recall for each model
precision_full <- cm_full$byClass['Precision']
recall_full <- cm_full$byClass['Recall']
f1_score_full <- 2 * (precision_full * recall_full) / (precision_full + recall_full)

precision_reduced <- cm_reduced$byClass['Precision']
recall_reduced <- cm_reduced$byClass['Recall']
f1_score_reduced <- 2 * (precision_reduced * recall_reduced) / (precision_reduced + recall_reduced)

precision_downsized <- cm_downsized$byClass['Precision']
recall_downsized <- cm_downsized$byClass['Recall']
f1_score_downsized <- 2 * (precision_downsized * recall_downsized) / (precision_downsized + recall_downsized)

# Print the F1 scores and precision for each model
cat("F1 Score full model:", f1_score_full, "Precision:", precision_full, "\n")
cat("F1 Score reduced model:", f1_score_reduced, "Precision:", precision_reduced, "\n")
cat("F1 Score for downsized model:", f1_score_downsized, "Precision:", precision_downsized, "\n")

#################  The GAM model has a better precision in prediction and a higher AUC value than the GLM method.

########### Visualizing the result from the model
# Now we can use the cropped raster for predictions and visualizations
predictors_df <- predict(cropped_raster, gams_reduced, type = 'response')
plot(predictors_df)


# Define breaks and corresponding colors
breaks1 <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
colors1 <- c("snow2", "orange1", "yellow2", "green3", 'green4')

# Plot the ensemble raster
plot(predictors_df, main="Current Distribution of Vitis", 
     breaks=breaks1, col=colors1, axis.args=list(at=breaks1, labels=breaks1))


# Add a legend to the plot
legend("bottomright",
       legend=c("Very Low", "Low", "Moderate", "High", 'Very High'),
       fill=colors_u4, title="Variability")



################ Variable IMportance plot for GAM


# Create a data frame with variable names and their importance measures
variable_importance <- data.frame(
  Variable = c("Bio_1", "Bio_2", "Bio_4", 
               "Bio_5", "Bio_10", "Bio_11", 
               "Bio_12", "Bio_13", "Bio_15", 
               "Bio_19"),
  Chi_sq = c(9.146, 12.945, 2.221, 24.410, 5.277, 1.811, 15.729, 1.087, 12.254, 25.645),
  p_value = c(0.069225, 0.000321, 0.136184, 0.003046, 0.192996, 0.178440, 0.026548, 0.585399, 0.041122, 0.002355)
)

# Order the data frame by Chi-squared values
variable_importance <- variable_importance[order(-variable_importance$Chi_sq),]

# Create the bar chart
ggplot(variable_importance, aes(x=reorder(Variable, Chi_sq), y=Chi_sq)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +  # Flip the coordinates to make the chart horizontal
  labs(title="Variable Importance in Predicting Species Presence",
       x="Chi-squared Value", y="Environmental Variable") +
  theme_minimal()









#### Create a variable that contains only important environmental variable required from combined_data
selected_combined_data <- combined_data %>%  select(c(3, 4, 5, 6, 7, 9, 13, 14, 16, 17))
names(selected_combined_data)

##### FUTURE ######################

# Create a RasterStack from the multi-band raster file
bioclim_stack <- stack("C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\Dissertation/wc2.1_2.5m_bioc_ACCESS-CM2_ssp126_2081-2100.tif")

# Check the RasterStack to confirm it contains all 19 bioclimatic variables
print(bioclim_stack)
plot(bioclim_stack)

# Crop the global raster to the study extent
future_raster <- crop(bioclim_stack, map_extent)
plot(future_raster)

names(future_raster)



#################  ENSEMBLE #####################


##### For ssp126
# Define the directory where your raster files are located
raster_dir <- "C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\Dissertation\\ssp126"

# List all files that end with '2081-2100.tif'
raster_files <- list.files(path = raster_dir, pattern = "ssp126_2081-2100.tif$", full.names = TRUE)

# Initialize a list to store prediction rasters
prediction_rasters <- list()

# Loop over each file
for (raster_file in raster_files) {
  # Load the raster file
  bioclim_stack <- stack(raster_file)
  
  # Crop the global raster to the study extent
  future_raster <- crop(bioclim_stack, map_extent)
  
  # Create a data frame for all raster cells
  all_cells_df <- as.data.frame(future_raster, xy=TRUE)
  
  # Select the necessary columns for prediction
  selected_future <- c(3, 4, 6, 7, 12, 13, 14, 15, 17, 21)
  future_data_prediction <- all_cells_df[, selected_future]
  names(future_data_prediction) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_4", 
                                     "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", 
                                     "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", 
                                     "wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_19")
  ##### Remove NA  ################
  #future_data_prediction <- na.omit(future_data_prediction)
  
  # Predict for the full spatial extent
  all_predictions <- predict(gams_reduced, newdata=future_data_prediction, type = 'response')
  
  # Convert the predictions to a RasterLayer
  full_prediction_raster <- raster(future_raster)
  
  # Convert all_predictions to a simple numeric vector
  all_predictions_vector <- as.vector(all_predictions)
  
  # Assign the predictions to the raster
  values(full_prediction_raster) <- all_predictions_vector
  
  # Store the prediction raster in the list
  prediction_rasters[[basename(raster_file)]] <- full_prediction_raster
  
  print(names(prediction_rasters))
}

# After the loop, create an ensemble raster by calculating the mean of all prediction rasters
ensemble_raster <- stack(prediction_rasters) |> calc(fun = mean, na.rm = TRUE)

# Define breaks and corresponding colors
breaks1 <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
colors1 <- c("snow2", "orange1", "yellow2", "green3", 'green4')

# Plot the ensemble raster
plot(ensemble_raster, main="Ensemble Predicted Species Distribution Presence (SSP126)", 
     breaks=breaks1, col=colors1, axis.args=list(at=breaks1, labels=breaks1))

# Add a legend to the plot
legend("bottomright", 
       legend=c("Very Low", "Low", "Moderate", "High", 'Very High'),
       fill=colors1, 
       title="Probability")


################# For 245 #################################
# Define the directory where your raster files are located
raster_dir_2 <- "C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\Dissertation\\ssp245"

# List all files that end with '2081-2100.tif'
raster_files2 <- list.files(path = raster_dir_2, pattern = "ssp245_2081-2100.tif$", full.names = TRUE)

# Initialize a list to store prediction rasters
prediction_rasters_2 <- list()

# Loop over each file
for (raster_file in raster_files2) {
  # Load the raster file
  bioclim_stack2 <- stack(raster_file)
  
  # Crop the global raster to the study extent
  future_raster2 <- crop(bioclim_stack2, map_extent)
  
  # Create a data frame for all raster cells
  all_cells_df_2 <- as.data.frame(future_raster2, xy=TRUE)
  
  # Select the necessary columns for prediction
  selected_future <- c(3, 4, 6, 7, 12, 13, 14, 15, 17, 21)
  future_data_prediction2 <- all_cells_df_2[, selected_future]
  names(future_data_prediction2) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_4", 
                                      "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", 
                                      "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", 
                                      "wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_19")
  
  ##### Remove NA  ################
  #future_data_prediction2 <- na.omit(future_data_prediction2)
  
  # Predict for the full spatial extent
  all_predictions2 <- predict(gams_reduced, newdata=future_data_prediction2, type = 'response')
  
  # Convert the predictions to a RasterLayer
  full_prediction_raster2 <- raster(future_raster2)
  
  # Convert all_predictions to a simple numeric vector
  all_predictions_vector2 <- as.vector(all_predictions2)
  
  # Assign the predictions to the raster
  values(full_prediction_raster2) <- all_predictions_vector2
  
  # Store the prediction raster in the list
  prediction_rasters_2[[basename(raster_file)]] <- full_prediction_raster2
  
  #print(names(prediction_rasters_2))
}

# After the loop, create an ensemble raster by calculating the mean of all prediction rasters
ensemble_raster2 <- stack(prediction_rasters_2) |> calc(fun = mean, na.rm = TRUE)

# Define breaks and corresponding colors
breaks2 <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
colors2 <- c("snow2", "orange1", "yellow2", "green3", 'green4')

# Plot the ensemble raster
plot(ensemble_raster2, main="Ensemble Predicted Species Distribution Presence (SSP245)", 
     breaks=breaks2, col=colors2, axis.args=list(at=breaks2, labels=breaks2))

# Add a legend to the plot
legend("bottomright", 
       legend=c("Very Low", "Low", "Moderate", "High"),
       fill=colors2, 
       title="Probability")

# Plot the ensemble raster
plot(ensemble_raster2, main="Ensemble Predicted Species Distribution foe SSP245")





##### For ssp370

# Define the directory where your raster files are located
raster_dir_3 <- "C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\Dissertation\\ssp370"

# List all files that end with '2081-2100.tif'
raster_files3 <- list.files(path = raster_dir_3, pattern = "ssp370_2081-2100.tif$", full.names = TRUE)

# Initialize a list to store prediction rasters
prediction_rasters_3 <- list()

# Loop over each file
for (raster_file in raster_files3) {
  # Load the raster file
  bioclim_stack3 <- stack(raster_file)
  
  # Crop the global raster to the study extent
  future_raster3 <- crop(bioclim_stack3, map_extent)
  
  # Create a data frame for all raster cells
  all_cells_df_3 <- as.data.frame(future_raster3, xy=TRUE)
  
  # Select the necessary columns for prediction
  selected_future <- c(3, 4, 6, 7, 12, 13, 14, 15, 17, 21)
  future_data_prediction3 <- all_cells_df_3[, selected_future]
  names(future_data_prediction3) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_4", 
                                      "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", 
                                      "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", 
                                      "wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_19")
  ##### Remove NA  ################
  #future_data_prediction3 <- na.omit(future_data_prediction3)
  
  # Predict for the full spatial extent
  all_predictions3 <- predict(gams_reduced, newdata=future_data_prediction3, type = 'response')
  
  # Convert the predictions to a RasterLayer
  full_prediction_raster3 <- raster(future_raster3)
  
  # Convert all_predictions to a simple numeric vector
  all_predictions_vector3 <- as.vector(all_predictions3)
  
  # Assign the predictions to the raster
  values(full_prediction_raster3) <- all_predictions_vector3
  
  # Store the prediction raster in the list
  prediction_rasters_3[[basename(raster_file)]] <- full_prediction_raster3
  
  #print(names(prediction_rasters_2))
}

# After the loop, create an ensemble raster by calculating the mean of all prediction rasters
ensemble_raster3 <- stack(prediction_rasters_3) |> calc(fun = mean, na.rm = TRUE)

# Define breaks and corresponding colors
breaks3 <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
colors3 <- c("snow2", "orange1", "yellow2", "green2", 'green4')

# Plot the ensemble raster
plot(ensemble_raster3, main="Ensemble Predicted Species Distribution Presence (SSP370)", 
     breaks=breaks3, col=colors3, axis.args=list(at=breaks3, labels=breaks3))

# Add a legend to the plot
legend("bottomleft", 
       legend=c("Very Low", "Low", "Moderate", "High", "Very High"),
       fill=colors3, 
       title="Probability")

# Plot the ensemble raster
plot(ensemble_raster3, main="Ensemble Predicted Species Distribution FOR SSP370")




####  FOR SSP 585
# Define the directory where your raster files are located
raster_dir_4 <- "C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\Dissertation\\ssp585"

# List all files that end with '2081-2100.tif'
raster_files4 <- list.files(path = raster_dir_4, pattern = "ssp585_2081-2100.tif$", full.names = TRUE)

# Initialize a list to store prediction rasters
prediction_rasters_4 <- list()

# Loop over each file
for (raster_file in raster_files4) {
  # Load the raster file
  bioclim_stack4 <- stack(raster_file)
  
  # Crop the global raster to the study extent
  future_raster4 <- crop(bioclim_stack4, map_extent)
  
  # Create a data frame for all raster cells
  all_cells_df_4 <- as.data.frame(future_raster4, xy=TRUE)
  
  # Select the necessary columns for prediction
  selected_future <- c(3, 4, 6, 7, 12, 13, 14, 15, 17, 21)
  future_data_prediction4 <- all_cells_df_4[, selected_future]
  names(future_data_prediction4) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_4", 
                                      "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", 
                                      "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", 
                                      "wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_19")
  ##### Remove NA  ################
  #future_data_prediction4 <- na.omit(future_data_prediction4)
  
  # Predict for the full spatial extent
  all_predictions4 <- predict(gams_reduced, newdata=future_data_prediction4, type = 'response')
  
  # Convert the predictions to a RasterLayer
  full_prediction_raster4 <- raster(future_raster4)
  
  # Convert all_predictions to a simple numeric vector
  all_predictions_vector4 <- as.vector(all_predictions4)
  
  # Assign the predictions to the raster
  values(full_prediction_raster4) <- all_predictions_vector4
  
  # Store the prediction raster in the list
  prediction_rasters_4[[basename(raster_file)]] <- full_prediction_raster4
  
  #print(names(prediction_rasters_2))
}

# After the loop, create an ensemble raster by calculating the mean of all prediction rasters
ensemble_raster4 <- stack(prediction_rasters_4) |> calc(fun = mean, na.rm = TRUE)

# Plot the ensemble raster
plot(ensemble_raster4, main="Ensemble Predicted Species Distribution for SSP585")


# Define breaks and corresponding colors
breaks4 <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
colors4 <- c("snow2", "orange1", "yellow2", "green2", 'green4')

# Plot the ensemble raster
plot(ensemble_raster4, main="Ensemble Predicted Species Distribution Presence (SSP585)", 
     breaks=breaks4, col=colors4, axis.args=list(at=breaks4, labels=breaks4))

# Add a legend to the plot
par(xpd = TRUE)
legend("bottomright", inset = c(-0.2, 0),
       legend=c("Very Low", "Low", "Moderate", "High", "Very High"),
       fill=colors4, 
       title="Probability")

###################################################################
#########################################################

###########################################################################
#################################################################
#### Another way to make the plot

# Define the number of colors you want to use
number_of_colors <- 5

# Generate a color gradient using terrain.colors
color_gradient <- rev(terrain.colors(number_of_colors))

# Plot the raster with the color gradient
plot(ensemble_raster4, 
     main="Ensemble Predicted Species Distribution for SSP126",
     col=color_gradient,
     breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
     axis.args=list(at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2)))

# Add a legend to the plot
legend("bottomright", 
       legend=c("Very Low", "Low", "Moderate", "High", "Very High"),
       fill=color_gradient, 
       title="Probability",
       border=NA)

############################################################
#################################################


par(mfrow = c(1, 1))
par(mar = c(5, 4, 4, 2) + 0.1) ## Default


#########################################################
#### CHECK IF THERE ARE ZERO IN THE PREDICTED VALUES

zero_values_indices <- which(values(predictors_df) == 0)

#####################################################################
#### Check if there is any cell that has same value in both raster and can lead to 0###

difference_raster_c <- overlay(predictors_df, ensemble_raster, fun = function(x, y) {
  # Check if the values are the same
  ifelse(x == y, 1, NA)  # Returns 1 for no change, NA for change
})

unchanged_indices <- which(values(difference_raster_c) == 1)



############### FOR 126 #############
# Calculate the difference
difference_raster <- overlay(predictors_df, ensemble_raster, fun = function(x, y) {
  return(x - y)
})

# Plot the difference
plot(difference_raster, main="Diff between Presence/Absence and Ensemble Prediction SSP126", )


summary(difference_raster)

################  FOR 126  #################

########################################
########### BREAKS 126 ############

####################################
# Define the number of colors you want to use
number_of_colors_e <- 4

# Generate a color gradient using terrain.colors
color_gradient_e <- terrain.colors(number_of_colors_e)

# Plot the raster with the color gradient
plot(difference_raster, 
     main="Diff between Presence/Absence and Ensemble Prediction SSP126",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))

# Add a legend to the plot
legend("bottomright", 
       legend=c("Significant Increase", "Increase", "Decrease", "Significant Decrease"),
       fill=color_gradient_e, 
       title="Probable Presence",
       border=NA)



##########################3


###############################################################


################  FOR 245  #################
# Calculate the difference
difference_raster2 <- overlay(predictors_df, ensemble_raster2, fun = function(x, y) {
  return(x - y)
})

# Plot the difference
plot(difference_raster2, main="Difference between Presence/Absence and Ensemble2 Prediction")


# Define the number of colors you want to use
number_of_colors_e <- 4

# Generate a color gradient using terrain.colors
color_gradient_e <- terrain.colors(number_of_colors_e)

# Plot the raster with the color gradient
plot(difference_raster2, 
     main="Diff between Presence/Absence and Ensemble Prediction SSP126",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))

# Add a legend to the plot
legend("bottomright", 
       legend=c("Significant Increase", "Increase", "Decrease", "Significant Decrease"),
       fill=color_gradient_e, 
       title="Probable Presence",
       border=NA)

########################################


############# FOR 370 #################
# Calculate the difference
difference_raster3 <- overlay(predictors_df, ensemble_raster3, fun = function(x, y) {
  return(x - y)
})

# Plot the difference
plot(difference_raster3, main="Difference between Presence/Absence and Ensemble3 Prediction")


# Define the number of colors you want to use
number_of_colors_e <- 4

# Generate a color gradient using terrain.colors
color_gradient_e <- terrain.colors(number_of_colors_e)

# Plot the raster with the color gradient
plot(difference_raster3, 
     main="Diff between Presence/Absence and Ensemble Prediction SSP370",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))

# Add a legend to the plot
legend("bottomright", 
       legend=c("Significant Increase", "Increase", "Decrease", "Significant Decrease"),
       fill=color_gradient_e, 
       title="Probable Presence",
       border=NA)

########################################
#
########### FOR 585 ###############
# Calculate the difference
difference_raster4 <- overlay(predictors_df, ensemble_raster4, fun = function(x, y) {
  return(x - y)
})

# Plot the difference
plot(difference_raster4, main="Difference between Presence/Absence and Ensemble4 Prediction")


# Define the number of colors you want to use
number_of_colors_e <- 4

# Generate a color gradient using terrain.colors
color_gradient_e <- terrain.colors(number_of_colors_e)

# Plot the raster with the color gradient
plot(difference_raster4, 
     main="Diff between Presence/Absence and Ensemble Prediction SSP585",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))

# Add a legend to the plot
legend("bottomright", 
       legend=c("Significant Increase", "Increase", "Decrease", "Significant Decrease"),
       fill=color_gradient_e, 
       title="Probable Presence",
       border=NA)



########################################
########### BREAKS 585 ############

# Define breaks for the difference values, ensuring a range around zero
breaks585 <- c(-1, -0.5, -0.05, 0.05, 0.5, 1)

# Define colors for each interval, including a range around zero for 'white'
colors585 <- c("blue", "lightblue", "white", "pink", "red")
plot(difference_raster4, main="Diff between Presence/Absence and Ensemble Prediction SSP585",
     breaks=breaks585, col=colors585, legend=TRUE)

# Add a custom legend to describe the intervals
legend("bottomright", legend=c("Increase (-1 to -0.5)", "Slight increase (-0.5 to -0.05)", 
                               "Neglible/No change (-0.05 to 0.05)", "Slight decrease (0.05 to 0.5)", 
                               "decrease (0.5 to 1)"), fill=colors585, bty="n")






############ Roc curves

roc_ssp126 <- roc(current_bin, bin1)
roc_ssp245 <- roc(current_bin, bin2)
roc_ssp370 <- roc(current_bin, bin3)
roc_ssp585 <- roc(current_bin, bin4)


# Plotting the ROC curves for each SSP scenario
plot(roc_ssp126, main="ROC Curves for Different SSP Scenarios", col="blue", lwd=2)
plot(roc_ssp245, add=TRUE, col="red", lwd=2)
plot(roc_ssp370, add=TRUE, col="green", lwd=2)
plot(roc_ssp585, add=TRUE, col="purple", lwd=2)
legend("bottomright", legend=c("SSP126", "SSP245", "SSP370", "SSP585"),
       col=c("blue", "red", "green", "purple"), lwd=2)





############################################################
############### SUITABILITY ACROSS ALL SSP ##################
###########################################################
par(mfrow = c(1, 1))

# Define the suitability threshold
threshold <- 0.48

# Identify suitable areas across all SSPs
suitable_areas_all_ssp <- overlay(ensemble_raster, ensemble_raster2, ensemble_raster3, ensemble_raster4, fun = function(x, y, z, w) {
  return((x > threshold) & (y > threshold) & (z > threshold) & (w > threshold))
})

# Plot the suitable areas
plot(suitable_areas_all_ssp, main="Areas suitable across all SSPs")



# Define breaks for the difference values, ensuring a range around zero
breaks_s1 <- c(0, 0.48, 1)

# Define colors for each interval, including a range around zero for 'white'
colors_s1 <- c("orange", "blue")
plot(suitable_areas_all_ssp, main="Areas suitable across all SSPs", xlab = 'Longitude', ylab = 'Latitude',
     breaks=breaks_s1, col=colors_s1, legend=TRUE)

# Add a custom legend to describe the intervals
legend("bottomleft", legend=c("Absence", "Presence"), fill=colors_s1, bty="n")



# Define breaks for the difference values, ensuring a range around zero
breaks_p <- c(0, 0.48, 1)

# Define colors for each interval, including a range around zero for 'white'
colors_p <- c("orange", "blue")
plot(binary_raster, main="Current Suitability of Vitis Vinifera", xlab = 'Longitude', ylab = 'Latitude',
     breaks=breaks_p, col=colors_p, legend=TRUE)

# Add a custom legend to describe the intervals
legend("bottomleft", legend=c("Absence", "Presence"), fill=colors_p, bty="n")




###############################################################################
#################### HISTOGRAN OF SUITABILITY SCORE #################
######################### CONTINUOUS DATA ###########################
# Assuming 'all_predictions' contains the model's prediction scores
hist(all_pv, main="Histogram of Suitability Scores SSP126",
     xlab="Suitability Score", ylab="Frequency", col="steelblue", border="black")


# Assuming 'all_predictions' contains the model's prediction scores
hist(all_pv2, main="Histogram of Suitability Scores SSP245",
     xlab="Suitability Score", ylab="Frequency", col="steelblue", border="black")

# Assuming 'all_predictions' contains the model's prediction scores
hist(all_pv3, main="Histogram of Suitability Scores SSP 370",
     xlab="Suitability Score", ylab="Frequency", col="steelblue", border="black")

# Assuming 'all_predictions' contains the model's prediction scores
hist(all_pv4, main="Histogram of Suitability Scores SSP585",
     xlab="Suitability Score", ylab="Frequency", col="steelblue", border="black")



##################################################
########### Determining Threshold ##################
############################################

library(pROC)

# Calculate the ROC curve
roc_obj <- roc(test_data$species, predicted_probs_reduced)

# Find the coordinates of the optimal point
coords12 <- coords(roc_obj, "best", ret="threshold", best.method="youden")

# Optimal threshold based on Youden's Index
optimal_threshold <- coords12$threshold

# Apply the optimal threshold to your predictions
predicted_classes <- ifelse(predicted_probs_reduced > optimal_threshold, 1, 0)

# Now you can create your confusion matrix or any other performance metric with this new threshold



#####Check coLLINEARITY


install.packages("corrplot")
install.packages("Hmisc")  # For calculating correlations with significance
library(corrplot)
library(Hmisc)

cor <- combined_data[-c(1, 2, 22)]

names(cor) <- c('BIO1', 'BIO10', 'BIO11', 'BIO12', 'BIO13', 'BIO14', 'BIO15', 'BIO16', 'BIO17',
                'BIO18', 'BIO19', 'BIO2', 'BIO3', 'BIO4', 'BIO5', 'BIO6', 'BIO7', 'BIO8', 'BIO9')

# Calculate the correlation matrix
cor_matrix <- rcorr(as.matrix(cor))

# Set the file name and dimensions
png(filename = "correlation.png", width = 3000, height = 2000, res = 300)

# Plot the correlation matrix with colors
corr_plot <- corrplot(cor_matrix$r, method = "color", type = "upper", 
                      tl.col = "black", tl.srt = 45, 
                      addCoef.col = "black", number.cex = 0.7,
                      col = colorRampPalette(c("blue", "white", "red"))(200))


dev.off()
################################################################
################################################################

##### Combining collinear variables


combination <- combined_data
combination$temperature_index <- (combination$wc2.1_2.5m_bio_1 + combination$wc2.1_2.5m_bio_10 + combination$wc2.1_2.5m_bio_11) / 3



combination$precipitation_index <- (combination$wc2.1_2.5m_bio_12 + combination$wc2.1_2.5m_bio_13) / 2


combination$seasonality_index <- (combination$wc2.1_2.5m_bio_4 + combination$wc2.1_2.5m_bio_15) / 2

cor(combination)
# Fit a new model using the indices
mode <- glm(species ~ temperature_index + precipitation_index + seasonality_index + wc2.1_2.5m_bio_2 + wc2.1_2.5m_bio_3 + wc2.1_2.5m_bio_5 + wc2.1_2.5m_bio_6 + wc2.1_2.5m_bio_7 + wc2.1_2.5m_bio_8 + wc2.1_2.5m_bio_9 + wc2.1_2.5m_bio_14 + wc2.1_2.5m_bio_16 + wc2.1_2.5m_bio_17 + wc2.1_2.5m_bio_18 + wc2.1_2.5m_bio_19, data = combination[, -c(1, 2)], family = binomial(link = "logit"))

summary(mode)

# Perform stepwise model selection based on AIC
steps <- step(mode, direction="both")
summary(steps)

BIC(steps, reduced_gam)



crops <- crop(biovars, map_extent)


names(crops)
# Create a temperature index
temperature_index <- (crops[[1]] + crops[[2]] + crops[[3]]) / 3

# Create a precipitation index
precipitation_index <- (crops[[4]] + crops[[5]]) / 2

# Create a seasonality index
seasonality_index <- (crops[[14]] + crops[[7]]) / 2


cr <- crops

names(crops) <- c("wc2.1_2.5m_bio_14", "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", 
                  "wc2.1_2.5m_bio_18", "wc2.1_2.5m_bio_19", "wc2.1_2.5m_bio_2",
                  "wc2.1_2.5m_bio_3", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6", "wc2.1_2.5m_bio_7",
                  "wc2.1_2.5m_bio_8", "wc2.1_2.5m_bio_9", "temperature_index", "precipitation_index", "seasonality_index")
# Update the RasterBrick with the new indices
crops <- addLayer(crops, temperature_index)
crops <- addLayer(crops, precipitation_index)
crops <- addLayer(crops, seasonality_index)

# Remove the original layers that contributed to the indices
crops <- dropLayer(crops, c(1, 2, 3, 4, 5, 14, 7))


# Make predictions using the modified RasterBrick
predict_g <- predict(crops, steps, type = 'response')

values(predict_g)

plot(predict_g)



set.seed(12)

# Fitting a GAM with all variables as smooth terms (excluding longitude and latitude)
gam1 <- gam(species ~ s(wc2.1_2.5m_bio_14) + s(wc2.1_2.5m_bio_16) + s(wc2.1_2.5m_bio_17) + 
              s(wc2.1_2.5m_bio_18) + s(wc2.1_2.5m_bio_19) + s(wc2.1_2.5m_bio_2) + 
              s(wc2.1_2.5m_bio_3) + s(wc2.1_2.5m_bio_5) + s(wc2.1_2.5m_bio_6) + 
              s(wc2.1_2.5m_bio_7) + s(wc2.1_2.5m_bio_8) + s(wc2.1_2.5m_bio_9) + 
              s(seasonality_index) + s(temperature_index) + s(precipitation_index), 
            data=combination, 
            family=binomial(link="logit"), method="REML", select=TRUE)

summary(gam1)

########## For the reduced version of the above model (Using only statistically significant variables)

gam2 <- gam(species ~ s(wc2.1_2.5m_bio_2) + s(wc2.1_2.5m_bio_5) +
              s(wc2.1_2.5m_bio_19) + s(temperature_index) + s(precipitation_index), 
            data=combination, 
            family=binomial(link="logit"))

summary(gam2)



BIC(gam2)
BIC(gams_reduced)



#####  Testing the model
set.seed(123) # for reproducibility
training <- createDataPartition(combination$species, p = 0.8, list = FALSE)
train <- combination[training, ]
tests <- combination[-training, ]

# Fit the models on the training data
gam_t1 <- gam(species ~ s(wc2.1_2.5m_bio_14) + s(wc2.1_2.5m_bio_16) + s(wc2.1_2.5m_bio_17) + 
                s(wc2.1_2.5m_bio_18) + s(wc2.1_2.5m_bio_19) + s(wc2.1_2.5m_bio_2) + 
                s(wc2.1_2.5m_bio_3) + s(wc2.1_2.5m_bio_5) + s(wc2.1_2.5m_bio_6) + 
                s(wc2.1_2.5m_bio_7) + s(wc2.1_2.5m_bio_8) + s(wc2.1_2.5m_bio_9) + 
                s(seasonality_index) + s(temperature_index) + s(precipitation_index), 
              data=trains, 
              family=binomial(link="logit"))

### The reduced version
gams_t2 <- gam(species ~ s(wc2.1_2.5m_bio_2) + s(wc2.1_2.5m_bio_5) +
                 s(wc2.1_2.5m_bio_19) + s(temperature_index) + s(precipitation_index), 
               data=trains, 
               family=binomial(link="logit"))
# Calculate predicted probabilities for the test set
predicted_probs_1 <- predict(gam_t1, newdata=tests, type="response")
predicted_probs_2 <- predict(gams_t2, newdata=tests, type="response")

# Calculate ROC curves
roc_curve_1 <- roc(tests$species, predicted_probs_1)
roc_curve_2 <- roc(tests$species, predicted_probs_2)

# Calculate AUC for each model
auc_1 <- auc(roc_curve_1)
auc_2 <- auc(roc_curve_2)

# Print AUC values
cat("AUC combined model", auc_1, "\n")
cat("AUC for combined model:", auc_2, "\n")

# Calculate the confusion matrix, specifying the positive class
cm_1 <- confusionMatrix(as.factor(round(predicted_probs_1)), tests$species, positive = "1")
cm_2 <- confusionMatrix(as.factor(round(predicted_probs_2)), tests$species, positive = "1")


precision_2 <- cm_2$byClass['Precision']
recall_2 <- cm_2$byClass['Recall']
f1_score_2 <- 2 * (precision_2 * recall_2) / (precision_2 + recall_2)

precision_1 <- cm_1$byClass['Precision']
recall_1 <- cm_1$byClass['Recall']
f1_score_1 <- 2 * (precision_1 * recall_1) / (precision_1 + recall_1)

# Print the F1 scores and precision for each model
cat("F1 Score combined Full model:", f1_score_1, "Precision:", precision_1, "\n")
cat("F1 Score combined model:", f1_score_2, "Precision:", precision_2, "\n")
#########################################################
###########################################################





##############################################
############## Principal component Analysis  ####################

com <- combined_data[, -c(1, 2, 22)]

combined_pca <- princomp(com, cor = TRUE)

summary(combined_pca)
plot(combined_pca)


names(combined_pca)

# Extract the scores (principal components) from the PCA object
pca_scores <- as.data.frame(combined_pca$scores)

# Add the longitude and latitude columns to the PCA scores data frame
pca_df <- pca_scores %>%
  mutate(decimalLongitude = combined_data$decimalLongitude, 
         decimalLatitude = combined_data$decimalLatitude, species = combined_data$species)

# If you want the longitude and latitude columns to be the first two columns
pca_df <- pca_df %>%
  select(decimalLongitude, decimalLatitude, everything())


# Mean of each variable
means <- combined_pca$center

# Standard deviation of each variable
sds <- combined_pca$scale

loadings <- combined_pca$loadings

##################################################


########################
pca_scores$species <- response




standardized_layers <- list()
for(i in 1:nlayers(crops)) {
  standardized_layer <- (crops[[i]] - means[i]) / sds[i]
  standardized_layers[[i]] <- standardized_layer
}
standardized_crops <- brick(standardized_layers)

# Convert standardized raster to matrix and apply loadings
standardized_matrix <- as.matrix(stack(standardized_crops))
pca_scores_matrix <- standardized_matrix %*% loadings



# Assuming 'pca_scores_matrix' is the matrix of PCA scores
# And 'selected_crops' is the RasterBrick of selected environmental variables

# Initialize an empty list to store the RasterLayers
pca_layers <- list()

# Loop through each column in the PCA scores matrix
for(i in 1:ncol(pca_scores_matrix)) {
  # Convert the column to a RasterLayer
  pca_layer <- raster(crops[[1]])  # Use the first layer of selected_crops as a template
  pca_layer <- setValues(pca_layer, pca_scores_matrix[, i])
  
  # Add the RasterLayer to the list
  pca_layers[[i]] <- pca_layer
}

# Combine the RasterLayers into a RasterBrick
pca_scores_raster <- brick(pca_layers)

# Set the names of the layers in the RasterBrick to match the PCA components
names(pca_scores_raster) <- paste0("Comp.", 1:ncol(pca_scores_matrix))

# Assuming 'pca_scores_raster' is your RasterBrick with PCA scores
#selected_pca_scores <- subset(pca_scores_raster, c("Comp.2", "Comp.3", "Comp.5", "Comp.6"))

# Use the selected PCA scores for prediction
predictors_pca <- predict(selected_pca_scores, model, type = 'response')

##### PCA GLM  #########

# Assuming 'species_presence' is a binary vector indicating species presence/absence
# and 'pca_scores' is a matrix of PCA scores for the retained components
model <- glm(response ~ Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6, data = pca_scores, family = binomial(link = "logit"))

summary(model)

# Perform stepwise model selection based on AIC
stepwise <- step(model, direction="both")
summary(stepwise)

s_model <- glm(response ~ Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.6, data = pca_scores, family = binomial(link = "logit"))

summary(s_model)

predictors_pca <- predict(pca_scores_raster, s_model, type = 'response')



values(predictors_pca)
plot(predictors_pca)


############### COMPARING PREDICTIVE POWER  ######################
set.seed(123) # Ensure reproducibility
splitIndexs <- createDataPartition(pca_scores$species, p = 0.8, list = FALSE)
glm_train_d <- pca_scores[splitIndexs, ]
glm_test_d <- pca_scores[-splitIndexs, ]

# Full model with all bioclimatic variables
full_mode <- glm(species ~ Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.5 + Comp.6, family=binomial(link="logit"), data=glm_train_d)

# Stepwise reduced model
step_mode <- glm(species ~  Comp.1 + Comp.2 + Comp.3 + Comp.4 + Comp.6, family = binomial(link = "logit"), data=glm_train_d)


# Predicted probabilities for the full model
predicted_full <- predict(full_mode, newdata=glm_test_d, type="response")

# Predicted probabilities for the stepwise reduced model
predicted_step <- predict(step_mode, newdata=glm_test_d, type="response")

# Confusion matrix for the full model
cm_f <- confusionMatrix(as.factor(round(predicted_full)), glm_test_d$species, positive = '1')

# Confusion matrix for the stepwise reduced model
cm_s <- confusionMatrix(as.factor(round(predicted_step)), glm_test_d$species, positive = '1')

# Calculate AUC for both models
roc_f <- roc(glm_test_d$species, predicted_full)
roc_s <- roc(glm_test_d$species, predicted_step)
auc_f <- auc(roc_f)
auc_s <- auc(roc_s)

# Print AUC values
cat("AUC for Full Model:", auc_f, "\n")
cat("AUC for Stepwise Reduced Model:", auc_s, "\n")
##### The result indicates that the mod model is best but because of ecologicial importance and uncertainty about future climate pattern
############################

############## PCA GAM ############
set.seed(11)

# Fitting a GAM with all variables as smooth terms (excluding longitude and latitude)
gam_m <- gam(species ~ s(Comp.1) + s(Comp.2) + s(Comp.3) + 
               s(Comp.4) + s(Comp.5) + s(Comp.6), data=pca_scores, 
             family=binomial(link="logit"), method="REML", select=TRUE)


summary(gam_m)

gam_O <- gam(species ~ s(Comp.1) + s(Comp.2) + s(Comp.3) + 
               s(Comp.4) + s(Comp.5) + s(Comp.6), data=pca_scores, 
             family=binomial(link="logit"))

summary(gam_O)


BIC(gams_red, gam_O)

#### REDUCED
gams_red <- gam(species ~ s(Comp.1) + s(Comp.2) + s(Comp.3) + 
                  s(Comp.4) + s(Comp.6), 
                data=pca_scores, family=binomial(link="logit"))

summary(gams_red)
AIC(gams_red)

set.seed(123) # for reproducibility
training_indexs <- createDataPartition(pca_scores$species, p = 0.8, list = FALSE)
trains <- pca_scores[training_indexs, ]
tests <- pca_scores[-training_indexs, ]

# Fit your models on the training data
gam_O <- gam(species ~ s(Comp.1) + s(Comp.2) + s(Comp.3) + 
               s(Comp.4) + s(Comp.5) + s(Comp.6), data=trains, 
             family=binomial(link="logit"))

gams_red <- gam(species ~ s(Comp.1) + s(Comp.2) + s(Comp.3) + 
                  s(Comp.4) + s(Comp.6), 
                data=trains, family=binomial(link="logit"))

# Calculate predicted probabilities for the test set
predicted_probs_O <- predict(gam_O, newdata=tests, type="response")
predicted_probs_red <- predict(gams_red, newdata=tests, type="response")

# Calculate ROC curves
roc_curve_O <- roc(tests$species, predicted_probs_O)
roc_curve_red <- roc(tests$species, predicted_probs_red)

# Calculate AUC for each model
auc_O <- auc(roc_curve_O)
auc_red <- auc(roc_curve_red)

# Print AUC values
cat("AUC PCA model", auc_O, "\n")
cat("AUC for downsized model:", auc_red, "\n")

# Calculate the confusion matrix, specifying the positive class
cm_O <- confusionMatrix(as.factor(round(predicted_probs_O)), tests$species, positive = "1")
cm_red <- confusionMatrix(as.factor(round(predicted_probs_red)), tests$species, positive = "1")

# Now the 'Positive' class should be correctly considered as '1'


# Extract precision and recall for each model
precision_red <- cm_red$byClass['Precision']
recall_red <- cm_red$byClass['Recall']
f1_score_red <- 2 * (precision_red * recall_red) / (precision_red + recall_red)

precision_O <- cm_O$byClass['Precision']
recall_O <- cm_O$byClass['Recall']
f1_score_O <- 2 * (precision_O * recall_O) / (precision_O + recall_O)

# Print the F1 scores and precision for each model
cat("F1 Score PCA:", f1_score_red, "Precision:", precision_red, "\n")
cat("F1 Score for FULL model:", f1_score_O, "Precision:", precision_O, "\n")

##### gams_red is the better model

##### Prediction
predictors_gam <- predict(pca_scores_raster, gams_red, type = 'response')
values(predictors_gam) <- na.omit(values(predictors_gam))
values(predictors_gam)
plot(predictors_gam)



###### Using the combined_data earlier created
predictors <- combined_data[, !(names(combined_data) %in% c("decimalLongitude", "decimalLatitude", "species"))]
response <- combined_data$species



##### FUTURE ######################

# Create a RasterStack from the multi-band raster file
bioclim_stacks <- stack("C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\Dissertation/wc2.1_2.5m_bioc_ACCESS-CM2_ssp126_2081-2100.tif")

# Check the RasterStack to confirm it contains all 19 bioclimatic variables
print(bioclim_stacks)
plot(bioclim_stack)

# Crop the global raster to the study extent
future_rast <- crop(bioclim_stack, map_extent)
plot(future_rast)

names(future_rast) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_3", 
                        "wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6", 
                        "wc2.1_2.5m_bio_7", "wc2.1_2.5m_bio_8", 
                        "wc2.1_2.5m_bio_9", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", 
                        "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_14", "
                        wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", 
                        "wc2.1_2.5m_bio_18", "wc2.1_2.5m_bio_19")


# Get the names of the layers in the same order as they appear in 'future_rast'
ordered_layer_names <- names(future_rast)

# Reorder 'means' and 'sds' to match the order of 'future_rast' layers
new_means <- means[ordered_layer_names]
new_sds <- sds[ordered_layer_names]

# Find the order of rows in 'loadings' that matches 'ordered_layer_names'
correct_order <- match(ordered_layer_names, rownames(loadings))

# Reorder the rows of 'loadings' to match the order of variables in 'future_rast'
new_loadings <- loadings[correct_order, ]

print(rownames(loadings))

# Now you can standardize the 'future_rast' using the reordered 'ordered_means' and 'ordered_sds'
standardized_future_layers <- list()
for(i in 1:nlayers(future_rast)) {
  standardized_future_layer <- (future_rast[[i]] - new_means[i]) / new_sds[i]
  standardized_future_layers[[i]] <- standardized_future_layer
}
standardized_future <- brick(standardized_future_layers)

# Convert standardized raster to matrix and apply loadings
standardized_future_matrix <- as.matrix(stack(standardized_future))
pca_future_scores_matrix <- standardized_future_matrix %*% new_loadings


# Initialize an empty list to store the RasterLayers
pca_future_layers <- list()

# Loop through each column in the PCA scores matrix
for(i in 1:ncol(pca_future_scores_matrix)) {
  # Convert the column to a RasterLayer
  pca_future_layer <- raster(future_rast[[1]])  # Use the first layer of selected_crops as a template
  pca_future_layer <- setValues(pca_future_layer, pca_future_scores_matrix[, i])
  
  # Add the RasterLayer to the list
  pca_future_layers[[i]] <- pca_future_layer
}

# Combine the RasterLayers into a RasterBrick
pca_future_scores_raster <- brick(pca_future_layers)

# Set the names of the layers in the RasterBrick to match the PCA components
names(pca_future_scores_raster) <- paste0("Comp.", 1:ncol(pca_future_scores_matrix))


############# For the future raster #################

# Create a data frame for all raster cells
all_df <- as.data.frame(pca_future_scores_raster, xy=TRUE)
names(all_df)
selected <- c(3, 4, 5, 6, 8)
future_data_predict <- select(all_df, all_of(selected))


# Predict for the full spatial extent
all_preds <- predict(gams_red, newdata=future_data_predict, type = 'response')

# Convert the predictions to a RasterLayer
full_pred_raster <- raster(pca_future_scores_raster)

# Convert all_predictions to a simple numeric vector
all_pred_vector <- as.vector(all_preds)

values(full_pred_raster) <- all_pred_vector

# Plot the predictions
plot(full_pred_raster, main="Predicted Species Distribution")




# Define the directory where your raster files are located
raster_n1 <- "C:\\Users\\user\\OneDrive\\Documents\\Dissertation\\Dissertation\\ssp126"

# List all files that end with '2081-2100.tif'
raster_files_n1 <- list.files(path = raster_n1, pattern = "ssp126_2081-2100.tif$", full.names = TRUE)

# Initialize a list to store prediction rasters
prediction_rasters_n1 <- list()

# Loop over each file
for (raster_file in raster_files_n1) {
  # Load the raster file
  bioclim_stack_n1 <- stack(raster_file)
  
  # Crop the global raster to the study extent
  future_rast <- crop(bioclim_stack_n1, map_extent)
  
  names(future_rast) <- c("wc2.1_2.5m_bio_1", "wc2.1_2.5m_bio_2", "wc2.1_2.5m_bio_3", 
                          "wc2.1_2.5m_bio_4", "wc2.1_2.5m_bio_5", "wc2.1_2.5m_bio_6", 
                          "wc2.1_2.5m_bio_7", "wc2.1_2.5m_bio_8", 
                          "wc2.1_2.5m_bio_9", "wc2.1_2.5m_bio_10", "wc2.1_2.5m_bio_11", 
                          "wc2.1_2.5m_bio_12", "wc2.1_2.5m_bio_13", "wc2.1_2.5m_bio_14", "
                        wc2.1_2.5m_bio_15", "wc2.1_2.5m_bio_16", "wc2.1_2.5m_bio_17", 
                          "wc2.1_2.5m_bio_18", "wc2.1_2.5m_bio_19")
  # Now you can standardize the 'future_rast' using the reordered 'ordered_means' and 'ordered_sds'
  standardized_future_layers <- list()
  for(i in 1:nlayers(future_rast)) {
    standardized_future_layer <- (future_rast[[i]] - new_means[i]) / new_sds[i]
    standardized_future_layers[[i]] <- standardized_future_layer
  }
  standardized_future <- brick(standardized_future_layers)
  
  # Convert standardized raster to matrix and apply loadings
  standardized_future_matrix <- as.matrix(stack(standardized_future))
  pca_future_scores_matrix <- standardized_future_matrix %*% new_loadings
  
  
  # Initialize an empty list to store the RasterLayers
  pca_future_layers <- list()
  
  # Loop through each column in the PCA scores matrix
  for(i in 1:ncol(pca_future_scores_matrix)) {
    # Convert the column to a RasterLayer
    pca_future_layer <- raster(future_rast[[1]])  # Use the first layer of selected_crops as a template
    pca_future_layer <- setValues(pca_future_layer, pca_future_scores_matrix[, i])
    
    # Add the RasterLayer to the list
    pca_future_layers[[i]] <- pca_future_layer
  }
  
  # Combine the RasterLayers into a RasterBrick
  pca_future_scores_raster <- brick(pca_future_layers)
  
  # Set the names of the layers in the RasterBrick to match the PCA components
  names(pca_future_scores_raster) <- paste0("Comp.", 1:ncol(pca_future_scores_matrix))
  
  # Create a data frame for all raster cells
  all_df <- as.data.frame(pca_future_scores_raster, xy=TRUE)
  
  # Select the necessary columns for prediction
  selected <- c(3, 4, 5, 6, 8)
  future_data_predict <- select(all_df, all_of(selected))
  
  ##### Remove NA  ################
  #future_data_predict <- na.omit(future_data_predict)
  
  # Predict for the full spatial extent
  all_preds <- predict(gams_red, newdata=future_data_predict, type = 'response')
  
  # Convert the predictions to a RasterLayer
  full_pred_raster <- raster(pca_future_scores_raster)
  
  # Convert all_predictions to a simple numeric vector
  all_pred_vector <- as.vector(all_preds)
  
  # Assign the predictions to the raster
  values(full_pred_raster) <- all_pred_vector
  
  # Store the prediction raster in the list
  prediction_rasters_n1[[basename(raster_file)]] <- full_pred_raster
  
  #print(names(prediction_rasters_2))
}

# After the loop, create an ensemble raster by calculating the mean of all prediction rasters
ensemble_raster_n1 <- stack(prediction_rasters_n1) |> calc(fun = mean, na.rm = TRUE)

plot(ensemble_raster_n1)





#############################################################
######### LATITUDE AND LONGITUDE SHIFT FOR VITIS ###################
# Define a threshold for presence, e.g., 0.5 or another value based on your model'soutput
presence_threshold <- 0.48
# Create binary rasters for presence (1) and absence (0)
# Assuming 'current_raster' and 'future_rast' are RasterLayer objects with presence probabilities
current_raster <- predictors_df 
future_rast <- ensemble_raster
future_rast2 <- ensemble_raster2
future_rast3 <- ensemble_raster3
future_rast4 <- ensemble_raster4
current_presence <- current_raster > presence_threshold
future_presence <- future_rast > presence_threshold
future_presence2 <- future_rast2 > presence_threshold
future_presence3 <- future_rast3 > presence_threshold
future_presence4 <- future_rast4 > presence_threshold
# Convert binary rasters to SpatialPolygonsDataFrame for presence areas
current_spdf <- rasterToPolygons(current_presence, fun=function(x){x==1}, dissolve=TRUE)
future_spdf <- rasterToPolygons(future_presence, fun=function(x){x==1}, dissolve=TRUE)
future_spdf2 <- rasterToPolygons(future_presence2, fun=function(x){x==1}, dissolve=TRUE)
future_spdf3 <- rasterToPolygons(future_presence3, fun=function(x){x==1}, dissolve=TRUE)
future_spdf4 <- rasterToPolygons(future_presence4, fun=function(x){x==1}, dissolve=TRUE)
# Calculate centroids of presence areas
current_centroid <- coordinates(sp::geometry(current_spdf))
future_centroid <- coordinates(sp::geometry(future_spdf))
future_centroid2 <- coordinates(sp::geometry(future_spdf2))
future_centroid3 <- coordinates(sp::geometry(future_spdf3))
future_centroid4 <- coordinates(sp::geometry(future_spdf4))
# Compute shifts
lat_shift <- future_centroid[2] - current_centroid[2]
lon_shift <- future_centroid[1] - current_centroid[1]
# Compute shifts for 2nd
lat_shift2 <- future_centroid2[2] - current_centroid[2]
lon_shift2 <- future_centroid2[1] - current_centroid[1]
# Compute shifts for 3rd
lat_shift3 <- future_centroid3[2] - current_centroid[2]
lon_shift3 <- future_centroid3[1] - current_centroid[1]
# Compute shifts for 4th
lat_shift4 <- future_centroid4[2] - current_centroid[2]
lon_shift4 <- future_centroid4[1] - current_centroid[1]
# Output the results
cat("Latitude Shift:", lat_shift, "degrees\n")
cat("Longitude Shift:", lon_shift, "degrees\n")
cat("Latitude Shift:", lat_shift2, "degrees\n")
cat("Longitude Shift:", lon_shift2, "degrees\n")
cat("Latitude Shift:", lat_shift3, "degrees\n")
cat("Longitude Shift:", lon_shift3, "degrees\n")
cat("Latitude Shift:", lat_shift4, "degrees\n")
cat("Longitude Shift:", lon_shift4, "degrees\n")
current_presence_cells
current_presence_count
current_presence
# Calculate the number of presence cells in each raster
current_presence_count <- sum(getValues(current_presence) == 1, na.rm = TRUE)
future_presence_count <- sum(getValues(future_presence) == 1, na.rm = TRUE)
future_presence_count2 <- sum(getValues(future_presence2) == 1, na.rm = TRUE)
future_presence_count3 <- sum(getValues(future_presence3) == 1, na.rm = TRUE)
future_presence_count4 <- sum(getValues(future_presence4) == 1, na.rm = TRUE)
current_presence_count
future_presence_count
future_presence_count2
future_presence_count3
future_presence_count4
current_df
#### Put all predictions in a vector format
all_pv <- na.omit(values(ensemble_raster))
all_pv2 <- na.omit(values(ensemble_raster2))
all_pv3 <- na.omit(values(ensemble_raster3))
all_pv4 <- na.omit(values(ensemble_raster4))
# If the data are not normally distributed or have unequal variances, use the Wilcoxon rank-sum test
wilcox_test <- wilcox.test(current_df, all_pv, alternative = "two.sided")
wilcox_test2 <- wilcox.test(current_df, all_pv2, alternative = "two.sided")
wilcox_test3 <- wilcox.test(current_df, all_pv3, alternative = "two.sided")
wilcox_test4 <- wilcox.test(current_df, all_pv4, alternative = "two.sided")
wilcox_test
wilcox_test2
wilcox_test3
wilcox_test4
colnames(data)
data$year
min(data$year)
min(na.omit(data$year))
# Define the latitude shifts for each SSP scenario
lat_shift_ssp126 <- -0.9315917
lat_shift_ssp245 <- -0.5178151
lat_shift_ssp370 <- -0.5651763
lat_shift_ssp585 <- 1.007641
# Define the time period over which the change occurs (from 2023 to the midpoint of 2081-2100)
time_period_years <- 2090 - 2024
number_of_decades <- time_period_years / 10
# Calculate the average rate of latitude change per decade for each scenario
rate_per_decade_ssp126 <- lat_shift_ssp126 / number_of_decades
rate_per_decade_ssp245 <- lat_shift_ssp245 / number_of_decades
rate_per_decade_ssp370 <- lat_shift_ssp370 / number_of_decades
rate_per_decade_ssp585 <- lat_shift_ssp585 / number_of_decades
# Output the results
cat("Average Rate of Latitude Change per Decade for SSP126:", rate_per_decade_ssp126, "degrees/decade\n")
cat("Average Rate of Latitude Change per Decade for SSP245:", rate_per_decade_ssp245, "degrees/decade\n")
cat("Average Rate of Latitude Change per Decade for SSP370:", rate_per_decade_ssp370, "degrees/decade\n")
cat("Average Rate of Latitude Change per Decade for SSP585:", rate_per_decade_ssp585, "degrees/decade\n")
##################################
rate_per_decade_ssp126 <- -0.1411503
rate_per_decade_ssp245 <- -0.07845683
rate_per_decade_ssp370 <- -0.08563277
rate_per_decade_ssp585 <- 0.1526729
# Conversion factor: 1 degree of latitude is approximately 111 kilometers
conversion_factor <- 111
# Convert degrees to kilometers
rate_km_per_decade_ssp126 <- rate_per_decade_ssp126 * conversion_factor
rate_km_per_decade_ssp245 <- rate_per_decade_ssp245 * conversion_factor
rate_km_per_decade_ssp370 <- rate_per_decade_ssp370 * conversion_factor
rate_km_per_decade_ssp585 <- rate_per_decade_ssp585 * conversion_factor
# Output the results
cat("Average Rate of Latitude Change per Decade for SSP126:", rate_km_per_decade_ssp126, "km/decade\n")
cat("Average Rate of Latitude Change per Decade for SSP245:", rate_km_per_decade_ssp245, "km/decade\n")
cat("Average Rate of Latitude Change per Decade for SSP370:", rate_km_per_decade_ssp370, "km/decade\n")
cat("Average Rate of Latitude Change per Decade for SSP585:", rate_km_per_decade_ssp585, "km/decade\n")




###  Normality Test ############################################333
# Create a histogram for the data
######### Convert current prediction to a Vector

par(mfrow = c(5, 2))

hist(current_df, breaks = 50, main = "Histogram of Current Predictions", xlab = "Predicted Values", col ='steelblue')
qqnorm(current_df)
qqline(current_df, col = "red", lwd = 2)

hist(all_pv, breaks = 50, main = "Histogram of SSP126 Predictions", xlab = "Predicted Values", col = 'steelblue')
qqnorm(all_pv)
qqline(all_pv, col = "red", lwd = 2)

hist(all_pv2, breaks = 50, main = "Histogram of SSP245 Predictions", xlab = "Predicted Values", col ='steelblue')
qqnorm(all_pv2)
qqline(all_pv2, col = "red", lwd = 2)

hist(all_pv3, breaks = 50, main = "Histogram of SSP370 Predictions", xlab = "Predicted Values", col = 'steelblue')
qqnorm(all_pv3)
qqline(all_pv3, col = "red", lwd = 2)

hist(all_pv4, breaks = 50, main = "Histogram of SSP585 Predictions", xlab = "Predicted Values", col ='steelblue')
qqnorm(all_pv4)
qqline(all_pv4, col = "red", lwd = 2)



# Wilcox test
wilcox_test_result <- wilcox.test(current_df, all_predictions_vector)
# Print the results
print(wilcox_test_result)
# Q-Q Plot for current raster values
qqnorm(current_df)
qqline(current_df, col = "red", lwd = 2)
qqnorm(all_pv)
qqline(all_pv, col = "red", lwd = 2)
# If the data are not normally distributed or have unequal variances, use the Wilcoxon rank-sum test  


#################
##########################  UNCERTAINTY ###################

# Stack the rasters to create a multi-layer raster object
ensemble_stack <- stack(prediction_rasters)
# Calculate the standard deviation across the ensemble models
uncertainty_sd1 <- calc(ensemble_stack, fun=sd)
# Define breaks and corresponding colors
breaks_u1 <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
colors_u1 <- c("#d73027", "#fc8d59", "#fee08b", '#91cf60', '#1a9850', 'blue')
# Plot the ensemble raster
plot(uncertainty_sd1, main="SDM Uncertainty Map - Standard Deviation (SSP126)",
     breaks=breaks_u1, col=colors_u1, axis.args=list(at=breaks_u1, labels=breaks_u1))
legend("bottomright",legend=c("None", "Low", "Moderate", "High", 'Very High'),
       fill=colors_u1, title="Variability")
# Create a histogram with breaks at intervals of 0.1 up to 0.5
hist_data <- hist(values(uncertainty_sd1), breaks=seq(0, 0.6, by=0.1), plot=FALSE)
plot(hist_data)
# Extract the counts for each interval
interval_counts <- hist_data$counts
interval_counts
# Assuming 'prediction_rasters' is a list of raster layers from your ensemble models



# Stack the rasters to create a multi-layer raster object
ensemble_stack2 <- stack(prediction_rasters_2)
# Calculate the standard deviation across the ensemble models
uncertainty_sd2 <- calc(ensemble_stack2, fun=sd)
# Define breaks and corresponding colors
breaks_u1 <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
colors_u1 <- c("#d73027", "#fc8d59", "#fee08b", '#91cf60', '#1a9850', 'blue')
# Plot the ensemble raster
plot(uncertainty_sd2, main="SDM Uncertainty Map - Standard Deviation (SSP245)", breaks=breaks_u1, col=colors_u1, axis.args=list(at=breaks_u1, labels=breaks_u1))
# Add a legend to the plot
legend("bottomright",legend=c("None", "Low", "Moderate", "High", 'Very High'),
       fill=colors_u1, title="Variability")
# Create a histogram with breaks at intervals of 0.1 up to 0.5
hist_data2 <- hist(values(uncertainty_sd2), breaks=seq(0, 0.6, by=0.1), plot=FALSE)
plot(hist_data2)
# Extract the counts for each interval
interval_counts2 <- hist_data2$counts
interval_counts2



# Stack the rasters to create a multi-layer raster object
ensemble_stack3 <- stack(prediction_rasters_3)
# Calculate the standard deviation across the ensemble models
uncertainty_sd3 <- calc(ensemble_stack3, fun=sd)
# Define breaks and corresponding colors
breaks_u1 <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
colors_u1 <- c("#d73027", "#fc8d59", "#fee08b", '#91cf60', '#1a9850', 'blue')
# Plot the ensemble raster
plot(uncertainty_sd3, main="SDM Uncertainty Map - Standard Deviation (SSP370)", 
     breaks=breaks_u1, col=colors_u1, axis.args=list(at=breaks_u1, labels=breaks_u1))
# Add a legend to the plot
legend("bottomright", legend=c("None", "Low", "Moderate", "High", 'Very High'), 
       fill=colors_u1, title="Variability")
# Create a histogram with breaks at intervals of 0.1 up to 0.5
hist_data3 <- hist(values(uncertainty_sd3), breaks=seq(0, 0.6, by=0.1), plot=FALSE)
plot(hist_data3)
# Extract the counts for each interval
interval_counts3 <- hist_data3$counts
interval_counts3


ensemble_stack4 <- stack(prediction_rasters_4)
# Calculate the standard deviation across the ensemble models
uncertainty_sd4 <- calc(ensemble_stack4, fun=sd)
# Define breaks and corresponding colors
breaks_u4 <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
colors_u4 <- c("#d73027", "#fc8d59", "#fee08b", '#91cf60', '#1a9850', 'blue')
# Plot the ensemble raster
plot(uncertainty_sd4, main="SDM Uncertainty Map - Standard Deviation (SSP585)",
     breaks=breaks_u4, col=colors_u4, axis.args=list(at=breaks_u4, labels=breaks_u4))
# Add a legend to the plot
legend("bottomright",
       legend=c("Very Low", "Low", "Moderate", "High", 'Very High', 'Extremely High'),
       fill=colors_u4, title="Variability")
# Create a histogram with breaks at intervals of 0.1 up to 0.5
hist_data4 <- hist(values(uncertainty_sd4), breaks=seq(0, 0.6, by=0.1), plot=FALSE)
plot(hist_data4)
# Extract the counts for each interval
interval_counts4 <- hist_data4$counts
interval_counts4

##############################################################

# Set up a 2x2 plotting area with adjusted margins
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1) + 0.1, oma = c(0, 0, 2, 0))

# Custom function to format y-axis labels in thousands without "K"
format_thousands <- function(x) {
  formatC(x / 1e3, format = "f", digits = 0)
}

# Plot the first histogram
plot(hist_data, main = "SSP126", col="green", border="black",
     xlab = "", ylab = "Frequency (Thousands)", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1, yaxt = "n")
axis(2, las = 1, at = axTicks(2), labels = format_thousands(axTicks(2)))

# Plot the second histogram
plot(hist_data2, main = "SSP245", col="red", border="black",
     xlab = "", ylab = "", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1, yaxt = "n")
axis(2, las = 1, at = axTicks(2), labels = format_thousands(axTicks(2)))

# Plot the third histogram
plot(hist_data3, main = "SSP370", col="lightblue", border="black",
     xlab = "Uncertainty", ylab = "Frequency (Thousands)", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1, yaxt = "n")
axis(2, las = 1, at = axTicks(2), labels = format_thousands(axTicks(2)))

# Plot the fourth histogram
plot(hist_data4, main = "SSP585", col="purple", border="black",
     xlab = "Uncertainty", ylab = "", 
     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.1, yaxt = "n")
axis(2, las = 1, at = axTicks(2), labels = format_thousands(axTicks(2)))

# Add a common title
mtext("Histograms of Uncertainty for Different SSP Scenarios", outer = TRUE, cex = 1.4, font = 2)




######################

####################################################################### 

par(mfrow = c(2, 2))

# Set up a 2x2 plotting area
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1) + 0.1, oma = c(0, 0, 2, 0))

# Define the number of colors you want to use
number_of_colors_e <- 4

# Generate a color gradient using terrain.colors
color_gradient_e <- terrain.colors(number_of_colors_e)

# Plot the raster with the color gradient
plot(difference_raster, 
     main="SSP126",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))


# Plot the raster with the color gradient
plot(difference_raster2, 
     main="SSP245",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))


# Plot the raster with the color gradient
plot(difference_raster3, 
     main="SSP370",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))

# Plot the raster with the color gradient
plot(difference_raster4, 
     main="SSP585",
     col=color_gradient_e,
     breaks=c(-1, -0.5, 0, 0.5, 1),
     axis.args=list(at=seq(-1, 1, by=0.5), labels=seq(-1, 1, by=0.5)))

# Add a common title
mtext("Differece Between Currrent and Ensemble for Different SSP Scenarios", outer = TRUE, cex = 1.1, font = 2)






#################################################################



######################################################
# Calculate the area of each cell
cell_areas <- area(ensemble_raster)

# Get the average area of the cells
average_cell_area <- mean(values(cell_areas), na.rm = TRUE)
print(average_cell_area)

# Calculate the total area
total_area_km2 <- presence_count * average_cell_area
print(total_area_km2)

# Example presence count
presence_count <- 95412


###########################################################


