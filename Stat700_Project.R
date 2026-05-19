library(nlme)
library(ggplot2)
data <- read.table("ACTG315.dat", header=TRUE)
colnames(data) <- c("No1", "No2", "patid", "t", "V", "CD4")
set.seed(123)
#only include first 3 months
data_sub <- subset(data, t <= 90)

#Linear mixed effect model
model_lme <- lme(V ~ t + t:CD4, 
               random = ~ 1 | patid, 
               data = data_sub, 
               method = "ML")
summary(model_lme) #AIC = 729.20

#Random coefficient mixed effect model
model_rcm <- lme(V ~ t + t:CD4, 
               random = ~ t + t:CD4 | patid, 
               data = data_sub, 
               method = "ML",
               control = lmeControl(opt = "optim"))
summary(model_rcm) #AIC = 727.63


#Nonlinear mixed effect
model_nlm <- nlme(V ~ b1*exp(-b2*t) + b3*exp(-b4*CD4),
                data = data_sub,
                fixed = b1 + b2 + b3 + b4 ~ 1,
                random = b1 ~ 1 | patid,
                start = c(b1=2.5, b2=0.05, b3=2.5, b4=0.001))
summary(model_nlm) #AIC = 602.04

#Nonlinear random coefficient mixed effect model --look at initial start values
model_nlrc <- nlme(V ~ b1*exp(-b2*t) + b3*exp(-b4*CD4),
                data = data_sub,
                fixed = b1 + b2 + b3 + b4 ~ 1,
                random = b1 + b3 ~ 1 | patid,
                start = c(b1=2.5, b2=0.05, b3=2.5, b4=0.001))
summary(model_nlrc) #AIC = 494.95

#Compare AIC values
aic_values <- AIC(model_lme, model_rcm, model_nlm, model_nlrc)
print(aic_values) #NLRC with clear lowest AIC

##Model Visualizations##

#Predictions for LME
#Level 1 = Subject-Specific -- blue line
#Level 0 = Population Average -- red line
data_sub$pred_indiv_lme <- predict(model_lme, level = 1)
data_sub$pred_pop_lme   <- predict(model_lme, level = 0)
#Plot LME
lme_plot <- ggplot(data_sub, aes(x = t, y = V, group = patid)) +
  # Raw data points
  geom_point(alpha = 0.4, color = "grey30") + 
  #Subject-specific predictions (Individualized fit)
  geom_line(aes(y = pred_indiv_lme), color = "blue", alpha = 0.4) +
  #Population-level prediction (The true average of your model)
  geom_line(aes(y = pred_pop_lme, group = NULL), color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "LME Model Fit", 
    subtitle = "Blue: Subject-Specific | Red: Population Average",
    y = "Viral Load (V)", 
    x = "Days (t)"
  )

#Predictions for RCM
data_sub$pred_indiv_rcm <- predict(model_rcm, level = 1)
data_sub$pred_pop_rcm   <- predict(model_rcm, level = 0)

#RCM plot
rcm_plot <- ggplot(data_sub, aes(x = t, y = V, group = patid)) +
  #Raw data points
  geom_point(alpha = 0.4, color = "grey30") + 
  #Subject-specific predictions
  geom_line(aes(y = pred_indiv_rcm), color = "blue", alpha = 0.4) +
  #Population-level prediction
  geom_line(aes(y = pred_pop_rcm, group = NULL), color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "RCM Model Fit", 
    subtitle = "Blue: Subject-Specific | Red: Population Average",
    y = "Viral Load (V)", 
    x = "Days (t)"
  )

#Predictions for NLM
data_sub$pred_indiv_nlm <- predict(model_nlm, level = 1)
data_sub$pred_pop_nlm   <- predict(model_nlm, level = 0)

#NLM plot
nlm_plot <- ggplot(data_sub, aes(x = t, y = V, group = patid)) +
  # Raw data points
  geom_point(alpha = 0.4, color = "grey30") + 
  # Subject-specific predictions
  geom_line(aes(y = pred_indiv_nlm), color = "blue", alpha = 0.4) +
  # Population-level prediction
  geom_line(aes(y = pred_pop_nlm, group = NULL), color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "NLM Model Fit", 
    subtitle = "Blue: Subject-Specific | Red: Population Average",
    y = "Viral Load (V)", 
    x = "Days (t)"
  )

#Predictions for NLRC
data_sub$pred_indiv_nlrc <- predict(model_nlrc, level = 1)
data_sub$pred_pop_nlrc   <- predict(model_nlrc, level = 0)

#NLRC plot
nlrc_plot <- ggplot(data_sub, aes(x = t, y = V, group = patid)) +
  # Raw data points
  geom_point(alpha = 0.4, color = "grey30") + 
  # Subject-specific predictions
  geom_line(aes(y = pred_indiv_nlrc), color = "blue", alpha = 0.4) +
  # Population-level prediction
  geom_line(aes(y = pred_pop_nlrc, group = NULL), color = "red", linewidth = 1.2) +
  theme_minimal() +
  labs(
    title = "NLRC Model Fit", 
    subtitle = "Blue: Subject-Specific | Red: Population Average",
    y = "Viral Load (V)", 
    x = "Days (t)"
  )


# Initialize a list of models
model_list <- list("LME" = model_lme, "RCM" = model_rcm, "NLM" = model_nlm, 
                   "NLRC" = model_nlrc)

# Define the parameter names as they should appear in the table
param_names <- c("Beta 1", "Beta 2", "Beta 3", "Beta 4", "Sigma^2", "AIC")

# Create an empty matrix to store results
results_matrix <- matrix(NA, nrow = length(param_names), ncol = 4)
colnames(results_matrix) <- c("LME", "RCM", "NLM", "NLRC")
rownames(results_matrix) <- param_names

# Loop through models and manually extract values
for (i in 1:length(model_list)) {
  mod <- model_list[[i]]
  f_eff <- fixef(mod)
  
  # Fill Fixed Effects
  results_matrix[1, i] <- round(f_eff[1], 4) # Beta 1
  results_matrix[2, i] <- round(f_eff[2], 4) # Beta 2
  results_matrix[3, i] <- round(f_eff[3], 4) # Beta 3
  
  # Models C and D have a 4th beta
  if (length(f_eff) >= 4) {
    results_matrix[4, i] <- round(f_eff[4], 4)
  }
  
  # Extract Sigma^2 (Residual Variance)
  results_matrix[5, i] <- round(mod$sigma^2, 4)
  
  # Extract AIC
  results_matrix[6, i] <- round(AIC(mod), 2)
}

flipped_matrix <- t(results_matrix)

# Convert to data frame for easier handling
final_table_flipped <- as.data.frame(flipped_matrix)

##Extract subject-specific predictions##
for (m_name in names(model_list)) {
  mod <- model_list[[m_name]]
  
  cat("\n--- Subject-Specific Predictions for Model", m_name, "---\n")
  
  # Obtain predicted values (V-hat_ij) for each observation
  # these use the individual-level coefficients
  subject_predictions <- predict(mod, level = 1)
  
  # Combine with original data for a clear output
  output_data <- data.frame(
    patid = data_sub$patid,
    Time  = data_sub$t,
    Observed_V = data_sub$V,
    Predicted_V = subject_predictions
  )
  
  print(head(output_data))
}

#Combines all results into one data frame
all_results <- data.frame(
  Observed = data_sub$V,
  Pred_LME   = predict(model_lme, level = 1),
  Pred_RCM   = predict(model_rcm, level = 1),
  Pred_NLM   = predict(model_nlm, level = 1),
  Pred_NLRC   = predict(model_nlrc, level = 1)
)

#Absolute error of each model
all_results$Error_LME <- abs(all_results$Observed - all_results$Pred_LME)
all_results$Error_RCM <- abs(all_results$Observed - all_results$Pred_RCM)
all_results$Error_NLM <- abs(all_results$Observed - all_results$Pred_NLM)
all_results$Error_NLRC <- abs(all_results$Observed - all_results$Pred_NLRC)

#Identify which models have the smallest errors
error_cols <- all_results[, c("Error_LME", "Error_RCM", "Error_NLM", "Error_NLRC")]
all_results$Best_Model_Index <- apply(error_cols, 1, which.min)

#Index model names
model_names <- c("Model LME", "Model RCM", "Model NLM", "Model NLRC")
all_results$Best_Model_Name <- model_names[all_results$Best_Model_Index]

#Prediction counts based on smallest errors
prediction_counts <- table(all_results$Best_Model_Name)
print(prediction_counts)

#Percentages from the prediction counts
prediction_percentages <- prop.table(prediction_counts) * 100
print(prediction_percentages)

#Obtain beta coefficients for each patient
fixef_summary <- list(
  LME  = fixef(model_lme),
  RCM  = fixef(model_rcm),
  NLM  = fixef(model_nlm),
  NLRC = fixef(model_nlrc)
)

fixef_summary

#Extract BLUPs for Linear Mixed Model
mu_lme <- ranef(model_lme)
head(mu_lme)

# Extract BLUPs for the Random Coefficient Model
mu_rcm <- ranef(model_rcm)
head(mu_rcm)

# Extract BLUPs for the Nonlinear Mixed Model
mu_nlm <- ranef(model_nlm)
head(mu_nlm)

# Extract BLUPs for the Nonlinear Random Coefficient Model
mu_nlrc <- ranef(model_nlrc)
head(mu_nlrc)

#Obtain beta_i for each model
beta_i_lme  <- coef(model_lme)
beta_i_rcm  <- coef(model_rcm)
beta_i_nlm  <- coef(model_nlm)
beta_i_nlrc <- coef(model_nlrc)

#View results for beta_i for each model
head(beta_i_lme)
head(beta_i_rcm)
head(beta_i_nlm)
head(beta_i_nlrc)









