####################################################################################################
### R script to run PLSR models to predict traits on leaves subject to alcohol/drying treatments.###
### This script requires the completion of "Spec_prep_trait_analysis.R" to load in               ###
### required dataframes and process reflectance spectra.                                         ###

# Set Working directory
setwd("C:/Users/rfortier/Dropbox/MBG Postdoc/Specimen prep/Data analysis")
#setwd("~/Library/CloudStorage/Dropbox/MBG Postdoc/Specimen prep/Data analysis")

# Load more libraries
library(caret)
library(pls)
library(MLmetrics)
library(prospectr)
library(plsVarSel)
library(parallel)
library(pbmcapply)
library(future.apply)
library(progressr) 

# If not done already, run trait script first to create required dataframes
#source("Spec_prep_trait_analysis.R")

### Model traits from spectra for all data. We will use only dry spectra first.
# Define traits to model
traits_list <- c("dry_thickness", "dry_LMA") 

### Global model ###
# Number of iterations
niterations <- 100

# Storage for global model results
global_predictions <- list()
global_stats <- list()

for(iteration in 1:niterations) {
  
  cat("Running global model iteration", iteration, "\n")
  start_time <- Sys.time()
  
  # Split data 60/40 by specimen (voucher + treatment combination)
  unique_specimens <- unique(dry_leaf_deriv$specimen_id)
  n_specimens <- length(unique_specimens)
  training_n <- round(n_specimens * 0.6)
  
  training_specimens <- sample(unique_specimens, training_n)
  validation_specimens <- setdiff(unique_specimens, training_specimens)
  
  calibrationData <- dry_leaf_deriv[dry_leaf_deriv$specimen_id %in% training_specimens, ]
  validationData <- dry_leaf_deriv[dry_leaf_deriv$specimen_id %in% validation_specimens, ]
  
  # Prepare spectral data
  X_train <- calibrationData %>% select(all_of(deriv_wavelength_cols))
  X_val <- validationData %>% select(all_of(deriv_wavelength_cols))
  
  # Store results for this iteration
  iter_predictions <- list()
  iter_stats <- data.frame()
  
  # Model each trait
  for(trait in traits_list) {
    
    Y_train <- calibrationData[[trait]]
    Y_val <- validationData[[trait]]
    
    # Build PLSR model
    model <- plsr(Y_train ~ as.matrix(X_train),
                  ncomp = 30,
                  method = "oscorespls",
                  validation = "CV",
                  segments = 10)
    
    ncomp <- selectNcomp(model, method = "onesigma", plot = FALSE)
    
    # Predict on validation data
    val_pred <- as.vector(predict(model, newdata = as.matrix(X_val), ncomp = ncomp)[,,1])
    val_fit <- lm(Y_val ~ val_pred)
    
    # Store predictions
    iter_predictions[[trait]] <- data.frame(
      iteration = iteration,
      trait = trait,
      voucher_no = validationData$voucher_no,
      treatment_code = validationData$treatment_code,
      measured = Y_val,
      predicted = val_pred
    )
    
    # Store statistics
    iter_stats <- rbind(iter_stats, data.frame(
      iteration = iteration,
      trait = trait,
      R2 = summary(val_fit)$r.squared,
      RMSE = RMSD(Y_val, val_pred),
      perRMSE = percentRMSD(Y_val, val_pred, 0.025, 0.975),
      bias = mean(val_pred, na.rm = TRUE) - mean(Y_val, na.rm = TRUE),
      ncomp = ncomp
    ))
  }
  
  global_predictions[[iteration]] <- do.call(rbind, iter_predictions)
  global_stats[[iteration]] <- iter_stats
  
  end_time <- Sys.time()
  cat("  Completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
}

# Combine results
global_predictions_df <- do.call(rbind, global_predictions)
global_stats_df <- do.call(rbind, global_stats)

# Summarize global model performance
global_summary <- global_stats_df %>%
  group_by(trait) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    mean_perRMSE = mean(perRMSE),
    sd_perRMSE = sd(perRMSE),
    mean_bias = mean(bias),
    sd_bias = sd(bias),
    mean_ncomp = mean(ncomp),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

global_stats_df <- global_stats_df %>%
  mutate(trait = recode(trait,
                        dry_thickness = "Thickness",
                        dry_LMA = "LMA"))

# Plot model performance
ggplot(global_stats_df, aes(x = trait, y = R2)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_classic(base_size = 14) +
  labs(x = "Trait",
       y = "R²") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot predicted vs measured values
# First summarize each leaf across all iterations
global_pred_summary <- global_predictions_df %>%
  group_by(trait, voucher_no, treatment_code, measured) %>%
  summarise(
    mean_predicted = mean(predicted, na.rm = TRUE),
    sd_predicted = sd(predicted, na.rm = TRUE),
    n_iterations = n(),
    .groups = "drop"
  ) %>%
  mutate(trait = recode(trait,
                        dry_LMA =  "LMA",
                        dry_thickness = "Thickness"))

facet_stats <- global_stats_df %>%
  group_by(trait) %>%
  summarise(
    R2 = mean(R2, na.rm = TRUE),
    perRMSE = mean(perRMSE, na.rm = TRUE),
    .groups = "drop"
  )

global_plots <- ggplot(global_pred_summary, aes(x = mean_predicted, y = measured)) +
  geom_errorbar(aes(xmin = mean_predicted - sd_predicted,
                    xmax = mean_predicted + sd_predicted),
                alpha = .3, colour = "grey50") +
  geom_point(alpha = .6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "red", linewidth = 1) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "blue", linewidth = 1) +
  geom_text(data = facet_stats,
            aes(x = -Inf, y = Inf,
                label = paste0("Mean R² = ", round(R2, 2),
                               "\nMean %RMSE = ", round(perRMSE*100, 2))),
            hjust = -.05, vjust = 1.1, size = 4,
            inherit.aes = FALSE) +
  facet_wrap(~trait, scales = "free") +
  labs(x = "", y = "",
       title = "Global model") +
  theme_bw()

# Make supplementary table of PLSR model results
#write.csv(global_summary, file = "TableS7.csv", row.names = F)

### Global model rarefied down to 5 leaves per tree - "mixed-treatment" model ###
# Add a code for each leaf
dry_leaf_deriv$leaf_code <- paste0(dry_leaf_deriv$voucher_no, "_", dry_leaf_deriv$treatment_code, "_", dry_leaf_deriv$leaf)
# For each iteration, randomly select five leaves per voucher_no before splitting.
global_mix_predictions <- list()
global_mix_stats <- list()

for(iteration in 1:niterations) {
  
  cat("Running mixed-treatment model", iteration, "\n")
  start_time <- Sys.time()
  
  # Sample five leaves per voucher_no
  leaf_ids <- dry_leaf_deriv %>%
    distinct(voucher_no, leaf_code) %>%
    group_by(voucher_no) %>%
    slice_sample(n = 5) %>%
    ungroup()
  
  dry_leaf_mix <- dry_leaf_deriv %>%
    semi_join(leaf_ids, by = c("voucher_no", "leaf_code"))
  
  # Split data 60/40 stratified by voucher_no
  leaf_split <- dry_leaf_mix %>%
    distinct(voucher_no, leaf_code) %>%
    group_by(voucher_no) %>%
    slice_sample(prop = 0.6) %>%
    ungroup()
  
  calibrationData <- dry_leaf_mix %>%
    semi_join(leaf_split, by = c("voucher_no", "leaf_code"))
  validationData <- dry_leaf_mix %>%
    anti_join(leaf_split, by = c("voucher_no", "leaf_code"))
  
  # Prepare spectral data
  X_train <- calibrationData %>% select(all_of(deriv_wavelength_cols))
  X_val   <- validationData   %>% select(all_of(deriv_wavelength_cols))
  
  iter_predictions <- list()
  iter_stats <- data.frame()
  
  for(trait in traits_list) {
    
    Y_train <- calibrationData[[trait]]
    Y_val   <- validationData[[trait]]
    
    model <- plsr(Y_train ~ as.matrix(X_train),
                  ncomp = 30,
                  method = "oscorespls",
                  validation = "CV",
                  segments = 10)
    
    ncomp <- selectNcomp(model, method = "onesigma", plot = FALSE)
    
    val_pred <- as.vector(predict(model, newdata = as.matrix(X_val), ncomp = ncomp)[,,1])
    val_fit  <- lm(Y_val ~ val_pred)
    
    iter_predictions[[trait]] <- data.frame(
      iteration = iteration,
      trait = trait,
      voucher_no = validationData$voucher_no,
      treatment_code = validationData$treatment_code,
      measured = Y_val,
      predicted = val_pred
    )
    
    iter_stats <- rbind(iter_stats, data.frame(
      iteration = iteration,
      trait = trait,
      R2 = summary(val_fit)$r.squared,
      RMSE = RMSD(Y_val, val_pred),
      perRMSE = percentRMSD(Y_val, val_pred, 0.025, 0.975),
      bias = mean(val_pred, na.rm = TRUE) - mean(Y_val, na.rm = TRUE),
      ncomp = ncomp
    ))
  }
  
  global_mix_predictions[[iteration]] <- do.call(rbind, iter_predictions)
  global_mix_stats[[iteration]] <- iter_stats
  
  end_time <- Sys.time()
  cat("  Completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
}

# Combine results
global_mix_predictions_df <- do.call(rbind, global_mix_predictions)
global_mix_stats_df       <- do.call(rbind, global_mix_stats)

# Summarize
global_mix_summary <- global_mix_stats_df %>%
  group_by(trait) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    mean_perRMSE = mean(perRMSE),
    sd_perRMSE = sd(perRMSE),
    mean_bias = mean(bias),
    sd_bias = sd(bias),
    mean_ncomp = mean(ncomp),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

global_mix_stats_df <- global_mix_stats_df %>%
  mutate(trait = recode(trait,
                        dry_thickness = "Thickness",
                        dry_LMA = "LMA"),
         treatment = "mixed")

# Plot model performance
ggplot(global_mix_stats_df, aes(x = trait, y = R2)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_classic(base_size = 14) +
  labs(x = "Trait",
       y = "R²",
       title = "Mixed-treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot predicted vs measured
global_mix_pred_summary <- global_mix_predictions_df %>%
  group_by(trait, voucher_no, treatment_code, measured) %>%
  summarise(
    mean_predicted = mean(predicted, na.rm = TRUE),
    sd_predicted   = sd(predicted, na.rm = TRUE),
    n_iterations   = n(),
    .groups = "drop"
  ) %>%
  mutate(trait = recode(trait,
                        dry_LMA = "LMA",
                        dry_thickness = "Thickness"))

facet_stats_mix <- global_mix_stats_df %>%
  group_by(trait) %>%
  summarise(
    R2 = mean(R2, na.rm = TRUE),
    perRMSE = mean(perRMSE, na.rm = TRUE),
    .groups = "drop"
  )
facet_stats_mix$treatment <- "mixed"

mix_plots <- ggplot(global_mix_pred_summary, aes(x = mean_predicted, y = measured)) +
  geom_errorbar(aes(xmin = mean_predicted - sd_predicted,
                    xmax = mean_predicted + sd_predicted),
                alpha = .3, colour = "grey50") +
  geom_point(alpha = .6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "red", linewidth = 1) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "blue", linewidth = 1) +
  geom_text(data = facet_stats_mix,
            aes(x = -Inf, y = Inf,
                label = paste0("Mean R² = ", round(R2, 2),
                               "\nMean %RMSE = ", round(perRMSE * 100, 2))),
            hjust = -.05, vjust = 1.1, size = 4,
            inherit.aes = FALSE) +
  facet_wrap(~trait, scales = "free") +
  labs(x = "", y = "",
       title = "Mixed-treatment") +
  theme_bw()

global_mix_plots <- ggarrange(plotlist = c(global_plots, mix_plots), ncol = 1, nrow = 2) 
annotate_figure(p = global_mix_plots, left = text_grob("Measured value", rot = 90, vjust = 1, size = 14),
                bottom = text_grob("Predicted value",vjust = -0.5, size = 14))


### Treatment specific models
# Combine the two alcohol treatments for all subsequent analyses
dry_leaf_deriv2 <- dry_leaf_deriv %>%
  mutate(alcohol_treatment = ifelse(alcohol_treatment == "control", "control", "alcohol"))

# Get unique treatment combinations
treatment_combos <- dry_leaf_deriv2 %>%
  distinct(dry_treatment, alcohol_treatment) %>%
  mutate(treatment_label = paste(dry_treatment, alcohol_treatment, sep = "_")) %>%
  mutate(dry_treatment = as.character(dry_treatment),
         alcohol_treatment = as.character(alcohol_treatment))

treatment_predictions <- list()
treatment_stats <- list()

for(t in 1:nrow(treatment_combos)) {
  
  current_dry <- treatment_combos$dry_treatment[t]
  current_alc <- treatment_combos$alcohol_treatment[t]
  treatment_label <- treatment_combos$treatment_label[t]
  
  cat("\nProcessing treatment:", treatment_label, "\n")
  
  # Filter data for this treatment
  treatment_data_full <- dry_leaf_deriv2 %>%
    filter(dry_treatment == current_dry & alcohol_treatment == current_alc)
  
  # If this is an alcohol treatment, subset to match control sample size
  if(current_alc == "alcohol") {
    control_colls <- dry_leaf_deriv2 %>%
      filter(dry_treatment == current_dry & alcohol_treatment == "control") %>%
      pull(image_name) %>%
      unique()
    
    alcohol_colls <- treatment_data_full %>%
      pull(image_name) %>%
      unique()
    
    if(length(alcohol_colls) > length(control_colls)) {
      set.seed(789 + t)
      
      # Sample 50% of vouchers per species (at least 1 per species)
      selected_colls <- treatment_data_full %>%
        distinct(Species, image_name) %>%
        group_by(Species) %>%
        slice_sample(prop = 0.5) %>%
        ungroup() %>%
        pull(image_name)
      
      treatment_data <- treatment_data_full %>%
        filter(image_name %in% selected_colls)
      
      cat("  Subsampled alcohol data:", length(alcohol_colls), "->",
          length(selected_colls), "vouchers\n")
    }
  } else {
    treatment_data <- treatment_data_full
  }
  
  treat_iter_preds <- list()
  treat_iter_stats <- list()
  
  for(iteration in 1:niterations) {
    
    # Split by specimen
    training_leaves <- treatment_data %>%
      distinct(specimen_id, leaf_code) %>%
      group_by(specimen_id) %>%
      slice_sample(prop = 0.6) %>%
      ungroup() %>%
      pull(leaf_code)
    
    calibrationData <- treatment_data %>% filter(leaf_code %in% training_leaves)
    validationData  <- treatment_data %>% filter(!leaf_code %in% training_leaves)
    
    # Prepare spectral data
    X_train <- calibrationData %>% select(all_of(deriv_wavelength_cols))
    X_val <- validationData %>% select(all_of(deriv_wavelength_cols))
    
    iter_preds <- list()
    iter_stats <- data.frame()
    
    for(trait in traits_list) {
      
      Y_train <- calibrationData[[trait]]
      Y_val <- validationData[[trait]]
      
      # Build model
      model <- plsr(Y_train ~ as.matrix(X_train),
                    ncomp = 30,
                    method = "oscorespls",
                    validation = "CV",
                    segments = 10)
      
      ncomp <- selectNcomp(model, method = "onesigma", plot = FALSE)
      
      # Predict
      val_pred <- as.vector(predict(model, newdata = as.matrix(X_val), ncomp = ncomp)[,,1])
      val_fit <- lm(Y_val ~ val_pred)
      
      # Store results
      iter_preds[[trait]] <- data.frame(
        iteration = iteration,
        treatment = treatment_label,
        trait = trait,
        voucher_no = validationData$voucher_no,
        measured = Y_val,
        predicted = val_pred
      )
      
      iter_stats <- rbind(iter_stats, data.frame(
        iteration = iteration,
        treatment = treatment_label,
        trait = trait,
        R2 = summary(val_fit)$r.squared,
        RMSE = RMSD(Y_val, val_pred),
        perRMSE = percentRMSD(Y_val, val_pred, 0.025, 0.975),
        bias = mean(val_pred, na.rm = TRUE) - mean(Y_val, na.rm = TRUE),
        ncomp = ncomp
      ))
    }
    
    treat_iter_preds[[iteration]] <- do.call(rbind, iter_preds)
    treat_iter_stats[[iteration]] <- iter_stats
  }
  
  treatment_predictions[[treatment_label]] <- do.call(rbind, treat_iter_preds)
  treatment_stats[[treatment_label]] <- do.call(rbind, treat_iter_stats)
  
  cat("  Completed", niterations, "iterations for", treatment_label, "\n")
}

# Combine treatment results
treatment_stats_df <- do.call(rbind, treatment_stats)
treatment_predictions_df <- do.call(rbind, treatment_predictions)

# Add mixed treatment
treatment_stats_df <- rbind(treatment_stats_df, global_mix_stats_df)
treatment_stats_df <- treatment_stats_df %>%
  mutate(treatment = factor(treatment, levels = c("mixed", "control_control", "extra dry_control", "control_alcohol", "extra dry_alcohol")))

# Summarize treatment-specific performance
treatment_summary <- treatment_stats_df %>%
  group_by(treatment, trait) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    mean_perRMSE = mean(perRMSE),
    sd_perRMSE = sd(perRMSE),
    mean_bias = mean(bias),
    sd_bias = sd(bias),
    mean_ncomp = mean(ncomp),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) 
  
#write.csv(treatment_summary, file = "TableS7b.csv", row.names = F)

# Create a plot for each treatment combination
# Summarize predictions: mean and SD for each specimen across iterations
treatment_pred_summary <- treatment_predictions_df %>%
  group_by(treatment, trait, voucher_no, measured) %>%
  summarise(
    mean_predicted = mean(predicted, na.rm = TRUE),
    sd_predicted = sd(predicted, na.rm = TRUE),
    n_iterations = n(),
    .groups = "drop") %>%
  mutate(treatment = factor(treatment, levels = c("control_control", "extra dry_control", "control_alcohol", "extra dry_alcohol"), ordered = TRUE))

treatment_labels <- c(
  "control_control"       = "DC - AC",
  "extra dry_control"     = "DE - AC",
  "control_alcohol"       = "DC - A50/80",
  "extra dry_alcohol"     = "DE - A50/80")
treatment_plot_list <- list()

for(treatment_name in levels(treatment_pred_summary$treatment)) {
  
  treatment_data <- treatment_pred_summary %>%
    filter(treatment == treatment_name)
  
  for(trait in traits_list) {
    
    trait_data <- treatment_data %>% filter(trait == !!trait)
    
    plot_stats <- treatment_summary %>%
      filter(treatment == treatment_name & trait == !!trait)
    r2_val   <- plot_stats$mean_R2
    per_rmse <- plot_stats$mean_perRMSE
    
    p <- ggplot(trait_data, aes(x = mean_predicted, y = measured)) +
      geom_errorbarh(aes(xmin = mean_predicted - sd_predicted,
                         xmax = mean_predicted + sd_predicted),
                     alpha = 0.5, color = "gray50") +
      geom_point(alpha = 0.5, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "red", linewidth = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
      annotate("text", x = -Inf, y = Inf,
               label = paste0(treatment_labels[[treatment_name]],
                              "\nR² = ", r2_val,
                              "\n%RMSE = ", per_rmse * 100),
               hjust = -0.1, vjust = 1.1, size = 4) +
      labs(x = "", y = "") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    
    treatment_plot_list[[paste(treatment_name, trait, sep = "_")]] <- p
  }
}

treatment_plots <- ggarrange(plotlist = treatment_plot_list, ncol = 2, nrow = 4) 
annotate_figure(p = treatment_plots, left = text_grob("Measured value", rot = 90, vjust = 1, size = 14),
                bottom = text_grob("Predicted value",vjust = -0.5, size = 14),
                top = text_grob("Thickness (mm)                         LMA (g m-2)", size = 14))


# ANOVAs on model results to test for treatment effect
r2_cld_list <- list()

treatment_stats_df$trait <- ifelse(treatment_stats_df$trait == "dry_LMA", "LMA", 
                                   ifelse(treatment_stats_df$trait == "dry_thickness", "Thickness", treatment_stats_df$trait))

for(tr in unique(treatment_stats_df$trait)) {
  trait_data <- treatment_stats_df %>% filter(trait == tr)
  
  anova_model <- aov(R2 ~ treatment, data = trait_data)
  tukey_result <- TukeyHSD(anova_model)
  tukey_cld <- multcompLetters4(anova_model, tukey_result)
  
  global_max <- max(trait_data$R2, na.rm = TRUE)
  
  r2_cld_list[[tr]] <- data.frame(
    treatment = names(tukey_cld$treatment$Letters),
    Letters   = tukey_cld$treatment$Letters
  ) %>%
    mutate(trait = tr,
           y_pos = global_max)
}

cld_r2 <- do.call(rbind, r2_cld_list)
cld_r2 <- cld_r2 %>%
  mutate(treatment = factor(treatment, levels = c("mixed", "control_control", "extra dry_control", "control_alcohol", "extra dry_alcohol")))


a <- ggplot(treatment_stats_df, aes(x = trait, y = R2, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = cld_r2, aes(x = trait, y = y_pos, label = Letters, group = treatment),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "R²", fill = "Treatment") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_discrete(labels = c("Mixed", "Dc - Ac", "De - Ac", "Dc - A50/80", "De - A50/80")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

rmse_cld_list <- list()

for(tr in unique(treatment_stats_df$trait)) {
  trait_data <- treatment_stats_df %>% filter(trait == tr)
  
  anova_model <- aov(perRMSE ~ treatment, data = trait_data)
  tukey_result <- TukeyHSD(anova_model)
  tukey_cld <- multcompLetters4(anova_model, tukey_result)
  
  global_max <- max(trait_data$perRMSE, na.rm = TRUE)
  
  rmse_cld_list[[tr]] <- data.frame(
    treatment = names(tukey_cld$treatment$Letters),
    Letters   = tukey_cld$treatment$Letters
  ) %>%
    mutate(trait = tr,
           y_pos = global_max)
}

cld_rmse <- do.call(rbind, rmse_cld_list)

b <- ggplot(treatment_stats_df, aes(x = trait, y = perRMSE, fill = treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = cld_rmse, aes(x = trait, y = y_pos, label = Letters, group = treatment),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Trait", y = "%RMSE", fill = "Treatment") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_discrete(labels = c("Mixed", "Dc - Ac", "De - Ac", "Dc - A50/80", "De - A50/80")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggarrange(plotlist = c(a,b), nrow = 2, ncol = 1, common.legend = TRUE)

## treatment transfer models##
# Define  transfer scenarios
transfer_results <- list()
transfer_predictions <- list()

# Define transfer scenarios
transfer_scenarios <- list(
  list(train = "control", test = "alcohol"),
  list(train = "control", test = "extra dry"),
  list(train = "alcohol", test = "control"),
  list(train = "extra dry", test = "control")
)


for(scenario_idx in seq_along(transfer_scenarios)) {
  
  scenario <- transfer_scenarios[[scenario_idx]]
  train_treatment <- scenario$train
  test_treatment <- scenario$test
  
  cat("\nScenario", scenario_idx, ": Training on", train_treatment, 
      "-> Testing on", test_treatment, "\n")
  
  # Get training and testing data based on simplified treatment groups
  if(train_treatment == "control") {
    # All non-alcohol leaves (both drying treatments)
    train_data <- dry_leaf_deriv2 %>%
      filter(alcohol_treatment == "control")
  } else if(train_treatment == "alcohol") {
    # Combined alcohol treatments (both drying treatments)
    train_data_full <- dry_leaf_deriv2 %>%
      filter(alcohol_treatment == "alcohol")
    
    # Subset by voucher to match control sample size
    control_vouchers <- dry_leaf_deriv2 %>% 
      filter(alcohol_treatment == "control") %>%
      pull(voucher_no) %>%
      unique()
    alcohol_vouchers <- train_data_full %>%
      pull(voucher_no) %>%
      unique()
    
    if(length(alcohol_vouchers) > length(control_vouchers)) {
      set.seed(123 + scenario_idx)
      selected_vouchers <- sample(alcohol_vouchers, length(control_vouchers))
      train_data <- train_data_full %>%
        filter(voucher_no %in% selected_vouchers)
      cat("  Training data subsampled:", length(alcohol_vouchers), "->", 
          length(control_vouchers), "vouchers\n")
    } else {
      train_data <- train_data_full
    }
  } else if(train_treatment == "extra dry") {
    # All extra dry leaves (all alcohol treatments)
    train_data <- dry_leaf_deriv2 %>%
      filter(dry_treatment == "extra dry")
  }
  
  if(test_treatment == "control") {
    # All non-alcohol leaves (both drying treatments)
    test_data <- dry_leaf_deriv2 %>%
      filter(alcohol_treatment == "control")
  } else if(test_treatment == "alcohol") {
    # Combined alcohol treatments (both drying treatments)
    test_data_full <- dry_leaf_deriv2 %>%
      filter(alcohol_treatment == "alcohol")
    
    # Subset by voucher to match control sample size
    control_vouchers <- dry_leaf_deriv2 %>% 
      filter(alcohol_treatment == "control") %>%
      pull(voucher_no) %>%
      unique()
    alcohol_vouchers <- test_data_full %>%
      pull(voucher_no) %>%
      unique()
    
    if(length(alcohol_vouchers) > length(control_vouchers)) {
      set.seed(456 + scenario_idx)  # Different seed for test data
      selected_vouchers <- sample(alcohol_vouchers, length(control_vouchers))
      test_data <- test_data_full %>%
        filter(voucher_no %in% selected_vouchers)
      cat("  Testing data subsampled:", length(alcohol_vouchers), "->", 
          length(control_vouchers), "vouchers\n")
    } else {
      test_data <- test_data_full
    }
  } else if(test_treatment == "extra dry") {
    # All extra dry leaves (all alcohol treatments)
    test_data <- dry_leaf_deriv2 %>%
      filter(dry_treatment == "extra dry")
  }
  
  # Check sufficient data
  if(nrow(train_data) < 10 | nrow(test_data) < 5) {
    cat("  Insufficient data, skipping...\n")
    next
  }
  
  scenario_stats <- list()
  scenario_preds <- list()
  
  for(iteration in 1:niterations) {
    
    # Sample up to 5 leaves per tree (voucher_no) for training and validation datasets
    train_data_iter <- train_data %>%
      group_by(voucher_no) %>%
      slice_sample(n = 5) %>%
      ungroup()
    
    test_data_iter <- test_data %>%
      group_by(voucher_no) %>%
      slice_sample(n = 5) %>%
      ungroup()
    
    # Prepare spectral data
    X_train <- train_data_iter %>% select(all_of(deriv_wavelength_cols))
    X_test  <- test_data_iter  %>% select(all_of(deriv_wavelength_cols))
    
    iter_stats <- data.frame()
    iter_preds <- list()
    
    for(trait in traits_list) {
      
      Y_train <- train_data_iter[[trait]]
      Y_test  <- test_data_iter[[trait]]
      
      # Build model on training treatment
      model <- plsr(Y_train ~ as.matrix(X_train),
                    ncomp = 30,
                    method = "oscorespls",
                    validation = "CV",
                    segments = 10)
      
      ncomp <- selectNcomp(model, method = "onesigma", plot = FALSE)
      
      # Test on different treatment
      test_pred <- as.vector(predict(model, newdata = as.matrix(X_test), ncomp = ncomp)[,,1])
      test_fit <- lm(Y_test ~ test_pred)
      
      # Store predictions
      iter_preds[[trait]] <- data.frame(
        iteration = iteration,
        train_treatment = train_treatment,
        test_treatment = test_treatment,
        trait = trait,
        specimen_id = test_data_iter$specimen_id,
        voucher_no = test_data_iter$voucher_no,
        measured = Y_test,
        predicted = test_pred
      )
      
      # Store statistics
      iter_stats <- rbind(iter_stats, data.frame(
        iteration = iteration,
        train_treatment = train_treatment,
        test_treatment = test_treatment,
        trait = trait,
        R2 = summary(test_fit)$r.squared,
        RMSE = RMSD(Y_test, test_pred),
        perRMSE = percentRMSD(Y_test, test_pred, 0.025, 0.975),
        bias = mean(test_pred, na.rm = TRUE) - mean(Y_test, na.rm = TRUE),
        ncomp = ncomp
      ))
    }
    
    scenario_stats[[iteration]] <- iter_stats
    scenario_preds[[iteration]] <- do.call(rbind, iter_preds)
  }
  
  transfer_results[[paste(train_treatment, "to", test_treatment)]] <- 
    do.call(rbind, scenario_stats)
  transfer_predictions[[paste(train_treatment, "to", test_treatment)]] <- 
    do.call(rbind, scenario_preds)
  
  cat("  Completed", niterations, "iterations\n")
}

# Combine transfer results
transfer_stats_df <- do.call(rbind, transfer_results)
transfer_predictions_df <- do.call(rbind, transfer_predictions)

# Summarize transfer performance
transfer_summary <- transfer_stats_df %>%
  group_by(train_treatment, test_treatment, trait) %>%
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2),
    mean_RMSE = mean(RMSE),
    sd_RMSE = sd(RMSE),
    mean_perRMSE = mean(perRMSE),
    sd_perRMSE = sd(perRMSE),
    mean_bias = mean(bias),
    sd_bias = sd(bias),
    mean_ncomp = mean(ncomp),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
transfer_summary <- transfer_summary %>%
  mutate(train_treatment = factor(train_treatment, levels = c("control", "extra dry", "alcohol")),
         test_treatment = factor(test_treatment, levels = c("control", "extra dry", "alcohol")))
transfer_summary <- transfer_summary %>%
  arrange(train_treatment, test_treatment)
#write.csv(transfer_summary, file = "TableS8.csv", row.names = F)

# Summarize predictions: mean and SD for each specimen across iterations
transfer_pred_summary <- transfer_predictions_df %>%
  group_by(train_treatment, test_treatment, trait, specimen_id, voucher_no, measured) %>%
  summarise(
    mean_predicted = mean(predicted, na.rm = TRUE),
    sd_predicted = sd(predicted, na.rm = TRUE),
    n_iterations = n(),
    .groups = "drop"
  )

# Create a plot for each transfer scenario and trait combination
transfer_plot_list <- list()

for(scenario_name in names(transfer_results)) {
  
  scenario_data <- transfer_pred_summary %>%
    filter(paste(train_treatment, "to", test_treatment) == scenario_name)
  
  if(nrow(scenario_data) == 0) next
  
  # Get training and testing treatment labels
  train_treat <- unique(scenario_data$train_treatment)
  test_treat <- unique(scenario_data$test_treatment)
  
  for(trait in traits_list) {
    
    trait_data <- scenario_data %>% filter(trait == !!trait)
    
    if(nrow(trait_data) == 0) next
    
    # Store R2 and RMSE for the plot
    plot_stats <- transfer_summary %>%
      filter(train_treatment == train_treat & test_treatment == test_treat & trait == !!trait)
    r2_val <- plot_stats$mean_R2
    per_rmse <- plot_stats$mean_perRMSE
    
    # Create plot
    p <- ggplot(trait_data, aes(x = mean_predicted, y = measured)) +
      geom_errorbarh(aes(xmin = mean_predicted - sd_predicted, 
                         xmax = mean_predicted + sd_predicted),
                     alpha = 0.5, color = "gray50") +
      geom_point(alpha = 0.5, size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                  color = "red", linewidth = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
      annotate("text", x = -Inf, y = Inf, 
               label = paste0("Train: ", train_treat, 
                              "\nTest: ", test_treat, 
                              "\nR² = ", r2_val, 
                              "\n%RMSE = ", per_rmse*100),
               hjust = 0, vjust = 1, size = 3) +
      labs(x = "", 
           y = "") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
    
    plot_name <- paste(scenario_name, trait, sep = "_")
    transfer_plot_list[[plot_name]] <- p
  }
}

transfer_plots <- ggarrange(plotlist = transfer_plot_list, ncol = 2, nrow = 4) 
annotate_figure(p = transfer_plots, left = text_grob("Measured value", rot = 90, vjust = 1, size = 14),
                bottom = text_grob("Predicted value",vjust = -0.5, size = 14),
                top = text_grob("Thickness (mm)                         LMA (g m-2)", size = 14))

# ANOVAs on transfer model results to test for scenario effect
transfer_stats_df <- transfer_stats_df %>%
  mutate(scenario = paste(train_treatment, "to", test_treatment),
         trait = recode(trait,
                        dry_thickness = "Thickness",
                        dry_LMA = "LMA"))

# Add mixed treatment
global_mix_stats_df2 <- global_mix_stats_df %>%
  mutate(scenario = "Mixed",
         train_treatment = NA,
         test_treatment = NA)
global_mix_stats_df2$treatment <- NULL

transfer_stats_df <- rbind(transfer_stats_df, global_mix_stats_df2)


# Define level order and labels
scenario_levels <- c("Mixed", "control to extra dry", "extra dry to control", 
                     "control to alcohol",   "alcohol to control")
scenario_labels <- c("Mixed", "Dc to De", "De to Dc", "Ac to A50/80", "A50/80 to Ac")

transfer_stats_df <- transfer_stats_df %>%
  mutate(scenario = factor(scenario, levels = scenario_levels))

# R2 ANOVA
transfer_r2_cld_list <- list()
for(tr in unique(transfer_stats_df$trait)) {
  trait_data <- transfer_stats_df %>% filter(trait == tr)
  
  anova_model  <- aov(R2 ~ scenario, data = trait_data)
  tukey_result <- TukeyHSD(anova_model)
  tukey_cld    <- multcompLetters4(anova_model, tukey_result)
  global_max   <- max(trait_data$R2, na.rm = TRUE)
  
  transfer_r2_cld_list[[tr]] <- data.frame(
    scenario = names(tukey_cld$scenario$Letters),
    Letters  = tukey_cld$scenario$Letters
  ) %>%
    mutate(
      trait    = tr,
      y_pos    = global_max,
      scenario = factor(scenario, levels = scenario_levels)  # <-- match factor order
    )
}
transfer_cld_r2 <- do.call(rbind, transfer_r2_cld_list)

tc_a <- ggplot(transfer_stats_df, aes(x = trait, y = R2, fill = scenario)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = transfer_cld_r2,
            aes(x = trait, y = y_pos, label = Letters, group = scenario),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "R²", fill = "Scenario") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_discrete(labels = scenario_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# %RMSE ANOVA
transfer_rmse_cld_list <- list()
for(tr in unique(transfer_stats_df$trait)) {
  trait_data <- transfer_stats_df %>% filter(trait == tr)
  
  anova_model  <- aov(perRMSE ~ scenario, data = trait_data)
  tukey_result <- TukeyHSD(anova_model)
  tukey_cld    <- multcompLetters4(anova_model, tukey_result)
  global_max   <- max(trait_data$perRMSE, na.rm = TRUE)
  
  transfer_rmse_cld_list[[tr]] <- data.frame(
    scenario = names(tukey_cld$scenario$Letters),
    Letters  = tukey_cld$scenario$Letters
  ) %>%
    mutate(
      trait    = tr,
      y_pos    = global_max,
      scenario = factor(scenario, levels = scenario_levels)  # <-- match factor order
    )
}
transfer_cld_rmse <- do.call(rbind, transfer_rmse_cld_list)

tc_b <- ggplot(transfer_stats_df, aes(x = trait, y = perRMSE, fill = scenario)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = transfer_cld_rmse,
            aes(x = trait, y = y_pos, label = Letters, group = scenario),
            position = position_dodge(width = 0.75),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Trait", y = "%RMSE", fill = "Scenario") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_discrete(labels = scenario_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggarrange(plotlist = list(tc_a, tc_b), nrow = 2, ncol = 1, common.legend = TRUE)









