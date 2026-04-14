####################################################################################################
### R script to run PLS-DA models to discriminate species from leaves subject to alcohol/drying  ###
### treatments. This script requires the completion of "Spec_prep_trait_analysis.R" to load in   ###
### required dataframes and pre-process reflectance spectra.                                     ###

# set working directory
setwd("C:/Users/rfortier/Dropbox/MBG Postdoc/Specimen prep/Data analysis")
#setwd("~/Library/CloudStorage/Dropbox/MBG Postdoc/Specimen prep/Data analysis")

# Run initial script if not already run 
#source("Spec_prep_trait_analysis.R")

# Load libraries
library(caret)
library(pls)

### PLS-DA for Species Discrimination ###

# Prepare data for PLS-DA - ensure Species is a factor and replace spaces with periods
dry_leaf_deriv <- dry_leaf_deriv %>%
  mutate(Species = as.factor(Species))
dry_leaf_deriv$Species <- gsub(" ", ".", dry_leaf_deriv$Species)

# Check species distribution
print(table(dry_leaf_deriv$Species))

# Make species list
species_list <- unique(dry_leaf_deriv$Species)

# Number of iterations for PLS-DA models
niterations_plsda <- 25

## Global PLS-DA Model ##
# Initialize storage
global_plsda_confMatrix_table <- data.frame()
global_plsda_predictions_perScan <- data.frame()
global_plsda_performanceMetrics <- data.frame()

# Enable parallel processing
library(doParallel)

# Detect available cores (leave one free)
cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

for(iteration in 1:niterations_plsda) {
  
  cat("Running global PLS-DA iteration", iteration, "\n")
  
  # Initialize empty calibration and validation datasets
  calibrationData <- data.frame()
  validationData <- data.frame()
  
  # Get list of all species
  species_list <- unique(dry_leaf_deriv$Species)
  
  # Split by species to ensure each species is represented in both training and validation
  for(sp in species_list) {
    
    # Subset to current species
    current_species_data <- dry_leaf_deriv[dry_leaf_deriv$Species == sp, ]
    
    # Get unique specimens for this species
    current_specimens <- unique(current_species_data$specimen_id)
    current_n <- length(current_specimens)
    
    # Calculate 60/40 split
    training_n <- round(current_n * 0.6)
    
    # Randomly select 60% of specimens for training
    training_specimens <- sample(current_specimens, training_n)
    validation_specimens <- setdiff(current_specimens, training_specimens)
    
    # Subset data for this species
    current_training_data <- current_species_data[current_species_data$specimen_id %in% training_specimens, ]
    current_validation_data <- current_species_data[current_species_data$specimen_id %in% validation_specimens, ]
    
    # Add to overall datasets
    calibrationData <- rbind(calibrationData, current_training_data)
    validationData <- rbind(validationData, current_validation_data)
  }
  
  # Downsample training data to balance classes
  downsampledData <- calibrationData %>%
    group_by(Species) %>%
    group_split() %>%
    map_df(~ sample_n(.x, min(table(calibrationData$Species))))
  
  # Original trainControl model
  train_ctrl <- trainControl(
    method = "repeatedcv", 
    number = 10, 
    repeats = 10, # reduce from 10 to speed things up
    sampling = "down", 
    classProbs = TRUE,
    summaryFunction = multiClassSummary,
    allowParallel = TRUE
  )
  
  # Train model on downsampled data to find optimal components
  X_down <- downsampledData %>% select(all_of(deriv_wavelength_cols))
  Y_down <- downsampledData$Species
  
  plsFit_down <- train(
    x = X_down,
    y = Y_down,
    method = "pls", 
    tuneLength = 30, 
    trControl = train_ctrl,
    metric = "Kappa",
    preProc = c("center", "scale")
  )
  
  best_ncomp <- plsFit_down$bestTune$ncomp
  
  # Upsample training data
  upsampledData <- calibrationData %>%
    group_by(Species) %>%
    group_split() %>%
    map_df(~ sample_n(.x, max(table(calibrationData$Species)), replace = TRUE))
  
  # Ensure species is a factor
  upsampledData$Species <- as.factor(upsampledData$Species)
  
  X_up <- upsampledData %>% select(all_of(deriv_wavelength_cols))
  Y_up <- upsampledData$Species
  
  # Train final model with upsampled data
  plsFit_up <- train(
    x = X_up,
    y = Y_up,
    method = "pls",
    tuneGrid = expand.grid(ncomp = 1:best_ncomp),
    trControl = train_ctrl,
    metric = "Kappa",
    preProc = c("center", "scale")
  )
  
  # Predict on validation data
  X_val <- validationData %>% select(all_of(deriv_wavelength_cols))
  Y_val <- validationData$Species
  
  predictions <- predict(plsFit_up, newdata = X_val)
  probabilities <- predict(plsFit_up, newdata = X_val, type = "prob")
  predicted_probs <- apply(probabilities, 1, max)
  
  # Ensure Y_val has same levels as predictions
  Y_val <- factor(Y_val, levels = levels(predictions))
  
  # Confusion matrix
  confMatrix <- confusionMatrix(predictions, Y_val)
  
  # Store confusion matrix table
  current_table_df <- as.data.frame(confMatrix$table) %>%
    mutate(iteration = iteration)
  global_plsda_confMatrix_table <- rbind(global_plsda_confMatrix_table, current_table_df)
  
  # Store predictions per scan
  current_predictionsPerScan <- data.frame(
    iteration = iteration,
    specimen_id = validationData$specimen_id,
    voucher_no = validationData$voucher_no,
    treatment_code = validationData$treatment_code,
    dry_treatment = validationData$dry_treatment,
    alcohol_treatment = validationData$alcohol_treatment,
    true_species = Y_val,
    predicted_species = predictions,
    predicted_prob = predicted_probs,
    correct = predictions == Y_val
  )
  global_plsda_predictions_perScan <- rbind(global_plsda_predictions_perScan, current_predictionsPerScan)
  
  # Store performance metrics
  newPerformanceRow <- data.frame(
    iteration = iteration,
    accuracy = confMatrix$overall['Accuracy'],
    kappa = confMatrix$overall['Kappa'],
    ncomp = plsFit_up$bestTune$ncomp
  )
  global_plsda_performanceMetrics <- rbind(global_plsda_performanceMetrics, newPerformanceRow)
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

# Summarize
global_plsda_summary <- global_plsda_performanceMetrics %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy   = sd(accuracy),
    mean_kappa    = mean(kappa),
    sd_kappa      = sd(kappa),
    mean_ncomp    = mean(ncomp)
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
global_plsda_summary$treatment_label <- "global"

# Calculate confusion matrix
global_confMatrix_summed <- global_plsda_confMatrix_table %>%
  group_by(Prediction, Reference) %>%
  summarise(Freq = sum(Freq), .groups = "drop") %>%
  group_by(Reference) %>%
  mutate(Percent = Freq / sum(Freq) * 100) %>%
  ungroup()

# Fix spp names
global_confMatrix_summed <- global_confMatrix_summed %>%
  mutate(Reference = as.character(Reference),
         Prediction = as.character(Prediction)) 
global_confMatrix_summed$Prediction <- gsub(".", " ", global_confMatrix_summed$Prediction, fixed = TRUE)
global_confMatrix_summed$Reference <- gsub(".", " ", global_confMatrix_summed$Reference, fixed = TRUE)

# Plot conf matrix
global_plsda_plot <- ggplot(global_confMatrix_summed, aes(x = Prediction, y = Reference, fill = Percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(round(Percent) == 0, "", round(Percent))), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Global",
       x = "",
       y = "",
       fill = "Percentage") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "italic"))  # Italicize y-axis text

#############################
### mixed-treatment ###
#############################

# Add a code for each leaf
dry_leaf_deriv$leaf_code <- paste0(dry_leaf_deriv$voucher_no, "_", dry_leaf_deriv$treatment_code, "_", dry_leaf_deriv$leaf)

global_mix_plsda_confMatrix_table    <- data.frame()
global_mix_plsda_predictions_perScan <- data.frame()
global_mix_plsda_performanceMetrics  <- data.frame()

# Initialize parallel processing
cl <- makeCluster(cores)
registerDoParallel(cl)

for(iteration in 1:niterations_plsda) {
  
  cat("Running mixed-treatment iteration", iteration, "\n")
  start_time <- Sys.time()
  
  # Sample five leaves per voucher_no (randomly each iteration)
  leaf_ids <- dry_leaf_deriv %>%
    distinct(voucher_no, leaf_code) %>%
    group_by(voucher_no) %>%
    slice_sample(n = 5) %>%
    ungroup()
  
  dry_leaf_mix <- dry_leaf_deriv %>%
    semi_join(leaf_ids, by = c("voucher_no", "leaf_code"))
  
  # Stratified 60/40 split by leaf_code within each species,
  # ensuring every species is represented in both calibration and validation
  calibrationData <- data.frame()
  validationData  <- data.frame()
  
  for(sp in unique(dry_leaf_mix$Species)) {
    sp_data      <- dry_leaf_mix[dry_leaf_mix$Species == sp, ]
    sp_specimens <- unique(sp_data$leaf_code)
    training_n   <- max(1, round(length(sp_specimens) * 0.6))
    
    # Guard: if only 1 specimen, put it in calibration only
    if(length(sp_specimens) == 1) {
      calibrationData <- rbind(calibrationData, sp_data)
      next
    }
    
    train_sp <- sample(sp_specimens, training_n)
    val_sp   <- setdiff(sp_specimens, train_sp)
    
    calibrationData <- rbind(calibrationData, sp_data[sp_data$leaf_code %in% train_sp, ])
    validationData  <- rbind(validationData,  sp_data[sp_data$leaf_code %in% val_sp, ])
  }
  
  # Downsample
  downsampledData <- calibrationData %>%
    group_by(Species) %>%
    group_split() %>%
    map_df(~ sample_n(.x, min(table(calibrationData$Species))))
  
  # Ensure species is a factor
  downsampledData$Species <- as.factor(downsampledData$Species)
  
  X_down <- downsampledData %>% select(all_of(deriv_wavelength_cols))
  Y_down <- downsampledData$Species
  
  plsFit_down <- train(
    x = X_down,
    y = Y_down,
    method = "pls",
    tuneLength = 20,
    trControl = train_ctrl,
    metric = "Kappa",
    preProc = c("center", "scale")
  )
  
  best_ncomp <- plsFit_down$bestTune$ncomp
  
  # Upsample
  upsampledData <- calibrationData %>%
    group_by(Species) %>%
    group_split() %>%
    map_df(~ sample_n(.x, max(table(calibrationData$Species)), replace = TRUE))
  
  # Ensure species is a factor
  upsampledData$Species <- as.factor(upsampledData$Species)
  
  X_up <- upsampledData %>% select(all_of(deriv_wavelength_cols))
  Y_up <- upsampledData$Species
  
  plsFit_up <- train(
    x = X_up,
    y = Y_up,
    method = "pls",
    tuneGrid = expand.grid(ncomp = 1:best_ncomp),
    trControl = train_ctrl,
    metric = "Kappa",
    preProc = c("center", "scale")
  )
  
  X_val <- validationData %>% select(all_of(deriv_wavelength_cols))
  Y_val <- factor(validationData$Species, levels = levels(Y_up))
  
  predictions     <- predict(plsFit_up, newdata = X_val)
  probabilities   <- predict(plsFit_up, newdata = X_val, type = "prob")
  predicted_probs <- apply(probabilities, 1, max)
  
  confMatrix <- confusionMatrix(predictions, Y_val)
  
  global_mix_plsda_confMatrix_table <- rbind(
    global_mix_plsda_confMatrix_table,
    as.data.frame(confMatrix$table) %>% mutate(iteration = iteration)
  )
  
  global_mix_plsda_predictions_perScan <- rbind(
    global_mix_plsda_predictions_perScan,
    data.frame(
      iteration         = iteration,
      specimen_id       = validationData$specimen_id,
      voucher_no        = validationData$voucher_no,
      treatment_code    = validationData$treatment_code,
      dry_treatment     = validationData$dry_treatment,
      alcohol_treatment = validationData$alcohol_treatment,
      true_species      = Y_val,
      predicted_species = predictions,
      predicted_prob    = predicted_probs,
      correct           = predictions == Y_val
    )
  )
  
  global_mix_plsda_performanceMetrics <- rbind(
    global_mix_plsda_performanceMetrics,
    data.frame(
      iteration = iteration,
      accuracy  = as.numeric(confMatrix$overall['Accuracy']),
      kappa     = as.numeric(confMatrix$overall['Kappa']),
      ncomp     = plsFit_up$bestTune$ncomp
    )
  )
  
  end_time <- Sys.time()
  cat("  Completed in", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

global_mix_plsda_performanceMetrics$treatment_label <- "mixed"

# Summarize
global_mix_plsda_summary <- global_mix_plsda_performanceMetrics %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy   = sd(accuracy),
    mean_kappa    = mean(kappa),
    sd_kappa      = sd(kappa),
    mean_ncomp    = mean(ncomp)
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
global_mix_plsda_summary$treatment_label <- "mixed"

global_mix_confMatrix_summed <- global_mix_plsda_confMatrix_table %>%
  group_by(Prediction, Reference) %>%
  summarise(Freq = sum(Freq), .groups = "drop") %>%
  group_by(Reference) %>%
  mutate(Percent = Freq / sum(Freq) * 100) %>%
  ungroup()

# Fix spp names
global_mix_confMatrix_summed <- global_mix_confMatrix_summed %>%
  mutate(Reference = as.character(Reference),
         Prediction = as.character(Prediction)) 
global_mix_confMatrix_summed$Prediction <- gsub(".", " ", global_mix_confMatrix_summed$Prediction, fixed = TRUE)
global_mix_confMatrix_summed$Reference <- gsub(".", " ", global_mix_confMatrix_summed$Reference, fixed = TRUE)

mix_plsda_plot <- ggplot(global_mix_confMatrix_summed, aes(x = Prediction, y = Reference, fill = Percent)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(round(Percent) == 0, "", paste0(round(Percent), "%"))), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Mixed-treatment",
       x = "",
       y = "",
       fill = "Percentage") +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "italic"))  # Italicize y-axis text

global_mix_plsda_plots <- ggarrange(plotlist = c(global_plsda_plot, mix_plsda_plot), ncol = 1, nrow = 2) 
annotate_figure(p = global_mix_plsda_plots, left = text_grob("True species", rot = 90, vjust = 1, size = 14),
                bottom = text_grob("Predicted species",vjust = -0.5, size = 14))


################################################################################
### Within-treatment PLS-DA models ###
################################################################################

# Combine the two alcohol treatments for all subsequent analyses
dry_leaf_deriv2 <- dry_leaf_deriv %>%
  mutate(alcohol_treatment = ifelse(alcohol_treatment == "control", "control", "alcohol"))

# Get unique treatment combinations
treatment_combos <- dry_leaf_deriv2 %>%
  distinct(dry_treatment, alcohol_treatment) %>%
  mutate(treatment_label = paste(dry_treatment, alcohol_treatment, sep = "_")) %>%
  mutate(dry_treatment    = as.character(dry_treatment),
         alcohol_treatment = as.character(alcohol_treatment))

treatment_plsda_confMatrix_table    <- data.frame()
treatment_plsda_predictions_perScan <- data.frame()
treatment_plsda_performanceMetrics  <- data.frame()

# Initialize parallel processing
cl <- makeCluster(cores)
registerDoParallel(cl)

for(t in 1:nrow(treatment_combos)) {
  
  current_dry    <- treatment_combos$dry_treatment[t]
  current_alc    <- treatment_combos$alcohol_treatment[t]
  treatment_label <- treatment_combos$treatment_label[t]
  
  cat("\nProcessing treatment:", treatment_label, "\n")
  
  # Filter to current treatment
  treatment_data_full <- dry_leaf_deriv2 %>%
    filter(dry_treatment == current_dry & alcohol_treatment == current_alc)
  
  # If alcohol treatment, subsample to match control sample size:
  # sample 50% of image_names per species (matching PLSR within-treatment logic)
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
      
      selected_colls <- treatment_data_full %>%
        distinct(Species, image_name) %>%
        group_by(Species) %>%
        slice_sample(prop = 0.5) %>%
        ungroup() %>%
        pull(image_name)
      
      treatment_data <- treatment_data_full %>%
        filter(image_name %in% selected_colls)
      
      cat("  Subsampled alcohol data:", length(alcohol_colls), "->",
          length(selected_colls), "collections\n")
    } else {
      treatment_data <- treatment_data_full
    }
  } else {
    treatment_data <- treatment_data_full
  }
  
  # Create leaf_code once before the iteration loop
  treatment_data <- treatment_data %>%
    mutate(leaf_code = paste0(specimen_id, "_", leaf))
  
  for(iteration in 1:niterations_plsda) {
    
    # Split 60% train / 40% validation by leaf_code within specimen_id
    training_leaves <- treatment_data %>%
      distinct(specimen_id, leaf_code) %>%
      group_by(specimen_id) %>%
      slice_sample(prop = 0.6) %>%
      ungroup() %>%
      pull(leaf_code)
    
    calibrationData <- treatment_data %>% filter(leaf_code %in% training_leaves)
    validationData  <- treatment_data %>% filter(!leaf_code %in% training_leaves)
    
    # Downsample
    downsampledData <- calibrationData %>%
      group_by(Species) %>%
      group_split() %>%
      map_df(~ sample_n(.x, min(table(calibrationData$Species))))
    
    # Ensure species is a factor
    downsampledData$Species <- as.factor(downsampledData$Species)
    
    X_down <- downsampledData %>% select(all_of(deriv_wavelength_cols))
    Y_down <- downsampledData$Species
    
    plsFit_down <- train(
      x = X_down,
      y = Y_down,
      method = "pls",
      tuneLength = min(20, nrow(downsampledData) - 1),
      trControl = train_ctrl,
      metric = "Kappa",
      preProc = c("center", "scale")
    )
    
    best_ncomp <- plsFit_down$bestTune$ncomp
    
    # Upsample
    upsampledData <- calibrationData %>%
      group_by(Species) %>%
      group_split() %>%
      map_df(~ sample_n(.x, max(table(calibrationData$Species)), replace = TRUE))
    # Ensure species is a factor
    upsampledData$Species <- as.factor(upsampledData$Species)
    
    X_up <- upsampledData %>% select(all_of(deriv_wavelength_cols))
    Y_up <- upsampledData$Species
    
    plsFit_up <- train(
      x = X_up,
      y = Y_up,
      method = "pls",
      tuneGrid = expand.grid(ncomp = 1:best_ncomp),
      trControl = train_ctrl,
      metric = "Kappa",
      preProc = c("center", "scale")
    )
    
    X_val <- validationData %>% select(all_of(deriv_wavelength_cols))
    Y_val <- factor(validationData$Species, levels = levels(Y_up))
    
    predictions     <- predict(plsFit_up, newdata = X_val)
    probabilities   <- predict(plsFit_up, newdata = X_val, type = "prob")
    predicted_probs <- apply(probabilities, 1, max)
    
    confMatrix <- confusionMatrix(predictions, Y_val)
    
    treatment_plsda_confMatrix_table <- rbind(
      treatment_plsda_confMatrix_table,
      as.data.frame(confMatrix$table) %>%
        mutate(iteration = iteration, treatment_label = treatment_label)
    )
    
    treatment_plsda_predictions_perScan <- rbind(
      treatment_plsda_predictions_perScan,
      data.frame(
        iteration         = iteration,
        treatment_label   = treatment_label,
        specimen_id       = validationData$specimen_id,
        voucher_no        = validationData$voucher_no,
        true_species      = Y_val,
        predicted_species = predictions,
        predicted_prob    = predicted_probs,
        correct           = predictions == Y_val
      )
    )
    
    treatment_plsda_performanceMetrics <- rbind(
      treatment_plsda_performanceMetrics,
      data.frame(
        iteration       = iteration,
        treatment_label = treatment_label,
        accuracy        = as.numeric(confMatrix$overall['Accuracy']),
        kappa           = as.numeric(confMatrix$overall['Kappa']),
        ncomp           = plsFit_up$bestTune$ncomp
      )
    )
  }
  
  cat("  Completed", niterations_plsda, "iterations for", treatment_label, "\n")
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

# Summarize treatment-specific performance
treatment_plsda_summary <- treatment_plsda_performanceMetrics %>%
  group_by(treatment_label) %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy   = sd(accuracy),
    mean_kappa    = mean(kappa),
    sd_kappa      = sd(kappa),
    mean_ncomp    = mean(ncomp),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

# Combine with mixed-treatment
treatment_plsda_performanceMetrics <- rbind(treatment_plsda_performanceMetrics, global_mix_plsda_performanceMetrics)
treatment_plsda_performanceMetrics <- treatment_plsda_performanceMetrics %>%
  mutate(treatment_label = factor(treatment_label, levels = c("mixed", "control_control", "extra dry_control", "control_alcohol", "extra dry_alcohol")))

treatment_plsda_summary <- rbind(treatment_plsda_summary, global_mix_plsda_summary)
treatment_plsda_summary <- treatment_plsda_summary %>%
  mutate(treatment_label = factor(treatment_label, levels = c("mixed", "control_control", "extra dry_control", "control_alcohol", "extra dry_alcohol")))

# Treatment-specific confusion matrices
treatment_confmat_plots <- list()

for(trt in c("control_control", "control_alcohol", "extra dry_control", "extra dry_alcohol")) {
  
  confmat_data <- treatment_plsda_confMatrix_table %>%
    filter(treatment_label == trt) %>%
    group_by(Prediction, Reference) %>%
    summarise(Freq = sum(Freq), .groups = "drop") %>%
    group_by(Reference) %>%
    mutate(Percent = Freq / sum(Freq) * 100) %>%
    ungroup()
  
  p <- ggplot(confmat_data, aes(x = gsub("."," ", Prediction, fixed = TRUE), y = gsub(".", " ", Reference, fixed = TRUE), fill = Percent)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(round(Percent) == 0, "", paste0(round(Percent), "%"))), size = 2.5) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = trt, x = "", y = "") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          axis.text.y = element_text(face = "italic"),
          legend.position = "none")
  
  treatment_confmat_plots[[trt]] <- p
}

treatment_confmat_combined <- ggarrange(plotlist = treatment_confmat_plots, ncol = 2, nrow = 2)
annotate_figure(treatment_confmat_combined,
                left = text_grob("True",face = "bold", size = 14, rot = 90),
                bottom = text_grob("Prediction", face = "bold", size = 14))

# ANOVAs on within-treatment model results to test for treatment effect
acc_cld_list <- list()

for(met in "accuracy") {
  anova_model  <- aov(accuracy ~ treatment_label, data = treatment_plsda_performanceMetrics)
  tukey_result <- TukeyHSD(anova_model)
  tukey_cld    <- multcompLetters4(anova_model, tukey_result)
  
  global_max <- max(treatment_plsda_performanceMetrics$accuracy, na.rm = TRUE)
  
  acc_cld_list[["accuracy"]] <- data.frame(
    treatment_label = names(tukey_cld$treatment_label$Letters),
    Letters         = tukey_cld$treatment_label$Letters,
    y_pos           = global_max
  )
}

cld_acc <- acc_cld_list[["accuracy"]]
cld_acc <- cld_acc %>%
  mutate(treatment_label = factor(treatment_label, levels = c("mixed", "control_control", "extra dry_control", "control_alcohol", "extra dry_alcohol")))


kap_cld_list <- list()

for(met in "kappa") {
  anova_model  <- aov(kappa ~ treatment_label, data = treatment_plsda_performanceMetrics)
  tukey_result <- TukeyHSD(anova_model)
  tukey_cld    <- multcompLetters4(anova_model, tukey_result)
  
  global_max <- max(treatment_plsda_performanceMetrics$kappa, na.rm = TRUE)
  
  kap_cld_list[["kappa"]] <- data.frame(
    treatment_label = names(tukey_cld$treatment_label$Letters),
    Letters         = tukey_cld$treatment_label$Letters,
    y_pos           = global_max
  )
}

cld_kap <- kap_cld_list[["kappa"]]
cld_kap <- cld_kap %>%
  mutate(treatment_label = factor(treatment_label, levels = c("mixed", "control_control", "extra dry_control", "control_alcohol", "extra dry_alcohol")))

ta <- ggplot(treatment_plsda_performanceMetrics,
             aes(x = treatment_label, y = accuracy, fill = treatment_label)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = cld_acc,
            aes(x = treatment_label, y = y_pos, label = Letters),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "Accuracy", fill = "Treatment") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(labels = c("Mixed", "Dc - Ac", "De - Ac", "Dc - A50/80", "De - A50/80")) +
  #scale_fill_discrete(labels = c("Mixed", "Dc - Ac", "De - Ac", "Dc - A50/80", "De - A50/80")) +
  scale_fill_discrete(guide = "none") +
  theme(axis.text.x = element_text(),
        legend.position = "bottom")

tb <- ggplot(treatment_plsda_performanceMetrics,
             aes(x = treatment_label, y = kappa, fill = treatment_label)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = cld_kap,
            aes(x = treatment_label, y = y_pos, label = Letters),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Treatment", y = "Kappa", fill = "Treatment") +
  scale_x_discrete(labels = c("Mixed", "Dc - Ac", "De - Ac", "Dc - A50/80", "De - A50/80")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_discrete(guide = "none") +
  theme(axis.text.x = element_text(),
        legend.position = "bottom")

ggarrange(plotlist = list(ta, tb), nrow = 2, ncol = 1, common.legend = TRUE)

# Make supplemental table
treatment_plsda_summary <- rbind(treatment_plsda_summary, global_plsda_summary)
#write.csv(treatment_plsda_summary, file = "TableS9.csv", row.names = F)

################################################################################
### Treatment transfer PLS-DA models ###
################################################################################

transfer_plsda_confMatrix_table    <- data.frame()
transfer_plsda_predictions_perScan <- data.frame()
transfer_plsda_performanceMetrics  <- data.frame()

# Define transfer scenarios (matching PLSR)
transfer_scenarios <- list(
  list(train = "control",   test = "alcohol"),
  list(train = "control",   test = "extra dry"),
  list(train = "alcohol",   test = "control"),
  list(train = "extra dry", test = "control")
)

# Initialize parallel processing
cl <- makeCluster(cores)
registerDoParallel(cl)

for(scenario_idx in seq_along(transfer_scenarios)) {
  
  train_treatment <- transfer_scenarios[[scenario_idx]]$train
  test_treatment  <- transfer_scenarios[[scenario_idx]]$test
  
  cat("\nScenario", scenario_idx, ": Training on", train_treatment,
      "-> Testing on", test_treatment, "\n")
  
  # Determine which experimental axis is being used
  using_dry_axis <- (train_treatment %in% c("control", "extra dry") &
                       test_treatment %in% c("control", "extra dry"))
  
  # Get training data
  if (using_dry_axis) {
    if (train_treatment == "control") {
      train_data <- dry_leaf_deriv2 %>%
        filter(dry_treatment == "control")
    } else if (train_treatment == "extra dry") {
      train_data <- dry_leaf_deriv2 %>%
        filter(dry_treatment == "extra dry")
    }
  } else {
    # Alcohol axis: control vs alcohol
    if (train_treatment == "control") {
      train_data <- dry_leaf_deriv2 %>%
        filter(alcohol_treatment == "control")
    } else if (train_treatment == "alcohol") {
      train_data_full <- dry_leaf_deriv2 %>%
        filter(alcohol_treatment == "alcohol")
      
      control_vouchers <- dry_leaf_deriv2 %>%
        filter(alcohol_treatment == "control") %>%
        pull(voucher_no) %>% unique()
      alcohol_vouchers <- train_data_full %>%
        pull(voucher_no) %>% unique()
      
      if (length(alcohol_vouchers) > length(control_vouchers)) {
        set.seed(123 + scenario_idx)
        selected_vouchers <- sample(alcohol_vouchers, length(control_vouchers))
        train_data <- train_data_full %>% filter(voucher_no %in% selected_vouchers)
        cat("  Training data subsampled:", length(alcohol_vouchers), "->",
            length(control_vouchers), "vouchers\n")
      } else {
        train_data <- train_data_full
      }
    }
  }
  
  # Get testing data
  if (using_dry_axis) {
    if (test_treatment == "control") {
      test_data <- dry_leaf_deriv2 %>%
        filter(dry_treatment == "control")
    } else if (test_treatment == "extra dry") {
      test_data <- dry_leaf_deriv2 %>%
        filter(dry_treatment == "extra dry")
    }
  } else {
    # Alcohol axis: control vs alcohol
    if (test_treatment == "control") {
      test_data <- dry_leaf_deriv2 %>%
        filter(alcohol_treatment == "control")
    } else if (test_treatment == "alcohol") {
      test_data_full <- dry_leaf_deriv2 %>%
        filter(alcohol_treatment == "alcohol")
      
      control_vouchers <- dry_leaf_deriv2 %>%
        filter(alcohol_treatment == "control") %>%
        pull(voucher_no) %>% unique()
      alcohol_vouchers <- test_data_full %>%
        pull(voucher_no) %>% unique()
      
      if (length(alcohol_vouchers) > length(control_vouchers)) {
        set.seed(456 + scenario_idx)
        selected_vouchers <- sample(alcohol_vouchers, length(control_vouchers))
        test_data <- test_data_full %>% filter(voucher_no %in% selected_vouchers)
        cat("  Testing data subsampled:", length(alcohol_vouchers), "->",
            length(control_vouchers), "vouchers\n")
      } else {
        test_data <- test_data_full
      }
    }
  }
  cat("  Training samples:", nrow(train_data), "| Testing samples:", nrow(test_data), "\n")
  
  for(iteration in 1:niterations_plsda) {
    
    # Sample up to 5 leaves per tree for training and test sets each iteration
    train_data_iter <- train_data %>%
      group_by(voucher_no) %>%
      slice_sample(n = 5) %>%
      ungroup()
    
    test_data_iter <- test_data %>%
      group_by(voucher_no) %>%
      slice_sample(n = 5) %>%
      ungroup()
    
    # Make a leaf code
    train_data_iter$leaf_code <- paste0(train_data_iter$specimen_id,"_", train_data_iter$leaf)
    test_data_iter$leaf_code <- paste0(test_data_iter$specimen_id,"_", test_data_iter$leaf)
    
    # Stratified 60/40 split by leaf_code within each species,
    # ensuring every species is represented in both calibration and validation
    calibrationData <- data.frame()
    validationData  <- data.frame()
    
    for(sp in unique(train_data_iter$Species)) {
      train_sp_data <- train_data_iter[train_data_iter$Species == sp, ]
      test_sp_data <- test_data_iter[test_data_iter$Species == sp, ]
      train_sp_specimens <- unique(train_sp_data$leaf_code)
      test_sp_specimens <- unique(test_sp_data$leaf_code)
      
      training_n   <- max(1, round(length(train_sp_specimens) * 0.6))
      
      train_sp <- sample(train_sp_specimens, training_n)
      val_sp   <- setdiff(test_sp_specimens, train_sp)
      
      calibrationData <- rbind(calibrationData, train_sp_data[train_sp_data$leaf_code %in% train_sp, ])
      validationData  <- rbind(validationData,  test_sp_data[test_sp_data$leaf_code %in% val_sp, ])
    }
    
    # Downsample
    downsampledData <- calibrationData %>%
      group_by(Species) %>%
      group_split() %>%
      map_df(~ sample_n(.x, min(table(calibrationData$Species))))
    
    # Ensure species is a factor
    downsampledData$Species <- as.factor(downsampledData$Species)
    
    X_down <- downsampledData %>% select(all_of(deriv_wavelength_cols))
    Y_down <- downsampledData$Species
    
    plsFit_down <- train(
      x = X_down,
      y = Y_down,
      method = "pls",
      tuneLength = min(20, nrow(downsampledData) - 1),
      trControl = train_ctrl,
      metric = "Kappa",
      preProc = c("center", "scale")
    )
    
    best_ncomp <- plsFit_down$bestTune$ncomp
    
    # Upsample
    upsampledData <- calibrationData %>%
      group_by(Species) %>%
      group_split() %>%
      map_df(~ sample_n(.x, max(table(calibrationData$Species)), replace = TRUE))
    
    # Ensure species is a factor
    upsampledData$Species <- as.factor(upsampledData$Species)
    
    X_up <- upsampledData %>% select(all_of(deriv_wavelength_cols))
    Y_up <- upsampledData$Species
    
    plsFit_up <- train(
      x = X_up,
      y = Y_up,
      method = "pls",
      tuneGrid = expand.grid(ncomp = 1:best_ncomp),
      trControl = train_ctrl,
      metric = "Kappa",
      preProc = c("center", "scale")
    )
    
    X_val <- validationData %>% select(all_of(deriv_wavelength_cols))
    Y_val <- factor(validationData$Species, levels = levels(Y_up))
    
    predictions     <- predict(plsFit_up, newdata = X_val)
    probabilities   <- predict(plsFit_up, newdata = X_val, type = "prob")
    predicted_probs <- apply(probabilities, 1, max)
    
    confMatrix <- confusionMatrix(predictions, Y_val)
    
    transfer_plsda_confMatrix_table <- rbind(
      transfer_plsda_confMatrix_table,
      as.data.frame(confMatrix$table) %>%
        mutate(iteration       = iteration,
               train_treatment = train_treatment,
               test_treatment  = test_treatment)
    )
    
    transfer_plsda_predictions_perScan <- rbind(
      transfer_plsda_predictions_perScan,
      data.frame(
        iteration         = iteration,
        train_treatment   = train_treatment,
        test_treatment    = test_treatment,
        specimen_id       = test_data_iter$specimen_id,
        voucher_no        = test_data_iter$voucher_no,
        true_species      = Y_val,
        predicted_species = predictions,
        predicted_prob    = predicted_probs,
        correct           = predictions == Y_val
      )
    )
    
    transfer_plsda_performanceMetrics <- rbind(
      transfer_plsda_performanceMetrics,
      data.frame(
        iteration       = iteration,
        train_treatment = train_treatment,
        test_treatment  = test_treatment,
        accuracy        = as.numeric(confMatrix$overall['Accuracy']),
        kappa           = as.numeric(confMatrix$overall['Kappa']),
        ncomp           = plsFit_up$bestTune$ncomp
      )
    )
  }
  
  cat("  Completed", niterations_plsda, "iterations\n")
}

# Stop parallel processing
stopCluster(cl)
registerDoSEQ()

# Summarize transfer performance
transfer_plsda_summary <- transfer_plsda_performanceMetrics %>%
  group_by(train_treatment, test_treatment) %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy   = sd(accuracy),
    mean_kappa    = mean(kappa),
    sd_kappa      = sd(kappa),
    mean_ncomp    = mean(ncomp),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))
transfer_plsda_summary <- transfer_plsda_summary %>%
  mutate(train_treatment = factor(train_treatment, levels = c("control", "extra dry", "alcohol")),
         test_treatment = factor(test_treatment, levels = c("control", "extra dry", "alcohol")))
transfer_plsda_summary <- transfer_plsda_summary %>%
  arrange(train_treatment, test_treatment)

# Make supp table
#write.csv(transfer_plsda_summary, file = "TableS9.csv", row.names = F)

# Transfer confusion matrices
transfer_confmat_plots <- list()

for(scenario in unique(paste(transfer_plsda_confMatrix_table$train_treatment,
                             "to",
                             transfer_plsda_confMatrix_table$test_treatment))) {
  
  confmat_data <- transfer_plsda_confMatrix_table %>%
    filter(paste(train_treatment, "to", test_treatment) == scenario) %>%
    group_by(Prediction, Reference) %>%
    summarise(Freq = sum(Freq), .groups = "drop") %>%
    group_by(Reference) %>%
    mutate(Percent = Freq / sum(Freq) * 100) %>%
    ungroup()
  
  if(nrow(confmat_data) == 0) next
  
  p <- ggplot(confmat_data, aes(x = gsub(".", " ", Prediction, fixed = TRUE), y = gsub(".", " ", Reference, fixed = TRUE), fill = Percent)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(round(Percent) == 0, "", paste0(round(Percent), "%"))), size = 2) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = scenario, x = "", y = "") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
          axis.text.y = element_text(face = "italic"),
          legend.position = "none")
  
  transfer_confmat_plots[[scenario]] <- p
}


transfer_confmat_combined <- ggarrange(plotlist = transfer_confmat_plots, ncol = 2, nrow = 2)
annotate_figure(transfer_confmat_combined, left = text_grob("True species",face = "bold", size = 14, rot = 90),
                bottom = text_grob("Predicted species", face = "bold", size = 14))

# ANOVAs on transfer model results to test for scenario effect
transfer_plsda_performanceMetrics <- transfer_plsda_performanceMetrics %>%
  mutate(scenario = paste(train_treatment, "to", test_treatment))

# Add mixed treatment
global_mix_plsda_performanceMetrics2 <- global_mix_plsda_performanceMetrics %>%
  mutate(scenario = "Mixed",
         train_treatment = NA,
         test_treatment  = NA,
         treatment_label = NULL)

transfer_plsda_performanceMetrics <- rbind(transfer_plsda_performanceMetrics, global_mix_plsda_performanceMetrics2)

# Define level order and labels
scenario_levels <- c("Mixed", "control to extra dry", "extra dry to control", 
                     "control to alcohol",   "alcohol to control")
scenario_labels <- c("Mixed", "Dc to De", "De to Dc", "Ac to A50/80", "A50/80 to Ac")

transfer_plsda_performanceMetrics <- transfer_plsda_performanceMetrics %>%
  mutate(scenario = factor(scenario, levels = scenario_levels))

transfer_acc_cld_list <- list()

anova_model  <- aov(accuracy ~ scenario, data = transfer_plsda_performanceMetrics)
tukey_result <- TukeyHSD(anova_model)
tukey_cld    <- multcompLetters4(anova_model, tukey_result)
global_max   <- max(transfer_plsda_performanceMetrics$accuracy, na.rm = TRUE)

transfer_acc_cld_list[["accuracy"]] <- data.frame(
  scenario = names(tukey_cld$scenario$Letters),
  Letters  = tukey_cld$scenario$Letters,
  y_pos    = global_max
)

transfer_cld_acc <- transfer_acc_cld_list[["accuracy"]]

tc_a <- ggplot(transfer_plsda_performanceMetrics,
               aes(x = scenario, y = accuracy, fill = scenario)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = transfer_cld_acc,
            aes(x = scenario, y = y_pos, label = Letters),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "Accuracy", fill = "Scenario") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(labels = scenario_labels) +
  scale_fill_discrete(guide = "none") +
  theme(axis.text.x = element_text(),
        legend.position = "bottom")

transfer_kap_cld_list <- list()

anova_model  <- aov(kappa ~ scenario, data = transfer_plsda_performanceMetrics)
tukey_result <- TukeyHSD(anova_model)
tukey_cld    <- multcompLetters4(anova_model, tukey_result)
global_max   <- max(transfer_plsda_performanceMetrics$kappa, na.rm = TRUE)

transfer_kap_cld_list[["kappa"]] <- data.frame(
  scenario = names(tukey_cld$scenario$Letters),
  Letters  = tukey_cld$scenario$Letters,
  y_pos    = global_max
)

transfer_cld_kap <- transfer_kap_cld_list[["kappa"]]

tc_b <- ggplot(transfer_plsda_performanceMetrics,
               aes(x = scenario, y = kappa, fill = scenario)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(data = transfer_cld_kap,
            aes(x = scenario, y = y_pos, label = Letters),
            vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  labs(x = "Scenario", y = "Kappa", fill = "Scenario") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(labels = scenario_labels) +
  scale_fill_discrete(guide = "none") +
  theme(axis.text.x = element_text(),
        legend.position = "bottom")

ggarrange(plotlist = list(tc_a, tc_b), nrow = 2, ncol = 1, common.legend = TRUE)


# Model comparison figure
model_summaries_plot <- treatment_plsda_summary %>%
  filter(treatment_label != "mixed" & treatment_label != "global") %>%
  separate_wider_delim(treatment_label, delim = "_", names = c("dry_treatment", "alcohol_treatment"))

library(grid) 

model_summaries_plot <- model_summaries_plot %>%
  mutate(
    dry_treatment = factor(dry_treatment, levels = c("extra dry", "control")),
    alcohol_treatment = factor(alcohol_treatment, levels = c("control", "alcohol")),
    x = as.numeric(alcohol_treatment),
    y = as.numeric(dry_treatment)
  )

# Make arrows
arrows_df <- tibble::tibble(
  x     = c(1.2, 1.8, 0.2, 0.1),
  xend  = c(1.8, 1.2, 0.2, 0.1),
  y     = c(2.7, 2.8, 1.2, 1.8),
  yend  = c(2.7, 2.8, 1.8, 1.2),
  mean_accuracy = c(
    transfer_plsda_summary$mean_accuracy[2], # control to alcohol
    transfer_plsda_summary$mean_accuracy[4], # alcohol to control
    transfer_plsda_summary$mean_accuracy[1], # control to extra dry
    transfer_plsda_summary$mean_accuracy[3]  # extra dry to control
  )
)

# Choose gradient colors
cols <- c("#EF7C12FF", "#E9E29CFF", "#89d06c")
acc_breaks  <- c(0,0.5,1)

# Offset for the labels so they don't sit on top of the arrowhead
x_offset <- 0.05
y_offset <- 0.05

# Read in PLSR model stats 
treatment_plsr_summary <- read.csv("TableS5b.csv")
treatment_plsr_summary <- treatment_plsr_summary %>%
  filter(treatment != "Mixed") %>%
  separate_wider_delim(treatment, delim = "-", names = c("dry_treatment", "alcohol_treatment")) %>%
  mutate(dry_treatment = tolower(dry_treatment),
         alcohol_treatment = tolower(alcohol_treatment),
         dry_treatment = factor(dry_treatment, levels = c("extra dry", "control")),
         alcohol_treatment = factor(alcohol_treatment, levels = c("control", "alcohol")),
         x = as.numeric(alcohol_treatment),
         y = as.numeric(dry_treatment))
treatment_plsr_thick <- treatment_plsr_summary %>%
  filter(trait == "Thickness")
treatment_plsr_lma <- treatment_plsr_summary %>%
  filter(trait == "LMA")

transfer_plsr_summary <- read.csv("TableS8.csv")
transfer_plsr_thick <- transfer_plsr_summary %>%
  filter(trait == "dry_thickness")
transfer_plsr_lma <- transfer_plsr_summary %>%
  filter(trait == "dry_LMA")

# Make arrows
arrows_df_thick <- tibble::tibble(
  x     = c(1.2, 1.8, 0.2, 0.1),
  xend  = c(1.8, 1.2, 0.2, 0.1),
  y     = c(2.7, 2.8, 1.2, 1.8),
  yend  = c(2.7, 2.8, 1.8, 1.2),
  mean_R2 = c(
    transfer_plsr_thick$mean_R2[2], # control to alcohol
    transfer_plsr_thick$mean_R2[4], # alcohol to control
    transfer_plsr_thick$mean_R2[1], # control to extra dry
    transfer_plsr_thick$mean_R2[3]  # extra dry to control
  )
)

arrows_df_lma <- tibble::tibble(
  x     = c(1.2, 1.8, 0.2, 0.1),
  xend  = c(1.8, 1.2, 0.2, 0.1),
  y     = c(2.7, 2.8, 1.2, 1.8),
  yend  = c(2.7, 2.8, 1.8, 1.2),
  mean_R2 = c(
    transfer_plsr_lma$mean_R2[2], # control to alcohol
    transfer_plsr_lma$mean_R2[4], # alcohol to control
    transfer_plsr_lma$mean_R2[1], # control to extra dry
    transfer_plsr_lma$mean_R2[3]  # extra dry to control
  )
)

# Plot PLS-DA summary
a <- ggplot(model_summaries_plot, aes(x = x, y = y, fill = mean_accuracy)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_accuracy, 2)), size = 4) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("AC", "A50/80"),
    position = "top") +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("DE", "DC")) +
  scale_fill_gradientn(
    name = "Model performance",
    colors = cols,
    values = acc_breaks,
    limits = c(0,1)) +
  scale_color_gradientn(
    name = "Model performance",
    colors = cols,
    values = acc_breaks,
    limits = c(0,1)) +
  coord_cartesian(clip = "off") +
  # draw arrows colored by mean_accuracy
  geom_segment(
    data = arrows_df,
    aes(x = x, y = y, xend = xend, yend = yend, color = mean_accuracy),
    arrow = arrow(length = unit(0.3, "cm")),
    linewidth = 1.3
  ) +
  geom_text(
    data = arrows_df,
    aes(x = ifelse(y > 2, (x + xend) / 2, ifelse(x == 0.2, x + 0.15, x - 0.15)),
        y = ifelse(x < 1, (y + yend) / 2, ifelse(y == 2.7, y - 0.1, y + 0.1)),
        label = round(mean_accuracy, 2)),
    size = 4,
    color = "black") +
  labs(title = "C) Species classification accuracy",
       x = "",
       y = "",
       fill = "Mean Accuracy",
       color = "Transfer Accuracy") +
  theme_void() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        axis.text.x  = element_text(size = 16, vjust = -12),
        axis.text.y = element_text(size = 16, hjust = 7))

# Plot thickness
b <- ggplot(treatment_plsr_thick, aes(x = x, y = y, fill = mean_R2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_R2, 2)), size = 4) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("AC", "A50/80"),
    position = "top") +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("DE", "DC")) +
  scale_fill_gradientn(
    name = "Model performance",
    colors = cols,
    values = acc_breaks,
    limits = c(0,1)) +
  scale_color_gradientn(
    name = "Model performance",
    colors = cols,
    values = acc_breaks,
    limits = c(0,1)) +
  coord_cartesian(clip = "off") +
  # draw arrows colored by mean_R2
  geom_segment(
    data = arrows_df_thick,
    aes(x = x, y = y, xend = xend, yend = yend, color = mean_R2),
    arrow = arrow(length = unit(0.3, "cm")),
    linewidth = 1.3,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = arrows_df_thick,
    aes(x = ifelse(y > 2, (x + xend) / 2, ifelse(x == 0.2, x + 0.15, x - 0.15)),
        y = ifelse(x < 1, (y + yend) / 2, ifelse(y == 2.7, y - 0.1, y + 0.1)),
        label = round(mean_R2, 2)),
    size = 4,
    color = "black") +
  labs(title = "A) Leaf thickness prediction R2",
       x = "",
       y = "",
       fill = "Mean R2",
       color = "Transfer R2") +
  theme_void() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        axis.text.x  = element_text(size = 16, vjust = -12),
        axis.text.y = element_text(size = 16, hjust = 7))

# Plot lma
c <- ggplot(treatment_plsr_lma, aes(x = x, y = y, fill = mean_R2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_R2, 2)), size = 4) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("AC", "A50/80"),
    position = "top") +
  scale_y_continuous(
    breaks = c(1, 2),
    labels = c("DE", "DC")) +
  scale_fill_gradientn(
    name = "Model performance",
    colors = cols,
    values = acc_breaks,
    limits = c(0,1)) +
  scale_color_gradientn(
    name = "Model performance",
    colors = cols,
    values = acc_breaks,
    limits = c(0,1)) +
  coord_cartesian(clip = "off") +
  # draw arrows colored by mean_R2
  geom_segment(
    data = arrows_df_lma,
    aes(x = x, y = y, xend = xend, yend = yend, color = mean_R2),
    arrow = arrow(length = unit(0.3, "cm")),
    linewidth = 1.3
  ) +
  geom_text(
    data = arrows_df_lma,
    aes(x = ifelse(y > 2, (x + xend) / 2, ifelse(x == 0.2, x + 0.15, x - 0.15)),
        y = ifelse(x < 1, (y + yend) / 2, ifelse(y == 2.7, y - 0.1, y + 0.1)),
        label = round(mean_R2, 2)),
    size = 4,
    color = "black") +
  labs(title = "B) LMA prediction R2",
    x = "",
       y = "",
       fill = "Mean R2",
       color = "Transfer R2") +
  theme_void() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(1, "cm"),
        axis.text.x  = element_text(size = 16, vjust = -12),
        axis.text.y = element_text(size = 16, hjust = 7))


summary_figure <- ggarrange(b, c, a, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
annotate_figure(summary_figure)







