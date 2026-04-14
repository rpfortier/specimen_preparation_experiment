###################################################################################
### R script to clean leaf trait data and spectra, and analyze changes in both  ### 
### due to alcohol/drying treatments. 

# set working directory
setwd("C:/Users/rfortier/Dropbox/MBG Postdoc/Specimen prep/Data analysis")
#setwd("~/Library/CloudStorage/Dropbox/MBG Postdoc/Specimen prep/Data analysis")

# Load libraries
library(readxl)
library(ggplot2)
library(ggspectra)
library(grDevices)
library(ggpubr)
library(ggpattern)
library(tibble)
library(purrr)
library(yardstick)
library(data.table)
library(emmeans)
library(tidyverse)
library(lme4)
library(lmerTest)
library(MuMIn)
library(multcompView)
library(broom)
library(paletteer)
library(vegan)

# Read in some functions for later
source("spectral_functions.R")

# Read in leaf trait data
leaf_traits <- read.csv("specimen_prep_leaf_data.csv")

# Change voucher_no to character
leaf_traits$voucher_no <- as.factor(leaf_traits$voucher_no)

# Convert mass to grams
leaf_traits$mass <- leaf_traits$mass/1000 #convert from mg to g

# Calculate leaf mass per area (LMA)
leaf_traits$fresh_LMA <- (leaf_traits$mass)/leaf_traits$fresh_area*1000 
leaf_traits$dry_LMA <- (leaf_traits$mass)/leaf_traits$dry_area*1000 

# Calculate other traits and proportional change
leaf_traits <- leaf_traits %>%
  mutate(dry_treatment = factor(dry_treatment),
         alcohol_treatment = factor(alcohol_treatment),
         treatment_code = gsub("-", "_", treatment_code),
         fresh_thickness = (rowMeans(select(., fresh_thickness_1, fresh_thickness_2, fresh_thickness_3))) / 1000,
         dry_thickness = (rowMeans(select(., dry_thickness_1, dry_thickness_2, dry_thickness_3))) / 1000,
         fresh_ELW = fresh_circle_perimeter/pi,
         dry_ELW = dry_circle_perimeter/pi,
         prop_thickness = dry_thickness/fresh_thickness,
         prop_ELW = dry_ELW/fresh_ELW,
         prop_area = dry_area/fresh_area,
         prop_length = dry_length/fresh_length,
         prop_LMA = dry_LMA/fresh_LMA)

# Set factor levels in the desired order
leaf_traits$dry_treatment <- factor(
  leaf_traits$dry_treatment,
  levels = c("control", "extra dry")
)
leaf_traits$alcohol_treatment <- factor(
  leaf_traits$alcohol_treatment,
  levels = c("control", "alcohol", "extra alcohol")
)

# Summarize the fresh trait values
leaf_traits_long <- leaf_traits %>%
  pivot_longer(
    cols = c("fresh_thickness", "fresh_LMA", "fresh_ELW", "fresh_length", "fresh_area",
             "dry_thickness", "dry_LMA", "dry_ELW", "dry_length", "dry_area"),
    names_to = "trait",
    values_to = "value"
  ) %>%
  separate(trait, into = c("type", "trait"), sep = "_", extra = "merge")

trait_table <- leaf_traits_long %>%
  summarise(mean = round(mean(value, na.rm = TRUE), 2),
            min  = round(min(value, na.rm = TRUE), 2),
            max  = round(max(value, na.rm = TRUE), 2),
            .by = c(trait, type))

# Make histograms of traits
ggplot(leaf_traits_long, aes(x = value, fill = type)) +
  geom_histogram(bins = 20, color = "black", alpha = 0.6, position = "identity") +
  facet_wrap(~ trait, scales = "free_x") +
  scale_fill_manual(values = c("fresh" = "steelblue", "dry" = "tomato")) +
  theme_bw() +
  labs(x = "Value", y = "Count", fill = "Type")

# Combine trait data with TRY data for supplementary table
TRY_traits <- read.csv("TRY_trait_range.csv")
trait_summary <- rbind(trait_table, TRY_traits)
trait_summary <- trait_summary %>%
  mutate(summary = paste0(round(mean, 2), " (", round(min, 2), " - ", round(max, 2), ")")) %>%
  select(trait, type, summary) %>%
  pivot_wider(names_from = type, values_from = summary)
#write.csv(trait_summary, file = "TableS2.csv")

# Read in tree data
tree_data <- read.csv("trees_collected.csv")
tree_data$voucher_no <- as.factor(tree_data$voucher_no)

# Combine leaf with tree data
leaf_traits <- left_join(leaf_traits, tree_data, by = "voucher_no")
leaf_traits$Species <- as.factor(leaf_traits$Species)

# Set thickness of Melastoma beccarianum to NA. Thickness measurements for this species are not reliable.
leaf_traits$fresh_thickness <- ifelse(leaf_traits$Species == "Melastoma beccarianum", NA, leaf_traits$fresh_thickness)
leaf_traits$dry_thickness <- ifelse(leaf_traits$Species == "Melastoma beccarianum", NA, leaf_traits$dry_thickness)
leaf_traits$prop_thickness <- ifelse(leaf_traits$Species == "Melastoma beccarianum", NA, leaf_traits$prop_thickness)

# Make species table
tables1 <- leaf_traits %>%
  group_by(Family, Species) %>%
  summarize(n_ind = length(unique(voucher_no)),
            thickness = paste0(round(mean(fresh_thickness),2)," (", round(min(fresh_thickness),2), "-", round(max(fresh_thickness),2),")"),
            LMA = paste0(round(mean(fresh_LMA),2)," (", round(min(fresh_LMA),2), "-", round(max(fresh_LMA),2),")"),
            ELW = paste0(round(mean(fresh_ELW),2)," (", round(min(fresh_ELW),2), "-", round(max(fresh_ELW),2),")"),
            area = paste0(round(mean(fresh_area),2)," (", round(min(fresh_area),2), "-", round(max(fresh_area),2),")"),
            length = paste0(round(mean(fresh_length),2)," (", round(min(fresh_length),2), "-", round(max(fresh_length),2),")")) 
#write.csv(tables1, "TableS1.csv", row.names = FALSE)


### ANOVA and Tukey HSD to test for treatment effects on trait proportional change
### List of traits:
traits <- c("thickness", "LMA", "ELW", "area", "length")

anova_models <- list()
anova_plots <- list()

for (trait in traits) {

  cat("Trait:", trait, "\n")
  
  prop_trait <- paste0("prop_", trait)
  
  anova_formula <- as.formula(
    paste(prop_trait, "~ dry_treatment * alcohol_treatment")
  )
  
  anova_model <- aov(anova_formula, data = leaf_traits)
  anova_models[[trait]] <- anova_model
  
  print(summary(anova_model))
  print(shapiro.test(residuals(anova_model)))
  
  tukey_result <- TukeyHSD(anova_model)
  tukey_cld <- multcompLetters4(anova_model, tukey_result)
  
  cld_df <- data.frame(
    Letters = tukey_cld$`dry_treatment:alcohol_treatment`$Letters)
  cld_df$dry_treatment <- sapply(strsplit(rownames(cld_df), ":"), `[`, 1)
  cld_df$alcohol_treatment <- sapply(strsplit(rownames(cld_df), ":"), `[`, 2)
  
  # Find the global max for this trait to place all labels at the same height
  global_max <- max(leaf_traits[[prop_trait]], na.rm = TRUE)
  
  label_positions <- leaf_traits %>%
    group_by(dry_treatment, alcohol_treatment) %>%
    summarise(ymax = global_max, .groups = "drop") %>%
    left_join(cld_df, by = c("dry_treatment", "alcohol_treatment"))
  
  label_positions$dry_treatment <- factor(label_positions$dry_treatment,
                                          levels = c("control", "extra dry"))
  label_positions$alcohol_treatment <- factor(label_positions$alcohol_treatment,
                                              levels = c("control", "alcohol", "extra alcohol"))
  
  # Boxplot with Tukey letters
  p <- ggplot(leaf_traits, aes(x = dry_treatment, y = .data[[prop_trait]], fill = alcohol_treatment)) +
    geom_boxplot() +
    geom_text(data = label_positions, aes(label = Letters, y = ymax, group = alcohol_treatment),
              position = position_dodge(width = 0.75),
              vjust = -0.5,      # small gap above the global max line
              size = 5,
              fontface = "bold") +
    labs(title = trait) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text()) +
    scale_fill_paletteer_d("ggthemes::Tableau_10") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    scale_x_discrete(labels = c("DC", "DE")) 
  
  anova_plots[[trait]] <- p
}

anova_plot <- ggarrange(plotlist = anova_plots, ncol = 3, nrow = 2, common.legend = TRUE, legend = "right") 
annotate_figure(p = anova_plot, left = text_grob("Ratio (dry/fresh)", rot = 90, vjust = 1, size = 14),
                  bottom = text_grob("Drying treatment",vjust = -0.5, size = 14))

# Summarize ANOVA model results
anova_summary <- map_dfr(traits, function(trait) {
  tidy(anova_models[[trait]]) %>%
    mutate(trait = trait)
}) %>%
  select(trait, term, df, sumsq, meansq, statistic, p.value) %>%
  mutate(
    across(c(sumsq, meansq), ~ round(.x, 3)),
    across(c(statistic), ~ round(.x, 2)),
    p.value = case_when(
      p.value < 0.001 ~ "<0.001",
      TRUE ~ as.character(round(p.value, 3))),
    term = recode(term,
                  "dry_treatment"                    = "Dry",
                  "alcohol_treatment"                = "Alcohol",
                  "dry_treatment:alcohol_treatment"  = "Dry x Alcohol"
    ),
    trait = str_to_title(trait)) %>%
  rename(
    Trait     = trait,
    Term      = term,
    df        = df,
    SS        = sumsq,
    MS        = meansq,
    F         = statistic,
    `p-value` = p.value)

#write.csv(anova_summary, file = "TableS3.csv", row.names = FALSE)

### GLMMs ###
# First a global model
global_plots  <- list()
global_models <- list()

for (trait in traits) {
  
  dry_trait   <- paste0("dry_", trait)
  fresh_trait <- paste0("fresh_", trait)
  
  ### Correlation between dry and fresh measurements with treatments as fixed effects
  m <- lmer(as.formula(paste(dry_trait, "~", fresh_trait, " + dry_treatment + alcohol_treatment ",
                             "+", fresh_trait, " : dry_treatment", "+", fresh_trait,  " : alcohol_treatment", 
                             "+ (1 | Species/voucher_no)")),
            data = leaf_traits)
  summary <- summary(m)
  p_val <- coef(summary)[fresh_trait, "Pr(>|t|)"]
  r2 <- r.squaredGLMM(m)
  
  # Scatter plot: fresh vs dry
  p <- ggplot(leaf_traits, aes(x = .data[[dry_trait]], y = .data[[fresh_trait]], color = alcohol_treatment, linetype = dry_treatment, shape = dry_treatment)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1, size = 3.5,
             label = paste0("p-value = ", format.pval(p_val, digits = 2), "\n",
                            "R² = ", round(r2[2], 3))) +
    labs(title = trait) +
    xlab("") + 
    ylab("") +
    scale_color_paletteer_d("ggthemes::Tableau_10") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")
  
  global_plots[[trait]] <- p
  global_models[[trait]] <- m
}

traits_plot <- ggarrange(plotlist = global_plots, ncol = 2, nrow = 3, common.legend = TRUE, legend = "right") 
annotate_figure(p = traits_plot, left = text_grob("Dry measurement", rot = 90, vjust = 1, size = 14),
                bottom = text_grob("Fresh measurement", vjust = -0.5, size = 14))


# Extract summary tables for each trait
extract_lmer_table <- function(model_list) {
  
  results <- data.frame()
  
  for (trait in names(model_list)) {
    coefs <- as.data.frame(coef(summary(model_list[[trait]])))
    coefs$term  <- rownames(coefs)
    coefs$trait <- trait
    rownames(coefs) <- NULL
    results <- rbind(results, coefs)
  }
  
  results <- results[, c("trait", "term", "Estimate", "Std. Error", "df", "t value", "Pr(>|t|)")]
  names(results) <- c("trait", "term", "estimate", "std_error", "df", "t_value", "p_value")
  results[, c("estimate", "std_error", "df", "t_value", "p_value")] <-
    round(results[, c("estimate", "std_error", "df", "t_value", "p_value")], 3)
  
  return(results)
}

results_table <- extract_lmer_table(global_models)
#write.csv(results_table, file = "TableS4.csv", row.names = F)

### Calculate how much variance is introduced by alcohol and drying vs intraspecific variance
### Variance Component Decomposition ###

variance_models <- list()
vc_results    <- list()

for (trait in traits) {
  
  dry_trait   <- paste0("dry_", trait)
  fresh_trait <- paste0("fresh_", trait)
  
  # Random intercepts for treatments and species to allow variance decomposition
  m <- lmer(as.formula(paste(dry_trait, "~", fresh_trait,
                             "+ (1 | dry_treatment)",
                             "+ (1 | alcohol_treatment)",
                             "+ (1 | Species/voucher_no)")),
            data = leaf_traits)
  
  summary <- summary(m)
  p_val   <- coef(summary)[fresh_trait, "Pr(>|t|)"]
  r2      <- r.squaredGLMM(m)
  
  # Variance component decomposition
  vc <- as.data.frame(VarCorr(m))
  vc <- vc[is.na(vc$var2), ]          # keep variance rows only (drop covariances)
  
  total_var  <- sum(vc$vcov)
  vc$trait   <- trait
  vc$pct_var <- round((vc$vcov / total_var) * 100, 2)
  
  vc_results[[trait]] <- vc[, c("trait", "grp", "vcov", "pct_var")]
  variance_models[[trait]] <- m
}

# Compile variance components across all traits
vc_table <- do.call(rbind, vc_results)
names(vc_table) <- c("trait", "grouping_factor", "variance", "pct_total_var")
rownames(vc_table) <- NULL

# Reorder grouping factors 
vc_table$grouping_factor <- factor(vc_table$grouping_factor, levels = c("alcohol_treatment", "dry_treatment", "Species", "voucher_no:Species", "Residual"), ordered = T)
vc_table$trait <- factor(vc_table$trait, levels = c("thickness", "LMA", "ELW", "area", "length"), ordered = T)
# Multiply variance by 1000 for thickness to put it on a larger scale to be able to actually interpret the values
vc_table$variance <- ifelse(vc_table$trait == "thickness", vc_table$variance * 1000, vc_table$variance)
vc_table$variance <- round(vc_table$variance, 3)
#write.csv(vc_table, "TableS5.csv", row.names = FALSE)

#Plot variance components
ggplot(vc_table, aes(x = trait, y = pct_total_var, fill = grouping_factor)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_paletteer_d("LaCroixColoR::KeyLime", name = "Random effect") +
  labs(y = "% of Total Variance",
       x = "Leaf Trait") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))


#########################
### Spectral analyses ###
#########################

# Read in spectra
dry_spectra <- read.csv("dry_hyperspectral_data.csv")
dry_spectra <- dry_spectra %>%
  mutate(treatment_code = gsub("-", "_", treatment_code))
dry_spectra$voucher_no <- as.factor(dry_spectra$voucher_no)

# Trim spectra
dry_spectra <- dry_spectra[,-c(6:105, 2057:2156)]

# Add hyperspectral data to leaf traits
dry_leaf_all <- left_join(leaf_traits, dry_spectra, by = join_by(voucher_no, treatment_code, leaf))

# Calculate mean reflectance between AB and AD for each leaf
wavelength_cols_trimmed <- c(paste0("Wavelength_", 450:2400))

# Calculate mean and SD spectra for each treatment combination
spectral_data_long <- dry_leaf_all %>%
  pivot_longer(cols = starts_with("Wavelength_"),
               names_to = "wavelength_col",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(str_remove(wavelength_col, "Wavelength_")))

spectral_data_long <- spectral_data_long %>%
  mutate(alcohol_treatment = factor(alcohol_treatment, 
                               levels = c("control", "alcohol", "extra alcohol"),
                               ordered = FALSE))

spectral_summary <- spectral_data_long %>%
  group_by(dry_treatment, alcohol_treatment, wavelength) %>%
  summarise(
    mean_reflectance = mean(reflectance, na.rm = TRUE),
    sd_reflectance = sd(reflectance, na.rm = TRUE),
    n = n(),
    se_reflectance = sd_reflectance / sqrt(n),
    .groups = "drop"
  )

ggplot(spectral_summary, aes(x = wavelength, y = mean_reflectance, color = alcohol_treatment)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_reflectance - sd_reflectance,
                  ymax = mean_reflectance + sd_reflectance,
                  fill = alcohol_treatment), alpha = 0.15, colour = NA) +
  scale_x_continuous(breaks = seq(450, 2400, by = 250)) +
  facet_wrap(vars(dry_treatment), nrow = 2) +
  labs(
       x = "Wavelength (nm)",
       y = "Mean reflectance",
       colour = "Alcohol treatment",
       fill   = "Alcohol treatment") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top")

### Test for differences in mean reflectance at each wavelength.
# Run an anova at each wavelength to test for a treatment effect on mean reflectance
wavelength_mean_tests <- spectral_data_long %>%
  group_by(wavelength) %>%
  do({
    model <- aov(reflectance ~ dry_treatment * alcohol_treatment, data = .)
    aov_result <- anova(model)
    tibble(
      drying_F = aov_result["dry_treatment", "F value"],
      drying_p = aov_result["dry_treatment", "Pr(>F)"],
      alcohol_F = aov_result["alcohol_treatment", "F value"],
      alcohol_p = aov_result["alcohol_treatment", "Pr(>F)"],
      interaction_F = aov_result["dry_treatment:alcohol_treatment", "F value"],
      interaction_p = aov_result["dry_treatment:alcohol_treatment", "Pr(>F)"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    drying_p_adj = p.adjust(drying_p, method = "fdr"),
    alcohol_p_adj = p.adjust(alcohol_p, method = "fdr"),
    interaction_p_adj = p.adjust(interaction_p, method = "fdr")
  )

anova_plot_data <- wavelength_mean_tests %>%
  select(wavelength,
         drying_p_adj,
         alcohol_p_adj,
         interaction_p_adj) %>%
  pivot_longer(cols = -wavelength,
               names_to = "term",
               values_to = "p_adj") %>%
  mutate(term = recode(term,
                       drying_p_adj = "Drying main effect",
                       alcohol_p_adj = "Alcohol main effect",
                       interaction_p_adj = "Drying x Alcohol")) 

anova_plot_data$term <- factor(anova_plot_data$term, levels = c("Drying main effect", "Alcohol main effect", "Drying x Alcohol"))

# Create wavelength gradient
wavelengths <- seq(400, 700, length.out = 300)

spec_colors <- alpha(
  colorRampPalette(c("purple", "blue", "cyan", "green", "yellow", "orange", "red"))(length(wavelengths)), 0.7)

gradient_raster <- matrix(spec_colors, nrow = 1)

ggplot(anova_plot_data, aes(wavelength, -log10(p_adj))) +
  annotation_raster(raster = gradient_raster, xmin = 400, xmax = 700, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray95", alpha = 0.7, xmin = 700, xmax = 1100, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray75", alpha = 0.7, xmin = 1100, xmax = 2000, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray55", alpha = 0.7, xmin = 2000, xmax = 2500, ymin = -Inf, ymax = Inf) +
  annotate("text", x = 550, y = 40, label = "VIS", size = 4) +
  annotate("text", x = 900, y = 40, label = "NIR", size = 4) +
  annotate("text", x = 1550, y = 40, label = "SWIR 1", size = 4) +
  annotate("text", x = 2250, y = 40, label = "SWIR 2", size = 4) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ term, ncol = 1) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "red") +
  labs(x = "Wavelength (nm)",
       y = expression(-log[10]("adjusted p-value"))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))

sig_regions <- wavelength_mean_tests %>%
  filter(drying_p_adj < 0.05 | alcohol_p_adj < 0.05 | interaction_p_adj < 0.05) %>%
  mutate(
    effect_type = case_when(
      drying_p_adj < 0.05 & alcohol_p_adj < 0.05 & interaction_p_adj < 0.05 ~ "All three",
      drying_p_adj < 0.05 & alcohol_p_adj < 0.05 ~ "Drying + Alcohol",
      drying_p_adj < 0.05 & interaction_p_adj < 0.05 ~ "Drying + Interaction",
      alcohol_p_adj < 0.05 & interaction_p_adj < 0.05 ~ "Alcohol + Interaction",
      drying_p_adj < 0.05 ~ "Drying only",
      alcohol_p_adj < 0.05 ~ "Alcohol only",
      interaction_p_adj < 0.05 ~ "Interaction only"
    )
  )

sig_ribbons <- wavelength_mean_tests %>%
  mutate(
    sig_drying = ifelse(drying_p_adj < 0.05, 1, 0),
    sig_alcohol = ifelse(alcohol_p_adj < 0.05, 1, 0),
    sig_interaction = ifelse(interaction_p_adj < 0.05, 1, 0)
  )

ggplot(wavelength_mean_tests, aes(x = wavelength)) +
  # Add shaded regions for significance
  geom_ribbon(data = sig_ribbons, 
              aes(ymin = 0, ymax = ifelse(sig_drying == 1, 40, 0)),
              fill = "blue", alpha = 0.2) +
  geom_ribbon(data = sig_ribbons, 
              aes(ymin = 0, ymax = ifelse(sig_alcohol == 1, 40, 0)),
              fill = "red", alpha = 0.2) +
  geom_ribbon(data = sig_ribbons, 
              aes(ymin = 0, ymax = ifelse(sig_interaction == 1, 40, 0)),
              fill = "purple", alpha = 0.4) +
  # Plot -log10(p) values
  geom_line(aes(y = -log10(drying_p_adj), color = "Drying"), linewidth = 1) +
  geom_line(aes(y = -log10(alcohol_p_adj), color = "Alcohol"), linewidth = 1) +
  geom_line(aes(y = -log10(interaction_p_adj), color = "Interaction"), linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(x = "Wavelength (nm)", 
       y = "-log10(adjusted p-value)",
       color = "Treatment") +
  scale_color_manual(values = c("Drying" = "blue", 
                                "Alcohol" = "red", 
                                "Interaction" = "purple")) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

# Now pairwise contrasts to see if different concentrations of alcohol affect spectra
contrast_results <- spectral_data_long %>%
  group_by(wavelength) %>%
  do({
    fit <- aov(reflectance ~ dry_treatment + alcohol_treatment, data = .)
    
    ## Estimated marginal means
    em_alc <- emmeans(fit, ~ alcohol_treatment)
    
    ## Contrasts
    ct_alc <- contrast(
      em_alc,
      method = list(
        "alc50 - control" = c(-1,  1,  0),
        "alc80 - control" = c(-1,  0,  1),
        "alc80 - alc50" = c(0, -1, 1)))
    
    bind_rows(tidy(ct_alc))
    
  }) %>%
  ungroup()

# p-value adjustment for multiple-testing correction
contrast_results <- contrast_results %>%
  group_by(contrast) %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr")) %>%
  mutate(contrast = recode(contrast,
                           `alc50 - control` = "A50 - Control",
                           `alc80 - control` = "A80 - Control",
                           `alc80 - alc50` = "A80 - A50"))

contrast_results$contrast <- factor(contrast_results$contrast, levels = c("A50 - Control", "A80 - Control", "A80 - A50"))

# Now plot the contrast results
ggplot(contrast_results, aes(wavelength, -log10(p_adj))) +
  annotation_raster(raster = gradient_raster, xmin = 400, xmax = 700, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray95", alpha = 0.6, xmin = 700, xmax = 1100, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray75", alpha = 0.6, xmin = 1100, xmax = 2000, ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray55", alpha = 0.6, xmin = 2000, xmax = 2500, ymin = -Inf, ymax = Inf) +
  annotate("text", x = 550, y = 40, label = "VIS", size = 4) +
  annotate("text", x = 900, y = 40, label = "NIR", size = 4) +
  annotate("text", x = 1550, y = 40, label = "SWIR 1", size = 4) +
  annotate("text", x = 2250, y = 40, label = "SWIR 2", size = 4) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  facet_wrap(~ contrast, ncol = 1) +
  labs(x = "Wavelength (nm)",
       y = expression(-log[10]("adjusted p-value"))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))


#############################
########## PermANOVA ########

# Make spectral matrix
spec_matrix <- dry_leaf_all %>% select(all_of(wavelength_cols_trimmed))

# Compute bray-curtis distance on the matrix
dist_bc <- vegdist(spec_matrix, method = "bray")

perm_dry <- adonis2(
  dist_bc ~ dry_treatment * alcohol_treatment,
  data         = dry_leaf_all,
  permutations = 999,
  by           = "terms",
  strata       = dry_leaf_all$voucher_no)
summary(perm_dry)
#write.csv(perm_dry, file = "TableS6.csv")

nmds_dry <- metaMDS(dist_bc, k = 2, trymax = 100, trace = FALSE)
nmds_dry$stress

nmds_scores <- as.data.frame(scores(nmds_dry, display = "sites")) %>%
  bind_cols(dry_leaf_all %>% select(dry_treatment, alcohol_treatment, Species))

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2,
                        color = alcohol_treatment,
                        shape = dry_treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = interaction(dry_treatment, alcohol_treatment), linetype = dry_treatment),
              linewidth = 1, level = 0.95) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  labs(color = "Alcohol treatment",
       shape = "Drying treatment",
       linetype = "Drying treatment") +
  theme_classic(base_size = 14) +
  theme(legend.background = element_blank(),
    legend.box.background = element_rect(color = "black")) 

############################################
# Now calculate first and second derivatives
dry_deriv1 <- as.data.frame(t(apply(dry_leaf_all[, grep("Wavelength_", colnames(dry_leaf_all))], 1, diff)))
colnames(dry_deriv1) <- paste0("Deriv_Wavelength_", 451:2400)

dry_deriv2 <- as.data.frame(t(apply(dry_deriv1, 1, diff)))
colnames(dry_deriv2) <- paste0("Deriv2_Wavelength_", 452:2400)

# Combine with original data
dry_leaf_deriv <- cbind(dry_leaf_all, dry_deriv1, dry_deriv2)

# Create unique specimen ID (voucher + treatment)
dry_leaf_deriv <- dry_leaf_deriv %>%
  mutate(specimen_id = paste(voucher_no, treatment_code, sep = "_"))

# Create vector of derivative wavelength column names
deriv_wavelength_cols <- c(
  paste0("Deriv_Wavelength_", 451:2400),
  paste0("Deriv2_Wavelength_", 452:2400)
)


























