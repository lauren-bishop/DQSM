# Clinical ~ Age + Sex + QSM * ROI + (1 | Subject)

library(readxl)
library(dplyr)
library(lme4)
library(lmerTest)   # gives Satterthwaite p-values for lmer via emmeans
library(emmeans)
library(ggplot2)

# Paths

DATA_DIR <- "/Users/laurenbishop/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Research Assistant - AZ.RW/QSM Dipole Seperation/results"
ROI_DIR  <- file.path(DATA_DIR, "3Altas")
OUTDIR   <- file.path(DATA_DIR, "patient_QSM_clinical_associations_bilateral")

DIAG_FILE <- file.path(DATA_DIR, "subject_diagnosis_cognitive_merged_with_LBD.xlsx")

ROI_FILES <- c(
  "compS_from_3AtlasCompositeBilat.csv",
  "DCS_from_3AtlasCompositeBilat.csv",
  "PCS_from_3AtlasCompositeBilat.csv"
)

ALPHA <- 0.05

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)


# Variables

PATIENT_DIAGNOSES <- c("PD", "PDD", "PD-MCI", "LB-MCI", "DLB")

CLINICAL_COLUMNS <- c(
  UPDRS_Total        = "UPDRS.Total",
  UPDRS_III          = "UPDRS_III",
  MoCA               = "MoCA.x",
  MMSE               = "MMSE.x",
  MiniACE            = "Mini.ACE",
  CognitiveComposite = "Cognitive.Composite"
)

EXCLUDE_IDS <- c(
  842, 926, 1049, 1093, 1099, 1101, 1103, 1108,
  1120, 1122, 1127, 792, 793, 1109, 832, 836,
  1027, 2013
)


# Helper

extract_id <- function(x) {
  x <- as.character(x)
  as.integer(sub(".*?(\\d+).*", "\\1", x))
}


# Clinical data

diag <- read_excel(DIAG_FILE) |>
  as.data.frame()

required_cols <- c(
  "Subject", "Sex", "Age", "Diagnosis",
  unname(CLINICAL_COLUMNS)
)

missing_cols <- setdiff(required_cols, names(diag))

if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse = ", "))
}

names(diag)[match(
  unname(CLINICAL_COLUMNS),
  names(diag)
)] <- names(CLINICAL_COLUMNS)

# Convert numeric variables
for (v in c("Age", names(CLINICAL_COLUMNS))) {
  diag[[v]] <- suppressWarnings(as.numeric(diag[[v]]))
}

# Clean sex variable
diag$Sex <- toupper(trimws(as.character(diag$Sex)))
diag$Sex[diag$Sex %in% c("FEMALE", "2")] <- "F"
diag$Sex[diag$Sex %in% c("MALE", "1")] <- "M"
diag$Sex[!diag$Sex %in% c("F", "M")] <- NA
diag$Sex <- factor(diag$Sex, levels = c("F", "M"))

# Keep patients and exclude specified subjects

diag <- diag |>
  filter(Diagnosis %in% PATIENT_DIAGNOSES) |>
  mutate(SubjectID = extract_id(Subject)) |>
  filter(!is.na(SubjectID), !SubjectID %in% EXCLUDE_IDS) |>
  select(-SubjectID)


# QSM data

master_long_list <- list()

for (f in ROI_FILES) {
  
  message("Reading: ", f)
  
  path <- file.path(ROI_DIR, f)
  
  if (!file.exists(path)) {
    warning("File does not exist: ", path)
    next
  }
  
  df <- read.csv(path, stringsAsFactors = FALSE)
  
  if (!"Subject" %in% names(df)) {
    warning("No Subject column in ", f)
    next
  }
  
  ids <- extract_id(df$Subject)
  
  df <- df[
    !is.na(ids) & !ids %in% EXCLUDE_IDS,
    ,
    drop = FALSE
  ]
  
  mean_cols <- grep("_mean$", names(df), value = TRUE)
  
  if (!length(mean_cols)) {
    warning("No ROI mean columns found in ", f)
    next
  }
  
  message("  Bilateral ROI columns: ", length(mean_cols))
  
  master_long_list[[f]] <- data.frame(
    Subject = rep(df$Subject, times = length(mean_cols)),
    ROI = rep(sub("_mean$", "", mean_cols), each = nrow(df)),
    QSM = as.numeric(unlist(
      df[, mean_cols, drop = FALSE],
      use.names = FALSE
    )),
    Measure = sub("_from_.*$", "", f),
    stringsAsFactors = FALSE
  )
}

master_long <- bind_rows(master_long_list)

if (!nrow(master_long)) {
  stop("No bilateral ROI data were successfully loaded.")
}


# Merge clinical and QSM data

master_long <- master_long |>
  inner_join(diag, by = "Subject") |>
  mutate(
    ROI = factor(ROI),
    Subject = factor(Subject)   # needed as a grouping factor for (1 | Subject)
  )


# Analysis
#
# Model direction swapped per review comment:
#   clinical_var ~ Age + Sex + QSM * ROI + (1 | Subject)
# Each subject contributes one row per bilateral ROI within a given Measure,
# so those rows are NOT independent -> a per-subject random intercept is
# required instead of a plain lm(). The quantity of interest is now the
# slope of QSM on the clinical score, estimated separately per ROI
# (emtrends on QSM, grouped by ROI), adjusted for Age and Sex.

run_clinical_analysis <- function(data, clinical_var) {
  
  d <- data |>
    filter(
      !is.na(.data[[clinical_var]]),
      !is.na(QSM),
      !is.na(Age),
      !is.na(Sex)
    )
  
  message(
    clinical_var, ": ",
    length(unique(d$Subject)),
    " complete-case subjects"
  )
  
  if (length(unique(d$Subject)) < 10) {
    message("Too few subjects - skipping")
    return(NULL)
  }
  
  results <- list()
  
  for (m in unique(d$Measure)) {
    
    dm <- d |>
      filter(Measure == m) |>
      droplevels()
    
    if (length(unique(dm[[clinical_var]])) < 2) {
      message("  ", m, ": insufficient clinical score variation")
      next
    }
    
    # Need multiple ROIs per subject for the random intercept to make sense;
    # if a subject only ever has one row here, (1 | Subject) is not doing
    # anything useful but it's harmless.
    formula <- as.formula(
      paste0(clinical_var, " ~ Age + Sex + QSM * ROI + (1 | Subject)")
    )
    
    mdl <- tryCatch(
      lmerTest::lmer(
        formula,
        data = dm,
        REML = TRUE,
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 500000)
        )
      ),
      error = function(e) {
        message("  Model failed: ", e$message)
        NULL
      },
      warning = function(w) {
        message("  Model warning: ", w$message)
        # still try to fit and keep the model rather than dropping it
        suppressWarnings(
          lmerTest::lmer(
            formula,
            data = dm,
            REML = TRUE,
            control = lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 500000)
            )
          )
        )
      }
    )
    
    if (is.null(mdl)) next
    
    if (isSingular(mdl, tol = 1e-4)) {
      message("  ", m, ": singular fit (random-effect variance ~0) - interpret with caution")
    }
    
    trends <- tryCatch(
      emtrends(mdl, ~ ROI, var = "QSM"),
      error = function(e) {
        message("  emtrends failed: ", e$message)
        NULL
      }
    )
    
    if (is.null(trends)) next
    
    trend_col <- "QSM.trend"
    
    res <- as.data.frame(summary(trends, infer = TRUE))
    
    if (!trend_col %in% names(res)) {
      message("  Trend column not found: ", trend_col)
      next
    }
    
    res <- res |>
      select(ROI, all_of(trend_col), SE, df, t.ratio, p.value) |>
      rename(
        slope = all_of(trend_col),
        p_value = p.value
      ) |>
      mutate(
        Measure = m,
        ClinicalVariable = clinical_var,
        ROI_global = paste0(m, "_", ROI)
      )
    
    results[[m]] <- res
  }
  
  if (!length(results)) {
    return(NULL)
  }
  
  results <- bind_rows(results)
  
  out_file <- file.path(
    OUTDIR,
    paste0(clinical_var, "_BILATERAL_QSM_associations.csv")
  )
  
  write.csv(results, out_file, row.names = FALSE)
  
  results
}


# Run analyses

all_results <- lapply(
  names(CLINICAL_COLUMNS),
  function(x) run_clinical_analysis(master_long, x)
)

names(all_results) <- names(CLINICAL_COLUMNS)

all_results_combined <- bind_rows(all_results)


# Global FDR

all_results_combined <- all_results_combined |>
  mutate(
    p_fdr_GLOBAL_ALL = p.adjust(p_value, method = "BH"),
    significant_GLOBAL_ALL = p_fdr_GLOBAL_ALL < ALPHA
  )

write.csv(
  all_results_combined,
  file.path(
    OUTDIR,
    "ALL_BILATERAL_CLINICAL_QSM_ASSOCIATIONS_GLOBAL_FDR.csv"
  ),
  row.names = FALSE
)

sig_global <- all_results_combined |>
  filter(significant_GLOBAL_ALL) |>
  arrange(p_fdr_GLOBAL_ALL)

print(
  sig_global |>
    select(
      Measure,
      ROI,
      ClinicalVariable,
      slope,
      SE,
      df,
      t.ratio,
      p_value,
      p_fdr_GLOBAL_ALL
    )
)


# Helper: fixed-effects-only predictions with 95% CI for a merMod, ignoring
# random-effect variance (i.e. a population-average "confidence" band,
# analogous to predict.lm(..., interval = "confidence")).

predict_fixed_ci <- function(mdl, newdata, df_for_ci = NULL) {
  
  fe_formula <- lme4::nobars(formula(mdl))[-2]   # drop random terms & response
  X <- model.matrix(fe_formula, data = newdata)
  
  beta <- lme4::fixef(mdl)
  X <- X[, names(beta), drop = FALSE]
  
  fit <- as.numeric(X %*% beta)
  
  V <- as.matrix(vcov(mdl))
  se_fit <- sqrt(rowSums((X %*% V) * X))
  
  # crude but standard: use large-sample Normal quantile unless a df is given
  crit <- if (!is.null(df_for_ci)) qt(0.975, df_for_ci) else qnorm(0.975)
  
  data.frame(
    fit = fit,
    lwr = fit - crit * se_fit,
    upr = fit + crit * se_fit
  )
}


# Significant association plots
# Axes are flipped relative to the old QSM ~ clinical model: QSM is now the
# x-axis (predictor) and the clinical score is the y-axis (outcome).

PLOT_DIR <- file.path(OUTDIR, "significant_association_plots")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!nrow(sig_global)) {
  
  message("No associations survived global FDR.")
  
} else {
  
  for (i in seq_len(nrow(sig_global))) {
    
    row_i <- sig_global[i, ]
    
    clinical_var <- as.character(row_i$ClinicalVariable)
    measure <- as.character(row_i$Measure)
    roi <- as.character(row_i$ROI)
    
    plot_data <- master_long |>
      filter(
        Measure == measure,
        ROI == roi,
        !is.na(QSM),
        !is.na(.data[[clinical_var]]),
        !is.na(Age),
        !is.na(Sex)
      )
    
    if (nrow(plot_data) < 3) next
    
    model_data <- master_long |>
      filter(
        Measure == measure,
        !is.na(QSM),
        !is.na(.data[[clinical_var]]),
        !is.na(Age),
        !is.na(Sex)
      ) |>
      mutate(
        ROI = droplevels(factor(ROI)),
        Sex = droplevels(factor(Sex)),
        Subject = droplevels(factor(Subject))
      )
    
    formula <- as.formula(
      paste0(clinical_var, " ~ Age + Sex + QSM * ROI + (1 | Subject)")
    )
    
    mdl <- tryCatch(
      suppressWarnings(
        lmerTest::lmer(
          formula,
          data = model_data,
          REML = TRUE,
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 500000)
          )
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(mdl)) next
    
    x_values <- seq(
      min(plot_data$QSM, na.rm = TRUE),
      max(plot_data$QSM, na.rm = TRUE),
      length.out = 100
    )
    
    prediction_grid <- data.frame(
      Age = mean(model_data$Age, na.rm = TRUE),
      Sex = factor(
        levels(model_data$Sex)[1],
        levels = levels(model_data$Sex)
      ),
      ROI = factor(
        roi,
        levels = levels(model_data$ROI)
      )
    )
    
    prediction_grid <- prediction_grid[
      rep(1, 100),
      ,
      drop = FALSE
    ]
    
    prediction_grid$QSM <- x_values
    
    ci <- predict_fixed_ci(
      mdl,
      prediction_grid,
      df_for_ci = as.numeric(row_i$df)
    )
    
    prediction_grid$fit <- ci$fit
    prediction_grid$lwr <- ci$lwr
    prediction_grid$upr <- ci$upr
    
    p_global <- as.numeric(row_i$p_fdr_GLOBAL_ALL)
    
    p_label <- if (p_global < 0.001) {
      "Global FDR < 0.001"
    } else {
      paste0("Global FDR = ", signif(p_global, 3))
    }
    
    slope_label <- paste0(
      "Adjusted slope (", clinical_var, " per unit QSM) = ",
      signif(as.numeric(row_i$slope), 3)
    )
    
    gp <- ggplot(
      plot_data,
      aes(
        x = QSM,
        y = .data[[clinical_var]]
      )
    ) +
      geom_point(alpha = 0.65, size = 2) +
      geom_ribbon(
        data = prediction_grid,
        aes(
          x = QSM,
          ymin = lwr,
          ymax = upr
        ),
        inherit.aes = FALSE,
        alpha = 0.2
      ) +
      geom_line(
        data = prediction_grid,
        aes(
          x = QSM,
          y = fit
        ),
        inherit.aes = FALSE,
        linewidth = 1
      ) +
      labs(
        title = paste0(measure, " - ", roi, " (bilateral)"),
        subtitle = paste(
          clinical_var,
          slope_label,
          p_label,
          sep = " | "
        ),
        x = "QSM (bilateral mean ROI value)",
        y = clinical_var
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11)
      )
    
    filename <- paste0(
      gsub("[^A-Za-z0-9]+", "-", measure), "_",
      gsub("[^A-Za-z0-9]+", "-", roi), "_",
      gsub("[^A-Za-z0-9]+", "-", clinical_var),
      "_bilateral_significant.png"
    )
    
    ggsave(
      file.path(PLOT_DIR, filename),
      gp,
      width = 6,
      height = 5.5,
      dpi = 300
    )
  }
}


# Significant association summary plot

sig_plot <- all_results_combined |>
  filter(significant_GLOBAL_ALL) |>
  mutate(
    Direction = case_when(
      slope > 0 ~ "Positive",
      slope < 0 ~ "Negative",
      TRUE ~ NA_character_
    ),
    Direction = factor(
      Direction,
      levels = c("Negative", "Positive")
    ),
    EffectStrength = abs(t.ratio),
    ClinicalVariable = factor(
      ClinicalVariable,
      levels = names(CLINICAL_COLUMNS)
    ),
    Measure = factor(
      Measure,
      levels = c("compS", "DCS", "PCS")
    )
  )

if (nrow(sig_plot)) {
  
  clinical_labels <- c(
    UPDRS_Total = "UPDRS Total",
    UPDRS_III = "UPDRS III",
    MoCA = "MoCA",
    MMSE = "MMSE",
    MiniACE = "Mini-ACE",
    CognitiveComposite = "Cognitive Composite"
  )
  
  roi_plot <- ggplot(
    sig_plot,
    aes(
      x = ClinicalVariable,
      y = ROI
    )
  ) +
    geom_point(
      aes(
        size = EffectStrength,
        colour = Direction
      ),
      alpha = 0.85
    ) +
    facet_wrap(
      ~ Measure,
      ncol = 1,
      scales = "free_y"
    ) +
    scale_colour_manual(
      values = c(
        Negative = "red",
        Positive = "blue"
      ),
      name = "Association"
    ) +
    scale_size_continuous(
      name = "|t|",
      range = c(2.5, 10)
    ) +
    scale_x_discrete(labels = clinical_labels) +
    labs(
      title = "Globally Significant Bilateral QSM-Clinical Associations",
      subtitle = paste(
        "BH FDR < 0.05 across clinical variables, MRI measures and bilateral ROIs",
        "Bubble size represents |t|; colour represents direction",
        sep = "\n"
      ),
      x = "Clinical Measure",
      y = "Bilateral ROI"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(
        face = "bold",
        size = 14
      ),
      panel.grid.minor = element_blank()
    )
  
  print(roi_plot)
  
  ggsave(
    file.path(
      OUTDIR,
      "GLOBAL_FDR_SIGNIFICANT_ROI_CLINICAL_BUBBLE_PLOT.png"
    ),
    roi_plot,
    width = 11,
    height = 16,
    dpi = 300
  )
}


# Save significant results

write.csv(
  sig_global,
  file.path(
    OUTDIR,
    "SIGNIFICANT_BILATERAL_QSM_CLINICAL_ASSOCIATIONS_ONLY.csv"
  ),
  row.names = FALSE
)

message(
  "Analysis complete. ",
  nrow(sig_global),
  " associations survived global FDR."
)

