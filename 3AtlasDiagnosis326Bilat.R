DATA_DIR <- "/Users/laurenbishop/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Research Assistant - AZ.RW/QSM Dipole Seperation/results"

ROI_DIR <- file.path(DATA_DIR, "3Altas")

DIAG_FILE <- file.path(
  DATA_DIR,
  "subject_diagnosis_cognitive_merged_with_LBD.xlsx"
)

ROI_FILES <- c(
  "compS_from_3AtlasCompositeBilat.csv",
  "DCS_from_3AtlasCompositeBilat.csv",
  "PCS_from_3AtlasCompositeBilat.csv"
)

DIAG_LEVELS <- c("HC", "PD", "PDD", "LB-MCI", "DLB", "PD-MCI")
COGNITIVE_LEVELS <- c("PDD", "LB-MCI", "DLB", "PD-MCI")

EXCLUDE_IDS <- c(
  842, 926, 1049, 1093, 1099, 1101, 1103, 1108,
  1120, 1122, 1127, 792, 793, 1109, 832, 836,
  1027, 2013
)

OUTDIR <- file.path(DATA_DIR, "simple_stats_results_per_file_Bilat")
ALPHA <- 0.05

if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)
}

# Packages
library(emmeans)
library(lme4)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Check that files exist
if (!file.exists(DIAG_FILE)) {
  stop("Diagnosis file not found: ", DIAG_FILE)
}

for (f in ROI_FILES) {
  path <- file.path(ROI_DIR, f)
  
  if (!file.exists(path)) {
    stop("ROI file not found: ", path)
  }
}

# Extract numeric subject ID
extract_id <- function(x) {
  sapply(as.character(x), function(xx) {
    
    m <- regexpr("\\d+", xx)
    
    if (m == -1) {
      return(NA_integer_)
    }
    
    as.integer(
      substr(
        xx,
        m,
        m + attr(m, "match.length") - 1
      )
    )
  })
}

# Read Excel or CSV file
read_table <- function(path) {
  
  ext <- tolower(tools::file_ext(path))
  
  if (ext %in% c("xlsx", "xls")) {
    readxl::read_excel(path)
  } else if (ext == "csv") {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file type: ", ext)
  }
}

# Read diagnosis file
diag_raw <- read_table(DIAG_FILE)

col_lower <- tolower(names(diag_raw))

get_col <- function(pattern) {
  i <- grep(pattern, col_lower)
  
  if (length(i) == 0) {
    return(NA_character_)
  }
  
  names(diag_raw)[i[1]]
}

col_subj <- get_col("^subject$")
col_diag <- get_col("^diagnos")
col_age <- get_col("^age$")
col_sex <- get_col("^sex$")

if (
  is.na(col_subj) ||
  is.na(col_diag) ||
  is.na(col_age) ||
  is.na(col_sex)
) {
  stop("Diagnosis file must contain Subject, Diagnosis, Age and Sex columns.")
}

diag <- diag_raw[, c(
  col_subj,
  col_diag,
  col_age,
  col_sex
)]

names(diag) <- c(
  "Subject",
  "Diagnosis",
  "Age",
  "Sex"
)

# Clean sex
sex_raw <- trimws(toupper(as.character(diag$Sex)))

diag$Sex <- case_when(
  sex_raw %in% c("F", "FEMALE", "2") ~ "F",
  sex_raw %in% c("M", "MALE", "1") ~ "M",
  TRUE ~ NA_character_
)

diag <- diag[!is.na(diag$Sex), ]

# Clean diagnosis and remove excluded subjects
diag <- diag[diag$Diagnosis %in% DIAG_LEVELS, ]

diag$SubjectID <- extract_id(diag$Subject)

diag <- diag[
  !is.na(diag$SubjectID) &
    !diag$SubjectID %in% EXCLUDE_IDS,
]

diag$SubjectID <- NULL
rownames(diag) <- NULL

diag$Diagnosis <- factor(
  diag$Diagnosis,
  levels = DIAG_LEVELS
)

# Create three-group variable
diag$Group3 <- case_when(
  diag$Diagnosis == "HC" ~ "HC",
  diag$Diagnosis == "PD" ~ "PD",
  diag$Diagnosis %in% COGNITIVE_LEVELS ~ "Cognitive",
  TRUE ~ NA_character_
)

diag$Group3 <- factor(
  diag$Group3,
  levels = c("HC", "PD", "Cognitive")
)

# Create two-group variable
diag$Group2 <- ifelse(
  diag$Diagnosis == "HC",
  "Control",
  "Affected"
)

diag$Group2 <- factor(
  diag$Group2,
  levels = c("Control", "Affected")
)

# Diagnosis summary
message("Diagnosis loaded: ", nrow(diag), " subjects.")

message("\nFull 6-group counts:")
print(table(diag$Diagnosis))

message("\nThree-group counts:")
print(table(diag$Group3))

message("\nTwo-group counts:")
print(table(diag$Group2))

# Build master ROI table
master_long_list <- list()

for (f in ROI_FILES) {
  
  path <- file.path(ROI_DIR, f)
  
  df <- read_table(path)
  
  if ("SubjectID" %in% names(df)) {
    df <- df[!df$SubjectID %in% EXCLUDE_IDS, ]
  }
  
  if (!"Subject" %in% names(df)) {
    next
  }
  
  subject_ids <- extract_id(df$Subject)
  
  df <- df[
    !is.na(subject_ids) &
      !subject_ids %in% EXCLUDE_IDS,
  ]
  
  mean_cols <- grep(
    "_mean$",
    names(df),
    value = TRUE
  )
  
  if (length(mean_cols) == 0) {
    next
  }
  
  nsub <- nrow(df)
  
  long <- data.frame(
    Subject = rep(
      df$Subject,
      times = length(mean_cols)
    ),
    ROI = rep(
      sub("_mean$", "", mean_cols),
      each = nsub
    ),
    QSM = as.numeric(
      unlist(
        df[mean_cols],
        use.names = FALSE
      )
    ),
    Measure = sub(
      "_from_.*$",
      "",
      f
    ),
    stringsAsFactors = FALSE
  ) %>%
    filter(!is.na(QSM))
  
  master_long_list[[f]] <- long
}

master_long <- bind_rows(master_long_list)

master_long <- merge(
  master_long,
  diag,
  by = "Subject"
)

master_long$ROI <- factor(master_long$ROI)
master_long$Age <- as.numeric(master_long$Age)

# Main analysis function
run_group_analysis <- function(
    data,
    group_col,
    group_levels,
    ref_level,
    analysis_name,
    outdir,
    alpha
) {
  
  sub_outdir <- file.path(
    outdir,
    analysis_name
  )
  
  plot_dir <- file.path(
    sub_outdir,
    "plots"
  )
  
  dir.create(
    plot_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  data <- data[
    !is.na(data[[group_col]]),
  ]
  
  data$Group <- factor(
    data[[group_col]],
    levels = group_levels
  )
  
  data$Sex <- factor(data$Sex)
  data$Subject <- factor(data$Subject)
  
  results <- list()
  
  for (m in unique(data$Measure)) {
    
    message("  Measure: ", m)
    
    dsub <- data[data$Measure == m, ]
    
    dsub$Group <- droplevels(dsub$Group)
    dsub$ROI <- droplevels(factor(dsub$ROI))
    dsub$Sex <- factor(dsub$Sex)
    dsub$Subject <- factor(dsub$Subject)
    
    if (nlevels(dsub$Group) < 2) {
      next
    }
    
    # Descriptive statistics
    desc_stats <- dsub %>%
      group_by(ROI, Group) %>%
      summarise(
        mean_QSM = mean(QSM, na.rm = TRUE),
        sd_QSM = sd(QSM, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = Group,
        values_from = c(mean_QSM, sd_QSM, n)
      )
    
    # Mixed model with random intercept by subject
    model <- lmer(
      QSM ~ Age + Sex + Group * ROI + (1 | Subject),
      data = dsub,
      REML = FALSE
    )
    
    # Estimated marginal means
    emm <- emmeans(
      model,
      ~ Group | ROI,
      cov.reduce = mean,
      lmer.df = "asymptotic"
    )
    
    # Compare groups against reference group
    ref_here <- if (ref_level %in% levels(dsub$Group)) {
      ref_level
    } else {
      levels(dsub$Group)[1]
    }
    
    contrasts <- emmeans::contrast(
      emm,
      method = "trt.vs.ctrl",
      ref = ref_here,
      adjust = "none"
    )
    
    res <- as.data.frame(
      summary(
        contrasts,
        infer = TRUE,
        lmer.df = "asymptotic"
      )
    )
    
    if (nrow(res) == 0) {
      next
    }
    
    # emmeans with asymptotic df often uses z.ratio rather than t.ratio
    ratio_col <- intersect(
      c("z.ratio", "t.ratio"),
      names(res)
    )
    
    if (length(ratio_col) == 0) {
      stop(
        "Could not find z.ratio or t.ratio in emmeans output. Available columns: ",
        paste(names(res), collapse = ", ")
      )
    }
    
    ratio_col <- ratio_col[1]
    
    required_cols <- c(
      "ROI",
      "contrast",
      "estimate",
      "SE",
      "p.value"
    )
    
    missing_cols <- setdiff(required_cols, names(res))
    
    if (length(missing_cols) > 0) {
      stop(
        "Missing expected columns from emmeans output: ",
        paste(missing_cols, collapse = ", "),
        ". Available columns are: ",
        paste(names(res), collapse = ", ")
      )
    }
    
    res$test_statistic <- res[[ratio_col]]
    
    if (!"df" %in% names(res)) {
      res$df <- Inf
    }
    
    res <- res[
      ,
      c(
        "ROI",
        "contrast",
        "estimate",
        "SE",
        "df",
        "test_statistic",
        "p.value"
      )
    ]
    
    names(res)[names(res) == "p.value"] <- "p_value"
    
    res <- merge(
      res,
      desc_stats,
      by = "ROI",
      all.x = TRUE
    )
    
    res$Measure <- m
    
    write.csv(
      res,
      file.path(
        sub_outdir,
        paste0(
          m,
          "_emmeans_contrasts.csv"
        )
      ),
      row.names = FALSE
    )
    
    results[[m]] <- res
  }
  
  if (length(results) == 0) {
    message("No results produced.")
    return(invisible(NULL))
  }
  
  # Combine results
  all_results <- bind_rows(results)
  
  all_results$ROI_global <- paste0(
    all_results$Measure,
    "_",
    all_results$ROI
  )
  
  # FDR across ROIs within each Measure and contrast combination
  all_results <- all_results %>%
    group_by(Measure, contrast) %>%
    mutate(
      p_fdr_roi = p.adjust(
        p_value,
        method = "BH"
      ),
      significant_roi_fdr = p_fdr_roi < alpha
    ) %>%
    ungroup()
  
  # Save ROI-level FDR results
  write.csv(
    all_results,
    file.path(
      sub_outdir,
      paste0(
        "ROI_FDR_",
        analysis_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  # Significant results
  sig <- all_results %>%
    filter(significant_roi_fdr)
  
  if (nrow(sig) == 0) {
    
    message(
      "No ROIs survive ROI-level FDR for ",
      analysis_name,
      "."
    )
    
  } else {
    
    message(
      "\nSignificant results for ",
      analysis_name,
      ":"
    )
    
    print(
      sig[
        ,
        c(
          "ROI_global",
          "contrast",
          "estimate",
          "p_value",
          "p_fdr_roi"
        )
      ]
    )
    
    # Plots
    for (i in seq_len(nrow(sig))) {
      
      result <- sig[i, ]
      
      plot_data <- data[
        data$Measure == result$Measure &
          data$ROI == result$ROI,
      ]
      
      plot_data$Group <- factor(
        plot_data[[group_col]],
        levels = group_levels
      )
      
      p_label <- if (
        result$p_fdr_roi < 0.001
      ) {
        "p_fdr < 0.001"
      } else {
        paste0(
          "p_fdr = ",
          signif(result$p_fdr_roi, 3)
        )
      }
      
      plot <- ggplot(
        plot_data,
        aes(
          x = Group,
          y = QSM,
          fill = Group
        )
      ) +
        geom_boxplot(
          outlier.shape = NA,
          alpha = 0.6,
          width = 0.5
        ) +
        geom_jitter(
          width = 0.12,
          alpha = 0.6,
          size = 1.8
        ) +
        labs(
          title = paste0(
            result$Measure,
            " - ",
            result$ROI
          ),
          subtitle = paste0(
            result$contrast,
            " | ",
            p_label
          ),
          x = NULL,
          y = "QSM (mean ROI value)"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "bold")
        )
      
      filename <- paste0(
        result$Measure,
        "_",
        gsub(
          "[^A-Za-z0-9]+",
          "-",
          result$ROI
        ),
        "_",
        gsub(
          "[^A-Za-z0-9]+",
          "-",
          result$contrast
        ),
        ".png"
      )
      
      ggsave(
        file.path(
          plot_dir,
          filename
        ),
        plot,
        width = 5,
        height = 5,
        dpi = 300
      )
    }
  }
  
  invisible(all_results)
}

# Run the three analyses

res_full6 <- run_group_analysis(
  master_long,
  group_col = "Diagnosis",
  group_levels = DIAG_LEVELS,
  ref_level = "HC",
  analysis_name = "Full6Group",
  outdir = OUTDIR,
  alpha = ALPHA
)

res_three <- run_group_analysis(
  master_long,
  group_col = "Group3",
  group_levels = c("HC", "PD", "Cognitive"),
  ref_level = "HC",
  analysis_name = "ThreeGroup",
  outdir = OUTDIR,
  alpha = ALPHA
)

res_two <- run_group_analysis(
  master_long,
  group_col = "Group2",
  group_levels = c("Control", "Affected"),
  ref_level = "Control",
  analysis_name = "TwoGroup",
  outdir = OUTDIR,
  alpha = ALPHA
)

# Demographic summaries
write_demo_summary <- function(
    data,
    group_col,
    out_path
) {
  
  summary <- data %>%
    group_by(.data[[group_col]]) %>%
    summarise(
      n = n(),
      age_mean = mean(
        as.numeric(Age),
        na.rm = TRUE
      ),
      age_sd = sd(
        as.numeric(Age),
        na.rm = TRUE
      ),
      n_female = sum(
        Sex == "F",
        na.rm = TRUE
      ),
      n_male = sum(
        Sex == "M",
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  write.csv(
    summary,
    out_path,
    row.names = FALSE
  )
  
  print(summary)
  
  summary
}

write_demo_summary(
  diag,
  "Diagnosis",
  file.path(
    OUTDIR,
    "demographic_summary_Full6Group.csv"
  )
)

write_demo_summary(
  diag,
  "Group3",
  file.path(
    OUTDIR,
    "demographic_summary_ThreeGroup.csv"
  )
)

write_demo_summary(
  diag,
  "Group2",
  file.path(
    OUTDIR,
    "demographic_summary_TwoGroup.csv"
  )
)

summarise_significant <- function(results) {
  
  if (is.null(results) || nrow(results) == 0) {
    return(
      data.frame(
        Measure = character(),
        contrast = character(),
        direction = character(),
        n_significant = integer()
      )
    )
  }
  
  results %>%
    filter(significant_roi_fdr) %>%
    mutate(
      direction = ifelse(
        estimate > 0,
        "Higher",
        "Lower"
      )
    ) %>%
    group_by(Measure, contrast, direction) %>%
    summarise(
      n_significant = n(),
      .groups = "drop"
    )
}

summary_full6 <- summarise_significant(res_full6)
summary_three <- summarise_significant(res_three)
summary_two <- summarise_significant(res_two)

write.csv(
  summary_full6,
  file.path(OUTDIR, "summary_significant_Full6.csv"),
  row.names = FALSE
)

write.csv(
  summary_three,
  file.path(OUTDIR, "summary_significant_ThreeGroup.csv"),
  row.names = FALSE
)

write.csv(
  summary_two,
  file.path(OUTDIR, "summary_significant_TwoGroup.csv"),
  row.names = FALSE
)

print(summary_full6)
print(summary_three)
print(summary_two)
