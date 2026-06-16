# 0. Environment Initialization -----------------------------------------------------------
rm(list = ls())
gc()

setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code//")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")

dir.create("./23_Multi_Model_Prognostic_Analysis/", showWarnings = FALSE, recursive = TRUE)


# TCGA training cohort; GSE165656 external validation cohort


############################################################
## Multi-model survival modeling for TCGA-LAML and GSE165656
## Corrected version:
## - TCGA as training cohort
## - GSE165656 as external validation cohort
## - No ComBat leakage
## - Validation set standardized by TCGA mean/sd
## - Multiple feature selection + modeling combinations
## - Best model selected by GSE165656 C-index
############################################################



############################################################
## 0. Load Required Packages
############################################################

pkg_list <- c(
  "survival",
  "survminer",
  "glmnet",
  "randomForestSRC",
  "CoxBoost",
  "gbm",
  "survivalsvm",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "ggplot2",
  "pheatmap"
)

for (pkg in pkg_list) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(survival)
library(survminer)
library(glmnet)
library(randomForestSRC)
library(CoxBoost)
library(gbm)
library(survivalsvm)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(pheatmap)

set.seed(123)

############################################################
## 1. File Path Configuration
############################################################

tcga_file <- "../TCGA_data/TCGA_LAML_879_feature_expr_with_ostime_new_136_samples.csv"
gse_file  <- "../TCGA_data/GSE165656_LAML_879_feature_expr_with_ostime_new.csv"

univ_cox_file <- "../Outdata/7.Prognostic_analysis/2_univariate_cox_significant_genes.csv"

outdata_dir <- "./23_Multi_Model_Prognostic_Analysis/outdata"
outplot_dir <- "./23_Multi_Model_Prognostic_Analysis//outplot"

dir.create(outdata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outplot_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 2. Import and Clean Raw Data
############################################################

cat("\n==============================\n")
cat("Step 1: Reading data\n")
cat("==============================\n")

tcga <- read.csv(
  tcga_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

gse <- read.csv(
  gse_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

colnames(tcga)[colnames(tcga) == "sample"] <- "ID"
colnames(gse)[colnames(gse) == "sample"] <- "ID"

if (!all(c("ID", "OS.time", "OS") %in% colnames(tcga))) {
  stop("TCGA dataset must contain three columns: ID, OS.time, OS.")
}

if (!all(c("ID", "OS.time", "OS") %in% colnames(gse))) {
  stop("GSE165656 dataset must contain three columns: ID, OS.time, OS.")
}

tcga$OS.time <- as.numeric(as.character(tcga$OS.time))
tcga$OS      <- as.numeric(as.character(tcga$OS))

gse$OS.time <- as.numeric(as.character(gse$OS.time))
gse$OS      <- as.numeric(as.character(gse$OS))

## Note:
## If OS.time in GSE165656 file is already in days, comment out the next line.
## If OS.time in GSE165656 file is in weeks, retain the next line.
gse$OS.time <- gse$OS.time * 7

tcga <- tcga[!is.na(tcga$OS.time) & !is.na(tcga$OS), ]
gse  <- gse[!is.na(gse$OS.time) & !is.na(gse$OS), ]

tcga$OS.time[tcga$OS.time <= 0] <- 1
gse$OS.time[gse$OS.time <= 0]   <- 1

tcga$OS <- ifelse(tcga$OS > 0, 1, 0)
gse$OS  <- ifelse(gse$OS > 0, 1, 0)

cat("TCGA samples:", nrow(tcga), "\n")
cat("GSE165656 samples:", nrow(gse), "\n")

############################################################
## 3. Identify Shared Genes and Candidate Predictive Genes
############################################################

cat("\n==============================\n")
cat("Step 2: Preparing gene list\n")
cat("==============================\n")

clin_cols <- c("ID", "OS.time", "OS")

tcga_genes <- setdiff(colnames(tcga), clin_cols)
gse_genes  <- setdiff(colnames(gse), clin_cols)

common_genes <- intersect(tcga_genes, gse_genes)

cat("TCGA gene count:", length(tcga_genes), "\n")
cat("GSE gene count:", length(gse_genes), "\n")
cat("Shared gene count:", length(common_genes), "\n")

## If univariate Cox significant gene file exists, further filter candidates;
## If not, use all shared genes as candidates.
if (file.exists(univ_cox_file)) {
  
  univ_genes_df <- read.csv(
    univ_cox_file,
    header = TRUE,
    check.names = FALSE
  )
  
  if ("gene" %in% colnames(univ_genes_df)) {
    univ_genes <- unique(univ_genes_df$gene)
  } else if ("gene_name" %in% colnames(univ_genes_df)) {
    univ_genes <- unique(univ_genes_df$gene_name)
  } else {
    stop("Column 'gene' or 'gene_name' not found in univariate Cox result file.")
  }
  
  model_genes <- intersect(common_genes, univ_genes)
  
  cat("Using genes with significant univariate Cox P-values.\n")
  cat("Candidate modeling gene count:", length(model_genes), "\n")
  
} else {
  
  model_genes <- common_genes
  
  cat("Univariate Cox result file missing. All shared genes will be used.\n")
  cat("Candidate modeling gene count:", length(model_genes), "\n")
}

if (length(model_genes) < 2) {
  stop("Fewer than 2 candidate modeling genes detected, please verify gene names.")
}

tcga <- tcga[, c(clin_cols, model_genes)]
gse  <- gse[, c(clin_cols, model_genes)]

tcga[, model_genes] <- lapply(
  tcga[, model_genes],
  function(x) as.numeric(as.character(x))
)

gse[, model_genes] <- lapply(
  gse[, model_genes],
  function(x) as.numeric(as.character(x))
)

na_genes <- model_genes[
  colSums(is.na(tcga[, model_genes, drop = FALSE])) > 0 |
    colSums(is.na(gse[, model_genes, drop = FALSE])) > 0
]

if (length(na_genes) > 0) {
  cat("Remove genes containing NA values:", length(na_genes), "\n")
  model_genes <- setdiff(model_genes, na_genes)
}

if (length(model_genes) < 2) {
  stop("Fewer than 2 candidate genes remaining after removing genes with NA values.")
}

tcga <- tcga[, c(clin_cols, model_genes)]
gse  <- gse[, c(clin_cols, model_genes)]

############################################################
## 4. Standardization Using Training Cohort Mean & SD
############################################################

cat("\n==============================\n")
cat("Step 3: Feature standardization\n")
cat("==============================\n")

x_train_raw <- as.matrix(tcga[, model_genes, drop = FALSE])
x_valid_raw <- as.matrix(gse[, model_genes, drop = FALSE])

train_mean <- apply(x_train_raw, 2, mean, na.rm = TRUE)
train_sd   <- apply(x_train_raw, 2, sd, na.rm = TRUE)

keep_genes <- names(train_sd)[train_sd > 0 & !is.na(train_sd)]

x_train_raw <- x_train_raw[, keep_genes, drop = FALSE]
x_valid_raw <- x_valid_raw[, keep_genes, drop = FALSE]

train_mean <- train_mean[keep_genes]
train_sd   <- train_sd[keep_genes]

x_train <- scale(x_train_raw, center = train_mean, scale = train_sd)
x_valid <- scale(x_valid_raw, center = train_mean, scale = train_sd)

est_dd <- data.frame(
  OS.time = tcga$OS.time,
  OS = tcga$OS,
  x_train,
  check.names = FALSE
)

val_dd_list <- list(
  TCGA = est_dd,
  GSE165656 = data.frame(
    OS.time = gse$OS.time,
    OS = gse$OS,
    x_valid,
    check.names = FALSE
  )
)

pre_var <- colnames(est_dd)[-c(1, 2)]

cat("Final modeling gene count:", length(pre_var), "\n")

write_tsv(
  data.frame(gene = pre_var),
  file.path(outdata_dir, "0_modeling_candidate_genes.txt")
)

save(
  train_mean,
  train_sd,
  pre_var,
  file = file.path(outdata_dir, "0_train_standardization_parameters.RData")
)

############################################################
## 5. Utility Helper Functions
############################################################

cat("\n==============================\n")
cat("Step 4: Define utility functions\n")
cat("==============================\n")

safe_model <- function(model_name, expr) {
  
  cat("\nRunning:", model_name, "\n")
  
  warn_msg <- character()
  
  out <- tryCatch(
    {
      withCallingHandlers(
        expr,
        warning = function(w) {
          warn_msg <<- c(warn_msg, conditionMessage(w))
          cat("Warning triggered in:", model_name, "\n")
          cat("Warning message:", conditionMessage(w), "\n")
          invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) {
      cat("Model failed:", model_name, "\n")
      cat("Error message:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  return(out)
}

calc_cindex <- function(rs_df) {
  
  rs_df <- rs_df[!is.na(rs_df$RS), ]
  
  if (nrow(rs_df) < 5) return(NA)
  if (length(unique(rs_df$RS)) <= 1) return(NA)
  if (length(unique(rs_df$OS)) <= 1) return(NA)
  
  fit <- tryCatch(
    coxph(Surv(OS.time, OS) ~ RS, data = rs_df),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NA)
  
  as.numeric(summary(fit)$concordance[1])
}

add_result <- function(result, model_name, rs_list) {
  
  tmp <- data.frame(
    ID = names(rs_list),
    Cindex = sapply(rs_list, calc_cindex),
    Model = model_name,
    stringsAsFactors = FALSE
  )
  
  result <- rbind(result, tmp)
  
  return(result)
}

############################################################
## 6. Feature Selection Functions
############################################################

############################################################
## Modification 1: Maximum number of genes retained after filtering
############################################################

max_feature_genes <- 30

select_all <- function(est_dd, pre_var) {
  pre_var
}

############################################################
## Modification 2: RSF feature selection
## Retain only top max_feature_genes genes ranked by VIMP
############################################################

select_rsf <- function(est_dd, pre_var, seed = 123,
                       rf_nodesize = 5,
                       max_feature_genes = 30) {
  
  set.seed(seed)
  
  fit <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = est_dd[, c("OS.time", "OS", pre_var), drop = FALSE],
    ntree = 1000,
    nodesize = rf_nodesize,
    splitrule = "logrank",
    importance = TRUE,
    forest = TRUE,
    seed = seed
  )
  
  vimp <- fit$importance
  
  vimp <- vimp[!is.na(vimp)]
  
  if (length(vimp) < 2) {
    stop("Fewer than 2 genes with valid RSF importance scores.")
  }
  
  ## Prioritize genes with positive VIMP values
  vimp_pos <- vimp[vimp > 0]
  
  if (length(vimp_pos) >= 2) {
    vimp_use <- sort(vimp_pos, decreasing = TRUE)
  } else {
    vimp_use <- sort(vimp, decreasing = TRUE)
  }
  
  rid <- names(vimp_use)[1:min(max_feature_genes, length(vimp_use))]
  
  rid <- rid[!is.na(rid)]
  
  if (length(rid) < 2) {
    stop("RSF selected fewer than 2 genes.")
  }
  
  return(rid)
}

############################################################
## Modification 3: LASSO feature selection
## Retain only top max_feature_genes genes ranked by absolute coefficient
############################################################

select_lasso <- function(est_dd, pre_var, seed = 123,
                         max_feature_genes = 30) {
  
  x <- as.matrix(est_dd[, pre_var, drop = FALSE])
  y <- Surv(est_dd$OS.time, est_dd$OS)
  
  set.seed(seed)
  
  fit <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    nfolds = 10,
    type.measure = "deviance",
    maxit = 1e6
  )
  
  coef_mat <- coef(fit, s = "lambda.min")
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  
  coef_use <- coef_vec[coef_vec != 0]
  
  if (length(coef_use) < 2) {
    
    coef_mat <- coef(fit, s = "lambda.1se")
    coef_vec <- as.numeric(coef_mat)
    names(coef_vec) <- rownames(coef_mat)
    
    coef_use <- coef_vec[coef_vec != 0]
  }
  
  if (length(coef_use) < 2) {
    stop("LASSO selected fewer than 2 genes.")
  }
  
  coef_use <- sort(abs(coef_use), decreasing = TRUE)
  
  rid <- names(coef_use)[1:min(max_feature_genes, length(coef_use))]
  
  return(rid)
}

############################################################
## Modification 4: CoxBoost feature selection
## Retain only top max_feature_genes genes ranked by absolute coefficient
############################################################

select_coxboost <- function(est_dd, pre_var, seed = 123,
                            max_feature_genes = 30) {
  
  set.seed(seed)
  
  x <- as.matrix(est_dd[, pre_var, drop = FALSE])
  
  pen <- optimCoxBoostPenalty(
    time = est_dd$OS.time,
    status = est_dd$OS,
    x = x,
    trace = FALSE,
    start.penalty = 500,
    parallel = FALSE
  )
  
  cv.res <- cv.CoxBoost(
    time = est_dd$OS.time,
    status = est_dd$OS,
    x = x,
    maxstepno = 500,
    K = 10,
    type = "verweij",
    penalty = pen$penalty
  )
  
  fit <- CoxBoost(
    time = est_dd$OS.time,
    status = est_dd$OS,
    x = x,
    stepno = cv.res$optimal.step,
    penalty = pen$penalty
  )
  
  beta <- as.numeric(coef(fit))
  names(beta) <- pre_var
  
  beta_use <- beta[beta != 0]
  
  if (length(beta_use) < 2) {
    stop("CoxBoost selected fewer than 2 genes.")
  }
  
  beta_use <- sort(abs(beta_use), decreasing = TRUE)
  
  rid <- names(beta_use)[1:min(max_feature_genes, length(beta_use))]
  
  return(rid)
}

select_stepcox_backward <- function(est_dd, pre_var, max_genes_for_step = 80) {
  
  if (length(pre_var) > max_genes_for_step) {
    stop(
      paste0(
        "StepCox skipped: total gene count exceeds ",
        max_genes_for_step,
        ". Apply feature filtering prior to StepCox."
      )
    )
  }
  
  dd <- est_dd[, c("OS.time", "OS", pre_var), drop = FALSE]
  
  fit0 <- coxph(Surv(OS.time, OS) ~ ., data = dd)
  
  fit <- step(
    fit0,
    direction = "backward",
    trace = 0
  )
  
  rid <- names(coef(fit))
  
  if (length(rid) < 2) {
    stop("StepCox selected fewer than 2 genes.")
  }
  
  return(rid)
}

############################################################
## 7. Model Fitting Functions
############################################################

############################################################
## Modification 1: Standard Cox model fitting function
## Constraints:
## 1. Raw full gene set (All) cannot be directly input to Cox
## 2. Only filtered small gene sets are allowed for Cox regression
############################################################

fit_cox_model <- function(est_dd, val_dd_list, rid,
                          max_genes_for_cox = 50,
                          epv_cutoff = 2) {
  
  event_num <- sum(est_dd$OS == 1, na.rm = TRUE)
  
  if (length(rid) < 2) {
    stop("Skip standard Cox: fewer than 2 candidate genes.")
  }
  
  if (length(rid) > max_genes_for_cox) {
    stop(
      paste0(
        "Skip standard Cox: filtered gene count = ", length(rid),
        ", exceeds max_genes_for_cox = ", max_genes_for_cox
      )
    )
  }
  
  if (length(rid) >= floor(event_num / epv_cutoff)) {
    stop(
      paste0(
        "Skip standard Cox: filtered gene count = ", length(rid),
        ", total events = ", event_num,
        ", variable count is too high relative to events."
      )
    )
  }
  
  dd <- est_dd[, c("OS.time", "OS", rid), drop = FALSE]
  
  fit <- coxph(
    Surv(OS.time, OS) ~ .,
    data = dd,
    control = coxph.control(iter.max = 100)
  )
  
  if (any(is.na(coef(fit)))) {
    stop("Skip standard Cox: NA coefficients detected after model fitting.")
  }
  
  rs <- lapply(val_dd_list, function(dat) {
    
    newdat <- dat[, c("OS.time", "OS", rid), drop = FALSE]
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(
        predict(
          fit,
          newdata = newdat,
          type = "lp"
        )
      )
    )
  })
  
  list(
    fit = fit,
    rs = rs,
    genes = names(coef(fit)),
    model_type = "Cox"
  )
}

############################################################
## Modification 2: Stepwise Cox model fitting function
## Constraints:
## 1. Raw full gene set (All) cannot be directly input to StepCox
## 2. Only filtered small gene sets are allowed for StepCox
############################################################

fit_stepcox_model <- function(est_dd, val_dd_list, rid,
                              direction = "backward",
                              max_genes_for_stepcox = 50,
                              epv_cutoff = 2) {
  
  event_num <- sum(est_dd$OS == 1, na.rm = TRUE)
  
  if (length(rid) < 2) {
    stop("Skip StepCox: fewer than 2 candidate genes.")
  }
  
  if (length(rid) > max_genes_for_stepcox) {
    stop(
      paste0(
        "Skip StepCox: filtered gene count = ", length(rid),
        ", exceeds max_genes_for_stepcox = ", max_genes_for_stepcox
      )
    )
  }
  
  if (length(rid) >= floor(event_num / epv_cutoff)) {
    stop(
      paste0(
        "Skip StepCox: filtered gene count = ", length(rid),
        ", total events = ", event_num,
        ", variable count is too high relative to events."
      )
    )
  }
  
  dd <- est_dd[, c("OS.time", "OS", rid), drop = FALSE]
  
  fit0 <- coxph(
    Surv(OS.time, OS) ~ .,
    data = dd,
    control = coxph.control(iter.max = 100)
  )
  
  if (any(is.na(coef(fit0)))) {
    stop("Skip StepCox: NA coefficients detected in initial Cox model.")
  }
  
  fit <- step(
    fit0,
    direction = direction,
    trace = 0
  )
  
  final_genes <- names(coef(fit))
  
  if (length(final_genes) < 1) {
    stop("StepCox final model contains no predictive genes.")
  }
  
  if (any(is.na(coef(fit)))) {
    stop("Skip StepCox: NA coefficients detected in final model.")
  }
  
  rs <- lapply(val_dd_list, function(dat) {
    
    newdat <- dat[, c("OS.time", "OS", final_genes), drop = FALSE]
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(
        predict(
          fit,
          newdata = newdat,
          type = "lp"
        )
      )
    )
  })
  
  list(
    fit = fit,
    rs = rs,
    genes = final_genes,
    model_type = "StepCox"
  )
}

fit_lasso_model <- function(est_dd, val_dd_list, rid, seed = 123) {
  
  x <- as.matrix(est_dd[, rid, drop = FALSE])
  y <- Surv(est_dd$OS.time, est_dd$OS)
  
  set.seed(seed)
  
  fit <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 1,
    nfolds = 10,
    type.measure = "deviance"
  )
  
  rs <- lapply(val_dd_list, function(dat) {
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(
        predict(
          fit,
          newx = as.matrix(dat[, rid, drop = FALSE]),
          s = "lambda.min",
          type = "link"
        )
      )
    )
  })
  
  coef_mat <- coef(fit, s = "lambda.min")
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  final_genes <- names(coef_vec)[coef_vec != 0]
  
  list(fit = fit, rs = rs, genes = final_genes, model_type = "Lasso")
}

fit_ridge_model <- function(est_dd, val_dd_list, rid, seed = 123) {
  
  x <- as.matrix(est_dd[, rid, drop = FALSE])
  y <- Surv(est_dd$OS.time, est_dd$OS)
  
  set.seed(seed)
  
  fit <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = 0,
    nfolds = 10,
    type.measure = "deviance"
  )
  
  rs <- lapply(val_dd_list, function(dat) {
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(
        predict(
          fit,
          newx = as.matrix(dat[, rid, drop = FALSE]),
          s = "lambda.min",
          type = "link"
        )
      )
    )
  })
  
  list(fit = fit, rs = rs, genes = rid, model_type = "Ridge")
}

fit_enet_model <- function(est_dd, val_dd_list, rid, alpha = 0.5, seed = 123) {
  
  x <- as.matrix(est_dd[, rid, drop = FALSE])
  y <- Surv(est_dd$OS.time, est_dd$OS)
  
  set.seed(seed)
  
  fit <- cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = alpha,
    nfolds = 10,
    type.measure = "deviance"
  )
  
  rs <- lapply(val_dd_list, function(dat) {
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(
        predict(
          fit,
          newx = as.matrix(dat[, rid, drop = FALSE]),
          s = "lambda.min",
          type = "link"
        )
      )
    )
  })
  
  coef_mat <- coef(fit, s = "lambda.min")
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  final_genes <- names(coef_vec)[coef_vec != 0]
  
  list(
    fit = fit,
    rs = rs,
    genes = final_genes,
    model_type = paste0("Enet_alpha_", alpha)
  )
}

fit_coxboost_model <- function(est_dd, val_dd_list, rid, seed = 123) {
  
  set.seed(seed)
  
  x <- as.matrix(est_dd[, rid, drop = FALSE])
  
  pen <- optimCoxBoostPenalty(
    time = est_dd$OS.time,
    status = est_dd$OS,
    x = x,
    trace = FALSE,
    start.penalty = 500,
    parallel = FALSE
  )
  
  cv.res <- cv.CoxBoost(
    time = est_dd$OS.time,
    status = est_dd$OS,
    x = x,
    maxstepno = 500,
    K = 10,
    type = "verweij",
    penalty = pen$penalty
  )
  
  fit <- CoxBoost(
    time = est_dd$OS.time,
    status = est_dd$OS,
    x = x,
    stepno = cv.res$optimal.step,
    penalty = pen$penalty
  )
  
  rs <- lapply(val_dd_list, function(dat) {
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(
        predict(
          fit,
          newdata = as.matrix(dat[, rid, drop = FALSE]),
          newtime = dat$OS.time,
          newstatus = dat$OS,
          type = "lp"
        )
      )
    )
  })
  
  beta <- as.numeric(coef(fit))
  names(beta) <- rid
  final_genes <- names(beta)[beta != 0]
  
  list(fit = fit, rs = rs, genes = final_genes, model_type = "CoxBoost")
}

fit_rsf_model <- function(est_dd, val_dd_list, rid, seed = 123, rf_nodesize = 5) {
  
  set.seed(seed)
  
  fit <- rfsrc(
    Surv(OS.time, OS) ~ .,
    data = est_dd[, c("OS.time", "OS", rid), drop = FALSE],
    ntree = 1000,
    nodesize = rf_nodesize,
    splitrule = "logrank",
    importance = TRUE,
    forest = TRUE,
    seed = seed
  )
  
  rs <- lapply(val_dd_list, function(dat) {
    
    newdat <- dat[, c("OS.time", "OS", rid), drop = FALSE]
    
    pred <- predict(fit, newdata = newdat)
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(pred$predicted)
    )
  })
  
  list(fit = fit, rs = rs, genes = rid, model_type = "RSF")
}

fit_gbm_model <- function(est_dd, val_dd_list, rid, seed = 123) {
  
  dd <- est_dd[, c("OS.time", "OS", rid), drop = FALSE]
  
  set.seed(seed)
  
  fit0 <- gbm(
    formula = Surv(OS.time, OS) ~ .,
    data = dd,
    distribution = "coxph",
    n.trees = 10000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10,
    n.cores = 1,
    verbose = FALSE
  )
  
  best_trees <- gbm.perf(fit0, method = "cv", plot.it = FALSE)
  
  if (is.na(best_trees) || best_trees < 1) {
    best_trees <- which.min(fit0$cv.error)
  }
  
  rs <- lapply(val_dd_list, function(dat) {
    
    newdat <- dat[, c("OS.time", "OS", rid), drop = FALSE]
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(
        predict(
          fit0,
          newdata = newdat,
          n.trees = best_trees,
          type = "link"
        )
      )
    )
  })
  
  list(
    fit = fit0,
    rs = rs,
    genes = rid,
    model_type = "GBM",
    best_trees = best_trees
  )
}

fit_svm_model <- function(est_dd, val_dd_list, rid, seed = 123) {
  
  dd <- est_dd[, c("OS.time", "OS", rid), drop = FALSE]
  
  set.seed(seed)
  
  fit <- survivalsvm(
    Surv(OS.time, OS) ~ .,
    data = dd,
    gamma.mu = 1
  )
  
  rs <- lapply(val_dd_list, function(dat) {
    
    newdat <- dat[, c("OS.time", "OS", rid), drop = FALSE]
    
    pred <- predict(fit, newdat)
    
    data.frame(
      OS.time = dat$OS.time,
      OS = dat$OS,
      RS = as.numeric(pred$predicted)
    )
  })
  
  list(fit = fit, rs = rs, genes = rid, model_type = "survivalSVM")
}

############################################################
## 8. Define Feature Selection & Model Training Pipelines
############################################################

cat("\n==============================\n")
cat("Step 5: Define all model combinations\n")
cat("==============================\n")

############################################################
## Modification 5: All feature selection pipelines limited to max_feature_genes genes
############################################################

max_feature_genes <- 30

feature_selectors <- list(
  
  All = function() {
    select_all(est_dd, pre_var)
  },
  
  RSF = function() {
    select_rsf(
      est_dd,
      pre_var,
      seed = 123,
      rf_nodesize = 5,
      max_feature_genes = max_feature_genes
    )
  },
  
  Lasso = function() {
    select_lasso(
      est_dd,
      pre_var,
      seed = 123,
      max_feature_genes = max_feature_genes
    )
  },
  
  CoxBoost = function() {
    select_coxboost(
      est_dd,
      pre_var,
      seed = 123,
      max_feature_genes = max_feature_genes
    )
  }
)

## Enable StepCox backward selection only if total candidate genes ≤ 80
if (length(pre_var) <= 80) {
  feature_selectors$StepCox_backward <- function() {
    select_stepcox_backward(est_dd, pre_var, max_genes_for_step = 80)
  }
}

model_builders <- list(
  ############################################################
  ## Modification 3: Cox / StepCox only accept pre-filtered small gene sets
  ############################################################
  
  Cox = function(rid) {
    fit_cox_model(
      est_dd,
      val_dd_list,
      rid,
      max_genes_for_cox = 50,
      epv_cutoff = 2
    )
  },
  
  StepCox_backward = function(rid) {
    fit_stepcox_model(
      est_dd,
      val_dd_list,
      rid,
      direction = "backward",
      max_genes_for_stepcox = 50,
      epv_cutoff = 2
    )
  },
  
  Lasso = function(rid) {
    fit_lasso_model(est_dd, val_dd_list, rid, seed = 123)
  },
  Ridge = function(rid) {
    fit_ridge_model(est_dd, val_dd_list, rid, seed = 123)
  },
  CoxBoost = function(rid) {
    fit_coxboost_model(est_dd, val_dd_list, rid, seed = 123)
  },
  RSF = function(rid) {
    fit_rsf_model(est_dd, val_dd_list, rid, seed = 123, rf_nodesize = 5)
  },
  GBM = function(rid) {
    fit_gbm_model(est_dd, val_dd_list, rid, seed = 123)
  },
  survivalSVM = function(rid) {
    fit_svm_model(est_dd, val_dd_list, rid, seed = 123)
  }
)

## Add Elastic Net models with alpha ranging 0.1 ~ 0.9
for (a in seq(0.1, 0.9, by = 0.1)) {
  model_builders[[paste0("Enet_alpha_", a)]] <- local({
    alpha_value <- a
    function(rid) {
      fit_enet_model(
        est_dd,
        val_dd_list,
        rid,
        alpha = alpha_value,
        seed = 123
      )
    }
  })
}

############################################################
## 9. Execute All Feature Selection + Model Combinations
############################################################

############################################################
## Modification 4: Run all model combinations with filtering rules
## Core logic:
## 1. All + Cox directly skipped
## 2. All + StepCox directly skipped
## 3. Cox/StepCox only allowed after RSF/Lasso/CoxBoost feature filtering
############################################################

cat("\n==============================\n")
cat("Step 6: Execute all model combinations\n")
cat("==============================\n")

result <- data.frame()
model_store <- list()
gene_store <- list()

max_genes_for_cox_models <- 50
event_num <- sum(est_dd$OS == 1, na.rm = TRUE)

for (fs_name in names(feature_selectors)) {
  
  cat("\n###################################################\n")
  cat("Feature selection pipeline:", fs_name, "\n")
  cat("###################################################\n")
  
  rid <- safe_model(
    paste0("FeatureSelect_", fs_name),
    feature_selectors[fs_name]]()
  )
  
  if (is.null(rid)) next
  
  rid <- intersect(rid, pre_var)
  rid <- unique(rid)
  
  if (length(rid) < 2) {
    cat("Skipping pipeline", fs_name, ": fewer than 2 selected genes\n")
    next
  }
  
  cat("Number of selected genes:", length(rid), "\n")
  
  gene_store[[paste0("FeatureSelect_", fs_name)]] <- rid
  
  write_tsv(
    data.frame(gene = rid),
    file.path(outdata_dir, paste0("genes_selected_by_", fs_name, ".txt"))
  )
  
  for (model_name in names(model_builders)) {
    
    full_model_name <- paste0(fs_name, " + ", model_name)
    
    if (fs_name == "All") {
      full_model_name <- model_name
    }
    
    ################################################
    ## Modification 4.1:
    ## "All" is not a formal feature selection method.
    ## Skip All + Cox and All + StepCox entirely.
    ################################################
    
    if (fs_name == "All" &&
        model_name %in% c("Cox", "StepCox_backward")) {
      
      cat(
        "Skip model:", full_model_name,
        "| Reason: Cox/StepCox requires prior feature filtering, cannot use full gene set directly.\n"
      )
      
      next
    }
    
    ################################################
    ## Modification 4.2:
    ## Skip Cox/StepCox even after filtering if gene count exceeds threshold
    ################################################
    
    if (model_name %in% c("Cox", "StepCox_backward") &&
        length(rid) > max_genes_for_cox_models) {
      
      cat(
        "Skip model:", full_model_name,
        "| Selected gene count =", length(rid),
        "> max_genes_for_cox_models =", max_genes_for_cox_models,
        "\n"
      )
      
      next
    }
    
    ################################################
    ## Modification 4.3:
    ## Skip Cox/StepCox if variable count is too high relative to event count
    ################################################
    
    if (model_name %in% c("Cox", "StepCox_backward") &&
        length(rid) >= floor(event_num / 2)) {
      
      cat(
        "Skip model:", full_model_name,
        "| Selected gene count =", length(rid),
        ">= event_num/2 =", floor(event_num / 2),
        "\n"
      )
      
      next
    }
    
    fit_obj <- safe_model(
      full_model_name,
      model_builders[[model_name]](rid)
    )
    
    if (is.null(fit_obj)) next
    
    result <- add_result(
      result = result,
      model_name = full_model_name,
      rs_list = fit_obj$rs
    )
    
    model_store[[full_model_name]] <- fit_obj
    
    cat("Completed model:", full_model_name, "\n")
  }
}

if (nrow(result) == 0) {
  stop("All model runs failed. Please check input data and parameters.")
}

write.csv(
  result,
  file.path(outdata_dir, "0_all_model_Cindex_long_corrected.csv"),
  row.names = FALSE
)

save(
  result,
  model_store,
  gene_store,
  file = file.path(outdata_dir, "0_all_model_objects_corrected.RData")
)

if (nrow(result) == 0) {
  stop("All model runs failed. Please check packages and input data.")
}

write.csv(
  result,
  file.path(outdata_dir, "0_all_model_Cindex_long_after_feature_selection_Cox.csv"),
  row.names = FALSE
)

save(
  result,
  model_store,
  gene_store,
  file = file.path(outdata_dir, "0_all_model_objects_after_feature_selection_Cox.RData")
)

############################################################
## 10. Summarize C-index Results
############################################################

cat("\n==============================\n")
cat("Step 7: Summarize C-index performance\n")
cat("==============================\n")

result2 <- result %>%
  group_by(Model, ID) %>%
  summarise(Cindex = mean(Cindex, na.rm = TRUE), .groups = "drop")

cindex_wide <- result2 %>%
  pivot_wider(names_from = ID, values_from = Cindex) %>%
  as.data.frame()

num_cols <- setdiff(colnames(cindex_wide), "Model")
cindex_wide[, num_cols] <- lapply(cindex_wide[, num_cols, drop = FALSE], as.numeric)

if (!"TCGA" %in% colnames(cindex_wide)) {
  cindex_wide$TCGA <- NA
}

if (!"GSE165656" %in% colnames(cindex_wide)) {
  cindex_wide$GSE165656 <- NA
}

cindex_wide$Mean_Cindex <- rowMeans(
  cindex_wide[, c("TCGA", "GSE165656"), drop = FALSE],
  na.rm = TRUE
)

cindex_wide <- cindex_wide[
  order(cindex_wide$GSE165656, cindex_wide$Mean_Cindex, decreasing = TRUE),
]

write.csv(
  cindex_wide,
  file.path(outdata_dir, "1_output_C_index_corrected.csv"),
  row.names = FALSE
)

cat("\nTop 20 models ranked by GSE165656 C-index:\n")
print(head(cindex_wide, 20))

############################################################
## 11. Generate C-index Heatmap
############################################################

cat("\n==============================\n")
cat("Step 8: Plot C-index heatmap\n")
cat("==============================\n")

cindex_plot <- cindex_wide[
  !is.na(cindex_wide$TCGA) & !is.na(cindex_wide$GSE165656),
]

############################################################
## Heatmap modification details:
## 1. Add Average_Cindex = mean(TCGA C-index, GSE165656 C-index)
## 2. Sort models descending by Average_Cindex
## 3. Remove redundant self-matching combinations (e.g. RSF + RSF, Lasso + Lasso)
## 4. Apply SCI journal-friendly color palette
############################################################

cindex_plot <- subset(
  cindex_plot,
  GSE165656 >= 0.50
)
## Remove identical feature-selection + model pairs, e.g. RSF + RSF
same_model_rows <- grepl("^(.+) \\+ \\1$", cindex_plot$Model)

if (any(same_model_rows)) {
  cat("Removing redundant self-matching model combinations:\n")
  print(cindex_plot$Model[same_model_rows])
}

cindex_plot <- cindex_plot[!same_model_rows, ]

## Calculate average C-index across two cohorts
cindex_plot$Average_Cindex <- rowMeans(
  cindex_plot[, c("TCGA", "GSE165656"), drop = FALSE],
  na.rm = TRUE
)

## Sort models from highest to lowest average C-index
cindex_plot <- cindex_plot[
  order(cindex_plot$Average_Cindex, decreasing = TRUE),
]

## Export heatmap source data
write.csv(
  cindex_plot,
  file.path(outdata_dir, "1_Cindex_heatmap_data_sorted_by_average.csv"),
  row.names = FALSE
)

if (nrow(cindex_plot) > 1) {
  
  dt <- as.data.frame(
    cindex_plot[, c("TCGA", "GSE165656", "Average_Cindex"), drop = FALSE]
  )
  
  rownames(dt) <- cindex_plot$Model
  
  ## Clean column labels for figure
  colnames(dt) <- c("TCGA", "GSE165656", "Average")
  
  ## SCI standard color scheme: blue-white-red (low value blue, high value red)
  heat_cols <- colorRampPalette(
    c("#3953A4", "#FFFFFF", "#ED2123")
  )(100)
  
  ## Dynamically adjust figure height based on total model count
  plot_height <- max(6, nrow(dt) * 0.20)
  
  pdf(
    file.path(outplot_dir, "1_Cindex_heatmap_sorted_by_average_SCI2.pdf"),
    width = 5,
    height = plot_height
  )
  
  pheatmap(
    dt,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.3f",
    fontsize = 9,
    fontsize_row = 7,
    fontsize_col = 10,
    angle_col = 45,
    border_color = "grey90",
    color = heat_cols,
    breaks = seq(
      min(dt, na.rm = TRUE),
      max(dt, na.rm = TRUE),
      length.out = 101
    ),
    main = "C-index of Survival Prediction Models"
  )
  
  dev.off()
}

############################################################
## 12. Identify Optimal Model
############################################################

cat("\n==============================\n")
cat("Step 9: Select best-performing model\n")
cat("==============================\n")

valid_rank <- cindex_wide[!is.na(cindex_wide$GSE165656), ]

if (nrow(valid_rank) == 0) {
  stop("No valid model produced calculable GSE165656 C-index.")
}

best_model_name <- valid_rank$Model[1]

cat("Best model selected based on external GSE165656 C-index:", best_model_name, "\n")

writeLines(
  best_model_name,
  file.path(outdata_dir, "2_best_model_name.txt")
)

best_obj <- model_store[[best_model_name]]

if (is.null(best_obj)) {
  stop("Failed to retrieve best model object from stored model list.")
}

save(
  best_obj,
  best_model_name,
  file = file.path(outdata_dir, "2_best_model_object.RData")
)

write_tsv(
  data.frame(gene = best_obj$genes),
  file.path(outdata_dir, "2_best_model_genes.txt")
)

############################################################
## 13. Export Model Coefficients (If Supported)
############################################################

cat("\n==============================\n")
cat("Step 10: Export model coefficients where available\n")
cat("==============================\n")

coef_df <- NULL

if (best_obj$model_type %in% c("Cox", "StepCox")) {
  
  fit <- best_obj$fit
  
  coef_df <- data.frame(
    gene = names(coef(fit)),
    beta = as.numeric(coef(fit))
  )
  
} else if (best_obj$model_type == "Lasso") {
  
  fit <- best_obj$fit
  
  coef_mat <- coef(fit, s = "lambda.min")
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  
  coef_vec <- coef_vec[coef_vec != 0]
  
  coef_df <- data.frame(
    gene = names(coef_vec),
    beta = as.numeric(coef_vec)
  )
  
} else if (grepl("^Enet_alpha_", best_obj$model_type)) {
  
  fit <- best_obj$fit
  
  coef_mat <- coef(fit, s = "lambda.min")
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  
  coef_vec <- coef_vec[coef_vec != 0]
  
  coef_df <- data.frame(
    gene = names(coef_vec),
    beta = as.numeric(coef_vec)
  )
  
} else if (best_obj$model_type == "Ridge") {
  
  fit <- best_obj$fit
  
  coef_mat <- coef(fit, s = "lambda.min")
  coef_vec <- as.numeric(coef_mat)
  names(coef_vec) <- rownames(coef_mat)
  
  coef_df <- data.frame(
    gene = names(coef_vec),
    beta = as.numeric(coef_vec)
  )
  
} else if (best_obj$model_type == "CoxBoost") {
  
  fit <- best_obj$fit
  
  beta <- as.numeric(coef(fit))
  names(beta) <- best_obj$genes
  
  coef_df <- data.frame(
    gene = names(beta),
    beta = beta
  )
  
  coef_df <- coef_df[coef_df$beta != 0, ]
}

if (!is.null(coef_df) && nrow(coef_df) > 0) {
  
  write_tsv(
    coef_df,
    file.path(outdata_dir, "2_best_model_coefficients.txt")
  )
  
  cat("Coefficient table exported successfully.\n")
  
} else {
  
  cat("Best model type does not output interpretable linear coefficients.\n")
}

############################################################
## 14. Generate Risk Score & Risk Group Stratification
############################################################

cat("\n==============================\n")
cat("Step 11: Calculate risk scores and stratify risk groups\n")
cat("==============================\n")

train_rs <- best_obj$rs$TCGA
valid_rs <- best_obj$rs$GSE165656

train_df <- data.frame(
  ID = tcga$ID,
  OS.time = train_rs$OS.time,
  OS = train_rs$OS,
  riskScore = train_rs$RS,
  cohort = "TCGA",
  stringsAsFactors = FALSE
)

valid_df <- data.frame(
  ID = gse$ID,
  OS.time = valid_rs$OS.time,
  OS = valid_rs$OS,
  riskScore = valid_rs$RS,
  cohort = "GSE165656",
  stringsAsFactors = FALSE
)

## Risk cutoff derived exclusively from training cohort median
risk_cutoff <- median(train_df$riskScore, na.rm = TRUE)

train_df$risk_group <- ifelse(train_df$riskScore >= risk_cutoff
