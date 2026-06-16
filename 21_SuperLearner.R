# ===================== Super Learner Diagnostic Model =====================

# ===================== 0. Environment Initialization =====================

rm(list = ls())
gc()

{
# Set working directory (modify according to your actual path)
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.Plot_Code/")
cat("[Initialization] Current working directory: ", getwd(), "\n\n")
source("/home/weili/Project/AML/human/AML_combined_analyse/0.Enviroment.R")
dir.create("./21_SuperLearner/", showWarnings = FALSE, recursive = TRUE)

outdata_dir <- "./21_SuperLearner/Outdata/"
outplot_dir <- "./21_SuperLearner/Outplot/"

if (!dir.exists(outdata_dir)) {
  dir.create(outdata_dir, recursive = TRUE)
}

if (!dir.exists(outplot_dir)) {
  dir.create(outplot_dir, recursive = TRUE)
}

cat("[Output Directories]\n")
cat("Data output directory: ", outdata_dir, "\n")
cat("Figure output directory: ", outplot_dir, "\n\n")


# ===================== 1. Install & Load Required R Packages =====================

need_packages <- c(
  "data.table",
  "SuperLearner",
  "glmnet",
  "ranger",
  "randomForest",
  "xgboost",
  "pROC",
  "caret",
  "dplyr",
  "tidyr",
  "ggplot2",
  "e1071",      # SVM, Naive Bayes
  "gbm",        # Gradient Boosting Machine
  "nnet"        # Neural Network
)

for (pkg in need_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(data.table)
library(SuperLearner)
library(glmnet)
library(ranger)
library(randomForest)
library(xgboost)
library(pROC)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(e1071)
library(gbm)
library(nnet)


cat("[All R packages loaded successfully]\n\n")


# ===================== 2. Parameter Settings =====================

# ------------------------------------------------------------#
# Rules for diagnostic model labeling:
# 1 = Diseased / AML / Tumor
# 0 = Healthy / Normal
#
# By default, DISEASE_CODE = 1
# If your dataset uses 1 as disease label, keep DISEASE_CODE = 1
# ------------------------------------------------------------

DISEASE_CODE <- 1

set.seed(123)

train_file <- "../Outdata/5.all_data_harmony/5.1_all_expr_train.csv"
test_file  <- "../Outdata/5.all_data_harmony/5.2_all_expr_test.csv"
val_file   <- "../Outdata/5.all_data_harmony/5.3_all_expr_val.csv"

sample_info_file <- "../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv"

cat("[Label Configuration]\n")
cat("In original Group column, value ", DISEASE_CODE, " will be defined as positive disease label (model label = 1)\n\n")


# ===================== 3. Load Datasets =====================

cat("[Loading expression matrices]\n")

train_expr <- data.table::fread(train_file, data.table = FALSE, check.names = FALSE)
test_expr  <- data.table::fread(test_file,  data.table = FALSE, check.names = FALSE)
val_expr   <- data.table::fread(val_file,   data.table = FALSE, check.names = FALSE)

sample_info <- read.csv(
  sample_info_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("train_expr dimension: ", nrow(train_expr), " × ", ncol(train_expr), "\n")
cat("test_expr dimension: ", nrow(test_expr), " × ", ncol(test_expr), "\n")
cat("val_expr dimension: ", nrow(val_expr), " × ", ncol(val_expr), "\n\n")

cat("First 10 columns of train_expr:\n")
print(head(colnames(train_expr), 10))
cat("\n")


# ===================== 4. Data Preprocessing Functions =====================

# Clean feature names
clean_feature_names <- function(x) {
  x <- gsub("[[:space:]]+", "", x)
  x <- gsub("-", ".", x)
  x <- make.names(x, unique = TRUE)
  return(x)
}


# Clean raw expression table loaded by fread
# fread warning explanation:
# Detected 1050 column names but data has 1051 columns, auto-added column V1
# This usually indicates the first column stores sample IDs / row names.
process_expr_table <- function(dat, dataset_name = "train") {
  
  dat <- as.data.frame(dat, check.names = FALSE)
  
  # Process auto-generated V1 column from fread
  if ("V1" %in% colnames(dat)) {
    
    # If Sample column exists, V1 is duplicate row ID in most cases
    if ("Sample" %in% colnames(dat)) {
      dat$V1 <- as.character(dat$V1)
      dat$Sample <- as.character(dat$Sample)
      
      same_rate <- mean(dat$V1 == dat$Sample, na.rm = TRUE)
      
      if (!is.na(same_rate) && same_rate > 0.8) {
        dat$V1 <- NULL
        cat("[", dataset_name, "] Detected V1 highly consistent with Sample column, V1 removed.\n", sep = "")
      } else {
        colnames(dat)[colnames(dat) == "V1"] <- "RowID"
        cat("[", dataset_name, "] Detected V1 inconsistent with Sample column, renamed V1 to RowID.\n", sep = "")
      }
      
    } else {
      colnames(dat)[colnames(dat) == "V1"] <- "Sample"
      cat("[", dataset_name, "] No Sample column found, V1 renamed to Sample.\n", sep = "")
    }
  }
  
  # Check mandatory Sample column
  if (!"Sample" %in% colnames(dat)) {
    stop(paste0("[", dataset_name, "] Missing required Sample column, please check input file."))
  }
  
  # Check mandatory Group column
  if (!"Group" %in% colnames(dat)) {
    stop(paste0("[", dataset_name, "] Missing required Group column, please check input file."))
  }
  
  # Convert sample ID to character format
  dat$Sample <- as.character(dat$Sample)
  
  # Convert group label to character to avoid factor bugs
  dat$Group <- as.character(dat$Group)
  
  # Clean all feature column names (exclude metadata columns)
  feature_cols <- setdiff(colnames(dat), c("Sample", "Group", "RowID"))
  cleaned_feature_cols <- clean_feature_names(feature_cols)
  
  colnames(dat)[match(feature_cols, colnames(dat))] <- cleaned_feature_cols
  
  # Tag dataset source
  dat$Dataset <- dataset_name
  
  return(dat)
}


train_dat <- process_expr_table(train_expr, "train")
test_dat  <- process_expr_table(test_expr,  "test")
val_dat   <- process_expr_table(val_expr,   "validation")

cat("\n[Dimensions after preprocessing]\n")
cat("train_dat: ", nrow(train_dat), " × ", ncol(train_dat), "\n")
cat("test_dat: ", nrow(test_dat), " × ", ncol(test_dat), "\n")
cat("val_dat: ", nrow(val_dat), " × ", ncol(val_dat), "\n\n")


# ===================== 5. Recode Group to Binary Label =====================

# Output label rule:
# 1 = Diseased / AML / Tumor
# 0 = Healthy / Normal
recode_group_to_label <- function(group, disease_code = 0) {
  group <- as.character(group)
  disease_code <- as.character(disease_code)
  label <- ifelse(group == disease_code, 1, 0)
  label <- as.numeric(label)
  return(label)
}

train_dat$label <- recode_group_to_label(train_dat$Group, DISEASE_CODE)
test_dat$label  <- recode_group_to_label(test_dat$Group,  DISEASE_CODE)
val_dat$label   <- recode_group_to_label(val_dat$Group,   DISEASE_CODE)

cat("[Original Group Distribution]\n")
cat("train:\n")
print(table(train_dat$Group, useNA = "ifany"))
cat("test:\n")
print(table(test_dat$Group, useNA = "ifany"))
cat("validation:\n")
print(table(val_dat$Group, useNA = "ifany"))

cat("\n[Binary label distribution after recoding: 0=Normal, 1=Disease]\n")
cat("train:\n")
print(table(train_dat$label, useNA = "ifany"))
cat("test:\n")
print(table(test_dat$label, useNA = "ifany"))
cat("validation:\n")
print(table(val_dat$label, useNA = "ifany"))
cat("\n")


# ===================== 6. Feature Alignment Across Datasets =====================

non_feature_cols <- c(
  "Sample",
  "Group",
  "label",
  "Dataset",
  "RowID"
)

train_features <- setdiff(colnames(train_dat), non_feature_cols)
test_features  <- setdiff(colnames(test_dat),  non_feature_cols)
val_features   <- setdiff(colnames(val_dat),   non_feature_cols)

common_features <- Reduce(intersect, list(train_features, test_features, val_features))

cat("[Feature Summary]\n")
cat("Feature count in train set: ", length(train_features), "\n")
cat("Feature count in test set: ", length(test_features), "\n")
cat("Feature count in validation set: ", length(val_features), "\n")
cat("Shared features across all sets: ", length(common_features), "\n\n")

if (length(common_features) < 2) {
  stop("Too few shared features detected, verify consistent gene names across train/test/validation files.")
}


# ===================== 7. Construct Numeric X & Y Matrices =====================

X_train_raw <- train_dat[, common_features, drop = FALSE]
X_test_raw  <- test_dat[,  common_features, drop = FALSE]
X_val_raw   <- val_dat[,   common_features, drop = FALSE]

Y_train <- train_dat$label
Y_test  <- test_dat$label
Y_val   <- val_dat$label

# Convert all expression values to numeric
to_numeric_df <- function(df) {
  df <- as.data.frame(df, check.names = FALSE)
  df[] <- lapply(df, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })
  return(df)
}

X_train_raw <- to_numeric_df(X_train_raw)
X_test_raw  <- to_numeric_df(X_test_raw)
X_val_raw   <- to_numeric_df(X_val_raw)

cat("[Missing value count after numeric conversion]\n")
cat("train NA total: ", sum(is.na(X_train_raw)), "\n")
cat("test NA total: ", sum(is.na(X_test_raw)), "\n")
cat("validation NA total: ", sum(is.na(X_val_raw)), "\n\n")


# ===================== 8. Feature Filtering: All-NA, High Missing Rate, Near-Zero Variance =====================

# All filtering thresholds derived solely from training set to avoid data leakage

# 1. Drop features with all NA values in training set
all_na_features <- sapply(X_train_raw, function(x) all(is.na(x)))

if (sum(all_na_features) > 0) {
  cat("Removed features with full NA values in training set: ", sum(all_na_features), "\n")
  X_train_raw <- X_train_raw[, !all_na_features, drop = FALSE]
  X_test_raw  <- X_test_raw[,  colnames(X_train_raw), drop = FALSE]
  X_val_raw   <- X_val_raw[,   colnames(X_train_raw), drop = FALSE]
}

# 2. Drop features with missing rate > 30% in training set
na_rate <- sapply(X_train_raw, function(x) mean(is.na(x)))
high_na_features <- na_rate > 0.3

if (sum(high_na_features) > 0) {
  cat("Removed features with missing rate >30% in training set: ", sum(high_na_features), "\n")
  X_train_raw <- X_train_raw[, !high_na_features, drop = FALSE]
  X_test_raw  <- X_test_raw[,  colnames(X_train_raw), drop = FALSE]
  X_val_raw   <- X_val_raw[,   colnames(X_train_raw), drop = FALSE]
}

# 3. Median imputation for remaining missing values (median calculated from training set only)
train_medians <- apply(X_train_raw, 2, median, na.rm = TRUE)

for (j in seq_along(X_train_raw)) {
  feature_name <- colnames(X_train_raw)[j]
  
  if (any(is.na(X_train_raw[[j]]))) {
    X_train_raw[[j]][is.na(X_train_raw[[j]])] <- train_medians[feature_name]
  }
  
  if (any(is.na(X_test_raw[[j]]))) {
    X_test_raw[[j]][is.na(X_test_raw[[j]])] <- train_medians[feature_name]
  }
  
  if (any(is.na(X_val_raw[[j]]))) {
    X_val_raw[[j]][is.na(X_val_raw[[j]])] <- train_medians[feature_name]
  }
}

cat("[Missing value count after median imputation]\n")
cat("train NA total: ", sum(is.na(X_train_raw)), "\n")
cat("test NA total: ", sum(is.na(X_test_raw)), "\n")
cat("validation NA total: ", sum(is.na(X_val_raw)), "\n\n")

# 4. Drop near-zero variance features
nzv_index <- caret::nearZeroVar(X_train_raw)

if (length(nzv_index) > 0) {
  cat("Removed near-zero variance features from training set: ", length(nzv_index), "\n")
  X_train_raw <- X_train_raw[, -nzv_index, drop = FALSE]
  X_test_raw  <- X_test_raw[,  colnames(X_train_raw), drop = FALSE]
  X_val_raw   <- X_val_raw[,   colnames(X_train_raw), drop = FALSE]
}

final_features <- colnames(X_train_raw)

cat("[Final feature count for model training]: ", length(final_features), "\n\n")


# ===================== 9. Feature Standardization =====================

# Standardization parameters (mean & SD) calculated only from training data
train_means <- apply(X_train_raw, 2, mean, na.rm = TRUE)
train_sds   <- apply(X_train_raw, 2, sd, na.rm = TRUE)

train_sds[is.na(train_sds) | train_sds == 0] <- 1

X_train <- as.data.frame(scale(X_train_raw, center = train_means, scale = train_sds))
X_test  <- as.data.frame(scale(X_test_raw,  center = train_means, scale = train_sds))
X_val   <- as.data.frame(scale(X_val_raw,   center = train_means, scale = train_sds))

# Force identical column order across all datasets
X_test <- X_test[, colnames(X_train), drop = FALSE]
X_val  <- X_val[,  colnames(X_train), drop = FALSE]

cat("[Final standardized matrix dimensions]\n")
cat("X_train: ", nrow(X_train), " × ", ncol(X_train), "\n")
cat("X_test: ", nrow(X_test), " × ", ncol(X_test), "\n")
cat("X_val: ", nrow(X_val), " × ", ncol(X_val), "\n\n")


# ===================== 10. Custom SuperLearner Base Learner Wrappers =====================

# Custom stable wrappers for multiple algorithms, including:
# 1. LASSO Logistic Regression
# 2. Ridge Logistic Regression
# 3. Elastic Net Logistic Regression
# 4. Ranger Random Forest
# 5. Shallow XGBoost
# 6. Deep XGBoost


SL.glmnet_lasso_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = Y,
    family = "binomial",
    alpha = 1,
    weights = obsWeights,
    nfolds = 5,
    type.measure = "auc"
  )
  
  pred <- as.numeric(
    predict(fit, newx = as.matrix(newX), s = "lambda.min", type = "response")
  )
  
  fit <- list(object = fit)
  class(fit) <- "SL.glmnet_lasso_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.glmnet_lasso_custom <- function(object, newdata, ...) {
  pred <- as.numeric(
    predict(object$object, newx = as.matrix(newdata), s = "lambda.min", type = "response")
  )
  return(pred)
}


SL.glmnet_ridge_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = Y,
    family = "binomial",
    alpha = 0,
    weights = obsWeights,
    nfolds = 5,
    type.measure = "auc"
  )
  
  pred <- as.numeric(
    predict(fit, newx = as.matrix(newX), s = "lambda.min", type = "response")
  )
  
  fit <- list(object = fit)
  class(fit) <- "SL.glmnet_ridge_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.glmnet_ridge_custom <- function(object, newdata, ...) {
  pred <- as.numeric(
    predict(object$object, newx = as.matrix(newdata), s = "lambda.min", type = "response")
  )
  return(pred)
}


SL.glmnet_enet_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  fit <- glmnet::cv.glmnet(
    x = as.matrix(X),
    y = Y,
    family = "binomial",
    alpha = 0.5,
    weights = obsWeights,
    nfolds = 5,
    type.measure = "auc"
  )
  
  pred <- as.numeric(
    predict(fit, newx = as.matrix(newX), s = "lambda.min", type = "response")
  )
  
  fit <- list(object = fit)
  class(fit) <- "SL.glmnet_enet_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.glmnet_enet_custom <- function(object, newdata, ...) {
  pred <- as.numeric(
    predict(object$object, newx = as.matrix(newdata), s = "lambda.min", type = "response")
  )
  return(pred)
}


SL.ranger_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dat <- as.data.frame(X)
  dat$Y <- factor(Y, levels = c(0, 1))
  
  fit <- ranger::ranger(
    Y ~ .,
    data = dat,
    probability = TRUE,
    num.trees = 500,
    mtry = max(1, floor(sqrt(ncol(X)))),
    min.node.size = 5,
    importance = "impurity",
    case.weights = obsWeights,
    seed = 123
  )
  
  pred <- predict(fit, data = as.data.frame(newX))$predictions[, "1"]
  
  fit <- list(object = fit)
  class(fit) <- "SL.ranger_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.ranger_custom <- function(object, newdata, ...) {
  pred <- predict(object$object, data = as.data.frame(newdata))$predictions[, "1"]
  return(pred)
}


SL.xgb_shallow_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dtrain <- xgboost::xgb.DMatrix(
    data = as.matrix(X),
    label = Y,
    weight = obsWeights
  )
  
  fit <- xgboost::xgb.train(
    data = dtrain,
    objective = "binary:logistic",
    eval_metric = "auc",
    nrounds = 80,
    max_depth = 2,
    eta = 0.05,
    subsample = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 3,
    lambda = 1,
    alpha = 0.5,
    verbose = 0
  )
  
  pred <- predict(fit, newdata = as.matrix(newX))
  
  fit <- list(object = fit)
  class(fit) <- "SL.xgb_shallow_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.xgb_shallow_custom <- function(object, newdata, ...) {
  pred <- predict(object$object, newdata = as.matrix(newdata))
  return(pred)
}


SL.xgb_deep_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dtrain <- xgboost::xgb.DMatrix(
    data = as.matrix(X),
    label = Y,
    weight = obsWeights
  )
  
  fit <- xgboost::xgb.train(
    data = dtrain,
    objective = "binary:logistic",
    eval_metric = "auc",
    nrounds = 100,
    max_depth = 3,
    eta = 0.03,
    subsample = 0.8,
    colsample_bytree = 0.7,
    min_child_weight = 5,
    lambda = 2,
    alpha = 1,
    verbose = 0
  )
  
  pred <- predict(fit, newdata = as.matrix(newX))
  
  fit <- list(object = fit)
  class(fit) <- "SL.xgb_deep_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.xgb_deep_custom <- function(object, newdata, ...) {
  pred <- predict(object$object, newdata = as.matrix(newdata))
  return(pred)
}


SL.logistic_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dat <- as.data.frame(X)
  dat$Y <- Y
  
  fit_obj <- tryCatch({
    glm(
      Y ~ .,
      data = dat,
      family = binomial(),
      weights = obsWeights,
      maxit = 100
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred <- tryCatch({
      as.numeric(predict(fit_obj, newdata = as.data.frame(newX), type = "response"))
    }, error = function(e) {
      rep(mean(Y), nrow(newX))
    })
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.logistic_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.logistic_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred <- tryCatch({
    as.numeric(predict(object$object, newdata = as.data.frame(newdata), type = "response"))
  }, error = function(e) {
    rep(object$fallback, nrow(newdata))
  })
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}
# ===================== 12. LASSO / Ridge / Elastic Net Wrappers =====================

# ---------- LASSO Regression (alpha = 1) ----------
SL.lasso_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  fit_obj <- tryCatch({
    glmnet::cv.glmnet(
      x = as.matrix(X),
      y = Y,
      family = "binomial",
      alpha = 1,
      weights = obsWeights,
      nfolds = 5,
      type.measure = "auc"
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred <- as.numeric(
      predict(fit_obj, newx = as.matrix(newX), s = "lambda.min", type = "response")
    )
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.lasso_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.lasso_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred <- as.numeric(
    predict(object$object, newx = as.matrix(newdata), s = "lambda.min", type = "response")
  )
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}


# ---------- Ridge Regression (alpha = 0) ----------
SL.ridge_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  fit_obj <- tryCatch({
    glmnet::cv.glmnet(
      x = as.matrix(X),
      y = Y,
      family = "binomial",
      alpha = 0,
      weights = obsWeights,
      nfolds = 5,
      type.measure = "auc"
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred <- as.numeric(
      predict(fit_obj, newx = as.matrix(newX), s = "lambda.min", type = "response")
    )
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.ridge_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.ridge_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred <- as.numeric(
    predict(object$object, newx = as.matrix(newdata), s = "lambda.min", type = "response")
  )
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}


# ---------- Elastic Net (alpha = 0.5) ----------
SL.elasticnet_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  fit_obj <- tryCatch({
    glmnet::cv.glmnet(
      x = as.matrix(X),
      y = Y,
      family = "binomial",
      alpha = 0.5,
      weights = obsWeights,
      nfolds = 5,
      type.measure = "auc"
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred <- as.numeric(
      predict(fit_obj, newx = as.matrix(newX), s = "lambda.min", type = "response")
    )
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.elasticnet_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.elasticnet_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred <- as.numeric(
    predict(object$object, newx = as.matrix(newdata), s = "lambda.min", type = "response")
  )
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}

# ===================== 13. Random Forest Wrapper =====================

SL.rf_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dat <- as.data.frame(X)
  dat$Y <- factor(Y, levels = c(0, 1))
  
  fit_obj <- tryCatch({
    ranger::ranger(
      Y ~ .,
      data = dat,
      probability = TRUE,
      num.trees = 500,
      mtry = max(1, floor(sqrt(ncol(X)))),
      min.node.size = 5,
      importance = "impurity",
      case.weights = obsWeights,
      seed = 123
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred_mat <- predict(fit_obj, data = as.data.frame(newX))$predictions
    pred <- pred_mat[, "1"]
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.rf_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.rf_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred_mat <- predict(object$object, data = as.data.frame(newdata))$predictions
  pred <- pred_mat[, "1"]
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}
# ===================== 14. SVM Wrapper =====================

SL.svm_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dat <- as.data.frame(X)
  dat$Y <- factor(Y, levels = c(0, 1))
  
  fit_obj <- tryCatch({
    e1071::svm(
      Y ~ .,
      data = dat,
      kernel = "radial",
      probability = TRUE,
      cost = 1,
      gamma = 1 / ncol(X),
      scale = FALSE
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred_class <- predict(
      fit_obj,
      newdata = as.data.frame(newX),
      probability = TRUE
    )
    
    prob_mat <- attr(pred_class, "probabilities")
    
    if ("1" %in% colnames(prob_mat)) {
      pred <- prob_mat[, "1"]
    } else {
      pred <- prob_mat[, ncol(prob_mat)]
    }
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.svm_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.svm_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred_class <- predict(
    object$object,
    newdata = as.data.frame(newdata),
    probability = TRUE
  )
  
  prob_mat <- attr(pred_class, "probabilities")
  
  if ("1" %in% colnames(prob_mat)) {
    pred <- prob_mat[, "1"]
  } else {
    pred <- prob_mat[, ncol(prob_mat)]
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}
# ===================== 15. XGBoost Wrapper =====================

SL.xgboost_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dtrain <- xgboost::xgb.DMatrix(
    data = as.matrix(X),
    label = Y,
    weight = obsWeights
  )
  
  fit_obj <- tryCatch({
    xgboost::xgb.train(
      data = dtrain,
      objective = "binary:logistic",
      eval_metric = "auc",
      nrounds = 100,
      max_depth = 3,
      eta = 0.03,
      subsample = 0.8,
      colsample_bytree = 0.8,
      min_child_weight = 3,
      lambda = 2,
      alpha = 0.5,
      verbose = 0
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred <- predict(fit_obj, newdata = as.matrix(newX))
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.xgboost_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.xgboost_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred <- predict(object$object, newdata = as.matrix(newdata))
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}
# ===================== 16. Gradient Boosting Machine Wrapper =====================

SL.gbm_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dat <- as.data.frame(X)
  dat$Y <- Y
  
  fit_obj <- tryCatch({
    gbm::gbm(
      formula = Y ~ .,
      data = dat,
      distribution = "bernoulli",
      weights = obsWeights,
      n.trees = 300,
      interaction.depth = 3,
      shrinkage = 0.03,
      n.minobsinnode = 5,
      bag.fraction = 0.8,
      train.fraction = 1,
      verbose = FALSE
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred <- as.numeric(
      predict(
        fit_obj,
        newdata = as.data.frame(newX),
        n.trees = 300,
        type = "response"
      )
    )
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.gbm_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.gbm_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred <- as.numeric(
    predict(
      object$object,
      newdata = as.data.frame(newdata),
      n.trees = 300,
      type = "response"
    )
  )
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}
# ===================== 17. Neural Network Wrapper =====================

SL.nnet_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dat <- as.data.frame(X)
  dat$Y <- Y
  
  fit_obj <- tryCatch({
    nnet::nnet(
      Y ~ .,
      data = dat,
      size = 5,
      decay = 0.01,
      maxit = 300,
      trace = FALSE,
      MaxNWts = 100000
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    pred <- as.numeric(
      predict(
        fit_obj,
        newdata = as.data.frame(newX),
        type = "raw"
      )
    )
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.nnet_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.nnet_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  pred <- as.numeric(
    predict(
      object$object,
      newdata = as.data.frame(newdata),
      type = "raw"
    )
  )
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}
# ===================== 18. Naive Bayes Wrapper =====================

SL.naivebayes_custom <- function(Y, X, newX, family, obsWeights, id, ...) {
  
  dat <- as.data.frame(X)
  dat$Y <- factor(Y, levels = c(0, 1))
  
  fit_obj <- tryCatch({
    e1071::naiveBayes(
      Y ~ .,
      data = dat,
      laplace = 1
    )
  }, error = function(e) {
    NULL
  })
  
  if (is.null(fit_obj)) {
    pred <- rep(mean(Y), nrow(newX))
  } else {
    prob_mat <- predict(
      fit_obj,
      newdata = as.data.frame(newX),
      type = "raw"
    )
    
    if ("1" %in% colnames(prob_mat)) {
      pred <- prob_mat[, "1"]
    } else {
      pred <- prob_mat[, ncol(prob_mat)]
    }
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  fit <- list(object = fit_obj, fallback = mean(Y))
  class(fit) <- "SL.naivebayes_custom"
  
  return(list(pred = pred, fit = fit))
}

predict.SL.naivebayes_custom <- function(object, newdata, ...) {
  
  if (is.null(object$object)) {
    return(rep(object$fallback, nrow(newdata)))
  }
  
  prob_mat <- predict(
    object$object,
    newdata = as.data.frame(newdata),
    type = "raw"
  )
  
  if ("1" %in% colnames(prob_mat)) {
    pred <- prob_mat[, "1"]
  } else {
    pred <- prob_mat[, ncol(prob_mat)]
  }
  
  pred[pred < 0] <- 0
  pred[pred > 1] <- 1
  
  return(pred)
}

# ===================== 19. Initialize SuperLearner Base Learner Library =====================

SL.library <- c(
  "SL.rf_custom",
  "SL.svm_custom",
  "SL.logistic_custom",
  "SL.lasso_custom",
  "SL.ridge_custom",
  "SL.elasticnet_custom",
  "SL.xgboost_custom",
  "SL.gbm_custom",
  "SL.nnet_custom",
  "SL.naivebayes_custom"
)

cat("[Base learner library loaded for SuperLearner]\n")
print(SL.library)
cat("\n")

}

library(SuperLearner)
library(pROC)
library(caret)

#===================== 20. Screen Optimal Meta-Learner =====================
# Iterate mainstream meta-learning algorithms, evaluate full classification metrics (AUC/ACC/SEN/SPE/PRE/F1/BACC + optimal threshold) across train/validation/test sets to select the best performer
library(SuperLearner)
library(pROC)
library(caret)
set.seed(123)

# Candidate meta-learning algorithms
meta_methods <- c(
  "method.NNLS",
  "method.NNLS2",
  "method.NNloglik",
  "method.CC_LS",
  "method.CC_nloglik",
  "method.AUC"
)
result_df <- data.frame()

# Metric calculation function (defined outside loop for efficiency)
calc_metric <- function(y_true, prob){
  roc_obj <- roc(y_true, prob)
  youden <- roc_obj$sensitivities + roc_obj$specificities - 1
  best_idx <- which.max(youden)
  best_thres <- roc_obj$thresholds[best_idx]
  y_pred <- as.integer(prob >= best_thres)
  
  cm <- confusionMatrix(factor(y_pred), factor(y_true), positive = "1")
  auc_val <- as.numeric(auc(roc_obj))
  acc <- cm$overall["Accuracy"]
  sen <- cm$byClass["Sensitivity"]
  spe <- cm$byClass["Specificity"]
  pre <- cm$byClass["Precision"]
  f1 <- cm$byClass["F1"]
  bacc <- cm$byClass["Balanced Accuracy"]
  
  return(c(
    AUC=auc_val,
    Accuracy=acc,
    Sensitivity=sen,
    Specificity=spe,
    Precision=pre,
    F1=f1,
    Balanced_Accuracy=bacc,
    Threshold=best_thres
  ))
}

total_meta <- length(meta_methods)
for(i in seq_along(meta_methods)){
  met <- meta_methods[i]
  progress_pct <- round(i / total_meta * 100, 1)
  cat(sprintf("\n[Progress: %s/%s | %.1f%%] Training meta-learner: %s\n", i, total_meta, progress_pct, met))
  
  # Try-catch error handling to skip failed algorithms
  one_res <- tryCatch({
    # Train SuperLearner ensemble
    sl_fit <- SuperLearner(
      Y = Y_train,
      X = X_train,
      family = binomial(),
      SL.library = SL.library,
      method = met,
      cvControl = list(V = 5),
      verbose = FALSE
    )
    
    # Generate predictions on all three datasets
    pred_tr <- predict(sl_fit, newdata = X_train)$pred[,1]
    pred_vl <- predict(sl_fit, newdata = X_val)$pred[,1]
    pred_ts <- predict(sl_fit, newdata = X_test)$pred[,1]
    
    # Calculate performance metrics
    met_tr <- calc_metric(Y_train, pred_tr)
    met_vl <- calc_metric(Y_val, pred_vl)
    met_ts <- calc_metric(Y_test, pred_ts)
    
    # Assemble result row
    row_tmp <- data.frame(
      Meta_Learner = met,
      Train_AUC=met_tr["AUC"],Train_Acc=met_tr["Accuracy"],Train_Sen=met_tr["Sensitivity"],
      Train_Spe=met_tr["Specificity"],Train_Pre=met_tr["Precision"],Train_F1=met_tr["F1"],
      Train_BACC=met_tr["Balanced_Accuracy"],Train_Thres=met_tr["Threshold"],
      
      Val_AUC=met_vl["AUC"],Val_Acc=met_vl["Accuracy"],Val_Sen=met_vl["Sensitivity"],
      Val_Spe=met_vl["Specificity"],Val_Pre=met_vl["Precision"],Val_F1=met_vl["F1"],
      Val_BACC=met_vl["Balanced_Accuracy"],Val_Thres=met_vl["Threshold"],
      
      Test_AUC=met_ts["AUC"],Test_Acc=met_ts["Accuracy"],Test_Sen=met_ts["Sensitivity"],
      Test_Spe=met_ts["Specificity"],Test_Pre=met_ts["Precision"],Test_F1=met_ts["F1"],
      Test_BACC=met_ts["Balanced_Accuracy"],Test_Thres=met_ts["Threshold"]
    )
    row_tmp
  }, error = function(e){
    # Print error message and skip failed meta-learner
    cat(sprintf("❌ %s failed, error message: %s\n", met, e$message))
    return(NULL)
  })
  
  # Append valid results only
  if(!is.null(one_res)){
    result_df <- rbind(result_df, one_res)
    cat(sprintf("✅ %s training complete\n",met))
  }
}

# Output full summary after iteration completes
cat("\n===== Iteration finished, valid meta-learner performance summary =====\n")
# print(round(result_df,4))

# Save results to CSV for backup
write.csv(result_df,"./MetaLearner_Result.csv",row.names = F)
cat("Meta-learner screening results saved to MetaLearner_Result.csv\n")

# Select optimal meta-learner if valid results exist
if(nrow(result_df) > 0){
  best_meta <- result_df$Meta_Learner[which.max(result_df$Test_AUC)]
  cat("\nOptimal meta-learner selected by test set AUC: ",best_meta,"\n")
  # Retrain full SuperLearner with optimal meta-algorithm
  final_fit <- SuperLearner(
    Y = Y_train,
    X = X_train,
    family = binomial(),
    SL.library = SL.library,
    method = best_meta,
    cvControl = list(V = 5),
    verbose = FALSE
  )
} else {
  cat("No valid meta-learner training results generated!\n")
}
# Train_AUC	Train set ROC-AUC: Near-perfect training performance often indicates overfitting
# Train_Acc	Overall Accuracy: Total proportion of correctly classified samples
# Train_Sen	Sensitivity (True Positive Rate): Proportion of all disease samples correctly identified
# Train_Spe	Specificity (True Negative Rate): Proportion of all healthy samples correctly identified
# Train_Pre	Precision: Ratio of true positives among all predicted disease samples
# Train_F1	F1-score: Harmonic mean of precision and sensitivity
# Train_BACC	Balanced Accuracy = (Sensitivity + Specificity) / 2: Balanced metric for imbalanced datasets
# Train_Thres	Optimal classification threshold via Youden Index: Samples with predicted probability ≥ threshold labeled as disease, threshold varies per meta-learner

# Meta_Learner Train_AUC Train_Acc Train_Sen Train_Spe Train_Pre Train_F1
# AUC        method.NNLS         1         1         1         1         1        1
# AUC1      method.NNLS2         1         1         1         1         1        1
# AUC2   method.NNloglik         1         1         1         1         1        1
# AUC3      method.CC_LS         1         1         1         1         1        1
# AUC4 method.CC_nloglik         1         1         1         1         1        1
# AUC5        method.AUC         1         1         1         1         1        1
# Train_BACC Train_Thres Val_AUC Val_Acc Val_Sen Val_Spe Val_Pre Val_F1 Val_BACC
# AUC           1   0.4685196       1       1       1       1       1      1        1
# AUC1          1   0.4693690       1       1       1       1       1      1        1
# AUC2          1   0.5206494       1       1       1       1       1      1        1
# AUC3          1   0.4664958       1       1       1       1       1      1        1
# AUC4          1   0.4838597       1       1       1       1       1      1        1
# AUC5          1   0.5254218       1       1       1       1       1      1        1
# Val_Thres  Test_AUC  Test_Acc  Test_Sen  Test_Spe  Test_Pre   Test_F1 Test_BACC
# AUC  0.7633073 0.9919095 0.9629630 0.9473684 0.9752066 0.9677419 0.9574468 0.9612875
# AUC1 0.7537984 0.9916485 0.9583333 0.9473684 0.9669421 0.9574468 0.9523810 0.9571553
# AUC2 0.6720362 0.9941714 0.9768519 0.9789474 0.9752066 0.9687500 0.9738220 0.9770770
# AUC3 0.7614412 0.9926925 0.9583333 0.9789474 0.9421488 0.9300000 0.9538462 0.9605481
# AUC4 0.7161939 0.9915615 0.9675926 0.9789474 0.9586777 0.9489796 0.9637306 0.9688125
# AUC5 0.6694324 0.9958243 0.9768519 0.9578947 0
