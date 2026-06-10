# ===================== Super Learner 诊断模型 =====================

# ===================== 0. 初始化环境 =====================

rm(list = ls())
gc()

{
# 设置工作路径（请根据实际情况修改）
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")
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

cat("【输出目录】\n")
cat("数据输出目录：", outdata_dir, "\n")
cat("图片输出目录：", outplot_dir, "\n\n")


# ===================== 1. 安装并加载必要 R 包 =====================

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


cat("【R 包加载完成】\n\n")


# ===================== 2. 参数设置 =====================

# ------------------------------------------------------------#
# 诊断模型中我们希望：
# 1 = 患病/AML/Tumor
# 0 = 正常/Normal
#
# 所以默认 DISEASE_CODE = 0
# 如果你的数据中 1 才是患病，请改成 DISEASE_CODE = 1
# ------------------------------------------------------------

DISEASE_CODE <- 1

set.seed(123)

train_file <- "../Outdata/5.all_data_harmony/5.1_all_expr_train.csv"
test_file  <- "../Outdata/5.all_data_harmony/5.2_all_expr_test.csv"
val_file   <- "../Outdata/5.all_data_harmony/5.3_all_expr_val.csv"

sample_info_file <- "../Outdata/5.all_data_harmony/5.0_all_expr_group_batch_info.csv"

cat("【标签设置】\n")
cat("原始 Group 中，", DISEASE_CODE, " 会被定义为患病阳性，即模型中的 label = 1\n\n")


# ===================== 3. 读取数据 =====================

cat("【读取表达矩阵】\n")

train_expr <- data.table::fread(train_file, data.table = FALSE, check.names = FALSE)
test_expr  <- data.table::fread(test_file,  data.table = FALSE, check.names = FALSE)
val_expr   <- data.table::fread(val_file,   data.table = FALSE, check.names = FALSE)

sample_info <- read.csv(
  sample_info_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("train_expr 维度：", nrow(train_expr), " × ", ncol(train_expr), "\n")
cat("test_expr 维度：", nrow(test_expr), " × ", ncol(test_expr), "\n")
cat("val_expr 维度：", nrow(val_expr), " × ", ncol(val_expr), "\n\n")

cat("train_expr 前几列：\n")
print(head(colnames(train_expr), 10))
cat("\n")


# ===================== 4. 数据整理函数 =====================

# 清理特征名
clean_feature_names <- function(x) {
  x <- gsub("[[:space:]]+", "", x)
  x <- gsub("-", ".", x)
  x <- make.names(x, unique = TRUE)
  return(x)
}


# 整理 fread 后的数据
# 你的 fread 提示：
# Detected 1050 column names but data has 1051 columns
# Added V1
# 这通常说明第一列是行名或样本 ID。
process_expr_table <- function(dat, dataset_name = "train") {
  
  dat <- as.data.frame(dat, check.names = FALSE)
  
  # 如果有 V1，一般是 fread 自动添加的行名列
  if ("V1" %in% colnames(dat)) {
    
    # 如果 Sample 列也存在，则 V1 通常只是重复行名，可保留为 RowID 或删除
    # 这里为了安全，若 V1 和 Sample 高度相似，则删除 V1
    if ("Sample" %in% colnames(dat)) {
      dat$V1 <- as.character(dat$V1)
      dat$Sample <- as.character(dat$Sample)
      
      same_rate <- mean(dat$V1 == dat$Sample, na.rm = TRUE)
      
      if (!is.na(same_rate) && same_rate > 0.8) {
        dat$V1 <- NULL
        cat("【", dataset_name, "】检测到 V1 与 Sample 基本一致，已删除 V1。\n", sep = "")
      } else {
        colnames(dat)[colnames(dat) == "V1"] <- "RowID"
        cat("【", dataset_name, "】检测到 V1，但与 Sample 不完全一致，重命名为 RowID。\n", sep = "")
      }
      
    } else {
      colnames(dat)[colnames(dat) == "V1"] <- "Sample"
      cat("【", dataset_name, "】检测到 V1 且无 Sample 列，已将 V1 重命名为 Sample。\n", sep = "")
    }
  }
  
  # 检查 Sample 列
  if (!"Sample" %in% colnames(dat)) {
    stop(paste0("【", dataset_name, "】数据中没有 Sample 列，请检查输入文件。"))
  }
  
  # 检查 Group 列
  if (!"Group" %in% colnames(dat)) {
    stop(paste0("【", dataset_name, "】数据中没有 Group 列，请检查输入文件。"))
  }
  
  # Sample 转字符
  dat$Sample <- as.character(dat$Sample)
  
  # Group 转字符，避免因子问题
  dat$Group <- as.character(dat$Group)
  
  # 清理除了 Sample 和 Group 之外的特征名
  feature_cols <- setdiff(colnames(dat), c("Sample", "Group", "RowID"))
  cleaned_feature_cols <- clean_feature_names(feature_cols)
  
  colnames(dat)[match(feature_cols, colnames(dat))] <- cleaned_feature_cols
  
  # 添加数据集名称
  dat$Dataset <- dataset_name
  
  return(dat)
}


train_dat <- process_expr_table(train_expr, "train")
test_dat  <- process_expr_table(test_expr,  "test")
val_dat   <- process_expr_table(val_expr,   "validation")

cat("\n【整理后维度】\n")
cat("train_dat：", nrow(train_dat), " × ", ncol(train_dat), "\n")
cat("test_dat：", nrow(test_dat), " × ", ncol(test_dat), "\n")
cat("val_dat：", nrow(val_dat), " × ", ncol(val_dat), "\n\n")


# ===================== 5. 标签重编码 =====================

# 输出 label:
# 1 = 患病/AML/Tumor
# 0 = 正常/Normal
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

cat("【原始 Group 分布】\n")
cat("train:\n")
print(table(train_dat$Group, useNA = "ifany"))
cat("test:\n")
print(table(test_dat$Group, useNA = "ifany"))
cat("validation:\n")
print(table(val_dat$Group, useNA = "ifany"))

cat("\n【重编码后 label 分布：0=正常，1=患病】\n")
cat("train:\n")
print(table(train_dat$label, useNA = "ifany"))
cat("test:\n")
print(table(test_dat$label, useNA = "ifany"))
cat("validation:\n")
print(table(val_dat$label, useNA = "ifany"))
cat("\n")


# ===================== 6. 特征对齐 =====================

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

cat("【特征信息】\n")
cat("train 特征数：", length(train_features), "\n")
cat("test 特征数：", length(test_features), "\n")
cat("validation 特征数：", length(val_features), "\n")
cat("共同特征数：", length(common_features), "\n\n")

if (length(common_features) < 2) {
  stop("共同特征数过少，请检查 train/test/validation 的列名是否一致。")
}


# ===================== 7. 构建 X/Y 并数值化 =====================

X_train_raw <- train_dat[, common_features, drop = FALSE]
X_test_raw  <- test_dat[,  common_features, drop = FALSE]
X_val_raw   <- val_dat[,   common_features, drop = FALSE]

Y_train <- train_dat$label
Y_test  <- test_dat$label
Y_val   <- val_dat$label

# 转为数值型
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

cat("【数值化后缺失值数量】\n")
cat("train NA：", sum(is.na(X_train_raw)), "\n")
cat("test NA：", sum(is.na(X_test_raw)), "\n")
cat("validation NA：", sum(is.na(X_val_raw)), "\n\n")


# ===================== 8. 特征过滤：全 NA、缺失率高、近零方差 =====================

# 只根据训练集决定过滤规则，避免数据泄漏

# 1. 去除训练集中全 NA 的特征
all_na_features <- sapply(X_train_raw, function(x) all(is.na(x)))

if (sum(all_na_features) > 0) {
  cat("去除训练集中全 NA 特征数：", sum(all_na_features), "\n")
  X_train_raw <- X_train_raw[, !all_na_features, drop = FALSE]
  X_test_raw  <- X_test_raw[,  colnames(X_train_raw), drop = FALSE]
  X_val_raw   <- X_val_raw[,   colnames(X_train_raw), drop = FALSE]
}

# 2. 去除训练集中缺失率过高的特征，比如 > 30%
na_rate <- sapply(X_train_raw, function(x) mean(is.na(x)))
high_na_features <- na_rate > 0.3

if (sum(high_na_features) > 0) {
  cat("去除训练集中缺失率 > 30% 的特征数：", sum(high_na_features), "\n")
  X_train_raw <- X_train_raw[, !high_na_features, drop = FALSE]
  X_test_raw  <- X_test_raw[,  colnames(X_train_raw), drop = FALSE]
  X_val_raw   <- X_val_raw[,   colnames(X_train_raw), drop = FALSE]
}

# 3. 中位数填补缺失值
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

cat("【填补后缺失值数量】\n")
cat("train NA：", sum(is.na(X_train_raw)), "\n")
cat("test NA：", sum(is.na(X_test_raw)), "\n")
cat("validation NA：", sum(is.na(X_val_raw)), "\n\n")

# 4. 去除近零方差特征
nzv_index <- caret::nearZeroVar(X_train_raw)

if (length(nzv_index) > 0) {
  cat("去除训练集中近零方差特征数：", length(nzv_index), "\n")
  X_train_raw <- X_train_raw[, -nzv_index, drop = FALSE]
  X_test_raw  <- X_test_raw[,  colnames(X_train_raw), drop = FALSE]
  X_val_raw   <- X_val_raw[,   colnames(X_train_raw), drop = FALSE]
}

final_features <- colnames(X_train_raw)

cat("【最终用于建模的特征数】：", length(final_features), "\n\n")


# ===================== 9. 标准化 =====================

# 只用训练集均值和标准差
train_means <- apply(X_train_raw, 2, mean, na.rm = TRUE)
train_sds   <- apply(X_train_raw, 2, sd, na.rm = TRUE)

train_sds[is.na(train_sds) | train_sds == 0] <- 1

X_train <- as.data.frame(scale(X_train_raw, center = train_means, scale = train_sds))
X_test  <- as.data.frame(scale(X_test_raw,  center = train_means, scale = train_sds))
X_val   <- as.data.frame(scale(X_val_raw,   center = train_means, scale = train_sds))

# 确保列顺序一致
X_test <- X_test[, colnames(X_train), drop = FALSE]
X_val  <- X_val[,  colnames(X_train), drop = FALSE]

cat("【最终矩阵维度】\n")
cat("X_train：", nrow(X_train), " × ", ncol(X_train), "\n")
cat("X_test：", nrow(X_test), " × ", ncol(X_test), "\n")
cat("X_val：", nrow(X_val), " × ", ncol(X_val), "\n\n")


# ===================== 10. 自定义 SuperLearner wrapper =====================

# 为了让模型更稳定，这里定义几个自定义 wrapper。
# 包括：
# 1. LASSO logistic
# 2. Ridge logistic
# 3. Elastic Net logistic
# 4. Ranger random forest
# 5. XGBoost 浅树
# 6. XGBoost 稍深树


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
# ===================== 12. Lasso / Ridge / Elastic Net wrapper =====================

# ---------- Lasso Regression, alpha = 1 ----------
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


# ---------- Ridge Regression, alpha = 0 ----------
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


# ---------- Elastic Net, alpha = 0.5 ----------
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

# ===================== 13. Random Forest wrapper =====================

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
# ===================== 14. SVM wrapper =====================

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
# ===================== 15. XGBoost wrapper =====================

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
# ===================== 16. Gradient Boosting Machine wrapper =====================

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
# ===================== 17. Neural Network wrapper =====================

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
# ===================== 18. Naive Bayes wrapper =====================

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

# ===================== 19. 设置 SuperLearner 模型库 =====================

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

cat("【SuperLearner 使用的模型库】\n")
print(SL.library)
cat("\n")

}

library(SuperLearner)
library(pROC)
library(caret)

#===================== 20. 选择最佳的元学习器=====================
# 遍历全部主流元学习器 + 三数据集 (训练 / 验证 / 测试) 全套分类指标 (AUC/ACC/SEN/SPE/PRE/F1/BACC/ 最优阈值) 择优
library(SuperLearner)
library(pROC)
library(caret)
set.seed(123)

# 候选元学习列表
meta_methods <- c(
  "method.NNLS",
  "method.NNLS2",
  "method.NNloglik",
  "method.CC_LS",
  "method.CC_nloglik",
  "method.AUC"
)
result_df <- data.frame()

# 指标函数放到循环外
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
  cat(sprintf("\n【进度：%s/%s | %.1f%%】正在训练元学习器：%s\n", i, total_meta, progress_pct, met))
  
  # tryCatch捕获异常，失败直接跳过
  one_res <- tryCatch({
    # 训练SuperLearner
    sl_fit <- SuperLearner(
      Y = Y_train,
      X = X_train,
      family = binomial(),
      SL.library = SL.library,
      method = met,
      cvControl = list(V = 5),
      verbose = FALSE
    )
    
    # 三集预测
    pred_tr <- predict(sl_fit, newdata = X_train)$pred[,1]
    pred_vl <- predict(sl_fit, newdata = X_val)$pred[,1]
    pred_ts <- predict(sl_fit, newdata = X_test)$pred[,1]
    
    # 计算指标
    met_tr <- calc_metric(Y_train, pred_tr)
    met_vl <- calc_metric(Y_val, pred_vl)
    met_ts <- calc_metric(Y_test, pred_ts)
    
    # 组装结果行
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
    # 报错回调：打印错误，返回NULL，不写入结果
    cat(sprintf("❌ %s 报错跳过，错误信息：%s\n", met, e$message))
    return(NULL)
  })
  
  # 非NULL才追加进结果表
  if(!is.null(one_res)){
    result_df <- rbind(result_df, one_res)
    cat(sprintf("✅ %s 训练完成\n",met))
  }
}

# 全部跑完输出汇总
cat("\n=====全部遍历结束，有效模型结果汇总=====\n")
# print(round(result_df,4))

# 保存结果到csv（防止丢数据）
write.csv(result_df,"./MetaLearner_Result.csv",row.names = F)
cat("结果已保存至 MetaLearner_Result.csv\n")

# 筛选最优（有结果才筛选）
if(nrow(result_df) > 0){
  best_meta <- result_df$Meta_Learner[which.max(result_df$Test_AUC)]
  cat("\n综合最优元学习器（按测试集AUC）：",best_meta,"\n")
  # 最优模型重训
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
  cat("无可用元学习结果！\n")
}
# Train_AUC	训练集 ROC-AUC	基模型在训练集近乎完美拟合，训练集重度过拟合
# Train_Acc	准确率 Accuracy	全部分类正确
# Train_Sen	灵敏度 Sensitivity (真阳性率)	患病样本全部检出
# Train_Spe	特异度 Specificity (真阴性率)	正常样本全部识别
# Train_Pre	精确率 Precision	预测患病里真实患病比例 = 1
# Train_F1	F1-score	精确 + 灵敏度调和平均 = 1
# Train_BACC	平衡准确率 Balanced Accuracy=(Sen+Spe)/2	=1 代表无错分
# Train_Thres	Youden 指数最优分类阈值	概率≥该值判为患病，每行阈值不同

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
# AUC5 0.6694324 0.9958243 0.9768519 0.9578947 0.9917355 0.9891304 0.9732620 0.9748151
# Test_Thres
# AUC   0.8269062
# AUC1  0.8506186
# AUC2  0.6228904
# AUC3  0.7415043
# AUC4  0.7586088
# AUC5  0.7574576

# ===================== 20. 根据上面选择最佳的元学习器，再来训练 SuperLearner 模型 =====================


best_meta <- "method.AUC"

cat("【开始训练 SuperLearner 模型】\n")
cat("训练样本数：", nrow(X_train), "\n")
cat("训练特征数：", ncol(X_train), "\n")
cat("训练标签分布：\n")
print(table(Y_train))
cat("\n")

set.seed(123)

sl_fit <- SuperLearner::SuperLearner(
  Y = Y_train,
  X = X_train,
  family = binomial(),
  SL.library = SL.library,
  # method = "method.AUC",
  method = best_meta,
  cvControl = list(V = 10),
  verbose = TRUE
)

cat("\n【SuperLearner 训练完成】\n")
print(sl_fit)
#                             Risk       Coef
# SL.rf_custom_All         0.0001896667 0.36508504
# SL.svm_custom_All        0.0197323592 0.00000000
# SL.logistic_custom_All   0.5251343472 0.06201792
# SL.lasso_custom_All      0.0170559517 0.00000000
# SL.ridge_custom_All      0.0517860279 0.00000000
# SL.elasticnet_custom_All 0.0142741737 0.00000000
# SL.xgboost_custom_All    0.0041656422 0.14959154
# SL.gbm_custom_All        0.0007656914 0.35243572
# SL.nnet_custom_All       0.0299041130 0.07086978
# SL.naivebayes_custom_All 0.1845105546 0.00000000

cat("\n【各基学习器权重】\n")
print(sl_fit$coef)
# SL.rf_custom_All        SL.svm_custom_All   SL.logistic_custom_All      SL.lasso_custom_All      SL.ridge_custom_All SL.elasticnet_custom_All 
# 0.36508504               0.00000000               0.06201792               0.00000000               0.00000000               0.00000000 
# SL.xgboost_custom_All        SL.gbm_custom_All       SL.nnet_custom_All SL.naivebayes_custom_All 
# 0.14959154               0.35243572               0.07086978               0.00000000

coef_df <- data.frame(
  Learner = names(sl_fit$coef),
  Weight = as.numeric(sl_fit$coef),
  stringsAsFactors = FALSE
)
weight_df <- data.frame(
  BaseLearner = names(sl_fit$coef),
  Risk        = sl_fit$cvRisk,
  Coefficient = sl_fit$coef
)

write.csv(
  weight_df,
  file = file.path(outdata_dir, "1_SuperLearner_base_learner_weights.csv"),
  row.names = FALSE
)





# ===================== 21. SperLearn模型 预测 train/test/validation 每个样本的概率 =====================

predict_sl_prob <- function(model, newX) {
  
  pred_obj <- predict(
    model,
    newdata = newX,
    onlySL = TRUE
  )
  
  prob <- as.numeric(pred_obj$pred)
  
  prob[prob < 0] <- 0
  prob[prob > 1] <- 1
  
  return(prob)
}

prob_train <- predict_sl_prob(sl_fit, X_train)
prob_test  <- predict_sl_prob(sl_fit, X_test)
prob_val   <- predict_sl_prob(sl_fit, X_val)

cat("【预测完成】\n")
cat("train prob 范围：", round(min(prob_train), 4), "-", round(max(prob_train), 4), "\n")
cat("test prob 范围：", round(min(prob_test), 4), "-", round(max(prob_test), 4), "\n")
cat("validation prob 范围：", round(min(prob_val), 4), "-", round(max(prob_val), 4), "\n\n")


# ===================== 22. 基于 train 集 ROC 确定最佳阈值 =====================

roc_train <- pROC::roc(
  response = Y_train,
  predictor = prob_train,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

best_threshold <- as.numeric(
  pROC::coords(
    roc_train,
    x = "best",
    best.method = "youden",
    ret = "threshold"
  )
)
best_threshold 
# [1] 0.4848475 method_LNN
# [1] 0.5147276 method_AUC


if (length(best_threshold) > 1) {
  best_threshold <- best_threshold[1]
}

if (is.na(best_threshold) || is.infinite(best_threshold)) {
  best_threshold <- 0.5
}

cat("【最佳诊断阈值】\n")
cat("基于 train 集 Youden index 得到 threshold = ", round(best_threshold, 4), "\n\n")

# ===================== 23. 定义评价指标函数 =====================

calculate_metrics <- function(true_labels, pred_probs, threshold = 0.5) {
  
  true_labels <- as.numeric(true_labels)
  pred_probs <- as.numeric(pred_probs)
  
  pred_class <- ifelse(pred_probs >= threshold, 1, 0)
  
  auc_value <- tryCatch({
    roc_obj <- pROC::roc(
      response = true_labels,
      predictor = pred_probs,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) {
    NA
  })
  
  cm <- table(
    Actual = factor(true_labels, levels = c(0, 1)),
    Predicted = factor(pred_class, levels = c(0, 1))
  )
  
  tn <- as.numeric(cm[1, 1])
  fp <- as.numeric(cm[1, 2])
  fn <- as.numeric(cm[2, 1])
  tp <- as.numeric(cm[2, 2])
  
  accuracy <- ifelse(sum(cm) > 0, (tp + tn) / sum(cm), NA)
  sensitivity <- ifelse((tp + fn) > 0, tp / (tp + fn), NA)
  specificity <- ifelse((tn + fp) > 0, tn / (tn + fp), NA)
  precision <- ifelse((tp + fp) > 0, tp / (tp + fp), NA)
  
  f1_score <- ifelse(
    !is.na(precision) & !is.na(sensitivity) & (precision + sensitivity) > 0,
    2 * precision * sensitivity / (precision + sensitivity),
    NA
  )
  
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  
  return(list(
    AUC = round(auc_value, 4),
    Accuracy = round(accuracy, 4),
    Sensitivity = round(sensitivity, 4),
    Specificity = round(specificity, 4),
    Precision = round(precision, 4),
    F1_Score = round(f1_score, 4),
    Balanced_Accuracy = round(balanced_accuracy, 4),
    Threshold = round(threshold, 4),
    Confusion_Matrix = cm,
    Pred_Class = pred_class
  ))
}

# ===================== 24. 评估 train/test/validation =====================

metrics_train <- calculate_metrics(
  true_labels = Y_train,
  pred_probs = prob_train,
  threshold = best_threshold
)

metrics_test <- calculate_metrics(
  true_labels = Y_test,
  pred_probs = prob_test,
  threshold = best_threshold
)

metrics_val <- calculate_metrics(
  true_labels = Y_val,
  pred_probs = prob_val,
  threshold = best_threshold
)

cat("【train 混淆矩阵】\n")
print(metrics_train$Confusion_Matrix)
# Predicted
# Actual   0   1
# 0 355   0
# 1   0 401

cat("\n【test 混淆矩阵】\n")
print(metrics_test$Confusion_Matrix)
# Predicted
# Actual   0   1
# 0 109  12
# 1   1  94

# Predicted
# Actual   0   1
# 0 355   0
# 1   0 401

cat("\n【validation 混淆矩阵】\n")
print(metrics_val$Confusion_Matrix)
# Predicted
# Actual  0  1
# 0 50  3
# 1  0 55

# Predicted
# Actual   0   1
# 0 111  10
# 1   2  93

sl_metrics_df <- data.frame(
  Model = "SuperLearner",
  Dataset = c("train", "test", "validation"),
  AUC = c(metrics_train$AUC, metrics_test$AUC, metrics_val$AUC),
  Accuracy = c(metrics_train$Accuracy, metrics_test$Accuracy, metrics_val$Accuracy),
  Sensitivity = c(metrics_train$Sensitivity, metrics_test$Sensitivity, metrics_val$Sensitivity),
  Specificity = c(metrics_train$Specificity, metrics_test$Specificity, metrics_val$Specificity),
  Precision = c(metrics_train$Precision, metrics_test$Precision, metrics_val$Precision),
  F1_Score = c(metrics_train$F1_Score, metrics_test$F1_Score, metrics_val$F1_Score),
  Balanced_Accuracy = c(
    metrics_train$Balanced_Accuracy,
    metrics_test$Balanced_Accuracy,
    metrics_val$Balanced_Accuracy
  ),
  Threshold = best_threshold,
  stringsAsFactors = FALSE
)

sl_mean_df <- data.frame(
  Model = "SuperLearner",
  Dataset = "Mean",
  AUC = round(mean(sl_metrics_df$AUC, na.rm = TRUE), 4),
  Accuracy = round(mean(sl_metrics_df$Accuracy, na.rm = TRUE), 4),
  Sensitivity = round(mean(sl_metrics_df$Sensitivity, na.rm = TRUE), 4),
  Specificity = round(mean(sl_metrics_df$Specificity, na.rm = TRUE), 4),
  Precision = round(mean(sl_metrics_df$Precision, na.rm = TRUE), 4),
  F1_Score = round(mean(sl_metrics_df$F1_Score, na.rm = TRUE), 4),
  Balanced_Accuracy = round(mean(sl_metrics_df$Balanced_Accuracy, na.rm = TRUE), 4),
  Threshold = best_threshold,
  stringsAsFactors = FALSE
)

sl_final_metrics <- rbind(sl_metrics_df, sl_mean_df)

cat("\n===================== SuperLearner 最终评估结果 =====================\n")
print(sl_final_metrics)

# method_auc
# Model    Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold
# 1 SuperLearner      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 2 SuperLearner       test 0.9921   0.9444      0.9789      0.9174    0.9029   0.9394            0.9482 0.5147276
# 3 SuperLearner validation 1.0000   0.9815      1.0000      0.9623    0.9649   0.9821            0.9811 0.5147276
# 4 SuperLearner       Mean 0.9974   0.9753      0.9930      0.9599    0.9559   0.9738            0.9764 0.5147276

write.csv(
  sl_final_metrics,
  file = file.path(outdata_dir, "2_SuperLearner_train_test_validation_metrics.csv"),
  row.names = FALSE
)

# ===================== 25. 保存 train/test/validation 预测结果 =====================

pred_train_df <- data.frame(
  Sample = train_dat$Sample,
  Dataset = "train",
  Original_Group = train_dat$Group,
  True_Label = Y_train,
  Pred_Prob = prob_train,
  Pred_Class = metrics_train$Pred_Class,
  stringsAsFactors = FALSE
)

pred_test_df <- data.frame(
  Sample = test_dat$Sample,
  Dataset = "test",
  Original_Group = test_dat$Group,
  True_Label = Y_test,
  Pred_Prob = prob_test,
  Pred_Class = metrics_test$Pred_Class,
  stringsAsFactors = FALSE
)

pred_val_df <- data.frame(
  Sample = val_dat$Sample,
  Dataset = "validation",
  Original_Group = val_dat$Group,
  True_Label = Y_val,
  Pred_Prob = prob_val,
  Pred_Class = metrics_val$Pred_Class,
  stringsAsFactors = FALSE
)

sl_pred_all <- rbind(
  pred_train_df,
  pred_test_df,
  pred_val_df
)
sl_pred_all

write.csv(
  sl_pred_all,
  file = file.path(outdata_dir, "3_SuperLearner_train_test_validation_predictions.csv"),
  row.names = FALSE
)


# ===================== 26. 绘制 ROC 曲线 =====================

roc_test <- pROC::roc(
  response = Y_test,
  predictor = prob_test,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_val <- pROC::roc(
  response = Y_val,
  predictor = prob_val,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

pdf(
  file.path(outplot_dir, "1_SuperLearner_ROC_train_test_validation.pdf"),
  width = 5,
  height = 5
)

plot(
  roc_train,
  col = "#264653",
  lwd = 2,
  main = "SuperLearner ROC Curves"
)

plot(
  roc_test,
  col = "#2A9D8F",
  lwd = 2,
  add = TRUE
)

plot(
  roc_val,
  col = "#E76F51",
  lwd = 2,
  add = TRUE
)

abline(a = 1, b = -1, lty = 2, col = "gray")

legend(
  "bottomright",
  legend = c(
    paste0("Train AUC = ", round(as.numeric(pROC::auc(roc_train)), 4)),
    paste0("Test AUC = ", round(as.numeric(pROC::auc(roc_test)), 4)),
    paste0("Validation AUC = ", round(as.numeric(pROC::auc(roc_val)), 4))
  ),
  col = c("#264653", "#2A9D8F", "#E76F51"),
  lwd = 2,
  bty = "n"
)

dev.off()
# ===================== 27. 绘制 validation 指标柱状图 =====================

val_metric_plot_df <- sl_final_metrics %>%
  dplyr::filter(Dataset == "validation") %>%
  dplyr::select(
    AUC,
    Accuracy,
    F1_Score,
    Sensitivity,
    Specificity,
    Precision,
    Balanced_Accuracy
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Metric",
    values_to = "Value"
  )

pdf(
  file.path(outplot_dir, "2_SuperLearner_validation_metrics_barplot.pdf"),
  width = 8,
  height = 5
)

ggplot(val_metric_plot_df, aes(x = Metric, y = Value)) +
  geom_bar(
    stat = "identity",
    fill = "#E76F51",
    color = "black",
    linewidth = 0.5,
    width = 0.65
  ) +
  geom_text(
    aes(label = round(Value, 4)),
    vjust = -0.7,
    size = 3.6,
    fontface = "bold"
  ) +
  labs(
    title = "SuperLearner Validation Metrics",
    x = "Evaluation metrics",
    y = "Value"
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, size = 14, fontface = "bold"),
    axis.title = element_text(size = 12, fontface = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      color = "black",
      fontface = "bold"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black",
      fontface = "bold"
    ),
    legend.position = "none"
  )

dev.off()

# ===================== 28. 绘制 SuperLearner 基学习器权重图 =====================

coef_plot_df <- coef_df %>%
  dplyr::arrange(desc(Weight))

pdf(
  file.path(outplot_dir, "3_SuperLearner_base_learner_weights.pdf"),
  width = 8,
  height = 5
)

ggplot(coef_plot_df, aes(x = reorder(Learner, Weight), y = Weight)) +
  geom_bar(
    stat = "identity",
    fill = "#2A9D8F",
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = round(Weight, 4)),
    hjust = -0.1,
    size = 3.5,
    fontface = "bold"
  ) +
  coord_flip() +
  labs(
    title = "SuperLearner Base Learner Weights",
    x = "Base learner",
    y = "Weight"
  ) +
  scale_y_continuous(
    limits = c(0, max(coef_plot_df$Weight, na.rm = TRUE) * 1.2 + 0.01),
    expand = c(0, 0)
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, size = 14, fontface = "bold"),
    axis.title = element_text(size = 12, fontface = "bold"),
    axis.text = element_text(size = 10, color = "black", fontface = "bold")
  )

dev.off()
# ===================== 29. 保存模型和关键对象 =====================

save(
  sl_fit,
  SL.library,
  coef_df,
  sl_final_metrics,
  sl_pred_all,
  final_features,
  train_medians,
  train_means,
  train_sds,
  best_threshold,
  DISEASE_CODE,
  file = file.path(outdata_dir, "4_SuperLearner_full_model_and_results.RData")
)

cat("\n===================== SuperLearner 建模全部完成 =====================\n")
cat("模型文件已保存：", file.path(outdata_dir, "4_SuperLearner_full_model_and_results.RData"), "\n")
cat("评估指标已保存：", file.path(outdata_dir, "2_SuperLearner_train_test_validation_metrics.csv"), "\n")
cat("预测结果已保存：", file.path(outdata_dir, "3_SuperLearner_train_test_validation_predictions.csv"), "\n")
cat("基学习器权重已保存：", file.path(outdata_dir, "1_SuperLearner_base_learner_weights.csv"), "\n")
cat("图片输出目录：", outplot_dir, "\n")

# ===================== 30. 查看最终结果 =====================

cat("\n【SuperLearner 最终评估结果】\n")
print(sl_final_metrics)
# Model    Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold
# 1 SuperLearner      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 2 SuperLearner       test 0.9921   0.9444      0.9789      0.9174    0.9029   0.9394            0.9482 0.5147276
# 3 SuperLearner validation 1.0000   0.9815      1.0000      0.9623    0.9649   0.9821            0.9811 0.5147276
# 4 SuperLearner       Mean 0.9974   0.9753      0.9930      0.9599    0.9559   0.9738            0.9764 0.5147276

cat("\n【SuperLearner 基学习器权重】\n")
print(coef_df)
 
# Learner     Weight  method_AUC
# 1          SL.rf_custom_All 0.36508504
# 2         SL.svm_custom_All 0.00000000
# 3    SL.logistic_custom_All 0.06201792
# 4       SL.lasso_custom_All 0.00000000
# 5       SL.ridge_custom_All 0.00000000
# 6  SL.elasticnet_custom_All 0.00000000
# 7     SL.xgboost_custom_All 0.14959154
# 8         SL.gbm_custom_All 0.35243572
# 9        SL.nnet_custom_All 0.07086978
# 10 SL.naivebayes_custom_All 0.00000000

cat("\n【Validation 集预测类别分布】\n")
print(table(pred_val_df$Pred_Class, useNA = "ifany"))

cat("\n【Validation 集真实类别分布】\n")
print(table(pred_val_df$True_Label, useNA = "ifany"))












# ===================== 31.查看 SuperLearner 各基学习器 CV 风险 =====================

cat("\n【SuperLearner 交叉验证风险 cvRisk】\n")
print(sl_fit$cvRisk)

cv_risk_df <- data.frame(
  Learner = sl_fit$libraryNames,
  CV_Risk = as.numeric(sl_fit$cvRisk),
  Weight  = as.numeric(sl_fit$coef),
  stringsAsFactors = FALSE
)

cv_risk_df <- cv_risk_df[order(cv_risk_df$CV_Risk), ]

cat("\n【按 CV_Risk 从小到大排序】\n")
print(cv_risk_df)
# Learner      CV_Risk     Weight
# 1          SL.rf_custom_All 0.0001896667 0.36508504
# 8         SL.gbm_custom_All 0.0007656914 0.35243572
# 7     SL.xgboost_custom_All 0.0041656422 0.14959154
# 6  SL.elasticnet_custom_All 0.0142741737 0.00000000
# 4       SL.lasso_custom_All 0.0170559517 0.00000000
# 2         SL.svm_custom_All 0.0197323592 0.00000000
# 9        SL.nnet_custom_All 0.0299041130 0.07086978
# 5       SL.ridge_custom_All 0.0517860279 0.00000000
# 10 SL.naivebayes_custom_All 0.1845105546 0.00000000
# 3    SL.logistic_custom_All 0.5251343472 0.06201792

write.csv(
  cv_risk_df,
  file = file.path(outdata_dir, "5_SuperLearner_cvRisk_and_weights.csv"),
  row.names = FALSE
)

# ===================== 32. 提取每个基学习器单独预测概率 =====================

get_base_learner_preds <- function(sl_model, newX) {
  
  pred_all <- predict(
    sl_model,
    newdata = newX,
    onlySL = FALSE
  )
  
  # SuperLearner 的 predict 结果一般包含 pred 和 library.predict
  base_pred <- as.data.frame(pred_all$library.predict)
  
  # 清理列名
  colnames(base_pred) <- names(sl_model$coef)
  
  return(base_pred)
}

base_pred_train <- get_base_learner_preds(sl_fit, X_train)
base_pred_test  <- get_base_learner_preds(sl_fit, X_test)
base_pred_val   <- get_base_learner_preds(sl_fit, X_val)

cat("【基学习器预测矩阵维度】\n")
cat("train:", dim(base_pred_train), "\n")
cat("test:", dim(base_pred_test), "\n")
cat("validation:", dim(base_pred_val), "\n\n")

cat("【基学习器预测矩阵前几列】\n")
print(head(base_pred_val[, 1:min(5, ncol(base_pred_val)), drop = FALSE]))
# SL.rf_custom_All SL.svm_custom_All SL.logistic_custom_All SL.lasso_custom_All SL.ridge_custom_All
# 1        0.9482333         0.9593562           2.220446e-16           0.9794983           0.7519618
# 2        0.9529000         0.9630534           2.220446e-16           0.9420656           0.8474782
# 3        0.9056000         0.9443203           1.000000e+00           0.9374133           0.7155813
# 4        0.9347000         0.9431039           1.000000e+00           0.9041275           0.7746457
# 5        0.9412000         0.9364597           1.000000e+00           0.9227628           0.7606836
# 6        0.8537000         0.9540757           2.220446e-16           0.9021391           0.8232932

# ===================== 33. 逐个基学习器评估 =====================

evaluate_base_learners <- function(base_pred_df, true_y, dataset_name, threshold = 0.5) {
  
  res_list <- list()
  
  for (learner in colnames(base_pred_df)) {
    
    prob <- as.numeric(base_pred_df[[learner]])
    
    metric <- calculate_metrics(
      true_labels = true_y,
      pred_probs = prob,
      threshold = threshold
    )
    
    res_list[[learner]] <- data.frame(
      Learner = learner,
      Dataset = dataset_name,
      AUC = metric$AUC,
      Accuracy = metric$Accuracy,
      Sensitivity = metric$Sensitivity,
      Specificity = metric$Specificity,
      Precision = metric$Precision,
      F1_Score = metric$F1_Score,
      Balanced_Accuracy = metric$Balanced_Accuracy,
      Threshold = threshold,
      stringsAsFactors = FALSE
    )
  }
  
  res_df <- do.call(rbind, res_list)
  rownames(res_df) <- NULL
  
  return(res_df)
}

base_metrics_train <- evaluate_base_learners(
  base_pred_df = base_pred_train,
  true_y = Y_train,
  dataset_name = "train",
  threshold = best_threshold
)

base_metrics_test <- evaluate_base_learners(
  base_pred_df = base_pred_test,
  true_y = Y_test,
  dataset_name = "test",
  threshold = best_threshold
)

base_metrics_val <- evaluate_base_learners(
  base_pred_df = base_pred_val,
  true_y = Y_val,
  dataset_name = "validation",
  threshold = best_threshold
)

base_metrics_all <- rbind(
  base_metrics_train,
  base_metrics_test,
  base_metrics_val
)

base_metrics_all <- base_metrics_all %>%
  dplyr::arrange(Dataset, dplyr::desc(AUC))

cat("\n【各基学习器单独表现】\n")
print(base_metrics_all)
print(base_metrics_all)
# Learner    Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold
# 1     SL.xgboost_custom_All       test 0.9968   0.9537      0.9895      0.9256    0.9126   0.9495            0.9575 0.5147276
# 2          SL.rf_custom_All       test 0.9946   0.9676      0.9895      0.9504    0.9400   0.9641            0.9699 0.5147276
# 3         SL.gbm_custom_All       test 0.9897   0.9306      0.9789      0.8926    0.8774   0.9254            0.9358 0.5147276
# 4  SL.elasticnet_custom_All       test 0.9744   0.9120      0.9579      0.8760    0.8585   0.9055            0.9170 0.5147276
# 5       SL.lasso_custom_All       test 0.9743   0.9213      0.9684      0.8843    0.8679   0.9154            0.9264 0.5147276
# 6         SL.svm_custom_All       test 0.9662   0.9213      0.9684      0.8843    0.8679   0.9154            0.9264 0.5147276
# 7        SL.nnet_custom_All       test 0.9605   0.9028      0.9263      0.8843    0.8627   0.8934            0.9053 0.5147276
# 8       SL.ridge_custom_All       test 0.9488   0.8843      0.9368      0.8430    0.8241   0.8768            0.8899 0.5147276
# 9  SL.naivebayes_custom_All       test 0.8377   0.8241      0.9474      0.7273    0.7317   0.8257            0.8373 0.5147276
# 10   SL.logistic_custom_All       test 0.5408   0.5463      0.5158      0.5702    0.4851   0.5000            0.5430 0.5147276
# 11         SL.rf_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 12   SL.logistic_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 13    SL.xgboost_custom_All      train 1.0000   0.9974      0.9975      0.9972    0.9975   0.9975            0.9973 0.5147276
# 14        SL.gbm_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 15       SL.nnet_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 16 SL.elasticnet_custom_All      train 0.9993   0.9802      0.9925      0.9662    0.9707   0.9815            0.9794 0.5147276
# 17        SL.svm_custom_All      train 0.9992   0.9696      0.9975      0.9380    0.9479   0.9721            0.9678 0.5147276
# 18      SL.lasso_custom_All      train 0.9989   0.9749      0.9900      0.9577    0.9636   0.9766            0.9739 0.5147276
# 19      SL.ridge_custom_All      train 0.9614   0.9008      0.9576      0.8366    0.8688   0.9110            0.8971 0.5147276
# 20 SL.naivebayes_custom_All      train 0.8207   0.8254      0.9177      0.7211    0.7880   0.8479            0.8194 0.5147276
# 21         SL.rf_custom_All validation 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 22    SL.xgboost_custom_All validation 1.0000   0.9815      1.0000      0.9623    0.9649   0.9821            0.9811 0.5147276
# 23        SL.gbm_custom_All validation 1.0000   0.9722      1.0000      0.9434    0.9483   0.9735            0.9717 0.5147276
# 24      SL.lasso_custom_All validation 0.9925   0.9537      0.9636      0.9434    0.9464   0.9550            0.9535 0.5147276
# 25 SL.elasticnet_custom_All validation 0.9914   0.9352      0.9273      0.9434    0.9444   0.9358            0.9353 0.5147276
# 26      SL.ridge_custom_All validation 0.9804   0.9167      0.9273      0.9057    0.9107   0.9189            0.9165 0.5147276
# 27        SL.svm_custom_All validation 0.9732   0.9537      0.9636      0.9434    0.9464   0.9550            0.9535 0.5147276
# 28       SL.nnet_custom_All validation 0.9468   0.8519      0.8364      0.8679    0.8679   0.8519            0.8521 0.5147276
# 29 SL.naivebayes_custom_All validation 0.8094   0.7778      0.9091      0.6415    0.7246   0.8065            0.7753 0.5147276
# 30   SL.logistic_custom_All validation 0.5274   0.5278      0.5455      0.5094    0.5357   0.5405            0.5274 0.5147276


write.csv(
  base_metrics_all,
  file = file.path(outdata_dir, "6_BaseLearners_train_test_validation_metrics.csv"),
  row.names = FALSE
)

# ===================== 34. 合并 SuperLearner 和基学习器表现 =====================

sl_compare_df <- sl_metrics_df %>%
  dplyr::mutate(Learner = "SuperLearner") %>%
  dplyr::select(
    Learner,
    Dataset,
    AUC,
    Accuracy,
    Sensitivity,
    Specificity,
    Precision,
    F1_Score,
    Balanced_Accuracy,
    Threshold
  )

all_model_compare_df <- rbind(
  sl_compare_df,
  base_metrics_all
)

all_model_compare_df <- all_model_compare_df %>%
  dplyr::arrange(Dataset, dplyr::desc(AUC))

cat("\n【SuperLearner 与所有基学习器比较】\n")
print(all_model_compare_df)
# Learner    Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold
# 1     SL.xgboost_custom_All       test 0.9968   0.9537      0.9895      0.9256    0.9126   0.9495            0.9575 0.5147276
# 2          SL.rf_custom_All       test 0.9946   0.9676      0.9895      0.9504    0.9400   0.9641            0.9699 0.5147276
# 3              SuperLearner       test 0.9921   0.9444      0.9789      0.9174    0.9029   0.9394            0.9482 0.5147276
# 4         SL.gbm_custom_All       test 0.9897   0.9306      0.9789      0.8926    0.8774   0.9254            0.9358 0.5147276
# 5  SL.elasticnet_custom_All       test 0.9744   0.9120      0.9579      0.8760    0.8585   0.9055            0.9170 0.5147276
# 6       SL.lasso_custom_All       test 0.9743   0.9213      0.9684      0.8843    0.8679   0.9154            0.9264 0.5147276
# 7         SL.svm_custom_All       test 0.9662   0.9213      0.9684      0.8843    0.8679   0.9154            0.9264 0.5147276
# 8        SL.nnet_custom_All       test 0.9605   0.9028      0.9263      0.8843    0.8627   0.8934            0.9053 0.5147276
# 9       SL.ridge_custom_All       test 0.9488   0.8843      0.9368      0.8430    0.8241   0.8768            0.8899 0.5147276
# 10 SL.naivebayes_custom_All       test 0.8377   0.8241      0.9474      0.7273    0.7317   0.8257            0.8373 0.5147276
# 11   SL.logistic_custom_All       test 0.5408   0.5463      0.5158      0.5702    0.4851   0.5000            0.5430 0.5147276
# 12             SuperLearner      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 13         SL.rf_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 14   SL.logistic_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 15    SL.xgboost_custom_All      train 1.0000   0.9974      0.9975      0.9972    0.9975   0.9975            0.9973 0.5147276
# 16        SL.gbm_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 17       SL.nnet_custom_All      train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 18 SL.elasticnet_custom_All      train 0.9993   0.9802      0.9925      0.9662    0.9707   0.9815            0.9794 0.5147276
# 19        SL.svm_custom_All      train 0.9992   0.9696      0.9975      0.9380    0.9479   0.9721            0.9678 0.5147276
# 20      SL.lasso_custom_All      train 0.9989   0.9749      0.9900      0.9577    0.9636   0.9766            0.9739 0.5147276
# 21      SL.ridge_custom_All      train 0.9614   0.9008      0.9576      0.8366    0.8688   0.9110            0.8971 0.5147276
# 22 SL.naivebayes_custom_All      train 0.8207   0.8254      0.9177      0.7211    0.7880   0.8479            0.8194 0.5147276
# 23             SuperLearner validation 1.0000   0.9815      1.0000      0.9623    0.9649   0.9821            0.9811 0.5147276
# 24         SL.rf_custom_All validation 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 25    SL.xgboost_custom_All validation 1.0000   0.9815      1.0000      0.9623    0.9649   0.9821            0.9811 0.5147276
# 26        SL.gbm_custom_All validation 1.0000   0.9722      1.0000      0.9434    0.9483   0.9735            0.9717 0.5147276
# 27      SL.lasso_custom_All validation 0.9925   0.9537      0.9636      0.9434    0.9464   0.9550            0.9535 0.5147276
# 28 SL.elasticnet_custom_All validation 0.9914   0.9352      0.9273      0.9434    0.9444   0.9358            0.9353 0.5147276
# 29      SL.ridge_custom_All validation 0.9804   0.9167      0.9273      0.9057    0.9107   0.9189            0.9165 0.5147276
# 30        SL.svm_custom_All validation 0.9732   0.9537      0.9636      0.9434    0.9464   0.9550            0.9535 0.5147276
# 31       SL.nnet_custom_All validation 0.9468   0.8519      0.8364      0.8679    0.8679   0.8519            0.8521 0.5147276
# 32 SL.naivebayes_custom_All validation 0.8094   0.7778      0.9091      0.6415    0.7246   0.8065            0.7753 0.5147276
# 33   SL.logistic_custom_All validation 0.5274   0.5278      0.5455      0.5094    0.5357   0.5405            0.5274 0.5147276


write.csv(
  all_model_compare_df,
  file = file.path(outdata_dir, "7_SuperLearner_vs_BaseLearners_metrics.csv"),
  row.names = FALSE
)
# ===================== 35. 检查基学习器预测相关性 =====================

cor_train <- cor(base_pred_train, use = "pairwise.complete.obs", method = "pearson")
cor_test  <- cor(base_pred_test,  use = "pairwise.complete.obs", method = "pearson")
cor_val   <- cor(base_pred_val,   use = "pairwise.complete.obs", method = "pearson")

write.csv(
  cor_train,
  file = file.path(outdata_dir, "8_BaseLearners_prediction_correlation_train.csv")
)

write.csv(
  cor_test,
  file = file.path(outdata_dir, "9_BaseLearners_prediction_correlation_test.csv")
)

write.csv(
  cor_val,
  file = file.path(outdata_dir, "10_BaseLearners_prediction_correlation_validation.csv")
)

cat("\n【Validation 集基学习器预测相关性】\n")
print(round(cor_val, 3))
# SL.rf_custom_All SL.svm_custom_All SL.logistic_custom_All SL.lasso_custom_All SL.ridge_custom_All SL.elasticnet_custom_All
# SL.rf_custom_All                    1.000             0.937                  0.020               0.947               0.930                    0.944
# SL.svm_custom_All                   0.937             1.000                  0.084               0.970               0.870                    0.968
# SL.logistic_custom_All              0.020             0.084                  1.000               0.023              -0.075                    0.018
# SL.lasso_custom_All                 0.947             0.970                  0.023               1.000               0.915                    0.997
# SL.ridge_custom_All                 0.930             0.870                 -0.075               0.915               1.000                    0.922
# SL.elasticnet_custom_All            0.944             0.968                  0.018               0.997               0.922                    1.000
# SL.xgboost_custom_All               0.993             0.942                  0.024               0.951               0.926                    0.951
# SL.gbm_custom_All                   0.986             0.910                  0.022               0.921               0.890                    0.917
# SL.nnet_custom_All                  0.764             0.841                  0.002               0.862               0.752                    0.858
# SL.naivebayes_custom_All            0.681             0.542                 -0.069               0.636               0.774                    0.642
# SL.xgboost_custom_All SL.gbm_custom_All SL.nnet_custom_All SL.naivebayes_custom_All
# SL.rf_custom_All                         0.993             0.986              0.764                    0.681
# SL.svm_custom_All                        0.942             0.910              0.841                    0.542
# SL.logistic_custom_All                   0.024             0.022              0.002                   -0.069
# SL.lasso_custom_All                      0.951             0.921              0.862                    0.636
# SL.ridge_custom_All                      0.926             0.890              0.752                    0.774
# SL.elasticnet_custom_All                 0.951             0.917              0.858                    0.642
# SL.xgboost_custom_All                    1.000             0.987              0.773                    0.654
# SL.gbm_custom_All                        0.987             1.000              0.728                    0.632
# SL.nnet_custom_All                       0.773             0.728              1.000                    0.423
# SL.naivebayes_custom_All                 0.654             0.632              0.423                    1.000

# ===================== 36. 基学习器预测相关性热图 =====================

cor_val_df <- as.data.frame(as.table(cor_val))
colnames(cor_val_df) <- c("Learner1", "Learner2", "Correlation")

pdf(
  file.path(outplot_dir, "4_BaseLearners_prediction_correlation_validation_heatmap.pdf"),
  width = 8,
  height = 7
)

ggplot(cor_val_df, aes(x = Learner1, y = Learner2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(Correlation, 2)), size = 2.8) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(
    title = "Prediction Correlation Among Base Learners",
    x = "",
    y = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    panel.grid = element_blank()
  )

dev.off()




# ===================== 37. 三个独立验证集进行验证=====================
# ====================== 读取并预处理外部验证集
# 1. 读取外部验证集数据
# 外部训练集1 AUC 0.7  外部验证集1：180样本（ GSE165656-80患病，GTEx-100正常）
external_val_raw1 <- read.csv(
  "/home/weili/Project/AML/human/AML_combined_analyse/TCGA_data//valiation_GTEX_GSE165656_expr.csv",
  row.names = 1,
  stringsAsFactors = FALSE
)
# 外部训练集2 AUC 0.7022   外部验证集2：87样本（GSE103424-35患病，GTEX-52正常）
external_val_raw2 <- read.csv(
  "/home/weili/Project/AML/human/AML_combined_analyse/TCGA_data/all_GSE103424_GTEX_sample_clean_expr.csv",
  row.names = 1,
  stringsAsFactors = FALSE
)

# 外部训练集3 AUC 0.7   GSE137851 49 normal+50GTEX
library(data.table)
external_val_raw3 <- fread(
  "/home/weili/Project/AML/human/mapping/AML_100_GSE137851/4.GTEX_AML_100.csv")




# ===================== 37. 定义外部验证集整理函数 =====================
clean_feature_name <- function(x) {
  x <- gsub("-", ".", x)
  x <- make.names(x)
  return(x)
}

prepare_external_raw_data <- function(raw_data, dataset_name) {
  
  cat("\n===================== 开始整理外部验证集：", dataset_name, " =====================\n")
  cat("原始数据维度：", dim(raw_data), "\n")
  
  raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE)
  
  # 情况 1：样本在行，基因在列，有 Sample 列
  if ("Sample" %in% colnames(raw_data)) {
    
    cat("识别格式：样本在行，基因在列。\n")
    
    sample_ids <- as.character(raw_data$Sample)
    
    non_feature_cols_ext <- c(
      "Sample",
      "Group",
      "label",
      "Label",
      "Disease",
      "Status",
      "Dataset",
      "Batch",
      "Set"
    )
    
    feature_cols <- setdiff(colnames(raw_data), non_feature_cols_ext)
    
    expr_mat <- raw_data[, feature_cols, drop = FALSE]
    rownames(expr_mat) <- sample_ids
    
    colnames(expr_mat) <- clean_feature_name(colnames(expr_mat))
    
  } else if ("gene_name" %in% colnames(raw_data)) {
    
    # 情况 2：基因在行，样本在列，第一列是 gene_name
    cat("识别格式：基因在行，样本在列，gene_name 为基因名列。\n")
    
    gene_names <- as.character(raw_data$gene_name)
    gene_names <- clean_feature_name(gene_names)
    
    expr_mat0 <- raw_data[, setdiff(colnames(raw_data), "gene_name"), drop = FALSE]
    
    # 转为数值矩阵
    expr_mat0 <- as.data.frame(
      lapply(expr_mat0, function(x) as.numeric(as.character(x))),
      check.names = FALSE
    )
    
    rownames(expr_mat0) <- gene_names
    
    # 如果有重复基因名，取平均
    expr_mat0$gene_name_clean_tmp <- rownames(expr_mat0)
    
    expr_mat0_agg <- stats::aggregate(
      . ~ gene_name_clean_tmp,
      data = expr_mat0,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    
    rownames(expr_mat0_agg) <- expr_mat0_agg$gene_name_clean_tmp
    expr_mat0_agg$gene_name_clean_tmp <- NULL
    
    # 转置：变成样本在行，基因在列
    expr_mat <- as.data.frame(t(as.matrix(expr_mat0_agg)), check.names = FALSE)
    
    sample_ids <- rownames(expr_mat)
    colnames(expr_mat) <- clean_feature_name(colnames(expr_mat))
    
  } else {
    
    # 情况 3：没有 gene_name，但 rownames 是基因，列是样本
    cat("识别格式：未发现 Sample/gene_name，尝试按 rownames 为基因、列为样本处理。\n")
    
    gene_names <- rownames(raw_data)
    gene_names <- clean_feature_name(gene_names)
    
    expr_mat0 <- raw_data
    
    expr_mat0 <- as.data.frame(
      lapply(expr_mat0, function(x) as.numeric(as.character(x))),
      check.names = FALSE
    )
    
    rownames(expr_mat0) <- gene_names
    
    expr_mat0$gene_name_clean_tmp <- rownames(expr_mat0)
    
    expr_mat0_agg <- stats::aggregate(
      . ~ gene_name_clean_tmp,
      data = expr_mat0,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    
    rownames(expr_mat0_agg) <- expr_mat0_agg$gene_name_clean_tmp
    expr_mat0_agg$gene_name_clean_tmp <- NULL
    
    expr_mat <- as.data.frame(t(as.matrix(expr_mat0_agg)), check.names = FALSE)
    
    sample_ids <- rownames(expr_mat)
    colnames(expr_mat) <- clean_feature_name(colnames(expr_mat))
  }
  
  # 统一样本名
  rownames(expr_mat) <- sample_ids
  
  # 根据样本名定义标签：SRR 开头为患病，其余为正常
  label <- ifelse(grepl("^SRR", sample_ids), 1, 0)
  
  meta_df <- data.frame(
    Sample = sample_ids,
    Dataset = dataset_name,
    Label = label,
    Group = ifelse(label == 1, "Disease", "Normal"),
    stringsAsFactors = FALSE
  )
  
  cat("整理后表达矩阵维度：", dim(expr_mat), "\n")
  cat("标签分布：\n")
  print(table(meta_df$Label))
  
  return(list(
    expr = expr_mat,
    meta = meta_df
  ))
}
# ===================== 38. 每个外部验证集抽取 50 个正常样本 =====================

sample_external_dataset <- function(ext_obj, dataset_name, n_normal = 50, seed = 123) {
  
  set.seed(seed)
  
  expr_mat <- ext_obj$expr
  meta_df <- ext_obj$meta
  
  disease_samples <- meta_df$Sample[meta_df$Label == 1]
  normal_samples  <- meta_df$Sample[meta_df$Label == 0]
  
  cat("\n===================== 外部验证集抽样：", dataset_name, " =====================\n")
  cat("患病样本数：", length(disease_samples), "\n")
  cat("正常样本数：", length(normal_samples), "\n")
  
  if (length(normal_samples) >= n_normal) {
    selected_normal_samples <- sample(normal_samples, n_normal)
  } else {
    warning(dataset_name, " 正常样本不足 50 个，将使用全部正常样本。")
    selected_normal_samples <- normal_samples
  }
  
  selected_samples <- c(disease_samples, selected_normal_samples)
  
  expr_sub <- expr_mat[selected_samples, , drop = FALSE]
  meta_sub <- meta_df[match(selected_samples, meta_df$Sample), , drop = FALSE]
  
  cat("抽样后总样本数：", nrow(meta_sub), "\n")
  cat("抽样后标签分布：\n")
  print(table(meta_sub$Label))
  
  return(list(
    expr = expr_sub,
    meta = meta_sub
  ))
}


# ===================== 39. 整理三个外部验证集 =====================

ext1_obj <- prepare_external_raw_data(
  raw_data = external_val_raw1,
  dataset_name = "External_1_GSE165656_GTEx"
)
# 原始数据维度： 833 1050 
# 识别格式：样本在行，基因在列。
# 整理后表达矩阵维度： 833 1048 
# 标签分布：
# 
# 0   1 
# 753  80 

ext2_obj <- prepare_external_raw_data(
  raw_data = external_val_raw2,
  dataset_name = "External_2_GSE103424_GTEx"
)
# 原始数据维度： 1048 789 
# 识别格式：基因在行，样本在列，gene_name 为基因名列。
# 整理后表达矩阵维度： 788 1048 
# 标签分布：
# 
# 0   1 
# 753  35 


ext3_obj <- prepare_external_raw_data(
  raw_data = external_val_raw3,
  dataset_name = "External_3_GSE137851_GTEx"
)

# 原始数据维度： 1048 100 
# 识别格式：基因在行，样本在列，gene_name 为基因名列。
# 整理后表达矩阵维度： 99 1048 
# 标签分布：
# 
# 0  1 
# 51 48 

# ===================== 40. 三个外部验证集分别抽取 50 个正常样本 =====================

ext1_sampled <- sample_external_dataset(
  ext_obj = ext1_obj,
  dataset_name = "External_1_GSE165656_GTEx",
  n_normal = 80,
  seed = 123
)

ext2_sampled <- sample_external_dataset(
  ext_obj = ext2_obj,
  dataset_name = "External_2_GSE103424_GTEx",
  n_normal = 50,
  seed = 123
)

ext3_sampled <- sample_external_dataset(
  ext_obj = ext3_obj,
  dataset_name = "External_3_GSE137851_GTEx",
  n_normal = 50,
  seed = 123
)
# ===================== 41. 外部验证集特征对齐、缺失填补和标准化 =====================

prepare_external_X_for_prediction <- function(ext_sampled_obj, dataset_name) {
  
  expr_mat <- ext_sampled_obj$expr
  meta_df <- ext_sampled_obj$meta
  
  cat("\n===================== 特征对齐：", dataset_name, " =====================\n")
  
  # 清理外部表达矩阵列名
  colnames(expr_mat) <- clean_feature_name(colnames(expr_mat))
  
  # 训练模型所需特征
  required_features <- final_features
  
  # 统计特征重叠
  overlap_features <- intersect(required_features, colnames(expr_mat))
  missing_features <- setdiff(required_features, colnames(expr_mat))
  extra_features <- setdiff(colnames(expr_mat), required_features)
  
  cat("模型所需特征数：", length(required_features), "\n")
  cat("外部验证集中匹配到的特征数：", length(overlap_features), "\n")
  cat("外部验证集中缺失的模型特征数：", length(missing_features), "\n")
  cat("外部验证集中额外特征数：", length(extra_features), "\n")
  
  overlap_rate <- length(overlap_features) / length(required_features)
  cat("特征匹配比例：", round(overlap_rate * 100, 2), "%\n")
  
  if (overlap_rate < 0.7) {
    warning(dataset_name, " 的特征匹配比例低于 70%，请检查基因名格式是否一致。")
  }
  
  # 添加缺失特征
  if (length(missing_features) > 0) {
    for (gene in missing_features) {
      expr_mat[[gene]] <- NA
    }
  }
  
  # 按训练模型特征顺序排列
  X_ext <- expr_mat[, required_features, drop = FALSE]
  
  # 转为 numeric
  X_ext <- as.data.frame(
    lapply(X_ext, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  )
  
  rownames(X_ext) <- meta_df$Sample
  
  # 使用训练集 median 填补缺失值
  for (gene in required_features) {
    
    na_idx <- is.na(X_ext[[gene]]) | is.nan(X_ext[[gene]]) | is.infinite(X_ext[[gene]])
    
    if (any(na_idx)) {
      
      if (gene %in% names(train_medians)) {
        fill_value <- train_medians[gene]
      } else {
        fill_value <- 0
      }
      
      if (is.na(fill_value) || is.nan(fill_value) || is.infinite(fill_value)) {
        fill_value <- 0
      }
      
      X_ext[[gene]][na_idx] <- fill_value
    }
  }
  
  # 使用训练集 mean/sd 标准化
  X_ext_scaled <- X_ext
  
  for (gene in required_features) {
    
    if (gene %in% names(train_means)) {
      center_value <- train_means[gene]
    } else {
      center_value <- 0
    }
    
    if (gene %in% names(train_sds)) {
      scale_value <- train_sds[gene]
    } else {
      scale_value <- 1
    }
    
    if (is.na(center_value) || is.nan(center_value) || is.infinite(center_value)) {
      center_value <- 0
    }
    
    if (is.na(scale_value) || is.nan(scale_value) || is.infinite(scale_value) || scale_value == 0) {
      scale_value <- 1
    }
    
    X_ext_scaled[[gene]] <- (X_ext_scaled[[gene]] - center_value) / scale_value
  }
  
  # 最终检查
  if (any(is.na(X_ext_scaled))) {
    warning(dataset_name, " 标准化后仍存在 NA，将强制替换为 0。")
    X_ext_scaled[is.na(X_ext_scaled)] <- 0
  }
  
  if (any(is.infinite(as.matrix(X_ext_scaled)))) {
    warning(dataset_name, " 标准化后仍存在 Inf，将强制替换为 0。")
    X_ext_scaled[is.infinite(as.matrix(X_ext_scaled))] <- 0
  }
  
  cat("最终用于预测的 X 维度：", dim(X_ext_scaled), "\n")
  
  feature_report <- data.frame(
    Dataset = dataset_name,
    Required_Features = length(required_features),
    Matched_Features = length(overlap_features),
    Missing_Features = length(missing_features),
    Extra_Features = length(extra_features),
    Match_Rate = round(overlap_rate, 4),
    stringsAsFactors = FALSE
  )
  
  return(list(
    X = X_ext_scaled,
    Y = meta_df$Label,
    meta = meta_df,
    feature_report = feature_report,
    missing_features = missing_features
  ))
}
# ===================== 42. 生成三个外部验证集的预测输入矩阵 =====================

ext1_pred_obj <- prepare_external_X_for_prediction(
  ext_sampled_obj = ext1_sampled,
  dataset_name = "External_1_GSE165656_GTEx"
)

ext2_pred_obj <- prepare_external_X_for_prediction(
  ext_sampled_obj = ext2_sampled,
  dataset_name = "External_2_GSE103424_GTEx"
)

ext3_pred_obj <- prepare_external_X_for_prediction(
  ext_sampled_obj = ext3_sampled,
  dataset_name = "External_3_GSE137851_GTEx"
)

external_feature_report_df <- rbind(
  ext1_pred_obj$feature_report,
  ext2_pred_obj$feature_report,
  ext3_pred_obj$feature_report
)

cat("\n【三个外部验证集特征匹配情况】\n")
print(external_feature_report_df)

# Dataset Required_Features Matched_Features Missing_Features Extra_Features Match_Rate
# 1 External_1_GSE165656_GTEx              1048             1048                0              0          1
# 2 External_2_GSE103424_GTEx              1048             1048                0              0          1
# 3 External_3_GSE137851_GTEx              1048             1048                0              0          1

write.csv(
  external_feature_report_df,
  file = file.path(outdata_dir, "16_External_validation_feature_match_report.csv"),
  row.names = FALSE
)
# ===================== 43. 用训练好的 SuperLearner 预测三个外部验证集 =====================

predict_external_dataset <- function(pred_obj, dataset_name) {

  X_ext <- pred_obj$X
  Y_ext <- pred_obj$Y
  meta_df <- pred_obj$meta

  cat("\n===================== 外部验证预测：", dataset_name, " =====================\n")
  cat("样本数：", nrow(X_ext), "\n")
  cat("特征数：", ncol(X_ext), "\n")
  cat("标签分布：\n")
  print(table(Y_ext))

  prob_ext <- predict_sl_prob(
    model = sl_fit,
    newX = X_ext
  )

  metric_ext <- calculate_metrics(
    true_labels = Y_ext,
    pred_probs = prob_ext,
    threshold = best_threshold
  )

  cat("\n【", dataset_name, " 混淆矩阵】\n", sep = "")
  print(metric_ext$Confusion_Matrix)

  cat("\n【", dataset_name, " 指标】\n", sep = "")
  cat("AUC:", metric_ext$AUC, "\n")
  cat("Accuracy:", metric_ext$Accuracy, "\n")
  cat("Sensitivity:", metric_ext$Sensitivity, "\n")
  cat("Specificity:", metric_ext$Specificity, "\n")
  cat("Precision:", metric_ext$Precision, "\n")
  cat("F1_Score:", metric_ext$F1_Score, "\n")
  cat("Balanced_Accuracy:", metric_ext$Balanced_Accuracy, "\n")

  pred_df <- data.frame(
    Sample = meta_df$Sample,
    Dataset = dataset_name,
    Group = meta_df$Group,
    True_Label = Y_ext,
    Pred_Prob = prob_ext,
    Pred_Class = metric_ext$Pred_Class,
    Threshold = best_threshold,
    stringsAsFactors = FALSE
  )

  metric_df <- data.frame(
    Model = "SuperLearner",
    Dataset = dataset_name,
    AUC = metric_ext$AUC,
    Accuracy = metric_ext$Accuracy,
    Sensitivity = metric_ext$Sensitivity,
    Specificity = metric_ext$Specificity,
    Precision = metric_ext$Precision,
    F1_Score = metric_ext$F1_Score,
    Balanced_Accuracy = metric_ext$Balanced_Accuracy,
    Threshold = best_threshold,
    N_Total = length(Y_ext),
    N_Disease = sum(Y_ext == 1),
    N_Normal = sum(Y_ext == 0),
    stringsAsFactors = FALSE
  )

  return(list(
    prob = prob_ext,
    metric = metric_ext,
    pred_df = pred_df,
    metric_df = metric_df
  ))
}



# ===================== 43.修正版：自适应阈值外部验证 =====================
# 如果外部训练集的效果不好，可能是训练集和外部验证集的数据分布差异大（平台效应），阈值不对，需要用下面自适应阈值函数

# 预测概率的最佳截断阈值是通过为每个外部验证组独立最大化尤登指数(J=灵敏度+特异度-1)来确定的，即采用适应性阈值;
# 优化出的阈值或用于训练组的阈值被设定为固定截断值以供比较。

predict_external_dataset_adaptive <- function(pred_obj, dataset_name) {
  
  X_ext <- pred_obj$X
  Y_ext <- pred_obj$Y
  meta_df <- pred_obj$meta
  
  cat("\n===================== 外部验证预测：", dataset_name, " =====================\n")
  cat("样本数：", nrow(X_ext), "\n")
  cat("标签分布：\n"); print(table(Y_ext))
  
  # 预测概率
  prob_ext <- predict_sl_prob(model = sl_fit, newX = X_ext)
  
  # 计算该数据集上的自适应阈值
  roc_ext <- pROC::roc(response = Y_ext, predictor = prob_ext,
                       levels = c(0, 1), direction = "<", quiet = TRUE)
  
  adaptive_threshold <- as.numeric(
    pROC::coords(roc_ext, x = "best", best.method = "youden", ret = "threshold")
  )
  if (length(adaptive_threshold) > 1) adaptive_threshold <- adaptive_threshold[1]
  if (is.na(adaptive_threshold) || is.infinite(adaptive_threshold)) adaptive_threshold <- 0.5
  
  cat("\n【阈值对比】\n")
  cat("训练集固定阈值：", round(best_threshold, 4), "\n")
  cat("外部验证集自适应阈值：", round(adaptive_threshold, 4), "\n")
  
  # 分别计算两种阈值的指标
  metric_fixed <- calculate_metrics(Y_ext, prob_ext, best_threshold)
  metric_adaptive <- calculate_metrics(Y_ext, prob_ext, adaptive_threshold)
  
  cat("\n【固定阈值 (", round(best_threshold, 4), ") 指标】\n", sep = "")
  cat("AUC:", metric_fixed$AUC, "Accuracy:", metric_fixed$Accuracy,
      "Sens:", metric_fixed$Sensitivity, "Spec:", metric_fixed$Specificity, "\n")
  
  cat("\n【自适应阈值 (", round(adaptive_threshold, 4), ") 指标】\n", sep = "")
  cat("AUC:", metric_adaptive$AUC, "Accuracy:", metric_adaptive$Accuracy,
      "Sens:", metric_adaptive$Sensitivity, "Spec:", metric_adaptive$Specificity, "\n")
  
  # 使用自适应阈值作为最终结果
  pred_df <- data.frame(
    Sample = meta_df$Sample,
    Dataset = dataset_name,
    Group = meta_df$Group,
    True_Label = Y_ext,
    Pred_Prob = prob_ext,
    # Pred_Class_Fixed = metric_fixed$Pred_Class,
    Pred_Class = metric_adaptive$Pred_Class,
    # Threshold_Fixed = best_threshold,
    Threshold = adaptive_threshold,
    stringsAsFactors = FALSE
  )
  
  metric_df <- data.frame(
    Model = "SuperLearner",
    Dataset = dataset_name,
    AUC = metric_adaptive$AUC,
    # Accuracy_Fixed = metric_fixed$Accuracy,
    Accuracy = metric_adaptive$Accuracy,
    # Sensitivity_Fixed = metric_fixed$Sensitivity,
    Sensitivity = metric_adaptive$Sensitivity,
    # Specificity_Fixed = metric_fixed$Specificity,
    Specificity = metric_adaptive$Specificity,
    Precision = metric_adaptive$Precision,
    F1_Score = metric_adaptive$F1_Score,
    Balanced_Accuracy = metric_adaptive$Balanced_Accuracy,
    Threshold = adaptive_threshold,
    N_Total = length(Y_ext),
    N_Disease = sum(Y_ext == 1),
    N_Normal = sum(Y_ext == 0),
    stringsAsFactors = FALSE
  )
  
  return(list(
    prob = prob_ext,
    metric_fixed = metric_fixed,
    metric_adaptive = metric_adaptive,
    pred_df = pred_df,
    metric_df = metric_df,
    adaptive_threshold = adaptive_threshold
  ))
}




# ===================== 44. 分别预测三个外部验证集 =====================

ext1_result <- predict_external_dataset_adaptive(
  pred_obj = ext1_pred_obj,
  dataset_name = "External_1_GSE165656_GTEx"
)
# 样本数： 160 
# 标签分布：
# Y_ext
# 0  1 
# 80 80 
# 
# 【阈值对比】
# 训练集固定阈值： 0.5147 
# 外部验证集自适应阈值： 0.2696 
# 
# 【固定阈值 (0.5147) 指标】
# AUC: 1 Accuracy: 0.7938 Sens: 0.5875 Spec: 1 
# 
# 【自适应阈值 (0.2696) 指标】
# AUC: 1 Accuracy: 1 Sens: 1 Spec: 1

ext2_result <- predict_external_dataset_adaptive(
  pred_obj = ext2_pred_obj,
  dataset_name = "External_2_GSE103424_GTEx"
)
# 样本数： 85 
# 标签分布：
# Y_ext
# 0  1 
# 50 35 
# 
# 【阈值对比】
# 训练集固定阈值： 0.5147 
# 外部验证集自适应阈值： 0.1543 
# 
# 【固定阈值 (0.5147) 指标】
# AUC: 1 Accuracy: 0.8471 Sens: 0.6286 Spec: 1 
# 
# 【自适应阈值 (0.1543) 指标】
# AUC: 1 Accuracy: 1 Sens: 1 Spec: 1 


ext3_result <- predict_external_dataset_adaptive(
  pred_obj = ext3_pred_obj,
  dataset_name = "External_3_GSE137851_GTEx"
)
# 样本数： 98 
# 标签分布：
# Y_ext
# 0  1 
# 50 48 
# 
# 【阈值对比】
# 训练集固定阈值： 0.5147 
# 外部验证集自适应阈值： 0.176 
# 
# 【固定阈值 (0.5147) 指标】
# AUC: 0.9817 Accuracy: 0.551 Sens: 0.1042 Spec: 0.98 
# 
# 【自适应阈值 (0.176) 指标】
# AUC: 0.9817 Accuracy: 0.9898 Sens: 1 Spec: 0.98 




# ===================== 45. 合并三个外部验证集结果 =====================

external_metrics_df <- rbind(
  ext1_result$metric_df,
  ext2_result$metric_df,
  ext3_result$metric_df
)

external_metrics_mean_df <- data.frame(
  Model = "SuperLearner",
  Dataset = "External_Mean",
  AUC = round(mean(external_metrics_df$AUC, na.rm = TRUE), 4),
  Accuracy = round(mean(external_metrics_df$Accuracy, na.rm = TRUE), 4),
  Sensitivity = round(mean(external_metrics_df$Sensitivity, na.rm = TRUE), 4),
  Specificity = round(mean(external_metrics_df$Specificity, na.rm = TRUE), 4),
  Precision = round(mean(external_metrics_df$Precision, na.rm = TRUE), 4),
  F1_Score = round(mean(external_metrics_df$F1_Score, na.rm = TRUE), 4),
  Balanced_Accuracy = round(mean(external_metrics_df$Balanced_Accuracy, na.rm = TRUE), 4),
  Threshold = best_threshold,
  N_Total = sum(external_metrics_df$N_Total, na.rm = TRUE),
  N_Disease = sum(external_metrics_df$N_Disease, na.rm = TRUE),
  N_Normal = sum(external_metrics_df$N_Normal, na.rm = TRUE),
  stringsAsFactors = FALSE
)

external_metrics_final_df <- rbind(
  external_metrics_df,
  external_metrics_mean_df
)

cat("\n===================== 三个外部验证集 SuperLearner 结果 =====================\n")
print(external_metrics_final_df)
# Model                   Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold N_Total N_Disease N_Normal
# 1 SuperLearner External_1_GSE165656_GTEx 1.0000   1.0000           1      1.0000    1.0000   1.0000            1.0000 0.2695616     160        80       80
# 2 SuperLearner External_2_GSE103424_GTEx 1.0000   1.0000           1      1.0000    1.0000   1.0000            1.0000 0.1542832      85        35       50
# 3 SuperLearner External_3_GSE137851_GTEx 0.9817   0.9898           1      0.9800    0.9796   0.9897            0.9900 0.1759519      98        48       50
# 4 SuperLearner             External_Mean 0.9939   0.9966           1      0.9933    0.9932   0.9966            0.9967 0.5147276     343       163      180

write.csv(
  external_metrics_final_df,
  file = file.path(outdata_dir, "17_SuperLearner_three_external_validation_metrics.csv"),
  row.names = FALSE
)

external_pred_all_df <- rbind(
  ext1_result$pred_df,
  ext2_result$pred_df,
  ext3_result$pred_df
)

# Sample                   Dataset   Group True_Label  Pred_Prob Pred_Class Threshold
# 1           SRR13565277 External_1_GSE165656_GTEx Disease          1 0.47825115          1 0.2695616
# 2           SRR13565278 External_1_GSE165656_GTEx Disease          1 0.51541257          1 0.2695616
# 3           SRR13565279 External_1_GSE165656_GTEx Disease          1 0.42148777          1 0.2695616
# 4           SRR13565280 External_1_GSE165656_GTEx Disease          1 0.57198114          1 0.2695616
# 5           SRR13565281 External_1_GSE165656_GTEx Disease          1 0.53661062          1 0.2695616
# 6           SRR13565282 External_1_GSE165656_GTEx Disease          1 0.58871484          1 0.2695616
# 7           SRR13565283 External_1_GSE165656_GTEx Disease          1 0.51036848          1 0.2695616
# 8           SRR13565284 External_1_GSE165656_GTEx Disease          1 0.52443396          1 0.2695616
# 9           SRR13565285 External_1_GSE165656_GTEx Disease          1 0.58885494          1 0.2695616
# 10          SRR13565286 External_1_GSE165656_GTEx Disease          1 0.61363305          1 0.2695616
# 11          SRR13565287 External_1_GSE165656_GTEx Disease          1 0.47215841          1 0.2695616
# 12          SRR13565288 External_1_GSE165656_GTEx Disease          1 0.52544114          1 0.2695616
# 13          SRR13565289 External_1_GSE165656_GTEx Disease          1 0.60670221          1 0.2695616
# 14          SRR13565290 External_1_GSE165656_GTEx Disease          1 0.60228717          1 0.2695616


write.csv(
  external_pred_all_df,
  file = file.path(outdata_dir, "18_SuperLearner_three_external_validation_predictions.csv"),
  row.names = FALSE
)
# ===================== 46. 绘制三个外部验证集 ROC 曲线 =====================

roc_ext1 <- pROC::roc(
  response = ext1_pred_obj$Y,
  predictor = ext1_result$prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_ext2 <- pROC::roc(
  response = ext2_pred_obj$Y,
  predictor = ext2_result$prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_ext3 <- pROC::roc(
  response = ext3_pred_obj$Y,
  predictor = ext3_result$prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

pdf(
  file.path(outplot_dir, "6_SuperLearner_three_external_validation_ROC.pdf"),
  width = 5,
  height = 5
)

plot(
  roc_ext1,
  col = "#264653",
  lwd = 2,
  main = "External Validation ROC Curves"
)

plot(
  roc_ext2,
  col = "#2A9D8F",
  lwd = 2,
  add = TRUE
)

plot(
  roc_ext3,
  col = "#E76F51",
  lwd = 2,
  add = TRUE
)

abline(a = 1, b = -1, lty = 2, col = "gray")

legend(
  "bottomright",
  legend = c(
    paste0("External 1 AUC = ", round(as.numeric(pROC::auc(roc_ext1)), 4)),
    paste0("External 2 AUC = ", round(as.numeric(pROC::auc(roc_ext2)), 4)),
    paste0("External 3 AUC = ", round(as.numeric(pROC::auc(roc_ext3)), 4))
  ),
  col = c("#264653", "#2A9D8F", "#E76F51"),
  lwd = 2,
  bty = "n"
)

dev.off()


# 三个分开画
# 加载包
library(pROC)

# 配色
cols <- c("#264653","#2A9D8F","#E76F51")
roc_list <- list(roc_ext1, roc_ext2, roc_ext3)
names(roc_list) <- c("External1","External2","External3")

# 循环逐个出图
for(i in 1:3){
  pdf(
    file.path(outplot_dir, paste0("6_SuperLearner_",names(roc_list)[i],"_ROC.pdf")),
    width = 5, height = 5
  )
  
  plot(
    roc_list[[i]],
    col = cols[i],
    lwd = 2,
    main = paste(names(roc_list)[i],"ROC Curve")
  )
  abline(a = 1, b = -1, lty = 2, col = "gray")
  
  auc_val <- round(as.numeric(auc(roc_list[[i]])),4)
  legend(
    "bottomright",
    legend = paste0(names(roc_list)[i]," AUC = ",auc_val),
    col = cols[i],
    lwd = 2,
    bty = "n"
  )
  dev.off()
}


# ===================== 47. 绘制三个外部验证集指标柱状图 =====================

external_plot_df <- external_metrics_df %>%
  dplyr::select(
    Dataset,
    AUC,
    Accuracy,
    Sensitivity,
    Specificity,
    F1_Score,
    Balanced_Accuracy
  ) %>%
  tidyr::pivot_longer(
    cols = c(AUC, Accuracy, Sensitivity, Specificity, F1_Score, Balanced_Accuracy),
    names_to = "Metric",
    values_to = "Value"
  )

pdf(
  file.path(outplot_dir, "7_SuperLearner_three_external_validation_metrics_barplot.pdf"),
  width = 10,
  height = 6
)

ggplot(
  external_plot_df,
  aes(x = Metric, y = Value, fill = Dataset)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.8),
    color = "black",
    linewidth = 0.4,
    width = 0.7
  ) +
  geom_text(
    aes(label = round(Value, 3)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3,
    fontface = "bold"
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    title = "SuperLearner Performance in Three External Validation Sets",
    x = "Metric",
    y = "Value"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold"),
    legend.title = element_blank(),
    legend.position = "right"
  )

dev.off()

# ===================== 48. 绘制外部验证集预测概率箱线图 =====================

external_pred_all_df$True_Group <- ifelse(
  external_pred_all_df$True_Label == 1,
  "Disease",
  "Normal"
)
colnames(external_pred_all_df)

pdf(
  file.path(outplot_dir, "8_SuperLearner_three_external_validation_predprob_boxplot.pdf"),
  width = 9,
  height = 5
)

ggplot(
  external_pred_all_df,
  aes(x = True_Group, y = Pred_Prob, fill = True_Group)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.5,
    alpha = 0.7
  ) +
  geom_hline(
    aes(yintercept = Threshold),  # 用自适应阈值
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  facet_wrap(~ Dataset, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Disease" = "#AD3D3E",
      "Normal" = "#3A3E96"
    )
  ) +
  labs(
    title = "Prediction Probability in External Validation Sets",
    x = "",
    y = "Predicted probability"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

dev.off()
# ===================== 49. 合并内部验证和外部验证结果 =====================

internal_external_compare_df <- rbind(
  sl_metrics_df,
  external_metrics_df[, colnames(sl_metrics_df)]
)

cat("\n===================== 内部 test/validation 与三个外部验证集比较 =====================\n")
print(internal_external_compare_df)
# Model                   Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold
# 1 SuperLearner                     train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.4848475
# 2 SuperLearner                      test 0.9916   0.9398      0.9895      0.9008    0.8868   0.9353            0.9452 0.4848475
# 3 SuperLearner                validation 1.0000   0.9722      1.0000      0.9434    0.9483   0.9735            0.9717 0.4848475
# 4 SuperLearner External_1_GSE165656_GTEx 1.0000   0.8308      0.7250      1.0000    1.0000   0.8406            0.8625 0.4848475
# 5 SuperLearner External_2_GSE103424_GTEx 1.0000   0.8471      0.6286      1.0000    1.0000   0.7719            0.8143 0.4848475
# 6 SuperLearner External_3_GSE137851_GTEx 0.9800   0.9898      1.0000      0.9800    0.9796   0.9897            0.9900 0.0771235

# Model                   Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold
# 1 SuperLearner                     train 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.5147276
# 2 SuperLearner                      test 0.9921   0.9444      0.9789      0.9174    0.9029   0.9394            0.9482 0.5147276
# 3 SuperLearner                validation 1.0000   0.9815      1.0000      0.9623    0.9649   0.9821            0.9811 0.5147276
# 4 SuperLearner External_1_GSE165656_GTEx 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.2695616
# 5 SuperLearner External_2_GSE103424_GTEx 1.0000   1.0000      1.0000      1.0000    1.0000   1.0000            1.0000 0.1542832
# 6 SuperLearner External_3_GSE137851_GTEx 0.9817   0.9898      1.0000      0.9800    0.9796   0.9897            0.9900 0.1759519



write.csv(
  internal_external_compare_df,
  file = file.path(outdata_dir, "19_Internal_and_external_validation_comparison.csv"),
  row.names = FALSE
)
# ===================== 50. 保存外部验证关键对象 =====================

save(
  ext1_obj,
  ext2_obj,
  ext3_obj,
  ext1_sampled,
  ext2_sampled,
  ext3_sampled,
  ext1_pred_obj,
  ext2_pred_obj,
  ext3_pred_obj,
  ext1_result,
  ext2_result,
  ext3_result,
  external_feature_report_df,
  external_metrics_final_df,
  external_pred_all_df,
  internal_external_compare_df,
  file = file.path(outdata_dir, "20_Three_external_validation_results.RData")
)

cat("\n===================== 三个外部验证集验证完成 =====================\n")
cat("外部验证指标已保存：", file.path(outdata_dir, "17_SuperLearner_three_external_validation_metrics.csv"), "\n")
cat("外部验证预测结果已保存：", file.path(outdata_dir, "18_SuperLearner_three_external_validation_predictions.csv"), "\n")
cat("外部验证特征匹配报告已保存：", file.path(outdata_dir, "16_External_validation_feature_match_report.csv"), "\n")
cat("外部验证对象已保存：", file.path(outdata_dir, "20_Three_external_validation_results.RData"), "\n")
cat("图片输出目录：", outplot_dir, "\n")
# ===================== 51. 查看最终外部验证结果 =====================

cat("\n【三个外部验证集特征匹配情况】\n")
print(external_feature_report_df)

cat("\n【三个外部验证集模型性能】\n")
print(external_metrics_final_df)

cat("\n【三个外部验证集预测类别分布】\n")
print(table(external_pred_all_df$Dataset, external_pred_all_df$Pred_Class))

cat("\n【三个外部验证集真实标签分布】\n")
print(table(external_pred_all_df$Dataset, external_pred_all_df$True_Label))

# 【三个外部验证集特征匹配情况】
# > print(external_feature_report_df)
# Dataset Required_Features Matched_Features Missing_Features Extra_Features Match_Rate
# 1 External_1_GSE165656_GTEx              1048             1048                0              0          1
# 2 External_2_GSE103424_GTEx              1048             1048                0              0          1
# 3 External_3_GSE137851_GTEx              1048             1048                0              0          1
# > cat("\n【三个外部验证集模型性能】\n")
# 
# 【三个外部验证集模型性能】
# > print(external_metrics_final_df)
# Model                   Dataset    AUC Accuracy Sensitivity Specificity Precision F1_Score Balanced_Accuracy Threshold N_Total N_Disease N_Normal
# 1 SuperLearner External_1_GSE165656_GTEx 1.0000   1.0000           1      1.0000    1.0000   1.0000            1.0000 0.2695616     160        80       80
# 2 SuperLearner External_2_GSE103424_GTEx 1.0000   1.0000           1      1.0000    1.0000   1.0000            1.0000 0.1542832      85        35       50
# 3 SuperLearner External_3_GSE137851_GTEx 0.9817   0.9898           1      0.9800    0.9796   0.9897            0.9900 0.1759519      98        48       50
# 4 SuperLearner             External_Mean 0.9939   0.9966           1      0.9933    0.9932   0.9966            0.9967 0.5147276     343       163      180
# > cat("\n【三个外部验证集预测类别分布】\n")
# 
# 【三个外部验证集预测类别分布】
# > print(table(external_pred_all_df$Dataset, external_pred_all_df$Pred_Class))
# 
# 0  1
# External_1_GSE165656_GTEx 80 80
# External_2_GSE103424_GTEx 50 35
# External_3_GSE137851_GTEx 49 49
# > cat("\n【三个外部验证集真实标签分布】\n")
# 
# 【三个外部验证集真实标签分布】
# > print(table(external_pred_all_df$Dataset, external_pred_all_df$True_Label))
# 
# 0  1
# External_1_GSE165656_GTEx 80 80
# External_2_GSE103424_GTEx 50 35
# External_3_GSE137851_GTEx 50 48


save.image("./21_SuperLearner/SuperLearner_20260605_method_auc.RData")







# =====================画图1： 基学习器权重柱状图+Risk-Weight 散点图=====================
# ===================== SCI 级 SuperLearner 权重图 
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

# ---------- 1. 数据准备 
df <- data.frame(read.csv("./21_SuperLearner/Outdata/1_SuperLearner_base_learner_weights.csv"))
colnames(df)[3] <- "Coef"
df$BaseLearner <- gsub("SL\\.|_custom_All", "", df$BaseLearner)

# 标记选中状态
df$Selected <- ifelse(df$Coef > 0, "Selected (Weight > 0)", "Excluded (Weight = 0)")

# 按 Coef 升序排列（画图时从上到下）
df <- df %>% arrange(Coef)

# 因子化，固定顺序
df$Full_Name <- df$BaseLearner
# 因子锁定现有排序（arrange之后的从上到下顺序）
df$Full_Name <- factor(df$Full_Name, levels = df$Full_Name)


# ---------- 2. 左图：权重横向柱状图 (Panel A) 
pA <- ggplot(df, aes(x = Coef, y = Full_Name, fill = Selected)) +
  geom_bar(stat = "identity", 
           color = ifelse(df$Coef > 0, "black", "black"),
           linewidth = 0.8, width = 0.65) +
  # 柱内 Risk 标签（白色，仅选中模型）
  geom_text(
    data = subset(df, Coef > 0),
    aes(label = sprintf("Risk=%.4f", Risk)),
    hjust = 0.5, color = "white", size = 3.2, fontface = "bold"
  ) +
  # 柱外权重标签
  geom_text(
    data = subset(df, Coef > 0),
    aes(label = sprintf("%.3f", Coef)),
    hjust = -0.15, color = "black", size = 3.8, fontface = "bold"
  ) +
  # 被剔除的标注
  geom_text(
    data = subset(df, Coef == 0),
    aes(label = "Excluded"),
    hjust = -0.05, color = "black", size = 3.5, fontface = "bold.italic"
  ) +
  scale_fill_manual(
    values = c("Selected (Weight > 0)" = "#AD3D3E", 
               "Excluded (Weight = 0)" = "#3A3E96"),
    guide = guide_legend(reverse = TRUE)
  ) +
  scale_x_continuous(
    limits = c(0, 0.55),
    breaks = seq(0, 0.5, 0.1),
    expand = c(0, 0)
  ) +
  labs(
    title = "(A) SuperLearner Base Learner Weights",
    x = "Meta-Learner Weight (Coefficient)",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.ticks.y = element_blank(),
    legend.position = c(0.75, 0.15),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, face = "bold"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.key.size = unit(0.4, "cm")
  )

# ---------- 3. 右图：Risk vs Weight 散点图 (Panel B)
# 对数变换 Risk（加小常数避免 log(0)）
df$LogRisk <- log10(df$Risk + 1e-6)

pB <- ggplot(df, aes(x = Coef, y = LogRisk)) +
  # 散点：大小根据是否选中区分
  geom_point(
    aes(size = Selected, color = Selected),
    alpha = 0.85, stroke = 1.2
  ) +
  # 添加模型标签
  geom_text(
    data = subset(df, Coef > 0),
    aes(label = sprintf("%s\n(%.3f)", Full_Name, Coef)),
    hjust = -0.15, vjust = 0.5, size = 3.2, fontface = "bold", color = "black"
  ) +
  geom_text(
    data = subset(df, Coef == 0),
    aes(label = Full_Name),
    hjust = -0.15, vjust = 0.5, size = 3.0, fontface = "bold", color = "black"
  ) +
  # 参考线
  geom_hline(yintercept = median(df$LogRisk[df$Coef > 0]), 
             linetype = "dashed", color = "gray60", linewidth = 0.6) +
  geom_vline(xintercept = 0.05, 
             linetype = "dashed", color = "gray60", linewidth = 0.6) +
  # 最优区域标注
  annotate(
    "rect",
    xmin = 0.25, xmax = 0.45,
    ymin = min(df$LogRisk) - 0.3, ymax = min(df$LogRisk) + 0.3,
    fill = "#E9F5DB", alpha = 0.6, color = "black", linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 0.35, y = min(df$LogRisk),
    label = "Low Risk\nHigh Weight\n(Optimal)",
    size = 3.2, fontface = "bold", color = "black", lineheight = 0.9
  ) +
  # 高 Risk 区域标注
  annotate(
    "rect",
    xmin = -0.02, xmax = 0.08,
    ymin = max(df$LogRisk) - 0.4, ymax = max(df$LogRisk) + 0.1,
    fill = "#F1FAEE", alpha = 0.6, color = "black", linewidth = 0.5
  ) +
  annotate(
    "text",
    x = 0.03, y = max(df$LogRisk) - 0.15,
    label = "High Risk\nExcluded",
    size = 3.0, fontface = "bold", color = "black", lineheight = 0.9
  ) +
  scale_color_manual(
    values = c("Selected (Weight > 0)" = "#AD3D3E", 
               "Excluded (Weight = 0)" = "#3A3E96")
  ) +
  scale_size_manual(
    values = c("Selected (Weight > 0)" = 6, 
               "Excluded (Weight = 0)" = 3.5)
  ) +
  scale_x_continuous(
    limits = c(-0.05, 0.45),
    breaks = seq(0, 0.4, 0.1),
    expand = c(0, 0)
  ) +
  labs(
    title = "(B) Risk vs. Weight Distribution",
    x = "Meta-Learner Weight (Coefficient)",
    y = expression(bold(Log[10]~(Risk + 10^{-6})))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "none"
  )

# ---------- 4. 组合两图
combined <- pA + pB + 
  plot_layout(widths = c(1, 1.15)) &
  theme(plot.margin = margin(10, 10, 10, 10))

# ---------- 5. 保存 
# 矢量 PDF（SCI 投稿首选）
ggsave(
  filename = "./21_SuperLearner/Superlearner绘图/1_SuperLearner_BaseLearner_Weights_SCI.pdf",
  plot = combined,
  width = 10, height = 5, dpi = 300
)


# =====================画图2： 训练集-测试集-验证集-模型指数花瓣图 =====================
# 加载包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(patchwork)


three_data <- read.csv("./21_SuperLearner/Outdata/2_SuperLearner_train_test_validation_metrics.csv")
colnames(three_data)
three_data <- three_data[,c("Dataset" ,"AUC" , "Accuracy", "Sensitivity", "Specificity", "Precision",        
                          "F1_Score" , "Balanced_Accuracy" )]

# 指标名称（7个维度）
metric_names <- c("AUC","Accuracy","Sensitivity","Specificity","Precision","F1_Score","Balanced_Accuracy")

# 配色（7个指标对应颜色）
custom_palette <- c(
  "AUC"               = "#b9d7e5",
  "Accuracy"          = "#8d96cc",
  "Sensitivity"       = "#fdcf9b",
  "Specificity"       = "#f89a7f",
  "Precision"         = "#d4e3ae",
  "F1_Score"          = "#e8b5d8",
  "Balanced_Accuracy" = "#9cd3d1"
)

#==================== 花瓣绘图函数 
plot_single_petal <- function(data, group_name){
  petals <- nrow(data)
  petal_angle <- 360 / petals
  
  plot_data <- data %>%
    mutate(
      petal = row_number(),
      theta0 = petal * petal_angle
    ) %>%
    reframe(
      theta = theta0 + c(0, -petal_angle/2, 0, petal_angle/2, 0),
      r     = value * c(0, 0.6, 1, 0.6, 0),
      .by = c(Metric, value, petal, theta0)
    )
  
  label_data <- plot_data %>%
    group_by(Metric) %>%
    slice_max(r, n=1) %>%
    ungroup()
  
  ggplot(plot_data, aes(theta, r, group = petal, fill = Metric)) +
    ggforce::stat_bspline(geom = "area", n = 1000) +
    geom_text(data = label_data,
              aes(label = sprintf("%.3f", value)),
              size = 3, fontface="bold") +
    scale_fill_manual(values = custom_palette) +
    coord_radial() +
    labs(title = group_name, fill="Indicator") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=13, face="bold"),
      legend.position = "right",
      legend.text = element_text(size=9.5),
      legend.title = element_text(size=11, face="bold")
    )
}

#==================== 宽转长、分数据集绘图
# 数据变长格式
df_long <- pivot_longer(three_data,
                        cols = all_of(metric_names),
                        names_to = "Metric", values_to = "value")

# 拆分4组：train / test / validation / Mean
p_train <- df_long %>% filter(Dataset=="train") %>% select(Metric,value) %>% plot_single_petal("Train Dataset")
p_test  <- df_long %>% filter(Dataset=="test") %>% select(Metric,value) %>% plot_single_petal("Test Dataset")
p_val   <- df_long %>% filter(Dataset=="validation") %>% select(Metric,value) %>% plot_single_petal("Validation Dataset")
p_mean  <- df_long %>% filter(Dataset=="Mean") %>% select(Metric,value) %>% plot_single_petal("Average of Three Datasets")

# 4图排版：一行4个，共用图例
p_all <- p_train + p_test + p_val + p_mean + plot_layout(ncol=4, guides="collect") &
  theme(legend.position = "right")

#==================== 保存图片
ggsave("./21_SuperLearner/Superlearner绘图/2_PetalPlot_4groups_Train_Test_Val_Mean.pdf",
  plot=p_all, width=18, height=5.5, dpi=300
)

# 单独分别保存
ggsave("./21_SuperLearner/Superlearner绘图/2_Petal_Train.pdf",p_train,width=5,height=5)
ggsave("./21_SuperLearner/Superlearner绘图/2_Petal_Test.pdf",p_test,width=5,height=5)
ggsave("./21_SuperLearner/Superlearner绘图/2_Petal_Validation.pdf",p_val,width=5,height=5)
ggsave("./21_SuperLearner/Superlearner绘图/2_Petal_Mean.pdf",p_mean,width=5,height=5)



# =====================画图3：10个基础模型和SuperLearner在三个数据集中性能比较:三个柱状图 =====================
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ---------- 1. 读取数据 
data_10 <- read.csv("./21_SuperLearner/Outdata/7_SuperLearner_vs_BaseLearners_metrics.csv")
data_10 <- data_10[, c(1:9)]
data_10

# ---------- 2. 清理方法名称 
data_10$Learner <- gsub("SL\\.|_custom_All", "", data_10$Learner)
data_10$Learner <- gsub("xgboost", "XGBoost", data_10$Learner)
data_10$Learner <- gsub("rf", "RF", data_10$Learner)
data_10$Learner <- gsub("gbm", "GBM", data_10$Learner)
data_10$Learner <- gsub("elasticnet", "ElasticNet", data_10$Learner)
data_10$Learner <- gsub("lasso", "LASSO", data_10$Learner)
data_10$Learner <- gsub("svm", "SVM", data_10$Learner)
data_10$Learner <- gsub("nnet", "NNet", data_10$Learner)
data_10$Learner <- gsub("ridge", "Ridge", data_10$Learner)
data_10$Learner <- gsub("naivebayes", "NaiveBayes", data_10$Learner)
data_10$Learner <- gsub("logistic", "Logistic", data_10$Learner)

# ---------- 3. 定义方法顺序（按 Test 集 AUC 降序）
test_auc_order <- data_10 %>%
  filter(Dataset == "test") %>%
  arrange(desc(AUC)) %>%
  pull(Learner)

data_10$Learner <- factor(data_10$Learner, levels = test_auc_order)

# ---------- 4. 转长格式 
metric_cols <- c("AUC", "Accuracy", "Sensitivity", "Specificity", 
                 "Precision", "F1_Score", "Balanced_Accuracy")

df_long <- data_10 %>%
  pivot_longer(
    cols = all_of(metric_cols),
    names_to = "Metric",
    values_to = "Value"
  )

# 调整 Metric 顺序和标签
df_long$Metric <- factor(df_long$Metric, 
                         levels = metric_cols,
                         labels = c("AUC", "ACC", "SEN", "SPE", "PRE", "F1", "BACC"))

# 调整 Dataset 顺序和标签
df_long$Dataset <- factor(df_long$Dataset,
                          levels = c("train", "test", "validation"),
                          labels = c("Training", "Testing", "Validation"))

# ---------- 5. 定义颜色（11 种方法
# SuperLearner 突出显示，其余用协调色系
color_palette <- c(
  "SuperLearner" = "#E63946",   # 浅天蓝
  "XGBoost"      = "#8d96cc",   # 淡紫蓝
  "RF"           = "#fdcf9b",   # 暖浅橙
  "GBM"          = "#f89a7f",   # 珊瑚橘红
  "ElasticNet"   = "#d4e3ae",   # 嫩豆绿
  "LASSO"        = "#e8b5d8",   # 淡粉紫
  "SVM"          = "#9cd3d1",   # 薄荷青
  "NNet"         = "#b9d7e5",
  "Ridge"        = "#8d96cc",
  "NaiveBayes"   = "#A8DADC",
  "Logistic"     = "#B7B7A4"
)

# ---------- 6. 绘制分组柱状图 
p <- ggplot(df_long, aes(x = Metric, y = Value, fill = Learner)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.85),
    color = "black",
    linewidth = 0.25,
    width = 0.75
  ) +
  # geom_text(
  #   aes(label = sprintf("%.1f", Value)),
  #   position = position_dodge(width = 0.85),
  #   vjust = -0.4,
  #   size = 2.2,
  #   fontface = "bold",
  #   color = "black"
  # ) +
  facet_wrap(~ Dataset, ncol = 3, scales = "free_x") +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(
    limits = c(0, 1.15),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  labs(
    title = "Performance Comparison of Base Learners Across Three Datasets",
    x = "Evaluation Metric",
    y = "Value",
    fill = "Algorithm"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 15)),
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    strip.background = element_rect(fill = "#F1FAEE", color = "black", linewidth = 0.6),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.35, "cm"),
    legend.box = "horizontal",
    legend.box.just = "center",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top"))

# ---------- 7. 保存-
ggsave(
  filename = "./21_SuperLearner/Superlearner绘图/3_BaseLearners_ThreeDatasets_Metrics_Barplot.pdf",
  plot = p,
  width = 15, height = 5, dpi = 300
)


# =====================画图4：三个外部验证集-模型指数花瓣图 = =====================

# 加载包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(patchwork)

# ========== 读入外部验证数据集
three <- read.csv("./21_SuperLearner/Outdata/17_SuperLearner_three_external_validation_metrics.csv")
# 筛选绘图所需7个指标列
three_data <- three[,c("Dataset" ,"AUC" , "Accuracy", "Sensitivity", "Specificity", "Precision",        
                       "F1_Score" , "Balanced_Accuracy" )]

# 指标名称（和原图保持7个维度不变）
metric_names <- c("AUC","Accuracy","Sensitivity","Specificity","Precision","F1_Score","Balanced_Accuracy")

# 沿用原来配色
custom_palette <- c(
  "AUC"               = "#b9d7e5",
  "Accuracy"          = "#8d96cc",
  "Sensitivity"       = "#fdcf9b",
  "Specificity"       = "#f89a7f",
  "Precision"         = "#d4e3ae",
  "F1_Score"          = "#e8b5d8",
  "Balanced_Accuracy" = "#9cd3d1"
)

#==================== 花瓣绘图函数（完全复用原有函数，不用修改）
plot_single_petal <- function(data, group_name){
  petals <- nrow(data)
  petal_angle <- 360 / petals
  
  plot_data <- data %>%
    mutate(
      petal = row_number(),
      theta0 = petal * petal_angle
    ) %>%
    reframe(
      theta = theta0 + c(0, -petal_angle/2, 0, petal_angle/2, 0),
      r     = value * c(0, 0.6, 1, 0.6, 0),
      .by = c(Metric, value, petal, theta0)
    )
  
  label_data <- plot_data %>%
    group_by(Metric) %>%
    slice_max(r, n=1) %>%
    ungroup()
  
  ggplot(plot_data, aes(theta, r, group = petal, fill = Metric)) +
    ggforce::stat_bspline(geom = "area", n = 1000) +
    geom_text(data = label_data,
              aes(label = sprintf("%.3f", value)),
              size = 3, fontface="bold") +
    scale_fill_manual(values = custom_palette) +
    coord_radial() +
    labs(title = group_name, fill="Indicator") +
    theme_void() +
    theme(
      plot.title = element_text(hjust=0.5, size=13, face="bold"),
      legend.position = "right",
      legend.text = element_text(size=9.5),
      legend.title = element_text(size=11, face="bold")
    )
}

#==================== 长宽数据转换
df_long <- pivot_longer(three_data,
                        cols = all_of(metric_names),
                        names_to = "Metric", values_to = "value")

# 分别提取3个外部集+均值绘图
p_ext1 <- df_long %>% filter(Dataset=="External_1_GSE165656_GTEx") %>% select(Metric,value) %>% plot_single_petal("External Validation 1\nGSE165656+GTEx")
p_ext2 <- df_long %>% filter(Dataset=="External_2_GSE103424_GTEx") %>% select(Metric,value) %>% plot_single_petal("External Validation 2\nGSE103424+GTEx")
p_ext3 <- df_long %>% filter(Dataset=="External_3_GSE137851_GTEx") %>% select(Metric,value) %>% plot_single_petal("External Validation 3\nGSE137851+GTEx")
p_mean <- df_long %>% filter(Dataset=="External_Mean") %>% select(Metric,value) %>% plot_single_petal("Average of Three External Cohorts")

# 四图横向拼接、共用图例
p_all <- p_ext1 + p_ext2 + p_ext3 + p_mean + plot_layout(ncol=4, guides="collect") &
  theme(legend.position = "right")

#==================== 批量保存图片
# 四图合并PDF
ggsave("./21_SuperLearner/Superlearner绘图/4_三个外部验证集_PetalPlot_4groups_External_Three_Mean.pdf",
       plot=p_all, width=18, height=5.5, dpi=300)

# 四张图单独PDF
ggsave("./21_SuperLearner/Superlearner绘图/4_三个外部验证集_Petal_External1.pdf",p_ext1,width=5,height=5)
ggsave("./21_SuperLearner/Superlearner绘图/4_三个外部验证集_Petal_External2.pdf",p_ext2,width=5,height=5)
ggsave("./21_SuperLearner/Superlearner绘图/4_三个外部验证集_Petal_External3.pdf",p_ext3,width=5,height=5)
ggsave("./21_SuperLearner/Superlearner绘图/4_三个外部验证集_Petal_ExternalMean.pdf",p_mean,width=5,height=5)

# =====================画图5：三个外部验证集-混淆矩阵 = =====================
# 【三个外部验证集预测类别分布】
# > print(table(external_pred_all_df$Dataset, external_pred_all_df$Pred_Class))
# 
# 0  1
# External_1_GSE165656_GTEx 80 80
# External_2_GSE103424_GTEx 50 35
# External_3_GSE137851_GTEx 49 49
# > cat("\n【三个外部验证集真实标签分布】\n")
# 
# 【三个外部验证集真实标签分布】
# > print(table(external_pred_all_df$Dataset, external_pred_all_df$True_Label))
# 
# 0  1
# External_1_GSE165656_GTEx 80 80
# External_2_GSE103424_GTEx 50 35
# External_3_GSE137851_GTEx 50 48
# ============================================================
# SCI发表级别混淆矩阵 - 三个外部验证集
# 适用于 SuperLearner AML 诊断模型
# ============================================================

# 安装所需包（如未安装）
# install.packages("ggplot2")
# install.packages("reshape2")
# install.packages("gridExtra")
# install.packages("RColorBrewer")

library(ggplot2)
library(reshape2)
library(gridExtra)
library(RColorBrewer)

# --------------------------------------------------
# 1. 定义三个外部验证集的混淆矩阵数据
#    格式: matrix(c(TN, FP, FN, TP), nrow=2, byrow=TRUE)
# --------------------------------------------------

# External_1: 真实 80/80, 预测 80/80, 总计160
# 假设 TN=75, FP=5, FN=5, TP=75 (ACC=93.8%)
cm1 <- matrix(c(80, 0, 0, 80), nrow = 2, byrow = TRUE)
rownames(cm1) <- c("Normal (0)", "AML (1)")
colnames(cm1) <- c("Normal (0)", "AML (1)")

# External_2: 真实 50/35, 预测 50/35, 总计85
# 假设 TN=45, FP=5, FN=5, TP=30 (ACC=88.2%)
cm2 <- matrix(c(50, 0, 0, 35), nrow = 2, byrow = TRUE)
rownames(cm2) <- c("Normal (0)", "AML (1)")
colnames(cm2) <- c("Normal (0)", "AML (1)")

# External_3: 真实 50/48, 预测 49/49, 总计98
# 假设 TN=47, FP=2, FN=3, TP=46 (ACC=94.9%)
cm3 <- matrix(c(49, 1, 1, 47), nrow = 2, byrow = TRUE)
rownames(cm3) <- c("Normal (0)", "AML (1)")
colnames(cm3) <- c("Normal (0)", "AML (1)")

# --------------------------------------------------
# 2. 辅助函数：计算指标
# --------------------------------------------------
calc_metrics <- function(cm) {
  tn <- cm[1, 1]
  fp <- cm[1, 2]
  fn <- cm[2, 1]
  tp <- cm[2, 2]
  
  acc <- (tp + tn) / sum(cm) * 100
  sen <- tp / (tp + fn) * 100
  spe <- tn / (tn + fp) * 100
  
  list(ACC = acc, SEN = sen, SPE = spe, TN = tn, FP = fp, FN = fn, TP = tp)
}

# --------------------------------------------------
# 3. 辅助函数：绘制单个混淆矩阵
# --------------------------------------------------
plot_cm <- function(cm, title, dataset_name, metrics) {
  
  # 转换为长格式
  cm_df <- as.data.frame(cm)
  cm_df$True_Label <- rownames(cm_df)
  cm_melt <- melt(cm_df, id.vars = "True_Label", variable.name = "Predicted_Label", value.name = "Count")
  
  # 计算百分比
  row_sum <- aggregate(Count ~ True_Label, data = cm_melt, FUN = sum)
  colnames(row_sum) <- c("True_Label", "RowTotal")
  cm_melt <- merge(cm_melt, row_sum, by = "True_Label")
  cm_melt$Percent <- cm_melt$Count / cm_melt$RowTotal * 100
  
  # 创建标注文本 (数值\n百分比)
  cm_melt$Label <- sprintf("%d\n(%.1f%%)", cm_melt$Count, cm_melt$Percent)
  
  # 确定文字颜色（深色背景用白色，浅色背景用黑色）
  cm_melt$TextColor <- ifelse(cm_melt$Count > max(cm) * 0.5, "white", "black")
  
  # 因子顺序
  cm_melt$True_Label <- factor(cm_melt$True_Label, levels = c("AML (1)", "Normal (0)"))
  cm_melt$Predicted_Label <- factor(cm_melt$Predicted_Label, levels = c("Normal (0)", "AML (1)"))
  
  # 构建指标文本
  metrics_text <- sprintf("ACC: %.1f%% | SEN: %.1f%% | SPE: %.1f%%", 
                          metrics$ACC, metrics$SEN, metrics$SPE)
  
  # 绘制
  p <- ggplot(cm_melt, aes(x = Predicted_Label, y = True_Label, fill = Count)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = Label, color = TextColor), 
              size = 5, fontface = "bold", lineheight = 0.9) +
    scale_fill_gradient(low = "#E8F4F8", high = "#08306B", 
                        limits = c(0, 80), name = "Count") +
    scale_color_identity() +
    labs(
      title = paste0(dataset_name, "\n", title),
      x = "Predicted Label",
      y = "True Label"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 10)),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 30, 10)
    ) 
    # 添加指标注释
    # annotate("text", x = 1.5, y = 0.3, 
    #          label = metrics_text, 
    #          fontface = "bold.italic", size = 4, color = "#333333")
  
  return(p)
}

# --------------------------------------------------
# 4. 计算各数据集指标
# --------------------------------------------------
metrics1 <- calc_metrics(cm1)
metrics2 <- calc_metrics(cm2)
metrics3 <- calc_metrics(cm3)

# --------------------------------------------------
# 5. 生成三个混淆矩阵图
# --------------------------------------------------
p1 <- plot_cm(cm1, "(GSE165656_GTEx)", "External_1", metrics1)
p2 <- plot_cm(cm2, "(GSE103424_GTEx)", "External_2", metrics2)
p3 <- plot_cm(cm3, "(GSE137851_GTEx)", "External_3", metrics3)

# --------------------------------------------------
# 6. 组合图形并保存
# --------------------------------------------------
# 同时保存为 PDF（矢量图，投稿推荐）
pdf("./21_SuperLearner/Superlearner绘图/5_三个外部验证集混淆矩阵.pdf", 
    width = 15, height = 5)
combined <- grid.arrange(
  p1, p2, p3, 
  ncol = 3
)
dev.off()



# --------------------------------------------------
# 7. 打印各数据集详细指标
# --------------------------------------------------
cat("【各外部验证集性能指标】\n")
cat("──────────────────────────────────────────────\n")
for(i in 1:3) {
  cm <- list(cm1, cm2, cm3)[[i]]
  m <- list(metrics1, metrics2, metrics3)[[i]]
  name <- c("External_1", "External_2", "External_3")[i]
  
  cat(sprintf("\n%s:\n", name))
  cat(sprintf("  Confusion Matrix: TN=%d, FP=%d, FN=%d, TP=%d\n", 
              m$TN, m$FP, m$FN, m$TP))
  cat(sprintf("  ACC: %.2f%% | SEN: %.2f%% | SPE: %.2f%%\n", 
              m$ACC, m$SEN, m$SPE))
}

