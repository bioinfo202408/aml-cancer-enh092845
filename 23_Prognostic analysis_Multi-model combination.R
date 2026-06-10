# 0. 环境初始化 -----------------------------------------------------------
rm(list = ls())
gc()

setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码//")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

dir.create("./23_预后分析_多模型组合/", showWarnings = FALSE, recursive = TRUE)


# TCGA训练集；GSE165656 测试集


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
## 0. Packages
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
## 1. Paths
############################################################

tcga_file <- "../TCGA_data/TCGA_LAML_879_feature_expr_with_ostime_new_136_samples.csv"
gse_file  <- "../TCGA_data/GSE165656_LAML_879_feature_expr_with_ostime_new.csv"

univ_cox_file <- "../Outdata/7.Prognostic analyse/2_univariate_cox_significant_genes.csv"

outdata_dir <- "./23_预后分析_多模型组合/outdata"
outplot_dir <- "./23_预后分析_多模型组合//outplot"

dir.create(outdata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outplot_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 2. Read and clean data
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
  stop("TCGA 数据必须包含 ID, OS.time, OS 三列。")
}

if (!all(c("ID", "OS.time", "OS") %in% colnames(gse))) {
  stop("GSE165656 数据必须包含 ID, OS.time, OS 三列。")
}

tcga$OS.time <- as.numeric(as.character(tcga$OS.time))
tcga$OS      <- as.numeric(as.character(tcga$OS))

gse$OS.time <- as.numeric(as.character(gse$OS.time))
gse$OS      <- as.numeric(as.character(gse$OS))

## 注意：
## 如果你的 GSE165656 文件里的 OS.time 已经是 day，请把下一行注释掉。
## 如果 GSE165656 文件里的 OS.time 是 week，则保留下一行。
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
## 3. Common genes and candidate genes
############################################################

cat("\n==============================\n")
cat("Step 2: Preparing genes\n")
cat("==============================\n")

clin_cols <- c("ID", "OS.time", "OS")

tcga_genes <- setdiff(colnames(tcga), clin_cols)
gse_genes  <- setdiff(colnames(gse), clin_cols)

common_genes <- intersect(tcga_genes, gse_genes)

cat("TCGA gene number:", length(tcga_genes), "\n")
cat("GSE gene number:", length(gse_genes), "\n")
cat("Common gene number:", length(common_genes), "\n")

## 如果存在单因素 Cox 显著基因文件，则使用该文件进一步筛选；
## 如果不存在，则使用所有共同基因。
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
    stop("单因素 Cox 文件中未找到 gene 或 gene_name 列。")
  }
  
  model_genes <- intersect(common_genes, univ_genes)
  
  cat("Using univariate Cox significant genes.\n")
  cat("Candidate gene number:", length(model_genes), "\n")
  
} else {
  
  model_genes <- common_genes
  
  cat("Univariate Cox file not found. Using all common genes.\n")
  cat("Candidate gene number:", length(model_genes), "\n")
}

if (length(model_genes) < 2) {
  stop("候选建模基因少于 2 个，请检查基因名。")
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
  cat("Remove genes with NA:", length(na_genes), "\n")
  model_genes <- setdiff(model_genes, na_genes)
}

if (length(model_genes) < 2) {
  stop("去除 NA 基因后，候选基因少于 2 个。")
}

tcga <- tcga[, c(clin_cols, model_genes)]
gse  <- gse[, c(clin_cols, model_genes)]

############################################################
## 4. Standardization using training parameters
############################################################

cat("\n==============================\n")
cat("Step 3: Standardization\n")
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

cat("Final modeling gene number:", length(pre_var), "\n")

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
## 5. Utility functions
############################################################

cat("\n==============================\n")
cat("Step 4: Defining functions\n")
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
          cat("Warning in:", model_name, "\n")
          cat("Warning:", conditionMessage(w), "\n")
          invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) {
      cat("Failed:", model_name, "\n")
      cat("Error:", conditionMessage(e), "\n")
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
## 6. Feature selection functions
############################################################

############################################################
## 修改 1：设置筛选后最多保留多少个基因
############################################################

max_feature_genes <- 30

select_all <- function(est_dd, pre_var) {
  pre_var
}

############################################################
## 修改 2：RSF 筛基因
## 只保留 VIMP 排名前 max_feature_genes 的基因
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
    stop("RSF importance 结果少于 2 个基因。")
  }
  
  ## 优先保留 VIMP > 0 的基因
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
## 修改 3：LASSO 筛基因
## 只保留绝对系数排名前 max_feature_genes 的基因
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
## 修改 4：CoxBoost 筛基因
## 只保留绝对系数排名前 max_feature_genes 的基因
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
        "StepCox is skipped because gene number > ",
        max_genes_for_step,
        ". Use prior feature selection before StepCox."
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
## 7. Modeling functions
############################################################

############################################################
## 修改 1：普通 Cox 建模函数
## 目的：
## 1. 不允许 All 137 genes 直接 Cox
## 2. 只允许筛选后少量基因进入 Cox
############################################################

fit_cox_model <- function(est_dd, val_dd_list, rid,
                          max_genes_for_cox = 50,
                          epv_cutoff = 2) {
  
  event_num <- sum(est_dd$OS == 1, na.rm = TRUE)
  
  if (length(rid) < 2) {
    stop("普通 Cox 跳过：基因数少于 2。")
  }
  
  if (length(rid) > max_genes_for_cox) {
    stop(
      paste0(
        "普通 Cox 跳过：筛选后基因数 = ", length(rid),
        "，超过 max_genes_for_cox = ", max_genes_for_cox
      )
    )
  }
  
  if (length(rid) >= floor(event_num / epv_cutoff)) {
    stop(
      paste0(
        "普通 Cox 跳过：筛选后基因数 = ", length(rid),
        "，事件数 = ", event_num,
        "，变量数仍然偏多。"
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
    stop("普通 Cox 拟合后存在 NA 系数，跳过。")
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
## 修改 2：StepCox 建模函数
## 目的：
## 1. 不允许 All 137 genes 直接 StepCox
## 2. 只允许筛选后少量基因进入 StepCox
############################################################

fit_stepcox_model <- function(est_dd, val_dd_list, rid,
                              direction = "backward",
                              max_genes_for_stepcox = 50,
                              epv_cutoff = 2) {
  
  event_num <- sum(est_dd$OS == 1, na.rm = TRUE)
  
  if (length(rid) < 2) {
    stop("StepCox 跳过：基因数少于 2。")
  }
  
  if (length(rid) > max_genes_for_stepcox) {
    stop(
      paste0(
        "StepCox 跳过：筛选后基因数 = ", length(rid),
        "，超过 max_genes_for_stepcox = ", max_genes_for_stepcox
      )
    )
  }
  
  if (length(rid) >= floor(event_num / epv_cutoff)) {
    stop(
      paste0(
        "StepCox 跳过：筛选后基因数 = ", length(rid),
        "，事件数 = ", event_num,
        "，变量数仍然偏多。"
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
    stop("StepCox 初始 Cox 模型存在 NA 系数，跳过。")
  }
  
  fit <- step(
    fit0,
    direction = direction,
    trace = 0
  )
  
  final_genes <- names(coef(fit))
  
  if (length(final_genes) < 1) {
    stop("StepCox final model has no gene.")
  }
  
  if (any(is.na(coef(fit)))) {
    stop("StepCox 最终模型存在 NA 系数，跳过。")
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
## 8. Define feature selectors and model builders
############################################################

cat("\n==============================\n")
cat("Step 5: Defining model combinations\n")
cat("==============================\n")

############################################################
## 修改 5：所有筛选方法最多保留 max_feature_genes 个基因
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

## 如果候选基因数量不太多，也可以启用 StepCox 作为特征筛选。
## 这里设置不超过 80 个基因才运行。
if (length(pre_var) <= 80) {
  feature_selectors$StepCox_backward <- function() {
    select_stepcox_backward(est_dd, pre_var, max_genes_for_step = 80)
  }
}

model_builders <- list(
  ############################################################
  ## 修改 3：Cox / StepCox 只允许筛选后少量基因进入
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

## Add Elastic Net models
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
## 9. Run all combinations
############################################################

############################################################
## 修改 4：重新运行所有模型组合
## 核心逻辑：
## 1. All + Cox 直接跳过
## 2. All + StepCox 直接跳过
## 3. 只有 RSF/Lasso/CoxBoost 等筛基因后，才允许 Cox/StepCox
############################################################

cat("\n==============================\n")
cat("Step 6: Running all model combinations\n")
cat("==============================\n")

result <- data.frame()
model_store <- list()
gene_store <- list()

max_genes_for_cox_models <- 50
event_num <- sum(est_dd$OS == 1, na.rm = TRUE)

for (fs_name in names(feature_selectors)) {
  
  cat("\n###################################################\n")
  cat("Feature selection:", fs_name, "\n")
  cat("###################################################\n")
  
  rid <- safe_model(
    paste0("FeatureSelect_", fs_name),
    feature_selectors[[fs_name]]()
  )
  
  if (is.null(rid)) next
  
  rid <- intersect(rid, pre_var)
  rid <- unique(rid)
  
  if (length(rid) < 2) {
    cat("Selected gene number < 2. Skip:", fs_name, "\n")
    next
  }
  
  cat("Selected gene number:", length(rid), "\n")
  
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
    ## 修改 4.1：
    ## All 不算真正的特征筛选。
    ## 所以 All + Cox 和 All + StepCox 直接跳过。
    ################################################
    
    if (fs_name == "All" &&
        model_name %in% c("Cox", "StepCox_backward")) {
      
      cat(
        "Skip:", full_model_name,
        "| Reason: Cox/StepCox must be performed after feature selection, not All genes.\n"
      )
      
      next
    }
    
    ################################################
    ## 修改 4.2：
    ## 即使经过筛选，如果基因数仍然太多，
    ## Cox 和 StepCox 也跳过。
    ################################################
    
    if (model_name %in% c("Cox", "StepCox_backward") &&
        length(rid) > max_genes_for_cox_models) {
      
      cat(
        "Skip:", full_model_name,
        "| selected gene number =", length(rid),
        "> max_genes_for_cox_models =", max_genes_for_cox_models,
        "\n"
      )
      
      next
    }
    
    ################################################
    ## 修改 4.3：
    ## 如果变量数相对事件数仍然偏多，
    ## Cox 和 StepCox 也跳过。
    ################################################
    
    if (model_name %in% c("Cox", "StepCox_backward") &&
        length(rid) >= floor(event_num / 2)) {
      
      cat(
        "Skip:", full_model_name,
        "| selected gene number =", length(rid),
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
    
    cat("Finished:", full_model_name, "\n")
  }
}

if (nrow(result) == 0) {
  stop("所有模型均运行失败，请检查数据和参数。")
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
  stop("所有模型均运行失败，请检查数据和包。")
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
## 10. Summarize C-index
############################################################

cat("\n==============================\n")
cat("Step 7: Summarizing C-index\n")
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

cat("\nTop 20 models by GSE165656 C-index:\n")
print(head(cindex_wide, 20))

############################################################
## 11. C-index heatmap
############################################################

cat("\n==============================\n")
cat("Step 8: Plotting C-index heatmap\n")
cat("==============================\n")

############################################################
## 11. C-index heatmap
############################################################

cat("\n==============================\n")
cat("Step 8: Plotting C-index heatmap\n")
cat("==============================\n")

cindex_plot <- cindex_wide[
  !is.na(cindex_wide$TCGA) & !is.na(cindex_wide$GSE165656),
]

############################################################
## 11. C-index heatmap
## 修改：
## 1. 新增 Average_Cindex = mean(TCGA, GSE165656)
## 2. 按 Average_Cindex 从高到低排序
## 3. 去除相同模型组合，如 RSF + RSF、Lasso + Lasso 等
## 4. 使用更适合 SCI 的配色
############################################################

cat("\n==============================\n")
cat("Step 8: Plotting C-index heatmap\n")
cat("==============================\n")

cindex_plot <- cindex_wide[
  !is.na(cindex_wide$TCGA) & !is.na(cindex_wide$GSE165656),
]
cindex_plot <- subset(
  cindex_plot,
  GSE165656 >= 0.50
)
## 去除相同模型组合，例如 RSF + RSF, Lasso + Lasso, CoxBoost + CoxBoost
same_model_rows <- grepl("^(.+) \\+ \\1$", cindex_plot$Model)

if (any(same_model_rows)) {
  cat("Remove same-model combinations:\n")
  print(cindex_plot$Model[same_model_rows])
}

cindex_plot <- cindex_plot[!same_model_rows, ]

## 新增平均 C-index
cindex_plot$Average_Cindex <- rowMeans(
  cindex_plot[, c("TCGA", "GSE165656"), drop = FALSE],
  na.rm = TRUE
)

## 按平均 C-index 从高到低排序
cindex_plot <- cindex_plot[
  order(cindex_plot$Average_Cindex, decreasing = TRUE),
]

## 保存用于画图的数据
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
  
  ## 列名美化
  colnames(dt) <- c("TCGA", "GSE165656", "Average")
  
  ## SCI风格配色：蓝-白-红，低值蓝，高值红
  heat_cols <- colorRampPalette(
    c("#3953A4", "#FFFFFF", "#ED2123")
  )(100)
  
  ## 根据模型数量自动调整高度
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
    main = "C-index of Survival Models"
  )
  
  dev.off()
}

############################################################
## 12. Select best model
############################################################

cat("\n==============================\n")
cat("Step 9: Selecting best model\n")
cat("==============================\n")

valid_rank <- cindex_wide[!is.na(cindex_wide$GSE165656), ]

if (nrow(valid_rank) == 0) {
  stop("没有任何模型得到 GSE165656 C-index。")
}

best_model_name <- valid_rank$Model[1]

cat("Best model selected by GSE165656 C-index:", best_model_name, "\n")

writeLines(
  best_model_name,
  file.path(outdata_dir, "2_best_model_name.txt")
)

best_obj <- model_store[[best_model_name]]

if (is.null(best_obj)) {
  stop("无法从 model_store 中找到最佳模型对象。")
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
## 13. Save best model coefficients if available
############################################################

cat("\n==============================\n")
cat("Step 10: Saving coefficients if available\n")
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
  
  cat("Coefficient file saved.\n")
  
} else {
  
  cat("Best model has no simple linear coefficients or coefficients unavailable.\n")
}

############################################################
## 14. Best model riskScore and risk group
############################################################

cat("\n==============================\n")
cat("Step 11: Risk score and risk groups\n")
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

## cutoff 只能来自训练集
risk_cutoff <- median(train_df$riskScore, na.rm = TRUE)

train_df$risk_group <- ifelse(train_df$riskScore >= risk_cutoff, "High", "Low")
valid_df$risk_group <- ifelse(valid_df$riskScore >= risk_cutoff, "High", "Low")

train_df$risk_group <- factor(train_df$risk_group, levels = c("Low", "High"))
valid_df$risk_group <- factor(valid_df$risk_group, levels = c("Low", "High"))

cat("Risk cutoff from TCGA median:", risk_cutoff, "\n")

cat("\nTCGA risk groups:\n")
print(table(train_df$risk_group))

cat("\nGSE165656 risk groups:\n")
print(table(valid_df$risk_group))

write_tsv(
  train_df,
  file.path(outdata_dir, "3_best_model_riskScore_TCGA.txt")
)

write_tsv(
  valid_df,
  file.path(outdata_dir, "3_best_model_riskScore_GSE165656.txt")
)

writeLines(
  as.character(risk_cutoff),
  file.path(outdata_dir, "3_best_model_cutoff_from_TCGA_median.txt")
)

############################################################
## 15. KM plotting function
############################################################

plot_km_best <- function(dat, title, outfile) {
  
  dat <- dat %>%
    mutate(
      OS.month = OS.time / 30,
      risk_group = factor(risk_group, levels = c("Low", "High"))
    )
  
  if (length(unique(dat$risk_group)) < 2) {
    cat("Only one risk group in ", title, ". Skip KM plot.\n")
    return(NULL)
  }
  
  fit <- survfit(Surv(OS.month, OS) ~ risk_group, data = dat)
  
  cox_fit <- tryCatch(
    coxph(Surv(OS.month, OS) ~ risk_group, data = dat),
    error = function(e) NULL
  )
  
  hr_label <- ""
  
  if (!is.null(cox_fit)) {
    
    hr <- round(exp(coef(cox_fit)), 2)
    ci <- round(exp(confint(cox_fit)), 2)
    pcox <- signif(summary(cox_fit)$waldtest["pvalue"], 3)
    
    hr_label <- paste0(
      "HR = ", hr,
      " (", ci[1], "-", ci[2], ")\n",
      "Cox p = ", pcox
    )
  }
  
  pub_theme <- theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", linewidth = 0.8),
      plot.title = element_text(hjust = 0.5, size = 13),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10, colour = "black"),
      legend.text = element_text(size = 10),
      legend.position = c(0.82, 0.85),
      legend.background = element_blank()
    )
  
  p <- ggsurvplot(
    fit,
    data = dat,
    title = title,
    xlab = "Time months",
    ylab = "Survival probability",
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    pval = TRUE,
    pval.size = 4.5,
    conf.int = FALSE,
    palette = c("#2E86AB", "#E63946"),
    risk.table.height = 0.22,
    legend.title = "",
    legend.labs = c("Low risk", "High risk"),
    ggtheme = pub_theme,
    surv.median.line = "hv"
  )
  
  if (hr_label != "") {
    p$plot <- p$plot +
      annotate(
        "text",
        x = max(dat$OS.month, na.rm = TRUE) * 0.55,
        y = 0.20,
        label = hr_label,
        size = 4
      )
  }
  
  pdf(outfile, width = 8, height = 6)
  print(p)
  dev.off()
  
  return(p)
}

############################################################
## 16. Plot KM curves
############################################################

cat("\n==============================\n")
cat("Step 12: Plotting KM curves\n")
cat("==============================\n")

plot_km_best(
  train_df,
  paste0("Overall survival in TCGA training cohort\n", best_model_name),
  file.path(outplot_dir, "2_best_model_TCGA_KM.pdf")
)

plot_km_best(
  valid_df,
  paste0("Overall survival in GSE165656 validation cohort\n", best_model_name),
  file.path(outplot_dir, "2_best_model_GSE165656_KM.pdf")
)

############################################################
## 17. Save final summary
############################################################

cat("\n==============================\n")
cat("Step 13: Final summary\n")
cat("==============================\n")

best_cindex <- cindex_wide[cindex_wide$Model == best_model_name, ]

summary_df <- data.frame(
  Best_Model = best_model_name,
  TCGA_Cindex = best_cindex$TCGA,
  GSE165656_Cindex = best_cindex$GSE165656,
  Mean_Cindex = best_cindex$Mean_Cindex,
  Cutoff_from_TCGA_median = risk_cutoff,
  Gene_number = length(best_obj$genes),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  file.path(outdata_dir, "4_final_best_model_summary.csv"),
  row.names = FALSE
)

print(summary_df)

cat("\nAll analyses completed successfully.\n")
cat("Output data directory:", outdata_dir, "\n")
cat("Output plot directory:", outplot_dir, "\n")
