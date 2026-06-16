# 参考：https://mp.weixin.qq.com/s?__biz=MzU5NjE2ODU0OQ==&mid=2247484972&idx=1&sn=dfa146cdecc9a1158a40e4c282937429&scene=21&poc_token=HHNeIWqjXu-ab4ymdTs5MrnF-8oNBt9tiiRw-uF4

rm(list = ls())
gc()

# 工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

# 创建输出目录
out_dir <- "./22_雷达图/"
dir.create(out_dir, showWarnings = F, recursive = T)

# 加载包
library(ggradar)
library(dplyr)
library(scales)
library(tibble)
library(plotthis)
library(ggplot2)
library(tidyr)

# 读入原始结果
df_raw <- read.csv("./21_SuperLearner/Outdata/0_最佳元模型选择_MetaLearner_Result.csv")
colnames(df_raw) <- c(
  "Meta_Learner", "Train_AUC", "Train_Acc", "Train_Sen", "Train_Spe", "Train_Pre", "Train_F1", "Train_BACC", "Train_Thres",
  "Test_AUC", "Test_Acc", "Test_Sen", "Test_Spe", "Test_Pre", "Test_F1", "Test_BACC", "Test_Thres",
  "Val_AUC", "Val_Acc", "Val_Sen", "Val_Spe", "Val_Pre", "Val_F1", "Val_BACC", "Val_Thres"
)

# 筛选绘图指标，剔除阈值列
radar_data <- df_raw %>%
  select(Meta_Learner, contains(c("Train_","Val_")), -contains("Thres"))
colnames(radar_data)[1] <- "group"

## 1、原始数据ggradar总图（保持原样输出不变）
p1 <- ggradar(
  plot.data = radar_data,
  base.size = 12,
  font.radar = "sans",
  values.radar = c("0%", "50%", "100%"),
  grid.min = 0, grid.mid = 0.5, grid.max = 1,
  plot.extent.x.sf = 1, plot.extent.y.sf = 1.1,
  label.centre.y = F,
  grid.line.width = 0.4,
  gridline.min.linetype = "longdash",
  gridline.mid.linetype = "longdash",
  gridline.max.linetype = "longdash",
  gridline.min.colour = "gray70",
  gridline.mid.colour = "#007A87",
  gridline.max.colour = "gray70",
  grid.label.size = 4,
  label.gridline.min = F,label.gridline.mid=F,label.gridline.max=F,
  axis.label.offset = 1.18,
  axis.label.size = 3.8,
  axis.line.colour = "gray60",
  group.line.width = 1.3,
  group.point.size = 3,
  group.colours = c("#E63946","#F77F00","#FCBF49","#06D6A0","#118AB2","#073B4C"),
  background.circle.colour = "#F5F5F4",
  background.circle.transparency = 0.2,
  legend.title = "Meta Learner",
  legend.text.size = 9,
  legend.position = "right",
  fill = TRUE, fill.alpha = 0.18,
  plot.title = "SuperLearner Meta-Learner Performance Radar"
)
ggsave("./22_雷达图/1_最佳元模型选择_SuperLearner_Nature_Radar_ggradar.pdf",p1,width=10,height=7,dpi=300)

## 2、分面蜘蛛图：log1p=ln(x+1)对数变换
radar_long <- pivot_longer(
  data = radar_data,
  cols = -group,
  names_to = "Index",
  values_to = "Score_raw"
)
colnames(radar_long) <- c("Meta_Learner","Index","Score_raw")

# 对数变换：y = log1p(x)=ln(x+1)
radar_long <- radar_long %>%
  mutate(
    Score_log = log1p(Score_raw), # 对数拉伸高分细微差异
    Lab_text = round(Score_raw,3)  # 标签保留原始真实数值
  )

# 固定配色
color_vec <- c(
  "method.NNLS"     = "#E63946",
  "method.NNLS2"    = "#F77F00",
  "method.NNloglik" = "#FCBF49",
  "method.CC_LS"    = "#06D6A0",
  "method.CC_nloglik"="#118AB2",
  "method.AUC"      = "#073B4C"
)
# 图片6色：灰、深绿、天蓝、橙、紫、浅粉
color_vec <- c(
  "method.NNLS"     = "#A0A0A4",   # 灰色
  "method.NNLS2"    = "#008000",   # 深绿色
  "method.NNloglik" = "#57A4FD",   # 亮蓝色
  "method.CC_LS"    = "#FF6000",   # 橙红色
  "method.CC_nloglik"="#AD07E3",   # 紫色
  "method.AUC"      = "#F7ABE8"    # 浅粉色
)



# 绘制基础蜘蛛图，径向使用对数变换后数据
spider_base <- SpiderPlot(
  data = radar_long,
  x = "Index",
  y = "Score_log",
  split_by = "Meta_Learner",
  group_by = "Meta_Learner",
  y_nbreaks = 5,
  fill = TRUE,
  alpha = 0.25,
  linewidth = 1.3,
  pt_size = 3,
  palcolor = color_vec,
  facet_ncol = 3,
  facet_scales = "fixed",
  legend.position = "none",
  title = "SuperLearner Performance (log1p transformed radius, raw value label)",
  aspect.ratio = 1,
  theme_args = list(plot.title = element_text(hjust=0.5,size=14))
)

# 顶点外侧添加原始数值标签，向外偏移0.03
final_spider <- spider_base +
  geom_text(
    aes(label = Lab_text, y = Score_log + 0.03),
    size = 2.8,
    check_overlap = TRUE
  )

# 导出PDF
ggsave(
  filename = paste0(out_dir,"2_log1p对数变换_分面异色_原始数值标注_Spider2.pdf"),
  plot = final_spider,
  width = 15,
  height = 8,
  dpi = 300
)

把中文都换成英文
