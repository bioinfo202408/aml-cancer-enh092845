# 画图参考：https://mp.weixin.qq.com/s/zI7VN5Q_yYDkzDE_KGhTXg
# 需要跑完GO/KEGG富集分析的结果


# 清空环境并加载必要的包
rm(list = ls())
gc() #清理内存中不再使用的对象

# 设置工作路径
setwd("/home/weili/Project/AML/human/AML_combined_analyse/0.画图代码/")
cat("【初始化】当前工作路径：", getwd(), "\n\n")

#查看临时目录
Sys.getenv("TMPDIR") 

dir.create("./9_GO_KEGG富集分析/")



# ==============================================
# 1. 加载包
# ==============================================
library(ggplot2)
library(dplyr)
library(patchwork)


# ==============================================
# ######KEGG富集分析####
# ==============================================
# ==============================================
# 2. 构建你的 KEGG 数据框（直接复制运行即可）
# ==============================================
kegg_data <- read.csv("./9_GO_KEGG富集分析/KEGG_Enh092845.csv")
kegg_data <- read.csv("./9_GO_KEGG富集分析/预后相关ncRNA相关mRNA_enriched_KEGG_terms.csv")


# # 保留感兴趣的通路
# keypathway <- c("Oxidative phosphorylation","Lysosome biogenesis","Phagosome","Adherens junction",
#                 "Insulin resistance","Thermogenesis","Salmonella infection",
#                 "Human cytomegalovirus infection","Vibrio cholerae infection","Parkinson disease")
# kegg_data <- kegg_data[kegg_data$Description %in% keypathway,]


# 取 Top10（按 Count 降序）
kegg_data <- kegg_data %>%
  arrange(desc(Count)) %>%
  head(10)


# 计算 -log10(p.adjust)
kegg_data <- kegg_data %>%
  mutate(logP = -log10(p.adjust))

# 因子化，让顺序从上到下和数据一致
kegg_data$Description <- factor(kegg_data$Description, levels = rev(kegg_data$Description))

# ==============================================
# 3. 配色（和你 GO 图保持一致）
# ==============================================
Color <- rev(c("#D65656", "#5FAB5F", "#DDB370", "#E1A99A", "#CFE6A1",
               "#72C3E3", "#D43B63", "#D796C3", "#9B3683", "#294B2E"))

# ==============================================
# 4. 左侧：棒棒糖图（Count）
# ==============================================
p1 <- ggplot(kegg_data, aes(x = Count, y = subcategory)) +
  geom_segment(aes(y = subcategory, yend = subcategory, x = 0, xend = Count),
               linewidth = 0.8, color = "grey80") +
  geom_point(aes(color = subcategory), size = 6) +
  scale_color_manual(values = Color) +
  scale_x_continuous(limits = c(0, max(kegg_data$Count) + 5), expand = c(0, 0)) +
  labs(title = "Cell marker numbers", x = NULL, y = NULL) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )


# ==============================================
# 5. 左侧：气泡图
# ==============================================
kegg_data$size <- as.numeric(sapply(strsplit(kegg_data$GeneRatio, "/"), "[", 1))
p3 <- ggplot(kegg_data, aes(x = 0, y = subcategory)) +  
  geom_point(aes(color = subcategory, size = size)) +  
  scale_color_manual(values = Color, guide = "none") +  
  scale_size_continuous(name = "Marker Num.", breaks = seq(0,100, 20), limits = c(0,100)) +  
  scale_x_discrete(expand = c(0,0)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_bw(base_size = 18) +  
  theme(    
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),  # 去掉加粗
    panel.border = element_blank(),    
    axis.text.x = element_text(color = "black", face = "plain"),       # 去掉加粗
    axis.text.y = element_text(hjust = 1,color = "black", face = "plain"), # 去掉加粗
    axis.ticks = element_blank(),    
    panel.grid = element_blank(),    
    legend.position = "right",
    legend.text = element_text(face = "plain"),  # 图例不加粗
    legend.title = element_text(face = "plain")  # 图例标题不加粗
  )


# 如果上面画的有问题，比如很多通路用一个subcategory，就用下面的画
# # ==============================================
# # 4. 左侧：棒棒糖图（每个条目一根棒！不合并！）
# # ==============================================
# # 关键：给每一行创建一个唯一ID，防止自动合并
# kegg_data$id <- 1:nrow(kegg_data)
# 
# p1 <- ggplot(kegg_data, aes(x = Count, y = id)) +
#   # 每一行都画一根棒棒
#   geom_segment(aes(y = id, yend = id, x = 0, xend = Count),
#                linewidth = 0.8, color = "grey80") +
#   # 每一行都画一个点，颜色按 subcategory
#   geom_point(aes(color = subcategory), size = 6) +
#   
#   # 重点：Y轴标签强制显示为 subcategory
#   scale_y_continuous(breaks = 1:nrow(kegg_data),
#                      labels = kegg_data$subcategory) +
#   
#   scale_color_manual(values = Color) +
#   scale_x_continuous(limits = c(0, max(kegg_data$Count) + 5), expand = c(0, 0)) +
#   labs(title = "Cell marker numbers", x = NULL, y = NULL) +
#   theme_bw(base_size = 18) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
#     axis.text.x = element_text(color = "black", face = "plain"),
#     axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
#     panel.grid = element_blank(),
#     legend.position = "none",
#     panel.border = element_rect(size = 0.5, fill = NA)
#   )
# # ==============================================
# # 5. 左侧：气泡图（每行独立，不合并，Y轴subcategory）
# # ==============================================
# kegg_data$size <- as.numeric(sapply(strsplit(kegg_data$GeneRatio, "/"), "[", 1))
# 
# # 关键：用唯一id做Y轴，强制每行显示
# kegg_data$id <- 1:nrow(kegg_data)
# 
# p3 <- ggplot(kegg_data, aes(x = 0, y = id)) +  
#   geom_point(aes(color = subcategory, size = size)) +  
#   scale_color_manual(values = Color, guide = "none") +  
#   scale_size_continuous(name = "Marker Num.", breaks = seq(0,100, 20), limits = c(0,100)) +  
#   scale_y_continuous(breaks = 1:nrow(kegg_data), labels = kegg_data$subcategory) + # 强制显示subcategory
#   scale_x_discrete(expand = c(0,0)) +  
#   labs(title = "", x = NULL, y = NULL) +  
#   theme_bw(base_size = 18) +  
#   theme(    
#     plot.title = element_text(hjust = 0, size = 18, face = "plain"),
#     panel.border = element_blank(),    
#     axis.text.x = element_blank(),    # 隐藏x轴文字
#     axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
#     axis.ticks = element_blank(),    
#     panel.grid = element_blank(),    
#     legend.position = "right",
#     legend.text = element_text(face = "plain"),
#     legend.title = element_text(face = "plain")
#   )




# ==============================================
# 6. 右侧：条形图（-log10(p.adjust)）
# ==============================================
p2 <- ggplot(kegg_data, aes(y = Description)) +
  geom_bar(aes(x = logP, fill = Description),
           stat = "identity", width = 0.5, color = "transparent", alpha = 0.7) +
  geom_text(aes(x = 0.1, label = Description), hjust = 0, size = 5.5, face = "plain") +
  labs(title = "KEGG enrichment item", x = "-log10(p.adjust)", y = NULL) +
  scale_fill_manual(values = Color) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )


# ==============================================
# 7. 合并并保存 PDF
# ==============================================
p <- p1 + p2 + plot_layout(widths = c(1, 1.8))
# 左侧棒棒图+右侧条形图
ggsave("./9_GO_KEGG富集分析/KEGG_棒棒_条形图_bar.pdf",
       plot = p,
       width = 10,
       height = 6,
       device = "pdf")

# 左侧气泡图+右侧条形图
p <- p3 + p2 + 
  plot_layout(widths = c(0.12, 1.2), guides = "collect")
# 左侧棒棒图+右侧条形图
ggsave("./9_GO_KEGG富集分析/KEGG_气泡_条形_bar.pdf",
       plot = p,
       width = 10,
       height = 6,
       device = "pdf")









# ==============================================
# ######GO富集分析######
# ==============================================
kegg_data <- read.csv("./9_GO_KEGG富集分析/预后相关ncRNA相关mRNA_enriched_GO_terms.csv")

# 计算 GeneRatio ：把 "33/739" 变成 0.04465
kegg_data$GeneRatio_value <- sapply(kegg_data$GeneRatio, function(x) eval(parse(text = x)))


# 取 Top10（按 GeneRatio_value 降序）
kegg_data <- kegg_data %>%
  arrange(desc(GeneRatio_value)) %>%
  head(10)


# 计算 -log10(p.adjust)
kegg_data <- kegg_data %>%
  mutate(logP = -log10(p.adjust))

# 因子化，让顺序从上到下和数据一致
kegg_data$Description <- factor(kegg_data$Description, levels = rev(kegg_data$Description))

# ==============================================
# 3. 配色（和你 GO 图保持一致）
# ==============================================
Color <- rev(c("#D65656", "#5FAB5F", "#DDB370", "#E1A99A", "#CFE6A1",
               "#72C3E3", "#D43B63", "#D796C3", "#9B3683", "#294B2E"))

# # ==============================================
# # 4. 左侧：棒棒糖图（Count）
# # ==============================================
# p1 <- ggplot(kegg_data, aes(x = Count, y = subcategory)) +
#   geom_segment(aes(y = subcategory, yend = subcategory, x = 0, xend = Count),
#                linewidth = 0.8, color = "grey80") +
#   geom_point(aes(color = subcategory), size = 6) +
#   scale_color_manual(values = Color) +
#   scale_x_continuous(limits = c(0, max(kegg_data$Count) + 5), expand = c(0, 0)) +
#   labs(title = "Cell marker numbers", x = NULL, y = NULL) +
#   theme_bw(base_size = 18) +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
#     axis.text.x = element_text(color = "black", face = "plain"),
#     axis.text.y = element_text(hjust = 1, color = "black", face = "plain"),
#     panel.grid = element_blank(),
#     legend.position = "none",
#     panel.border = element_rect(size = 0.5, fill = NA)
#   )

# ==============================================
# 5. 左侧：气泡图
# ==============================================
p3 <- ggplot(kegg_data, aes(x = 0, y = Description)) +  
  geom_point(aes(color = Description, size = Count)) +  
  scale_color_manual(values = Color, guide = "none") +  
  scale_size_continuous(name = "Marker Num.", breaks = seq(0,100, 20), limits = c(0,100)) +  
  scale_x_discrete(expand = c(0,0)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_bw(base_size = 18) +  
  theme(    
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),  # 去掉加粗
    panel.border = element_blank(),    
    # axis.text.x = element_text(color = "black", face = "plain"),       # 去掉加粗
    # axis.text.y = element_text(hjust = 1,color = "black", face = "plain"), # 去掉加粗
    axis.text.x = element_blank(),  # 也隐藏x轴文字
    axis.text.y = element_blank(),  # ✅ 隐藏y轴所有文字（你要的效果）
    axis.ticks = element_blank(),    
    panel.grid = element_blank(),    
    legend.position = "right",
    legend.text = element_text(face = "plain"),  # 图例不加粗
    legend.title = element_text(face = "plain")  # 图例标题不加粗
  )

# ==============================================
# 6. 右侧：条形图（-log10(p.adjust)）
# ==============================================
p2 <- ggplot(kegg_data, aes(y = Description)) +
  geom_bar(aes(x = logP, fill = Description),
           stat = "identity", width = 0.5, color = "transparent", alpha = 0.7) +
  geom_text(aes(x = 0.1, label = Description), hjust = 0, size = 5.5, face = "plain") +
  labs(title = "GO enrichment item", x = "-log10(p.adjust)", y = NULL) +
  scale_fill_manual(values = Color) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )


# ==============================================
# 7. 合并并保存 PDF
# ==============================================
# p <- p1 + p2 + plot_layout(widths = c(1, 1.8))
# # 左侧棒棒图+右侧条形图
# ggsave("./9_GO_KEGG富集分析/KEGG_棒棒_条形图_bar.pdf",
#        plot = p,
#        width = 10,
#        height = 6,
#        device = "pdf")

# 左侧气泡图+右侧条形图
p <- p3 + p2 + 
  plot_layout(widths = c(0.12, 1.2), guides = "collect")
# 左侧棒棒图+右侧条形图
ggsave("./9_GO_KEGG富集分析/GO_气泡_条形_bar.pdf",
       plot = p,
       width = 8,
       height = 6,
       device = "pdf")




# ==============================================
# ######IRF1 KO数据 扰动基因-GO富集分析######
# ==============================================
library(openxlsx)
library(tidyverse)
library(patchwork)

# 读取数据
kegg_data <- read.xlsx("./9_GO_KEGG富集分析/IRF1_KO_virtual_ko_go_result_ENRICH.xlsx")
kegg_data <- kegg_data[-1,]

# 从 "分子/分母" 格式中提取分子，作为 Count 列
kegg_data$Count <- as.integer(sapply(strsplit(kegg_data$InTerm_InList, "/"), function(x) x[1]))

# 计算 GeneRatio
# 安全计算 GeneRatio，自动跳过 13/- 这种脏数据
kegg_data$GeneRatio_value <- sapply(kegg_data$InTerm_InList, function(x) {
  # 只保留 数字/数字 这种格式
  if(grepl("^\\d+/\\d+$", x)) {
    eval(parse(text = x))
  } else {
    NA  # 脏数据直接设为NA，不报错
  }
})

# 去掉脏数据行（确保画图干净）
kegg_data <- kegg_data[!is.na(kegg_data$GeneRatio_value), ]


# 按 LogP 从大到小排序，取 TOP 10
# ======================================================
kegg_data <- kegg_data %>%
  # arrange(desc(LogP)) %>%   # 按显著性排序：从大到小
  arrange(LogP) %>%   # 从小到大（去掉 desc 就是升序）
  head(10)                  # 取前10

# 计算 logP（和你原来一致）
kegg_data <- kegg_data %>%
  mutate(logP = -LogP)

# 因子化，保证画图从上到下是 Top1 → Top10
kegg_data$Description <- factor(kegg_data$Description, levels = rev(kegg_data$Description))

# ==============================================
# 配色
# ==============================================
Color <- rev(c("#D65656", "#5FAB5F", "#DDB370", "#E1A99A", "#CFE6A1",
               "#72C3E3", "#D43B63", "#D796C3", "#9B3683", "#294B2E"))

# ==============================================
# 左侧：气泡图
# ==============================================
p3 <- ggplot(kegg_data, aes(x = 0, y = Term)) +    # 这里修改左边气泡要显示什么内容
  geom_point(aes(color = Term, size = Count)) +  
  scale_color_manual(values = Color, guide = "none") +  
  scale_size_continuous(name = "Marker Num.", breaks = seq(0,30, 10), limits = c(0,30)) +  
  scale_x_discrete(expand = c(0,0)) +  
  labs(title = "", x = NULL, y = NULL) +  
  theme_bw(base_size = 18) +  
  theme(    
    plot.title = element_text(hjust = 0, size = 18, face = "plain"),
    panel.border = element_blank(),    
    axis.text.x = element_blank(),
    # axis.text.y = element_blank(),  # 加了这个Y轴就不显示
    axis.ticks = element_blank(),    
    panel.grid = element_blank(),    
    legend.position = "right",
    legend.text = element_text(face = "plain"),
    legend.title = element_text(face = "plain")
  )

# ==============================================
# 右侧：条形图（-log10(p.adjust)）
# ==============================================
p2 <- ggplot(kegg_data, aes(y = Description)) +
  geom_bar(aes(x = logP, fill = Description),
           stat = "identity", width = 0.5, color = "transparent", alpha = 0.7) +
  geom_text(aes(x = 0.1, label = Description), hjust = 0, size = 5.5, face = "plain") +
  labs(title = "GO enrichment item", x = "-log10(p.adjust)", y = NULL) +
  scale_fill_manual(values = Color) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
    axis.text.x = element_text(color = "black", face = "plain"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    panel.border = element_rect(size = 0.5, fill = NA)
  )

# ==============================================
# 合并 + 出图
# ==============================================
p <- p3 + p2 + 
  plot_layout(widths = c(0.12, 1.2), guides = "collect")

ggsave("./9_GO_KEGG富集分析/IRF1_KO_GO_TOP10_LogP_气泡条形图.pdf",
       plot = p,
       width = 12,
       height = 7,
       device = "pdf",
       bg="white")

cat("✅ GO TOP10 按 LogP 排序图已生成！\n")
