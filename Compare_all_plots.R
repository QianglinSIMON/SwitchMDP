library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

create_plot <- function(data,  y_var = "RMSE", title, show_legend = T) {
  
  #data<-data_hu_DIM_rmse0; show_legend = FALSE; title="RMSE:0-burn-in"
  # 将 "rho0" 和 "m" 转换为因子变量，并确保 rho0 的水平一致
  data <- data %>%
    mutate_at(c("rho0", "m"), as_factor) %>%
    mutate(rho0 = factor(rho0, levels = c("0.3", "0.5", "0.7", "0.9")))  # 强制一致的 rho0 水平
  
  # 移除包含 NA 的行（确保数据完整）
  data <- na.omit(data)
  
  # 定义一个标签映射
  f_names <- list(
    "0.3" = "rho == 0.3",
    "0.5" = "rho == 0.5",
    "0.7" = "rho == 0.7",
    "0.9" = "rho == 0.9"
  )
  
  # 自定义标签器函数
  f_labeller <- function(variable, value) {
    return(f_names[value])
  }
  
  # 创建 ggplot 图形
  p <- ggplot(data, aes(x = n, y = .data[[y_var]], color = m)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    facet_grid(. ~ rho0, 
               labeller = labeller(rho0 = f_labeller,
                                   .default = label_parsed),  # 使用自定义的 labeller
               scales = "free") +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),  # 标题加粗加大
      axis.text = element_text(size = 20, face = "bold"),  # 坐标轴文本加粗
      axis.title.x = element_text(size = 20, face = "bold"),  # x轴标签加粗加大
      axis.title.y = element_text(size = 20, face = "bold"),  # y轴标签加粗加大
      axis.line = element_line(size = 1.5),  # 坐标轴线加粗
      axis.text.x = element_text(size = 20),
      legend.text = element_text(size = 20),  # 图例文字加大加粗
      legend.title = element_text(size = 20, face = "bold"),  # 图例标题加大加粗
      strip.text = element_text(size = 20, face = "bold"),  # facet标签加大加粗
      legend.position = ifelse(show_legend, "bottom", "none"),  # 控制图例显示与否
      legend.box = "horizontal",  # 图例横向排列
      legend.spacing.x = unit(0.5, "cm"),  # 调整图例之间的间距
      legend.direction = "horizontal"  # 确保图例水平排列
    )
  
  # 如果显示图例，则添加图例的具体设置
  if (show_legend) {
    p <- p + guides(color = guide_legend(ncol = length(unique(data$m)) + 1))  # 确保图例项不换行
  }
  
  return(p)
}

## burn in length= 0--##
data_hu_DIM_rmse0 <- read.csv("result_Hu_DIM_rmse_burn0.csv")[,-1]

plot_hu_DIM_rmse0 <- create_plot(data_hu_DIM_rmse0, y_var = "RMSE", title = "RMSE: 0-Burn-in")

ggsave("plot_hu_DIM_rmse0.png", plot_hu_DIM_rmse0, width = 12, height = 8, dpi = 300) 

data_hu_DIM_bias0 <- read.csv("result_Hu_DIM_bias_burn0.csv")[,-1]

plot_hu_DIM_bias0 <- create_plot(data_hu_DIM_bias0, y_var = "Bias", title = "Bias: 0-Burn-in")

ggsave("plot_hu_DIM_bias0.png", plot_hu_DIM_bias0, width = 12, height = 8, dpi = 300) 


data_hu_DIM_sd0 <- read.csv("result_Hu_DIM_SD_burn0.csv")[,-1]

plot_hu_DIM_sd0 <- create_plot(data_hu_DIM_sd0, y_var = "SD", title = "SD: 0-Burn-in")

ggsave("plot_hu_DIM_sd0.png", plot_hu_DIM_sd0, width = 12, height = 8, dpi = 300) 


##--burn length= 2--##

data_hu_DIM_rmse1 <- read.csv("result_Hu_DIM_rmse_burn1.csv")[,-1]

plot_hu_DIM_rmse1 <- create_plot(data_hu_DIM_rmse1, y_var = "RMSE", title = "RMSE: 2-Burn-in")

ggsave("plot_hu_DIM_rmse1.png", plot_hu_DIM_rmse1, width = 12, height = 8, dpi = 300) 

data_hu_DIM_bias1 <- read.csv("result_Hu_DIM_bias_burn1.csv")[,-1]

plot_hu_DIM_bias1 <- create_plot(data_hu_DIM_bias1, y_var = "Bias", title = "Bias: 2-Burn-in")

ggsave("plot_hu_DIM_bias1.png", plot_hu_DIM_bias1, width = 12, height = 8, dpi = 300) 


data_hu_DIM_sd1 <- read.csv("result_Hu_DIM_SD_burn1.csv")[,-1]

plot_hu_DIM_sd1 <- create_plot(data_hu_DIM_sd1, y_var = "SD", title = "SD: 2-Burn-in")

ggsave("plot_hu_DIM_sd1.png", plot_hu_DIM_sd1, width = 12, height = 8, dpi = 300) 



##--burn length= 3--##

data_hu_DIM_rmse2 <- read.csv("result_Hu_DIM_rmse_burn2.csv")[,-1]

plot_hu_DIM_rmse2 <- create_plot(data_hu_DIM_rmse2, y_var = "RMSE", title = "RMSE: 3-Burn-in")

ggsave("plot_hu_DIM_rmse2.png", plot_hu_DIM_rmse2, width = 12, height = 8, dpi = 300) 

data_hu_DIM_bias2 <- read.csv("result_Hu_DIM_bias_burn2.csv")[,-1]

plot_hu_DIM_bias2 <- create_plot(data_hu_DIM_bias2, y_var = "Bias", title = "Bias: 3-Burn-in")

ggsave("plot_hu_DIM_bias2.png", plot_hu_DIM_bias2, width = 12, height = 8, dpi = 300) 


data_hu_DIM_sd2 <- read.csv("result_Hu_DIM_SD_burn2.csv")[,-1]

plot_hu_DIM_sd2 <- create_plot(data_hu_DIM_sd2, y_var = "SD", title = "SD: 3-Burn-in")

ggsave("plot_hu_DIM_sd2.png", plot_hu_DIM_sd2, width = 12, height = 8, dpi = 300) 



##--burn length= 6--##

data_hu_DIM_rmse3 <- read.csv("result_Hu_DIM_rmse_burn3.csv")[,-1]

plot_hu_DIM_rmse3 <- create_plot(data_hu_DIM_rmse3, y_var = "RMSE", title = "RMSE: 6-Burn-in")

ggsave("plot_hu_DIM_rmse3.png", plot_hu_DIM_rmse3, width = 12, height = 8, dpi = 300) 

data_hu_DIM_bias3 <- read.csv("result_Hu_DIM_bias_burn3.csv")[,-1]

plot_hu_DIM_bias3 <- create_plot(data_hu_DIM_bias3, y_var = "Bias", title = "Bias: 6-Burn-in")

ggsave("plot_hu_DIM_bias3.png", plot_hu_DIM_bias3, width = 12, height = 8, dpi = 300) 


data_hu_DIM_sd3 <- read.csv("result_Hu_DIM_SD_burn3.csv")[,-1]

plot_hu_DIM_sd3 <- create_plot(data_hu_DIM_sd3, y_var = "SD", title = "SD: 6-Burn-in")

ggsave("plot_hu_DIM_sd3.png", plot_hu_DIM_sd3, width = 12, height = 8, dpi = 300) 


##--IS ---## 
data_IS_linear_mse<-read.csv("result_IS_rmse_linear.csv")[,-1]

plot_IS_linear_mse <- create_plot(data_IS_linear_mse, y_var = "RMSE", title = "IS-Linear DGP")

ggsave("plot_IS_linear_mse.png", plot_IS_linear_mse, width = 12, height = 8, dpi = 300) 


data_IS_linear_bias<-read.csv("result_IS_bias_linear.csv")[,-1]

plot_IS_linear_bias <- create_plot(data_IS_linear_bias, y_var = "Bias", title = "IS-Linear DGP")

ggsave("plot_IS_linear_bias.png", plot_IS_linear_bias, width = 12, height = 8, dpi = 300) 


data_IS_linear_sd<-read.csv("result_IS_SD_linear.csv")[,-1]

plot_IS_linear_sd <- create_plot(data_IS_linear_sd, y_var = "SD", title = "IS-Linear DGP")

ggsave("plot_IS_linear_sd.png", plot_IS_linear_sd, width = 12, height = 8, dpi = 300) 


# MSE/RMSE 结果
data_IS_nonlinear_mse <- read.csv("result_IS_rmse_nonlinear.csv")[,-1]
plot_IS_nonlinear_mse <- create_plot(data_IS_nonlinear_mse, y_var = "RMSE", title = "IS-Nonlinear DGP")
ggsave("plot_IS_nonlinear_mse.png", plot_IS_nonlinear_mse, width = 12, height = 8, dpi = 300) 

# Bias 结果
data_IS_nonlinear_bias <- read.csv("result_IS_bias_nonlinear.csv")[,-1]
plot_IS_nonlinear_bias <- create_plot(data_IS_nonlinear_bias, y_var = "Bias", title = "IS-Nonlinear DGP")
ggsave("plot_IS_nonlinear_bias.png", plot_IS_nonlinear_bias, width = 12, height = 8, dpi = 300) 

# SD 结果
data_IS_nonlinear_sd <- read.csv("result_IS_SD_nonlinear.csv")[,-1]
plot_IS_nonlinear_sd <- create_plot(data_IS_nonlinear_sd, y_var = "SD", title = "IS-Nonlinear DGP")
ggsave("plot_IS_nonlinear_sd.png", plot_IS_nonlinear_sd, width = 12, height = 8, dpi = 300) 

# --- MIS-HT Bojinov ---
data_MIS_linear_mse<-read.csv("result_HT_rmse_linear.csv")[,-1]

plot_MIS_linear_mse <- create_plot(data_MIS_linear_mse, y_var = "RMSE", title = "MIS-Linear DGP")

ggsave("plot_MIS_linear_mse.png", plot_MIS_linear_mse, width = 12, height = 8, dpi = 300) 


data_MIS_linear_bias<-read.csv("result_HT_bias_linear.csv")[,-1]

plot_MIS_linear_bias <- create_plot(data_MIS_linear_bias, y_var = "Bias", title = "MIS-Linear DGP")

ggsave("plot_MIS_linear_bias.png", plot_MIS_linear_bias, width = 12, height = 8, dpi = 300) 


data_MIS_linear_sd<-read.csv("result_HT_SD_linear.csv")[,-1]

plot_MIS_linear_sd <- create_plot(data_MIS_linear_sd, y_var = "SD", title = "MIS-Linear DGP")

ggsave("plot_MIS_linear_sd.png", plot_MIS_linear_sd, width = 12, height = 8, dpi = 300) 



# MSE/RMSE 结果 - MIS非线性版本
data_MIS_nonlinear_mse <- read.csv("result_HT_rmse_nonlinear.csv")[,-1]
plot_MIS_nonlinear_mse <- create_plot(data_MIS_nonlinear_mse, y_var = "RMSE", title = "MIS-Nonlinear DGP")
ggsave("plot_MIS_nonlinear_mse.png", plot_MIS_nonlinear_mse, width = 12, height = 8, dpi = 300) 

# Bias 结果 - MIS非线性版本
data_MIS_nonlinear_bias <- read.csv("result_HT_bias_nonlinear.csv")[,-1]
plot_MIS_nonlinear_bias <- create_plot(data_MIS_nonlinear_bias, y_var = "Bias", title = "MIS-Nonlinear DGP")
ggsave("plot_MIS_nonlinear_bias.png", plot_MIS_nonlinear_bias, width = 12, height = 8, dpi = 300) 

# SD 结果 - MIS非线性版本
data_MIS_nonlinear_sd <- read.csv("result_HT_SD_nonlinear.csv")[,-1]
plot_MIS_nonlinear_sd <- create_plot(data_MIS_nonlinear_sd, y_var = "SD", title = "MIS-Nonlinear DGP")
ggsave("plot_MIS_nonlinear_sd.png", plot_MIS_nonlinear_sd, width = 12, height = 8, dpi = 300)



##--- summary all methods---##
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(dplyr)

#---读取线性模型的结果---#
data_ate_rmse_linear_ols <- read.csv("result_OLS_rmse_linear.csv")[,-1]
data_ate_rmse_linear_lstd <- read.csv("result_MLSTD_rmse_linear.csv")[,-1]
data_ate_rmse_linear_drl <- read.csv("result_DRL_MLSTD_rmse_linear.csv")[,-1]
data_ate_rmse_linear_ht <- read.csv("result_HT_rmse_linear.csv")[,-1]  # 添加 HT 数据
data_ate_rmse_linear_Is <- read.csv("result_IS_rmse_linear.csv")[,-1]  # 添加 IS 数据
data_hu_DIM_rmse_linear <- read.csv("result_Hu_DIM_rmse_burn1_linear.csv")[,-1]


#----读取非线性模型的结果---#
data_ate_rmse_nonlinear_ols <- read.csv("result_OLS_rmse_nonlinear.csv")[,-1]
data_ate_rmse_nonlinear_lstd <- read.csv("result_MLSTD_rmse_nonlinear.csv")[,-1]
data_ate_rmse_nonlinear_drl <- read.csv("result_DRL_MLSTD_rmse_nonlinear.csv")[,-1]
data_ate_rmse_nonlinear_ht <- read.csv("result_HT_rmse_nonlinear.csv")[,-1]  # 添加 HT 数据
data_ate_rmse_nonlinear_Is <- read.csv("result_IS_rmse_nonlinear.csv")[,-1]  # 添加 IS 数据
data_hu_DIM_rmse_nonlinear <- read.csv("result_Hu_DIM_rmse_burn1_nonlinear.csv")[,-1]


# 添加方法和模型类型信息
data_ate_rmse_linear_ols$Methods_name <- "OLS"
data_ate_rmse_linear_lstd$Methods_name <- "LSTD"
data_ate_rmse_linear_drl$Methods_name <- "DRL"
data_ate_rmse_linear_ht$Methods_name <- "Bojinov et al.(2023)"  # 添加 HT 方法名
data_ate_rmse_linear_Is$Methods_name<-"Xiong et al.(2024)"
data_hu_DIM_rmse_linear$Methods_name<-"Hu et al.(2022)"



data_ate_rmse_linear_ols$DGP <- "Linear"
data_ate_rmse_linear_lstd$DGP <- "Linear"
data_ate_rmse_linear_drl$DGP <- "Linear"
data_ate_rmse_linear_ht$DGP <- "Linear"  # 添加 HT 模型类型
data_ate_rmse_linear_Is$DGP<-"Linear"
data_hu_DIM_rmse_linear$DGP<-"Linear"




data_ate_rmse_nonlinear_ols$Methods_name <- "OLS"
data_ate_rmse_nonlinear_lstd$Methods_name <- "LSTD"
data_ate_rmse_nonlinear_drl$Methods_name <- "DRL"
data_ate_rmse_nonlinear_ht$Methods_name <- "Bojinov et al.(2023)"  # 添加 HT 方法名
data_ate_rmse_nonlinear_Is$Methods_name <-"Xiong et al.(2024)"
data_hu_DIM_rmse_nonlinear$Methods_name<-"Hu et al.(2022)"

data_ate_rmse_nonlinear_ols$DGP <- "Nonlinear"
data_ate_rmse_nonlinear_lstd$DGP <- "Nonlinear"
data_ate_rmse_nonlinear_drl$DGP <- "Nonlinear"
data_ate_rmse_nonlinear_ht$DGP <- "Nonlinear"  # 添加 HT 模型类型
data_ate_rmse_nonlinear_Is$DGP<-"Nonlinear"
data_hu_DIM_rmse_nonlinear$DGP<-"Nonlinear"


# 合并所有数据
all_data <- rbind(
  data_ate_rmse_linear_ols, 
  data_ate_rmse_linear_lstd, 
  data_ate_rmse_linear_drl,
  data_ate_rmse_linear_ht,  # 添加 HT 数据
  data_ate_rmse_linear_Is,
  data_hu_DIM_rmse_linear,
  data_ate_rmse_nonlinear_ols, 
  data_ate_rmse_nonlinear_lstd, 
  data_ate_rmse_nonlinear_drl,
  data_ate_rmse_nonlinear_ht,  # 添加非线性 HT 数据
  data_ate_rmse_nonlinear_Is,
  data_hu_DIM_rmse_nonlinear
)


all_data <- as_tibble(all_data)
all_data<-mutate(all_data,log_mse=log(MSE))
# 将 "rho0" 和 "m" 转换为因子变量，并确保 rho0 的水平一致
# 将 "rho0" 和 "m" 转换为因子变量，并确保 rho0 和 Methods_name 的水平一致



# 确保Methods_name的顺序
all_data <- all_data %>%
  mutate_at(c("rho0", "m", "DGP", "Methods_name"), as_factor) %>%
  mutate(
    rho0 = factor(rho0, levels = c("0.3", "0.5", "0.7", "0.9")),
    Methods_name = factor(Methods_name, 
                          levels = c("DRL", "LSTD", "OLS", "Xiong et al.(2024)", 
                                     "Bojinov et al.(2023)", "Hu et al.(2022)"))
  )

# 自定义调色板 - 高区分度、无粉色、Hu et al.(2022)为黑色
method_colors <- c(
  "DRL" = "#1F78B4",     # 深蓝色
  "LSTD" = "#33A02C",    # 深绿色
  "OLS" = "#E31A1C",     # 红色
  "Xiong et al.(2024)" = "#FF7F00", # 橙色
  "Bojinov et al.(2023)" = "#6A3D9A", # 紫色
  "Hu et al.(2022)" = "#000000"     # 黑色
)

# 创建绘图
p2 <- ggplot(all_data, aes(x = m, y = log_mse, fill = Methods_name, color = Methods_name)) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA) +
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = method_colors) +
  #ggtitle("Comparing LogMSE") +
  facet_grid(. ~ DGP, scales = "free") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 34, face = "bold"),
    axis.text = element_text(size = 30, face = "bold"),
    axis.title.x = element_text(size = 30, face = "bold"),
    axis.title.y = element_text(size = 30, face = "bold"),
    axis.line = element_line(size = 1.5),
    axis.text.x = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.title = element_blank(),  # 去除图例标题
   # legend.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 30, face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.spacing.x = unit(0.5, "cm"),
    legend.direction = "horizontal"
  )

p2

ggsave("compare_LogMSE_all.png", p2, width = 12, height = 10, dpi = 300)  # 非线性 OLS 图





