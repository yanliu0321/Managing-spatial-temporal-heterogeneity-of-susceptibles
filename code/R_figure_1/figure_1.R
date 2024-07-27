
rm(list = ls())
# 安装必要的R包
#install.packages("sf")
#install.packages("ggplot2")
#install.packages("raster")
#install.packages("readr")
#install.packages("dplyr")
#install.packages("jsonlite")
#install.packages("stringr")

# 加载必要的R包
library(sf)
library(ggplot2)
library(raster)
library(readr)
library(dplyr)
library(jsonlite)
library(stringr) 

#县
json_data <- fromJSON("D:\\Rstudio_data\\map_figure1\\china_country.geojson")


# "黄浦区"   "徐汇区"   "长宁区"   "静安区"   "普陀区"   "虹口区"   "杨浦区"   "闵行区"
# "宝山区"   "嘉定区"   "浦东新区" "金山区"   "松江区"   "青浦区"   "奉贤区"   "崇明区"
# 下面n中的xxxx代表上面对应区域的行政区划代码，如需运行程序，请自行更改
n <- c(which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx"),
       which(json_data[["features"]][["properties"]][["gb"]]=="156xxxx")
)
coordinates <- json_data[["features"]][["geometry"]][n,]
data1 <- data.frame(lng= coordinates$coordinates[[1]][,,,1],lat=coordinates$coordinates[[1]][,,,2])
data2 <- data.frame(lng= coordinates$coordinates[[2]][,,,1],lat=coordinates$coordinates[[2]][,,,2])
data3 <- data.frame(lng= coordinates$coordinates[[3]][,,,1],lat=coordinates$coordinates[[3]][,,,2])
data4 <- data.frame(lng= coordinates$coordinates[[4]][,,,1],lat=coordinates$coordinates[[4]][,,,2])
data5 <- data.frame(lng= coordinates$coordinates[[5]][,,,1],lat=coordinates$coordinates[[5]][,,,2])
data6 <- data.frame(lng= coordinates$coordinates[[6]][,,,1],lat=coordinates$coordinates[[6]][,,,2])
data7 <- data.frame(lng= coordinates$coordinates[[7]][,,,1],lat=coordinates$coordinates[[7]][,,,2])
data8 <- data.frame(lng= coordinates$coordinates[[8]][,,,1],lat=coordinates$coordinates[[8]][,,,2])
data9 <- data.frame(lng= coordinates$coordinates[[9]][,,,1],lat=coordinates$coordinates[[9]][,,,2])
data10 <- data.frame(lng= coordinates$coordinates[[10]][,,,1],lat=coordinates$coordinates[[10]][,,,2])

data11_1 <- data.frame(lng= coordinates$coordinates[[11]][[1]][,,1],lat=coordinates$coordinates[[11]][[1]][,,2])
data11_2 <- data.frame(lng= coordinates$coordinates[[11]][[2]][,,1],lat=coordinates$coordinates[[11]][[2]][,,2])
data11_3 <- data.frame(lng= coordinates$coordinates[[11]][[3]][,,1],lat=coordinates$coordinates[[11]][[3]][,,2])

data12 <- data.frame(lng= coordinates$coordinates[[12]][,,,1],lat=coordinates$coordinates[[12]][,,,2])
data13 <- data.frame(lng= coordinates$coordinates[[13]][,,,1],lat=coordinates$coordinates[[13]][,,,2])
data14 <- data.frame(lng= coordinates$coordinates[[14]][,,,1],lat=coordinates$coordinates[[14]][,,,2])
data15 <- data.frame(lng= coordinates$coordinates[[15]][,,,1],lat=coordinates$coordinates[[15]][,,,2])

data16_1 <- data.frame(lng= coordinates$coordinates[[16]][[1]][,,1],lat=coordinates$coordinates[[16]][[1]][,,2])
data16_2 <- data.frame(lng= coordinates$coordinates[[16]][[2]][,,1],lat=coordinates$coordinates[[16]][[2]][,,2])
data16_3 <- data.frame(lng= coordinates$coordinates[[16]][[3]][,,1],lat=coordinates$coordinates[[16]][[3]][,,2])

############ 计算每个区域的质心，用于在地图中标记各个区名字 
# 建立数据框列表
data0 <-list(
  data1, data2, data3, data4, data5, data6, data7, data8, data9, data10,
  data11_3,
  data12,  data13, data14, data15,
  data16_3
)
# 函数：将数据框转换为多边形
convert_to_polygon <- function(df) {
  if (is.list(df)) {
    lapply(df, function(sub_df) st_polygon(list(as.matrix(sub_df))))
  } else {
    st_polygon(list(as.matrix(df)))
  }
}

multi_polygons0 <- convert_to_polygon(data0)

# 创建带有区名称的sf对象
geometry0 <- st_sfc(multi_polygons0, crs = 4326)
centroids <- st_centroid(geometry0)
# 将质心数据框转化为数据框
centroid_df <- as.data.frame(st_coordinates(centroids))
colnames(centroid_df) <- c('lng','lat')
name <-  c('Huangpu',	'Xuhui',	'Changning',	'Jingan',	'Putuo',	'Hongkou',
           'Yangpu',	'Minghang',	'Baoshan',	'Jiading',	'Pudong',	'Jinshan',
           'Songjiang',	'Qingpu',	'Fengxian',	'Chongming')
name1 <-  c('HP',	'XH',	'CN',	'JA',	'PT',	'HK',
            'YP',	'Minghang',	'Baoshan',	'Jiading',	'Pudong',	'Jinshan',
            'Songjiang',	'Qingpu',	'Fengxian',	'Chongming')
centroid_df$name1 <- name1

#################################################################################
# 建立数据框列表 data1 到 data16_3
data <- list(
  data1, data2, data3, data4, data5, data6, data7, data8, data9, data10,
  data11_1, data11_2, data11_3,
  data12,  data13, data14, data15,
  data16_1, data16_2, data16_3
)

multi_polygons <- convert_to_polygon(data)

# 创建带有区名称的sf对象
geometry <- st_sfc(multi_polygons, crs = 4326)
geometry <- st_sf(geometry = geometry, district = rep(name, times = c(1,1,1,1,1,1,1,1,1,1,3,1,1,1,1,3)))



###############################   设置各个区上色数据 
epidemic_data <-  read_csv("D:\\Rstudio_data\\map_figure1\\epidemic data.csv")
# 对表格里的每列数据进行累计
epidemic_data_cumulative <- epidemic_data %>%
  mutate(across(-Time, cumsum))

# 自定义颜色方案
color_scheme <- c("0" = "#ffffff", "1-3" = "#a9de00", "3-10" = "#3acf3d", 
                  "10-20" = "#00ae7b", "20-40" = "#006e8f", "40-100" = "#4a3785", ">100" = "#4f0057")


##  a
epidemic_df1 <- data.frame(name,reported_case=t(epidemic_data_cumulative[which(epidemic_data_cumulative$Time=='2022/3/1'),-1]))

# 合并geometry数据和epidemic_data数据,并创建区间标签
geometry1 <- geometry
geometry1 <- geometry1 %>%
  left_join(epidemic_df1, by = c("district" = "name" )) %>%
  mutate(case_category = cut(reported_case, 
                             breaks = c(0, 1, 3, 10, 20, 40, 100, Inf), 
                             labels = c("0", "1-3", "3-10", "10-20", "20-40", "40-100", ">100"),
                             right = FALSE))


# 绘制上海市地图和热图
ggplot(data = geometry1) +
  geom_sf(aes(fill = case_category)) +
  geom_text(data = centroid_df, aes(x = lng, y = lat, label = name1), size = 5, color = "black") +
  scale_fill_manual(values = color_scheme) +
  labs(title = "(a) March 1, 2022", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(family = "serif"),
        plot.title = element_text(size = 20, face = "bold"),  # 改变标题字体大小
        legend.position = "left",  # 将图例放置在右侧
        legend.title = element_text(size = 18),  # 图例标题字体大小
        legend.text = element_text(size = 18),  # 图例文本字体大小
        legend.box.background = element_rect(color = "black"),  # 图例边框颜色
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.text.x = element_blank(),  # 去掉横轴刻度
        axis.text.y = element_blank(),  # 去掉纵轴刻度
        axis.ticks = element_blank(),  # 去掉坐标轴的刻度
        axis.title.x = element_blank(),  # 去掉横轴标签
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(title = "HP: Huangpu\nXH: Xuhui\nCN: Changning\nJA: Jingan\nPT: Putuo\nHK: Hongkou\nYP: Yangpu\n\n\nCases"))  # 添加额外的说明文字



##  b
epidemic_df2 <- data.frame(name,reported_case=t(epidemic_data_cumulative[which(epidemic_data_cumulative$Time=='2022/3/3'),-1]))

# 合并geometry数据和epidemic_data数据,并创建区间标签
geometry2 <- geometry
geometry2 <- geometry2 %>%
  left_join(epidemic_df2, by = c("district" = "name" ))  %>%
  mutate(case_category = cut(reported_case, 
                             breaks = c(0, 1, 3, 10, 20, 40, 100, Inf), 
                             labels = c("0", "1-3", "3-10", "10-20", "20-40", "40-100", ">100"),
                             right = FALSE))

# 绘制上海市地图和热图
ggplot(data = geometry2) +
  geom_sf(aes(fill = case_category)) +
  geom_text(data = centroid_df, aes(x = lng, y = lat, label = name1), size = 5, color = "black") +
  scale_fill_manual(values = color_scheme) +
  labs(title = "(b) March 3, 2022", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(family = "serif"),
        plot.title = element_text(size = 20, face = "bold"),  # 改变标题字体大小
        legend.position = "left",  # 将图例放置在右侧
        legend.title = element_text(size = 18),  # 图例标题字体大小
        legend.text = element_text(size = 18),  # 图例文本字体大小
        legend.box.background = element_rect(color = "black"),  # 图例边框颜色
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.text.x = element_blank(),  # 去掉横轴刻度
        axis.text.y = element_blank(),  # 去掉纵轴刻度
        axis.ticks = element_blank(),  # 去掉坐标轴的刻度
        axis.title.x = element_blank(),  # 去掉横轴标签
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(title = "HP: Huangpu\nXH: Xuhui\nCN: Changning\nJA: Jingan\nPT: Putuo\nHK: Hongkou\nYP: Yangpu\n\n\nCases"))  # 添加额外的说明文字


##  c
epidemic_df3 <- data.frame(name,reported_case=t(epidemic_data_cumulative[which(epidemic_data_cumulative$Time=='2022/3/6'),-1]))

# 合并geometry数据和epidemic_data数据,并创建区间标签
geometry3 <- geometry
geometry3 <- geometry3 %>%
  left_join(epidemic_df3, by = c("district" = "name" ))  %>%
  mutate(case_category = cut(reported_case, 
                             breaks = c(0, 1, 3, 10, 20, 40, 100, Inf), 
                             labels = c("0", "1-3", "3-10", "10-20", "20-40", "40-100", ">100"),
                             right = FALSE))

# 绘制上海市地图和热图
ggplot(data = geometry3) +
  geom_sf(aes(fill = case_category)) +
  geom_text(data = centroid_df, aes(x = lng, y = lat, label = name1), size = 5, color = "black") +
  scale_fill_manual(values = color_scheme) +
  labs(title = "(c) March 6, 2022", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(family = "serif"),
        plot.title = element_text(size = 20, face = "bold"),  # 改变标题字体大小
        legend.position = "left",  # 将图例放置在右侧
        legend.title = element_text(size = 18),  # 图例标题字体大小
        legend.text = element_text(size = 18),  # 图例文本字体大小
        legend.box.background = element_rect(color = "black"),  # 图例边框颜色
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.text.x = element_blank(),  # 去掉横轴刻度
        axis.text.y = element_blank(),  # 去掉纵轴刻度
        axis.ticks = element_blank(),  # 去掉坐标轴的刻度
        axis.title.x = element_blank(),  # 去掉横轴标签
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(title = "HP: Huangpu\nXH: Xuhui\nCN: Changning\nJA: Jingan\nPT: Putuo\nHK: Hongkou\nYP: Yangpu\n\n\nCases"))  # 添加额外的说明文字


##  d
epidemic_df4 <- data.frame(name,reported_case=t(epidemic_data_cumulative[which(epidemic_data_cumulative$Time=='2022/3/8'),-1]))

# 合并geometry数据和epidemic_data数据,并创建区间标签
geometry4 <- geometry
geometry4 <- geometry4 %>%
  left_join(epidemic_df4, by = c("district" = "name" ))  %>%
  mutate(case_category = cut(reported_case, 
                             breaks = c(0, 1, 3, 10, 20, 40, 100, Inf), 
                             labels = c("0", "1-3", "3-10", "10-20", "20-40", "40-100", ">100"),
                             right = FALSE))

# 绘制上海市地图和热图
ggplot(data = geometry4) +
  geom_sf(aes(fill = case_category)) +
  geom_text(data = centroid_df, aes(x = lng, y = lat, label = name1), size = 5, color = "black") +
  scale_fill_manual(values = color_scheme) +
  labs(title = "(d) March 8, 2022", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(family = "serif"),
        plot.title = element_text(size = 20, face = "bold"),  # 改变标题字体大小
        legend.position = "left",  # 将图例放置在右侧
        legend.title = element_text(size = 18),  # 图例标题字体大小
        legend.text = element_text(size = 18),  # 图例文本字体大小
        legend.box.background = element_rect(color = "black"),  # 图例边框颜色
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.text.x = element_blank(),  # 去掉横轴刻度
        axis.text.y = element_blank(),  # 去掉纵轴刻度
        axis.ticks = element_blank(),  # 去掉坐标轴的刻度
        axis.title.x = element_blank(),  # 去掉横轴标签
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(title = "HP: Huangpu\nXH: Xuhui\nCN: Changning\nJA: Jingan\nPT: Putuo\nHK: Hongkou\nYP: Yangpu\n\n\nCases"))  # 添加额外的说明文字


##  e
epidemic_df5 <- data.frame(name,reported_case=t(epidemic_data_cumulative[which(epidemic_data_cumulative$Time=='2022/3/13'),-1]))

# 合并geometry数据和epidemic_data数据,并创建区间标签
geometry5 <- geometry
geometry5 <- geometry5 %>%
  left_join(epidemic_df5, by = c("district" = "name" ))  %>%
  mutate(case_category = cut(reported_case, 
                             breaks = c(0, 1, 3, 10, 20, 40, 100, Inf), 
                             labels = c("0", "1-3", "3-10", "10-20", "20-40", "40-100", ">100"),
                             right = FALSE))

# 绘制上海市地图和热图
ggplot(data = geometry5) +
  geom_sf(aes(fill = case_category)) +
  geom_text(data = centroid_df, aes(x = lng, y = lat, label = name1), size = 5, color = "black") +
  scale_fill_manual(values = color_scheme) +
  labs(title = "(e) March 13, 2022", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(family = "serif"),
        plot.title = element_text(size = 20, face = "bold"),  # 改变标题字体大小
        legend.position = "left",  # 将图例放置在右侧
        legend.title = element_text(size = 18),  # 图例标题字体大小
        legend.text = element_text(size = 18),  # 图例文本字体大小
        legend.box.background = element_rect(color = "black"),  # 图例边框颜色
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.text.x = element_blank(),  # 去掉横轴刻度
        axis.text.y = element_blank(),  # 去掉纵轴刻度
        axis.ticks = element_blank(),  # 去掉坐标轴的刻度
        axis.title.x = element_blank(),  # 去掉横轴标签
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(title = "HP: Huangpu\nXH: Xuhui\nCN: Changning\nJA: Jingan\nPT: Putuo\nHK: Hongkou\nYP: Yangpu\n\n\nCases"))  # 添加额外的说明文字


##  f
epidemic_df6 <- data.frame(name,reported_case=t(epidemic_data_cumulative[which(epidemic_data_cumulative$Time=='2022/3/17'),-1]))

# 合并geometry数据和epidemic_data数据,并创建区间标签
geometry6 <- geometry
geometry6 <- geometry6 %>%
  left_join(epidemic_df6, by = c("district" = "name" ))  %>%
  mutate(case_category = cut(reported_case, 
                             breaks = c(0, 1, 3, 10, 20, 40, 100, Inf), 
                             labels = c("0", "1-3", "3-10", "10-20", "20-40", "40-100", ">100"),
                             right = FALSE))

# 绘制上海市地图和热图
ggplot(data = geometry6) +
  geom_sf(aes(fill = case_category)) +
  geom_text(data = centroid_df, aes(x = lng, y = lat, label = name1), size = 5, color = "black") +
  scale_fill_manual(values = color_scheme) +
  labs(title = "(f) March 17, 2022", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(text = element_text(family = "serif"),
        plot.title = element_text(size = 20, face = "bold"),  # 改变标题字体大小
        legend.position = "left",  # 将图例放置在右侧
        legend.title = element_text(size = 18),  # 图例标题字体大小
        legend.text = element_text(size = 18),  # 图例文本字体大小
        legend.box.background = element_rect(color = "black"),  # 图例边框颜色
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),  # 去掉次网格线
        axis.text.x = element_blank(),  # 去掉横轴刻度
        axis.text.y = element_blank(),  # 去掉纵轴刻度
        axis.ticks = element_blank(),  # 去掉坐标轴的刻度
        axis.title.x = element_blank(),  # 去掉横轴标签
        axis.title.y = element_blank()) + 
  guides(fill = guide_legend(title = "HP: Huangpu\nXH: Xuhui\nCN: Changning\nJA: Jingan\nPT: Putuo\nHK: Hongkou\nYP: Yangpu\n\n\nCases"))  # 添加额外的说明文字



