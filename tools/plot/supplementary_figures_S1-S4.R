
#Making supplementary figures S1-4 for the codon bias paper


#install.packages('gcookbook')
#install.packages('cowplot')
library(ggplot2)
library(gcookbook)
library(scales)
library(grid)
library(cowplot)
library(stringr)



library(rmarkdown)
render("S4FGH_sqrtPPR_function_of_tAI_smallCDSdeltaG_sqrt.rmd")


setwd("C:/Users/Magiz/Desktop")


y_title<-str_wrap("Fixed effects slope", width=13)


############## Figure S1 ############## 
#these are the figures for deciding which of CAI/nCAI/tAI/ntAI and ensemble/mfe/bp_prob to use
#this is the sqrtPPR version of figure 1


#types<-c("CAI", "nlCAI", "tAI", "ntAI", "CAI", "nlCAI", "tAI", "ntAI") #syn only tAI confint had some warnings
#group<-c("all", "all", "all", "all", "synonymous only", "synonymous only", "synonymous only", "synonymous only")
#slopes<-c(17.6296, 17.2315, 39.277, 53.540, 13.094, 12.638, 43.918, 45.652)
#lower<-c(12.39258, 12.01506, 30.56147, 39.47122, 5.656207, 5.306476, 27.28993, 13.66286)
#upper<-c(22.98069, 22.56631, 48.15024, 67.71153, 20.78859, 20.22202, 62.48222, 77.50339)

# cb data
cb.data <- data.frame(measure = c("CAI","CAI","nlCAI","nlCAI","tAI","tAI","ntAI","ntAI"), 
                      set = str_wrap(c(rep(x = c("all (1620 genes)","synonymous only (185 genes)"), times = 4)), width = 16),
                      slope = c(17.6296, 13.094, 17.2315, 12.638, 39.277, 43.918, 53.540, 45.652),
                      down95ci = c(12.39258, 5.656207, 12.01506, 5.306476, 30.56147, 27.28993, 39.47122, 13.66286),
                      up95ci = c(22.98069, 20.78859, 22.56631, 20.22202, 48.15024, 62.48222, 67.71153, 77.50339)
)
cb.data$measure <- factor(cb.data$measure, levels = c("CAI","nlCAI","tAI","ntAI"), ordered = TRUE)
str(cb.data)



#types<-c("ensemble deltaG", "ensemble deltaG", "mfe deltaG", "mfe deltaG", "mean bp prob", "mean bp prob")
#group<-c("all", "synonymous only", "all", "synonymous only", "all", "synonymous only")
#slopes<-c(-13.089, -27.625, -11.039, -19.008, 4.956, 4.794)
#lower<-c(-16.81145, -44.51314, -14.6455, -34.30597, 1.006035, -3.798673)
#upper<-c(-9.371116, -11.08517, -7.436932, -3.623943, 8.918971, 13.93096)

# rna data
rna.data <- data.frame(measure = c("ensemble deltaG","ensemble deltaG","mfe deltaG","mfe deltaG","mean bp prob","mean bp prob"), 
                       set = str_wrap(c(rep(x = c("all (1458 genes)","synonymous only (176 genes)"), times = 3)), width = 16),
                       slope = c(-13.089, -27.625, -11.039, -19.008, 4.956, 4.794),
                       down95ci = c(-16.81145, -44.51314, -14.6455, -34.30597, 1.006035, -3.798673),
                       up95ci = c(-9.371116, -11.08517, -7.436932, -3.623943, 8.918971, 13.93096)
)
rna.data$measure <- factor(rna.data$measure, levels = c("mfe deltaG","ensemble deltaG","mean bp prob"), ordered = TRUE)
str(rna.data)


colorBlindGrey8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
tiff(file="figureS1.tiff", width = 13, height = 8, units = 'cm', res = 300)
pd <- position_dodge(1)
fig1a <- ggplot(cb.data, aes(x = measure, y = slope, fill = set)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab("Fixed effects slope") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="top" )+
  labs(fill = str_wrap("Gene Set", width = 4))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[3]))

fig1b <- ggplot(rna.data, aes(x = measure, y = slope, fill = set)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab("Fixed effects slope") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="top")+
  labs(fill = str_wrap("Gene Set", width = 4))+
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))+
  scale_fill_manual(values=c(colorBlindGrey8[7], colorBlindGrey8[2]))

plot_grid(fig1a, fig1b, labels = c('A', 'B'), label_size = 12, nrow = 1)
dev.off()




############## Figure S2 ############## 
#this is the sqrtPPR version of Figure 2
#I only did the models for A-c

# tAI data for high and low median detlaG
#types<-c("tAI for low deltaG genes (723)", "tAI for high deltaG genes (724)") 
#slopes<-c(66.250, 18.986)
#lower<-c(50.11935, 9.603971)
#upper<-c(82.75412, 28.52574)


tAI.data <- data.frame(
  set = str_wrap(c("tAI for high deltaG genes (724)","tAI for low deltaG genes (723)"), width = 10),
  slope = c(18.986, 66.250),
  down95ci = c(9.603971, 50.11935),
  up95ci = c(28.52574, 82.75412)
)
tAI.data$set <- factor(tAI.data$set, levels = str_wrap(c("tAI for low deltaG genes (723)","tAI for high deltaG genes (724)"), width = 10), ordered = TRUE)
str(tAI.data)

# deltaG data for high and low median tAI
#types<-c("deltaG for low tAI genes (723)", "deltaG for high tAI genes (724)") 
#slopes<-c(-6.0738, -22.028)
#lower<-c(-9.224195, -28.53739)
#upper<-c(-3.050558, -15.47559)

deltaG.data <- data.frame(
  set = str_wrap(c("deltaG for low tAI genes (723)","deltaG for high tAI genes (724)"), width = 10),
  slope = c(-6.0738, -22.028),
  down95ci = c(-9.224195, -28.53739),
  up95ci = c(-3.050558, -15.47559)
)
#tAI.data$set <- factor(tAI.data$set, levels = c("low mF (724 genes)","high mF (723 genes)"), ordered = TRUE)
str(deltaG.data)

# tAI, deltaG, and tAI:deltaG
##                  Fixed Effects Slope    2.5%          97.5%
## tAI               11.819               -0.1861483     23.84029   (R slope was not scaled)
## ensembleG         6.902                -3.51263       17.24587   (R slope was not scaled)
## tAI:ensembleG    -48.122               -72.83125     -23.30038   (R slope was not scaled)
interact.m1.data <- data.frame(
  term = c("tAI","deltaG","tAI:deltaG"),
  slope = c(11.819, 6.902, -48.122),
  down95ci = c(-0.1861483, -3.51263, -72.83125),
  up95ci = c(23.84029, 17.24587, -23.30038)
)
interact.m1.data$term <- factor(interact.m1.data$term, levels = c("tAI","deltaG","tAI:deltaG"), ordered = TRUE)
str(interact.m1.data)


# combined figure 2
tiff(file="figureS2.tiff", width = 13, height = 10, units = 'cm', res = 300)
pd <- position_dodge(1)
fig2a <- ggplot(tAI.data, aes(x = set, y = slope)) + 
  geom_bar(stat="identity", position=pd, fill = colorBlindGrey8[6]) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))
fig2b <- ggplot(deltaG.data, aes(x = set, y = slope)) + 
  geom_bar(stat="identity", position=pd, fill = colorBlindGrey8[7]) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))
fig2c <- ggplot(interact.m1.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 8))+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))

plot_grid(fig2a, fig2b, fig2c, labels = c('A', 'B', 'C'), label_size = 12, nrow = 2, hjust=-2, vjust=1)
dev.off()






############## Figure S3 ############## 
#This is the sqrtPPR version of figure 3


#sqrtPPR ~ L5tAI + ensembleG + L5tAI:ensembleG + (0+L5tAI|gene) + (1|gene)
##
##                  Fixed Effects Slope    2.5%          97.5%
## L5tAI             -3.404               -12.57192      6.031097   (R slope was not scaled)
## ensembleG         -2.075               -9.460978      5.279253   (R slope was not scaled)
## L5tAI:ensembleG   -16.174              -32.97105      0.7694393  (R slope was not scaled)
L5.tAI.deltaG.data <- data.frame(
  term = str_wrap(c("5\' + linker-tAI","deltaG","5\' + linker-tAI:deltaG"), width = 8),
  slope = c(-3.404, -2.075, -16.174),
  down95ci = c(-12.57192, -9.460978, -32.97105),
  up95ci = c(6.031097, 5.279253, 0.7694393)
)
L5.tAI.deltaG.data$term <- factor(L5.tAI.deltaG.data$term, levels = str_wrap(c("5\' + linker-tAI","deltaG","5\' + linker-tAI:deltaG"), width = 8), ordered = TRUE)
str(L5.tAI.deltaG.data)

#sqrtPPR ~ D3tAI + ensembleG + D3tAI:ensembleG + (0+D3tAI|gene) + (0+ensembleG|gene) + (0+D3tAI:ensembleG|gene) + (1|gene)
##
##                  Fixed Effects Slope    2.5%          97.5%
## d3tAI             7.227                -7.722094      22.16008   (R slope was not scaled)
## ensembleG         12.251                0.8926793     23.56571   (R slope was not scaled)
## d3tAI:ensembleG  -59.136               -85.71652     -32.46815   (R slope was not scaled)
D3.tAI.deltaG.data <- data.frame(
  term = str_wrap(c("domain + 3\'-tAI","deltaG","domain + 3\'-tAI:deltaG"), width = 8),
  slope = c(7.227, 12.251, -59.136),
  down95ci = c(-7.722094, 0.8926793, -85.71652),
  up95ci = c(22.16008, 23.56571, -32.46815)
)
D3.tAI.deltaG.data$term <- factor(D3.tAI.deltaG.data$term, levels = str_wrap(c("domain + 3\'-tAI","deltaG","domain + 3\'-tAI:deltaG"), width = 8), ordered = TRUE)
str(D3.tAI.deltaG.data)

#sqrtPPR ~ tAI + UTR5 + CDS + UTR3 + tAI:UTR5 + tAI:CDS + tAI:UTR3 + (0+tAI|gene) + (0+UTR3|gene)+(1|gene) ### 1312 genes
##
##                  Fixed Effects Slope    2.5%          97.5%
## tAI              -4.754                -17.76681      8.272454   (R slope was not scaled)
## UTR5              354.18                46.34468      662.2212   (R slope was scaled by 10)
## CDS               253.191               145.7415      360.6425   (R slope was not scaled)
## UTR3              460.85                122.7808      796.4774   (R slope was scaled by 10)
## tAI:UTR5         -1060.43              -1784.47      -336.888    (R slope was scaled by 10)
## tAI:CDS          -827.925              -1087.547     -567.896    (R slope was not scaled)
## tAI:UTR3         -1306.17              -2089.263     -516.8423   (R slope was scaled by 10)
cds.utr.deltaG.tAI.data <- data.frame(
  term = str_wrap(c("tAI","5\'UTR-deltaG","CDS-deltaG","3\'UTR-deltaG","tAI:5\'UTR-deltaG","tAI:CDS-deltaG","tAI:3\'UTR-deltaG"), width = 12),
  slope = c(-4.754, 354.18, 253.191, 460.85, -1060.43, -827.925, -1306.17),
  down95ci = c(-17.76681, 46.34468, 145.7415, 122.7808, -1784.47, -1087.547, -2089.263),
  up95ci = c(8.272454, 662.2212, 360.6425, 796.4774, -336.888, -567.896, -516.8423)
)
cds.utr.deltaG.tAI.data$term <- factor(cds.utr.deltaG.tAI.data$term, levels = str_wrap(c("tAI","5\'UTR-deltaG","CDS-deltaG","3\'UTR-deltaG","tAI:5\'UTR-deltaG","tAI:CDS-deltaG","tAI:3\'UTR-deltaG"), width = 12), ordered = TRUE)
str(cds.utr.deltaG.tAI.data)



# figure 3
colorBlindGrey8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
tiff(file="figureS3.tiff", width = 10, height = 10, units = 'cm', res = 300)
pd <- position_dodge(1)
fig3a <- ggplot(L5.tAI.deltaG.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-90,25)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))
fig3b <- ggplot(D3.tAI.deltaG.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-90,25)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))
fig3c <- ggplot(cds.utr.deltaG.tAI.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6],colorBlindGrey8[7],colorBlindGrey8[7],colorBlindGrey8[7],colorBlindGrey8[1],colorBlindGrey8[1],colorBlindGrey8[1]))

r1 <- plot_grid(fig3a, fig3b, labels = c('A', 'B'), label_size = 12, nrow = 1, hjust=-2, vjust=1)
r2 <- plot_grid(fig3c, labels = c('C'), label_size = 12, nrow = 1, hjust=-2, vjust=0)
plot_grid(r1, r2, nrow = 2)
dev.off()



############## Figure S4 ############## 
# fine-scale local mF

##### logPPR numbers #####
#logPPR ~ logPPR ~ tAI + cap + upA + dwA + dwS + tAI:cap + tAI:upA + tAI:dwA + tAI:dwS + (0+tAI|gene) + (0+cap|gene) + (0+upA|gene) + (0+dwA|gene) + (0+dwS|gene) + (1|gene) ### 774 genes
##
##                  Fixed Effects Slope    2.5%          97.5%
## tAI              1.17509                0.7697925    1.579418
## cap              -0.008028              -0.05284203  0.03665219  (R slope was scaled by 1/10)
## upA              0.024033               -0.008039338 0.05560816  (R slope was scaled by 1/10)
## dwA              0.007211               -0.03176275  0.04621319  (R slope was scaled by 1/10)
## dwS              0.010622               -0.02234939  0.04399511  (R slope was scaled by 1/10)
## tAI:cap          0.006916               -0.09598388  0.1100781   (R slope was scaled by 1/10)
## tAI:upA          -0.051185              -0.1244014   0.02302314  (R slope was scaled by 1/10)
## tAI:dwA          -0.004045              -0.09395915  0.08589954  (R slope was scaled by 1/10)
## tAI:dwS          -0.021283              -0.09959116  0.05611616  (R slope was scaled by 1/10)
local.deltaG.tAI.data <- data.frame(
  term = str_wrap(c("tAI","5\'Cap-deltaG","-9 to +3 Start-deltaG","+4 to +10 Start-deltaG","+1 to +18 Stop-deltaG","tAI:5\'Cap-deltaG","tAI:-9 to +3 Start-deltaG","tAI:+4 to +10 Start-deltaG","tAI:+1 to +18 Stop-deltaG"), width = 10),
  slope = c(1.17509,-0.008028,0.024033,0.007211,0.010622,0.006916,-0.051185,-0.004045,-0.021283),
  down95ci = c(0.7697925,-0.05284203,-0.008039338,-0.03176275,-0.02234939,-0.09598388,-0.1244014,-0.09395915,-0.09959116),
  up95ci = c(1.579418,0.03665219,0.05560816,0.04621319,0.04399511,0.1100781,0.02302314,0.08589954,0.05611616)
)
local.deltaG.tAI.data$term <- factor(local.deltaG.tAI.data$term, levels = str_wrap(c("tAI","5\'Cap-deltaG","-9 to +3 Start-deltaG","+4 to +10 Start-deltaG","+1 to +18 Stop-deltaG","tAI:5\'Cap-deltaG","tAI:-9 to +3 Start-deltaG","tAI:+4 to +10 Start-deltaG","tAI:+1 to +18 Stop-deltaG"), width = 10), ordered = TRUE)
str(local.deltaG.tAI.data)

# checking power to detect effect of mF in 40 bp regions in the CDS
#logPPR ~ tAI + onefour + tAI:onefour + (0+tAI|gene) + (0+onefour|gene) + (1|gene)  #1322 genes (1/4 CDS)
##
##                  Fixed Effects Slope    2.5%          97.5%
##tAI               1.23903                0.9259004     1.551467      (R slope was not scaled)
##onefour          -0.004366              -0.02154977    0.01280856    (R slope was scaled by 1/10)
##tAI:onefour       0.011834              -0.02797082    0.05165319    (R slope was scaled by 1/10)

#logPPR ~ tAI + onehalf + tAI:onehalf + (0+tAI|gene) + (0+onehalf|gene) + (1|gene)  #1313 genes (1/2 CDS)
##
##                  Fixed Effects Slope    2.5%          97.5%
##tAI               1.152326               0.8425403     1.461559      (R slope was not scaled)
##onehalf           0.0004739             -0.01787049    0.01865711    (R slope was scaled by 1/10)
##tAI:onehalf       0.0035978             -0.03855118    0.0461482     (R slope was scaled by 1/10)

#logPPR ~ tAI + threefour + tAI:threefour + (0+tAI|gene) + (0+threefour|gene) + (1|gene)  #1314 genes (3/4 CDS)
##
##                  Fixed Effects Slope    2.5%          97.5%
##tAI               1.067096               0.7561221     1.377538      (R slope was not scaled)
##threefour         0.0009134             -0.0142555     0.01597644    (R slope was scaled by 1/10)
##tAI:threefour     0.0016502             -0.03328128    0.03685957    (R slope was scaled by 1/10)

cds25.int.deltaG.tAI.data <- data.frame(
  term = str_wrap(c("tAI","CDS25-deltaG","tAI:CDS25-deltaG"), width = 12),
  slope = c(1.23903,-0.004366,0.011834),
  down95ci = c(0.9259004,-0.02154977,-0.02797082),
  up95ci = c(1.551467,0.01280856,0.05165319)
)
cds25.int.deltaG.tAI.data$term <- factor(cds25.int.deltaG.tAI.data$term, levels = str_wrap(c("tAI","CDS25-deltaG","tAI:CDS25-deltaG"), width = 12), ordered = TRUE)
str(cds25.int.deltaG.tAI.data)
cds50.int.deltaG.tAI.data <- data.frame(
  term = str_wrap(c("tAI","CDS50-deltaG","tAI:CDS50-deltaG"), width = 12),
  slope = c(1.152326,0.0004739,0.0035978),
  down95ci = c(0.8425403,-0.01787049,-0.03855118),
  up95ci = c(1.461559,0.01865711,0.0461482)
)
cds50.int.deltaG.tAI.data$term <- factor(cds50.int.deltaG.tAI.data$term, levels = str_wrap(c("tAI","CDS50-deltaG","tAI:CDS50-deltaG"), width = 12), ordered = TRUE)
str(cds50.int.deltaG.tAI.data)
cds75.int.deltaG.tAI.data <- data.frame(
  term = str_wrap(c("tAI","CDS75-deltaG","tAI:CDS75-deltaG"), width = 12),
  slope = c(1.067096,0.0009134,0.0016502),
  down95ci = c(0.7561221,-0.0142555,-0.03328128),
  up95ci = c(1.377538,0.01597644,0.03685957)
)
cds75.int.deltaG.tAI.data$term <- factor(cds75.int.deltaG.tAI.data$term, levels = str_wrap(c("tAI","CDS75-deltaG","tAI:CDS75-deltaG"), width = 12), ordered = TRUE)
str(cds75.int.deltaG.tAI.data)


##### sqrtPPR numbers #####
#sqrtPPR ~ tAI + cap + upA + dwA + dwS + tAI:cap + tAI:upA + tAI:dwA + tAI:dwS +(0+tAI|gene)+(0+cap|gene)+(0+upA|gene)+(0+dwA|gene)+(0+dwS|gene) + (1|gene) #774 genes
##
##                  Fixed Effects Slope    2.5%          97.5%
## tAI               25.5994              15.37151       35.79892   (R slope was not scaled)
## cap               0.1257              -0.902655       1.150817   (R slope was not scaled)
## upA               0.3599              -0.3546989      1.062611   (R slope was not scaled)
## dwA              -0.4381              -1.372198       0.4951409  (R slope was not scaled)
## dwS               0.1753              -0.5646391      0.9253702  (R slope was not scaled)
## tAI:cap          -0.4813              -2.842932       1.8865     (R slope was not scaled)
## tAI:upA          -0.7780              -2.407888       0.876483   (R slope was not scaled)
## tAI:dwA           1.2828              -0.8702399      3.439278   (R slope was not scaled)
## tAI:dwS          -0.3834              -2.144795       1.354494   (R slope was not scaled)
local.deltaG.tAI.data.sqrt <- data.frame(
  term = str_wrap(c("tAI","5\'Cap-deltaG","-9 to +3 Start-deltaG","+4 to +10 Start-deltaG","+1 to +18 Stop-deltaG","tAI:5\'Cap-deltaG","tAI:-9 to +3 Start-deltaG","tAI:+4 to +10 Start-deltaG","tAI:+1 to +18 Stop-deltaG"), width = 10),
  slope = c(25.5994,0.1257,0.3599,-0.4381,0.1753,-0.4813,-0.7780,1.2828,-0.3834),
  down95ci = c(15.37151,-0.902655,-0.3546989,-1.372198,-0.5646391,-2.842932,-2.407888,-0.8702399,-2.144795),
  up95ci = c(35.79892,1.150817,1.062611,0.4951409,0.9253702,1.8865,0.876483,3.439278,1.354494)
)
local.deltaG.tAI.data.sqrt$term <- factor(local.deltaG.tAI.data.sqrt$term, levels = str_wrap(c("tAI","5\'Cap-deltaG","-9 to +3 Start-deltaG","+4 to +10 Start-deltaG","+1 to +18 Stop-deltaG","tAI:5\'Cap-deltaG","tAI:-9 to +3 Start-deltaG","tAI:+4 to +10 Start-deltaG","tAI:+1 to +18 Stop-deltaG"), width = 10), ordered = TRUE)
str(local.deltaG.tAI.data.sqrt)


#sqrtPPR ~ tAI + onefour + tAI:onefour + (0+tAI|gene) + (0+onefour|gene) + (1|gene)  #1322 genes (1/4 CDS)
##
##                  Fixed Effects Slope    2.5%          97.5%
##tAI               27.91440               20.25226      35.55821     (R slope was not scaled)
##onefour          -0.03082               -0.4059275     0.3456071    (R slope was not scaled)
##tAI:onefour       0.09169               -0.7810869     0.9611812    (R slope was not scaled)

#sqrtPPR ~ tAI + onehalf + tAI:onehalf + (0+tAI|gene) + (0+onehalf|gene) + (1|gene)  #1313 genes (1/2 CDS)
##
##                  Fixed Effects Slope    2.5%          97.5%
##tAI               27.77645               20.1753       35.36042     (R slope was not scaled)
##onehalf          -0.01897               -0.4140938     0.3717285    (R slope was not scaled)
##tAI:onehalf       0.14477               -0.7617331     1.062494     (R slope was not scaled)

#sqrtPPR ~ tAI + threefour + tAI:threefour + (0+tAI|gene) + (0+threefour|gene) + (1|gene)  #1314 genes (3/4 CDS)
##
##                  Fixed Effects Slope    2.5%          97.5%
##tAI               25.44731               17.9043       32.97481     (R slope was not scaled)
##threefour        -0.02733               -0.3619499     0.3032954    (R slope was not scaled)
##tAI:threefour     0.15814               -0.608933      0.9356227    (R slope was not scaled)

cds25.int.deltaG.tAI.data.sqrt <- data.frame(
  term = str_wrap(c("tAI","CDS25-deltaG","tAI:CDS25-deltaG"), width = 12),
  slope = c(27.91440,-0.03082,0.09169),
  down95ci = c(20.25226,-0.4059275,-0.7810869),
  up95ci = c(35.55821,0.3456071,0.9611812)
)
cds25.int.deltaG.tAI.data.sqrt$term <- factor(cds25.int.deltaG.tAI.data.sqrt$term, levels = str_wrap(c("tAI","CDS25-deltaG","tAI:CDS25-deltaG"), width = 12), ordered = TRUE)
str(cds25.int.deltaG.tAI.data.sqrt)

cds50.int.deltaG.tAI.data.sqrt <- data.frame(
  term = str_wrap(c("tAI","CDS50-deltaG","tAI:CDS50-deltaG"), width = 12),
  slope = c(27.77645,-0.01897,0.14477),
  down95ci = c(20.1753,-0.4140938,-0.7617331),
  up95ci = c(35.36042,0.3717285,1.062494)
)
cds50.int.deltaG.tAI.data.sqrt$term <- factor(cds50.int.deltaG.tAI.data.sqrt$term, levels = str_wrap(c("tAI","CDS50-deltaG","tAI:CDS50-deltaG"), width = 12), ordered = TRUE)
str(cds50.int.deltaG.tAI.data.sqrt)
cds75.int.deltaG.tAI.data.sqrt <- data.frame(
  term = str_wrap(c("tAI","CDS75-deltaG","tAI:CDS75-deltaG"), width = 12),
  slope = c(25.44731,-0.02733,0.15814),
  down95ci = c(17.9043,-0.3619499,-0.608933),
  up95ci = c(32.97481,0.3032954,0.9356227)
)
cds75.int.deltaG.tAI.data.sqrt$term <- factor(cds75.int.deltaG.tAI.data.sqrt$term, levels = str_wrap(c("tAI","CDS75-deltaG","tAI:CDS75-deltaG"), width = 12), ordered = TRUE)
str(cds75.int.deltaG.tAI.data.sqrt)



colorBlindGrey8 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
tiff(file="figureS4.tiff", width = 10, height = 20, units = 'cm', res = 300)
pd <- position_dodge(1)
figS6a <- ggplot(local.deltaG.tAI.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-1,2)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7], colorBlindGrey8[7], colorBlindGrey8[7], colorBlindGrey8[7], colorBlindGrey8[1], colorBlindGrey8[1], colorBlindGrey8[1], colorBlindGrey8[1]))
figS6b <- ggplot(cds25.int.deltaG.tAI.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-0.5,2)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))
figS6c <- ggplot(cds50.int.deltaG.tAI.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-0.5,2)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))
figS6d <- ggplot(cds75.int.deltaG.tAI.data, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-0.5,2)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))

figS6e <- ggplot(local.deltaG.tAI.data.sqrt, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-5,40)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7], colorBlindGrey8[7], colorBlindGrey8[7], colorBlindGrey8[7], colorBlindGrey8[1], colorBlindGrey8[1], colorBlindGrey8[1], colorBlindGrey8[1]))
figS6f <- ggplot(cds25.int.deltaG.tAI.data.sqrt, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-1,40)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))
figS6g <- ggplot(cds50.int.deltaG.tAI.data.sqrt, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-1,40)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))
figS6h <- ggplot(cds75.int.deltaG.tAI.data.sqrt, aes(x = term, y = slope, fill = term)) + 
  geom_bar(stat="identity", position=pd) +
  geom_errorbar(aes(ymin=down95ci,ymax=up95ci),width=0.2,position=pd) +
  xlab("")+
  ylab(y_title) +
  ylim(-1,40)+
  theme_bw() +
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  theme(axis.ticks = element_line(color="black"))+
  theme(legend.position="none")+
  scale_fill_manual(values=c(colorBlindGrey8[6], colorBlindGrey8[7],colorBlindGrey8[1]))

r1 <- plot_grid(figS6a, labels = c('A'), label_size = 12, nrow = 1, hjust=-2, vjust=1)
r2 <- plot_grid(figS6b, figS6c, figS6d, labels = c('B','C','D'), label_size = 12, nrow = 1, hjust=-2, vjust=0)
r3<- plot_grid(figS6e, labels = c('E'), label_size = 12, nrow = 1, hjust=-2, vjust=1)
r4 <- plot_grid(figS6f, figS6g, figS6h, labels = c('F','G','H'), label_size = 12, nrow = 1, hjust=-2, vjust=0)
plot_grid(r1, r2, r3, r4, nrow = 4)
dev.off()



