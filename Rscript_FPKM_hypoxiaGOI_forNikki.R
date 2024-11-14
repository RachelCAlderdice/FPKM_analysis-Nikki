#FPKM analysis of hypoxia associated GOI, for Nikki's PhD Chapter Apr 2024

#load dependencies
library("ggplot2")

#load input files
fpkm.pre<-read.csv("./input_files/hawaii_fpkm_matrix_withfunctions_RA.csv", sep = ",", header = TRUE, row.names = 1)
map.pre<-read.csv("./input_files/metadata_RA.csv", sep = ",")

#format files
row.names(map.pre) = map.pre$Sample
fpkm.sor=fpkm.pre[ , order(names(fpkm.pre))]
map=map.pre[order(rownames(map.pre)), ]
rownames(map)=gsub("\\-", "_", rownames(map))

#confirm names in files match
all(colnames(fpkm.sor) %in% row.names(map))
all(colnames(fpkm.sor) == row.names(map))

#subset data by treatment groups
map.A=subset(map, map$Temperature == "30" & map$Treatment == "mangrove_reef")
map.B=subset(map, map$Temperature == "30" & map$Treatment == "reef_reef")
map.C=subset(map, map$Temperature == "30" & map$Treatment == "wild_mangrove")
map.D=subset(map, map$Temperature == "30" & map$Treatment == "wild_reef")

#select genes of interest
gene=read.delim("./input_files/hypoxia_GOI.txt", header=TRUE)
genes=subset(gene, gene$anno == "HIF1A") #change gene name accordingly

#determine mean and SE of FPKM count data for each treatment group
# (those transcripts with the same gene annotations are summed together)

#30 mangrove_reef
fpkm=fpkm.sor[,(names(fpkm.sor)[(names(fpkm.sor) %in% row.names(map.A))])]
fpkm.a=subset(fpkm, row.names(fpkm) %in% genes$GeneID)
fpkm.colsum <- apply(fpkm.a, 2, sum)
Mean <- mean(fpkm.colsum)
x <- sd(fpkm.colsum)
SE <- x/sqrt(ncol(fpkm.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("30",length(nrow(out)))
out$Treatment = rep("mangrove_reef",length(nrow(out)))
mr_30 <- out

#30 reef_reef
fpkm=fpkm.sor[,(names(fpkm.sor)[(names(fpkm.sor) %in% row.names(map.B))])]
fpkm.a=subset(fpkm, row.names(fpkm) %in% genes$GeneID)
fpkm.colsum <- apply(fpkm.a, 2, sum)
Mean <- mean(fpkm.colsum)
x <- sd(fpkm.colsum)
SE <- x/sqrt(ncol(fpkm.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("30",length(nrow(out)))
out$Treatment = rep("reef_reef",length(nrow(out)))
rr_30 <- out

#30 wildmangrove
fpkm=fpkm.sor[,(names(fpkm.sor)[(names(fpkm.sor) %in% row.names(map.C))])]
fpkm.a=subset(fpkm, row.names(fpkm) %in% genes$GeneID)
fpkm.colsum <- apply(fpkm.a, 2, sum)
Mean <- mean(fpkm.colsum)
x <- sd(fpkm.colsum)
SE <- x/sqrt(ncol(fpkm.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("30",length(nrow(out)))
out$Treatment = rep("wildmangrove",length(nrow(out)))
wm_30 <- out

#30 wildreef
fpkm=fpkm.sor[,(names(fpkm.sor)[(names(fpkm.sor) %in% row.names(map.D))])]
fpkm.a=subset(fpkm, row.names(fpkm) %in% genes$GeneID)
fpkm.colsum <- apply(fpkm.a, 2, sum)
Mean <- mean(fpkm.colsum)
x <- sd(fpkm.colsum)
SE <- x/sqrt(ncol(fpkm.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("30",length(nrow(out)))
out$Treatment = rep("wildreef",length(nrow(out)))
wr_30 <- out

#combine each treatment group results
out=rbind(mr_30, rr_30, wm_30, wr_30)

#add gene name label 
out$Type = rep("HIF1A",length(nrow(out))) #change gene name accordingly

#create table of results
write.table(out, 
            "./results/HIF1A_fpkm_temp30.txt", #change gene name accordingly
            sep = "\t", 
            quote = FALSE, 
            row.names= FALSE, 
            col.names = TRUE) 

#create scatterplot of fpkm expression pattern over timepoints for each GOI

#load input files
data <- read.delim("./results/HIF1A_fpkm_temp30.txt") #change file accordingly for different genes

#format data
data$Treatment <- as.character(data$Treatment)
data$Treatment <- factor(data$Treatment, levels=unique(data$Treatment))

#custom labels for plot legend
custom_labels <- c("mangrove_reef" = "M-R",
                   "reef_reef" = "R-R",
                   "wildmangrove" = "WM",
                   "wildreef" = "WR")

#plot data
#HIF1A gene as example here, change scale accordingly 
a<-ggplot(data, aes(x = Treatment, y = Mean, color = Treatment, group = Treatment)) + 
  ylim(0, 2200) + #change according to range in 'data' dataframe
  geom_bar(aes(fill = Treatment), colour= "black", stat = "identity") + 
  ylab(label="FPKMs") + 
  xlab("Treatment") + 
  scale_x_discrete(labels = custom_labels) + 
  theme_classic(base_size = 29) + 
  theme(aspect.ratio = 1, axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2, color= "black") + 
  scale_color_manual(values=c("#9933CC","orange", "#CCCCFF", "#FFFF99")) + 
  scale_fill_manual(values=c("#9933CC","orange", "#CCCCFF", "#FFFF99")) 

#add gene name label to plot
grob <- grobTree(textGrob("HIF1A", #change gene name accordingly
                          x=0.1,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=22, fontface="italic")))
#print plot
a + annotation_custom(grob)
#save plot
ggsave("./results/HIF1A_fpkms_barplot.pdf", width=9, height=10, dpi= 300)
#repeat for each GOI, then combine and finalize format of pdf plots in affinity designer
