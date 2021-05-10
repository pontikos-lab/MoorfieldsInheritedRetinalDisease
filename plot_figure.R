library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(cowplot)


d <- read('Supp_Table_1.csv')
names(d)[1]<-'Gene'


d2 <- read('genes_transcripts_mart_export.txt')
d2$Gene <- d2[,'Gene name']
merge(d,d2,by='Gene')->d3
X <- as.data.frame( do.call('rbind', by(d3, d3$Gene, function (x) {
gene.name <- x[,'Gene']
transcript.length <- max(x[ 'Transcript length (including UTRs and CDS)'])
unique(data.frame(gene.name=gene.name, transcript.length=transcript.length))
})))
rownames(X) <- X$gene.name
X <- X[d$Gene,]


d$MOI <- d[, "Possible modes of inheritance"]
d$families.affected.number <- d[,'Families affected (number)']
d$individuals.affected.number <- d[,'Individuals affected (number)']
d$longest.transcript.length <- d[, 'Length of longest transcript']
d$MOI2 <- d$MOI
d[which(d$MOI=='Dominant and Recessive'),'MOI2']<-'Recessive/Dominant'
d[which(d$MOI=='Recessive and Dominant'),'MOI2']<-'Recessive/Dominant'
d <- d[-which(d$MOI==''),]

plot.moi <- function (x) {
data <- subset(d,MOI2==x)
ggplot(data = data, aes(x=longest.transcript.length,y=families.affected.number,label=Gene)) +
scale_x_log10( breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
theme_bw() +
#geom_text_repel(aes(label=Gene),size=2, arrow=arrow(length=unit(0.03,'npc'),type='closed',ends='first'))+
geom_text_repel(aes(label=Gene),size=2, segment.size=0.2, segment.color="grey50")+
geom_point(color=ifelse(data$longest.transcript.length<36000,'green','red'))+
labs(title = "", subtitle = "", x = expression(Log[10]*" Length of Longest Transcript"), y = expression(Log[10]*" Number of Affected Families"))+
geom_smooth(method="lm", size=0.25) +
#facet_wrap(~ MOI2, scales="free") +
stat_cor(method = "pearson", label.x.npc='left', label.y=2.5)
}
p1 <- plot.moi('Recessive')
p2 <- plot.moi('X-linked')
p3 <- plot.moi('Dominant')
p4 <- plot.moi('Mitochondrial inheritance')
plot_grid(p1,p2,p3,p4,labels=c('A','B','C','D'),label_size=12)


d$Gene = with(d,reorder(Gene,families.affected.number,median))

ggplot(d,aes(x=Gene,y=families.affected.number))+geom_bar(stat='identity')



ggplot(data = d, aes(x=individuals.affected.number,y=families.affected.number,label=Gene)) +
#scale_y_continuous(trans='log10') +
#annotation_logticks() +
#scale_x_continuous(trans='log10') +
#scale_x_log10( breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) +
#stat_cor(method = "pearson", label.x = 3, label.y = 30)+
theme_bw() +
#geom_text_repel(aes(label=Gene),size=2, arrow=arrow(length=unit(0.03,'npc'),type='closed',ends='first'))+
geom_text_repel(aes(label=Gene),size=2, segment.size=0.2, segment.color="grey50")+
geom_point(color='red')+
#labs(title = "", subtitle = "", x = expression(Log[10]*" Length of Longest Transcript"), y = expression(Log[10]*" Number of Affected Families"))+
geom_smooth(method="lm", size=0.25) +
facet_wrap(~ MOI2, scales="free")
