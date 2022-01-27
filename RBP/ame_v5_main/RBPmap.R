## test RBmap

tmp_5utr <- read.table('Documents/UCLA/Ribo_New/RBP/ame_v5_main/summary_RBPmap_5utr.txt')
library(ggplot2)

tmp_5utr$V2 <- -tmp_5utr$V2

for (r in unique(tmp_5utr$V3)){
  for (g in unique(tmp_5utr$V1)){
    if (nrow(tmp_5utr[tmp_5utr$V1 == g & tmp_5utr$V3 == r,]) == 0){
      tmp_5utr<- rbind(tmp_5utr, data.frame(V1=g,V2=unique(tmp_5utr$V2[tmp_5utr$V1 == g]),V3=r,V4=0))
    }
  }
}


cc <- c()
pval <- c()
for (r in unique(tmp_5utr$V3)){
  tmp <- cor.test(tmp_5utr[tmp_5utr$V3==r,"V2"], tmp_5utr[tmp_5utr$V3==r,"V4"], method = "spearman")
  cc <- c(cc, tmp$estimate)
  pval <- c(pval, tmp$p.value)
}

names(cc) <- names(pval) <- unique(tmp_5utr$V3)

tmp <- cbind(cc,pval)
tmp[order(tmp[,"pval"]),]

p <- ggplot(tmp_5utr[tmp_5utr$V3 %in% row.names(tmp)[tmp[,"pval"] <= 0.05],], aes(x=V2,y=V4)) + geom_point() + facet_wrap(~V3, scales = "free")
p

tmp_3utr <- read.table('Documents/UCLA/Ribo_New/RBP/ame_v5_main/summary_RBPmap_3utr.txt')
library(ggplot2)

tmp_3utr$V2 <- -tmp_3utr$V2

for (r in unique(tmp_3utr$V3)){
  for (g in unique(tmp_3utr$V1)){
    if (nrow(tmp_5utr[tmp_3utr$V1 == g & tmp_3utr$V3 == r,]) == 0){
      tmp_5utr<- rbind(tmp_3utr, data.frame(V1=g,V2=unique(tmp_3utr$V2[tmp_3utr$V1 == g]),V3=r,V4=0))
    }
  }
}

cc <- c()
pval <- c()
for (r in unique(tmp_3utr$V3)){
  tmp <- cor.test(tmp_3utr[tmp_3utr$V3==r,"V2"], tmp_3utr[tmp_3utr$V3==r,"V4"], method = "spearman")
  cc <- c(cc, tmp$estimate)
  pval <- c(pval, tmp$p.value)
}

names(cc) <- names(pval) <- unique(tmp_3utr$V3)

tmp <- cbind(cc,pval)
tmp <- tmp[order(tmp[,"pval"]),]
tmp <- cbind(tmp, adjpval=p.adjust(tmp[,"pval"],"bonferroni"))


p <- ggplot(tmp_3utr[tmp_3utr$V3 %in% row.names(tmp)[tmp[,"pval"] <= 0.05],], aes(x=V2,y=V4)) + geom_point() 
p <- p + geom_smooth(method = "lm") 
p <- p + facet_wrap(~V3, scales = "free")
p
