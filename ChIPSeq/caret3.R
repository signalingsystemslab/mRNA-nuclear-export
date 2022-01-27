library(data.table)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(leaps)
library(gridExtra)
library(caret)
library(skimr)
library(RANN)
library(doSNOW)
library(plyr)
library(grid)

#Converts ensemble ID, prepares data set
processdata <- function(counts, annotation) {
    colnames(annotation)[1] <- 'peakid'
    setorder(annotation, peakid)
    setorder(counts, V4)
    gene <- annotation[, c(16,3,4)]
    allcount <- data.frame(counts[, 11:12])
    #Normalize by peak depth
    factor <- apply(allcount, 2, sum) / min(apply(allcount, 2, sum))
    for (i in 1:2) {
        allcount[, i] <- allcount[, i] / factor[i]
    }
    table <- data.frame(gene, log(allcount + 1, base = 2))
    colnames(table) <- c('gene', 'start', 'end', 'rep1', 'rep2')
    return(table)
}

#Load histone data Sets
H3K27ac <- processdata(fread("finalnaive/multicov_H3K27ac.bed"), fread("finalnaive/H3K27acannotate.txt"))
H3K36me3 <- processdata(fread("finalnaive/multicov_H3K36me3.bed"), fread("finalnaive/H3K36me3annotate.txt"))
H3K4me3 <- processdata(fread("finalnaive/multicov_H3K4me3.bed"), fread("finalnaive/H3K4me3annotate.txt"))
H3K79me2 <- processdata(fread("finalnaive/multicov_H3K79me2.bed"), fread("finalnaive/H3K79me2annotate.txt"))

#mRNA paramater load
param <- fread('parameter.csv')[, -1]
colnames(param) <- c('name', 'k1prime', 'k2', 'k2prime', 'kdeg', 'transport')
curated <- fread('curated.csv')
splicing <- readRDS('splicingprob.rds')
curated <- curated[match(splicing$Name, curated$Name),]
param <- param[match(splicing$Name, param$name),]

#gets gene length
length <- as.numeric(apply(curated[, 7], 1, function(x) gsub(',.*', '', x))) - as.numeric(apply(curated[, 6], 1, function(x) gsub(',.*', '', x)))
#strand information
strand <- apply(curated[, 4], 1, function(x) {
    if (x == '+')
        return(1)
    else
        return(-1)
    })
length <- length * strand
curated$length <- length
#append features to curated
curated$intron <- curated$`Nb Exons (Max)` - 1
curated$prob <- splicing[match(curated$Name, splicing$Name), 5]
curated$kdeg <- param[match(curated$Name, param$name)]$kdeg
#get diane tss
tss <- as.numeric(apply(curated[, 6], 1, function(x) gsub(',.*', '', x)))
filttss <- cbind(curated[, c(2:4)], tss, curated[, c(9:12)])
filttss <- filttss[!is.na(filttss$tss)]

H3K27ac <- H3K27ac %>% filter(gene %in% param$name) %>% arrange(gene)
H3K36me3 <- H3K36me3 %>% filter(gene %in% param$name) %>% arrange(gene)
H3K4me3 <- H3K4me3 %>% filter(gene %in% param$name) %>% arrange(gene)
H3K79me2 <- H3K79me2 %>% filter(gene %in% param$name) %>% arrange(gene)

#combine histone replicates
agg <- function(data) {
    tssmatch <- filttss[match(unlist(data$gene), unlist(filttss$Name)),]$tss
    return(data.frame('gene' = data$gene,
                         'distance' = data$start-tssmatch,
                         'width' = data$end - data$start,
                         'naive' = apply(data[, 4:5], 1, mean)))
}

aggH3K27ac <- agg(H3K27ac)
aggH3K36me3 <- agg(H3K36me3)
aggH3K4me3 <- agg(H3K4me3)
aggH3K79me2 <- agg(H3K79me2)

aggH3K27ac <- aggH3K27ac[!is.na(aggH3K27ac$distance),]
aggH3K36me3 <- aggH3K36me3[!is.na(aggH3K36me3$distance),]
aggH3K4me3 <- aggH3K4me3[!is.na(aggH3K4me3$distance),]
aggH3K79me2 <- aggH3K79me2[!is.na(aggH3K79me2$distance),]

#Sliding window by distance from TSS
slide <- function(data, start, end, bins) {
    startpoints <- seq(start, end, length.out=bins)
    df <- data.frame(matrix(nrow = length(unique(data$gene)), ncol = 0))
    for (point in c(1:(length(startpoints)-1))) {
        bool <- data$distance >= startpoints[point] & data$distance <= startpoints[point+1]
        temp <- data
        temp[!bool,]$naive <- 0
        temp$gene <- as.character(temp$gene)
        temp <- temp %>% group_by(gene) %>% dplyr::summarise(naive = list(naive))
        temp$naive <- apply(temp, 1, function(x) {
            y <- unlist(x[2])
            return(mean(y[y != 0]))
        })
        temp[is.na(temp[, 2]), 2] <- 0
        df <- cbind(df, temp$naive)
    }
    temp <- data %>% group_by(gene) %>% dplyr::summarise(naive = mean(naive))
    df$total <- unlist(temp[, 2])
    temp <- data %>% group_by(gene) %>% dplyr::summarise(pos = mean(naive[distance > 0]), neg = mean(naive[distance < 0]))
    temp[is.na(temp)] <- 0
    df$down <- unlist(temp[, 2])
    df$up <- unlist(temp[, 3])
    df$gene <- temp$gene
    return(df)
}

#Combines all the bins into input matrix
#enhancer
H3K27acbin <- slide(aggH3K27ac, -15000, 15000, 5)
#active
H3K36me3bin <- slide(aggH3K36me3, -25000, 25000, 9)
#promoter
H3K4me3bin <- slide(aggH3K4me3, -5000, 5000, 5)
#active
H3K79me2bin <- slide(aggH3K79me2, -25000, 25000, 7)

#Get common genes and combine features
common <- Reduce(union, list(H3K27acbin$gene, H3K36me3bin$gene, H3K4me3bin$gene, H3K79me2bin$gene))
combine <- data.frame(H3K27ac = H3K27acbin[match(common, H3K27acbin$gene), - ncol(H3K27acbin)],
                      H3K36me3 = H3K36me3bin[match(common, H3K36me3bin$gene), - ncol(H3K36me3bin)],
                      H3K4me3 = H3K4me3bin[match(common, H3K4me3bin$gene), - ncol(H3K4me3bin)],
                      H3K79me2 = H3K79me2bin[match(common, H3K79me2bin$gene), - ncol(H3K79me2bin)])
combine <- cbind(combine, filttss[match(common, filttss$Name), 5:7])
combine[is.na(combine)] <- 0

set.seed(123)
#Series of feature combinations
runindex <- list(c(2:35), c(36), c(37), c(38), c(36:38),c(2:38))
count <- 0

for (i in runindex[3]) {
    ### Machine Learning for k1prime
    combine2 <- cbind(param[match(common, param$name), 6], combine)
    colnames(combine2) <- c('param', paste("H3K27ac", 1:(ncol(H3K27acbin) - 4), sep = "_"), paste("H3K27ac", c('total', 'down', 'up'), sep = "_"),
                                     paste("H3K36me3", 1:(ncol(H3K36me3bin) - 4), sep = "_"), paste("H3K36me3", c('total', 'down', 'up'), sep = "_"),
                                     paste("H3K4me3", 1:(ncol(H3K4me3bin) - 4), sep = "_"), paste("H3K4me3", c('total', 'down', 'up'), sep = "_"),
                                     paste("H3K79me2", 1:(ncol(H3K79me2bin) - 4), sep = "_"), paste("H3K79me2", c('total', 'down', 'up'), sep = "_"),
                                     'length', 'intron', 'prob')

    # Step 1: Get row numbers for the training data
    trainRowNumbers <- createDataPartition(combine2$param, p = 0.8, list = FALSE)
    trainData <- combine2[trainRowNumbers,]

    preProcess_range_model <- preProcess(trainData[, 2:(ncol(trainData) - 3)], method='range')
    trainData[, 2:(ncol(trainData) - 3)] <- predict(preProcess_range_model, newdata = trainData[, 2:(ncol(trainData) - 3)])
    ind <- unlist(i)
    train_x <- trainData[, ..ind]
    train_y <- trainData$param

    train.control <- trainControl(method = 'repeatedcv',
                                    number = 5,
                                    repeats = 3,
                                    search = 'grid')

    cl <- makeCluster(2, type = 'SOCK')
    registerDoSNOW(cl)

    #Number of Iterations and the Learning Rate
    nrounds <- 500
    tune_grid <- expand.grid(
        nrounds = seq(from = 200, to = nrounds, by = 100),
        eta = c(0.025, 0.05, 0.1, 0.3),
        max_depth = c(2, 3, 4, 5),
        gamma = 0,
        colsample_bytree = 0.1,
        min_child_weight = 1,
        subsample = 1
    )

    #Number of Iterations and the Learning Rate
    xgb_tune <- train(x = train_x,
                        y = train_y,
                        method = 'xgbTree',
                        tuneGrid = tune_grid,
                        trControl = train.control,
                        verbose = TRUE)

    print("1 done")
    stopCluster(cl)
    cl <- makeCluster(2, type = 'SOCK')
    registerDoSNOW(cl)

    #Maximum Depth and Minimum Child Weight
    tune_grid2 <- expand.grid(
        nrounds = seq(from = 50, to = nrounds, by = 100),
        eta = xgb_tune$bestTune$eta,
        max_depth = ((xgb_tune$bestTune$max_depth - 1):(xgb_tune$bestTune$max_depth + 1)),
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = c(1, 2, 3),
        subsample = 1
    )

    xgb_tune2 <- train(x = train_x,
                        y = train_y,
                        method = 'xgbTree',
                        trControl = train.control,
                        verbose = TRUE)
    print("2 done")

    stopCluster(cl)
    cl <- makeCluster(2, type = 'SOCK')
    registerDoSNOW(cl)

    #Column and Row Sampling
    tune_grid3 <- expand.grid(
        nrounds = seq(from = 50, to = nrounds, by = 100),
        eta = xgb_tune$bestTune$eta,
        max_depth = xgb_tune2$bestTune$max_depth,
        gamma = 0,
        colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
        min_child_weight = xgb_tune2$bestTune$min_child_weight,
        subsample = c(0.5, 0.75, 1.0)
    )

    xgb_tune3 <- train(x = train_x,
                        y = train_y,
                        method = 'xgbTree',
                        tuneGrid = tune_grid3,
                        trControl = train.control,
                        verbose = TRUE)
    print("3 done")
    stopCluster(cl)
    cl <- makeCluster(2, type = 'SOCK')
    registerDoSNOW(cl)

    #Gamma
    tune_grid4 <- expand.grid(
      nrounds = seq(from = 50, to = nrounds, by = 100),
      eta = xgb_tune$bestTune$eta,
      max_depth = xgb_tune2$bestTune$max_depth,
      gamma = c(0.5, 1.5, 2, 3, 4),
      colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
      min_child_weight = xgb_tune2$bestTune$min_child_weight,
      subsample = xgb_tune3$bestTune$subsample
    )

    xgb_tune4 <- train(x = train_x,
                        y = train_y,
                        method = 'xgbTree',
                        tuneGrid = tune_grid4,
                        trControl = train.control,
                        verbose = TRUE)
    print("4 done")
    stopCluster(cl)
    cl <- makeCluster(2, type = 'SOCK')
    registerDoSNOW(cl)

    #Reducing learning rate
    tune_grid5 <- expand.grid(
        nrounds = seq(from = 100, to = 2500, by = 100),
        eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
        max_depth = xgb_tune2$bestTune$max_depth,
        gamma = xgb_tune4$bestTune$gamma,
        colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
        min_child_weight = xgb_tune2$bestTune$min_child_weight,
        subsample = xgb_tune3$bestTune$subsample
    )

    xgb_tune5 <- train(x = train_x,
                        y = train_y,
                        method = 'xgbTree',
                        tuneGrid = tune_grid5,
                        trControl = train.control,
                        verbose = TRUE)
    print("5 done")
    stopCluster(cl)
    cl <- makeCluster(2, type = 'SOCK')
    registerDoSNOW(cl)

    #Final Model
    final_grid <- expand.grid(
        nrounds = xgb_tune5$bestTune$nrounds,
        eta = xgb_tune5$bestTune$eta,
        max_depth = xgb_tune5$bestTune$max_depth,
        gamma = xgb_tune5$bestTune$gamma,
        colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
        min_child_weight = xgb_tune5$bestTune$min_child_weight,
        subsample = xgb_tune5$bestTune$subsample
    )

    xgb_model <- train(x = train_x,
                        y = train_y,
                        method = 'xgbTree',
                        tuneGrid = final_grid,
                        trControl = train.control,
                        verbose = TRUE)

    stopCluster(cl)

    assign(paste('trainmodel', count, sep = ''), xgb_model)
    count <- count+1
}

chip <- trainmodel0
length <- trainmodel1
intron <- trainmodel2
prob <- trainmodel3
feat <- trainmodel4
chipfeat <- trainmodel5

#Summarize resamples 
results <- resamples(list('chip' = chip, 'length' = length, 'intron' = intron, 'prob' = prob, 'feat' = feat, 'chip+feat' = chipfeat))
summary(results)
bwplot(results)
dotplot(results)

#Resample visualizations
resamples <- results$values[, grep('Rsquared', colnames(results$values))]
colnames(resamples) <- c('chip', 'length', 'intron', 'prob', 'feat', 'chipfeat')
resamples2 <- data.frame('rsquared' = c(resamples[, 1], resamples[, 2], resamples[, 3], resamples[, 4], resamples[, 5], resamples[, 6]), 'model' = rep(colnames(resamples), each = 30))

resamples3 <- data.frame('model' = colnames(resamples), 'mean' = apply(resamples, 2, mean), 'sd' = apply(resamples, 2, sd))
resamples3$model <- factor(resamples3$model, levels = resamples3$model)

#errorbar plot
points <- ggplot(resamples3, aes(x = model, y = mean)) + geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean - sd,
                      ymax = mean + sd),
                  width = 0.5, size = 0.5) + theme_bw() + labs(x = 'Models', y = bquote(R ^ 2)) + theme(
       axis.text = element_text(size = 14),
       axis.title = element_text(size = 16),
       ) + coord_flip() + ylim(-0.02, 1)

#Important variables
data.frame(varImp(chip, scale = F)[1])
features <- c('chip', 'length', 'prob', 'intron')
models <- c('chip', 'length','prob','intron','feat', 'chipfeat')
aggmap <- data.frame(matrix(nrow = length(models), ncol = length(features)))
rownames(aggmap) <- models
colnames(aggmap) <- features
aggmap[1,] <- c(1, NA,NA,NA,NA)
aggmap[2,] <- c(NA, 1, NA,NA,NA)
aggmap[3,] <- c(NA, NA, 1, NA, NA)
aggmap[4,] <- c(NA, NA, NA, 1, NA)
aggmap[5,] <- c(NA, unlist(data.frame(varImp(feat, scale = F)[1])))
aggmap[6,] <- c(1 - (0.281490 + 0.2450728 + 0.115440), 0.281490, 0.245072, 0.115440)
aggmap <- aggmap[,4:1]
aggmap <- t(aggmap)

#Color gradient for heatmap
distmat <- dist(t(aggmap))
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette) {
    stopifnot(length(colors) == 4)
    ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
    ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
    return(c(ramp1, ramp2))
}

#feature heatmap
svg("heatmap.svg")
paletteLength <- 50
myColor <- colorRampPalette(c("white", "darkgreen"))(paletteLength)
myBreaks <- seq(0, 1, 1 / 50)
important <- pheatmap(aggmap, cluster_cols = F, cluster_rows = F, legend = TRUE, main = 'Feature Importance for Models', color = myColor, fontsize = 10, cellheight = 40, cellwidth = 40, breaks = myBreaks)

grid.text(23, x = 0.565, y = 0.16)
grid.text(7, x = 0.565, y = 0.555)
grid.text(14, x = 0.565, y = 0.715)
grid.text(9, x = 0.565, y = 0.875)

dev.off()
ggsave(file = "points.svg", plot = points, width = 4, height = 4)
ggsave(file = 'box.svg', plot = box, width = 4, height = 4)

#P-values mann whitney
unique <- as.character(resamples3$model)
for (i in unique[1:length(unique)]) {
    temp <- wilcox.test(resamples2[which(resamples2$model == unique[4]), 1], resamples2[which(resamples2$model == i), 1])
    print(paste(unique[4], i))
    print(temp$p.value)

}

#Figure 1
filtparam <- param[match(common, param$name),]
totalind <- grep('total', names(combine))
downind <- grep('down', names(combine))
upind <- grep('up', names(combine))
totalmark <- data.frame(common, filtparam$transport, combine[, totalind])
downmark <- data.frame(common, filtparam$transport, combine[, downind])
upmark <- data.frame(common, filtparam$transport, combine[, upind])
colnames(totalmark) <- c('gene', 'transport', 'H3K27ac', 'H3K36me3', 'H3K4me3', 'H3K79me2')
colnames(downmark) <- c('gene', 'transport', 'H3K27ac', 'H3K36me3', 'H3K4me3', 'H3K79me2')
colnames(upmark) <- c('gene', 'transport', 'H3K27ac', 'H3K36me3', 'H3K4me3', 'H3K79me2')

plotmark <- function(x, mark) {
    ggplot(x, aes(x = log(mark + 1, base = 2), y = transport)) + geom_point(alpha = 0.6, size = 2) + theme_bw() + labs(x = expression('ChIP-seq Signal' ~ (Log[2] + 1)), y = 'Transport Rate') + theme(axis.title = element_text(size = 14), legend.text = element_text(size = 12), legend.title = element_text(size = 14))
}
svg('totalsignal.svg')
g1 <- plotmark(totalmark, totalmark$H3K27ac) + labs(title = 'H3K27ac Total ChIP-seq Signal') + stat_cor(label.x = 0.5, label.y = 0.25)
g2 <- plotmark(totalmark, totalmark$H3K36me3) + labs(title = 'H3K36me3 Total ChIP-seq Signal') + stat_cor(label.x = 0.5, label.y = 0.25)
g3 <- plotmark(totalmark, totalmark$H3K4me3) + labs(title = 'H3K4me3 Total ChIP-seq Signal') + stat_cor(label.x = 0.5, label.y = 0.25)
g4 <- plotmark(totalmark, totalmark$H3K79me2) + labs(title = 'H3K79me2 Total ChIP-seq Signal') + stat_cor(label.x = 0.5, label.y = 0.25)
grid.arrange(g1, g2, g3, g4, nrow = 2)
dev.off()

chipimp <- data.frame(rownames(as.data.frame(varImp(chip)[1])), unlist(varImp(chip, scale = F)[1]))
colnames(chipimp) <- c('feature', 'importance')
chipimp$feature <- factor(chipimp$feature, levels = chipimp$feature[28:1])
mark <- c()
for (i in 1:28) {
    if (length(grep('H3K27ac', chipimp[i, 1]))>0) {
        mark <- c(mark,'H3K27ac')
    } else if (length(grep('H3K4me3', chipimp[i, 1])) > 0) {
        mark <- c(mark, 'H3K4me3')
    } else if (length(grep('H3K79me2', chipimp[i, 1])) > 0) {
        mark <- c(mark, 'H3K79me2')
    } else if (length(grep('H3K36me3', chipimp[i, 1])) > 0) {
        mark <- c(mark, 'H3K36me3')
    }
}
chipimp$mark <- mark
ggplot(chipimp[1:23,], aes(x = feature, y = importance)) + geom_point(stat = 'identity', aes(col = mark), size = 6) + coord_flip() + theme_bw()+
    theme(axis.title.y = element_blank(), axis.text = element_text(size = 16), axis.title = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 16))