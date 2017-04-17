setwd("~/Documents/IRTData")
library('mirt')
library('RCurl')
library('MASS')
library('pwr')
diag.df <- read.table('diag.txt', header = T)
diag1 <- diag.df[,1:2]
sort(table(diag1[,2]),decreasing = T)
top3 <- c('DSM4:V71.09', 'DSM4:295.XX','DSM4:296.70A')
diag.top3 <- diag1[which(diag1[,2] %in% top3),]#10545
halluc.vars <- read.table("halluc_vars.txt")
#write.table(halluc.vars, file = "halluc_vars.txt")
everything.df <- read.table("everything_201703170901.tsv", sep ='\t', stringsAsFactors = FALSE, header = TRUE, na.strings = 
                              c("9", "u", "U", "","-1", "-9999", "-999999", "uuu", "99","999", "9999","998","9998", ".", NA), fill = TRUE, quote = "")
dim(everything.df)#48367  4414
names(everything.df)[1:2] <- c("ind_id", "table_set_id")

dup.rows <- everything.df[duplicated(everything.df[,1]),1]#duplicated rows
everything.df <-everything.df[!everything.df[,1]%in% dup.rows,]
dimnames(everything.df)[[1]] <- everything.df[,1]
match(tolower(halluc.vars[,1]), tolower(colnames(everything.df)))
halluc.df <- everything.df[,match(tolower(halluc.vars[,1]), tolower(colnames(everything.df)))[1:82]]
halluc.mat <- cleaning(halluc.df, rows = T)
write.table(halluc.mat, file = 'halluc_mat.txt')
dim(halluc.mat) # 24794    82
diag.top3 <- diag.top3[which(diag.top3[,1] %in% rownames(halluc.mat)),]
dim(diag.top3) #7966    2
halluc.mat <- halluc.mat[which(rownames(halluc.mat) %in% diag.top3[,1]),]
dim(halluc.mat)
halluc.mat <- cleaning(halluc.mat, rows = F)
dim(halluc.mat) #7966   82


group0 <- as.factor(diag.top3[,2])
group0 <- droplevels(group0)
length(group0)

dif.mirt1 <- mirt(halluc.mat, 1, technical = list(NCYCLES = 5000))#1 factor
dif.mirt1.fscores <- fscores(dif.mirt1, full.scores = TRUE, full.scores.SE = TRUE, method = 'EAP')
theta.est <- dif.mirt1.fscores[, "F1", drop = FALSE]
dif.re1 = matrix(NA, nrow = ncol(halluc.mat), ncol = 3)
for(i in 1:ncol(halluc.mat)){
  dif.stat <- DIF_logit(theta.est,halluc.mat[,i], group0)
  if(is.na(dif.stat)){next}
  dif.re1[i,] = c(dif.stat$stat, dif.stat$df, dif.stat$power)
}
boxplot(dif.re1[,3])
boxplot(-log(dif.re1[,3]))
dif.re2 = matrix(NA, nrow = ncol(halluc.mat), ncol = 3)
for(i in 1:ncol(halluc.mat)){
  dif.stat <- DIF_MH(theta.est,halluc.mat[,i], group0)
  if(is.na(dif.stat)){next}
  dif.re2[i,] = c(dif.stat$stat, dif.stat$df, dif.stat$power)
}
boxplot(dif.re1[,3], dif.re2[,3], main='Power',xlab= c('logit vs MH'))
sum(dif.re1[,3]<0.8,na.rm = T)
sum(dif.re2[,3]<0.8,na.rm = T)
plot(dif.re1[,3],dif.re2[,3])
abline(coef=c(0,1))
