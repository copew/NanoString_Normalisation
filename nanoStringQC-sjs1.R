library(NanoStringNorm)
library(ggplot2)
library(reshape2)
dev.off()

coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000', '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
coliClust <- data.frame(c("1", "2", "3", "4ER+", "4ER-", "5", "6", "7","8", "9", "10"),coliClust)
colnames(coliClust)[1] <- "Clust11"
colnames(coliClust)[2] <- "Colour Codes"
coliClust$`Colour Codes` <- as.character(coliClust$`Colour Codes`)

###This divides iC4 into iC4ER+ and iC4ER-. You might want to omit that for now, in that case do 
          
               
coliClust10 <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#8B0000','#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
coliClust10 <- data.frame(c("1", "2", "3", "4", "5", "6", "7","8", "9", "10"),coliClust10)
colnames(coliClust10)[1] <- "iClust"
colnames(coliClust10)[2] <- "Colour Codes"
coliClust10$`Colour Codes` <- as.character(coliClust10$`Colour Codes`)

save(coliClust10, file = "iC10_colour.RData")

#str(coliClust10)

##remember to change input and output folder names!!!!!!!!!!!!!!####

batchNumber <- "Combined"
metadataFolder <- "/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Nanostring/Other Info/"
nanostringDataFolder <- paste0("/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Nanostring/Digital Analyser/",batchNumber,"/")
outputFolder <- paste0("/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Nanostring/Nanostring QC/",batchNumber,"QC/")

dir.create(outputFolder,showWarnings = FALSE)

iC10genesFile <- paste0(metadataFolder,"for_nanostring.txt")
metabricExpressionData <- paste0(metadataFolder,"METABRIC_CURTIS_EXPRESSION.RData")


pdf(paste0(outputFolder,"metrics.pdf"), width=11, height=8)

nanoStringFiles <- list.files(nanostringDataFolder, pattern = "MB-AD", recursive = T, full.names = T)
nanoStringData <- lapply(nanoStringFiles, read.csv,
                         skip = 26, header = T, nrows = 773,quote=NULL,
                         stringsAsFactors = F)
#extract metabric id from file id
metabricID <- regmatches(nanoStringFiles, regexpr("MB-AD-[0-9]*", nanoStringFiles))

#rename count column to metabric ID
for (i in seq_along(nanoStringData)){
  colnames(nanoStringData[[i]])[colnames(nanoStringData[[i]])=="Count"] <-metabricID[[i]]
}

#convert list to dataframe
nanoStringData <- do.call("cbind", nanoStringData)
nanoStringData <- nanoStringData[, !duplicated(colnames(nanoStringData), fromLast = FALSE)] 

#dim(nanoStringData)

##exclude dodgy samples
#nanoStringData<-nanoStringData[,!colnames(nanoStringData) %in% c("MB-AD-0414", "MB-AD-0129")]
#head(nanoStringData)
#dim(nanoStringData)

librarySize=data.frame(sample=names(nanoStringData[,c(4:ncol(nanoStringData))]),size=colSums(nanoStringData[,c(4:ncol(nanoStringData))]))
ggplot(librarySize,aes(x=sample,y=size))+
  geom_bar(stat="identity",color="black",fill="gray80")+
  labs(x="Sample",y="Total counts",title="Distribution of raw counts across samples")+
  theme_bw(base_size=16)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




#normalisation using NanoStringNorm
normalisedNanoStringData <- NanoStringNorm(
  x = nanoStringData,
  anno = NA,
  CodeCount = 'sum',
  Background = 'mean',
  SampleContent = 'housekeeping.sum',
  round.values = FALSE,
  take.log = FALSE,
  return.matrix.of.endogenous.probes = TRUE
)


normalisedLibrarySize <- data.frame(sample=colnames(normalisedNanoStringData),size=colSums(normalisedNanoStringData))
meltednormalisedNanoStringData<-melt(normalisedNanoStringData)

# ggplot(normalisedLibrarySize,aes(x=sample,y=size))+
#   geom_bar(stat="identity",color="black",fill="gray80")+
#   labs(x="Sample",y="Normalised total counts",title="Distribution of normalised counts across samples")+
#   theme_bw(base_size=16)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# ggplot(meltednormalisedNanoStringData,aes(x=Var2,y=log(value+1.1)))+
#   geom_boxplot(color="black",fill="gray80")+
#   labs(x="Sample",y="Log Normalised total counts",title="Distribution of normalised counts within sample")+
#   scale_y_continuous(breaks=c(0,2,4,6,8,10))+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))


meltedMeans<-data.frame(gene=rownames(normalisedNanoStringData),mean=rowMeans(normalisedNanoStringData))
meltedMeans[meltedMeans$mean<=1,]
meltedMeans[meltedMeans$mean>1000,]

#distribution of all genes
# ggplot(meltedMeans,aes(x=(mean)))+
#   geom_histogram(binwidth = 0.15, color="black", fill="gray80")+
#   theme_bw(base_size=16)+
#   scale_x_log10(breaks=c(0,1,5,10,25,50,100,250,500,1000,2500,5000,10000))+
#   labs(x="Mean Normalised Nanostring Count across samples (log scale)",y="Number of genes",title="Distribution of gene counts")


#220 genes identfied by marcus
iCgenes220 <- scan(iC10genesFile,what = character())

dim(normalisedNanoStringData)

saveRDS(normalisedNanoStringData, "nanostring_total.RDS")

nanoStringIC<-normalisedNanoStringData[rownames(normalisedNanoStringData) %in% iCgenes220,]
meltedMeans<-data.frame(gene=rownames(nanoStringIC),mean=rowMeans(nanoStringIC))
# ggplot(meltedMeans,aes(x=(mean)))+
#   geom_histogram(binwidth = 0.25, color="black", fill="gray80")+
#   theme_bw()+
#   scale_x_log10(breaks=c(0,1,5,10,25,50,100,250,500,1000,2500,5000,10000))+
#   labs(x="Mean Normalised Nanostring Count across samples (log scale)",y="Number of genes",title="Distribution of iC10 classifier gene counts")

dev.off()


load(metabricExpressionData)

mb_expression <- Exp[which(Exp[,1] %in% rownames(nanoStringIC)),]

nanoStringIC=nanoStringIC[sort(rownames(nanoStringIC)),]
mb_expression[,1]<-as.character(mb_expression[,1])
mb_expression <- mb_expression[match(rownames(nanoStringIC), mb_expression[,1]),]
mb_expression=mb_expression[!is.na(mb_expression$SYMBOL),]
colnames(mb_expression) <- sub(".", "-", colnames(mb_expression), fixed=T)
colnames(mb_expression)

colnames(nanoStringIC) <- sub("-AD-", "-", colnames(nanoStringIC), fixed=T)
colnames(nanoStringIC)

rownames(mb_expression) <- mb_expression[,1]
mb_expression <- mb_expression[,which(colnames(mb_expression) %in% colnames(nanoStringIC))]

#dim(mb_expression)
#colnames(mb_expression)
nanoStringIC <- nanoStringIC[,colnames(mb_expression)]
nanoStringIC <- nanoStringIC[rownames(nanoStringIC)%in% rownames(mb_expression),]
dim(nanoStringIC)

mean(rownames(nanoStringIC) == rownames(mb_expression))
mean(colnames(nanoStringIC) == colnames(mb_expression))




###all the r


cors <- data.frame(
  cors=sapply(1:ncol(mb_expression),function(i){
    cor(nanoStringIC[,i], mb_expression[,i],method="spearman")
  }),
  MB.ID=colnames(mb_expression)
  )

 ggplot(cors, aes(x=cors)) + geom_histogram(binwidth=0.01)+theme_bw(base_size=16)+labs(x="Correlation of Samples", y="Frequency")

cors_asc <- cors[with(cors, order(cors)), ] 
head(cors_asc, 25) 


 #View(cors)


 
##plots per sample

pdf(paste0(outputFolder,"Plots_per_sample.pdf"), width=8.3, height=11.7)
par(mfrow=c(4, 4))
for (i in 1:ncol(mb_expression)) {
   plot(log2(nanoStringIC[,i]), mb_expression[,i], main=colnames(mb_expression)[i], xlab="Nanostring",
        ylab="Illumina")
 
   linMod<-lm(unlist(mb_expression[,i]) ~ log2(as.numeric(nanoStringIC[,i]+1)))
   modsum = summary(linMod)
   r2 = modsum$adj.r.squared
   my.p = modsum$coefficients[2,4]
   rp = vector('expression',2)
 rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
   rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                      list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
   abline(linMod,col="red")
   legend('topleft', legend = rp, bty = 'n')
}
dev.off()

##plots per gene

pdf(paste0(outputFolder,"Plots_per_gene.pdf"))
par(mfrow=c(2, 2))
for (i in 1:nrow(mb_expression)) {
  plot(log2(nanoStringIC[i,]+1), mb_expression[i,], main=rownames(mb_expression)[i], xlab="Nanostring",
       ylab="Illumina")
  linMod<-lm(unlist(mb_expression[i,]) ~ log2(as.numeric(nanoStringIC[i,]+1)))
  modsum = summary(linMod)
  r2 = modsum$adj.r.squared
  my.p = modsum$coefficients[2,4]
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE),
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE),
                     list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
  abline(linMod,col="red")
  legend('topleft', legend = rp, bty = 'n')

  }
dev.off()

#r values per gene
mb_expression_matrix <- data.matrix(mb_expression) 
cors_gene <- data.frame(
  cors=sapply(1:nrow(mb_expression_matrix),function(i){
    cor(nanoStringIC[i, ],mb_expression_matrix[i, ], method="spearman")
  }), 
  Genes = rownames(mb_expression_matrix)
  )

ggplot(cors_gene, aes(x=cors)) +
  geom_histogram(binwidth=0.02)+
  theme_bw(base_size = 16)+
  ylab("Frequency")+
  xlab("Correlation of Genes")

summary(cors_gene$cors)

cors_gene_asc <- cors_gene[with(cors_gene, order(cors)),]
head(cors_gene_asc, 15)
tail(cors_gene_asc, 15)


# r value per cluster

##number of samples per cluster
MB_intClust <- read.csv("/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Nanostring/Other info/MB intClust.csv")
colnames(MB_intClust)[3] <- "ER Status"
MB_intClust$Clust11 <- paste(MB_intClust$IntClust, MB_intClust$`ER Status`, sep="")
head(MB_intClust)

clust_corr11 <- merge(cors, MB_intClust, by='MB.ID')
head(clust_corr11)

clust_corr11$Clust11 <- as.factor(clust_corr11$Clust11)
clust_corr11 <- within(clust_corr11, Clust11 <- factor(Clust11, levels=c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10")))

clust_corr11 <- merge(clust_corr11, coliClust, by="Clust11")
clust_corr11$`Colour Codes` <- as.character(clust_corr11$`Colour Codes`)
ggplot(clust_corr11, aes(x=Clust11,y=cors)) + 
  theme_bw(base_size=16)+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.3, col=clust_corr11$`Colour Codes`, size=2)+
  xlab("Integrative Clusters")+
  ylab("Correlation")




#View(clust_corr11)
##number of samples
ggplot(clust_corr11, aes(x=Clust11, fill=Clust11))+
  geom_bar()+
  theme_bw(base_size = 20)+
  theme(legend.position = "none")+
  ylab("Number of Samples")+
  xlab("Integrative Clusters")+
  scale_fill_manual(values=coliClust$`Colour Codes`)





clust_corr <- merge(cors, MB_intClust, by='MB.ID')
clust_corr$IntClust <- as.factor(clust_corr$IntClust)
clust_corr <- merge(clust_corr, coliClust10, by="IntClust")
clust_corr$`Colour Codes` <- as.character(clust_corr$`Colour Codes`)
str(clust_corr)

colnames(clust_corr)
summary(clust_corr$cors)
save(clust_corr, file="clust_corr.RData")


View(clust_corr)
ggplot(clust_corr, aes(x=IntClust,y=cors)) + 
  theme_bw(base_size=16)+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width=0.3, col=clust_corr$`Colour Codes`, size=2)+
  xlab("Integrative Clusters")+
  ylab("Correlation")







IntClust4 <- clust_corr[clust_corr$IntClust=="4", ]
#colnames(IntClust4)
ggplot(IntClust4, aes(X4.hormone.status,cors)) + 
  theme_bw()+
  geom_boxplot(outlier.shape=NA, width=0.5)+
  geom_jitter(width=0.25,aes(col=X4.hormone.status))+
  xlab("Oestrogen Receptor Status")+
  ylab("Correlation")+
  theme(legend.title=element_blank())+ 
  expand_limits(y=c(0,1))
 
# ###comes with package
# # pdf('NanoStringNorm_Example_Plots_All.pdf'); 
# #  Plot.NanoStringNorm( x = normalisedNanoStringData, label.best.guess = TRUE,
# #    plot.type = 'all' )
# # dev.off();


# ##correlation by ER status
# Metabric_Clinical <- read.csv("~/OneDrive - University Of Cambridge/Documents/PhD/Nanostring/Other Info/Metabric.Clinical.csv")
# str(Metabric_Clinical)
# names(Metabric_Clinical)[names(Metabric_Clinical) == 'METABRIC_ID'] <- 'MB.ID'
# corr_clinical <- merge(clust_corr, Metabric_Clinical, by="MB.ID")
# ggplot(corr_clinical, aes(er_status,cors)) + 
#   theme_bw()+
#   geom_boxplot(outlier.shape=NA, width=0.5)+
#   geom_jitter(width=0.25,aes(col=er_status))+
#   xlab("Oestrogen Receptor Status")+
#   ylab("Correlation")+
#   theme(legend.title=element_blank())+
#   scale_color_hue(labels = c("Positive", "Negative", "NA"))+ 
#   expand_limits(y=c(0,1))+
#   scale_fill_manual(values=coliClust10)
# 
# ##compare by er status
# cor_ERpositive <- corr_clinical[corr_clinical$er_status=="pos", 3]
# cor_ERnegative <- corr_clinical[corr_clinical$er_status=="neg", 3]
# 
# t.test(cor_ERnegative, cor_ERpositive)
# 
# 
# ##compare by her2 status
# corr_clinical$HER2_status <- gsub("null", "0", corr_clinical$HER2_status)
# summary(aov(cors~HER2_status, corr_clinical))
# ggplot(corr_clinical, aes(HER2_status,cors)) + 
#   theme_bw()+
#   geom_boxplot(outlier.shape=NA, width=0.5)+
#   geom_jitter(width=0.25,aes(col=HER2_status))+
#   xlab("HER2 Status")+
#   ylab("Correlation")+
#   theme(legend.title=element_blank())+
#    expand_limits(y=c(0,1))
# 
# 
# summary(clust_corr)
# 
# 
# 
# 
#   
# ##count diffeerent between clusters?
# # corr_cluster <- lm(cors~IntClust, data=clust_corr)
# # summary(corr_cluster)
# # anova(corr_cluster)
# corr_anova <- aov(cors~IntClust, data=clust_corr)
# summary(corr_anova)


write.csv(nanoStringIC, "/Users/cope01/Documents/OneDrive - University Of Cambridge/Documents/PhD/Nanostring/Other Info/nanoStringIC.csv")
