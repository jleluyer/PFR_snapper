ls()
rm(list=ls())
ls()

sessionInfo()
#https://support.bioconductor.org/p/56568/
#setwd("/Users/jeremyleluyer/Documents/post-doc/projects/transgenic/transciptome/muscle/03_results/glm/03_results/")

## test Glm approahc with edger ##
library('edgeR')
setwd("~/Documents/Projets/Snapper/")
counts <- read.table("02_data/htseq_matrix.txt",header=T,na.strings="NA")
annotation<-read.csv("01_info_files/snapper.annotation.csv",header=T,sep=";",na.strings="NA")

vectorname<-row.names(targets)
vectorname

#counts<-counts[,vectorname]
# Create design
targets <- read.table("02_data/design.txt",header=T,na.strings="NA")
targets$genotype <- relevel(targets$genotype, ref=c("Domestic"))
#targets$genotype <- relevel(targets$genotype, ref=c("Domestic"))
targets


#reassignng column name
rownames(counts) <- counts[,1]
counts[,1]<- NULL
counts<-counts[,vectorname] #force order

dim(counts)
str(counts)
colSums(counts) / 1e06

#nb gene with low count
table( rowSums( counts ) )[ 1:30 ]

#create object
Group <-  factor(paste(targets$temperature,targets$genotype,sep="."))
Group
cds <- DGEList( counts , group = Group )
head(cds$counts)
cds$samples

levels(cds$samples$group)

#filtering low counts
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 2, ]
dim( cds )

#calcultate normalization factors
cds <- calcNormFactors( cds )
cds$samples

#library size
cds$samples$lib.size*cds$samples$norm.factors

#Create design (can be created manually)
design <- model.matrix(~genotype *temperature, data=targets)
colnames(design)
design
#library size
cds$samples$lib.size*cds$samples$norm.factors
cds$samples$norm.factors
#estimate dispersion
cds <-estimateDisp(cds,design) #estimate all dispersions in one run #add different prior values (prior.df = INT)
names( cds )
write.table(cds$counts, file="countTEMP.txt",quote=F)
#create a logcpm matrix
logcpm <- cpm(cds, prior.count=2, log=TRUE)
write.table(logcpm, file="03_results/logcpm_edger_snapper.txt", sep = "\t" , quote=F,row.names = TRUE)

# Testing comparison any treatment test
fit <- glmFit(cds, design)
colnames(fit)
lrt <- glmLRT(fit, coef=4) #see contrast explanation https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
lrt$comparison # which groups have been compared
topTags( lrt , n = 20 , sort.by = "PValue" )
topTags(lrt)


#### plot#### 
#plot MDS
plotMDS( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )
d<-as.data.frame(plotMDS( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) ))

meanVarPlot <- plotMeanVar( cds ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , #these are arguments about what is plotted
                            pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" ) #these arguments are to make it look prettier


#sort byt fold change
resultsByFC.tgw <- topTags( lrt , n = nrow( lrt$table ) , sort.by = "logFC" )$table

head( resultsByFC.tgw )

#store full toptags
resultsTbl.tgw <- topTags( lrt , n = nrow( lrt$table ) )$table
head(resultsTbl.tgw)

#extract significant with FDR < 0.05
resOrdered = resultsTbl.tgw[order(resultsTbl.tgw$FDR),]
resSig = subset(resOrdered, FDR<0.01)
#write table in a file
head(resSig)
dim(resSig)
write.table(resSig, file="03_results/interaction_temp_geno.txt",sep= "\t",quote=F )


# Heatmap interaction
library(circlize)
## heatmap top 20 tags

#source("https://bioconductor.org/biocLite.R")
#biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

log_cpm20topInter<-read.table("~/project_snapper/03_results/GO_enrichment/DE/interaction/topinteraction_cpm.txt",sep="\t",header=T)
#rownames(log_cpm20topInter) <- log_cpm20topInter[,1]
#log_cpm20topInter[,1] <- NULL
matrix_locpminter<-as.matrix(log_cpm20topInter)

#Complex heatmap
type = gsub("s\\d+_", "", colnames(matrix_locpminter))
ha = HeatmapAnnotation(df = data.frame(type = type))
tiff(file = "top_20_tags_DE_temperature.tiff", width = 12000, height = 9000, units = "px", res = 800)
Heatmap(matrix_locpminter,name = "Log CPM", km=1, col = colorRamp2(c(-1, 0, 4), c("dodgerblue4", "white", "yellow")),
        show_row_names = T, show_column_names = T,clustering_distance_rows = "euclidean")
dev.off()


### plot interaction
tiff(file = "redaction/Figure_geneinteraction.tiff",width = 15, height = 11, units = "cm",res=300)
library(ggplot2)
log_counts_interaction<-read.table("~/project_snapper/03_results/GO_enrichment/DE/interaction/transposed_interaction_logpcm.txt",sep="\t",header=T)
str(log_counts_interaction)
plot<- (ggplot(data = log_counts_interaction, aes(x=Temperature, y=TRINITY_DN82676_c0_g1_i2,group=Genotype)) 
        + geom_point(aes(fill=Genotype,shape=Genotype),size=5)
        + scale_shape_manual(values=c(21,22))
        + scale_fill_manual(values=c("grey20","white"))
        + scale_y_continuous(limits = c(-1,4))
        + stat_summary(fun.y = mean, geom="line",aes(linetype=Genotype))
) + theme_bw()+ 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(colour="black",size=14,face="plain"), 
        axis.text.x = element_text(colour="black",size=14,face="plain"),
        axis.title.x=element_text(colour="black",size=14,face="bold"),
        axis.title.y=element_text(colour="black",size=14,face="bold"),
        legend.title=element_blank()
        )+ ggtitle("")+labs(y="LogCPM")
plot
dev.off()
##########################################################################################################################
#################################### glm for additive traits snapper #####################################################
##########################################################################################################################


#### Carefull remove genes in interaction first

#see post
ls()
rm(list=ls())
ls()
#https://support.bioconductor.org/p/56568/
#setwd("/Users/jeremyleluyer/Documents/post-doc/projects/transgenic/transciptome/muscle/03_results/glm/03_results/")
## test Glm approahc with edger ##
library('edgeR')
setwd("~/Documents/Projets/Snapper/")

counts <- read.table("02_data/htseq_matrix.txt",header=T,na.strings="NA")



# Create design
targets <- read.table("02_data/design.txt",header=T,na.strings="NA")
targets


#reassignng column name
rownames(counts) <- counts[,1]
counts[,1]<- NULL
dim(counts)
str(counts)
colSums(counts) / 1e06

#nb gene with low count
table( rowSums( counts ) )[ 1:30 ]
#create object
Group <-  factor(paste(targets$temperature,targets$genotype,sep="."))
Group
cds <- DGEList( counts , group = Group )
head(cds$counts)
cds$samples
levels(cds$samples$group)

#filtering low counts
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 2, ]
dim( cds )

#calcultate normalization factors
cds <- calcNormFactors( cds )
cds$samples

#library size
cds$samples$lib.size*cds$samples$norm.factors

#Create design (can be created manually)
design.oneway <- model.matrix(~0+Group)
colnames(design.oneway) <- levels(Group)
design.oneway

#library size
cds$samples$lib.size*cds$samples$norm.factors
cds$samples$norm.factors
#estimate dispersion
cds <-estimateDisp(cds,design.oneway) #estimate all dispersions in one run #add different prior values (prior.df = INT)
names( cds )

#make contrast
#my.contrasts <- makeContrasts(
# cond.geno = (bassin.wild-bassin.trans) - (stream.wild-stream.trans),
#cond.enviro = (bassin.wild-stream.wild) - (bassin.trans-stream.trans),
#levels=design.oneway)
#my.contrasts

# Testing comparison any treatment test
fit <- glmFit(cds, design.oneway)
colnames(fit)
?glm()
#lrt <- glmLRT(fit, contrast=c(1,-1,1,-1)) #geno effecto#see contrast explanation https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
lrt <- glmLRT(fit, contrast=c(-1,-1,1,1)) # temperature effect


lrt$comparison # which groups have been compared
topTags( lrt , n = 20 , sort.by = "PValue" )
topTags(lrt)


#### plot#### 
#plot MDS
plotMDS( cds , main = "MDS Plot for Count Data", labels = colnames( cds$counts ) )

meanVarPlot <- plotMeanVar( cds ,show.raw.vars=TRUE , show.tagwise.vars=TRUE , show.binned.common.disp.vars=FALSE , show.ave.raw.vars=FALSE , dispersion.method = "qcml" , NBline = TRUE , nbins = 100 , #these are arguments about what is plotted
                            pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" ) #these arguments are to make it look prettier



#sort byt fold change
resultsByFC.tgw <- topTags( lrt , n = nrow( lrt$table ) , sort.by = "logFC" )$table

head( resultsByFC.tgw )

#store full toptags
resultsTbl.tgw <- topTags( lrt , n = nrow( lrt$table ) )$table
head(resultsTbl.tgw)

#extract significant with FDR < 0.05
resOrdered = resultsTbl.tgw[order(resultsTbl.tgw$FDR),]
resSig = subset(resOrdered, FDR<0.01)
reSig2=subset(resSig,abs(logFC)>2)

dim(reSig2)
#write table in a file
write.table(reSig2, file="03_results/DEGs/temperature_additive_deg.txt", sep = "\t" , row.names = TRUE,quote=F)

###########################################################################################################################
########################### Genotype effect ###############################################################################
##########################################################################################################################

colnames(fit)
lrt <- glmLRT(fit, contrast=c(1,-1,1,-1)) #geno effecto#see contrast explanation https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf


lrt$comparison # which groups have been compared
topTags( lrt , n = 20 , sort.by = "PValue" )
topTags(lrt)

#sort byt fold change
resultsByFC.tgw <- topTags( lrt , n = nrow( lrt$table ) , sort.by = "logFC" )$table

head( resultsByFC.tgw )

#store full toptags
resultsTbl.tgw <- topTags( lrt , n = nrow( lrt$table ) )$table
head(resultsTbl.tgw)

#extract significant with FDR < 0.05
resOrdered = resultsTbl.tgw[order(resultsTbl.tgw$FDR),]
resSig = subset(resOrdered, FDR<0.01)
reSig2.geno=subset(resSig,abs(logFC)>2)


dim(reSig2.geno)
#write table in a file
write.table(reSig2.geno, file="03_results/DEGs/genotype_additive_deg.txt", sep = "\t" , row.names = TRUE,quote=F)




###########################################################################################################################
########################### cluster additive gene #########################################################################
##########################################################################################################################
log_deg_additive<-read.table("03_results/DEGs/list_additive_logcpm.txt",header=T)
str(log_deg_additive)
#install.packages("pheatmap")
library(pheatmap)
mat <- as.matrix(log_deg_additive)
mat  <- mat - rowMeans(mat)
#anno <- as.data.frame(colData(rld)[, c("cell","dex")])
aka2 = data.frame(ID = factor(c("dh","dh","dh","dl","dl","dl","wh","wh","wl","wl","wl")))
subj<-colnames(mat)
rownames(aka2)<-subj
aka3 = list(ID = c(wh = "dodgerblue2", dh="tomato2", dl="darkgoldenrod2", wl="springgreen3"))

tiff(file = "03_results/heatmap_clusters.tiff", width = 15, height = 15, units = "cm", res = 300)
set.seed(1)
pheatmap(mat,kmeans_k=8,annotation_col = aka2, 
         annotation_legend = FALSE,
         annotation_names_col = FALSE,
         annotation_colors = aka3[1])
dev.off()

###########################################################################################################################
########################### Visualize interaction #########################################################################
##########################################################################################################################

library(ggplot2)
install.packages("ggplot")

#For interaction
log_counts_interaction <- as.data.frame(t(read.table("03_results/logcpm_edger_snapper.txt",header=T)))
write.table(log_counts_interaction, file="03_results/logcpm_edger_snapper_plot.txt", sep = "\t" , row.names = TRUE,quote=F)
log_counts_interaction <- as.data.frame(read.table("03_results/logcpm_edger_snapper_plot.txt",header=T))
plot<- (ggplot(data = log_counts_interaction, aes(x=genotype, y=TRINITY_DN82676_c0_g1_i2,group=temperature,shape=genotype)) 
+ geom_point(aes(colour=genotype)) 
+ scale_y_continuous(limits = c(-5, 15))
+ stat_summary(fun.y = mean, geom="line",aes(color=environment))
+ ggtitle("cluster 9")
)
plot
 


###########################################################################################################################
########################### RDA on gene expression #########################################################################
##########################################################################################################################


ls()
rm(list=ls())
ls()

setwd("~/Documents/Projets/Snapper/")

#charger matrice expression et la transposer (individus en objet et marqueur en variables)
genet=read.table("03_results/logcpm_edger_snapper.txt",header=TRUE)
genet_trans=t(genet)
row.names(genet_trans)
design=read.table("02_data/design.txt",header=T)
df<-merge(design, genet_trans, by=0, all=TRUE)
rownames(df)=df$Row.names
df$Row.names <- NULL
str(df)
#mettre dans l'ordre
df=df[ order(df$temperature, df$genotype), ]


#produire PCOA à partir d'une matrice euclidean
library(ape)
library(vegan)
library(cluster)

genet.pcoa=pcoa(daisy(df[,3:ncol(df)], metric="euclidean"))
genet.pcoa$values


write.table(genet.pcoa$vector, file="03_results/Genet_PC-factor.txt", sep="\t",quote=F)

#Produire sélection de variable de la db-RDA (var2 = sex, var3 = river, var4 = treatment)
ordistep(rda(genet.pcoa$vectors[,1:3]~df$temperature+df$genotype, scale=F))

#Produire db-RDA avec modèle sélectionné par ordistep, soit sex + river
rda_genet=rda(formula = genet.pcoa$vectors[, 1:3] ~ df$genotype + df$temperature, scale = F)

#test significance of global model
anova(rda_genet, step=1000)


#Test significance for each variable
anova(rda_genet, step=1000, by='margin')

#Calculé le pourcentage de variation expliqué par le modèle global      
RsquareAdj(rda_genet)

#partial db-RDA for river effect alone / control for Sex
rda_genet2=rda(genet.pcoa$vectors[,1:3],as.numeric(df$genotype),as.numeric(df$temperature), scale=F)

#Testé significativité du modèle
anova(rda_genet2, step=1000)

#Calculé le % de variance expliquée par River seul 
RsquareAdj(rda_genet2)

#produire partial db-RDA effet Sex seul / controlé par River
rda_genet3=rda(genet.pcoa$vectors[,1:3],as.numeric(df$temperature),as.numeric(df$genotype), scale=F)

#Testé significativité du modèle
anova(rda_genet3, step=1000)

RsquareAdj(rda_genet3)


# GRAPHIQUES DE LA RDA;
# -------------------;
library(ggplot2)
sommaire = summary(rda_genet)
sommaire
df1  <- data.frame(sommaire$sites[,1:2])       # PC1 and PC2
df2  <- data.frame(sommaire$species[,1:2])

# prepare df1 with group info
df1<-merge(df1, design, by=0, all=TRUE)
rownames(df1)=df1$Row.names
df1$Row.names <- NULL
df1$group <-  factor(paste(df1$temperature,df1$genotype,sep="-"))

# loadings for PC1 and PC2
rda.plot.rna <- ggplot(df1, aes(x=RDA1, y=RDA2,colour=group)) + 
  geom_point(size=4) +
  # stat_ellipse(df1,aes(x=RDA1, y=RDA2,color=group,group=group),type = "norm")+
  # geom_text(aes(label=rownames(df1)),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed()

df2subset<-df2[c("Axis.1","Axis.2"),]

rda.biplot.rna <- rda.plot.rna +
  geom_segment(data=df2subset, aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(0.01,"npc"))) + labs(x="RDA1 (50.1 %)",y="RDA2 (17.7%)") 

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black",size=15),
             axis.text.y=element_text(colour="black",size=15),
             axis.title.x=element_text(colour="black",size=15,face="bold"),
             axis.title.y=element_text(colour="black",size=15,face="bold"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"),
             legend.title=element_blank(),
             aspect.ratio=1,
             legend.background = element_rect(size=0.4, linetype="solid", 
                                              colour ="black"),
             legend.position = c(0.8,0.85),
             legend.text=element_text(size=15),
             legend.key = element_rect(fill=NA)) 
couleurs=c("tomato2","dodgerblue2","darkgoldenrod2","springgreen3") 
graph.rda<-rda.biplot.rna +theme + 
  annotate("text", x = c(8,-3), y = c(-1,5), label = c("Temperature***","Genotype**"),size=c(6,6),fontface =2) +
  scale_color_manual(values=couleurs)
graph.rda

tiff(file = "03_results/rdas_global.tiff", width = 25, height = 15, units = "cm", res = 300)
graph.rda
dev.off()
