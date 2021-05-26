
## Setting up packages and data
#required packages

deps = c("vegan","dae","plyr","scales","lme4","nlme","MASS","sets","phyloseq",
         "mgcv","indicspecies", "randomForest","mvabund", "MuMIn","paleotree",
         "TreeDist","phytools","MicEco","ecodist","mmod","poppr","psych",
         "ape","magrittr","pegas");


for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

#sample metadata

url <- "https://ndownloader.figshare.com/files/27898803"
download.file(url, destfile="all_sample.csv", 'libcurl')

url <- "https://ndownloader.figshare.com/files/27898749"
download.file(url, destfile="tick_biodiv.csv", 'libcurl')
tick <- read.csv(file="tick_biodiv.csv",sep=";")

#OTU numbers

url <- "https://ndownloader.figshare.com/files/27898776"
download.file(url, destfile="tick99.shared", 'libcurl')
otu <- read.delim(file="tick99.shared",header=TRUE)
row.names(otu)<- otu$Group
otu <- otu[,-c(1:3)]

#OTU names

url <- "https://ndownloader.figshare.com/files/27898767"
download.file(url, destfile="tick99.taxonomy", 'libcurl')
cons.tax <- read.table(file = "tick99.taxonomy", 
                       header = TRUE, row.names = 1)
otu.taxa <- gsub("\\(\\d*\\)", "", cons.tax[, "Taxonomy"])
otu.taxa <- gsub("\"", "", otu.taxa)
otu.taxa <- gsub("unclassified;", "", otu.taxa)
otu.taxa <- gsub(";$", "", otu.taxa)
otu.taxa <- gsub("^.*;", "", otu.taxa)
names(otu.taxa) <- rownames(cons.tax)

# metrics calculcated in mothur

url <- "https://ndownloader.figshare.com/files/27898791"
download.file(url, destfile="UF_weight.csv", 'libcurl')

url <- "https://ndownloader.figshare.com/files/27898770"
download.file(url, destfile="tick99.tre", 'libcurl')

#tick population genetics

url <- "https://ndownloader.figshare.com/files/28035786"
download.file(url, destfile="tick_allele_R2.csv", 'libcurl')
tick_pop <- read.genalex("tick_allele_R2.csv", ploidy=2, sep=";")

#calculating genetic distances

ticdist   <- provesti.dist(tick_pop)

## Descriptive statistics

otunum <- otu
otu$amplSum <- rowSums(otunum)
summary(otu$amplSum)

#Selection of samples
otu_tick <- subset(otu, amplSum > 500)
stopifnot(sum(otu_tick$group != tick$group) == 0) #does everything match?
otu_tick <- subset(otu_tick, select=c(-amplSum))

#Let's take out potential contaminating amplicons and rarefy the data
otu_rare <- rrarefy(otu_tick, 500)
dt <- otu_tick
dt[dt>0] <-1 
otu_rare <- otu_rare[, colSums(dt != 0) > 1]
otu_rare <- otu_rare[, colSums(otu_rare != 0) > 1]
otu_rare = as.data.frame(otu_rare)

# core analysis

nseqs <- apply(otu_rare, 1, sum)[1]
otu.relabund <- otu_rare / nseqs
otu.med <- aggregate(otu.relabund, by=list(tick$group), median)
otu.et <- otu.med$Group.1 #and format the table
otu.med <- otu.med[,-1]

total.n.otus <- ncol(otu.med)  #total number of OTUs across all samples
sobs <- rowSums(otu.med>0)
otu.distribution <- quantile(sobs, c(0.025, 0.50, 0.975))  #95% CI for OTUs

present90 <- apply(otu.med>0, 2, sum) >= 0.90*nrow(otu.med) #true of false for OTUs in 90% of samples
notus.present90 <- sum(present90) #number of OTUs in 90% 
otu.present90 <- colnames(otu.med)[present90] #what are the OTUs in 90% ?
otu_names.present90 <- levels(factor(otu.taxa[otu.present90])) #what taxa are represented?
otu.percent.present90 <- 100 * sum(otu.med[,present90])/nrow(otu.med) #what percentage of seqs do they represent?

presentAll <- apply(otu.med>0, 2, sum) == nrow(otu.med)
otu.presentAll<- colnames(otu.med)[presentAll] 
notus.presentAll <- sum(presentAll)
otu_names.presentAll <- levels(factor(otu.taxa[otu.presentAll]))
otu.percent.presentAll <- 100 * sum(otu.med[,presentAll])/nrow(otu.med)


##Alphadiversity

# Phylodiversity already in metadata from mothur 
# Inverse simpson
tick$simpson <- diversity(otu_rare, index="invsimpson")

# Mixed effects models
m1.nlme = lme(simpson ~ masl*tick_stage*month+(I(masl^2)), random= ~1|site,  data=tick)
m1.nlme = lme(simpson ~ masl*tick_stage*month, random= ~1|site,  data=tick)
m1.nlme = lme(simpson ~ masl+tick_stage+month+masl:tick_stage+masl:month+tick_stage:month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+masl:month+tick_stage*month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+masl:month+tick_stage+tick_stage:month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+tick_stage+tick_stage:month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+tick_stage, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ tick_stage, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ 1, random= ~1|site, data=tick)


summary(m1.nlme)
r.squaredGLMM(m1.nlme)

m2.nlme = lme(phylodiv ~ masl*tick_stage*month+(I(masl^2)), random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl*tick_stage*month, random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl+tick_stage+month+masl:tick_stage+masl:month+tick_stage:month, random= ~1|site, data=tick)
m2.nlme = lme(phylodiv ~ masl+masl:month+tick_stage*month, random= ~1|site, data=tick)
m2.nlme = lme(phylodiv ~ masl+month+tick_stage+tick_stage:month, random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl+tick_stage+ month , random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl+month, random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ 1, random= ~1|site, data=tick)

summary(m2.nlme)
r.squaredGLMM(m2.nlme)

##Figure 3
#Fig3a

tick_phylodiv <- tick[,c(7:9)]
tick_phylodiv$site <- rep("FE",nrow(tick_phylodiv))
pred_phylodiv <- predict(m2.nlme, newdata=tick_phylodiv)

p1 <- ggplot(tick, aes(x=masl, y=phylodiv, colour=as.factor(month))) +
  scale_colour_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  geom_point(size=3) +
  geom_line(aes(y=pred_phylodiv, group=month, size="Months")) +
  scale_size_manual(name="Months", values=c("Months"=0.5, "Population"=3)) +
  theme_bw(base_size=22) 
print(p1)

#Fig3b

p2 <- ggplot(tick, aes(x=tick_stage,y=simpson,color=tick_stage))+
  geom_jitter(position=position_jitter(0.1))+
  stat_summary(fun=median, geom="crossbar", ymax=-100, ymin=-100)+
  scale_colour_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  theme_bw(base_size=22) 
print(p2)


##Betadiversity

betad1 <- vegdist(otu_rare, method="bray") #Braycurtis

betad2 <- vegdist(otu_rare, method="bray", binary=TRUE) #Jaccard


unif_weig <- read.csv(file="UF_weight.csv",sep=";", dec=",") #UniFrac weighted
rownames(unif_weig) <- unif_weig$X
unif_weight <- unif_weig[,-1]
UFw <- as.dist(unif_weight, diag=FALSE, upper=FALSE)

#braycurtis

adonis(betad1~tick_stage*masl*site*month, tick, perm=200)
adonis(betad1~tick_stage+masl*site*month+tick_stage:masl:site:month+tick_stage:masl:month+tick_stage:site:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~tick_stage+masl*site*month+tick_stage:masl:month+tick_stage:site:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~tick_stage+masl*site*month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+month+tick_stage+site:month+masl:site+masl:site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+tick_stage+site:month+masl:site+masl:site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+tick_stage+site:month+masl:site+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+tick_stage+site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+tick_stage+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+tick_stage+tick_stage:masl:month+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+tick_stage+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(betad1~masl+site+tick_stage+tick_stage:site, tick, perm=200) # THIS ONE
adonis(betad1~masl+site+tick_stage, tick, perm=200)
adonis(betad1~site+tick_stage, tick, perm=200)
adonis(betad1~tick_stage, tick, perm=200)

#Jaccard

adonis(betad2~tick_stage*masl*site*month, tick, perm=200)
adonis(betad2~tick_stage+masl+site*month+tick_stage:masl:month+tick_stage:site:month+tick_stage:masl+tick_stage:month+tick_stage:site+masl:site:month+masl:tick_stage+masl:site, tick, perm=200)
adonis(betad2~tick_stage+masadonis(UFw ~tick_stage*masl*site*month, tick, perm=200)
l+site*month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site+masl:site:month+masl:tick_stage+masl:site, tick, perm=200)
adonis(betad2~tick_stage+masl+site+site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site+masl:site:month+masl:tick_stage+masl:site, tick, perm=200)
adonis(betad2~tick_stage+masl+site+site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site+masl:site:month+masl:tick_stage, tick, perm=200)
adonis(betad2~tick_stage+masl+site+site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site+masl:tick_stage, tick, perm=200)
adonis(betad2~tick_stage+masl+site+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site+masl:tick_stage, tick, perm=200)
adonis(betad2~tick_stage+masl+site+tick_stage:masl+tick_stage:month+tick_stage:site+masl:tick_stage, tick, perm=200)
adonis(betad2~tick_stage+masl+site+tick_stage:masl+tick_stage:site+masl:tick_stage, tick, perm=200)
adonis(betad2~tick_stage+masl+site+tick_stage:masl, tick, perm=200)
adonis(betad2~tick_stage+masl+masl:tick_stage, tick, perm=200)
adonis(betad2~tick_stage+masl, tick, perm=200)

#Uni-Frac

adonis(UFw~tick_stage*masl*site*month, tick, perm=200)
adonis(UFw~tick_stage+masl*site*month+tick_stage:masl:site:month+tick_stage:masl:month+tick_stage:site:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~tick_stage+masl*site*month+tick_stage:masl:month+tick_stage:site:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~tick_stage+masl*site*month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+month+tick_stage+site:month+masl:site+masl:site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+tick_stage+site:month+masl:site+masl:site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+tick_stage+site:month+masl:site+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+tick_stage+site:month+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+tick_stage+tick_stage:masl:month+tick_stage:masl+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+tick_stage+tick_stage:masl:month+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+tick_stage+tick_stage:month+tick_stage:site, tick, perm=200)
adonis(UFw~masl+site+tick_stage+tick_stage:site, tick, perm=200) # THIS ONE
adonis(UFw~masl+site+tick_stage, tick, perm=200)
adonis(UFw~site+tick_stage, tick, perm=200)
adonis(UFw~tick_stage, tick, perm=200)

#FIGURE 3

mod <- metaMDS(otu_rare)
MDS1 <- as.numeric(mod$points[,1])
MDS2 <- as.numeric(mod$points[,2])
MDS <- cbind(MDS1, MDS2)
row.names(MDS) <- row.names(tick)
MDS <- data.frame(MDS)
MDS$group1 <- tick$tick_stage
MDS$masl <- as.numeric(tick$masl)
MDS$group2 <- tick$elev

MDS_group1 <- as.factor(MDS$group1)
MDS_group2 <- as.factor(MDS$group2)

find_hull <- function(MDS) MDS[chull(MDS$MDS1,MDS$MDS2), ]

hulls1 <- ddply(MDS, "MDS_group1", find_hull)
hulls2 <- ddply(MDS, "MDS_group2", find_hull)
hulls1$group1 <- as.factor(hulls1$group1)
hulls2$group2 <- as.factor(hulls2$group2)

#Fig 3a

ggplot(NULL, aes())+
  geom_polygon(data=hulls2, aes(MDS1,MDS2,fill=factor(group2)), alpha=0.3)+
  geom_point(data=MDS, aes(MDS1,MDS2,colour=group2,shape=group2), alpha=1, size=4)+
  scale_shape_manual(values=c(15,16,18))+
  labs(x="Axis 1",y="Axis 2")+
  theme(axis.text=element_text(colour="black"),
        legend.title=element_blank(),
        panel.background=element_rect(fill="white", colour="black"))

#Fig 3b 

ggplot(NULL, aes())+
  geom_polygon(data=hulls1, aes(MDS1,MDS2,fill=factor(group1)), alpha=0.3)+
  geom_point(data=MDS, aes(MDS1,MDS2,colour=group1,shape=group1), alpha=1, size=4)+
  scale_shape_manual(values=c(15,16,18))+
  labs(x="Axis 1",y="Axis 2")+
  theme(axis.text=element_text(colour="black"),
        legend.title=element_blank(),
        panel.background=element_rect(fill="white", colour="black"))


#Group dispersion

p1 <- rep(NA, 4)
mod <- with(tick, betadisper(betad1,masl))
anova(mod)
p1[1] <-  anova(mod)$Pr
mod <- with(tick, betadisper(betad1,site))
anova(mod)
p1[2] <-  anova(mod)$Pr
mod <- with(tick, betadisper(betad1,tick_stage))
anova(mod)
p1[3] <-  anova(mod)$Pr
mod <- with(tick, betadisper(betad1,month))
anova(mod)
p1[4] <-  anova(mod)$Pr
p1ad <- p.adjust(p1, "hochberg")

p2 <- rep(NA, 4)
mod <- with(tick, betadisper(betad2,masl))
anova(mod)
p2[1] <-  anova(mod)$Pr
mod <- with(tick, betadisper(betad2,site))
anova(mod)
p2[2] <-  anova(mod)$Pr
mod <- with(tick, betadisper(betad2,tick_stage))
anova(mod)
p2[3] <-  anova(mod)$Pr
mod <- with(tick, betadisper(betad2,month))
anova(mod)
p2[4] <-  anova(mod)$Pr
mod <- with(tick, betadisper(UFw,masl))
anova(mod)
p2ad <- p.adjust(p2, "hochberg")

p3 <- rep(NA, 4)
p3[1] <-  anova(mod)$Pr
mod <- with(tick, betadisper(UFw,site))
anova(mod)
p3[2] <-  anova(mod)$Pr
mod <- with(tick, betadisper(UFw,tick_stage))
anova(mod)
p3[3] <-  anova(mod)$Pr
mod <- with(tick, betadisper(UFw,month))
anova(mod)
p3[4] <-  anova(mod)$Pr
p3ad <- p.adjust(p3, "hochberg")

# Random forest

tick <- droplevels(tick)
perm = 1e4
set.seed(1)
forest <- randomForest(x = as.matrix(otu.relabund), y = as.factor(tick$site), importance = TRUE, 
                       proximity = TRUE, ntree = perm)  #run random forest

confusion <- forest[[5]][, 1:2]
oob.error <- (confusion[1, 2] + confusion[2, 1])/sum(confusion) * 100
oob.error
print(forest)

#picante analysis
# First phyloseq 

# Assign variables for imported data
sharedfile = "tick99.shared"
taxfile = "tick99.taxonomy"
treefile <- "tick99.tre" 
allfile <- "all_sample.csv"

# Import mothur data
tick_phyloseq <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile,
                             mothur_tree_file = treefile)

colnames(tax_table(tick_phyloseq)) <- c("Kingdom", "Phylum", "Class", 
                                       "Order", "Family", "Genus")

all <- read.delim(allfile)
all <- sample_data(all)
rownames(all) <- all$group
tick_phylose <- merge_phyloseq(tick_phyloseq, all)

#remove samples removed in earlier analysis

tick_phylos = subset_samples(tick_phylose, group != "F10" & group != "F38" & group != "F36" & 
                                group != "F5" & group != "F50" & group != "F8" & group != "H15" & 
                                group != "M20" & group != "M25" & group != "M3" & 
                                group != "H32" & group != "M33" & group != "M8" & group != "M9")


tick_ps <- sample_data(tick)
rownames(tick_ps) <- tick_ps$group
tick_phylo <- merge_phyloseq(tick_phylos, tick_ps)
tick_phyl <- rarefy_even_depth(tick_phylo, sample.size=500)
tick_phy <- filter_taxa(tick_phyl, function(x) sum(x>0) > (0.02*length(x)), TRUE) 


#community and phylogeny files

picOTU <- function(physeq) {
  require("picante")
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

comm <- picOTU(tick_phy)
phy <- tick_phy@phy_tree

phydist <- cophenetic(phy)

#NTI and NRI

ses.mpd.result <- ses.mpd(comm, phydist, null.model="taxa.labels")
ses.mntd.result <- ses.mntd(comm, phydist, null.model="taxa.labels")

nti_nri <- ses.mntd.result$mntd.obs.z/ses.mpd.result$mpd.obs.z #Phylogenetic clustering!

alpha_phy <- data.frame(tick$elev,tick$masl,tick$tick_stage,tick$site, tick$month, ses.mpd.result$mpd.obs.z,
                                     ses.mntd.result$mntd.obs.z,nti_nri)
names(alpha_phy) <- c('elevation', 'masl','stage','site','month','mpd', 'mntd','nti_nri')

summary(lm(nti_nri~elevation+stage+month, data=alpha_phy))

m3.nlme = lme(mpd ~ masl*stage*month+(I(masl^2)), random= ~1|site,  data=alpha_phy)
m3.nlme = lme(mpd ~ masl+month+stage+masl:month+stage:month+ masl:stage+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ masl+month+stage+masl:month+ masl:stage+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ masl+month+stage+masl:stage+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ masl+stage+month+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ stage+month+(I(month^2)), random= ~1|site,  data=alpha_phy)
m3.nlme = lme(mpd ~ stage+month, random= ~1|site,  data=alpha_phy)
m3.nlme = lme(mpd ~ stage, random= ~1|site,  data=alpha_phy)
m3.nlme = lme(mpd ~ 1, random= ~1|site,  data=alpha_phy)

summary(m3.nlme)
r.squaredGLMM(m3.nlme)

m4.nlme = lme(mntd ~ masl*stage*month+(I(masl^2)), random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ masl*stage*month, random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ masl+month+stage+masl:month+stage:month+ masl:stage, random= ~1|site, data=alpha_phy)
m4.nlme = lme(mntd ~ masl+month+stage+masl:month+masl:stage, random= ~1|site, data=alpha_phy)
m4.nlme = lme(mntd ~ masl+month+stage+masl:stage, random= ~1|site, data=alpha_phy)
m4.nlme = lme(mntd ~ masl+stage+month, random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ stage+month, random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ stage, random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ 1, random= ~1|site,  data=alpha_phy)

summary(m4.nlme)
r.squaredGLMM(m4.nlme)


#FIGURE 3C-D
p <- ggplot(alpha_phy, aes(x=stage,y=mpd,color=stage))+
  geom_jitter(position=position_jitter(0.1))+
  stat_summary(fun=median, geom="crossbar", ymax=-100, ymin=-100)+
  scale_colour_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  theme_bw(base_size=22) 
print(p)

p <- ggplot(alpha_phy, aes(x=stage,y=mntd,color=stage))+
  geom_jitter(position=position_jitter(0.1))+
  stat_summary(fun=median, geom="crossbar", ymax=-100, ymin=-100)+
  scale_colour_manual(values=c("#66c2a5","#fc8d62","#8da0cb"))+
  theme_bw(base_size=22) 
print(p)

low <- tick$elev=="low"
middle <- tick$elev=="middle"
high <- tick$elev=="high"

nti_nri[high] 
nti_nri[middle]
nti_nri[low]

describe(nti_nri)
describe(nti_nri[low])
describe(nti_nri[middle])
describe(nti_nri[high])

##Phylosymbiosis

#LOESS

loesstick <- loess(UFw ~ticdist, span=0.2)
smoothed <- predict(loesstick)

plot(ticdist,UFw)
lines(smoothed, x=ticdist)

#pseudoR2

ss.dist <- sum(scale(UFw, scale=FALSE)^2)
ss.resid <- sum(resid(loesstick)^2)
1-ss.resid/ss.dist

# Mantel and host distances

mtest1 <- mantel.rtest(UFw, ticdist,nrepet = 9999) #No significance

# Trees and host genetics

tick_tree <- ticdist %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade

comm_tree <- UFw %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade


comm_tree$tip.label <- tick_tree$tip.label
MatchingSplitDistance(tick_tree, comm_tree)

distances <- integer(10000)


for (i in 1:10000) {
  random <- rtree(79, tip.label=tick_tree$tip.label)
  distances[i] <- MatchingSplitDistance(random,tick_tree)
}
hist(distances) #No significance



##Phylogenetic community turnover
# Phylogenetic turnover

#First, load up RC-code from Stegen et al. 2013

devtools::source_url("https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r?raw=TRUE")


#NB. These are computationally very heavy 

#sesbetaNTI <- ses.comdistnt(comm, phydist, null.model="taxa.labels", abundance.weighted=TRUE)
#RCbray <- raup_crick_abundance(otu_rare, plot_names_in_col1 = FALSE, reps=50)

#sesbetaNTI <- read.csv(file="sesBETA.csv",dec=".", stringsAsFactors = FALSE)
#rownames(sesbetaNTI) <- colnames(sesbetaNTI) 

#RCbray2 <- read.csv(file="RCbray2_csc.csv",dec=".")
#rownames(RCbray2) <- colnames(RCbray2) 

quantile(sesbetaNTI, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2, c(0.25, 0.50, 0.75), na.rm=TRUE) 

## All samples
## classification: 
# Variable selection      betaNTI > 2

class_VS <- sesbetaNTI>=2
sum(class_VS, na.rm=TRUE)/2 #59

# Homogenous selection:   betaNTI < -2

class_HS <- -2>=sesbetaNTI
sum(class_HS, na.rm=TRUE)/2 #273

# Dispersal limitation    betaNTI < 2, RCbray > 0.95

class_DLtemp <- ((sesbetaNTI<2 & sesbetaNTI>(-2)) + (RCbray2>=0.95))
class_DL <- (class_DLtemp==2)
sum(class_DL, na.rm=TRUE)/2 #280


# Homogenizing dispersal  betaNTI < 2, RCbray < -0.95

class_HDtemp <- ((sesbetaNTI<2 & sesbetaNTI>(-2)) + (-0.95>=RCbray2))
class_HD <- (class_HDtemp==2)
sum(class_HD, na.rm=TRUE)/2 #44

# Undominated             betaNTI < 2, RCbray < 0.95

class_UDtemp <- ((sesbetaNTI<2 & sesbetaNTI>(-2)) + ((RCbray2<0.95)&(RCbray2>(-0.95))))
class_UD <- (class_UDtemp==2)
sum(class_UD, na.rm=TRUE)/2 #2425




#SITES
#Select sites


FEL <- tick$elev=="low"&tick$site=="FE"
FEM <- tick$elev=="middle"&tick$site=="FE"
FLL <- tick$elev=="low"&tick$site=="FL"
FLM <- tick$elev=="middle"&tick$site=="FL"
PAL <- tick$elev=="low"&tick$site=="PA"
PAM <- tick$elev=="middle"&tick$site=="PA"

sesbetaNTI_FEL <- sesbetaNTI[FEL,FEL]
sesbetaNTI_FEM <- sesbetaNTI[FEM,FEM]
sesbetaNTI_FLL <- sesbetaNTI[FLL,FLL]
sesbetaNTI_FLM <- sesbetaNTI[FLM,FLM]
sesbetaNTI_PAL <- sesbetaNTI[PAL,PAL]
sesbetaNTI_PAM <- sesbetaNTI[PAM,PAM]

RCbray2_FEL <- RCbray2[FEL,FEL]
RCbray2_FEM <- RCbray2[FEM,FEM]
RCbray2_FLL <- RCbray2[FLL,FLL]
RCbray2_FLM <- RCbray2[FLM,FLM]
RCbray2_PAL <- RCbray2[PAL,PAL]
RCbray2_PAM <- RCbray2[PAM,PAM]



quantile(sesbetaNTI_FEL, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_FEL, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_FEM, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_FEM, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_FLL, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_FLL, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_FLM, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_FLM, c(0.25, 0.50, 0.75), na.rm=TRUE) 

quantile(sesbetaNTI_PAL, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_PAL, c(0.25, 0.50, 0.75), na.rm=TRUE) 

quantile(sesbetaNTI_PAM, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_PAM, c(0.25, 0.50, 0.75), na.rm=TRUE) 


## classification: 
# Variable selection      betaNTI > 2

class_VSFEL <- sesbetaNTI_FEL>=2
sum(class_VSFEL, na.rm=TRUE)/2 #6

class_VSFEM <- sesbetaNTI_FEM>=2
sum(class_VSFEM, na.rm=TRUE)/2 #5

class_VSFLL <- sesbetaNTI_FLL>=2
sum(class_VSFLL, na.rm=TRUE)/2 #8

class_VSFLM <- sesbetaNTI_FLM>=2
sum(class_VSFLM, na.rm=TRUE)/2 #0

class_VSPAL <- sesbetaNTI_PAL>=2
sum(class_VSPAL, na.rm=TRUE)/2 #0

class_VSPAM <- sesbetaNTI_PAM>=2
sum(class_VSPAM, na.rm=TRUE)/2 #0

# Homogenous selection:   betaNTI < -2

class_HSFEL <- sesbetaNTI_FEL<=-2
sum(class_HSFEL, na.rm=TRUE)/2 #6

class_HSFEM <- sesbetaNTI_FEM<=-2
sum(class_HSFEM, na.rm=TRUE)/2 #5

class_HSFLL <- sesbetaNTI_FLL<=-2
sum(class_HSFLL, na.rm=TRUE)/2 #11

class_HSFLM <- sesbetaNTI_FLM<=-2
sum(class_HSFLM, na.rm=TRUE)/2 #5

class_HSPAL <- sesbetaNTI_PAL<=-2
sum(class_HSPAL, na.rm=TRUE)/2 #6

class_HSPAM <- sesbetaNTI_PAM<=-2
sum(class_HSPAM, na.rm=TRUE)/2 #3

# Dispersal limitation    betaNTI < 2, RCbray > 0.95

class_DLtemp <- ((sesbetaNTI_FEL<2 & sesbetaNTI_FEL>(-2)) + (RCbray2_FEL>=0.95))
class_DLFEL <- (class_DLtemp==2)
sum(class_DLFEL, na.rm=TRUE)/2 #9

class_DLtemp <- ((sesbetaNTI_FEM<2 & sesbetaNTI_FEM>(-2)) + (RCbray2_FEM>=0.95))
class_DLFEM <- (class_DLtemp==2)
sum(class_DLFEM, na.rm=TRUE)/2 #4

class_DLtemp <- ((sesbetaNTI_FLL<2 & sesbetaNTI_FLL>(-2)) + (RCbray2_FLL>=0.95))
class_DLFLL <- (class_DLtemp==2)
sum(class_DLFLL, na.rm=TRUE)/2 #13

class_DLtemp <- ((sesbetaNTI_FLM<2 & sesbetaNTI_FLM>(-2)) + (RCbray2_FLM>=0.95))
class_DLFLM <- (class_DLtemp==2)
sum(class_DLFLM, na.rm=TRUE)/2 #10

class_DLtemp <- ((sesbetaNTI_PAL<2 & sesbetaNTI_PAL>(-2)) + (RCbray2_PAL>=0.95))
class_DLPAL <- (class_DLtemp==2)
sum(class_DLPAL, na.rm=TRUE)/2 #1

class_DLtemp <- ((sesbetaNTI_PAM<2 & sesbetaNTI_PAM>(-2)) + (RCbray2_PAM>=0.95))
class_DLPAM <- (class_DLtemp==2)
sum(class_DLPAM, na.rm=TRUE)/2 #0

# Homogenizing dispersal  betaNTI < 2, RCbray < -0.95

class_HDtemp <- ((sesbetaNTI_FEL<2 & sesbetaNTI_FEL>(-2)) + (-0.95>=RCbray2_FEL))
class_HDFEL <- (class_HDtemp==2)
sum(class_HDFEL, na.rm=TRUE)/2 #4

class_HDtemp <- ((sesbetaNTI_FEM<2 & sesbetaNTI_FEM>(-2)) + (-0.95>=RCbray2_FEM))
class_HDFEM <- (class_HDtemp==2)
sum(class_HDFEM, na.rm=TRUE)/2 #1

class_HDtemp <- ((sesbetaNTI_FLL<2 & sesbetaNTI_FLL>(-2)) + (-0.95>=RCbray2_FLL))
class_HDFLL <- (class_HDtemp==2)
sum(class_HDFLL, na.rm=TRUE)/2 #2

class_HDtemp <- ((sesbetaNTI_FLM<2 & sesbetaNTI_FLM>(-2)) + (-0.95>=RCbray2_FLM))
class_HDFLM <- (class_HDtemp==2)
sum(class_HDFLM, na.rm=TRUE)/2 #0

class_HDtemp <- ((sesbetaNTI_PAL<2 & sesbetaNTI_PAL>(-2)) + (-0.95>=RCbray2_PAL))
class_HDPAL <- (class_HDtemp==2)
sum(class_HDPAL, na.rm=TRUE)/2 #1

class_HDtemp <- ((sesbetaNTI_PAM<2 & sesbetaNTI_PAM>(-2)) + (-0.95>=RCbray2_PAM))
class_HDPAM <- (class_HDtemp==2)
sum(class_HDPAM, na.rm=TRUE)/2 #0

# Undominated             betaNTI < 2, RCbray < 0.95

class_UDtemp <- ((sesbetaNTI_FEL<2 & sesbetaNTI_FEL>(-2)) + ((RCbray2_FEL<0.95)&(RCbray2_FEL>(-0.95))))
class_UDFEL <- (class_UDtemp==2)
sum(class_UDFEL, na.rm=TRUE)/2 #72

class_UDtemp <- ((sesbetaNTI_FEM<2 & sesbetaNTI_FEM>(-2)) + ((RCbray2_FEM<0.95)&(RCbray2_FEM>(-0.95))))
class_UDFEM <- (class_UDtemp==2)
sum(class_UDFEM, na.rm=TRUE)/2 #53

class_UDtemp <- ((sesbetaNTI_FLL<2 & sesbetaNTI_FLL>(-2)) + ((RCbray2_FLL<0.95)&(RCbray2_FLL>(-0.95))))
class_UDFLL <- (class_UDtemp==2)
sum(class_UDFLL, na.rm=TRUE)/2 #176

class_UDtemp <- ((sesbetaNTI_FLM<2 & sesbetaNTI_FLM>(-2)) + ((RCbray2_FLM<0.95)&(RCbray2_FLM>(-0.95))))
class_UDFLM <- (class_UDtemp==2)
sum(class_UDFLM, na.rm=TRUE)/2 #40

class_UDtemp <- ((sesbetaNTI_PAL<2 & sesbetaNTI_PAL>(-2)) + ((RCbray2_PAL<0.95)&(RCbray2_PAL>(-0.95))))
class_UDPAL <- (class_UDtemp==2)
sum(class_UDPAL, na.rm=TRUE)/2 #37

class_UDtemp <- ((sesbetaNTI_PAM<2 & sesbetaNTI_PAM>(-2)) + ((RCbray2_PAM<0.95)&(RCbray2_PAM>(-0.95))))
class_UDPAM <- (class_UDtemp==2)
sum(class_UDPAM, na.rm=TRUE)/2 #12



##TICK STAGE by SITE


FELF <- tick$elev=="low"&tick$site=="FE"&tick$tick_stage=="F"
FELM <- tick$elev=="low"&tick$site=="FE"&tick$tick_stage=="M"

FEMF <- tick$elev=="middle"&tick$site=="FE"&tick$tick_stage=="F"
FEMM <- tick$elev=="middle"&tick$site=="FE"&tick$tick_stage=="M"
FEMN <- tick$elev=="middle"&tick$site=="FE"&tick$tick_stage=="N"

FLLF <- tick$elev=="low"&tick$site=="FL"&tick$tick_stage=="F"
FLLM <- tick$elev=="low"&tick$site=="FL"&tick$tick_stage=="M"

FLMF <- tick$elev=="middle"&tick$site=="FL"&tick$tick_stage=="F"
FLMM <- tick$elev=="middle"&tick$site=="FL"&tick$tick_stage=="M"
FLMN <- tick$elev=="middle"&tick$site=="FL"&tick$tick_stage=="N"

PALF <- tick$elev=="low"&tick$site=="PA"&tick$tick_stage=="F"
PALM <- tick$elev=="low"&tick$site=="PA"&tick$tick_stage=="M"

PAMF <- tick$elev=="middle"&tick$site=="PA"&tick$tick_stage=="F"
PAMM <- tick$elev=="middle"&tick$site=="PA"&tick$tick_stage=="M"


sesbetaNTI_FELM <- sesbetaNTI[FELM,FELM]
sesbetaNTI_FELF <- sesbetaNTI[FELF,FELF]

sesbetaNTI_FEMM <- sesbetaNTI[FEMM,FEMM]
sesbetaNTI_FEMF <- sesbetaNTI[FEMF,FEMF]
sesbetaNTI_FEMN <- sesbetaNTI[FEMN,FEMN]

sesbetaNTI_FLLM <- sesbetaNTI[FLLM,FLLM]
sesbetaNTI_FLLF <- sesbetaNTI[FLLF,FLLF]

sesbetaNTI_FLMM <- sesbetaNTI[FLMM,FLMM]
sesbetaNTI_FLMF <- sesbetaNTI[FLMF,FLMF]
sesbetaNTI_FLMN <- sesbetaNTI[FLMN,FLMN]

sesbetaNTI_PALM <- sesbetaNTI[PALM,PALM]
sesbetaNTI_PALF <- sesbetaNTI[PALF,PALF]

sesbetaNTI_PAMM <- sesbetaNTI[PAMM,PAMM]
sesbetaNTI_PAMF <- sesbetaNTI[PAMF,PAMF]

RCbray2_FELM <- RCbray2[FELM,FELM]
RCbray2_FELF <- RCbray2[FELF,FELF]

RCbray2_FEMM <- RCbray2[FEMM,FEMM]
RCbray2_FEMF <- RCbray2[FEMF,FEMF]
RCbray2_FEMN <- RCbray2[FEMN,FEMN]

RCbray2_FLLM <- RCbray2[FLLM,FLLM]
RCbray2_FLLF <- RCbray2[FLLF,FLLF]

RCbray2_FLMM <- RCbray2[FLMM,FLMM]
RCbray2_FLMF <- RCbray2[FLMF,FLMF]
RCbray2_FLMN <- RCbray2[FLMN,FLMN]

RCbray2_PALM <- RCbray2[PALM,PALM]
RCbray2_PALF <- RCbray2[PALF,PALF]

RCbray2_PAMM <- RCbray2[PAMM,PAMM]
RCbray2_PAMF <- RCbray2[PAMF,PAMF]

WITHIN <- FELF+FELM+FEMM+FEMN+FLLF+FLLM+FLMF+FLMM+FLMN+PALF+PALM+PAMF+PAMM
WITH=ifelse(WITHIN=="1",TRUE,FALSE)
sesbetaNTI_WITH <- sesbetaNTI[WITH,WITH]
RCbray2_WITH <- RCbray2[WITH,WITH]


#For between-site comparisons

class_VS_WITH <- sesbetaNTI_WITH>=2
sum(class_VS_WITH, na.rm=TRUE)/2 #49

class_HS_WITH <- -2>=sesbetaNTI_WITH
sum(class_HS_WITH, na.rm=TRUE)/2 #193

class_DLtemp <- ((sesbetaNTI_WITH<2 & sesbetaNTI_WITH>(-2)) + (RCbray2_WITH>=0.95))
class_DL_WITH <- (class_DLtemp==2)
sum(class_DL_WITH, na.rm=TRUE)/2 #206

class_HDtemp <- ((sesbetaNTI_WITH<2 & sesbetaNTI_WITH>(-2)) + (-0.95>=RCbray2_WITH))
class_HD_WITH <- (class_HDtemp==2)
sum(class_HD_WITH, na.rm=TRUE)/2 #35

class_UDtemp <- ((sesbetaNTI_WITH<2 & sesbetaNTI_WITH>(-2)) + ((RCbray2_WITH<0.95)&(RCbray2_WITH>(-0.95))))
class_UD_WITH <- (class_UDtemp==2)
sum(class_UD_WITH, na.rm=TRUE)/2 #1795

## classification: 
# Variable selection      betaNTI > 2

class_VS_FELM <- sesbetaNTI_FELM>=2
sum(class_VS_FELM, na.rm=TRUE)/2 #0 ... These were counted then manually - see males and females down below.


## On to total numbers and Figure 5/Table S3

#column All

all <- c(sum(class_VS, na.rm=TRUE)/2,sum(class_HS, na.rm=TRUE)/2, sum(class_DL, na.rm=TRUE)/2, 
         sum(class_HD, na.rm=TRUE)/2,sum(class_UD, na.rm=TRUE)/2)


#column Within-site

elevL<- c(sum(class_VSFEL, na.rm=TRUE)/2+sum(class_VSFLL, na.rm=TRUE)/2+sum(class_VSPAL, na.rm=TRUE)/2,
          sum(class_HSFEL, na.rm=TRUE)/2+sum(class_HSFLL, na.rm=TRUE)/2+sum(class_HSPAL, na.rm=TRUE)/2,
          sum(class_DLFEL, na.rm=TRUE)/2+sum(class_DLFLL, na.rm=TRUE)/2+sum(class_DLPAL, na.rm=TRUE)/2,
          sum(class_HDFEL, na.rm=TRUE)/2+sum(class_HDFLL, na.rm=TRUE)/2+sum(class_HDPAL, na.rm=TRUE)/2,
          sum(class_UDFEL, na.rm=TRUE)/2+sum(class_UDFLL, na.rm=TRUE)/2+sum(class_UDPAL, na.rm=TRUE)/2)
elevM<- c(sum(class_VSFEM, na.rm=TRUE)/2+sum(class_VSFLM, na.rm=TRUE)/2+sum(class_VSPAM, na.rm=TRUE)/2,
          sum(class_HSFEM, na.rm=TRUE)/2+sum(class_HSFLM, na.rm=TRUE)/2+sum(class_HSPAM, na.rm=TRUE)/2,
          sum(class_DLFEM, na.rm=TRUE)/2+sum(class_DLFLM, na.rm=TRUE)/2+sum(class_DLPAM, na.rm=TRUE)/2,
          sum(class_HDFEM, na.rm=TRUE)/2+sum(class_HDFLM, na.rm=TRUE)/2+sum(class_HDPAM, na.rm=TRUE)/2,
          sum(class_UDFEM, na.rm=TRUE)/2+sum(class_UDFLM, na.rm=TRUE)/2+sum(class_UDPAM, na.rm=TRUE)/2)
total<- elevL+elevM

#Column Between-site

betwsite <- c(sum(class_VS_WITH, na.rm=TRUE)/2-total[1],sum(class_HS_WITH, na.rm=TRUE)/2-total[2],
              sum(class_DL_WITH, na.rm=TRUE)/2-total[3],sum(class_HD_WITH, na.rm=TRUE)/2-total[4],
              sum(class_UD_WITH, na.rm=TRUE)/2-total[5])
 
#Tick sex

female <- c(2,13,4,3,178)
male <- c(4,2,2,2,63)

#Collating data

group <- c(rep("All" , 5), rep("Between Sites" , 5) , rep("Total" , 5) , rep("Low" , 5), 
           rep("Middle" , 5), rep("Female" , 5) , rep("Male" , 5))
group <- factor(group, levels = c("All","Between Sites","Total","Low","Middle","Female","Male"))


proc <- rep(c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),7)
proc <- factor(proc, levels = c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"))

value <- c(all,betwsite,total,elevL,elevM,female,male)

data <- data.frame(group,proc,value)


#Chi-squared tests

test <- as.table(cbind(c(value[6:10]),c(value[11:15])))
dimnames(test)<- list(factor= c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),
                      level= c("All","Between Sites","Total","Low","Middle","Female","Male"))
chisq.test(test)

test <- as.table(cbind(c(value[16:20]),c(value[21:25])))
dimnames(test)<- list(factor= c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),
                      level= c("All","Between Sites","Total","Low","Middle","Female","Male"))
chisq.test(test)

test <- as.table(cbind(c(value[26:30]),c(value[31:35])))
dimnames(test)<- list(factor= c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),
                      level= c("All","Between Sites","Total","Low","Middle","Female","Male"))
chisq.test(test)

#Figure 5

ggplot(data, aes(fill=proc, y=value, x=group)) + 
  geom_bar(position="fill", stat="identity")