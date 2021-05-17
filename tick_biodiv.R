
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
m1.nlme = lme(simpson ~ masl*tick_stage*month+(I(month^2)), random= ~1|site,  data=tick)
m1.nlme = lme(simpson ~ masl*tick_stage*month, random= ~1|site,  data=tick)
m1.nlme = lme(simpson ~ masl+tick_stage+month+masl:tick_stage+masl:month+tick_stage:month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+masl:month+tick_stage*month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+masl:month+tick_stage+tick_stage:month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+tick_stage+tick_stage:month, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ masl+tick_stage, random= ~1|site, data=tick)
m1.nlme = lme(simpson ~ tick_stage, random= ~1|site, data=tick)


summary(m1.nlme)
r.squaredGLMM(m1.nlme)

m2.nlme = lme(phylodiv ~ masl*tick_stage*month+(I(month^2)), random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl*tick_stage*month, random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl+tick_stage+month+masl:tick_stage+masl:month+tick_stage:month, random= ~1|site, data=tick)
m2.nlme = lme(phylodiv ~ masl+masl:month+tick_stage*month, random= ~1|site, data=tick)
m2.nlme = lme(phylodiv ~ masl+month+tick_stage+tick_stage:month, random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl+tick_stage+ month , random= ~1|site,  data=tick)
m2.nlme = lme(phylodiv ~ masl+month, random= ~1|site,  data=tick)

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

m3.nlme = lme(mpd ~ masl*stage*month+(I(month^2)), random= ~1|site,  data=alpha_phy)
m3.nlme = lme(mpd ~ masl+month+stage+masl:month+stage:month+ masl:stage+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ masl+month+stage+masl:month+ masl:stage+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ masl+month+stage+masl:stage+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ masl+stage+month+(I(month^2)), random= ~1|site, data=alpha_phy)
m3.nlme = lme(mpd ~ stage+month+(I(month^2)), random= ~1|site,  data=alpha_phy)
m3.nlme = lme(mpd ~ stage+month, random= ~1|site,  data=alpha_phy)
m3.nlme = lme(mpd ~ stage, random= ~1|site,  data=alpha_phy)
summary(m3.nlme)
r.squaredGLMM(m3.nlme)

m4.nlme = lme(mntd ~ masl*stage*month+(I(month^2)), random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ masl*stage*month, random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ masl+month+stage+masl:month+stage:month+ masl:stage, random= ~1|site, data=alpha_phy)
m4.nlme = lme(mntd ~ masl+month+stage+masl:month+masl:stage, random= ~1|site, data=alpha_phy)
m4.nlme = lme(mntd ~ masl+month+stage+masl:stage, random= ~1|site, data=alpha_phy)
m4.nlme = lme(mntd ~ masl+stage+month, random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ stage+month, random= ~1|site,  data=alpha_phy)
m4.nlme = lme(mntd ~ stage, random= ~1|site,  data=alpha_phy)
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

# On ELEVATIION
## select H, M and L groups

low <- tick$elev=="low"
middle <- tick$elev=="middle"
high <- tick$elev=="high"

sesbetaNTI_H <- sesbetaNTI[high,high]
sesbetaNTI_M <- sesbetaNTI[middle,middle]
sesbetaNTI_L <- sesbetaNTI[low,low]

RCbray2_H <- RCbray2[high,high]
RCbray2_M <- RCbray2[middle,middle]
RCbray2_L <- RCbray2[low,low]

quantile(sesbetaNTI_H, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_H, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_M, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_M, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_L, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_L, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_H, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_H, c(0.25, 0.50, 0.75), na.rm=TRUE) 

quantile(sesbetaNTI_M, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_M, c(0.25, 0.50, 0.75), na.rm=TRUE) 

quantile(sesbetaNTI_L, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_L, c(0.25, 0.50, 0.75), na.rm=TRUE) 


## classification: 
# Variable selection      betaNTI > 2

class_VSH <- sesbetaNTI_H>=2
sum(class_VSH, na.rm=TRUE)/2 #0

class_VSM <- sesbetaNTI_M>=2
sum(class_VSM, na.rm=TRUE)/2 #11

class_VSL <- sesbetaNTI_L>=2
sum(class_VSL, na.rm=TRUE)/2 #15

# Homogenous selection:   betaNTI < -2

class_HSH <- -2>=sesbetaNTI_H
sum(class_HSH, na.rm=TRUE)/2 #1

class_HSM <- -2>=sesbetaNTI_M
sum(class_HSM, na.rm=TRUE)/2 #47

class_HSL <- -2>=sesbetaNTI_L
sum(class_HSL, na.rm=TRUE)/2 #76

# Dispersal limitation    betaNTI < 2, RCbray > 0.95

class_DLtemp <- ((sesbetaNTI_H<2 & sesbetaNTI_H>(-2)) + (RCbray2_H>=0.95))
class_DLH <- (class_DLtemp==2)
sum(class_DLH, na.rm=TRUE)/2 #0

class_DLtemp <- ((sesbetaNTI_M<2 & sesbetaNTI_M>(-2)) + (RCbray2_M>=0.95))
class_DLM <- (class_DLtemp==2)
sum(class_DLM, na.rm=TRUE)/2 #44

class_DLtemp <- ((sesbetaNTI_L<2 & sesbetaNTI_L>(-2)) + (RCbray2_L>=0.95))
class_DLL <- (class_DLtemp==2)
sum(class_DLL, na.rm=TRUE)/2 #83


# Homogenizing dispersal  betaNTI < 2, RCbray < -0.95

class_HDtemp <- ((sesbetaNTI_H<2 & sesbetaNTI_H>(-2)) + (-0.95>=RCbray2_H))
class_HDH <- (class_HDtemp==2)
sum(class_HDH, na.rm=TRUE)/2 #0

class_HDtemp <- ((sesbetaNTI_M<2 & sesbetaNTI_M>(-2)) + (-0.95>=RCbray2_M))
class_HDM <- (class_HDtemp==2)
sum(class_HDM, na.rm=TRUE)/2 #4

class_HDtemp <- ((sesbetaNTI_L<2 & sesbetaNTI_L>(-2)) + (-0.95>=RCbray2_L))
class_HDL <- (class_HDtemp==2)
sum(class_HDL, na.rm=TRUE)/2 #20

# Undominated             betaNTI < 2, RCbray < 0.95

class_UDtemp <- ((sesbetaNTI_H<2 & sesbetaNTI_H>(-2)) + ((RCbray2_H<0.95)&(RCbray2_H>(-0.95))))
class_UDH <- (class_UDtemp==2)
sum(class_UDH, na.rm=TRUE)/2 #9

class_UDtemp <- ((sesbetaNTI_M<2 & sesbetaNTI_M>(-2)) + ((RCbray2_M<0.95)&(RCbray2_M>(-0.95))))
class_UDM <- (class_UDtemp==2)
sum(class_UDM, na.rm=TRUE)/2 #300

class_UDtemp <- ((sesbetaNTI_L<2 & sesbetaNTI_L>(-2)) + ((RCbray2_L<0.95)&(RCbray2_L>(-0.95))))
class_UDL <- (class_UDtemp==2)
sum(class_UDL, na.rm=TRUE)/2 #796

all <- c(sum(class_VS, na.rm=TRUE),sum(class_HS, na.rm=TRUE)/2, sum(class_DL, na.rm=TRUE)/2, 
         sum(class_HD, na.rm=TRUE)/2,sum(class_UD, na.rm=TRUE)/2)
elevL<- c(sum(class_VSL, na.rm=TRUE),sum(class_HSL, na.rm=TRUE)/2, sum(class_DLL, na.rm=TRUE)/2, 
          sum(class_HDL, na.rm=TRUE)/2,sum(class_UDL, na.rm=TRUE)/2)
elevM<- c(sum(class_VSM, na.rm=TRUE),sum(class_HSM, na.rm=TRUE)/2, sum(class_DLM, na.rm=TRUE)/2, 
          sum(class_HDM, na.rm=TRUE)/2,sum(class_UDM, na.rm=TRUE)/2)
elevH<- c(sum(class_VSH, na.rm=TRUE),sum(class_HSH, na.rm=TRUE)/2, sum(class_DLH, na.rm=TRUE)/2, 
          sum(class_HDH, na.rm=TRUE)/2,sum(class_UDH, na.rm=TRUE)/2)

elevation <- c(rep("all" , 5) , rep("low" , 5) , rep("middle" , 5) , rep("high" , 5) )
proc <- rep(c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),4)
value <- c(all,elevL,elevM,elevH)
data <- data.frame(elevation,proc,value)

#Figure in Supplement

ggplot(data, aes(fill=proc, y=value, x=elevation)) + 
  geom_bar(position="fill", stat="identity")
proc <- rep(c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),4)

# On TICK STAGE
## select M, F and N groups

M <- tick$tick_stage=="M"
F <- tick$tick_stage=="F"
N <- tick$tick_stage=="N"

sesbetaNTI_M <- sesbetaNTI[M,M]
sesbetaNTI_F <- sesbetaNTI[F,F]
sesbetaNTI_N <- sesbetaNTI[N,N]

RCbray2_M <- RCbray2[M,M]
RCbray2_F <- RCbray2[F,F]
RCbray2_N <- RCbray2[N,N]

quantile(sesbetaNTI_M, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_M, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_F, c(0.025, 0.50, 0.975), na.rm=TRUE) 
quantile(RCbray2_F, c(0.025, 0.50, 0.975), na.rm=TRUE) 

quantile(sesbetaNTI_N, c(0.025, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_N, c(0.025, 0.50, 0.75), na.rm=TRUE) 

quantile(sesbetaNTI_M, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_M, c(0.25, 0.50, 0.75), na.rm=TRUE) 

quantile(sesbetaNTI_F, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_F, c(0.25, 0.50, 0.75), na.rm=TRUE) 

quantile(sesbetaNTI_N, c(0.25, 0.50, 0.75), na.rm=TRUE) 
quantile(RCbray2_N, c(0.25, 0.50, 0.75), na.rm=TRUE) 



## classification: 
# Variable selection      betaNTI > 2

class_VS_M <- sesbetaNTI_M>=2
sum(class_VS_M, na.rm=TRUE)/2 #23

class_VS_F <- sesbetaNTI_F>=2
sum(class_VS_F, na.rm=TRUE)/2 #3

class_VS_N <- sesbetaNTI_N>=2
sum(class_VS_N, na.rm=TRUE)/2 #3

# Homogenous selection:   betaNTI < -2

class_HS_M <- -2>=sesbetaNTI_M
sum(class_HS_M, na.rm=TRUE)/2 #40

class_HS_F <- -2>=sesbetaNTI_F
sum(class_HS_F, na.rm=TRUE)/2 #98

class_HS_N <- -2>=sesbetaNTI_N
sum(class_HS_N, na.rm=TRUE)/2 #0

# Dispersal limitation    betaNTI < 2, RCbray > 0.95

class_DLtemp <- ((sesbetaNTI_M<2 & sesbetaNTI_M>(-2)) + (RCbray2_M>=0.95))
class_DL_M <- (class_DLtemp==2)
sum(class_DL_M, na.rm=TRUE)/2 #27

class_DLtemp <- ((sesbetaNTI_F<2 & sesbetaNTI_F>(-2)) + (RCbray2_F>=0.95))
class_DL_F <- (class_DLtemp==2)
sum(class_DL_F, na.rm=TRUE)/2 #28

class_DLtemp <- ((sesbetaNTI_N<2 & sesbetaNTI_N>(-2)) + (RCbray2_N>=0.95))
class_DL_N <- (class_DLtemp==2)
sum(class_DL_N, na.rm=TRUE)/2 #15


# Homogenizing dispersal  betaNTI < 2, RCbray < -0.95

class_HDtemp <- ((sesbetaNTI_M<2 & sesbetaNTI_M>(-2)) + (-0.95>=RCbray2_M))
class_HD_M <- (class_HDtemp==2)
sum(class_HD_M, na.rm=TRUE)/2 #6

class_HDtemp <- ((sesbetaNTI_F<2 & sesbetaNTI_F>(-2)) + (-0.95>=RCbray2_F))
class_HD_F <- (class_HDtemp==2)
sum(class_HD_F, na.rm=TRUE)/2 #25

class_HDtemp <- ((sesbetaNTI_N<2 & sesbetaNTI_N>(-2)) + (-0.95>=RCbray2_N))
class_HD_N <- (class_HDtemp==2)
sum(class_HD_N, na.rm=TRUE)/2 #1

# Undominated             betaNTI < 2, RCbray < 0.95

class_UDtemp <- ((sesbetaNTI_M<2 & sesbetaNTI_M>(-2)) + ((RCbray2_M<0.95)&(RCbray2_M>(-0.95))))
class_UD_M <- (class_UDtemp==2)
sum(class_UD_M, na.rm=TRUE)/2 #400

class_UDtemp <- ((sesbetaNTI_F<2 & sesbetaNTI_F>(-2)) + ((RCbray2_F<0.95)&(RCbray2_F>(-0.95))))
class_UD_F <- (class_UDtemp==2)
sum(class_UD_F, na.rm=TRUE)/2 #512

class_UDtemp <- ((sesbetaNTI_N<2 & sesbetaNTI_N>(-2)) + ((RCbray2_N<0.95)&(RCbray2_N>(-0.95))))
class_UD_N <- (class_UDtemp==2)
sum(class_UD_N, na.rm=TRUE)/2 #26

stageM<- c(sum(class_VS_M, na.rm=TRUE),sum(class_HS_M, na.rm=TRUE)/2, sum(class_DL_M, na.rm=TRUE)/2, 
          sum(class_HD_M, na.rm=TRUE)/2,sum(class_UD_M, na.rm=TRUE)/2)
stageF<- c(sum(class_VS_F, na.rm=TRUE),sum(class_HS_F, na.rm=TRUE)/2, sum(class_DL_F, na.rm=TRUE)/2, 
          sum(class_HD_F, na.rm=TRUE)/2,sum(class_UD_F, na.rm=TRUE)/2)
stageN<- c(sum(class_VS_N, na.rm=TRUE),sum(class_HS_N, na.rm=TRUE)/2, sum(class_DL_N, na.rm=TRUE)/2, 
          sum(class_HD_N, na.rm=TRUE)/2,sum(class_UD_N, na.rm=TRUE)/2)

stage <- c(rep("all" , 5) , rep("male" , 5) , rep("female" , 5) , rep("nymph" , 5) )
proc <- rep(c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),4)
value2 <- c(all,stageM,stageF,stageN)
data2 <- data.frame(stage,proc,value2)

#Figure in Supplement

ggplot(data2, aes(fill=proc, y=value2, x=stage)) + 
  geom_bar(position="fill", stat="identity")
proc <- rep(c("Variable selection","Homogenous selection","Dispersal limitation","Homogenizing dispersal","Undominated"),4)