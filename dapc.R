## Adapting workflow from Jessi Rick's DAPC script associated with Junker et al 2020 (sardines); pulled from Jessi's github
## git clone https://github.com/jessicarick/lake-tanganyika-sardines.git

source("packages_functions.R")

ccvcf <- read.vcfR("CC_alignedjune2021_miss0.4_mac3_Q20_DP3_thin90_maf001.recode.vcf")
##ccvcf <- read.vcfR("CCCS_cccsref0.9test_miss0.4_mac3_Q30_DP3_ind95_maf001_1snpcontig_sexinds.recode.vcf")
##ccvcf <- read.vcfR("CCCS_cccsref0.9test_miss0.2_mac3_Q20_knownsexinds_thin90_maf001.recode.vcf")


##convert the genotypes (encoded as ./., 0/0, 0/1, 1/1) into genotype values (0,1,2) in a "genlight" object

ccgen <- vcfR2genlight(ccvcf)
col.names <- unlist(strsplit(indNames(ccgen),"/project/rrg-emandevi/agstreams_2020/radsex_indivs/align2dace/bwa_assem2/"))
##col.names <- unlist(strsplit(indNames(ccgen),"/project/def-emandevi/agstreams_2020/bwa_assem_cccsref0.9test/"))
col.names.clean <- unlist(strsplit(col.names,".sorted.bam"))

## Pull in lists of males and females, merge in the order in the vcf
females.paths <- read.table("females_bampaths.txt")
females1 <- unlist(strsplit(as.character(females.paths[,1]), "/project/rrg-emandevi/agstreams_2020/radsex_indivs/align2dace/bwa_assem2/"))
females2 <- unlist(strsplit(females1,".sorted.bam"))
females <- data.frame(females2, rep("F", length(females2)))
colnames(females) <- c("ind", "sex")

males.paths <- read.table("males_bampaths.txt")
males1 <- unlist(strsplit(as.character(males.paths[,1]), "/project/rrg-emandevi/agstreams_2020/radsex_indivs/align2dace/bwa_assem2/"))
males2 <- unlist(strsplit(males1,".sorted.bam"))
males <- data.frame(males2, rep("M", length(males2)))
colnames(males) <- c("ind", "sex")

sexinfo <- rbind(females, males)

## Merge sexinfo with the names in VCF order

fishinfo <- merge(data.frame(col.names.clean), sexinfo, by.x="col.names.clean", by.y="ind", all.x=T, all.y=T)
colnames(fishinfo) <- c("ind", "sex")

## Do initial PCA
## I *think* that genotype matrix is supposed to be nind rows x nloci columns

ccpca <- do.pca(t(as.matrix(ccgen)))
pcSummarycc <- summary(ccpca)
scree <- plot(ccpca, type="lines")

## Do DAPC with sex as grouping variable
pop(ccgen) <- fishinfo$sex
dapc1 <- dapc(ccgen, n.pca=20, n.da=3)

scatter(dapc1, scree.da=TRUE,
        bg="white", pch=20, cell=0,
        cstar=0, solid=.4, cex=2,clab=0,
        leg=FALSE, col=viridis(5)[1:2])

loadingplot(dapc1$var.contr[,1],
            cex.lab=0.5,srt=90,byfac=FALSE)

## create null distribution by randomizing sexes ##
threshold.rand <- randomize.dapc(ccgen,pop=as.factor(fishinfo$sex),niter=100, return.all=TRUE, npca=20)
threshold <- quantile(threshold.rand,c(0.99),na.rm=TRUE)

## plot random and empirical loadings ##
loadings.cc <- as.data.frame(dapc1$var.contr)
loadings.cc$type <- "empirical"
rand.loadings.cc <- as.data.frame(threshold.rand)
rand.loadings.cc$type <- "random"
colnames(rand.loadings.cc)[1] <- "LD1"

all.loadings.cc <- rbind(loadings.cc,rand.loadings.cc)

p <- ggplot(all.loadings.cc) +
  geom_density(aes(x=LD1,group=type,fill=type),alpha=0.5) +
  #coord_cartesian(ylim=c(0,10000),xlim=c(0,0.0050)) +
  theme_minimal()+
  scale_fill_manual(values=c("gray20","gray90"),na.value = "grey90")+
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=14),
        axis.title = element_text(size=16)) +
  geom_vline(xintercept=threshold,linetype="dashed")
p


## plot loading plot with this threshold ##
loadingplot(dapc1$var.contr[,1],
            threshold=0.025,
            cex.lab=2,srt=90,
            byfac=FALSE,main="",
            xlab="SNP Location",ylab="SNP Loading",
            cex.axis=1.5)
abline(h=threshold,lty=2,col="gray20",lwd=2)

## pull out significant loci and plot by scaffold ##
sex.loci <- data.frame(SNP = locNames(ccgen)[dapc1$var.contr[,1] > threshold],
                       loading = dapc1$var.contr[dapc1$var.contr[,1] > threshold,])

write.csv(sex.loci, "CC_dacealign_miss0.4_sexloci_DAPC.csv", row.names=F, quote=F)
##write.csv(sex.loci, "CC_artificialref_miss0.2_sexloci_DAPC.csv", row.names=F, quote=F)
##write.csv(sex.loci, "CC_artificialref_miss0.4_sexloci_DAPC.csv", row.names=F, quote=F)

## make figure of dapc results
##layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))

##
pdf("dapc_cc_dacealign_miss0.4.pdf", width=8, height=4)
##pdf("dapc_cc_artificialref_miss0.2.pdf", width=8, height=4)
##pdf("dapc_cc_artificialref_miss0.4.pdf", width=8, height=4)
par(mfrow=c(1,2))
scatter(dapc1,
        scree.da=TRUE,
        bg="white",
        pch=20, cell=0, cstar=0,
        solid=0.4, cex=2,clab=0, leg=FALSE,
        col=viridis(5),
        cex.lab=2, cex.axis=1.5)

loadingplot(dapc1$var.contr[,1],
            threshold=0.025,
            cex.lab=2,srt=90,byfac=FALSE,main="",xlab="SNP Location",ylab="SNP Loading")
abline(h=threshold,lty=2,col="gray20",lwd=2)

dev.off()


## This part irrelevant for all except dace alignment, but that only IDs 5 scaffolds
## barplot(sort(table(sex.loci$scaffold),decreasing=TRUE),
##        las=2,
##        xlab="",
##        ylab="number of significant sex loci")


## CHeck heterozygosity for sex loci

ccgen_sex <- ccgen[,which(locNames(ccgen) %in% sex.loci$SNP)]
pop(ccgen_sex) <- fishinfo$sex

ccdiv <- gl.Ho(ccgen)
ccdiv_sex <- gl.Ho(ccgen_sex)

plot(ccdiv, xlab='Locus number',ylab='Observed heterozygosity')

#######################
## plotting sex loci ##
#######################
ccSexA <- ccgen[ccgen$pop == "F"]
ccSexB <- ccgen[ccgen$pop == "M"]

ccgenSexA_sex <- ccSexA[,which(locNames(ccSexA) %in% sex.loci$SNP)]
ccgenSexB_sex <- ccSexB[,which(locNames(ccSexB) %in% sex.loci$SNP)]

ccdiv_sexA <- gl.Ho(ccSexA)
ccdiv_sexB <- gl.Ho(ccSexB)
ccdiv_sexA_sexloci <- gl.Ho(ccgenSexA_sex)
ccdiv_sexB_sexloci <- gl.Ho(ccgenSexB_sex)

diff <- ccdiv_sexA_sexloci - ccdiv_sexB_sexloci

##pdf("sex_loci_Ho_artificialref0.4.pdf", width=8, height=4)
pdf("sex_loci_Ho_dacealign0.4.pdf", width=8, height=4)
##pdf("sex_loci_Ho_artificialref0.2.pdf", width=8, height=4)
par(mfrow=c(1,2))
plot(ccdiv_sexA,
     xlab='Locus',ylab='Observed heterozygosity',
     main="All loci, colored by group",xaxt='n',
     pch=19)
points(ccdiv_sexB,
       xlab='Locus',ylab='Observed heterozygosity',
       col="turquoise",
       pch=19)

plot(ccdiv_sexA_sexloci,
     xlab='Locus',ylab='Observed heterozygosity',
     ylim=c(0,1),
     main="Sex Loci, colored by group",
     xaxt='n', pch=19)
axis(1, at=1:length(ccdiv_sexA_sexloci),
     labels=row.names(as.data.frame(ccdiv_sexA_sexloci)),
     las=2)
points(ccdiv_sexB_sexloci,
       xlab='Locus number',ylab='Observed heterozygosity',
       col="turquoise",
       pch=19)

dev.off()

plot(diff,dapc1$var.contr[which(locNames(ccSexA) %in% sex.loci$SNP),1], xlab = "Ho diff F-M", ylab="DAPC loading")

plot(ccdiv_sexA_sexloci, ccdiv_sexB_sexloci)

hetsexloci <- cbind(ccdiv_sexA_sexloci, ccdiv_sexB_sexloci, diff, dapc1$var.contr[which(locNames(ccSexA) %in% sex.loci$SNP),1])
colnames(hetsexloci) <- c("hetF", "hetM", "diffFM", "loading")

##write.csv(hetsexloci, "het_sexloci_artificialref0.4.csv", quote=F)
##write.csv(hetsexloci, "het_sexloci_artificialref0.2.csv", quote=F)
write.csv(hetsexloci, "het_sexloci_dacealign0.4.csv", quote=F)

write.csv(data.frame(as.matrix(ccgen_sex), ccgen_sex$pop), "genotypes_sexloci_dacealign0.4.csv")
##write.csv(data.frame(as.matrix(ccgen_sex), ccgen_sex$pop), "genotypes_sexloci_artificialref0.4.csv")
##write.csv(data.frame(as.matrix(ccgen_sex), ccgen_sex$pop), "genotypes_sexloci_artificialref0.2.csv")
