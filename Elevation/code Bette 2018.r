if (!require(igraph)) {install.packages("igraph",repos = "http://cran.us.r-project.org"); require(effects)}
library("igraph")
if (!require(qgraph)) {install.packages("qgraph",repos = "http://cran.us.r-project.org"); require(qgraph)}
library("qgraph")
if (!require(foreign)) {install.packages("foreign",repos = "http://cran.us.r-project.org"); require(foreign)}
library("foreign")
if (!require(bootnet)) {install.packages("bootnet",repos = "http://cran.us.r-project.org"); require(bootnet)}
library("bootnet")
if (!require(dplyr)) {install.packages("dplyr",repos = "http://cran.us.r-project.org"); require(dplyr)}
library("dplyr")
if (!require(gridExtra)) {install.packages("gridExtra",repos = "http://cran.us.r-project.org"); require(gridExtra)}
library("gridExtra")
if (!require(gridBase)) {install.packages("gridBase",repos = "http://cran.us.r-project.org"); require(gridBase)}
library("gridBase")
source("http://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")
if (!require(devtools)) {install.packages("devtools",repos = "http://cran.us.r-project.org"); require(devtools)}
library("devtools")
if (!require(networktools)) {install.packages("networktools",repos = "http://cran.us.r-project.org"); require(networktools)}
library("networktools")
if (!require(CTT)) {install.packages("CTT",repos = "http://cran.us.r-project.org"); require(CTT)}
library("CTT")

sessionInfo()

fulldata <- read.spss("SmilesData.sav", to.data.frame=T)
Demograph <- fulldata[!is.na(fulldata$PTGI1),]

##Demographics Descriptives
dim(Demograph)
table(Demograph$gender)
table(Demograph$race)
mean(Demograph$age, na.rm = TRUE)
sd(Demograph$age, na.rm = TRUE)
table(Demograph$cause)
table(Demograph$relat)
table(Demograph$peduc)

##Univariate Descriptives

CGitems <- Demograph[,c(paste("ICG", c(4:5, 7, 9:11, 14:15, 17, 22, 25, 26, 28), sep = ""))]
PTGitems <-Demograph[,c(paste("PTGI", c(1:2, 5, 7:8, 10:11, 18:20), sep = ""))]

colnames(CGitems)<-c("noacc", "yearn", "anger", "shock", "notrst", "nocare", "discon", "mnglss", "numb", "futmng", "identi", "wldvw", "nctrl")
colnames(PTGitems)<-c("chgpri", "applif", "spirit", "nwpth", "close", "handle", "better",  "faith", "strgr", "people")
descriptCG <- cbind(apply(CGitems, 2, mean), apply(CGitems, 2, sd))
descriptPTG<-cbind(apply(PTGitems,2, mean), apply(PTGitems, 2, sd))
colnames(descriptCG)<-c("Mean", "SD")
colnames(descriptPTG) <-c("Mean", "SD")

##Reliability Analyses
itemAnalysis(CGitems) ## reliability of .91
itemAnalysis(PTGitems) ## reliability of .92

##Network analysis

###Making the data matrices
netdata3 <- fulldata[,c(paste("ICG", c(4:5, 7, 9:11, 14:15, 17, 22, 25, 26, 28), sep=""), paste("PTGI", c(1:2, 5, 7:8, 10:11, 18:20), sep=""))]
netdata3 <- na.omit(netdata3)
dim(netdata3)
netdata3<-data.matrix(netdata3)

dim(netdata3)
netdata3<-data.matrix(netdata3)
nodenamesCG <- c("noacc", "yearn", "anger", "shock", "notrst", "ncare", "discon", "mnglss", "numb", "futmng", "identi", "wldvw", "nctrl")
nodenamesPTGI <- c("chgpri", "applif", "spirit", "newpth", "close", "handle", "better",  "faith", "strngr", "people")
nodenamesALL <- c(nodenamesCG, nodenamesPTGI)
nodenamesCG_plot <- c("noacc", "yearn", "anger", "shock", "notrst", "ncare", "discon", "mnglss", "numb", "futmng", "IDENTI", "wldvw", "NCTRL")
nodenamesPTGI_plot <- c("chgpri", "applif", "spirit", "NWPTH", "close", "handle", "better",  "faith", "STRGR", "people")
nodenamesALL_plot <- c("noacc", "yearn", "anger", "shock", "notrst", "NCARE", "discon", "mnglss", "numb", "futmng", "identi", "WLDVW", "nctrl", "chgpri", "applif", "spirit", "nwpth", "close", "handle", "better",  "faith", "strgr", "people")

cg_ptg <- netdata3[,c(paste("ICG", c(4:5, 7, 9:11, 14:15, 17, 22, 25, 26, 28), sep=""), paste("PTGI", c(1:2, 5, 7:8, 10:11, 18:20), sep=""))]
cg_only <- netdata3[,c(paste("ICG", c(4:5, 7, 9:11, 14:15, 17, 22, 25, 26, 28), sep=""))]
ptg_only <- netdata3[,c(paste("PTGI", c(1:2, 5, 7:8, 10:11, 18:20), sep=""))]
cg_ptg_d <- data.frame(cg_ptg)
colnames(netdata3[,1:23]) <- nodenamesALL

colnames(cg_ptg_d) <- nodenamesALL
colnames(cg_only) <- nodenamesCG
colnames(cg_ptg) <- nodenamesALL
colnames(ptg_only) <- nodenamesPTGI

###CG Network

cor_cg_only <- invisible(cor_auto(cg_only))
glasso_cg_only <- EBICglasso(cor_cg_only, n=dim(cg_only)[1], gamma=0.5)
CGnet <- qgraph(glasso_cg_only, layout="spring", labels=nodenamesCG_plot, vsize=9,
                label.cex=c(rep(.7,10),.6,.7,.6),
                label.font=c(rep(1,10),2,1,2),
                label.color=c(rep(1, 10),"blue", 1,"blue"),
                label.scale=F, DoNotPlot=F)

pdf("CGn.pdf", width = 3.8, height = 4)
CGbw<-makeBW(CGnet)
dev.off()


###CG Expected Influence
EI_cg <- expectedInf(glasso_cg_only)
pdf("CGei.pdf", width = 5)
plot(EI_cg$step1, order="value", zscore=F, yaxt = "n")
dev.off()
###Edge Weight and EI stability
set.seed(123)
net <- estimateNetwork(cg_only, default="EBICglasso")

net_boot <- bootnet_flex(net, statistics=c("edge", "expectedInf"), nBoots=1000, type="case", caseN = 50)

CorStabCG<-corStability(net_boot, statistics=c("edge", "expectedInf"))

###PTG Network

cor_ptg_only <- invisible(cor_auto(ptg_only))
glasso_ptg_only <- EBICglasso(cor_ptg_only, n=dim(ptg_only)[1], gamma=0.5)
PTGInet<- qgraph(glasso_ptg_only, layout="spring", labels=nodenamesPTGI_plot, vsize=9,
                 label.cex=c(rep(.7,3),.6,rep(.7, 4),.6,.7),
                 label.font=c(rep(1,3),2,rep(1,4),2),
                 label.color=c(rep(1,3),"blue", rep(1,4),"blue"),
                 label.scale=F, DoNotPlot=F)

PTGbw<-makeBW(PTGInet)

pdf("PTGn.pdf", width = 3.8, height = 4)
PTGbw<-makeBW(PTGInet)
dev.off()


###PTG Expected Influence
EI_ptg <- expectedInf(glasso_ptg_only)
PTGIei <- plot(EI_ptg$step1, order="value", zscore=F)
PTGIei

pdf("PTGei.pdf", width = 5)
plot(EI_ptg$step1, order="value", zscore=F)
dev.off()

###Edge Weight and EI stability
net2 <- estimateNetwork(ptg_only, default="EBICglasso")

net_boot2 <- bootnet_flex(net2, statistics=c("edge", "expectedInf"), nBoots=1000, type="case", caseN = 50)

CorStabPTG <-corStability(net_boot2, statistics=c("edge", "expectedInf"))

###Combined Network
cor_netCombo <- cor_auto(cg_ptg)
glassoCombo <- EBICglasso(cor_netCombo, n=dim(cg_ptg)[1], gamma=0.5)
COMBOnet <- qgraph(glassoCombo, layout="spring", labels=nodenamesALL_plot, groups=list(CG=1:13, PTG=14:23),legend=F, 
                   colors= c("white", "lightblue"),
                   label.cex=c(rep(.6, 5),.5,rep(.6,5),.5, rep(.6,11)),
                   label.font=c(rep(1, 5),2, rep(1,5),2, rep(1,11)),
                   label.color=c(rep(1, 5),"darkred",rep(1,5),"blue", rep(1,11)),
                   label.scale=F, DoNotPlot=F)

write.csv(cor_netCombo, file = "S3.csv")
ComboBW <- makeBW(COMBOnet)

pdf("COMBOn.pdf", width = 5, height = 4)
COMBObw<-makeBW(COMBOnet)
dev.off()

bridgesCombo <- bridge(glassoCombo, communities=c(rep("CG", 13), rep("Growth", 10)))
BEI<-plot(bridgesCombo, order="value", zscore=F, include = "Bridge Expected Influence (1-step)")

pdf("BEI.pdf", width = 5)
plot(bridgesCombo, order="value", zscore=F, include = "Bridge Expected Influence (1-step)")
dev.off()

###spinglass community detection
spinglass.community(graph_from_adjacency_matrix(glassoCombo, mode = "lower", weighted = TRUE), spins = 2)

###Edge weight stability
net3 <- estimateNetwork(cg_ptg, default="EBICglasso")

net_boot3 <- bootnet_flex(net3, statistics=c("edge", "expectedInf"), nBoots=1000, type="case", caseN = 50)

CorStabPTG <-corStability(net_boot2, statistics=c("edge", "expectedInf"))

