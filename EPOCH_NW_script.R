########################################################################################
########################################################################################
###                 Network Analysis of the Chinese version of EPOCH                 ###
###                  Zeng, Guang; Peng, Kaiping; Hu, Chuan-Peng                      ### 
###               Email = hcp4715@gmail.com       twitter= @hcp4715                  ###
########################################################################################
########################################################################################

########################################################################################
############################# Acknowledgement ##########################################
###  This script is based on the R code from the following two papers ##################
###  Briganti,Eiko, LINKOWSKI, Braun,Kempenaers, (2018),                             ###
###        http://dx.doi.org/10.17632/b8d5n5h3gc.1                                   ###
###  Bellet, Jones, Neimeyer, & McNally. (2018).                                     ###
###        http://dx.doi.org/10.1177/2167702618777454                                ###
########################################################################################
########################################################################################

########################################################################################
### Preparing
###
rm(list = ls())     # remove all variables

curWD <- dirname(rstudioapi::getSourceEditorContext()$path) #Get the directory ofcurrent script


#required packages (install them first = install.packages(name of the package))
if (!require(qgraph)) {install.packages("qgraph",repos = "http://cran.us.r-project.org"); require(qgraph)}
library("qgraph")     # qgraph for network analysis

if (!require(bootnet)) {install.packages("bootnet",repos = "http://cran.us.r-project.org"); require(bootnet)}
library("bootnet")    # bootnet for network analysis

if (!require(tidyverse)) {install.packages("tidyverse",repos = "http://cran.us.r-project.org"); require(tidyverse)}
library("tidyverse")  # use tidyverse for data manipulation

if (!require(haven)) {install.packages("haven",repos = "http://cran.us.r-project.org"); require(haven)}
library("haven")      # read spss file

if (!require(mgm)) {install.packages("mgm",repos = "http://cran.us.r-project.org"); require(mgm)}
library("mgm")        # 


#library(qgraph)
#library(stats)
#library(readr)
#library(bootnet)
#library(igraph)
#library(dplyr)
#library(mgm)
#library(reshape2)
#library(data.table)

# Load the data in spss format 
data <- haven::read_spss("EPOCHreducedForR2.sav")  #  the short version data

# the colnames for test at T1
t1name <- c("E1.1", "E2.1", "E3.1", "E4.1",
            "P1.1", "P2.1", "P3.1", "P4.1", 
            "O1.1", "O2.1", "O3.1", "O4.1",
            "C1.1", "C2.1", "C3.1", "C4.1",
            "H1.1", "H2.1", "H3.1", "H4.1")

# the colnames for test at T2
t2name <- c("E1.2", "E2.2", "E3.2", "E4.2",
            "P1.2", "P2.2", "P3.2", "P4.2", 
            "O1.2", "O2.2", "O3.2", "O4.2",
            "C1.2", "C2.2", "C3.2", "C4.2",
            "H1.2", "H2.2", "H3.2", "H4.2")

# Name for each item
epochName <- c("E-forgetTime", "E-focus",   "E-flow",      "E-learn",
               "P-beginEnd",   "P-finish",  "P-plan",      "P-hardwork",
               "O-future",      "O-best",    "O-goodthing", "O-solution",
               "C-share",      "C-support", "C-caring",    "C-friends",
               "H-happy",      "H-fun",     "H-life",      "H-joy")

# select columns for the first testing
df1 <- data %>%              
      dplyr::select(t1name)

# select columns for the second testing
df2 <- data %>%              
      dplyr::select(t2name) %>%
      dplyr::filter(complete.cases(data))


# column names
#colnames(data) <- c(1:28)

#correlation matrix
cor_davis <- cor(df1, method="spearman")

#names
#names<- c("1FS", "2EC", "3PT_R", "4EC_R", "5FS", "6PD", "7FS_R", 
#   "8PT","9EC", "10PD", "11PT", "12FS_R", "13PD_R", "14EC_R", "15PT_R", 
#     "16FS", "17PD", "18EC_R", "19PD_R", "20EC", "21PT", "22EC", "23FS", 
#       "24PD", "25PT", "26FS", "27PD", "28PT")

#groups
#gr <- list(c(1, 5, 7, 12, 16, 23, 26), c(3, 8, 11, 15, 21, 25, 28),
#  c(2, 4, 9, 14, 18, 20, 22), c(6, 10, 13, 17, 19, 24, 27))

gr <- list(c(1, 2, 3, 4),
           c(5, 6, 7, 8),
           c(9, 10, 11, 12),
           c(13, 14, 15, 16),
           c(17, 18, 19, 20))

#gr <- list(c("E1.1", "E2.1", "E3.1", "E4.1"),
#      c("P1.1", "P2.1", "P3.1", "P4.1"),
#      c("O1.1", "O2.1", "O3.1", "O4.1"),
#      c("C1.1", "C2.1", "C3.1", "C4.1"),
#      c("H1.1", "H2.1", "H3.1", "H4.1"))

# estimate gaussian graphical model using spearman correlations
network1 <- estimateNetwork(df1, default="EBICglasso", corMethod = "cor", corArgs =
                    list(method = "spearman", use = "pairwise.complete.obs"))
# node predictability
type <- rep('g', 20)                 #g=gaussian, 20 = number of nodes in the network

fit1 <- mgm(df1,                     # estimate k-degree mixed graphical models
            type = type,
            level = rep(1, 20))

pred1<- predict(fit1, df1)
pred1$error
mean(pred1$error$R2) # shows the mean predictability of a node in the network = 0.4874

#visualize network with original data
graph1 <-plot(network1,
            labels = TRUE,
            nodeNames = epochName,
            pie = pred1$error$R2,
            layout = "spring",
            groups = gr,
            color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
            legend.cex = .35)

pdf("Figure1.pdf", width=10, height=7)
qgraph(graph1)
dev.off()

########################################################################################
####  Below is the original code from Briganti et al (2018),not changed yet          ###
########################################################################################


#graph1 weight matrix
graph1mat <- getWmat(graph1)
graph1mat #visualize weight matrix

#mean graph1 weight matrix
mean(graph1mat) #visualize the mean edge weight of the network = 0.025

#spinglass algorithm
#run the spinglass algorithm 100 times to get a mean of communities, max and min values
g = as.igraph(graph1, attributes=TRUE)
matrix_spinglass <- matrix(NA, nrow=1,ncol=100)
for (i in 1:100) {
  set.seed(i)
  spinglass <- spinglass.community(g)
  matrix_spinglass[1,i] <- max(spinglass$membership) 
}
mean(as.vector(matrix_spinglass)) # 4.2
max(as.vector(matrix_spinglass)) # 6
min(as.vector(matrix_spinglass)) # 4
median(as.vector(matrix_spinglass)) # 4
sgc <- spinglass.community(g) #4 communities identified!
sgc$membership #shows the membership of an item to a community

#walktrap algorithm
glasso.ebic <-EBICglasso(S=cor_davis, n = nrow(data)) #build a graph object for the algorithm
graph.glasso <-as.igraph(qgraph(glasso.ebic, layout = "spring", vsize = 3))
wc<- walktrap.community(graph.glasso) #shows the membership of an item to a community
n.dim <- max(wc$membership)  #visualize communities of items: 6, 10 and 17 form a 5th community!

#eigenvalue plot
plot(eigen(cor_davis)$values, type="b")
abline(h=1,col="red", lty = 3)

#Fig.2 Davis Network centrality plot, in the Supplementary Materials
graph1_cp <- centralityPlot(graph1) #using the graph from the data
pdf("Figure2bis.pdf", width=10, height=7) 
qgraph(graph1_cp) 
dev.off()

#centrality criteria 
graph1.c <- centrality(graph1)
graph1.c$InDegree
graph1.c$Closeness
graph1.c$Betwennness
cor(graph1.c$InDegree, graph1.c$Betweenness, method = "spearman") 
cor(graph1.c$InDegree, graph1.c$Closeness, method = "spearman") 
cor(graph1.c$Closeness, graph1.c$Betweenness, method = "spearman") 

cent <- as.data.frame(scale(centrality(graph1)$InDegree))
cent <- mutate(cent, id = rownames(cent))
colnames(cent) <- c("1", "IRI_Item")
cent_long <- melt(cent, id="IRI_Item")

pdf("Figure2.pdf", width=6, height=4, useDingbats=FALSE)
strengthplot <- ggplot(data=cent_long, aes(x=IRI_Item, y=value, group=1)) +
  geom_line() +
  geom_point(shape = 21, fill = "white", size = 1.5, stroke = 1) +
  xlab(" ") + ylab("Centrality") +
  scale_y_continuous(limits = c(-3, 3)) + 
  scale_x_discrete(limits=c(1:28)) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
strengthplot
# coord_flip()
dev.off()

# stabilility
boot1 <- bootnet(network1, ncores=7, nboots=2000)                
boot2 <- bootnet(network1, ncores=7, nboots=2000, type="case")
save(boot1, file = "boot1davis.Rdata")
save(boot2, file = "boot2davis.Rdata")

#load the boots: modify file path
load(file = "boot1.Rdata")
load(file = "boot2.Rdata")

# Fig3 - Edge weight bootstrap
fig3 <- plot(boot1, labels = FALSE, order = "sample")
pdf("Figure3.pdf", width=10, height=7) 
plot(boot1, labels = FALSE, order = "sample") 
dev.off()

#Fig3 - Edge weight difference: is edge X significantly larger than edge Y? Black=Y Gray=N 
boot3 <- plot(boot1, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
plot(boot3)
pdf("Figure4.pdf", width=10, height=7)
plot(boot3)
dev.off()

#Fig5 - Centrality Bootstrap
plot(boot2)
cs1 <- corStability(boot2)  
pdf("Figure5.pdf", width=10, height=7)
plot(boot2)
dev.off()

#Fig6 - Centrality difference: is node X significantly more central than node Y? Black=Yes, Gray=N
boot4 <- plot(boot1, "strength", order="sample", labels=TRUE) 
pdf("Figure6.pdf", width=10, height=7)
plot(boot4) 
dev.off()

