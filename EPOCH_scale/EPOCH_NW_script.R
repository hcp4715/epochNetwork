########################################################################################
###                                                                                  ###
###                 Network Analysis of the Chinese version of EPOCH                 ###
###                  Zeng, Guang; Peng, Kaiping; Hu, Chuan-Peng                      ### 
###               Email = hcp4715@gmail.com       twitter= @hcp4715                  ###
###                                                                                  ###
########################################################################################

########################################################################################
############################# Acknowledgement ##########################################
###  This script is based on the R code from the following papers:                   ###
###  Briganti,Eiko, LINKOWSKI, Braun,Kempenaers, (2018),                             ###
###        http://dx.doi.org/10.17632/b8d5n5h3gc.1                                   ###
###  Bellet, Jones, Neimeyer, & McNally. (2018).                                     ###
###        http://dx.doi.org/10.1177/2167702618777454                                ###
###  Fried et al. 2018, Clinical Psychological Science,                              ###
###        https://doi.org/10.1177/2167702617745092                                  ###
########################################################################################
########################################################################################

# ---------- Table of Contents ----------------------------------------------------------
# ---------- 1. Load Libraries ----------------------------------------------------------
# ---------- 2. Load & manipulate data --------------------------------------------------
# ---------- 3. Estimate networks for 1st and 2nd measurement ---------------------------
# ---------- 4. Stability and accuracy for estimated networks ---------------------------
# ---------- 5. Network comparison tests for 1st and 2nd measurement ---- ---------------

# ---------------------------------------------------------------------------------------
# ---------- 1. Load libraries ----------------------------------------------------------
# ---------------------------------------------------------------------------------------

rm(list = ls())     # remove all variables
curWD <- dirname(rstudioapi::getSourceEditorContext()$path) #Get the directory ofcurrent script
setwd(curWD)

#required packages (install them first = install.packages(name of the package))
if (!require(qgraph)) {install.packages("qgraph",repos = "http://cran.us.r-project.org"); require(qgraph)}
library("qgraph")     # qgraph for network analysis

if (!require(bootnet)) {install.packages("bootnet",repos = "http://cran.us.r-project.org"); require(bootnet)}
library("bootnet")    # bootnet for network analysis

if (!require(haven)) {install.packages("haven",repos = "http://cran.us.r-project.org"); require(haven)}
library("haven")      # read spss file

if (!require(mgm)) {install.packages("mgm",repos = "http://cran.us.r-project.org"); require(mgm)}
library("mgm")        # 

if (!require(tidyverse)) {install.packages("tidyverse",repos = "http://cran.us.r-project.org"); require(tidyverse)}
library("tidyverse")  # use tidyverse for data manipulation

#library(stats)
#library(data.table)


# ---------------------------------------------------------------------------------------
# ---------- 2. Load & manipulate data --------------------------------------------------
# ---------------------------------------------------------------------------------------

data <- haven::read_spss("EPOCH_Clean_no_id.sav")  #  the short version data

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
df1_all <- data %>%              
      dplyr::select(t1name)

# select columns for the second testing
df1_r <- data %>%                              # 1st measurement with re-test
      dplyr::select(c(t1name,t2name)) %>%
      dplyr::filter(complete.cases(data)) %>%
      dplyr::select(t1name)
df2 <- data %>%                                # re-test
      dplyr::select(c(t1name,t2name)) %>%
      dplyr::filter(complete.cases(data)) %>%
      dplyr::select(t1name)

# ---------------------------------------------------------------------------------------
# ---------- 3. Estimate networks individually ------------------------------------------
# ---------------------------------------------------------------------------------------

#correlation matrix
df1_all_corr <- qgraph::cor_auto(df1_all)
df1_r_corr   <- qgraph::cor_auto(df1_r)
df2_corr     <- qgraph::cor_auto(df2)

# corrr::rplot(df2_corr,shape = 15,colors = c('red','green')) # visualize the corelation

# group the items for later legend
gr <- list(c(1, 2, 3, 4),
           c(5, 6, 7, 8),
           c(9, 10, 11, 12),
           c(13, 14, 15, 16),
           c(17, 18, 19, 20))

# estimate gaussian graphical model using spearman correlations
# network1_all <- estimateNetwork(df1, default="EBICglasso", corMethod = "cor", corArgs =
#                    list(method = "spearman", use = "pairwise.complete.obs"))
network1_all <- bootnet::estimateNetwork(df1_all,default = "EBICglasso")
network1_r   <- bootnet::estimateNetwork(df1_r,  default = "EBICglasso")
network2     <- bootnet::estimateNetwork(df2,    default = "EBICglasso")

# using average layout late
Layout <- qgraph::averageLayout(network1_all,network1_r,network2)

### Follow previous studies, we analyzed the node predictability
# Re-estimate individual networks via mgm 
type <- rep('g', 20)                 #g=gaussian, 20 = number of nodes in the network

fit1_all  <- mgm::mgm(na.omit(df1_all), type=type, lev=rep(1,20)) 
pred1_all <- predict(fit1_all, na.omit(df1_all))

fit1_r  <- mgm::mgm(na.omit(df1_r), type=type, lev=rep(1,20)) 
pred1_r <- predict(fit1_r, na.omit(df1_r))

fit2    <- mgm::mgm(na.omit(df2), type=type, lev=rep(1,20)) 
pred2   <- predict(fit2, na.omit(df2))

# The average node predictability(the proportion of explained variance("R2"))
mean(pred1_all$error$R2) # .48725
mean(pred1_r$error$R2)   # .4968
mean(pred2$error$R2)     # .49675

### Plot networks
# network of all data
pdf("Figure1.pdf", width=10, height=7)
network1_all_G <- plot(network1_all,
                       title="All data", 
                       maximum=.47,
                       labels = TRUE,
                       nodeNames = epochName,
                       pie = pred1_all$error$R2,
                       layout = Layout,
                       groups = gr,
                       color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                       legend.cex = .35,
                       border.width=2, border.color='#555555')       # max=0.47
#qgraph(graph1)
dev.off()

pdf("Fig_S1_retest_T1_T2.pdf", width= 16, height= 9)
par(mfrow=c(1,2))
network1_r_G <- plot(network1_r,
                     title="Re-test data T1", 
                     maximum=.47,
                     labels = TRUE,
                     nodeNames = epochName,
                     pie = pred1_r$error$R2,
                     layout = Layout,
                     groups = gr,
                     color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                     legend.cex = .35,
                     border.width=2, border.color='#555555')
network2_G <- plot(network2,
                   title="Re-test data T2", 
                   maximum=.47,
                   labels = TRUE,
                   nodeNames = epochName,
                   pie = pred2$error$R2,
                   layout = Layout,
                   groups = gr,
                   color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                   legend.cex = .35,
                   border.width=2, border.color='#555555')
dev.off()

### Save output
# save(data1cor,data2cor,data3cor,data4cor, file="data_cormatrices.RData")
save(network1_all, network1_r, network2, file="networks_3_sep.RData")
# save(fit1,pred1,fit2,pred2,fit3,pred3,fit4,pred4, file="mgm.Rdata")

# ---------------------------------------------------------------------------------------
# ---------- 4. Community and centrality of the estimated networks ----------------------
# ---------------------------------------------------------------------------------------

# (1) used the Eigenvalue
# eigenvalue plot
plot(eigen(df1_all_corr)$values, type="b")
abline(h=1,col="red", lty = 3)

# (2) Using spinglass algorithm to get the community informtion
# run the spinglass algorithm 100 times to get a mean of communities, max and min values
g_all <- igraph::as.igraph(network1_all_G, attributes=TRUE)

max_comm  <- matrix(NA, nrow=1,ncol=100)   # create a empty matrix 
comm_memb <- matrix(NA, nrow=100,ncol=20)  # create a empty matrix for storing the membership
for (i in 1:100) {
      set.seed(i)
      spinglass <- igraph::spinglass.community(g_all)
      max_comm[1,i] <- max(spinglass$membership)
      comm_memb[i,] <- spinglass$membership
}

summary(as.vector(matrix_spinglass))      # min = 5; median = 5, mean = 5.09; max = 7
hist(as.vector(matrix_spinglass))         # view the distribution of the # of communities

set.seed(100)
sgc <- igraph::spinglass.community(g_all) # 5 communities identified
sgc$membership                            # shows the membership of an item to a community
# 4 5 4 4 5 5 5 5 3 1 1 1 2 2 2 2 3 3 3 3

# (3) using EGA package to detect the community (if not installed, use following two lines of code)
# library("devtools")
# devtools::install_github('hfgolino/EGA')
library("EGA")
ega <- EGA::EGA(df1_all, plot.EGA = TRUE)
ega$wc  #3 4 3 3 4 4 4 4 1 2 2 2 2 2 2 2 1 1 1 1

# walktrap algorithm for community (random walk approach)
glasso.ebic   <- qgraph::EBICglasso(S=df1_all_corr, n = nrow(df1_all)) #build a graph object for the algorithm
graph.glasso  <- igraph::as.igraph(qgraph(glasso.ebic, layout = "spring", vsize = 3))
walktrap_comm <- igraph::walktrap.community(graph.glasso) #shows the membership of an item to a community
walktrap_comm$membership

## plot two approach of communities and save as "Fig_S1_comm.pdf"
group.spinglass<- list(c(1,3,4), c(2,5:8), c(9,17:20), c(10:12), c(13:16))
group.ega      <- list(c(1,3,4), c(2,5:8), c(9,17:20), c(10:16))
pdf("Fig_S2_comm.pdf", width= 16, height= 9)
par(mfrow=c(2,2))
par(mfrow = c(1, 2))
graph_ega <- qgraph(df1_all_corr, graph="glasso", layout="spring", sampleSize = nrow(data),groups=group.spinglass, 
                    vsize=7, cut=0, maximum=.45, border.width=1.5,
                    color=c("red", "green", "blue", "orange","white"), title="igraph spinglass")

graph_sg <- qgraph(df1_all_corr, graph="glasso", layout="spring", sampleSize = nrow(data),groups=group.ega, 
                 vsize=7, cut=0, maximum=.45, border.width=1.5,
                 color=c("red", "green", "blue", "orange"), title="ega walktrap")
dev.off()

### Centrality
#Fig.2 Davis Network centrality plot, in the Supplementary Materials
All_cp <- centralityPlot(network1_all_G) #using the graph from the data
pdf("Figure_S3_centrality_all.pdf", width=10, height=7) 
All_cp 
dev.off()

# centrality criteria 
graph1.c <- centrality(network1_all_G)
graph1.c$InDegree
graph1.c$Closeness
graph1.c$Betweenness
cor(graph1.c$InDegree, graph1.c$Betweenness) # 0.646
cor(graph1.c$InDegree, graph1.c$Closeness)   # 0.722
cor(graph1.c$Closeness, graph1.c$Betweenness)# 0.5268 

cent <- as.data.frame(scale(centrality(network1_all_G)$InDegree))
cent <- dplyr::mutate(cent, id = rownames(cent))
colnames(cent) <- c("1", "EPOCH_Item")
cent_long <- reshape2::melt(cent, id="EPOCH_Item")

# plot Figure 2. 
strengthplot <- ggplot2::ggplot(data=cent_long, aes(x=EPOCH_Item, y=value, group=1)) +
      geom_line() +
      geom_point(shape = 21, fill = "white", size = 1.5, stroke = 1) +
      xlab(" ") + ylab("Centrality") +
      scale_y_continuous(limits = c(-3, 3)) + 
      scale_x_discrete(limits=c(1:20)) +
      theme_bw() +
      theme(panel.grid.minor=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
pdf("Figure2.pdf", width=6, height=4, useDingbats=FALSE)
strengthplot
# coord_flip()
dev.off()

# ---------------------------------------------------------------------------------------
# ---------- 4. Stability and accuracy for individually estimated networks --------------
# ---------------------------------------------------------------------------------------

### Estimate and save stability and accuracy
boot1a <- bootnet::bootnet(network1_all, nBoots = 1000, nCores = 8)
boot1b <- bootnet::bootnet(network1_all, nBoots = 1000, type = "case",  nCores = 8)

boot2a <- bootnet::bootnet(network2, nBoots = 1000, nCores = 8)
boot2b <- bootnet::bootnet(network2, nBoots = 1000, type = "case",  nCores = 8)

### Plot edge weight CI
pdf("Fig2_edge_weight_ci.pdf")
plot(boot1a, labels = FALSE, order = "sample") 
plot(boot2a, labels = FALSE, order = "sample") 
dev.off()

### Plot centrality stability
pdf("FigS3_centrality_stablity.pdf") 
plot(boot1b)
plot(boot2b)
dev.off()

### Centrality stability coefficient
cs1 <- bootnet::corStability(boot1b)  
cs2 <- bootnet::corStability(boot2b) 
cs <- matrix(NA,2,3)
cs[1,1:3] <- round(cs1, digits=2)
cs[2,1:3] <- round(cs2, digits=2)
colnames(cs) <- c("Betweenness", "Closeness", "Strength")
rownames(cs) <- c("All_data_T1", "Rep_Data_T2")
sink("TableS1.txt")
cs
sink()

### Edge weights diff test
pdf("Fig4_Edge_weight_diff.pdf")
plot(boot1a, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
plot(boot2a, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
dev.off()

### Centrality diff test
pdf("Fig5_centrality_diff.pdf")
plot(boot1a, "strength", order="sample", labels=FALSE) 
plot(boot2a, "strength", order="sample", labels=FALSE) 
dev.off()

### Conclusion: 
# all networks have CS coefficients for node strength > 0.5; looks good!

### test- retest network comparison
dft1_c1  <- cor(df1_r)
dft1_c1b <- qgraph::cor_auto(df1_r)
cor(dft1_c1[lower.tri(dft1_c1)], dft1_c1b[lower.tri(dft1_c1b)], 	method="spearman") #0.9936

dft2_c1  <- cor(df2)
dft2_c1b <- qgraph::cor_auto(df2)
cor(dft2_c1[lower.tri(dft2_c1)], dft2_c1b[lower.tri(dft2_c1b)], 	method="spearman") #0.9936

# since they are very similar, we can probably use the NetworkComparisonTest::NCT() 
# that is only validated for Pearson, not poylchoric correlations

### NCT topology
set.seed(100)
compare_12 <- NetworkComparisonTest::NCT(df1_r,df2, it=5000, binary.data=FALSE, test.edges=TRUE, edges='all', progressbar=TRUE)

compare_12$nwinv.pval     # output p-value: 1

### quantification of differences: count significantly different edges (total number )
sum(compare_12$einv.pvals$"p-value" < 0.05) # 0 of all

### NCT global strength
compare_12$glstrinv.pval  # 1

### the network are similar
