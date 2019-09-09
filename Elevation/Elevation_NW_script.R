########################################################################################
###                                                                                  ###
###                 Network Analysis of the Chinese version of Elevation             ###
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

# ---------------------------------------------------------------------------------------
# ---------- 2. Load & manipulate data --------------------------------------------------
# ---------------------------------------------------------------------------------------

#data <- haven::read_spss("EPOCH_Clean_no_id.sav")  #  the short version data

data <- read.csv("Elevation_Chengdu_4.csv")  #  the short version data

# clean data
data_v <- data %>%
   dplyr::filter(Age > 0 & Age < 50)

ElevName <- c("Eleva1", "Eleva2", "Eleva3", "Eleva4" , "Eleva5", "Eleva6")

PANASName <- c("Active", "Enthusiastic", "Happy", "Excited", "Proud", "Cheerful", "Delighted", "Strong", "Gratitude",
               "Ashamed", "Upset", "Horrified", "Scared", "Irritable", "Nervous","Jittery",  "Guilty", "Angry" )

EnvyName  <- c("Envy1","Envy2","Envy3" , "Envy4" ,"Envy5" ,"Envy6","Envy7","Envy8")

#################################
########  Elevation  ############
#################################
Data_Elev <- data_v %>%
   dplyr::select(ElevName) %>%
   tidyr::drop_na()


Data_Elev_m <- data_v %>%
   dplyr::filter(Gender == 1) %>%  # select male
   dplyr::select(ElevName) %>%     # only elevation 
   tidyr::drop_na()

Data_Elev_f <- data_v %>%
   dplyr::filter(Gender == 2) %>%  # select female
   dplyr::select(ElevName) %>%     # only elevation 
   tidyr::drop_na()

Data_Elev_corr    <- qgraph::cor_auto(Data_Elev)
Data_Elev_corr_m  <- qgraph::cor_auto(Data_Elev_m)
Data_Elev_corr_f  <- qgraph::cor_auto(Data_Elev_f)

network1_Elev     <- bootnet::estimateNetwork(Data_Elev,default = "EBICglasso")
network1_Elev_m   <- bootnet::estimateNetwork(Data_Elev_m,default = "EBICglasso")
network1_Elev_f   <- bootnet::estimateNetwork(Data_Elev_f,default = "EBICglasso")

type <- rep('g', 6)                 # g=gaussian, 6 = number of nodes in the network

# using average layout late
Layout <- qgraph::averageLayout(network1_Elev,network1_Elev_m,network1_Elev_f)

fit1_elev  <- mgm::mgm(na.omit(Data_Elev), type=type, lev=rep(1,6)) 
pred1_elev <- predict(fit1_elev, na.omit(Data_Elev))

fit1_elev_m  <- mgm::mgm(na.omit(Data_Elev_m), type=type, lev=rep(1,6)) 
pred1_elev_m <- predict(fit1_elev_m, na.omit(Data_Elev_m))

fit1_elev_f  <- mgm::mgm(na.omit(Data_Elev_f), type=type, lev=rep(1,6)) 
pred1_elev_f <- predict(fit1_elev_f, na.omit(Data_Elev_f))

# network of all data
pdf("Figure1_elevation.pdf", width= 8, height=6)
network1_Elev <- plot(network1_Elev,
                       title="Elevation", 
                       maximum=.47,
                       labels = TRUE,
                       nodeNames = ElevName,
                       pie = pred1_elev$error$R2,
                       layout = Layout,
                       #groups = gr,
                       #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                       legend.cex = .35,
                       border.width=2, border.color='#555555')       # max=0.47
#qgraph(graph1)
dev.off()

# network of Male
pdf("Figure1_elevation_Male.pdf", width= 8, height=6)
network1_Elev_m <- plot(network1_Elev_m,
                      title="Elevation", 
                      maximum=.47,
                      labels = TRUE,
                      nodeNames = ElevName,
                      pie = pred1_elev_m$error$R2,
                      layout = Layout,
                      #groups = gr,
                      #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                      legend.cex = .35,
                      border.width=2, border.color='#555555')       # max=0.47
#qgraph(graph1)
dev.off()

# network of Male
pdf("Figure1_elevation_Female.pdf", width= 8, height=6)
network1_Elev_f <- plot(network1_Elev_f,
                        title="Elevation", 
                        maximum=.47,
                        labels = TRUE,
                        nodeNames = ElevName,
                        pie = pred1_elev_f$error$R2,
                        layout = Layout,
                        #groups = gr,
                        #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                        legend.cex = .35,
                        border.width=2, border.color='#555555')       # max=0.47
#qgraph(graph1)
dev.off()

#################################
########  PANAS      ############
#################################

Data_PANAS <- data_v %>%
   dplyr::select(PANASName) %>%
   tidyr::drop_na()

Data_PANAS_corr   <- qgraph::cor_auto(Data_PANAS)
network1_PANAS    <- bootnet::estimateNetwork(Data_PANAS,default = "EBICglasso")

fit1_PANAS  <- mgm::mgm(na.omit(Data_PANAS), type=type_panas, lev=rep(1,18)) 
pred1_PANAS <- predict(fit1_PANAS, na.omit(Data_PANAS))

type_panas <- rep('g', 18)                 # g=gaussian, 6 = number of nodes in the network

g_PANAS <- list(c(1:9),c(10:18))

# network of all data of PANAS
pdf("Figure2_PANAS.pdf", width= 12, height=6)
network1_PANAS_p <- plot(network1_PANAS,
                      title="PANAS", 
                      maximum=.47,
                      labels = TRUE,
                      nodeNames = PANASName,
                      pie = pred1_PANAS$error$R2,
                      layout = 'spring',
                      groups = g_PANAS,
                      color = c("#d7191c", "#fdae61"),
                      #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                      legend.cex = .35,
                      border.width=2, border.color='#555555')       # max=0.47
#qgraph(graph1)
dev.off()

#################################
####  Elevation & PANAS  ########
#################################
Elev_PANAS_Name <- c(ElevName,PANASName)
Data_Elev_PANAS <- data_v %>%
   dplyr::select(Elev_PANAS_Name) %>%
   tidyr::drop_na()

Data_Elev_PANAS_corr   <- qgraph::cor_auto(Data_Elev_PANAS)
network1_Elev_PANAS    <- bootnet::estimateNetwork(Data_Elev_PANAS,default = "EBICglasso")

type_Elev_panas <- rep('g', 24)                 # g=gaussian, 6 = number of nodes in the network

fit1_Elev_PANAS  <- mgm::mgm(na.omit(Data_Elev_PANAS), type=type_Elev_panas, lev=rep(1,24)) 
pred1_Elev_PANAS <- predict(fit1_Elev_PANAS, na.omit(Data_Elev_PANAS))


g_Elev_PANAS <- list(c(1:6),c(7:15),c(16:24))

# network of all data of PANAS
pdf("Figure3_Elev_PANAS.pdf", width= 12, height=6)
network1_PANAS_p <- plot(network1_Elev_PANAS,
                         title="Elev_PANAS", 
                         maximum=.47,
                         labels = TRUE,
                         nodeNames = Elev_PANAS_Name,
                         pie = pred1_Elev_PANAS$error$R2,
                         layout = 'spring',
                         groups = g_Elev_PANAS,
                         color = c("#d7191c", "#fdae61","#abd9e9"),
                         #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                         legend.cex = .35,
                         border.width=2, border.color='#555555')       # max=0.47
#qgraph(graph1)
dev.off()


Elev_PANAS_Cent <- centralityPlot(network1_Elev_PANAS) #using the graph from the data

### plot figure s2: all centrality measures;
pdf("Figure_4_centrality_Elev_PANAS.pdf", width=6, height=10) 
Elev_PANAS_Cent 
dev.off()

#################################
####  Elevation & Envy   ########
#################################

Elev_Envy_Name <- c(ElevName,EnvyName)
Data_Elev_Envy <- data_v %>%
   dplyr::select(Elev_Envy_Name) %>%
   tidyr::drop_na()

Data_Elev_Envy_corr   <- qgraph::cor_auto(Data_Elev_Envy)
network1_Elev_Envy    <- bootnet::estimateNetwork(Data_Elev_Envy,default = "EBICglasso")

type_Elev_Envy <- rep('g', 14)                 # g=gaussian, 6 = number of nodes in the network

fit1_Elev_Envy  <- mgm::mgm(na.omit(Data_Elev_Envy), type=type_Elev_Envy, lev=rep(1,14)) 
pred1_Elev_Envy <- predict(fit1_Elev_Envy, na.omit(Data_Elev_Envy))


g_Elev_Envy <- list(c(1:6),c(7:14))

# network of all data of PANAS
pdf("Figure5_Elev_Envy.pdf", width= 12, height=6)
network1_Envy_p <- plot(network1_Elev_Envy,
                         title="Elev_Envy", 
                         maximum=.47,
                         labels = TRUE,
                         nodeNames = Elev_Envy_Name,
                         pie = pred1_Elev_Envy$error$R2,
                         layout = 'spring',
                         groups = g_Elev_Envy,
                         color = c("#d7191c", "#fdae61"),
                         #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
                         legend.cex = .35,
                         border.width=2, border.color='#555555')       # max=0.47
#qgraph(graph1)
dev.off()


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
      dplyr::select(t2name)

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
mean(pred2$error$R2)     # .48775

# the node predictability of the central items:
pred1_all$errors$R2[20]  # H4 0.654
pred1_all$errors$R2[2]   # E2 0.53
pred1_all$errors$R2[12]  # O4 0.52

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

### plot figure S1, graph of test and re-test graph
pdf("Fig_S1_retest_T1_T2.pdf", width= 24, height= 9)
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
# ---------- 4. Centrality of the estimated networks ------------------------------------
# ---------------------------------------------------------------------------------------

### Centrality

All_cp <- centralityPlot(network1_all_G) #using the graph from the data

### plot figure s2: all centrality measures;
pdf("Figure_S2_centrality_all.pdf", width=10, height=7) 
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

# transfom the node strength data
cent <- as.data.frame(scale(centrality(network1_all_G)$InDegree))
cent <- dplyr::mutate(cent, id = rownames(cent))
colnames(cent) <- c("1", "EPOCH_Item")
cent_long <- reshape2::melt(cent, id="EPOCH_Item")

# plot Figure 2: the node strength as . 
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
# edge weight accuracy for the whole network
boot1a <- bootnet::bootnet(network1_all, nBoots = 1000, nCores = 8)

# centrality stability for the whole network
boot1b <- bootnet::bootnet(network1_all, nBoots = 1000, type = "case",  nCores = 8)

# edge weight accuracy for the re-test network 
#boot2a <- bootnet::bootnet(network2, nBoots = 1000, nCores = 8)

# centrality stability for the re-test network
#boot2b <- bootnet::bootnet(network2, nBoots = 1000, type = "case",  nCores = 8)

### Plot edge weight CI as Figure S3
pdf("FigS3_edge_weight_ci.pdf")
plot(boot1a, labels = FALSE, order = "sample")   # all data
# plot(boot2a, labels = FALSE, order = "sample") # replication data
dev.off()

### Edge weights diff test, plotted as Figure S4
pdf("FigS4_Edge_weight_diff.pdf")
plot(boot1a, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
#plot(boot2a, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
dev.off()


### Plot centrality stability, as Figure S5
pdf("FigS5_centrality_stablity.pdf") 
plot(boot1b)
#plot(boot2b)
dev.off()

### Centrality stability coefficient (CS-coefficient)
cs1 <- bootnet::corStability(boot1b)  
# cs2 <- bootnet::corStability(boot2b) # re-test network
cs <- matrix(NA,1,3)
cs[1,1:3] <- round(cs1, digits=2)
# cs[2,1:3] <- round(cs2, digits=2)
colnames(cs) <- c("Betweenness", "Closeness", "Strength")
#rownames(cs) <- c("All_data_T1", "Rep_Data_T2")
rownames(cs) <- c("All_data_T1")
sink("TableS1.txt")
cs
sink()


### Centrality diff test as figure S6
pdf("FigS6_centrality_diff.pdf")
plot(boot1a, "strength", order="sample", labels=FALSE) 
#plot(boot2a, "strength", order="sample", labels=FALSE) 
dev.off()

### Conclusion: all networks have CS coefficients for node strength > 0.5; looks good!

### Save output
save(boot1a, file = "boot1a.Rdata")
save(boot1b, file = "boot1b.Rdata")
#save(boot2a, file = "boot2a.Rdata")
#save(boot2b, file = "boot2b.Rdata")

##### test- retest network comparison
#####
dft1_c1  <- cor(df1_r)
dft1_c1b <- qgraph::cor_auto(df1_r)
cor(dft1_c1[lower.tri(dft1_c1)], dft1_c1b[lower.tri(dft1_c1b)], 	method="spearman") #0.9936

dft2_c1  <- cor(df2)
dft2_c1b <- qgraph::cor_auto(df2)
cor(dft2_c1[lower.tri(dft2_c1)], dft2_c1b[lower.tri(dft2_c1b)], 	method="spearman") #0.99164

# since they are very similar, we can probably use the NetworkComparisonTest::NCT() 
# that is only validated for Pearson, not poylchoric correlations

### NCT topology
# the colnames for test at T2
t3name <- c("E1", "E2", "E3", "E4",
            "P1", "P2", "P3", "P4", 
            "O1", "O2", "O3", "O4",
            "C1", "C2", "C3", "C4",
            "H1", "H2", "H3", "H4")
colnames(df1_r) <- colnames(df2) <- t3name  # re-name the column name of these two data set
set.seed(100)
compare_12 <- NetworkComparisonTest::NCT(df1_r,df2, it=5000, binary.data=FALSE, 
                                         test.edges=TRUE, edges='all', progressbar=TRUE)

compare_12$nwinv.pval     # output p-value: 0.0526

### quantification of differences: count significantly different edges (total number )
sum(compare_12$einv.pvals$"p-value" < 0.05) # 1 of 190, < 0.5%

### NCT global strength
compare_12$glstrinv.pval  # 0.01

save(compare_12, file = "compare_test-retest.Rdata")

### end ###
