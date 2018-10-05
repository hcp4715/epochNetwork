
##############################################################################
################## R code for Network of Positive Tratis #####################
##############################################################################
#
# Author   Date(Y-M-D)  Log of change
# =======  ===========  ==============
# hcp      18.10.05     revised code based on Bellet et al 2018.
#
#
# 
##### Part 1: read the .sav file ##########

# remove previous variables in memory
Sys.setlocale("LC_ALL", "English")  # set local encoding to English
Sys.setenv(LANG = "en") # set the feedback language to English

rm(list = setdiff(ls(), lsf.str())) # remove all variables except functions

# load libraries
library(tidyverse)
library(haven)
library('bootnet')
library('qgraph')

if (!require(networktools)) {install.packages("networktools",repos = "http://cran.us.r-project.org"); require(networktools)}
library("networktools")

source('Bellet_Supplemental_Script_S2.R')

#spss2date <- function(x) as.Date(x/86400, origin = "1582-10-14")
#df.raw <- foreign::read.spss("EPOCHcombinedB.sav", reencode='utf-8',to.data.frame=TRUE)
#df.raw$name <- as.character(df.raw$name)
#df.raw$surveyDate <- spss2date(df.raw$surveyDate)
#str(df.raw)
#summary(df.raw)

# load the reduced dataset
df.s <- haven::read_spss("EPOCHreducedForR2.sav")  #  the short version data
df.l <- haven::read_spss("EPOCHcombinedB.sav")     #  the long version data

###### estimate network model for short version ####
df.s <- df.s[,1:50]    # remove the summary score

# the colnames for test at T1
t1name <- c("E1.1", "P1.1", "O1.1", "C1.1", "H1.1",
            "E2.1", "P2.1", "O2.1", "C2.1", "H2.1",
            "E3.1", "P3.1", "O3.1", "C3.1", "H3.1",
            "E4.1", "P4.1", "O4.1", "C4.1", "H4.1")

# the colnames for test at T2
t2name <- c("E1.2", "P1.2", "O1.2", "C1.2", "H1.2",
            "E2.2", "P2.2", "O2.2", "C2.2", "H2.2",
            "E3.2", "P3.2", "O3.2", "C3.2", "H3.2",
            "E4.2", "P4.2", "O4.2", "C4.2", "H4.2")

# reliability of two test
alpha.t1 <- psych::alpha(df.s[,t1name]) 
print(alpha.t1$total)  # 0.936

alpha.t2 <- psych::alpha(df.s[,t2name]) 
print(alpha.t2$total)  # 0.935

# estimate the network 
cor.epochNet.t1 <- invisible(qgraph::cor_auto(df.s[,t1name]))
glasso_epochNet.t1 <- qgraph::EBICglasso(cor.epochNet.t1, n=dim(df.s[,t1name])[1], gamma=0.5)

epochNet.t1 <- qgraph::qgraph(glasso_epochNet.t1, layout="spring", labels=t1name, vsize=9,
                label.cex=c(rep(.7,10),.6,.7,.6),
                label.font=c(rep(1,10),2,1,2),
                label.color=c(rep(1, 10),"blue", 1,"blue"),
                label.scale=F, DoNotPlot=F)

pdf("epochNet.t1.pdf", width = 3.8, height = 4)
epochbw <- makeBW(epochNet.t1)
dev.off()

# expected influence
EI_epoch.t1 <- networktools::expectedInf(glasso_epochNet.t1)
pdf("epoch.t1_EI.pdf", width = 5)
plot(EI_epoch.t1$step1, order="value", zscore=F, yaxt = "n")
dev.off()

####Edge Weight and EI stability
set.seed(123)
net <- bootnet::estimateNetwork(df.s[,t1name], default="EBICglasso")

net_boot <- bootnet_flex(net, statistics=c("edge", "expectedInf"), nBoots=1000, type="case", caseN = 50)

CorStabEpochT1 <- corStability(net_boot, statistics=c("edge", "expectedInf")) # correlation stability of the bootnet




epochNet.t1 <- df.s[,t1name] %>%
      estimateNetwork(default = "EBICglasso")

plot(epochNet.t1,layout = 'spring',labels = TRUE)
centralityPlot(epochNet.t1)

epochNet.t2 <- estimateNetwork(df[,t2name],default = "EBICglasso")
plot(epochNet.t2,layout = 'spring',labels = TRUE)
centralityPlot(epochNet.t2)

# clean the name variable, may need new code:
tmp <- df.raw$name %>%
        as.character()
tmp2 <- setNames(data.frame(matrix(ncol = 3, nrow = length(tmp) )), c("name", "school1",'unknown'))        

for (ii in 1:nrow(tmp2)){
        #a <- wk[tmp[ii]]
        new <- c()
        a_tmp <- strsplit(tmp[ii],' ')[[1]]
        for (jj in 1:length(a_tmp)){
                if(stringr::str_length(a_tmp[jj])>0){
                        new = c(new,a_tmp[jj])
                }
        }
        if(length(new) == 0){
                tmp2$name[ii] <- "NA"
                tmp2$school1[ii] <- NA
        } else if (length(new)>0 & (length(new)<=2)){
                tmp2$name[ii] <- new[1]
                tmp2$school1[ii] <- new[2]
        } else if (length(new) == 3){
                tmp2$name[ii] <- new[1]
                tmp2$school1[ii] <- new[2]
                tmp2$unknown[ii] <- new[3]
        }
        
}

# if all three columns of tmp2 is not na, paste them
tmp2$name[rowSums(!is.na(tmp2)) == 3] <- paste(tmp2$name[rowSums(!is.na(tmp2)) == 3],
                                               tmp2$school1[rowSums(!is.na(tmp2)) == 3],
                                               tmp2$unknown[rowSums(!is.na(tmp2)) == 3],
                                               sep = '')

tmp2[rowSums(!is.na(tmp2)) == 3,c('school1','unknown')] <- NA
#tmp2$unknown[rowSums(!is.na(tmp2)) == 3] <- NA

tmp3 <- tmp2[,1:2]
tmp3$len <- stringr::str_length(tmp3$name)
tmp4 <- tmp3[tmp3$len>12,]
tmp4$name <- strsplit(tmp4$name,'??')[[1]][1]
tmp4$name <-
summary(tmp3$len)
a <- df.raw$name[1]
new <- c()
a_tmp <- strsplit(a,' ')[[1]]
for (i in 1:length(a_tmp)){
        if(stringr::str_length(a_tmp[i])>0){
                print('>0')
                new = c(new,a_tmp[i])
        }
        
        print(c("length:", length(new)))
}
new


for (ii in 1:length(tmp)){
        #a <- wk[tmp[ii]]
        if (length(a) == 1){
                tmp2$name[ii] <- a[1]
        }else if (length(a) == 2){
                tmp2$school1[ii] <- a[2]
        }else if (length(a) >=3){
                tmp2$school1[ii] <- a[2]
                print(paste('the third string is:', a[3]))
                
        }
        
}

