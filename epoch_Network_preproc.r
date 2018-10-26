
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
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre7') # for 64-bit version
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
if (!require(xlsx)) {install.packages("xlsx",repos = "http://cran.us.r-project.org"); require(networktools)}
library("xlsx")

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

epochName <- c("E-forgetTime", "E-focus",   "E-flow",      "E-learn",
               "P-beginEnd",   "P-finish",  "P-plan",      "P-hardwork",
               "O-futur",      "O-best",    "O-goodthing", "O-solution",
               "C-share",      "C-support", "C-caring",    "C-friends",
               "H-happy",      "H-fun",     "H-life",      "H-joy")

# reliability of two test
alpha.t1 <- psych::alpha(df.s[,t1name]) 
print(alpha.t1$total)  # 0.936

alpha.t2 <- psych::alpha(df.s[,t2name]) 
print(alpha.t2$total)  # 0.935

# estimate the network for Time point 1
cor.epochNet.t1 <- invisible(qgraph::cor_auto(df.s[,t1name]))
glasso_epochNet.t1 <- qgraph::EBICglasso(cor.epochNet.t1, n=dim(df.s[,t1name])[1], gamma=0.5)

epochNet.t1 <- qgraph::qgraph(glasso_epochNet.t1, layout="spring", labels=epochName, vsize=9,
                label.cex = c(rep(.7,10),.6,.7,.6),
                label.font = c(rep(1,20)),
                label.color = c(rep('blue', 4),rep('red',4),rep(1,4),rep('#993300',4),rep("#009999",4)),
                label.scale = F, DoNotPlot=F)

pdf("epochNet.t1.pdf", width = 3.8, height = 4)
epochbw <- makeBW(epochNet.t1)
dev.off()

# estimate the network for Time point 2
cor.epochNet.t2 <- invisible(qgraph::cor_auto(df.s[,t2name]))
glasso_epochNet.t2 <- qgraph::EBICglasso(cor.epochNet.t2, n=dim(df.s[,t2name])[1], gamma=0.5)

epochNet.t2 <- qgraph::qgraph(glasso_epochNet.t2, layout="spring", labels=epochName, vsize=9,
                               label.cex = c(rep(.7,10),.6,.7,.6),
                              label.font = c(rep(1,20)),
                              label.color = c(rep('blue', 4),rep('red',4),rep(1,4),rep('#993300',4),rep("#009999",4)),
                              label.scale = F, DoNotPlot=F)

pdf("epochNet.t2.pdf", width = 3.8, height = 4)
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

for (ii in 1:nrow(tmp2)) {
      #a <- wk[tmp[ii]]
      new <- c()
      a_tmp <- strsplit(tmp[ii], ' ')[[1]]
      for (jj in 1:length(a_tmp)) {
            if (stringr::str_length(a_tmp[jj]) > 0) {
                  new = c(new, a_tmp[jj])
            }
      }
      if (length(new) == 0) {
            tmp2$name[ii] <- "NA"
            tmp2$school1[ii] <- NA
      } else if (length(new) > 0 & (length(new) <= 2)) {
            tmp2$name[ii] <- new[1]
            tmp2$school1[ii] <- new[2]
      } else if (length(new) == 3) {
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

tmp3 <- tmp2[, 1:2]
tmp3$len <- stringr::str_length(tmp3$name)
tmp4 <- tmp3[tmp3$len > 12, ]
tmp4$name <- strsplit(tmp4$name, '??')[[1]][1]
tmp4$name <- summary(tmp3$len)
a <- df.raw$name[1]
new <- c()
a_tmp <- strsplit(a,' ')[[1]]
for (i in 1:length(a_tmp)) {
      if (stringr::str_length(a_tmp[i]) > 0) {
            print('>0')
            new = c(new, a_tmp[i])
      }
      
      print(c("length:", length(new)))
}
new


for (ii in 1:length(tmp)) {
      #a <- wk[tmp[ii]]
      if (length(a) == 1) {
            tmp2$name[ii] <- a[1]
      } else if (length(a) == 2) {
            tmp2$school1[ii] <- a[2]
      } else if (length(a) >= 3) {
            tmp2$school1[ii] <- a[2]
            print(paste('the third string is:', a[3]))
            
      }
      
}


tmp <- with(df.l,table(name))
repeatedName <- data.frame(count = tmp[tmp >2])
df.l$name[which(nchar(df.l$name, type = "chars") == 19)] # " jiao   ming    hui"
df.l$name[which(nchar(df.l$name, type = "chars") == 17)] # "tian   zhi  xiang"
df.l$name[which(nchar(df.l$name, type = "chars") == 16)] # "huahgzhdhjuyhjfd" "wang  su    meng"
df.l$name[which(nchar(df.l$name, type = "chars") == 15)] # "            王紫漩" "....3.........." 
df.l$name[which(nchar(df.l$name, type = "chars") == 14)] # "陈思航。Jack Smith" "0zhangmengting"  
df.l$name[which(nchar(df.l$name, type = "chars") == 13)] 
# "yang  rui  xi"  "徐子玥  成都市 康河小学" "李文俊   成都市康河小学" "钟启伟   成都市康河小学"
# "马丹丹   成都市康河小学" "sunweixiangvb" 
df.l$name[which(nchar(df.l$name, type = "chars") == 12)] 
df.l$name[which(nchar(df.l$name, type = "chars") == 11)] 

schlName1 <- data.frame(with(df.l,table(school)))
#write.csv(schlName1, 'UniqueSchoolName.csv',row.names = F)
write.table(schlName1, 'UniqueSchoolName_1_raw.csv',quote = FALSE, sep = ',',row.names = F)
schlName2 <- data.frame(with(df.l,table(school_A)))
write.table(schlName2, 'Unique_School_A_Name_1_raw.csv',quote = FALSE, sep = ',',row.names = F)
#write.csv(schlName2, 'Unique_School_A_Name.csv',row.names = F)

# read the file with school ID
# read the coded file

UniqueSchoolName <- read.csv('UniqueSchoolName.csv', header = T,encoding = "UTF-8", stringsAsFactors=FALSE)
Unique_School_A_Name <- read.csv('Unique_School_A_Name.csv', header = T,encoding = "UTF-8", stringsAsFactors=FALSE)

Sys.setlocale(category = "LC_ALL", locale = "chs") #cht for traditional Chinese, etc.
# extract the useful columns
UniqueSchoolName <- UniqueSchoolName[,1:5]
# rename the columns
colnames(UniqueSchoolName) <- c("School","Freq","SchoolID","City","CityID")
colnames(Unique_School_A_Name) <- c("School","Freq","SchoolID","City","CityID")
# clear the strings by extracting the characters between " "
UniqueSchoolName$SchoolNew <- gsub(".*[\"]([^.]+)[\"].*", "\\1", UniqueSchoolName$School)
Unique_School_A_Name$SchoolNew <- gsub(".*[\"]([^.]+)[\"].*", "\\1", Unique_School_A_Name$School)

# Match the School name and asign the school ID
df.l$schoolID <-UniqueSchoolName[match(df.l$school, UniqueSchoolName$SchoolNew),3]
df.l$schoolID[is.na(df.l$schoolID)] <- Unique_School_A_Name[match(df.l$school_A[is.na(df.l$schoolID)], Unique_School_A_Name$SchoolNew),3]

# Get the unasigned school name
UniqueSchoolName2 <- data.frame(with(df.l[is.na(df.l$schoolID),],table(school)))
write.table(UniqueSchoolName2, 'UniqueSchoolName2.csv',quote = FALSE, sep = ',',row.names = F)
Unique_School_A_Name2 <- data.frame(with(df.l[is.na(df.l$schoolID),],table(school_A)))
write.table(Unique_School_A_Name2, 'Unique_School_A_Name2.csv',quote = FALSE, sep = ',',row.names = F)

# 2nd run
UniqueSchoolName2 <- read.csv('UniqueSchoolName2_coded.csv', sep = ',',header = T,encoding = "UTF-8", stringsAsFactors=FALSE)
colnames(UniqueSchoolName2) <- c("School","Freq","SchoolID","City","CityID")
df.l$schoolID[is.na(df.l$schoolID)] <- UniqueSchoolName2[match(df.l$school_A[is.na(df.l$schoolID)], UniqueSchoolName2$School),3]

Unique_School_A_Name2 <- read.csv('Unique_School_A_Name2_coded.csv', sep = ',',header = T,encoding = "UTF-8", stringsAsFactors=FALSE)
colnames(UniqueSchoolName2) <- c("School","Freq","SchoolID","City","CityID")
df.l$schoolID[is.na(df.l$schoolID)] <- UniqueSchoolName2[match(df.l$school_A[is.na(df.l$schoolID)], UniqueSchoolName2$School),3]
