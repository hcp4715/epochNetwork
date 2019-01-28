##############################################################################
###                                                                        ###
###                R code for Network of Positive Tratis                   ###
###                                                                        ###
##############################################################################
# Author   Date(Y-M-D)  Log of change
# =======  ===========  ==============
# hcp      18.10.05     revised code based on Bellet et al 2018.
#
#
# 
##### Part 1: read the .sav file ##########

# remove previous variables in memory
rm(list = ls())     # remove all variables
curWD <- dirname(rstudioapi::getSourceEditorContext()$path) #Get the directory ofcurrent script
setwd(curWD)

Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre7') # for 64-bit version
Sys.setlocale("LC_ALL", "English")  # set local encoding to English
Sys.setenv(LANG = "en") # set the feedback language to English

# load libraries
library(tidyverse)
library(haven)
library(bootnet)
library(qgraph)

if (!require(networktools)) {install.packages("networktools",repos = "http://cran.us.r-project.org"); require(networktools)}
library("networktools")
if (!require(xlsx)) {install.packages("xlsx",repos = "http://cran.us.r-project.org"); require(networktools)}
library("xlsx")

#source('Bellet_Supplemental_Script_S2.R')

#spss2date <- function(x) as.Date(x/86400, origin = "1582-10-14")
#df.raw <- foreign::read.spss("EPOCHcombinedB.sav", reencode='utf-8',to.data.frame=TRUE)
#df.raw$name <- as.character(df.raw$name)
#df.raw$surveyDate <- spss2date(df.raw$surveyDate)
#str(df.raw)
#summary(df.raw)

# load the reduced dataset
#df.s <- haven::read_spss("EPOCHreducedForR2.sav")  #  the short version data
df.l <- haven::read_spss("EPOCHcombinedB 10.18.sav")     #  the long version data

###### estimate network model for short version ####
# df.s <- df.s[,1:50]    # remove the summary score

# the colnames scales
epochName <- c(c(paste("E",1:4, sep="")),c(paste("P",1:4,sep = '')),
            c(paste("O",1:4,sep = '')),c(paste("C",1:4,sep = '')),c(paste("H",1:4,sep = '')))

# names for deppression
depresName <- c(paste("D",1:5,sep = ''))

# names for anxiety
anxName <- c(paste("A",1:5,sep = ''))

# names for health
healthName <- c(paste("Health",1:5,sep = ''))

# names for intro
introName <- c(paste("intro",1:5,sep = ''))

# names for exter
exterName <- c(paste("exter",1:5,sep = ''))

# names for Ident
IdentName <- c(paste("Ident",1:5,sep = ''))

# names for resilience
resNames <- df.l %>% 
  dplyr:: select(starts_with("res")) %>% colnames()
resName <- resName1[c(1:6)]

## data with all 6 item
df.res2 <- df.l %>%
  dplyr::filter(!is.na(res1) & 
                  !is.na(res2) &
                  !is.na(res3) &
                  !is.na(res4) &
                  !is.na(res5) &
                  !is.na(res6)) %>%            # eliminate the columns without resilience data
  dplyr::select_if(~sum(!is.na(.)) > 0) %>%  # eliminate columns that only have NAs
  dplyr::select(-c(epochE,epochP,epochO,epochC,epochH,epochAll,resilience))

df.res2 %>% 
  select_if(function(x) any(is.na(x))) %>% 
  summarise_each(funs(sum(is.na(.)))) -> extra_NA

df.res2_2 <- df.res2 %>%
  dplyr::filter(!is.na(age)) %>%             # has age info
  dplyr::select_if(~sum(!is.na(.)) > 0) %>%  # no all NA
  dplyr::filter(!(age <0 | age > 20)) %>%       # bizard age
  dplyr::select(-c(name,set,surveyDate,birthday,grade,gradeS,class,classS,major,major_A,school,
                   Schoolname, schooltype,provinces, date,id810,object,object_A,
                   depress,anxiety,Health,res2,res4,res6,
                   Intro, Ident,Integ,Exter,
                   GM2,GM4,GM6,GM8,GM))
corrplot::corrplot(cor(df.res2_2[,-c(1:4)]))
network1 <- bootnet::estimateNetwork(df.res2_2[,-c(1:4)],default = "EBICglasso")

plot(network1,
     #title="Resilience", 
     maximum=.47,
     #labels = TRUE,
     #nodeNames = resName1,
     #pie = pred1_all$error$R2,
     layout = 'spring')

summary(df.res2_2)

## data with all 4 item
df.res1 <- df.l %>%
  dplyr::filter(!is.na(res1) & 
                  !is.na(res2) &
                  !is.na(res3) &
                  !is.na(res4)) %>%            # eliminate the columns without resilience data
  dplyr::select_if(~sum(!is.na(.)) > 0) %>%  # eliminate columns that only have NAs
  dplyr::select(-c(epochE,epochP,epochO,epochC,epochH,epochAll,resilience))


df.res <- df.res2[,resName] %>%
  dplyr::mutate(res2_r = 6 - res2,
                  res4_r = 6 - res4,
                  res6_r = 6 - res6) %>%
  dplyr::select(res1,res2_r,res3,res4_r,res5,res6_r)
str(df.res)
summary(df.res)
psych::omega(df.res)
head(df.res)
corrplot::corrplot(cor(df.res))
df.res_corr <- qgraph::cor_auto(df.res)

network1_res <- bootnet::estimateNetwork(df.res,default = "EBICglasso")

plot(network1_res,
     title="Resilience", 
     maximum=.47,
     labels = TRUE,
     nodeNames = resName1,
     #pie = pred1_all$error$R2,
     layout = 'spring',
     #groups = gr,
     #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"),
     #color = c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6", "#aa00ff"),
     legend.cex = .35,
     border.width=2, border.color='#555555')       # max=0.47

# names for GM
GMName <- c(paste("GM",1:8,sep = ''))

# names for FR
FRName <- c(paste("FR",1:18,sep = ''))

# names for MR
MRName <- c(paste("MR",1:18,sep = ''))

# names scheng
colnames(df.l)[colnames(df.l) == 'SchEng1'] <- 'scheng1'
colnames(df.l)[colnames(df.l) == 'schEng5'] <- 'scheng5'
schengName <- c(paste("scheng",1:5,sep = ''))

# names for cope
copeName <- c(paste("cope",1:5,sep = ''))

# names for emp
empName <- c(paste("emp",1:3,sep = ''))

# names for SelAwar
selAwarName <- c(paste("SelAwar",1:5,sep = ''))

# names for belong
belongName <- c(paste("belong",1:4sep = ''))

# names for Grit
gritName <- c(paste("Grit",1:3,sep = ''))

# names for seficy
seficyName <- c(paste("seficy",1:10,sep = ''))

# names for MAP (??)
MAPName <- c(paste("MAP",1:4,sep = ''))
PAPName <- c(paste("PAP",1:4,sep = ''))
MAVName <- c(paste("MAV",1:4,sep = ''))
PAVOName <- c(paste("PAVO",1:3,sep = ''))

# names for Belongm
BelongmName <- c(paste("Belongm",1:4,sep = ''))

# names for Belongc
BelongcName <- c(paste("Belongc",1:4,sep = ''))

# names for Belongt
BelongtName <- c(paste("Belongt",1:4,sep = ''))

# names for Belongfri
BelongfriName <- c(paste("Belongfri",1:4,sep = ''))

# names for interperson
interpersonName <- c(paste("interperson",1:5,sep = ''))

# other names
relName <- c('schpf','TeaRel','ClaRel','ParRel','initiation','NegAss','exposure','ES','InterSk')


epochItems <- c("E-forgetTime", "E-focus",   "E-flow",      "E-learn",
               "P-beginEnd",   "P-finish",  "P-plan",      "P-hardwork",
               "O-futur",      "O-best",    "O-goodthing", "O-solution",
               "C-share",      "C-support", "C-caring",    "C-friends",
               "H-happy",      "H-fun",     "H-life",      "H-joy")

# clean the name variable, may need new code:
tmp <- df.l$name %>%
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


### restart to clean the data
### Basic information:
### dataset of cleaned data 
###  dataset ID: 249  810 1340 1410 1854 2293 2325 3588 4543 4843 5111 
### sample size: 234  778 1340 1279 1737 2129 1664 2271 2607 1322 2493

### dataset of uncleaned data:
###  dataset ID: 249  810 1340 1410 1854 2293 2325 3588 4543 4843 5111 
### sample size: 234  778 1340 1291 1762 2161 1716 3038 3826 4338 4321

### This means: dataset 249, 810 and 1340 are clean

### recode the school ID
# read the file with school ID
# read the coded file
Sys.setlocale(category = "LC_ALL", locale = "chs") #cht for traditional Chinese, etc.
# extract the useful columns
library(xlsx)
UniSchName  <- read.xlsx2('UniqueSchoolName1101.xlsx', 'UniqueSchoolName', header = T,stringsAsFactors=FALSE)
UniSchAName <- read.xlsx2('Unique_School_A_Name1101.xlsx', 'Unique_School_A_Name', header = T,stringsAsFactors=FALSE)
#UniSchName  <- read.csv('UniqueSchoolName1101.csv', sep = ',',header = T,encoding = "UTF-8", stringsAsFactors=FALSE)
#UniSchAName <- read.csv('Unique_School_A_Name1101.csv', sep = ',', quote = "",header = T,encoding = "UTF-8", stringsAsFactors=FALSE)

#UniqueSchoolName <- UniqueSchoolName[,1:5]
# rename the columns
#colnames(UniqueSchoolName) <- c("School","Freq","SchoolID","City","CityID")
#colnames(Unique_School_A_Name) <- c("School","Freq","SchoolID","City","CityID")
# clear the strings by extracting the characters between " "
#UniqueSchoolName$SchoolNew <- gsub(".*[\"]([^.]+)[\"].*", "\\1", UniqueSchoolName$School)
#Unique_School_A_Name$SchoolNew <- gsub(".*[\"]([^.]+)[\"].*", "\\1", Unique_School_A_Name$School)

# Match the School name and asign the school ID
# df.l$schoolID <- UniSchName[match(df.l$school, UniSchName$school),3]
#df.l$schoolID <- as.character(df.l$schoolID)
df.l$schoolID[is.na(df.l$schoolID)] <- 
      UniSchAName[match(df.l$school_A[is.na(df.l$schoolID)], UniSchAName$school),3]

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

### temporarily forget the problem
resilName <- c("res1","res3","res5", "res2_reverse", "res4_reverse","res6_reverse")
df.1 <- subset(df.l, set %in% c(249,810,1340))
resScale <- df.1 %>%
      select(epochName,resilName) #%>%
      #complete.cases()
