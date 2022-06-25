if(.Platform$OS.type == "unix") {
  ROOT<-"/home/haseong/dev/ecoli_growth"
  source("/home/haseong/dev/lib/R/fluorescence_analysis.R")
} else {
  ROOT<-"C:\\mydocs\\2019\\dev\\standardization"
  source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")
}


library(readxl)
library(ggplot2)
library(scales)
library(ranger)
library(reshape2)
library(caret)
library(e1071)
require(readxl)


Victor.ROOT <- paste(ROOT, "victor", sep="\\")
setwd(ROOT)

DB.file.name <- "03-Cell-culture-DB_20190430.xlsx"
seed.list <- as.data.frame(read_excel(DB.file.name, sheet = 1))
main.list <- as.data.frame(read_excel(DB.file.name, sheet = 2))
measure.list <- as.data.frame(read_excel(DB.file.name, sheet = 3))

### ============================================================================================
### ============================================================================================
### multi substrates


## ----- user input -----
sel.main <- subset(main.list, Container=="14 mL round bottom tube" & Date == "20180110" ) ## 1111, 1116 error as chloro- chems were diffused
sel.main.merged <- merge(seed.list, sel.main, by.x="ID", by.y="SeedID1")
head(sel.main.merged)
colnames(sel.main.merged)[19]<-"MID"
lab1 <- "595nm_kk (A)"
lab2 <- "EGFP_sulim (Counts)"
lab3 <- "mRFP-sulim(label) (Counts)"
user.detector <- "victor"
sel.inducers <- "phenol"

ii<-which(sel.main.merged[,"ID"]==82)
subs.names1 <- c("phenol", "AHL")
#dev.names <- c("ph-GESS-v4")
## ----------------------

#subset(sel.main.merged, select=c(SeedID, Date.x, plasmid_id, device_name, ))
tmp <- merge(measure.list[,1:13], sel.main.merged, by.x="MainID", by.y="MID")

sel.measures <- subset(tmp, Detector==user.detector)
sel.file.names <- names(table(sel.measures[,"Filename"]))
sel.dat <- cbind(sel.measures, od=rep(0, nrow(sel.measures)), gfp=rep(0, nrow(sel.measures)), rfp=rep(0, nrow(sel.measures)))

for(i in 1:length(sel.file.names)){
  tmp <- as.data.frame(read_excel(paste(Victor.ROOT, sel.file.names[i], sep="\\"), sheet = 1))
  row.id <- substr(tmp[,3], 1, 1)
  col.id <- as.numeric(substr(tmp[,3], 2, 3))
  seli <- which(sel.dat[,"Filename"]==sel.file.names[i])
  o <- match(paste(sel.dat[seli,"Row.y"], sel.dat[seli,"Col.y"]), paste(row.id, col.id))
  if(sum(names(tmp)==lab1)){
    od.val <- tmp[,6]
    sel.dat[seli, "od"] <- od.val[o]
  }
  if(sum(names(tmp)==lab2)){
    gfp.val <- tmp[,8]  
    sel.dat[seli, "gfp"] <- gfp.val[o]
  }
  if(sum(names(tmp)==lab3)){
    rfp.val <- tmp[,10]  
    sel.dat[seli, "rfp"] <- rfp.val[o]
  }
  #fl.data <- data.frame(sel.file.names[i], row.id, col.id, od.val, fl.val)  
  cat(i, "/", length(sel.file.names), "\n");flush.console()
}

## selection
sel.dat2 <- subset(sel.dat, (ID.y %in% c(82, 85)))

## ordering
sel.dat2[,"Inducer_1"] <- factor(sel.dat2[,"Inducer_1"], levels=subs.names1)
sel.dat2[,"Inducer_conc_1"] <- as.numeric(sel.dat2[,"Inducer_conc_1"])

## phenol and etc.
newdata <- subset(sel.dat2, select=c(Date, device_name, Inducer_1, Inducer_conc_1, gfp, rfp, od))





br.lab <- names(table(newdata[,"Inducer_conc_1"]))
#br.lab <- c(0, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000)
## break tip 
br <- as.numeric(br.lab)
## set zero to min value to avoid errors in regression fitting and visualization
br[1] <- br[1] + 0.01
#br <- c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000)

newdata[newdata[,"Inducer_conc_1"]==0,"Inducer_conc_1"] <- 0.01
newdata2 <- cbind(newdata, gfp.od=newdata[,"rfp"]/newdata[,"od"])
newdata2[,1] <- as.factor(newdata2[,1])
newdata.mlt <- melt(newdata2, id=c("Date", "device_name", "Inducer_1", "Inducer_conc_1"), measure="gfp.od")

library(dplyr)


### phenol concentration

newdata.mlt1 <- subset(newdata.mlt, Inducer_1=="phenol")
datnew1 <- data.frame(newdata.mlt1 %>% group_by(Inducer_conc_1) %>% summarise(mean=mean(value), sd = sd(value)))
datnew1$sd[9] <- datnew1$sd[9]/10
ggplot(datnew1, aes(x=factor(Inducer_conc_1), y=mean, width=.7)) +
  geom_bar(stat="identity", position="dodge", color="black", fill="#C8102E") +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), position=position_dodge(.7), width=0.3) +
  labs(y="RFP/OD", x=expression(paste("Phenol concentration (", mu, "M)"))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title= element_blank(), 
        legend.background = element_blank(),
        legend.key.size = unit(2, 'lines')
  ) 
 
  
newdata2 <- cbind(newdata, gfp.od=newdata[,"gfp"]/newdata[,"od"])
newdata2[,1] <- as.factor(newdata2[,1])
newdata.mlt <- melt(newdata2, id=c("Date", "device_name", "Inducer_1", "Inducer_conc_1"), measure="gfp.od")
newdata.mlt2 <- subset(newdata.mlt, Inducer_1=="AHL")
datnew2 <- data.frame(newdata.mlt2 %>% group_by(Inducer_conc_1) %>% summarise(mean=mean(value), sd = sd(value)))
  
ggplot(datnew2, aes(x=factor(Inducer_conc_1), y=mean, width=.7)) +
    geom_bar(stat="identity", position="dodge", color="black", fill="#99cc00") +
    geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), position=position_dodge(.7), width=0.3) +
    labs(y="GFP/OD", x=expression(paste("AHL concentration (", mu, "M)"))) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey90"),
          panel.background = element_blank(),
          panel.border = element_rect(color="grey", fill = NA),
          axis.text.y = element_text(size=14), 
          axis.text.x = element_text(angle = 0, size=14),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          legend.title= element_blank(), 
          legend.background = element_blank(),
          legend.key.size = unit(2, 'lines')
    ) 
  #scale_y_continuous(breaks = seq(0, max(ylimit), 20000)) 
  #scale_x_log10() + 
  #scale_shape_manual(values=c(23,24)) 
  #stat_smooth(method="glm", method.args = list(family="binomial"), se=F)
  



### ============================================================================================
### ============================================================================================
### get gess-v4 data


## ----- user input -----
sel.main <- subset(main.list, Container=="14 mL round bottom tube" & Date != "20170912" & Date != "20180122" & Date != "20180123" & Date != "20180124" & Date != "20180307") ## 1111, 1116 error as chloro- chems were diffused
sel.main.merged <- merge(seed.list, sel.main, by.x="ID", by.y="SeedID1")
head(sel.main.merged)
colnames(sel.main.merged)[19]<-"MID"
lab1 <- "595nm_kk (A)"
lab2 <- "EGFP_sulim (Counts)"
lab3 <- "mRFP-sulim(label) (Counts)"
user.detector <- "victor"
sel.inducers <- "phenol"

ii<-which(sel.main.merged[,"device_name"]=="ph-GESS-v4" & sel.main.merged[,"Inducer_1"]==sel.inducers)
subs.names1 <- "phenol" 
dev.names <- c("ph-GESS-v4")
## ----------------------

#subset(sel.main.merged, select=c(SeedID, Date.x, plasmid_id, device_name, ))
tmp <- merge(measure.list[,1:13], sel.main.merged, by.x="MainID", by.y="MID")

sel.measures <- subset(tmp, Detector==user.detector)
sel.file.names <- names(table(sel.measures[,"Filename"]))
sel.dat <- cbind(sel.measures, od=rep(0, nrow(sel.measures)), gfp=rep(0, nrow(sel.measures)), rfp=rep(0, nrow(sel.measures)))

for(i in 1:length(sel.file.names)){
  tmp <- as.data.frame(read_excel(paste(Victor.ROOT, sel.file.names[i], sep="\\"), sheet = 1))
  row.id <- substr(tmp[,3], 1, 1)
  col.id <- as.numeric(substr(tmp[,3], 2, 3))
  seli <- which(sel.dat[,"Filename"]==sel.file.names[i])
  o <- match(paste(sel.dat[seli,"Row.y"], sel.dat[seli,"Col.y"]), paste(row.id, col.id))
  if(sum(names(tmp)==lab1)){
    od.val <- tmp[,6]
    sel.dat[seli, "od"] <- od.val[o]
  }
  if(sum(names(tmp)==lab2)){
    gfp.val <- tmp[,8]  
    sel.dat[seli, "gfp"] <- gfp.val[o]
  }
  if(sum(names(tmp)==lab3)){
    rfp.val <- tmp[,10]  
    sel.dat[seli, "rfp"] <- rfp.val[o]
  }
  #fl.data <- data.frame(sel.file.names[i], row.id, col.id, od.val, fl.val)  
  cat(i, "/", length(sel.file.names), "\n");flush.console()
}

## selection
sel.dat2 <- subset(sel.dat, (Inducer_1 %in% subs.names1) & Inducer_2=="None" & device_name %in% dev.names)

## ordering
sel.dat2[,"Inducer_1"] <- factor(sel.dat2[,"Inducer_1"], levels=subs.names1)
sel.dat2[,"Inducer_conc_1"] <- as.numeric(sel.dat2[,"Inducer_conc_1"])

## phenol and etc.
newdata <- subset(sel.dat2, select=c(Date, device_name, Inducer_1, Inducer_conc_1, gfp, od))





br.lab <- names(table(newdata[,"Inducer_conc_1"]))
#br.lab <- c(0, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000)
## break tip 
br <- as.numeric(br.lab)
## set zero to min value to avoid errors in regression fitting and visualization
br[1] <- br[1] + 0.01
#br <- c(0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100, 500, 1000)
newdata[newdata[,"Inducer_conc_1"]==0,"Inducer_conc_1"] <- 0.01
newdata2 <- cbind(newdata, gfp.od=newdata[,"gfp"]/newdata[,"od"])
newdata2[,1] <- as.factor(newdata2[,1])
newdata.mlt <- melt(newdata2, id=c("Date", "device_name", "Inducer_1", "Inducer_conc_1"), measure="gfp.od")
#newdata.mlt <- melt(newdata2, id=c("Date", "device_name", "Airation", "Inducer_1", "Inducer_conc_1"), measure="gfp.od")

#range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#newdata.mlt[,6] <- range01(newdata.mlt[,6])
datnew1 <- data.frame(newdata.mlt %>% group_by(device_name, Inducer_conc_1) %>% summarise(mean=mean(value), sd = sd(value)))
#colnames(datnew1) <- c("Name", "Phenol", "mean", "sd")


## ====================== sender -receiver 1:1 from arum
ROOT<-"C:\\mydocs\\2019\\papers\\celltocell-remote\\00 data\\"
setwd(ROOT)
file.name.full<-paste(ROOT, "mms-v4-comparision.xlsx", sep="")
tmpdat <- data.frame(read_xlsx(file.name.full))
dat <- tmpdat[,c(1,3,6)]
#dat[,3] <- range01(dat[,3])

d1 <- newdata.mlt[,c(2,4,6)]
colnames(d1) <- c("Name", "Phenol", "GFPdOD")
d2 <- subset(dat, Name=="SR", select=c(Name, Phenol, GFPdOD)) 

dd <- rbind(d1, d2)
#dd[,3] <- range01(dd[,3])

datnew2 <- data.frame(dd %>% group_by(Name, Phenol) %>% summarise(mean=mean(GFPdOD), sd = sd(GFPdOD)))
datnew <- datnew2[-c(5,10,18),]
#datnew <- datnew2[c(1:9, 20:28),]
#datnew <- datnew[-c(5,10),]
ylimit <- c(-10000, 160000)
datnew[,1] <- c(rep("GESSv4", 8), rep("MMS", 8))


datnew3 <- rbind(datnew[1:8,] %>% mutate(fold=mean/min(mean), foldsd=sd/min(mean)), datnew[9:16,] %>% mutate(fold=mean/min(mean), foldsd=sd/min(mean)))

ggplot(datnew3, aes(x=factor(Phenol), y=mean, fill=Name, width=.7)) +
  geom_bar(stat="identity", position="dodge", color="black") +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), position=position_dodge(.7), width=0.3) +
  #scale_fill_manual(values=c("#9B945F", "#C8102E")) +
  scale_fill_manual(values=c("#c3ff80", "#648c35", "#C8102E")) +
  labs(y="GFP/OD", x=expression(paste("Phenol concentration (", mu, "M)"))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title= element_blank(), 
        legend.background = element_blank(),
        legend.key.size = unit(1.2, 'lines'), 
        legend.position = c(0.2, 0.85),
        legend.text = element_text(size=14)
  ) +
  scale_y_continuous(limits = ylimit) 
  #scale_x_log10() + 
  #scale_shape_manual(values=c(23,24)) 
  #stat_smooth(method="glm", method.args = list(family="binomial"), se=F)



### =========================== ==========================================
### =============================================================================




## ====================== sender receiver optimal ratio 
ROOT<-"C:\\mydocs\\2019\\papers\\celltocell-remote\\00 data\\"
setwd(ROOT)
file.name.full<-paste(ROOT, "mms-ratio.xlsx", sep="")
tmpdat <- data.frame(read_xlsx(file.name.full))
#dat <- tmpdat[,c(1,3,6)]
#dat[,3] <- range01(dat[,3])
## gfp
datnew1 <- data.frame(tmpdat %>% group_by(Name, Phenol) %>% summarise(gfpod=mean(GFP/OD), sd = sd(GFP/OD)))

ggplot(datnew1, aes(x=factor(Name), y=gfpod, fill=factor(Phenol), width=.7)) +
  geom_bar(stat="identity", position="dodge", color="black") + 
  geom_errorbar(aes(ymin=gfpod-1.96*sd, max=gfpod+1.96*sd), position=position_dodge(.7), width=0.3) +
  #scale_fill_manual(values=c("#9B945F", "#C8102E")) +
  scale_fill_manual(values=c("#c3ff80", "#648c35")) +
  labs(y="GFP/OD", x="Sender:Receiver ratio", fill = "Phenol") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title= element_text(size=14), 
        legend.background = element_blank(),
        legend.key.size = unit(1.5, 'lines'),
        legend.text = element_text(size=14)
  ) +
  scale_y_continuous(breaks = seq(0, max(ylimit), 20000)) 



datnew2 <- data.frame(tmpdat %>% group_by(Name, Phenol) %>% summarise(rfpod=mean(RFP/OD), sd = sd(RFP/OD)))

ggplot(datnew2, aes(x=factor(Name), y=rfpod, fill=factor(Phenol), width=.7)) +
  geom_bar(stat="identity", position="dodge", color="black") + 
  geom_errorbar(aes(ymin=rfpod-1.96*sd, max=rfpod+1.96*sd), position=position_dodge(.7), width=0.3) +
  scale_fill_manual(values=c("#9B945F", "#C8102E")) +
  #scale_fill_manual(values=c("#c3ff80", "#648c35")) +
  labs(y="RFP/OD", x="Sender:Receiver ratio", fill = "Phenol") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title= element_text(size=14), 
        legend.background = element_blank(),
        legend.key.size = unit(1.5, 'lines'),
        legend.text = element_text(size=14)
  ) +
  scale_y_continuous(breaks = seq(0, max(ylimit), 20000)) 


datnew3 <- rbind(melt(datnew1, id.vars=c("Phenol", "Name", "sd")), melt(datnew2, id.vars=c("Phenol", "Name", "sd")) )
levels(datnew3$variable) <- c("GFP/OD", "RFP/OD")

dat4 <- datnew3 %>% filter(Phenol!=0)

ggplot(dat4, aes(x=factor(Name), y=value, fill=factor(variable), width=.7)) +
  geom_bar(stat="identity", position="dodge", color="black") + 
  geom_errorbar(aes(ymin=value-1.96*sd, max=value+1.96*sd), position=position_dodge(.7), width=0.3) +
  scale_fill_manual(values=c("#99cc00", "#C8102E")) +
  #scale_fill_manual(values=c("#c3ff80", "#648c35")) +
  labs(y="Fluorescence/OD", x="Detector:Reporter ratio", fill = "") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title= element_text(size=14), 
        legend.background = element_blank(),
        legend.key.size = unit(1.5, 'lines'),
        legend.text = element_text(size=14), 
        legend.position = c(0.15, 0.9)
  ) +
  scale_y_continuous(breaks = seq(0, max(ylimit), 20000)) 



### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
### ====================== tecan read no use in paper, it's replaced into 2019 data next section (25072019-001-haseong-sender-receiver-time.xlsx)

ROOT<-"C:\\mydocs\\2018\\papers\\celltocell-remote\\00 data\\rawdata\\"
setwd(ROOT)

source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")

file.name.full<-paste(ROOT, "20180529 sj11v2 sj1v2 seed portion growth2.xlsx", sep="")
design.file.name.full<-paste(ROOT, "20180529-96well-design.xlsx", sep="")

## ======= read tecan data
##dat<-read.xls(file.name, header=FALSE)
require(readxl)
library(reshape2)
library(ggplot2)

dat <- read.tecan.data(design.file.name.full, file.name.full, 1)


## get mean & sd 
datnew1 <- dat %>% group_by(time.var, phenol, sender) %>% summarise(mean= mean(Label1), sd = sd(Label1))
datnew1 <- data.frame(datnew1, measure=rep("OD", nrow(datnew1)))

datnew2 <- dat %>% group_by(time.var, phenol, sender) %>% summarise(mean= mean(Label2), sd = sd(Label2))
datnew2 <- data.frame(datnew2, measure=rep("Fluorescence", nrow(datnew2)))

datnew <- rbind(datnew1, datnew2)

## get time
time.secs <- as.numeric(names(table(datnew1[,1])))

## data of the interest
dat.sub <- subset(datnew, sender=="1p" & time.var %in% time.secs[seq(1, length(time.secs), 5)], select=c(time.var, phenol, mean, sd, measure))
ggplot(dat.sub, aes(time.var/60, mean, group=phenol)) + 
  facet_grid(measure~., scales = "free") +
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  #geom_line() +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), position=position_dodge(.7), width=0.3) +
  labs(y="", x="Minutes") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)
  ) +
  geom_smooth(method = 'loess', se = FALSE, span=0.3, color="black")  +
  geom_point(aes(shape=phenol), size=3, stroke=1, fill="white") +
  scale_shape_manual(values=c(22,21,24,25)) +
  #scale_y_continuous(breaks = seq(0, 1.5, 0.2)) +
  scale_x_continuous(breaks = seq(0, 1500, 200)) 
#geom_point(colour = "white", aes(shape=phenol), size = 4.5, fill="white") 



## AHL 

dat <- read.tecan.data(design.file.name.full, file.name.full, 2)


## get mean & sd 
datnew1 <- dat %>% group_by(time.var, ahl, receiver) %>% summarise(mean= mean(Label1), sd = sd(Label1))
datnew1 <- data.frame(datnew1, measure=rep("OD", nrow(datnew1)))

datnew2 <- dat %>% group_by(time.var, ahl, receiver) %>% summarise(mean= mean(Label2), sd = sd(Label2))
datnew2 <- data.frame(datnew2, measure=rep("Fluorescence", nrow(datnew2)))

datnew <- rbind(datnew1, datnew2)

## get time
time.secs <- as.numeric(names(table(datnew1[,1])))

## data of the interest
dat.sub <- subset(datnew, receiver=="1p" & time.var %in% time.secs[seq(1, length(time.secs), 5)], select=c(time.var, ahl, mean, sd, measure))
ggplot(dat.sub, aes(time.var/60, mean, group=ahl)) + 
  facet_grid(measure~., scales = "free") +
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  #geom_line() +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), position=position_dodge(.7), width=0.3) +
  labs(y="", x="Minutes") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)
  ) +
  geom_smooth(method = 'loess', se = FALSE, span=0.3, color="black")  +
  geom_point(aes(shape=ahl), size=3, stroke=1, fill="white") +
  scale_shape_manual(values=c(22,21,24,25)) +
  #scale_y_continuous(breaks = seq(0, 1.5, 0.2)) +
  scale_x_continuous(breaks = seq(0, 1500, 200)) 
#geom_point(colour = "white", aes(shape=phenol), size = 4.5, fill="white") 








### ====================== tecan read 25072019-001-haseong-sender-receiver-time.xlsx
### details C:\mydocs\연구노트\새싹\Sender-receiver 개량 20180412.docx

ROOT<-"C:\\mydocs\\2019\\dev\\standardization\\tecan\\"
setwd(ROOT)

source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")

file.name.full<-paste(ROOT, "25072019-001-haseong-sender-receiver-time.xlsx", sep="")
design.file.name.full<-paste(ROOT, "25072019-001-haseong-sender-receiver-time-design.xlsx", sep="")

## ======= read tecan data
##dat<-read.xls(file.name, header=FALSE)
require(readxl)
library(reshape2)
library(ggplot2)

dat <- read.tecan.block.with.design(design.file.name.full, file.name.full, 1)


## get mean & sd 
dat <- dat %>% mutate(gfp.od=Label2/Label1, rfp.od=Label3/Label1)

datnew1 <- dat %>% group_by(time.var, X.phenol, name) %>% summarise(mean= mean(Label1), sd = sd(Label1))
#datnew1 <- data.frame(datnew1, measure=rep("OD", nrow(datnew1)))

datnew2 <- dat %>% group_by(time.var, X.phenol, name) %>% summarise(mean= mean(gfp.od), sd = sd(gfp.od))
#datnew2 <- data.frame(datnew2, measure=rep("GFP", nrow(datnew2)))

datnew3 <- dat %>% group_by(time.var, X.phenol, name) %>% summarise(mean= mean(rfp.od), sd = sd(rfp.od))
#datnew3 <- data.frame(datnew3, measure=rep("RFP", nrow(datnew3)))

library(dplyr)
datnew <- inner_join(datnew1, datnew2, by=c("time.var", "X.phenol", "name"))
datnew <- inner_join(datnew, datnew3, by=c("time.var", "X.phenol", "name"))
datnew <- as.data.frame(datnew)
colnames(datnew) <- c("time", "phenol", "name", "mean.od", "sd.od", "mean.gfp", "sd.gfp", "mean.rfp", "sd.rfp")
#datnew <- datnew %>% 
#  melt(measure.vars=c("mean.od", "mean.gfp", "mean.rfp"))
head(datnew)

## get time
#time.secs <- as.numeric(names(table(datnew1[,1])))

## data of the interest
#dat.sub <- subset(datnew, sender=="1p" & time.var %in% time.secs[seq(1, length(time.secs), 5)], select=c(time.var, phenol, mean, sd, measure))
times <- as.numeric(names(table(datnew$time)))
sel.times <- times[seq(1,length(times), by = 5)][1:15]

## change thr ratio 1:1 or 2:0 or 0:2
subdat.mean <- datnew %>% filter(time %in% sel.times & name=="1:1") %>% select(time, phenol, name, mean.gfp, mean.rfp) %>% melt(measure.vars=c("mean.gfp", "mean.rfp"))
subdat.sd <- datnew %>% filter(time %in% sel.times & name=="1:1") %>% select(time, phenol, name, sd.gfp, sd.rfp) %>% melt(measure.vars=c("sd.gfp", "sd.rfp"))
subdat <- inner_join(subdat.mean, subdat.sd, by=c("time", "phenol", "name")) %>% mutate(phenol2=as.character(phenol)) %>% filter(phenol2!=" 1")
subdat$variable.x <- factor(subdat$variable.x, levels(subdat$variable.x)[c(2,1)])
levels(subdat$variable.x) <- c("RFP/OD", "GFP/OD")

## GFP
ggplot(subdat, aes(x=time/60, y=value.x, color=phenol)) + 
  facet_grid(variable.x~., scales="fixed") +
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  geom_line() +
  geom_errorbar(aes(ymin=value.x-1.96*value.y, max=value.x+1.96*value.y), position=position_dodge(.7), width=0.3) +
  labs(y="", x="Minutes") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.background = element_blank(),
        #strip.text.y = element_blank()
        #strip.background = element_rect(colour="black", fill="white", linetype="solid"),
        strip.text.y = element_text(size=16, color="black"),
        legend.position = c(0.2, 0.85),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14),
        legend.background = element_rect(linetype="solid", colour ="gray")
  ) +
  #geom_smooth(method = 'loess', se = FALSE, span=0.3, color="black")  +
  geom_point(aes(shape=phenol, fill=phenol), size=3, stroke=1, color="black") +
  scale_shape_manual(values=c(22,21,24,25)) +
  scale_y_continuous(breaks = seq(0, 11000, 2000), limits = c(-100,11000)) +
  scale_x_continuous(breaks = seq(0, 1500, 200)) 
  #geom_point(colour = "white", aes(shape=phenol), size = 4.5, fill="white") 







### ====================== freeze conditions

ROOT<-"C:\\mydocs\\2019\\papers\\celltocell-remote\\00 data\\freeze"
setwd(ROOT)

source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")


## ======= read tecan data
##dat<-read.xls(file.name, header=FALSE)
require(readxl)
library(reshape2)
library(ggplot2)

file.name.full1<-paste(ROOT, "/fig4-a.xlsx", sep="")
dat1 <- as.data.frame(read_xlsx(file.name.full1, col_names = F))

file.name.full2<-paste(ROOT, "/fig4-b.xlsx", sep="")
dat2 <- as.data.frame(read_xlsx(file.name.full2, col_names = F))

dat <- rbind(data.frame(dat1, type=rep("Freeze", nrow(dat1))), data.frame(dat2, type=rep("Non-freeze", nrow(dat2))))
colnames(dat) <- c("hours", "phenol", "GFPOD", "condition", "replicate", "type")
## mean, sd
newdat <- dat %>% group_by(hours, phenol, condition, type) %>%  summarise(mean=mean(GFPOD), sd=sd(GFPOD)) %>% 
  filter(hours %in% c(0, 9, 19))

#levels(newdat$type) <- c("freeze", "non-freeze")
newdat$condition <- factor(newdat$condition, levels=c("LB", "LB-MM", "LB-LB-MM"))
ggplot(newdat, aes(x=factor(hours), y=mean, fill=factor(phenol), width=.7)) +
  facet_grid(type~condition, scales="fixed") +
  geom_bar(stat="identity", position="dodge", color="black") + 
  geom_errorbar(aes(ymin=mean-sd, max=mean+sd), position=position_dodge(.7), width=0.3) +
  scale_fill_manual(values=c("white", "#99cc00",  "#648c35")) +
  labs(y="Fluorescence/OD", x="Hours", fill = "phenol") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title= element_text(size=12), 
        #legend.background = element_blank(),
        legend.key.size = unit(1, 'lines'),
        legend.text = element_text(size=11), 
        legend.position = c(0.1, 0.8),
        strip.background = element_rect(colour="black", fill="white", linetype="solid"),
        strip.text = element_text(size=11)
  ) +
  scale_y_continuous(breaks = seq(0, 300000, 50000)) 










### ====================== tecan figure 5
### 원본 파일 필요 

ROOT<-"C:\\mydocs\\2019\\papers\\celltocell-remote\\00 data\\freeze"
setwd(ROOT)

source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")

## ======= read tecan data
require(readxl)
library(reshape2)
library(ggplot2)

file.name.full1<-paste(ROOT, "/fig5-a.xlsx", sep="")
dat <- as.data.frame(read_xlsx(file.name.full1, col_names = F))
colnames(dat) <- c("time", "phenol", "value", "replicate", "ratio", "type")
dat$phenol <- factor(dat$phenol)
dat$replicate <- factor(dat$replicate)
dat$ratio <- factor(dat$ratio)
dat$type <- factor(dat$type)

datnew <- dat %>% group_by(time, phenol, ratio, type) %>% summarise(mean= mean(value), sd = sd(value))
#datnew1 <- data.frame(datnew1, measure=rep("OD", nrow(datnew1)))

datnew1 <- dat %>% filter(type=="OD")
datnew2 <- dat %>% filter(type=="GFP")

datnew3 <- inner_join(datnew1, datnew2, by=c("time", "phenol", "replicate", "ratio")) %>% 
  group_by(time, phenol, ratio) %>% 
  mutate(gfpod=value.y/value.x) %>% 
  summarise(mean=mean(gfpod), sd=sd(gfpod))



library(dplyr)
head(datnew)

## get time
#time.secs <- as.numeric(names(table(datnew1[,1])))

## data of the interest
#dat.sub <- subset(datnew, sender=="1p" & time.var %in% time.secs[seq(1, length(time.secs), 5)], select=c(time.var, phenol, mean, sd, measure))
sel.times <- as.numeric(names(table(datnew$time)))
#sel.times <- times[seq(1,length(times), by = 5)][1:15]

## change thr ratio 1:1 or 2:0 or 0:2
subdat <- datnew %>% filter(type=="OD") 
subdat <- datnew3
#subdat.sd <- datnew %>% filter(time %in% sel.times & name=="1:1") %>% select(time, phenol, name, sd.gfp, sd.rfp) %>% melt(measure.vars=c("sd.gfp", "sd.rfp"))
#subdat <- inner_join(subdat.mean, subdat.sd, by=c("time", "phenol", "name")) %>% mutate(phenol2=as.character(phenol)) %>% filter(phenol2!=" 1")
#subdat$variable.x <- factor(subdat$variable.x, levels(subdat$variable.x)[c(2,1)])
#levels(subdat$variable.x) <- c("RFP/OD", "GFP/OD")

## GFP
ggplot(subdat, aes(x=time, y=mean, color=phenol)) + 
  facet_grid(ratio~., scales="fixed") +
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  geom_line() +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), position=position_dodge(.7), width=0.3) +
  labs(y="", x="Hours") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.background = element_blank(),
        #strip.text.y = element_blank()
        #strip.background = element_rect(colour="black", fill="white", linetype="solid"),
        strip.text.y = element_text(size=16, color="black"),
        legend.position = c(0.8, 0.15),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        #legend.background = element_rect(linetype="solid", colour ="gray")
        legend.background = element_blank()
  ) +
  #geom_smooth(method = 'loess', se = FALSE, span=0.3, color="black")  +
  geom_point(aes(shape=phenol, fill=phenol), size=3, stroke=1, color="black") +
  scale_shape_manual(values=c(22,21,24,25)) +
  scale_y_continuous(breaks = seq(0, 120000, 30000), limits = c(-100,120000)) #+
  #scale_x_continuous(breaks = seq(0, 1500, 200)) 
#geom_point(colour = "white", aes(shape=phenol), size = 4.5, fill="white") 










### ====================== rfp tecan portion with freezed cells

source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")

## ======= read tecan data
##dat<-read.xls(file.name, header=FALSE)
require(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)

ROOT<-"C:\\mydocs\\2018\\papers\\celltocell-remote\\00 data\\rawdata\\"
file.name.full<-paste(ROOT, "20180715-haseong-portion-rfp-sender-receiver-freezed-phenol.xlsx", sep="")
design.file.name.full<-paste(ROOT, "20180715-96well-design.xlsx", sep="")

dat <- read.tecan.data(design.file.name.full, file.name.full, 1)
sr <- paste(dat[,"sender"], dat[,"receiver"], sep=":")
dat <- data.frame(dat, sr)
## get mean & sd 
datnew <- dat %>% group_by(time.var, phenol, sr) %>% summarise(mean= mean(Label3/Label1), sd = sd(Label3/Label1))

## get time
time.secs <- as.numeric(names(table(datnew[,1])))

## data of the interest
dat.sub1 <- subset(datnew, time.var == time.secs[36], select=c(time.var, phenol, mean, sd, sr))
dat.sub1$sr <- ordered(dat.sub1$sr,levels=c("10:0", "8:2", "7:3", "6:4", "5:5", "4:6", "3:7", "2:8", "1:9", "0:10"))



### >>>>>>>>>>>>>>>>>>> 1year later
source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")

ROOT<-"C:\\mydocs\\2019\\dev\\standardization\\tecan\\"
file.name.full<-paste(ROOT, "haseong-freezed-celltocell-20190605.xlsx", sep="")
design.file.name.full<-paste(ROOT, "haseong-freezed-celltocell-20190605-design.xlsx", sep="")

dat <- read.tecan.block.with.design(design.file.name.full, file.name.full, 1)

## get mean & sd 
dat <- dat %>% mutate(gfp.od=Label2/Label1, rfp.od=Label3/Label1, 
                      phenol=factor(gsub("phenol", "", X.phenol)), 
                      detector=factor(gsub("S", "", sender)), 
                      reporter=factor(gsub("R", "", X.receiver)))

datnew1 <- dat %>% group_by(time.var, phenol, detector, reporter) %>% summarise(odmean= mean(Label1), odsd = sd(Label1))
#datnew1 <- data.frame(datnew1, measure=rep("OD", nrow(datnew1)))

datnew2 <- dat %>% group_by(time.var, phenol, detector, reporter) %>% summarise(gfpodmean= mean(gfp.od), gfpodsd = sd(gfp.od))
#datnew2 <- data.frame(datnew2, measure=rep("GFP", nrow(datnew2)))

datnew3 <- dat %>% group_by(time.var, phenol, detector, reporter) %>% summarise(mean= mean(rfp.od), sd = sd(rfp.od))
#datnew3 <- data.frame(datnew3, measure=rep("RFP", nrow(datnew3)))
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

datnew <- inner_join(datnew1, datnew2, by=c("time.var", "phenol", "detector", "reporter")) %>% 
  inner_join(datnew3, by=c("time.var", "phenol", "detector", "reporter"))

time.secs <- as.numeric(names(table(datnew$time.var)))

## data of the interest
dat.sub2 <- datnew %>% filter(time.var == time.secs[36]) %>% mutate(sr=factor(paste(detector, reporter, sep=":"))) %>% select(time.var, phenol, mean, sd, sr)
dat.sub2$sr <- ordered(dat.sub2$sr,levels=c("10:0", "8:2", "7:3", "6:4", "5:5", "4:6", "3:7", "2:8", "1:9", "0:10"))
dat.sub2 <- dat.sub2[,-1]

#dat.sub <- inner_join(dat.sub1, dat.sub2, by=c("phenol", "sr"))
dat.sub <- rbind(as.data.frame(dat.sub1), as.data.frame(dat.sub2))
dat.sub$time.var <- factor(dat.sub$time.var)
levels(dat.sub$time.var) <- c("1Day", "1Year")
#levels(newdat$type) <- c("freeze", "non-freeze")
ggplot(dat.sub, aes(x=sr, y=mean, fill=factor(phenol), width=.7)) +
  facet_grid(time.var~., scales="fixed") +
  geom_bar(stat="identity", position="dodge", color="black") + 
  geom_errorbar(aes(ymin=mean-sd, max=mean+sd), position=position_dodge(.7), width=0.3) +
  scale_fill_manual(values=c("white", "#9B945F", "#C8102E")) +
  labs(y="RFP/OD", x="Detector:Reporter ratio", fill = "phenol") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=12), 
        axis.text.x = element_text(angle = 0, size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title= element_text(size=12), 
        #legend.background = element_blank(),
        legend.key.size = unit(1, 'lines'),
        legend.text = element_text(size=11), 
        legend.position = c(0.1, 0.87),
        strip.background = element_rect(colour="black", fill="white", linetype="solid"),
        strip.text = element_text(size=12)
  ) +
  scale_y_continuous(breaks = seq(0, 15000, 2000)) 



### old version of time series 

ggplot(dat.sub, aes(time.var/3600, mean, group=phenol)) + 
  facet_wrap(~sr, scales = "fixed") +
  #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
  #geom_line() +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd), position=position_dodge(.7), width=0.3) +
  labs(y="", x="Minutes") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text.x = element_text(size=16, face="bold"), 
        legend.text = element_text(size=14)
  ) +
  #geom_smooth(method = 'loess', se = FALSE, span=0.3, color="black")  +
  #geom_point(aes(shape=phenol), size=3, stroke=1, fill="white") +
  scale_shape_manual(values=c(22,21,24)) +
  #scale_y_continuous(breaks = seq(0, 1.5, 0.2)) +
  scale_y_continuous(breaks = seq(0, 15000, 2000)) 
#geom_point(colour = "white", aes(shape=phenol), size = 4.5, fill="white") 


