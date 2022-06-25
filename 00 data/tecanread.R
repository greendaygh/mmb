library(readxl)
source("C:\\mydocs\\lib\\R\\fluorescence_analysis.R")

#### =========================================================================
####  =========== 
#### =========================================================================


file.name.full <- "C:\\mydocs\\2019\\dev\\standardization\\tecan\\20180715-haseong-portion-rfp-sender-receiver-freezed-phenol.xlsx"
#file.name.full<-paste(ROOT, "20180715-haseong-portion-rfp-sender-receiver-freezed-phenol.xlsx", sep="")
design.file.name.full<-"C:\\mydocs\\2018\\dev\\standardization\\04-96well-design-201807152.xlsx"

## ======= read tecan data
##dat<-read.xls(file.name, header=FALSE)
require(readxl)
library(reshape2)
library(ggplot2)
library(plyr)

dat <- read.tecan.data(design.file.name.full, file.name.full)

dat2 <- dat %>% mutate(specfl=Label3/Label1, type=factor(paste("S", sender, "R", receiver, sep=""), levels=c("S10R0", "S8R2", "S7R3", "S6R4", "S5R5", "S4R6", "S3R7", "S2R8", "S1R9", "S0R10"))) %>% 
  select(time=time.var, phenol, rep, type, specfl)  

dat2 <- ddply(dat2, c("time", "type", "phenol"), summarise,
              mean = mean(specfl),
              sd = sd(specfl),
              time = mean(time)
              )

time.secs <- as.numeric(names(table(dat$time.var)))
dat_0y <- subset(dat2, dat2$time %in% time.secs[seq(1, length(time.secs), 4)])

ggplot(datnew, aes(x=time.var/60, y=mean, col=phenol)) + 
  geom_point() +
  geom_errorbar(aes(ymin=mean-1.96*sd, max=mean+1.96*sd)) +
  facet_wrap(~type, nrow=2) +
  labs(y="Specific fluorescence", x="Minutes") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_blank(),
        panel.border = element_rect(color="grey", fill = NA),
        axis.text.y = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text.x = element_text(size = 14)
  ) +
  geom_smooth(method = 'loess', se = FALSE, span=0.2)  +
  geom_point(aes(shape=phenol), size=4.5, stroke=2, fill="white") +
  scale_shape_manual(values=c(22,23,24)) 






#### =========================================================================
####  =========== after one year
#### =========================================================================


data_file_name2 <- "C:\\mydocs\\2019\\dev\\standardization\\tecan\\haseong-freezed-celltocell-20190605.xlsx"

mydesign <- as.data.frame(read_excel(data_file_name2, sheet="EXPDESIGN", range="A1:L8", skip = 0, col_names=F))
mydata <- as.data.frame(read_excel(data_file_name2, sheet=1, col_names=F))
#head(mydesign)

#tmp <- read.tecan.data.columsample(data_file_name)

pos1 <- rep(LETTERS[1:8], time=12)
pos2 <- rep(sprintf("%02d", 1:12), each=8)
well_position_labels <- paste(pos1, pos2, sep="")
well_position_matrix <- matrix(well_position_labels, nrow=8, ncol=12)

sel_well_idx <- which(!is.na(mydesign))
sel_well_pos_names <- well_position_matrix[sel_well_idx]
sel_well_sample_names <- as.matrix(mydesign)[sel_well_idx]

label.idx<-grep(pattern="^Label[0-9]$", as.character(mydata[,1]))
num.labels<-length(label.idx)
label.types<-as.character(mydata[label.idx,1])

strt.idx <- which(as.character(mydata[,1])=="0s")
ed.idx <- c(c(strt.idx-1)[-1], (strt.idx[2]-strt.idx[1])*length(strt.idx)+strt.idx[1]-1)
well.names <- as.character(unlist(mydata[strt.idx[1]+1,-1]))

## label1
tmpdat <- mydata[strt.idx[1]:ed.idx[1],]
time.idx <- grep(".+s$", tmpdat[,1])
time.sec<-as.numeric(gsub(pattern="s", replacement = "", tmpdat[time.idx,1]))
num.times <- length(time.sec)

lab1mat <- matrix(0, nrow=num.times, ncol=length(sel_well_idx))
colnames(lab1mat) <- sel_well_sample_names
rownames(lab1mat) <- as.character(time.sec)
for(i in 1:num.times){
  dat_start <- time.idx[i]+2
  dat_end <- time.idx[i]+9
  valmat <- as.matrix(tmpdat[dat_start:dat_end,2:13])
  lab1mat[i,] <- valmat[sel_well_idx]
}


## label2
tmpdat <- mydata[strt.idx[2]:ed.idx[2],]
time.idx <- grep(".+s$", tmpdat[,1])
time.sec<-as.numeric(gsub(pattern="s", replacement = "", tmpdat[time.idx,1]))
num.times <- length(time.sec)

lab2mat <- matrix(0, nrow=num.times, ncol=length(sel_well_idx))
colnames(lab2mat) <- sel_well_sample_names
rownames(lab2mat) <- as.character(time.sec)
for(i in 1:num.times){
  dat_start <- time.idx[i]+2
  dat_end <- time.idx[i]+9
  valmat <- as.matrix(tmpdat[dat_start:dat_end,2:13])
  lab2mat[i,] <- valmat[sel_well_idx]
}


## label3
tmpdat <- mydata[strt.idx[3]:ed.idx[3],]
time.idx <- grep(".+s$", tmpdat[,1])
time.sec<-as.numeric(gsub(pattern="s", replacement = "", tmpdat[time.idx,1]))
num.times <- length(time.sec)

lab3mat <- matrix(0, nrow=num.times, ncol=length(sel_well_idx))
colnames(lab3mat) <- sel_well_sample_names
rownames(lab3mat) <- as.character(time.sec)
for(i in 1:num.times){
  dat_start <- time.idx[i]+2
  dat_end <- time.idx[i]+9
  valmat <- as.matrix(tmpdat[dat_start:dat_end,2:13])
  lab3mat[i,] <- valmat[sel_well_idx]
}

library(reshape2)

odmlt <- melt(lab1mat, value.name="od")
gfpmlt <- melt(lab2mat, value.name="gfp")
rfpmlt <- melt(lab3mat, value.name="rfp")

library(dplyr)

join_dat <- inner_join(odmlt, gfpmlt, by=c("Var1", "Var2"))
join_dat <- inner_join(join_dat, rfpmlt, by=c("Var1", "Var2"))


library(ggplot2)

ggplot(join_dat, aes(x=Var1, y=rfp, color=Var2)) +
  geom_point()


join_dat2 <- join_dat %>% filter(!is.na(Var2)) %>% mutate_if(is.factor, as.character)
condmat <- do.call(rbind.data.frame, strsplit(join_dat2$Var2, split=";"))
colnames(condmat) <- c("Var3", "Var4", "Rep", "Var5", "Var6")
join_dat3 <- cbind(join_dat2, condmat) %>% mutate(Var7 = paste(Var3, Var4, sep=""))

dat_1y <- join_dat3 %>% select(type=Var7, phenol=Var6, rfp=rfp, od=od, time=Var1, rep=Rep) %>% 
          mutate(specfl=rfp/od) %>% 
          ddply(c("type", "phenol", "time"), summarise, mean=mean(specfl), sd=sd(specfl))
          #group_by(time, phenol, type) %>% 

ggplot(join_dat3, aes(x=Var1, y=rfp/od, color=Var6)) +
  facet_grid(cols=vars(Var7)) +
  geom_point() +
  geom_smooth()


#### =========================================================================
####  =========== plot
#### =========================================================================

#dat_0y
#dat_1y

head(dat_0y)
head(dat_1y)

ggplot(dat_0y, aes(x=time, y=mean, color=phenol)) +
  facet_grid(cols=vars(type)) +
  geom_point() +
  geom_smooth()

dat_1y$type <- factor(dat_1y$type, levels=c("S10R0", "S8R2", "S7R3", "S6R4", "S5R5", "S4R6", "S3R7", "S2R8", "S1R9", "S0R10"))
ggplot(dat_0y, aes(x=time, y=mean, color=phenol)) +
  facet_grid(cols=vars(type)) +
  geom_point() +
  geom_smooth() +
  ylim(c(0, 15000))

ggplot(dat_1y, aes(x=time, y=mean, color=phenol)) +
  facet_grid(cols=vars(type)) +
  geom_point() +
  geom_smooth() +
  ylim(c(0, 15000))


# dat_0y





