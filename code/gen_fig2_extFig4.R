# This R code produces Fig. 2 and Extended Data Fig. 4 in the manuscript Gilbert, Huang et al. 
# For AMP, it analyzes least sensitive viruses, and analyzes the 703/081 and 704/085 trials pooled.

# The user needs to install the package sievePH available at CRAN
library(sievePH)

# Set the working directory:
# If Linux:
# setwd("/trials/vaccine/p704/analysis/manuscripts/NeutTiterBiomarker") 
# If Windows:
# setwd("T:/vaccine/p704/analysis/manuscripts/NeutTiterBiomarker")
# assuming this script is run from the 'code' directory
setwd("..")
  
# Define the inputs needed, for 
# 703+704 least sensitive viruses:

# Fig. 2
figureic80 <- "./figures/Fig2ID80_leastsensitiveBAMA_703704.pdf"
outpreic80 <- './data/outpreic80poolleastsens703704_VE'
inputfileic80pool <- './data/inputdataic80poolleastsens703704.dat'
# Set A data (all N=274 animals from Pegu et al. 2019 Cell Host Microbe):
datNHP80A <- read.csv("./data/NHPID80protectionresults_outliersremoved.csv")
# Set B data: Repeat restricting to CD4bs epitope bnAbs and removing SF162P3, and TZM-bl assay only:
datNHP80B <- read.csv("./data/NHPID80protectionresultsCD4bsonlyminusSF162P3.csv")
# Set C data: Repeat restricting to all epitopes excluding MPER and removing SF162P3, and TZM-bl assay only:
datNHP80C <- read.csv("./data/NHPID80protectionresultsallminusMPERandSF162P3.csv")

# Extended Data Fig. 4.
figureic50 <- "./figures/ExtFig4ID50_leastsensitiveBAMA_703704.pdf"
outpreic50 <- './data/outpreic50poolleastsens703704_VE'
inputfileic50pool <- './data/inputdataic50poolleastsens703704.dat'
datNHP50A <- read.csv("./data/NHPID50protectionresults.csv")
datNHP50B <- read.csv("./data/NHPID50protectionresultsCD4bsonlyminusSF162P3.csv")
datNHP50C <- read.csv("./data/NHPID50protectionresultsallminusMPERandSF162P3.csv")

#titleic80 <- "PE by ID80 of Least Sens. Var.: 703+704"
#titleic90 <- "PE by IC50 of Least Sens. Var.: 703+704"

outfilebandconc <- './data/bandwidthconc703704.dat'

# Completed file processing
######################################################

bandconc <- scan(outfilebandconc) 
hbandic50 <- round(bandconc[2],2)
hbandic80 <- round(bandconc[3],2)
medianmidpointconcpool704 <- bandconc[7] 
medianmidpointconc10704 <- bandconc[5]
medianmidpointconc30704 <- bandconc[6]
  
mat <- matrix(scan(inputfileic50pool),ncol=9,byrow=T)
trt <- mat[,4]
V <- mat[,8]  # IC50 most sensitive founder variant
keep <- mat[,7]==1
grp <- trt+1
#Original mark scale
VOrig <- mat[,9]

mat <- matrix(scan(inputfileic80pool),ncol=9,byrow=T)
V2 <- mat[,8] # IC80 most sensitive founder variant
keep2 <- mat[,7]==1
V2Orig <- mat[,9]

minic50 <- min(VOrig[keep])
maxminusminic50 <- max(VOrig[keep] - minic50)
minic80 <- min(V2Orig[keep2])
maxminusminic80 <- max(V2Orig[keep2]-minic80)

# getlab puts the plotting tick marks at the right place

getlab <- function(V,VOrig,keep,ic=TRUE) {
  xx <- sort(V[keep])
  yy <- sort(10^VOrig[keep])
  if (!ic) {
    yy <- sort(VOrig[keep]) }
  xxextend <- xx[1]
  yyextend <- yy[1]
  for (i in 1:(length(xx)-1)) {
    xxextend <- c(xxextend,seq(xx[i],xx[i+1],len=200))
    yyextend <- c(yyextend,seq(yy[i],yy[i+1],len=200))
  }
  
  
  labelans <- c(yy[1],0.3,1,3,5,max(yy)) 
  # later change to labelans <- c(yy[1],1,2,3,5,10,20,50,max(yy))
  if(!ic) { 
    labelans <- c(yy[1],1,2,3,4,max(yy)) }
  m <- length(labelans)
  atans <- rep(NA,m)
  # may be of interest to change to labelans <- c(yy[1],0.1,0.2,0.4,0.6,1,2,5,10,20,50,100)
  atans[1] <- xx[1]
  atans[m] <- xx[length(xx)]
  indices <- c(1:length(xxextend))
  
  for (i in 2:(m-1)) {
    ind <- min(indices[yyextend>=labelans[i]])
    atans[i] <- xxextend[ind] }
  
  return(list(atans,labelans))
}

ans <- getlab(V,VOrig,keep)
atansic50 <- ans[[1]]
labelansic50 <- round(ans[[2]],2)
ans <- getlab(V2,V2Orig,keep2)
atansic80 <- ans[[1]]
labelansic80 <- round(ans[[2]],2)

VACAMP_PE_A0<-matrix(scan(outpreic80,skip=2),ncol=9,byrow=T)

vAMP<-VACAMP_PE_A0[,1]
logRR <-VACAMP_PE_A0[,2]
logRRSE<-VACAMP_PE_A0[,3]
# 95% Confidence Interval of PE(v)
PE <- 1-exp(logRR)
PElow<-1-exp(logRR+1.96*logRRSE)
PEup<-1-exp(logRR-1.96*logRRSE)
mat <- cbind(vAMP,PE,PElow,PEup)
colnames(mat) <- c("mark","TE","LB","UB")
d <- data.frame(mat)

# Properties of the data frame:
# a univariate continuous mark is considered
# 4 columns with variables in this order
# column names don't need to match those above
# first column are mark values (e.g., on the log10 scale) on a plotting grid 
# (a grid of 200 values is sufficient)
# each row contains inference for a single value on the grid

# save 'd' as a list, name its first component, and change its class
d <- list(pe=d)
class(d) <- "summary.sievePH"

source("./code/plot.summary.sievePH.R")

# 3. use the 'plot' method for objects of class 'summary.sievePH' ---------

# with scatter and box plots of marks
# suppose that 'mark' and 'tx' are numeric vectors of mark values and treatment indicators
# for everyone in the analysis cohort (mark=NA for right-censored participants);
# provide directory path and file name in the 'pdf' function

vAMPOrig <- 10**(vAMP*maxminusminic80 + minic80)
marknew <- V2
marknew[V2Orig==99] <- NA

# Version of the PE-by-IC80 plot in terms of Predicted ID80 at median mid-point concentration.
# The only difference is the x-axis uses labels on the median mid-point concentration/IC80 scale:

mat <- cbind(vAMP,rev(PE),rev(PElow),rev(PEup))
colnames(mat) <- c("mark","TE","LB","UB")
d <- data.frame(mat)
d <- list(pe=d)
class(d) <- "summary.sievePH"

#BAMA 
conc10 <- 0.5417*9.3 + (1-0.5417)*12.2
conc30 <- 0.6316*28.4 + (1-0.6316)*32.1
cpoolBAMA <- (60/(60+47))*conc10 + (1-(60/(60+47)))*conc30

vAMPOrig <- 10**(vAMP*maxminusminic80 + minic80)

# BAMA:
vAMPOrig <- cpoolBAMA/vAMPOrig

# Add to the plot the main result of Pegu et al. 2019 (Protection by ID80)

# Use linear interpolation on the log10 scale.
# vAMP goes from 0.02 to 1.0 in increments of 0.01.
# vAMPOrig is ID80 for AMP defined at each of the 99 points

# vAMPscaledNHPid80titer <- rep(NA,length(datNHP$id80titer))
# Set the lowest ID80's of the NHP grid to equal the lowest value in the AMP grid:
log10NHPid80A <- datNHP80A$id80titer
log10NHPid80B <- datNHP80B$id80titer
log10NHPid80C <- datNHP80C$id80titer

# Only include the NHP titers that are >= the smallest AMP one and <= the largest AMP one
kpA <- log10NHPid80A >= min(log10(vAMPOrig)) & log10NHPid80A <= max(log10(vAMPOrig))
log10NHPid80A <- log10NHPid80A[kpA]
kpB <- log10NHPid80B >= min(log10(vAMPOrig)) & log10NHPid80B <= max(log10(vAMPOrig))
log10NHPid80B <- log10NHPid80B[kpB]
kpC <- log10NHPid80C >= min(log10(vAMPOrig)) & log10NHPid80C <= max(log10(vAMPOrig))
log10NHPid80C <- log10NHPid80C[kpC]


# Linearly interpolate the in-range values:
# Reverse the order of log10(ID80) in AMP for easier calculations:
log10AMPid80 <- sort(log10(vAMPOrig))
xaxisvalueNHPA <- rep(NA,length(kpA[kpA]))
xaxisvalueNHPB <- rep(NA,length(kpB[kpB]))
xaxisvalueNHPC <- rep(NA,length(kpC[kpC]))

for (i in 1:length(kpA[kpA])) {
# Locate the 2 points to interpolate between
lowerpt <- max(log10AMPid80[log10AMPid80 <= log10NHPid80A[i]])
upperpt <- min(log10AMPid80[log10AMPid80 >= log10NHPid80A[i]])
indlowerpt <- c(1:99)[log10AMPid80==lowerpt]
indupperpt <- c(1:99)[log10AMPid80==upperpt]
percway <- (log10NHPid80A[i] - lowerpt)/(upperpt - lowerpt)
xaxisvalueNHPA[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}

for (i in 1:length(kpB[kpB])) {
lowerpt <- max(log10AMPid80[log10AMPid80 <= log10NHPid80B[i]])
upperpt <- min(log10AMPid80[log10AMPid80 >= log10NHPid80B[i]])
indlowerpt <- c(1:99)[log10AMPid80==lowerpt]
indupperpt <- c(1:99)[log10AMPid80==upperpt]
percway <- (log10NHPid80B[i] - lowerpt)/(upperpt - lowerpt)
xaxisvalueNHPB[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}

for (i in 1:length(kpC[kpC])) {
lowerpt <- max(log10AMPid80[log10AMPid80 <= log10NHPid80C[i]])
upperpt <- min(log10AMPid80[log10AMPid80 >= log10NHPid80C[i]])
indlowerpt <- c(1:99)[log10AMPid80==lowerpt]
indupperpt <- c(1:99)[log10AMPid80==upperpt]
percway <- (log10NHPid80C[i] - lowerpt)/(upperpt - lowerpt)
xaxisvalueNHPC[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}

#}

# Scale ID80 in a similar way
#BAMA:
#xtickAtquant <- c(2,5,10,20,30,60,125,250)
# Change to a log2 scale (antilogged)
#xtickAtquant <- c(2^c(1:8))
xtickAtquant <- c(2^c(1:7),234)

newvAMP <- rep(NA,length(xtickAtquant))
revvAMPOrig <- sort(vAMPOrig)
for (i in 1:length(xtickAtquant)) {
# Locate the 2 points to interpolate between
lowerpt <- max(revvAMPOrig[revvAMPOrig <= xtickAtquant[i]])
upperpt <- min(revvAMPOrig[revvAMPOrig >= xtickAtquant[i]])
indlowerpt <- c(1:99)[revvAMPOrig==lowerpt]
indupperpt <- c(1:99)[revvAMPOrig==upperpt]
percway <- (xtickAtquant[i] - lowerpt)/(upperpt - lowerpt)
newvAMP[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}


prot <- rev(PE)

upperpt <- min(rev(log10(vAMPOrig))[prot >= 0.50])
lowerpt <- max(rev(log10(vAMPOrig))[prot <= 0.50])
percway <- (0.50 - max(prot[prot <= 0.50]))/(min(prot[prot >= 0.50]) - max(prot[prot <= 0.50]))
AMPid8050 <- 10^(lowerpt + percway*(upperpt-lowerpt))
AMPind50 <- min(c(1:length(vAMP))[prot >= 0.50])

upperpt <- min(rev(log10(vAMPOrig))[prot >= 0.75])
lowerpt <- max(rev(log10(vAMPOrig))[prot <= 0.75])
percway <- (0.75 - max(prot[prot <= 0.75]))/(min(prot[prot >= 0.75]) - max(prot[prot <= 0.75]))
AMPid8075 <- 10^(lowerpt + percway*(upperpt-lowerpt))
AMPind75 <- min(c(1:length(vAMP))[prot >= 0.75])

upperpt <- min(rev(log10(vAMPOrig))[prot >= 0.90])
lowerpt <- max(rev(log10(vAMPOrig))[prot <= 0.90])
percway <- (0.90 - max(prot[prot <= 0.90]))/(min(prot[prot >= 0.90]) - max(prot[prot <= 0.90]))
AMPid8090 <- 10^(lowerpt + percway*(upperpt-lowerpt))
AMPind90 <- min(c(1:length(vAMP))[prot >= 0.90])

# Compute values of NHP id80 titer at 50%, 75%, 90% PE:
protNHP80A <- datNHP80A$protectlevel[kpA]

upperpt <- min(log10NHPid80A[protNHP80A >= 50])
lowerpt <- max(log10NHPid80A[protNHP80A <= 50])
percway <- (50 - max(protNHP80A[protNHP80A <= 50]))/
           (min(protNHP80A[protNHP80A >= 50]) - max(protNHP80A[protNHP80A <= 50]))
NHPid80A50 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80A50 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80A50])

upperpt <- min(log10NHPid80A[protNHP80A >= 75])
lowerpt <- max(log10NHPid80A[protNHP80A <= 75])
percway <- (75 - max(protNHP80A[protNHP80A <= 75]))/
           (min(protNHP80A[protNHP80A >= 75]) - max(protNHP80A[protNHP80A <= 75]))
NHPid80A75 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80A75 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80A75])

upperpt <- min(log10NHPid80A[protNHP80A >= 90])
lowerpt <- max(log10NHPid80A[protNHP80A <= 90])
percway <- (90 - max(protNHP80A[protNHP80A <= 90]))/
           (min(protNHP80A[protNHP80A >= 90]) - max(protNHP80A[protNHP80A <= 90]))
NHPid80A90 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80A90 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80A90])


protNHP80B <- datNHP80B$protectlevel[kpB]

upperpt <- min(log10NHPid80B[protNHP80B >= 50])
lowerpt <- max(log10NHPid80B[protNHP80B <= 50])
percway <- (50 - max(protNHP80B[protNHP80B <= 50]))/
           (min(protNHP80B[protNHP80B >= 50]) - max(protNHP80B[protNHP80B <= 50]))
NHPid80B50 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80B50 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80B50])

upperpt <- min(log10NHPid80B[protNHP80B >= 75])
lowerpt <- max(log10NHPid80B[protNHP80B <= 75])
percway <- (75 - max(protNHP80B[protNHP80B <= 75]))/
           (min(protNHP80B[protNHP80B >= 75]) - max(protNHP80B[protNHP80B <= 75]))
NHPid80B75 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80B75 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80B75])

upperpt <- min(log10NHPid80B[protNHP80B >= 90])
lowerpt <- max(log10NHPid80B[protNHP80B <= 90])
percway <- (90 - max(protNHP80B[protNHP80B <= 90]))/
           (min(protNHP80B[protNHP80B >= 90]) - max(protNHP80B[protNHP80B <= 90]))
NHPid80B90 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80B90 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80B90])


protNHP80C <- datNHP80C$protectlevel[kpC]

upperpt <- min(log10NHPid80C[protNHP80C >= 50])
lowerpt <- max(log10NHPid80C[protNHP80C <= 50])
percway <- (50 - max(protNHP80C[protNHP80C <= 50]))/
           (min(protNHP80C[protNHP80C >= 50]) - max(protNHP80C[protNHP80C <= 50]))
NHPid80C50 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80C50 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80C50])

upperpt <- min(log10NHPid80C[protNHP80C >= 75])
lowerpt <- max(log10NHPid80C[protNHP80C <= 75])
percway <- (75 - max(protNHP80C[protNHP80C <= 75]))/
           (min(protNHP80C[protNHP80C >= 75]) - max(protNHP80C[protNHP80C <= 75]))
NHPid80C75 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80C75 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80C75])

upperpt <- min(log10NHPid80C[protNHP80C >= 90])
lowerpt <- max(log10NHPid80C[protNHP80C <= 90])
percway <- (90 - max(protNHP80C[protNHP80C <= 90]))/
           (min(protNHP80C[protNHP80C >= 90]) - max(protNHP80C[protNHP80C <= 90]))
NHPid80C90 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind80C90 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid80C90])


### The figure shows 3 lines for NHP for set A, B, C: 

protectlevelA <- datNHP80A$protectlevel
protectlevelB <- datNHP80B$protectlevel



protectlevelC <- datNHP80C$protectlevel

pdf(figureic80,width=0.95*7, height=0.95*6)
xlabs <- quantile(vAMPOrig,prob=seq(0,1,len=length(xtickAtquant)))
xlabs <- c(round(xlabs[1:4],1),round(xlabs[5:8],0))
xlabs <- xtickAtquant
plot(d, 
     mark=NULL, 
     tx=NULL, 
     ylim=c(-0.4, 1),
     xtickAt=quantile(vAMP,prob=seq(0,1,len=8)),
     xtickLab=xlabs,
     ytickAt=seq(-0.4, 1, by=0.20),
     ytickLab=seq(-40, 100, by=20),
     xlab=expression(paste("Predicted ", ID[80]," against autologous virus")),
     ylab="Prevention Efficacy (%)",
     txLab=c("Placebo", "VRC01"))
lines(xaxisvalueNHPA,protectlevelA[kpA]/100,lwd=6,lty=1,col="blue")
lines(xaxisvalueNHPB,protectlevelB[kpB]/100,lwd=6,lty=1,col="orange")
lines(xaxisvalueNHPC,protectlevelC[kpC]/100,lwd=6,lty=1,col="green")
abline(h=0.5,lty=3,lwd=1)
text(0.1,0.52,"50% Efficacy")
text(0.1,0.77,"75% Efficacy")
text(0.1,0.92,"90% Efficacy")

abline(h=0.75,lty=3,lwd=1)
abline(h=0.90,lty=3,lwd=1)
#mtext(as.character(paste("bandwidth",hbandic80)), side=1, line=-6, adj=1, cex=1.2)
#mtext("BAMA for AMP, ELISA for NHP study", side=1, line=-3, adj=1, cex=1.1)
legend(0.6,-0.05,legend=c("AMP","NHP Set A", "NHP Set B", "NHP Set C"),col=c("black","blue","orange","green"),
lty=c(1,1,1,1),lwd=4,cex=1.1)
#title(title8)
dev.off()


# Repeat for IC50:
VACAMP_PE_A0<-matrix(scan(outpreic50,skip=2),ncol=9,byrow=T)

vAMP<-VACAMP_PE_A0[,1]
logRR <-VACAMP_PE_A0[,2]
logRRSE<-VACAMP_PE_A0[,3]
# 95% Confidence Interval of PE(v)
PE <- 1-exp(logRR)
PElow<-1-exp(logRR+1.96*logRRSE)
PEup<-1-exp(logRR-1.96*logRRSE)
mat <- cbind(vAMP,PE,PElow,PEup)
colnames(mat) <- c("mark","TE","LB","UB")
d <- data.frame(mat)
d <- list(pe=d)
class(d) <- "summary.sievePH"

vAMPOrig <- 10**(vAMP*maxminusminic50 + minic50)
marknew <- V
marknew[VOrig==99] <- NA

mat <- cbind(vAMP,rev(PE),rev(PElow),rev(PEup))
colnames(mat) <- c("mark","TE","LB","UB")
d <- data.frame(mat)
d <- list(pe=d)
class(d) <- "summary.sievePH"

#BAMA 
conc10 <- 0.5417*9.3 + (1-0.5417)*12.2
conc30 <- 0.6316*28.4 + (1-0.6316)*32.1
cpoolBAMA <- (60/(60+47))*conc10 + (1-(60/(60+47)))*conc30

# IC50 scale
vAMPOrig <- 10**(vAMP*maxminusminic50 + minic50)
# ID50 scale

# BAMA
vAMPOrig <- cpoolBAMA/vAMPOrig

# Add to the plot the main result of Pegu et al. 2019 (Protection by ID50)

# Use linear interpolation on the log10 scale.
# vAMP goes from 0.02 to 1.0 in increments of 0.01.
# vAMPOrig is ID50 for AMP defined at each of the 99 points

# vAMPscaledNHPid50titer <- rep(NA,length(datNHP$id50titer))
# Set the lowest ID50's of the NHP grid to equal the lowest value in the AMP grid:
log10NHPid50A <- datNHP50A$id50titer
log10NHPid50B <- datNHP50B$id50titer
log10NHPid50C <- datNHP50C$id50titer

# Only include the NHP titers that are >= the smallest AMP one and <= the largest AMP one
kpA <- log10NHPid50A >= min(log10(vAMPOrig)) & log10NHPid50A <= max(log10(vAMPOrig))
log10NHPid50A <- log10NHPid50A[kpA]
kpB <- log10NHPid50B >= min(log10(vAMPOrig)) & log10NHPid50B <= max(log10(vAMPOrig))
log10NHPid50B <- log10NHPid50B[kpB]
kpC <- log10NHPid50C >= min(log10(vAMPOrig)) & log10NHPid50C <= max(log10(vAMPOrig))
log10NHPid50C <- log10NHPid50C[kpC]


# Linearly interpolate the in-range values:
# Reverse the order of log10(ID50) in AMP for easier calculations:
log10AMPid50 <- sort(log10(vAMPOrig))
xaxisvalueNHPA <- rep(NA,length(kpA[kpA]))
xaxisvalueNHPB <- rep(NA,length(kpB[kpB]))
xaxisvalueNHPC <- rep(NA,length(kpC[kpC]))


for (i in 1:length(kpA[kpA])) {
# Locate the 2 points to interpolate between
lowerpt <- max(log10AMPid50[log10AMPid50 <= log10NHPid50A[i]])
upperpt <- min(log10AMPid50[log10AMPid50 >= log10NHPid50A[i]])
indlowerpt <- c(1:99)[log10AMPid50==lowerpt]
indupperpt <- c(1:99)[log10AMPid50==upperpt]
percway <- (log10NHPid50A[i] - lowerpt)/(upperpt - lowerpt)
xaxisvalueNHPA[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}

for (i in 1:length(kpB[kpB])) {
lowerpt <- max(log10AMPid50[log10AMPid50 <= log10NHPid50B[i]])
upperpt <- min(log10AMPid50[log10AMPid50 >= log10NHPid50B[i]])
indlowerpt <- c(1:99)[log10AMPid50==lowerpt]
indupperpt <- c(1:99)[log10AMPid50==upperpt]
percway <- (log10NHPid50B[i] - lowerpt)/(upperpt - lowerpt)
xaxisvalueNHPB[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}

for (i in 1:length(kpC[kpC])) {
lowerpt <- max(log10AMPid50[log10AMPid50 <= log10NHPid50C[i]])
upperpt <- min(log10AMPid50[log10AMPid50 >= log10NHPid50C[i]])
indlowerpt <- c(1:99)[log10AMPid50==lowerpt]
indupperpt <- c(1:99)[log10AMPid50==upperpt]
percway <- (log10NHPid50C[i] - lowerpt)/(upperpt - lowerpt)
xaxisvalueNHPC[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}

#}


# Scale ID50 in a similar way
#BAMA:
#xtickAtquant <- c(2,5,10,25,60,150,400,800)
# Change to log2 scale (antilogged) except 800 at the top:
xtickAtquant <- c(2^c(1:9),800)

newvAMP <- rep(NA,length(xtickAtquant))
revvAMPOrig <- sort(vAMPOrig)

for (i in 1:length(xtickAtquant)) {
# Locate the 2 points to interpolate between
lowerpt <- max(revvAMPOrig[revvAMPOrig <= xtickAtquant[i]])
upperpt <- min(revvAMPOrig[revvAMPOrig >= xtickAtquant[i]])
indlowerpt <- c(1:99)[revvAMPOrig==lowerpt]
indupperpt <- c(1:99)[revvAMPOrig==upperpt]
percway <- (xtickAtquant[i] - lowerpt)/(upperpt - lowerpt)
newvAMP[i] <- vAMP[indlowerpt] + percway*(vAMP[indupperpt] - vAMP[indlowerpt])
}


prot <- rev(PE)

upperpt <- min(rev(log10(vAMPOrig))[prot >= 0.50])
lowerpt <- max(rev(log10(vAMPOrig))[prot <= 0.50])
percway <- (0.50 - max(prot[prot <= 0.50]))/(min(prot[prot >= 0.50]) - max(prot[prot <= 0.50]))
AMPid5050 <- 10^(lowerpt + percway*(upperpt-lowerpt))
AMPind50 <- min(c(1:length(vAMP))[prot >= 0.50])

upperpt <- min(rev(log10(vAMPOrig))[prot >= 0.75])
lowerpt <- max(rev(log10(vAMPOrig))[prot <= 0.75])
percway <- (0.75 - max(prot[prot <= 0.75]))/(min(prot[prot >= 0.75]) - max(prot[prot <= 0.75]))
AMPid5075 <- 10^(lowerpt + percway*(upperpt-lowerpt))
AMPind75 <- min(c(1:length(vAMP))[prot >= 0.75])

upperpt <- min(rev(log10(vAMPOrig))[prot >= 0.90])
lowerpt <- max(rev(log10(vAMPOrig))[prot <= 0.90])
percway <- (0.90 - max(prot[prot <= 0.90]))/(min(prot[prot >= 0.90]) - max(prot[prot <= 0.90]))
AMPid5090 <- 10^(lowerpt + percway*(upperpt-lowerpt))
AMPind90 <- min(c(1:length(vAMP))[prot >= 0.90])

# Compute values of NHP id50 titer at 50%, 75%, 90% PE:
protNHP50A <- datNHP50A$protectlevel[kpA]

upperpt <- min(log10NHPid50A[protNHP50A >= 50])
lowerpt <- max(log10NHPid50A[protNHP50A <= 50])
percway <- (50 - max(protNHP50A[protNHP50A <= 50]))/
           (min(protNHP50A[protNHP50A >= 50]) - max(protNHP50A[protNHP50A <= 50]))
NHPid50A50 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50A50 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50A50])

upperpt <- min(log10NHPid50A[protNHP50A >= 75])
lowerpt <- max(log10NHPid50A[protNHP50A <= 75])
percway <- (75 - max(protNHP50A[protNHP50A <= 75]))/
           (min(protNHP50A[protNHP50A >= 75]) - max(protNHP50A[protNHP50A <= 75]))
NHPid50A75 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50A75 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50A75])

upperpt <- min(log10NHPid50A[protNHP50A >= 90])
lowerpt <- max(log10NHPid50A[protNHP50A <= 90])
percway <- (90 - max(protNHP50A[protNHP50A <= 90]))/
           (min(protNHP50A[protNHP50A >= 90]) - max(protNHP50A[protNHP50A <= 90]))
NHPid50A90 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50A90 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50A90])


protNHP50B <- datNHP50B$protectlevel[kpB]

upperpt <- min(log10NHPid50B[protNHP50B >= 50])
lowerpt <- max(log10NHPid50B[protNHP50B <= 50])
percway <- (50 - max(protNHP50B[protNHP50B <= 50]))/
           (min(protNHP50B[protNHP50B >= 50]) - max(protNHP50B[protNHP50B <= 50]))
NHPid50B50 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50B50 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50B50])

upperpt <- min(log10NHPid50B[protNHP50B >= 75])
lowerpt <- max(log10NHPid50B[protNHP50B <= 75])
percway <- (75 - max(protNHP50B[protNHP50B <= 75]))/
           (min(protNHP50B[protNHP50B >= 75]) - max(protNHP50B[protNHP50B <= 75]))
NHPid50B75 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50B75 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50B75])

upperpt <- min(log10NHPid50B[protNHP50B >= 90])
lowerpt <- max(log10NHPid50B[protNHP50B <= 90])
percway <- (90 - max(protNHP50B[protNHP50B <= 90]))/
           (min(protNHP50B[protNHP50B >= 90]) - max(protNHP50B[protNHP50B <= 90]))
NHPid50B90 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50B90 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50B90])


protNHP50C <- datNHP50C$protectlevel[kpC]

upperpt <- min(log10NHPid50C[protNHP50C >= 50])
lowerpt <- max(log10NHPid50C[protNHP50C <= 50])
percway <- (50 - max(protNHP50C[protNHP50C <= 50]))/
           (min(protNHP50C[protNHP50C >= 50]) - max(protNHP50C[protNHP50C <= 50]))
NHPid50C50 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50C50 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50C50])

upperpt <- min(log10NHPid50C[protNHP50C >= 75])
lowerpt <- max(log10NHPid50C[protNHP50C <= 75])
percway <- (75 - max(protNHP50C[protNHP50C <= 75]))/
           (min(protNHP50C[protNHP50C >= 75]) - max(protNHP50C[protNHP50C <= 75]))
NHPid50C75 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50C75 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50C75])

upperpt <- min(log10NHPid50C[protNHP50C >= 90])
lowerpt <- max(log10NHPid50C[protNHP50C <= 90])
percway <- (90 - max(protNHP50C[protNHP50C <= 90]))/
           (min(protNHP50C[protNHP50C >= 90]) - max(protNHP50C[protNHP50C <= 90]))
NHPid50C90 <- 10^(lowerpt + percway*(upperpt - lowerpt))
NHPind50C90 <- min(c(1:length(vAMP))[rev(vAMPOrig) >= NHPid50C90])


### The figure shows 3 lines for NHP for set A, B, C: 

protectlevelA <- datNHP50A$protectlevel
protectlevelB <- datNHP50B$protectlevel
protectlevelC <- datNHP50C$protectlevel


pdf(figureic50,width=0.95*7, height=0.95*6)
xlabs <- quantile(vAMPOrig,prob=seq(0,1,len=8))
xlabs <- c(round(xlabs[1:4],1),round(xlabs[5:10],0))
xlabs <- xtickAtquant
plot(d, 
     mark=NULL, 
     tx=NULL, 
     ylim=c(-0.4, 1),
     xtickAt=quantile(vAMP,prob=seq(0,1,len=length(xtickAtquant))),
     xtickLab=xlabs,
     ytickAt=seq(-0.4, 1, by=0.20),
     ytickLab=seq(-40, 100, by=20),
     xlab=expression(paste("Predicted ", ID[50]," against autologous virus")),
     ylab="Prevention Efficacy (%)",
     txLab=c("Placebo", "VRC01"),cex.axis=0.2)
lines(xaxisvalueNHPA,protectlevelA[kpA]/100,lwd=6,lty=1,col="blue")
lines(xaxisvalueNHPB,protectlevelB[kpB]/100,lwd=6,lty=1,col="orange")
lines(xaxisvalueNHPC,protectlevelC[kpC]/100,lwd=6,lty=1,col="green")
abline(h=0.5,lty=3,lwd=1)
text(0.1,0.52,"50% Efficacy")
text(0.1,0.77,"75% Efficacy")
text(0.1,0.92,"90% Efficacy")

abline(h=0.75,lty=3,lwd=1)
abline(h=0.90,lty=3,lwd=1)
#mtext(as.character(paste("bandwidth",hbandic50)), side=1, line=-6, adj=1, cex=1.2)
#mtext("BAMA for AMP, ELISA for NHP study", side=1, line=-3, adj=1, cex=1.1)
legend(0.6,-0.05,legend=c("AMP","NHP Set A", "NHP Set B", "NHP Set C"),col=c("black","blue","orange","green"),
lty=c(1,1,1,1),lwd=4,cex=1.1)
#title(title9)
dev.off()

q(save='no')
