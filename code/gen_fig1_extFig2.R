### Generate Fig. 1 and Extended Data Fig. 2
library(magrittr)
library(dplyr)

source('CCC_time.R')

df = read.csv('../data/AMP_autologous_neutralization.csv')

# For display compute categories, truncate predicted titer at 0.001 and jitter measured titers that were below LOD
df$cat = 1 + ifelse(df$sensitivity < 1, 0, 1) + ifelse(df$sensitivity > 3, 1, 0)
df$pred.titer.trunc = pmax(0.001, df$pred.titer)
df$titer_j = df$titer + ifelse(df$titer == 5, runif(nrow(df), -0.3,0.3), 0)

# create 3 distinct colors (like default 3 color scheme in ggplot)
COL = c('green','blue', 'red')

xlim = range(c(df$titer))
ylim = range(df$pred.titer)
ylim = xlim

# When experimental titer is greater than the limit of detection the predicted titer
# is > 8.7
sss = subset(df, titer >= 10)
min(sss$pred.titer)

# When the experimental titer is at the LOD the predicted titer is <= 36
sss = subset(df, titer==5)
max(sss$pred.titer)

AT = c(0.001,0.01, 0.1, 1, 10, 100)
LEGEND = c(expression(paste("<1 ", mu, "g/mL")), expression(paste("1-3 ", mu, "g/mL")), expression(paste(">3 ", mu, "g/mL")))

set.seed(1)

for(pc in c(80, 50)) {
  
  if( pc == 80 ) {
    pdf('../figures/fig1.pdf')
  } else {
    pdf('../figures/extFig2.pdf')
  }
  ss = subset(df, poscrit==pc)
  sss = subset(ss, titer>=10, select=c('pub_id','isolate','visitno','titer','pred.titer'))
  
  # compute on isolate level
  sss = sss %>% group_by(isolate) %>%
    arrange(isolate, visitno) %>% 
    mutate(n=1, rep = cumsum(n))
  
  est = CCC_time(sss$titer, sss$pred.titer, sss$rep)
  
  t = table('Experimental' = ifelse(ss$titer>=10, '>=10', '<10'), 'Predicted'=ifelse(ss$pred.titer>=10, '>=10', '<10'))
  
  ylim = range(ss$pred.titer) 
  ylim[1] = max(0.001, ylim[1])
  
  plot(ss$titer_j, ss$pred.titer.trunc, pch=ss$cat+14, log='xy', col=COL[ss$cat],
       xlab=sprintf('Experimental Serum ID%d Titer',pc), ylab=sprintf('Predicted Serum ID%d Titer (PT%d)', pc, pc),
       xlim=xlim, ylim=ylim, axes=F)
  box()
  axis(1, at=c(1, 5, 10, 20, 50, 100, 200))
  axis(2, at=AT[-1], labels=as.character(AT[-1]), las=2)
  axis(2, at=0.001, labels=expression(""<=0.001), las=2)
  abline(0,1)
  abline(v=10, h=10, lty=3)
  legend(x='bottomright', bg='white', title=sprintf('Virus category (IC%d)', pc), legend=LEGEND, pch=15:17, col=COL[1:3])
  
  print("******************************************\n")
  print(sprintf("n=%d experimental IC%d values from %d infected VRC01 recipient cases and %d synthesized viruses", 
                nrow(ss), pc, length(unique(ss$pub_id)), length(unique(ss$isolate))))
  print(sprintf("CCC %0.2f (95%% CI %0.2f to %0.2f)", est$CCC, est$CCC.l, est$CCC.u))
  print(t)
  print(sprintf("%d of %d Experimental ID%d titers where below the limit of detection when the Predicted ID%d was greater than 10 (range %0.1f to %0.1f)",
                t[1,2], sum(t), pc, pc, min(subset(ss, titer==5 & pred.titer >=10)$pred.titer),  max(subset(ss, titer==5 & pred.titer >=10)$pred.titer)))
  print(sprintf("%d of %d Experimental ID%d titers where at or above the limit of detection when the Predicted ID%d was less than 10 (range %0.1f to %0.1f)",
                t[2,1], sum(t), pc, pc, min(subset(ss, titer>=10 & pred.titer < 10)$pred.titer),  max(subset(ss, titer>=10 & pred.titer<10)$pred.titer)))

  dev.off()
}

q(save='no')


