#' Plotting Mark-Specific Proportional Hazards Model Fits
#'
#' \code{plot} method for class \code{summary.sievePH}. For univariate marks, it plots point and interval estimates of the mark-specific treatment effect parameter specified by \code{contrast} in \code{\link{summary.sievePH}}, and,
#' optionally, scatter/box plots of the observed mark values by treatment. For bivariate marks, plotting is restricted to the point estimate, which is displayed as a surface. No plotting is provided for marks of higher dimensions.
#'
#' @param x an object returned by \code{\link{summary.sievePH}}
#' @param mark either a numeric vector specifying a univariate continuous mark or a data frame specifying a multivariate continuous mark.
#' For subjects with a right-censored time-to-event, the value(s) in \code{mark} should be set to \code{NA}.
#' @param tx a numeric vector indicating the treatment group (1 if treatment, 0 if placebo)
#' @param xlim a numeric vector of length 2 specifying the x-axis range (\code{NULL} by default)
#' @param ylim a numeric vector of length 2 specifying the y-axis range (\code{NULL} by default)
#' @param zlim a numeric vector of length 2 specifying the z-axis range in a 3-dimensional plot (\code{NULL} by default)
#' @param xtickAt a numeric vector specifing the position of x-axis tickmarks (\code{NULL} by default)
#' @param xtickLab a numeric vector specifying labels for tickmarks listed in \code{xtickAt}. If \code{NULL} (default), the labels are determined by \code{xtickAt}.
#' @param ytickAt a numeric vector specifing the position of y-axis tickmarks (\code{NULL} by default)
#' @param ytickLab a numeric vector specifying labels for tickmarks listed in \code{ytickAt}. If \code{NULL} (default), the labels are determined by \code{ytickAt}.
#' @param xlab a character string specifying the x-axis label (\code{NULL} by default)
#' @param ylab a character string specifying the y-axis label (\code{NULL} by default)
#' @param zlab a character string specifying the z-axis label in a 3-dimensional plot (\code{NULL} by default)
#' @param txLab a character vector of length 2 specifying the placebo and treatment labels (in this order). The default labels are \code{placebo} and \code{treatment}.
#' @param title a character string specifying the plot title (\code{NULL} by default)
#' @param ... other arguments to be passed to plotting functions
#'
#' @details
#' For bivariate marks, \code{markGrid} in \code{\link{summary.sievePH}} must have equally spaced values for each component.
#'
#' @return None. The function is called solely for plot generation.
#'
#' @examples
#' n <- 500
#' tx <- rep(0:1, each=n/2)
#' tm <- c(rexp(n/2, 0.2), rexp(n/2, 0.2 * exp(-0.4)))
#' cens <- runif(n, 0, 15)
#' eventTime <- pmin(tm, cens, 3)
#' eventInd <- as.numeric(tm <= pmin(cens, 3))
#' mark <- ifelse(eventInd==1, c(rbeta(n/2, 2, 5), rbeta(n/2, 2, 2)), NA)
#' markRng <- range(mark, na.rm=TRUE)
#'
#' # fit a model with a univariate mark
#' fit <- sievePH(eventTime, eventInd, mark, tx)
#' sfit <- summary(fit, markGrid=seq(markRng[1], markRng[2], length.out=10))
#' plot(sfit, mark, tx)
#'
#' @seealso \code{\link{sievePH}}, \code{\link{sievePHipw}}, \code{\link{sievePHaipw}} and \code{\link{summary.sievePH}}
#'
#' @export
plot.summary.sievePH <- function(x, mark=NULL, tx=NULL, xlim=NULL, ylim=NULL, zlim=NULL, xtickAt=NULL, xtickLab=NULL, ytickAt=NULL, ytickLab=NULL, 
                                 xlab=NULL, ylab=NULL, zlab=NULL, txLab=c("Placebo", "Treatment"), title=NULL, 
                                 compTE=NULL, sec.xtickAt=NULL, sec.xtickLab=NULL, sec.xlab=NULL,
                                 parMar=c(5, 6, 2, 1), ...){
  contrast <- names(x)[length(names(x))]
  
  cexAxis <- 1.2
  cexLab <- 1.4
  cexTitle <- 1.6
  cexText <- 1.2
  cexLegend <- 1.2
  
  par(mar=parMar, oma=rep(0, 4), cex.axis=cexAxis, cex.lab=cexLab, cex.main=cexTitle)
  
  # a 2-dimensional plot only when the mark is univariate
  if (NCOL(x[[contrast]])==4){
    if (is.null(xlim)){
      xlim <- range(x[[contrast]][, 1])
    }
    
    if (is.null(ylim)){
      ylim <- range(x[[contrast]][, -1], na.rm=TRUE)
    }
    ySplit <- ylim[2]
    
    # need extra room for box plots on the top
    if (!any(c(is.null(mark), is.null(tx)))){
      ylim <- c(ylim[1], ylim[2] + ifelse(is.null(sec.xtickLab), 0.35, 0.35) * (ylim[2] - ylim[1]))
    }
    
    if (is.null(xlab)){ xlab <- colnames(x[[contrast]])[1] }
    if (is.null(ylab)){ ylab <- switch(colnames(x[[contrast]])[2], TE="Treatment Efficacy", HR="Hazard Ratio", LogHR="Log Hazard Ratio") }
    
    plot(x[[contrast]][, 1], x[[contrast]][, 2], xlim=xlim, ylim=ylim, type="n", xlab="", ylab="", xaxt=ifelse(is.null(xtickAt), "s", "n"), yaxt="n", bty="l", ...)
    
    if (!is.null(xtickAt)){
      if (is.null(xtickLab)){ xtickLab <- xtickAt }
      axis(side=1, at=xtickAt[seq(1, length(xtickAt), by=2)], labels=xtickLab[seq(1, length(xtickLab), by=2)], cex.axis=cexAxis)
      axis(side=1, at=xtickAt[seq(2, length(xtickAt), by=2)], labels=xtickLab[seq(2, length(xtickLab), by=2)], cex.axis=cexAxis)
    }
    
    if (!is.null(sec.xtickAt) && !is.null(sec.xtickLab)){
      axis(side=1, line=4.5, at=sec.xtickAt[seq(1, length(sec.xtickAt), by=2)], labels=sec.xtickLab[seq(1, length(sec.xtickLab), by=2)], cex.axis=cexAxis, 
           col="magenta3", col.ticks="magenta3", col.axis="magenta3")
      axis(side=1, line=4.5, at=sec.xtickAt[seq(2, length(sec.xtickAt), by=2)], labels=sec.xtickLab[seq(2, length(sec.xtickLab), by=2)], cex.axis=cexAxis, 
           col="magenta3", col.ticks="magenta3", col.axis="magenta3")
    }
    
    if (!is.null(ytickAt)){
      if (is.null(ytickLab)){ ytickLab <- ytickAt }
      axis(side=2, at=ytickAt, labels=ytickLab, las=1, cex.axis=cexAxis)
    } else {
      # to avoid overlapping tickmarks
      if (!any(c(is.null(mark), is.null(tx)))){
        axis(side=2, at=axTicks(2)[axTicks(2) <= 1.1 * max(x[[contrast]][, 4], na.rm=TRUE)])
      }
    }
    
    if (is.null(sec.xlab)){
      mtext(xlab, side=1, line=3, cex=cexLab)
    } else { 
      mtext(xlab, side=1, line=2.7, cex=cexLab)
      mtext(sec.xlab, side=1, line=7, cex=cexLab, col="magenta3") 
    }
    
    mtext(ylab, side=2, line=3.5, las=3, cex=cexLab)
    
    if (!is.null(title)){ mtext(title, side=3, font=2, line=1, cex=cexTitle, outer=TRUE, at=0, adj=0) }
    
    abline(h=ifelse(colnames(x[[contrast]])[2]=="HR", 1, 0), col="gray70", lwd=2)
    
    ylimGap <- 0
    if (!is.null(compTE)){
      ylimGap <- 0.1
      # colCI <- "darkgoldenrod2"
      # colRGB <- c(col2rgb(colCI))
      # colRGB <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha=255*0.55, maxColorValue=255)
      # polygon(c(x[[contrast]][, 1], rev(x[[contrast]][, 1])), 
      #         c(ifelse(compTE$LB2 >= ylim[1] & compTE$LB2 <= ySplit, compTE$LB2, ylim[1]), rev(ifelse(compTE$UB2 <= 1, compTE$UB2, 1))), col=colRGB, border=NA)
      lines(x[[contrast]][, 1], ifelse(compTE$TE2 >= ylim[1] + ylimGap & compTE$TE2 <= ySplit, compTE$TE2, NA), lwd=4, col="magenta3")
    }
    
    lines(x[[contrast]][, 1], ifelse(x[[contrast]][, 2] >= ylim[1] + ylimGap & x[[contrast]][, 2] <= ySplit, x[[contrast]][, 2], NA), lwd=4)
    lines(x[[contrast]][, 1], ifelse(x[[contrast]][, 3] >= ylim[1] + ylimGap & x[[contrast]][, 3] <= ySplit, x[[contrast]][, 3], NA), lwd=3.5, lty="dashed")
    lines(x[[contrast]][, 1], ifelse(x[[contrast]][, 4] >= ylim[1] + ylimGap & x[[contrast]][, 4] <= ySplit, x[[contrast]][, 4], NA), lwd=3.5, lty="dashed")
    
    # text(min(out$v), -0.6, paste0("Marginal Sieve Test P ",ifelse(pMarginalSieve<0.001,"< 0.001",paste0("= ",format(pMarginalSieve,digits=2))),marginalSignifMark), pos=4, cex=cexText)
    
    #legend("bottomleft", fill=colCI, border=colCI, legend="95% Pointwise CI", cex=cexLegend, bty="n")
    #legend("bottomleft", lwd=c(3.5,2), lty=c("dashed","longdash"), legend=c("95% Pointwise CI","Overall Hazard-Ratio PE"), col=c("black","darkorange"), cex=cexLegend, bty="n")
    if (is.null(compTE)){
      legend(x=xlim[1], y=ifelse(colnames(x[[contrast]])[2]=="TE", ylim[1] + 0.05 * (ylim[2] - ylim[1]), ySplit), lwd=3.5, lty="dashed", 
             legend="95% Pointwise CI", col="black", cex=cexLegend, bty="n")  
    } else {
      legend(x=xlim[1], y=ifelse(colnames(x[[contrast]])[2]=="TE", ylim[1] + 0.12 * (ylim[2] - ylim[1]), ySplit), lwd=c(3.5, 4), lty=c("dashed", "solid"), 
             legend=c("95% Pointwise CI", expression("Est. PE by Measured" ~ IC[80] ~ "(Quantile-Matched)")), col=c("black", "magenta3"), cex=cexLegend, bty="n")
    }
    
    # add scatter/box plots of the observed mark values by treatment
    if (!any(c(is.null(mark), is.null(tx)))){
      if (is.null(sec.xtickLab)){
        par(fig=c(0,1,0.7,1), new=TRUE)  
      } else {
        par(fig=c(0,1,0.75,1), new=TRUE)  
      }
      
      data <- na.omit(cbind(mark, tx))
      plotMarkHoriz(data[, 1], data[, 2], parMar=c(0.5, parMar[-1]), yLim=xlim, txLab=txLab)
      
      par(fig=c(0,1,0,1), new=TRUE)
    }
    
    # a 3-dimensional plot (a surface) when the mark is bivariate
  } else if (NCOL(x[[contrast]])==5){
    if (is.null(xlim)){
      xlim <- range(x[[contrast]][, 1])
    }
    
    if (is.null(ylim)){
      ylim <- range(x[[contrast]][, 2])
    }
    
    if (is.null(zlim)){
      zlim <- range(x[[contrast]][, 3])
    }
    
    if (is.null(xlab)){ xlab <- colnames(x[[contrast]])[1] }
    if (is.null(ylab)){ ylab <- paste0("\n", colnames(x[[contrast]])[2]) }
    if (is.null(zlab)){ zlab <- switch(colnames(x[[contrast]])[3], TE="\n\nTreatment Efficacy", HR="\n\nHazard Ratio", LogHR="\n\nLog Hazard Ratio") }
    
    # the first two arguments must be vectors with equally spaced values in ascending order
    persp(sort(unique(x[[contrast]][, 1])), sort(unique(x[[contrast]][, 2])), getOuterProduct(x[[contrast]], zlim),
          xlab=xlab, ylab=ylab, zlab=zlab, col="lightgreen", theta=150, phi=20, ticktype="detailed",
          nticks=5, xlim=xlim, ylim=ylim, zlim=zlim, r=3, expand=0.8, main=title)
  } else {
    stop("Plotting of results is available for univariate and bivariate marks only.")
  }
}

plotMarkHoriz <- function(mark, tx, parMar, yLim, txLab=c("Placebo", "Treatment")){
  cexAxis <- 1.3
  cexLab <- 1.4
  
  par(mar=parMar, oma=c(0,0,0,0), cex.axis=cexAxis, cex.lab=cexLab)
  # vioplot(mark[tx==0], mark[tx==1], at=c(0.5, 1.5), names=NA, horizontal=TRUE, drawRect=TRUE, col="white", border=c("blue", "red3"), rectCol=c("blue", "red3"), 
  #         lineCol=c("blue", "red3"), frame.plot=FALSE, yaxt="n")
  boxplot(mark ~ as.factor(tx), at=c(0.5,1.5), xlim=c(0,2), ylim=c(yLim[1], yLim[2]), frame.plot=FALSE, xaxt="n", yaxt="n",
          xlab="", ylab="", boxwex=0.75, outline=FALSE, border="black", lwd=2.5, horizontal=TRUE)
  axis(side=2, at=c(0.5,1.5), labels=txLab, cex.axis=cexAxis, las=1)
  points(mark, jitter(tx + 0.5, factor=0.7), col=ifelse(tx==1, "red3", "blue"), pch=ifelse(tx==1, 24, 21), lwd=1.5, cex=ifelse(tx==1, 1, 1), bg="white")
  # points(mark, jitter(tx + 0.5, factor=0.6), col=ifelse(tx==1, "red3", "blue"), pch=20, lwd=2, cex=0.6)
}

getOuterProduct <- function(df, zlim){
  mark1 <- sort(unique(df[, 1]))
  mark2 <- sort(unique(df[, 2]))
  
  out <- matrix(NA, nrow=length(mark1), ncol=length(mark2))
  for (i in 1:length(mark1)){
    for (j in 1:length(mark2)){
      idx <- which(df[, 1]==mark1[i] & df[, 2]==mark2[j])
      if (length(idx) > 1){
        stop("There are replicates on the marker grid.")
      } else if (length(idx)==1){
        out[i, j] <- df[idx, 3]
      }
      
    }
  }
  out <- ifelse(out >= zlim[1] & out <= zlim[2], out, NA)
  return(out)
}
