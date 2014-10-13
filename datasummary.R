# TO RUN Script:
#      source('datasummary.R')
#
# install package for nifti file
#     install.packages("oro.nifti", repos="http://R-Forge.R-project.org")
# load nifti package
# library(oro.nifti)
library(car)
# jac2009 = readNIfTI("Brain_Stemlogjacobian2009.nii.gz", reorient=FALSE)

stopQuietly <- function(...)
  {
  blankMsg <- sprintf( "\r%s\r", paste( rep(" ", getOption( "width" ) - 1L ), collapse = " ") );
  stop( simpleError( blankMsg ) );
  } # stopQuietly()

#args <- commandArgs( trailingOnly = TRUE )
# manually pass in command line arguments for debugging
args <- c( "datasummaryL2_10sourceNewton46.txt","bestfit", ".3")
args <- c( "datasummary.txt","bestfit", ".7")
args <- c( "datasummaryL2_10sourceNewton49.txt","heating", ".3")
print(typeof(args))
print(args)

if( length( args ) < 3 )
  {
  cat( "Usage: Rscript datasummary.R datasummary.txt plotid dcethreshold", sep = "" )
  stopQuietly()
  }

# global vars

# read data
rawdata = read.table(args[1],header = TRUE, sep = ',')
PlotID=args[2]
dcetreshold= as.numeric( args[3] )
NumBreaks = 10

# for large objective functions and data on the boundary of optimization
# space most likely and error in the optimization
alphaconversion = 1.e7   # 1.e7 [cm^2 /s * 1e3] =  1 [m^2/s]
mueffconversion = 1.e-2  # 1.e-2 [1/cm]  = 1 [1/m]
lower_bound_alpha  = 1.19227e-07 * alphaconversion 
upper_bound_alpha  = 2.91296e-07 * alphaconversion 
lower_bound_mu_eff = 8.00000e-01 * mueffconversion 
upper_bound_mu_eff = 5.30000e+03 * mueffconversion 
lower_bound_robin  = 0.00000e+00 
upper_bound_robin  = 1.00000e+04 

rawdata$alpha  = rawdata$alpha  * alphaconversion
rawdata$mu_eff = rawdata$mu_eff * mueffconversion 
iterstats        =  subset(rawdata , obj<=5.e5   
                                   & alpha  > (1.34e-7*alphaconversion) 
                                   & mu_eff < (5.e3*mueffconversion)
                          )
                                   ##& alpha  < upper_bound_alpha
                                   ##& alpha  < (1.9e-7*alphaconversion)
                                   ##& alpha  > lower_bound_alpha
                                   ##& mu_eff < upper_bound_mu_eff 
                                   ##& mu_eff > lower_bound_mu_eff ) 
 
# view first 10 lines
NumDataPoints = length(iterstats$iddata)
print(head(iterstats,n=10))
print( paste(paste(paste("selected ", NumDataPoints ), " of "), length(rawdata$iddata)) )

# print columns
names(iterstats)

CrossOver=1.0
#color code
colcode <- character(nrow(iterstats))
## colcode[] <- "black"
## colcode[qoisubset$dirmethod == 'GR'] <- 1
## colcode[qoisubset$dirmethod == 'S2'] <- 3
## 
symcode <- integer(nrow(iterstats ))
## symcode[qoisubset$dirmethod == 'GR'] <- 1
## symcode[qoisubset$dirmethod == 'S2'] <- 3
## 
# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 2.5) )
    h <- hist(x, breaks=NumBreaks, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col='cyan', ...)
    meanvalue <- mean(x)
    localdigits <- 4
    avgtxt <- paste("avg=", format( meanvalue , digits = localdigits ))
    stdtxt <- paste("std=", format( sd(x)     , digits = localdigits ))
    mintxt <- paste("min=", format( min(x)    , digits = localdigits ))
    maxtxt <- paste("max=", format( max(x)    , digits = localdigits ))
    text(meanvalue , 1.2, mintxt )
    text(meanvalue , 1.4, avgtxt )
    text(meanvalue , 1.6, stdtxt )
    text(meanvalue , 1.8, maxtxt )
    if (do.legend) legend(43,2.5,c("GR","S2"), col = c(1,3), pch = c(1,3), bg="white")
    do.legend <<- FALSE
  
}
## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    #if(cor(x, y) <  0.0) txt <- paste("cor=-", txt)
    #if(cor(x, y) >= 0.0) txt <- paste("cor=" , txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    #text(0.5, 0.5, txt, cex = cex.cor * r)
    text(0.5, 0.5, txt, cex = 1.2* r)
}

# lexical scope on color code and symbol code
# http://www.johndcook.com/R_language_for_programmers.html#scope
panel.linear <- function(x, y)
{
        points(x,y)
        #points(x,y,col=colcode,pch=symcode)
	#abline(lm(y~x), col='red')
        abline(h=CrossOver,v=CrossOver,col='blue',lty=2)
}


# plot global summary plot
do.legend <- FALSE
pdf(paste(paste('datasummary',PlotID,sep=""),'.pdf',sep=""))
 #pairs(~obj+dice+alpha+mu_eff,data=rawdata,
 pairs(~obj+dice+alpha+mu_eff,data=iterstats,
        diag.panel  = panel.hist, 
        lower.panel = panel.linear,
        upper.panel = panel.cor,
        main=paste("N = "   ,NumDataPoints,sep="")
       )
dev.off()

parametersummary <- function(param, paramname,literature_value,unitdimension)
{
    pdf(paste(paste(paste('datasummary',paramname,sep=""),PlotID,sep=""),'.pdf',sep=""))
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 2.5) )
    #h <- hist(x, breaks=15, FALSE)
    histinfo <- hist(param, freq=FALSE, breaks=NumBreaks,
                     xlab=paste(paramname,unitdimension,sep=""), 
                     main=paste("N = "   ,NumDataPoints,sep=""), col="cyan")
    print("sum" )
    print(histinfo$breaks )
    print(histinfo$density )
    print(sum(diff(histinfo$breaks)*histinfo$density))
    #breaks <- h$breaks; nB <- length(breaks)
    #y <- h$counts; y <- y/max(y)
    #rect(breaks[-nB], 0, breaks[-1], y, col='cyan')
    curve(dnorm(x, mean=mean(param), sd=sd(param)), add=TRUE, col="darkblue", lwd=2)
    abline(v=literature_value,col='red',lty=2)
    meanvalue <- mean(param)
    localdigits <- 4
    avgstdtxt <- paste(paste(paste("avg=" , format( meanvalue , digits = localdigits )),
                                  " std="), format( sd(param) , digits = localdigits ))
    minmaxtxt <- paste(paste(paste("min=",  format( min(param), digits = localdigits )),
                                  " max="), format( max(param), digits = localdigits ))
    text(meanvalue , mean(histinfo$density)+1*  histinfo$density[1], minmaxtxt )
    text(meanvalue , mean(histinfo$density)+1.2*histinfo$density[1], avgstdtxt )
    ## if (do.legend) legend(43,2.5,c("GR","S2"), col = c(1,3), pch = c(1,3), bg="white")
    ## do.legend <<- FALSE
    dev.off()
  
}

# plot volume change 
initial_alpha  = 1.38546e-07  * alphaconversion 
initial_mu_eff = 1.80000e+02  * mueffconversion 
initial_robin  = 0.00000e+00 
parametersummary(iterstats$mu_eff ,'mueff',initial_mu_eff,' [1/cm]'   )
parametersummary(iterstats$alpha  ,'alpha',initial_alpha ,' [cm^2/s x 1.e3]' )

## plot(iterstats$alpha, iterstats$mu_eff, 
##   	 xlab="alpha", ylab="mu_eff", pch=3,cex=.5,
## abline(v=upper_bound_alpha ,h=upper_bound_mu_eff ,col='blue',lty=2)
## abline(v=lower_bound_alpha ,h=lower_bound_mu_eff ,col='blue',lty=2)
## abline(v=initial_alpha ,h=initial_mu_eff ,col='red',lty=2)
