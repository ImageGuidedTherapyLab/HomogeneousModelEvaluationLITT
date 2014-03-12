# TO RUN Script:
#      source('processhistory')
#
# install package for nifti file
#     install.packages("oro.nifti", repos="http://R-Forge.R-project.org")
# load nifti package
# library(oro.nifti)
library(car)
# jac2009 = readNIfTI("Brain_Stemlogjacobian2009.nii.gz", reorient=FALSE)
PlotID='007'
iterstats = read.table('workdir/Patient0006/007/opt/iterationhistory.txt',header = TRUE, sep = ',')



# view first 10 lines
print(head(iterstats,n=10))

# print columns
names(iterstats)

## # summary stats of patient data
## patientstats = read.table("../BrainCases.csv"  , header = FALSE, sep = ';',skip=1)
## print(paste("num=", format(length( patientstats$V4))))
## print(paste("avg=", format(mean(   patientstats$V4))))
## print(paste("med=", format(median( patientstats$V4))))
## print(paste("std=", format(sd(     patientstats$V4))))
## print(paste("min=", format(min(    patientstats$V4))))
## print(paste("max=", format(max(    patientstats$V4))))
## print(summary(patientstats$V4))
## 
## # extract subsets
## wholebrainlabel                         =  subset(iterstats , labeltype=='ICBM_FullBrain' & jacobian=='full' ) 
## Cerebrum				=  subset(wholebrainlabel, LabelID==1  ) 
## Brain_Stem                            	=  subset(wholebrainlabel, LabelID==58 ) 
## Cerebellum                          	=  subset(wholebrainlabel, LabelID==67 ) 
## Ventricle                              	=  subset(wholebrainlabel, LabelID==255) 
## 
## # rename columns
## names(Cerebrum)[  names(Cerebrum)   == 'Mean'] <- 'Cerebrum'
## names(Brain_Stem)[names(Brain_Stem) == 'Mean'] <- 'Brain_Stem'
## names(Cerebellum)[names(Cerebellum) == 'Mean'] <- 'Cerebellum'
## names(Ventricle )[names(Ventricle)  == 'Mean'] <- 'Ventricle'
## 
## # concatenate multiple data anatomy
## dataqoi = cbind(Cerebrum,  Brain_Stem,Cerebellum , Ventricle )
## qoisubset = subset(dataqoi , studytime <= 400.0  &mrn!=888    );PlotID='TimeZeroFull';
CrossOver=1.0
## 
## print(paste("TE num=", format(length( dataqoi$TE))))
## print(paste("TE avg=", format(mean(   dataqoi$TE))))
## print(paste("TE med=", format(median( dataqoi$TE))))
## print(paste("TE std=", format(sd(     dataqoi$TE))))
## print(paste("TE min=", format(min(    dataqoi$TE))))
## print(paste("TE max=", format(max(    dataqoi$TE))))
## 
## print(paste("TR num=", format(length( dataqoi$TR))))
## print(paste("TR avg=", format(mean(   dataqoi$TR))))
## print(paste("TR med=", format(median( dataqoi$TR))))
## print(paste("TR std=", format(sd(     dataqoi$TR))))
## print(paste("TR min=", format(min(    dataqoi$TR))))
## print(paste("TR max=", format(max(    dataqoi$TR))))
## 
## # count images columns
## counttotalimages                            =  subset(dataqoi , dirmethod =='GR' & metrictype =='CC' & LabelID==1 ) 
## countsubsetimages                           =  subset(qoisubset , dirmethod =='GR' & metrictype =='CC' & LabelID==1 ) 
## print(paste("total  images=", format(length( counttotalimages$Ventricle))))
## print(paste("subset images=", format(length( countsubsetimages$Ventricle))))
## 
## #qoisubset = subset(dataqoi , subset = PtId!=947329  & ImageBase =='timezero'  & Jac=='projection'  );PlotID='TimeZeroProj';CrossOver=1.0
## #qoisubset = subset(dataqoi , subset = PtId!=947329  & ImageBase =='template'  & Jac=='full'  );PlotID='TemplateFull';CrossOver=0.0
## #qoisubset = subset(dataqoi , subset = PtId!=947329  & ImageBase =='template'  & Jac=='projection'  );PlotID='TemplateProj';CrossOver=0.0
## 
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
    h <- hist(x, plot = FALSE)
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
    if(cor(x, y) <  0.0) txt <- paste("cor=-", txt)
    if(cor(x, y) >= 0.0) txt <- paste("cor=" , txt)
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
	abline(lm(y~x), col='red')
        abline(h=CrossOver,v=CrossOver,col='blue',lty=2)
}


# plot volume change 
do.legend <- FALSE
pdf(paste(paste('/var/www/fuentes/ipresults/processhistory',PlotID,sep=""),'.pdf',sep=""))
 pairs(~iter+mu_eff+k+dz+obj,data=iterstats,
        diag.panel  = panel.hist, 
        lower.panel = panel.linear,
        upper.panel = panel.cor,
       )
## scatterplotMatrix(~PtNum+Time+Lobes+White_Matter+CSF|Dir,data=dataqoi,
##          diagonal='density',
##          subset = PtId!=947329  & ImageBase =='timezero'  & Jac=='full' ,
##   	 main="Three Cylinder Options")
dev.off()
## #      source('stats.R')
## 
## do.legend <- TRUE
## pdf(paste(paste('/var/www/fuentes/CohortWholeTime',PlotID,sep=""),'.pdf',sep=""))
## # pairs(~PtNum+age+Cerebrum+Cerebellum+Brain_Stem+CSF,data=qoisubset,
##  pairs(~age+studytime+dose+Cerebrum+Ventricle,data=dataqoi,
##         diag.panel  = panel.hist, 
##         lower.panel = panel.linear,
##         upper.panel = panel.cor,
##        )
## dev.off()
