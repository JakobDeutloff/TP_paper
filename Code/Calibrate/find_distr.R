# Identifies parameters of the probability distributions using the rrisk distributins package

install.packages('rriskDistributions')
library('rriskDistributions')

# mean tip threshold values
tipPFTP <- 4.0
tipAMOC <- 4.0
tipGRIS <- 1.5
tipWAIS <- 1.5
tipAMAZ <- 3.5
tipBORF <- 4.0
tipAWSI <- 6.3
tipEAIS <- 7.5
tipEASB <- 3.0
tipGLCR <- 2.0
tipLABC <- 1.8
tipTUND <- 4.0
tipPFAT <- 1.5
tipBARI <- 1.6
tipREEF <- 1.5
tipSAHL <- 2.8
# max tip threshold
tipPFTP_max <- 6.0
tipAMOC_max <- 8.0
tipGRIS_max <- 3.0
tipWAIS_max <- 3.0
tipAMAZ_max <- 6.0
tipBORF_max <- 5.0
tipAWSI_max <- 8.7
tipEAIS_max <- 10.0
tipEASB_max <- 6.0
tipGLCR_max <- 3.0
tipLABC_max <- 3.8
tipTUND_max <- 7.2
tipPFAT_max <- 2.3
tipBARI_max <- 1.7
tipREEF_max <- 2.0
tipSAHL_max <- 3.5
# min tip threshold
tipPFTP_min <- 3.0
tipAMOC_min <- 1.4
tipGRIS_min <- 0.8
tipWAIS_min <- 1.0
tipAMAZ_min <- 2.0
tipBORF_min <- 1.4
tipAWSI_min <- 4.5
tipEAIS_min <- 5.0
tipEASB_min <- 2.0
tipGLCR_min <- 1.5
tipLABC_min <- 1.1
tipTUND_min <- 1.5
tipPFAT_min <- 1.0
tipBARI_min <- 1.5
tipREEF_min <- 1.0
tipSAHL_min <- 2.0

# mean timescales 
yAMAZ = 100
yPFTP = 50
yPFAT = 200
#max tiescales
yAMAZ_max = 200
yPFTP_max = 300
yPFAT_max = 300
# min timescales
yAMAZ_min = 50
yPFTP_min = 10
yPFAT_min = 100 

# Set which quantiles min, mean and max threshold correspond to
quantiles <- c(0.05, 0.5, 0.95)

## Tipping points
# Get log-normal distributions
get.lnorm.par(p=quantiles, q=c(tipPFTP_min, tipPFTP, tipPFTP_max))
get.lnorm.par(p=quantiles, q=c(tipAMOC_min, tipAMOC, tipAMOC_max))
get.lnorm.par(p=quantiles, q=c(tipGRIS_min, tipGRIS, tipGRIS_max))
get.lnorm.par(p=quantiles, q=c(tipWAIS_min, tipWAIS, tipWAIS_max))
get.lnorm.par(p=quantiles, q=c(tipAMAZ_min, tipAMAZ, tipAMAZ_max))
get.lnorm.par(p=quantiles, q=c(tipBORF_min, tipBORF, tipBORF_max))
get.lnorm.par(p=quantiles, q=c(tipAWSI_min, tipAWSI, tipAWSI_max))
get.lnorm.par(p=quantiles, q=c(tipEAIS_min, tipEAIS, tipEAIS_max))
get.lnorm.par(p=quantiles, q=c(tipEASB_min, tipEASB, tipEASB_max))
get.lnorm.par(p=quantiles, q=c(tipGLCR_min, tipGLCR, tipGLCR_max))
get.lnorm.par(p=quantiles, q=c(tipLABC_min, tipLABC, tipLABC_max))
get.lnorm.par(p=quantiles, q=c(tipTUND_min, tipTUND, tipTUND_max))
get.lnorm.par(p=quantiles, q=c(tipPFAT_min, tipPFAT, tipPFAT_max))
get.lnorm.par(p=quantiles, q=c(tipBARI_min, tipBARI, tipBARI_max))
get.lnorm.par(p=quantiles, q=c(tipREEF_min, tipREEF, tipREEF_max))
get.lnorm.par(p=quantiles, q=c(tipSAHL_min, tipSAHL, tipSAHL_max))

# Try other distributions
distPFTP = fit.perc(p=quantiles, q=c(tipPFTP_min, tipPFTP, tipPFTP_max))
distAMOC = fit.perc(p=quantiles, q=c(tipAMOC_min, tipAMOC, tipAMOC_max))
distGRIS = fit.perc(p=quantiles, q=c(tipGRIS_min, tipGRIS, tipGRIS_max))
distWAIS = fit.perc(p=quantiles, q=c(tipWAIS_min, tipWAIS, tipWAIS_max))
distAMAZ = fit.perc(p=quantiles, q=c(tipAMAZ_min, tipAMAZ, tipAMAZ_max))
distBORF = fit.perc(p=quantiles, q=c(tipBORF_min, tipBORF, tipBORF_max))
distAWSI = fit.perc(p=quantiles, q=c(tipAWSI_min, tipAWSI, tipAWSI_max))
distEAIS = fit.perc(p=quantiles, q=c(tipEAIS_min, tipEAIS, tipEAIS_max))
distEASB = fit.perc(p=quantiles, q=c(tipEASB_min, tipEASB, tipEASB_max))
distGLCR = fit.perc(p=quantiles, q=c(tipGLCR_min, tipGLCR, tipGLCR_max))
distLABC = fit.perc(p=quantiles, q=c(tipLABC_min, tipLABC, tipLABC_max))
distTUND = fit.perc(p=quantiles, q=c(tipTUND_min, tipTUND, tipTUND_max))
distPFAT = fit.perc(p=quantiles, q=c(tipPFAT_min, tipPFAT, tipPFAT_max))
distBARI = fit.perc(p=quantiles, q=c(tipBARI_min, tipBARI, tipBARI_max))
distREEF = fit.perc(p=quantiles, q=c(tipREEF_min, tipREEF, tipREEF_max))
distSAHL = fit.perc(p=quantiles, q=c(tipSAHL_min, tipSAHL, tipSAHL_max))

## Timescales
# Get distributions
dist_y_PFTP = fit.perc(p=quantiles, q=c(yPFTP_min, yPFTP, yPFTP_max))
dist_y_AMAZ = fit.perc(p=quantiles, q=c(yAMAZ_min, yAMAZ, yAMAZ_max))
dist_y_PFAT = fit.perc(p=quantiles, q=c(yPFAT_min, yPFAT, yPFAT_max))

