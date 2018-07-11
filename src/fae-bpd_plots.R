#-----------------------------------------------------------------------
# This analyzes and plots the behavioral part of the attribution study
# in patients with BPD and healthy controls
# PH, 1/7/17
#-----------------------------------------------------------------------
rm(list=ls())
library('plyr')
library('lme4')
library('lsmeans')

plotreg <- function(df, lty, lwd, pch) {
  # Linear regression plots
  #
  #   Args:
  #
  #   Returns:
  #  
  reg    <- lm(df$y ~ df$x)
  p      <- points(x=df$x, y=df$y, pch=pch)
  a      <- abline(reg, lwd=lwd, lty=lty)
  return(reg)
}

plotdot <- function(df, pch=19, lwd=2, cex=0.75, l=FALSE, lty=1) {
  # Bar plots with error bars
  #
  #   Args:
  #
  #   Returns:
  tmp  <- points(x=df$x, y=df$y, pch=pch, cex=cex)
  if (!is.na(df$eb[1])) {
    mebs    <- df$y - df$eb
    pebs    <- df$y + df$eb
    a       <- arrows(df$x, mebs, df$x, pebs, length=0.00,
                      angle=90, code=3, lwd=lwd)
  }
  if (l == TRUE) {
    tmp    <- lines(df$x, df$y, lty=lty)
  }
} 
    
plotbar <- function(df, cols, factors, xlim, ylim, axes=TRUE,
                    xpd=TRUE, abl=NULL) {
  # Bar plots with error bars
  #
  #   Args:
  #
  #   Returns:
  # create means and se
  # restrict dataframe to columns of interest
  #mxy         <- max(df$m)*1.45
  tabbedmeans <- tapply(df$m, factors, function(x) c(x=x))
  tabbedeb    <- tapply(df$eb, factors, function(x) c(x=x))
  barcenters  <- barplot(height = tabbedmeans, beside = TRUE,
                        las = 1, axes=axes, col=cols, space=c(0, 0.4),
                        cex.names=0.8, ylim=ylim, xpd=xpd)
  if (!is.na(df$eb[1])) {
    s           <- segments(barcenters, tabbedmeans-tabbedeb,
                            barcenters, tabbedmeans+tabbedeb,
                            lwd=1.0)
    a           <- arrows  (barcenters, tabbedmeans-tabbedeb,
                            barcenters, tabbedmeans+tabbedeb,
                            lwd=1.0, angle=90, code=3, length=0.02)
  }
  #axis(1, at=1.5, labels="")
  if (!is.null(abl)) abline(h=abl, lwd=2)
  #lines(0:6, rep(-0.01, 7), lwd=2)
  #box(bty="l")
  #box(bty="7", col="white")
} 

calcmlm <- function(df=NULL){
  # Tests main hypothesis
  # Args:
  # Returns:
  # Build a subset for this test
  require('lmerTest')
  sdf <- subset(df, qtypen==1 &
                    ((expversion==1 & run==2 & (trial==1 | trial==2)) |
                     (expversion==2 & run==1 & (trial==1 | trial==2)) |
                     (expversion==3 & run==2 & (trial==1 | trial==2)) |
                     (expversion==4 & run==1 & (trial==1 | trial==4))))
  sdf$id <- factor(sdf$id)
  sdf$expversion <- factor(sdf$expversion)
  sdf$group <- factor(sdf$group)
  sdf$valence <- factor(sdf$valence)
  sdf$trait.attributions.val <- sdf$trait.attributions.val * 100  
  m1 <- lmer(trait.attributions.val ~ group * valence +
               (1 |id), sdf)
  print(anova(m1))
  print(summary(m1))

  # produce custom contrast plots
  detach("package:lmerTest", unload=TRUE)
  fcts     <- c("group", "valence")
  cv       <- matrix(c( -1,  -1,  1,  1, 
                         1,  -1,  0,  0,
                         0,   0,  1, -1,
                        -1,   1, -1,  1),
                     nrow=4, ncol=4, byrow=TRUE)

  cl       <- list(
    "Pos minus Neg"       =cv[1, ],
    "HCC minus BPD, Neg"  =cv[2, ],
    "HCC minus BPD, Pos"  =cv[3, ],
    "HCC minus BPD"       =cv[4, ]
  )

  x2lab    <- ""
  x1lab    <- "% Trait Attributions\n(adjusted Means with 95% CI)"

  plotcontrasts(model=m1, fcts=fcts, cl=cl, x1lab=x1lab)

  require('lmerTest')
  m2 <- lmer(trait.attributions.val ~ group * valencen * empathy +
               (1 |id), sdf)
  print(anova(m2))

  m3 <- lmer(trait.attributions.val ~ group * valencen * paranoia +
               (1 |id), sdf)
  print(anova(m3))

  m4 <- lmer(trait.attributions.val ~ group * valencen +
               need.for.cognition + (1 |id), sdf)
  print(anova(m4))

  # Test group effect on like ratings
  m5a <- lmer(resp ~ group * valence  +
               (1 + valence |id), subset(df, qtype=="like"))
  print(anova(m5a))

  m5b <- lmer(resp ~ group * valence * ttype +
               (1 + valence * ttype|id), subset(df, qtype=="like"))
  print(anova(m5b))

  #m5b1 <- lmer(resp ~ valence * ttype +
  #             (1 + valence * ttype|id), subset(df, qtype=="like"))
  #print(anova(m5b1))
  
  
  detach("package:lmerTest", unload=TRUE)
  m5b <- lmer(resp ~ group * valence * ttype +
               (1 + valence * ttype|id), subset(df, qtype=="like"))
  lsm <- lsmeans(m5b, c("group", "valence", "ttype"))
  lsm.contr <- contrast(lsm, list("PH1"=c(-1, -1, 0, 0, 1, 1, 0, 0),
                                  "PH2"=c( 0,  0, 1, 1, 0, 0,-1,-1)))
  print(lsm.contr)

  # produce custom contrast plots
  fcts     <- c("group", "valence", "ttype")
  cv       <- matrix(c( -1,  -1,  1,  1,  0,  0,  0,  0,  
                         0,   0,  0,  0, -1, -1,  1,  1,
                        -1,  -1,  0,  0,  1,  1,  0,  0,
                         0,   0, -1, -1,  0,  0,  1,  1),
                     nrow=4, ncol=8, byrow=TRUE)

  cl       <- list(
    "Sit minus Disp, Pos"     =cv[4, ],
    "Sit minus Disp, Neg"     =cv[3, ],
    "Pos minus Neg, Sit"      =cv[2, ],
    "Pos minus Neg, Disp"     =cv[1, ]
  )

  x1lab    <- "Like Evaluations\n(adjusted Means with 95% CI)"
  plotcontrasts(model=m5b, fcts=fcts, cl=cl, xlim=c(-4, 4), x1lab=x1lab)

  
  require('lmerTest')
  # Test for order effect
  df$cueorder <- factor(df$cueorder)
  m6 <- lmer(resp ~ group * valence * qtype + cueorder +
               (1 + qtype * valence |id), df)
  print(anova(m6))

  # Test response times
  df$cueorder <- factor(df$cueorder)
  m7 <- lmer(likert ~ group * valence  + trial + run +
               (1 + valence + trial + run |id), subset(df, qtype=="like"))
  print(anova(m7))
  
  m8 <- lmer(attribrt ~ group * valence  + trial + run +
               (1 + valence + trial + run |id), subset(df, qtype=="attrib"))
  print(anova(m8))
  detach("package:lmerTest", unload=TRUE)
}

plotcontrasts <- function(model=NULL, fn=NULL, fcts=NULL, cl,
                          x1lab=NULL, x2lab=NULL, eb="ci",
                          xlim=c(-70, 70)) {
  # Plots custom contrasts
  #
  # Args:
  # Returns:

  par(mar=c(8, 15, 8, 2))
  lsm       <- lsmeans(model, fcts)
  tmp       <- print(lsm)
  lsm.contr <- contrast(lsm, cl)
  tmp       <- print(summary(lsm.contr))
  lsm.ci    <- confint(lsm.contr)
  lsmdf     <- lsm.ci[c("contrast", "estimate", "SE", "lower.CL",
                        "upper.CL")]

  # create dot plot
  if (is.null(x1lab)) {
    x1lab     <- "% Trait Attributions \n (adjusted Means with 95% CI)"
  }
  if (is.null(x2lab)) {
    #x2lab     <- "Differential valence \n (Positive minus Negative behavior)"
    x2lab     <- "Custom Contrasts"
  }
  tmp       <- plotcontrdots(lsmdf, xlim=xlim, x1lab=x1lab, x2lab=x2lab, eb)
  
  if(!is.null(fn)) dev.off()
}

plotcontrdots <- function(df=NULL, xlim=c(-0.5, 0.5),
                     x1lab=NULL, x2lab=NULL, eb="ci") {
  # Plots forest plot of contrasts
  # Args:
  # Returns:
  tmp <- dotchart(df$estimate, xlim=xlim, pch=19, font.lab=2)
  for (i in 1:nrow(df)) { 
    switch(eb,
           ci = { lines(x=c(df$lower.CL[i], df$upper.CL[i]),
                        y=c(i, i), lwd=2)},
           se = { lines(x=c(df$estimate[i] - df$SE[i],
                            df$estimate[i] + df$SE[i]),
                        y=c(i, i), lwd=2)})
  }
  mtext(x1lab, 1, at=0, font=2, line=3.5)
  axis(2, at=(1:nrow(df)), labels=paste(df$contrast), las=2)
  abline(v=0, lwd=1, lty=2)
  box(lwd=2)
  mtext(x2lab, 3, at=0, cex=1.1, font=2, line=0.1)
}

intcorr <- function(df=NULL) {
  # Calculates two different methods of showing differences in slopes
  # Args:
  # Returns:
  df   <- data.frame(y=c(0.4375, 0.4167, 0.5313, 0.4516, 0.5417, 0.5172,
                         0.1500, 0.5161, 0.5313, 0.5000, 0.4839, 0.3871,
                         0.3125, 0.5313, 0.4063, 0.5517, 0.3871, 0.7188,
                         0.7188, 0.5484, 0.9375, 0.5938, 0.4375, 0.8750,
                         0.9063, 0.6774, 0.5625, 0.5000, 0.5313),
                    x1=c("B", "B", "B", "B", "A", "A", "A", "B", "A",
                         "A", "B", "A", "B", "A", "A", "B", "B", "B",
                         "A", "B", "A", "A", "B", "B", "A", "A", "A",
                         "B", "A"),
                    x2=c(4.00, 3.63, 3.67, 3.63, 3.57, 3.47, 4.27,
                         2.17, 3.87, 3.60, 3.43, 4.30, 4.13, 4.67,
                         4.13, 3.37, 2.63, 2.33, 3.30, 2.33, 3.57, 3.73,
                         3.50, 3.63, 2.57, 3.43, 3.93, 2.89, 4.23))
  plot(y ~ x2, data=subset(df, x1=="A"))
  abline(lm(y ~ x2, data=subset(df, x1=="A")), lty=3)
  points(y ~ x2, data=subset(df, x1=="B"), pch=19)
  abline(lm(y ~ x2, data=subset(df, x1=="B")))
  l    <- lm(y ~ x1 * x2, data=df)
  summary(l)

  r1   <- sqrt(summary(lm(y ~ x2, data=subset(df, x1=="A")))$r.squared)
  r2   <- sqrt(summary(lm(y ~ x2, data=subset(df, x1=="B")))$r.squared)
  z1   <- (1/2) * (log(1+r1) - log(1-r1))
  n1   <- nrow(subset(df, x1=="A"))
  n2   <- nrow(subset(df, x1=="B"))
  z2   <- (1/2) * (log(1+r2) - log(1-r2))
  z    <- (z1 - z2) / (sqrt((1/(n1 - 3)) + (1/(n2 - 3))))
  pval <- 2 * pnorm(-abs(z))
}

#-----------------------------------------------------------------------
# Main program
# Prepare data frame
#-----------------------------------------------------------------------
df                   <- read.csv('../data/fae-bpd.csv')
df.tab0              <- subset(df, trial==1&run==1&qtypen==1)
df.tab0              <- df.tab0[ , c("id", "age", "gender", "group",
                                     "trait.attributions",
                                     "like.evaluations",
                                     "fmrimovparam",
                                     "beta_12_4248",
                                     "beta12_2050",
                                     "beta60_2830",
                                     "beta_62_2614",
                                     "beta_8436",
                                     "empathy",
                                     "empathy.suffering",
                                     "empathy.positivesharing",
                                     "empathy.responsivecrying",
                                     "empathy.emotionalattention",
                                     "empathy.feelforothers",
                                     "empathy.emotionalcontagion",
                                     "paranoia",
                                     "paranoia.frequency",
                                     "paranoia.conviction",
                                     "paranoia.distress",
                                     "iri.total",
                                     "iri.persp",
                                     "iri.fantasy",
                                     "iri.emp",
                                     "iri.distress",
                                     "need.for.cognition",
                                     "bdi", "stais", "stait")]

#-----------------------------------------------------------------------
# Create data frame for covariates at baseline
#-----------------------------------------------------------------------
covs                 <- c("Age",
                          "Trait.Attributions",
                          "Like.Evaluations",
                          "STAIT", "STAIS", "BDI", "IRI.total",
                          "beta_12_4248",
                          "beta12_2050",
                          "beta60_2830",
                          "beta_62_2614",
                          "beta_8436",
                          "IRI.persp", "IRI.fantasy", "IRI.emp",
                          "IRI.distress", "Need.for.Cognition",
                          "Empathy",
                          "Empathy.Suffering",
                          "Empathy.Positivesharing",
                          "Empathy.Responsivecrying",
                          "Empathy.Emotionalattention",
                          "Empathy.Feelforothers",
                          "Empathy.Emotionalcontagion",
                          "Paranoia",
                          "Paranoia.Frequency",
                          "Paranoia.Conviction",
                          "Paranoia.Distress",
                          "fmrimovparam")

df.tab1              <- data.frame(matrix(ncol=9, nrow=length(covs) + 2))
colnames(df.tab1)    <- c("Char","HCCN","HCCMean","HCCSD","BPDN",
                          "BPDMean","BPDSD","Stat","P")

pdf("../output/figures/fae-bpd_1_overview.pdf")
m                    <- NA
tmp                  <- calcmlm(df=df)
#detach("package:lmerTest", unload=TRUE)

#-----------------------------------------------------------------------
# Add gender and run chi-square test
#-----------------------------------------------------------------------
df.tmp               <- subset(df, group=="HCC"&gender=="M"
                               &trial==1&run==1&qtypen==1)
df.tab1$Char[1]      <- "Males"
df.tab1$HCCN[1]      <- length(df.tmp$id)
m[1]                 <- length(df.tmp$id)

df.tmp               <- subset(df, group=="HCC"&gender=="F"
                               &trial==1&run==1&qtypen==1)
df.tab1$Char[2]      <- "Females"
df.tab1$HCCN[2]      <- length(df.tmp$id)
m[2]                 <- length(df.tmp$id)

df.tmp               <- subset(df, group=="BPD"&gender=="M"
                               &trial==1&run==1&qtypen==1)
df.tab1$Char[1]      <- "Males"
df.tab1$BPDN[1]      <- length(df.tmp$id)
m[3]                 <- length(df.tmp$id)

df.tmp               <- subset(df, group=="BPD"&gender=="F"
                               &trial==1&run==1&qtypen==1)
df.tab1$Char[2]      <- "Females"
df.tab1$BPDN[2]      <- length(df.tmp$id)
m[4]                 <- length(df.tmp$id)

cq                   <- chisq.test(matrix(c(m), 2, 2, byrow=TRUE))
df.tab1$Stat[1]      <- round(cq$statistic, 1)
df.tab1$P[1]         <- round(cq$p.value, 2)
df.tab1$Stat[2]      <- round(cq$statistic, 1)
df.tab1$P[2]         <- round(cq$p.value, 2)
#-----------------------------------------------------------------------
# loop over numeric covariates
#-----------------------------------------------------------------------
for (i in 1:length(covs)) {
  df.tab1$Char[i+2]      <- covs[i] 
  p                      <- subset(df.tab0, group=="HCC")
  df.tab1$HCCN[i+2]      <- length(p[complete.cases(p[ , tolower(covs[i])]),
                                     tolower(covs[i])])

  df.tab1$HCCMean[i+2]   <- round(mean(p[complete.cases(p[ , tolower(covs[i])]),
                                         tolower(covs[i])]), 1)

  df.tab1$HCCSD[i+2]     <- round(sd(p[complete.cases(p[ , tolower(covs[i])]),
                                       tolower(covs[i])]), 1)

  p                      <- subset(df.tab0, group=="BPD")
  df.tab1$BPDN[i+2]      <- length(p[complete.cases(p[ , tolower(covs[i])]),
                                     tolower(covs[i])])

  df.tab1$BPDMean[i+2]   <- round(mean(p[complete.cases(p[ , tolower(covs[i])]),
                                         tolower(covs[i])]), 1)

  df.tab1$BPDSD[i+2]     <- round(sd(p[complete.cases(p[ , tolower(covs[i])]),
                                       tolower(covs[i])]), 1)

  t                      <- t.test(df.tab0[ , tolower(covs[i])] ~ df.tab0$group)
  df.tab1$Stat[i+2]      <- round(t$statistic, 1)
  df.tab1$P[i+2]         <- round(t$p.value, 2)
}
#-----------------------------------------------------------------------
# Adapt header and first column
#-----------------------------------------------------------------------
colnames(df.tab1)          <- c("Char", "HCC", "Mean", "SD",
                                "BPD", "Mean ","SD ","Stat","P")

df.tab1$Char               <- c("Males", "Females", "Age",
                                "Trait.Attributions",
                                "Like.Evaluations",
                                "STAIT", "STAIS", "BDI",
                                "beta_12_4248",
                                "beta12_2050",
                                "beta60_2830",
                                "beta_62_2614",
                                "beta_8436",
                                "IRI.total",
                                "IRI.persp",
                                "IRI.fantasy",
                                "IRI.emp",
                                "IRI.distress",
                                "Need.for.Cognition",
                                "Empathy",
                                "Empathy.Suffering",
                                "Empathy.Positivesharing",
                                "Empathy.Responsivecrying",
                                "Empathy.Emotionalattention",
                                "Empathy.Feelforothers",
                                "Empathy.Emotionalcontagion",
                                "Paranoia",
                                "Paranoia.Frequency",
                                "Paranoia.Conviction",
                                "Paranoia.Distress",
                                "fmrimovparam")

# Correct zeros in p-values
df.tab1$P[df.tab1$P==0]    <- 0.001

# Write the data frame to disk
write.csv(df.tab1, "../output/tables/fae-bpd_1_overview.csv",
          row.names=FALSE, na="")


# Plot regressions
rcovs                      <- c("empathy",
                                "paranoia",
                                "need.for.cognition")
rcovsl                     <- c("Empathy",
                                "Paranoia",
                                "Need for Cognition")

par(mar = c(6, 1, 1, 2), oma = c(9, 6, 2, 6))
lf                         <- layout(matrix(c(1:4), 2, 2, byrow = TRUE))
#layout.show(lf)
#-----------------------------------------------------------------------
# Covariates predicting trait attributions 
#-----------------------------------------------------------------------
for (i in 1:length(rcovs)) { 
  df.tab0.hcc   <- subset(df.tab0, group=="HCC")
  df.tab0.hcc$y <- df.tab0.hcc$trait.attributions * 100
  df.tab0.hcc$x <- df.tab0.hcc[ , rcovs[i]]

  tmp           <- plot(x=df.tab0.hcc[ , rcovs[i]],
                        y=df.tab0.hcc$trait.attributions * 100,
                        xlim=c((min(df.tab0[, rcovs[i]], na.rm=TRUE)),
                               (max(df.tab0[, rcovs[i]], na.rm=TRUE))),
                        ylim=c((min(df.tab0$trait.attributions * 100, na.rm=TRUE) * 0.9),
                               (max(df.tab0$trait.attributions * 100, na.rm=TRUE) * 1.1)),
                        col="white", xaxt="n", yaxt="n", xlab="",
                        ylab="", bty="l")

  reg           <- plotreg(df.tab0.hcc, lwd=2, lty=3, pch=1)
  rhcc          <- ifelse(reg$coef[2]<0, -1, 1) *
                   round(sqrt(summary(reg)$r.squared), 2)

  df.tab0.bpd   <- subset(df.tab0, group=="BPD")
  df.tab0.bpd$y <- df.tab0.bpd$trait.attributions * 100
  df.tab0.bpd$x <- df.tab0.bpd[ , rcovs[i]]
  reg           <- plotreg(df.tab0.bpd, lwd=2, lty=1, pch=19)
  rbpd          <- ifelse(reg$coef[2]<0, -1, 1) *
                   round(sqrt(summary(reg)$r.squared), 2)
  tmp           <- axis(1)
  tmp           <- mtext(rcovsl[i], 1, 3, font=2)
  if (!i == 2) {
    tmp         <- axis(2, las=1)
    tmp         <- mtext("% Trait Attributions", 2, 3, font=2)
  }
  #tmp           <- text(max(dfbquote(paste(r[HCC], "=", .(rhcc))), 3, adj=1, cex=0.3)
  l1            <- bquote(paste(r[HCC], .(ifelse(rhcc>0, " =   ", " = ")),
                                .(rhcc)))
  l2            <- bquote(paste(r[BPD], .(ifelse(rbpd>0, " =   ", " = ")),
                                .(rbpd)))
  #tmp           <- mtext(bquote(paste(r[BPD], "=", .(rbpd))), 3, adj=1, line=-1, cex=0.3)
  #tmp           <- legend("topright", legend=c(as.expression(l1),
  #                                             as.expression(l2)),
  #                        bty="n", cex=0.75, inset=c(0, -0.1), xpd=NA)
  # include analysis of covariance
  ac            <- lm(df.tab0$trait.attributions * 100
                      ~ df.tab0$group * df.tab0[ , rcovs[i]])
  tmp           <- print(anova(ac))
  sac           <- summary(ac)
  if (i == 1) {
      tab3        <- round(sac$coef, 3)
      #tab3.col1   <- rbind(rcovs[i], rep("", 3))
      tab3.col2   <- rbind("Intercept", "Group", rcovs[i],
                             paste("Group *", rcovs[i]))
      tab3        <- cbind(tab3.col2, tab3)
  } else {
      #tab3.col1   <- rbind(rcovs[i], rep("", 3))
      tab3.col2   <- rbind("Intercept", "Group", rcovs[i],
                             paste("Group *", rcovs[i]))
      sac$coef    <- cbind(tab3.col2, round(sac$coef, 3))
      tab3        <- rbind(tab3, sac$coef)
  } 
                        
}
plot.new()
legend("topleft", legend=c("HCC", "BPD"), pch=c(1, 19),
       lty=c(2, 1), bty="n", lwd=2)

# write the coefficients to disk
#tab3[ , c(2:5)]  <- round(tab3[c(2:12), c(2:5)], 3)
write.csv(tab3, "../output/tables/fae-bpd_2_ancova.csv",
          row.names=FALSE, na="")
#-----------------------------------------------------------------------
# Covariates predicting likability ratings 
#-----------------------------------------------------------------------
tmp <- par(mar = c(6, 1, 1, 2), oma = c(9, 6, 2, 6))
lf  <- layout(matrix(c(1:4), 2, 2, byrow = TRUE))
for (i in 1:length(rcovs)) { 
  dfm           <- ddply(subset(df, qtype=="like"),
                          c("id", "group"),
                          summarise, N=sum(!is.na(resp)),
                          mlike=mean(resp, na.rm=TRUE))

  df.tab0       <- merge(df.tab0, dfm, all=TRUE)
  df.tab0.hcc   <- subset(df.tab0, group=="HCC")
  df.tab0.hcc$y <- df.tab0.hcc$mlike
  df.tab0.hcc$x <- df.tab0.hcc[ , rcovs[i]]

  tmp           <- plot(x=df.tab0.hcc[ , rcovs[i]],
                        y=df.tab0.hcc$mlike,
                        xlim=c((min(df.tab0[, rcovs[i]], na.rm=TRUE)),
                        (max(df.tab0[, rcovs[i]], na.rm=TRUE))),
                        ylim=c((min(df.tab0$mlike, na.rm=TRUE) * 0.9),
                        (max(df.tab0$mlike, na.rm=TRUE) * 1.1)),
                        col="white", xaxt="n", yaxt="n", xlab="",
                        ylab="", bty="l")

  reg           <- plotreg(df.tab0.hcc, lwd=2, lty=3, pch=1)
  rhcc          <- ifelse(reg$coef[2]<0, -1, 1) *
                   round(sqrt(summary(reg)$r.squared), 2)

  df.tab0.bpd   <- subset(df.tab0, group=="BPD")
  df.tab0.bpd$y <- df.tab0.bpd$mlike
  df.tab0.bpd$x <- df.tab0.bpd[ , rcovs[i]]
  reg           <- plotreg(df.tab0.bpd, lwd=2, lty=1, pch=19)
  rbpd          <- ifelse(reg$coef[2]<0, -1, 1) *
                   round(sqrt(summary(reg)$r.squared), 2)
  l1            <- bquote(paste(r[HCC], .(ifelse(rhcc>0, " =   ", " = ")),
                                .(rhcc)))
  l2            <- bquote(paste(r[BPD], .(ifelse(rbpd>0, " =   ", " = ")),
                                .(rbpd)))
  #tmp           <- legend("topright", legend=c(as.expression(l1),
  #                                             as.expression(l2)),
  #                        bty="n", cex=0.75, inset=c(0, -0.1), xpd=NA)
  tmp           <- axis(1)
  tmp           <- mtext(rcovsl[i], 1, 3, font=2)
  if (!i == 2) {
    tmp <- axis(2, las=1)
    tmp <- mtext("Like Evaluations", 2, 3, font=2)
  }
  
  # include analysis of covariance
  lac          <- lm(df.tab0$like.evaluations
                     ~ df.tab0$group * df.tab0[ , rcovs[i]])
  tmp          <- print(anova(ac))
  slac         <- summary(lac)
  if (i == 1) {
      tab4        <- round(slac$coef, 3)
      #tab3.col1   <- rbind(rcovs[i], rep("", 3))
      tab4.col2   <- rbind("Intercept", "Group", rcovs[i],
                             paste("Group *", rcovs[i]))
      tab4        <- cbind(tab4.col2, tab4)
  } else {
      #tab3.col1   <- rbind(rcovs[i], rep("", 3))
      tab4.col2   <- rbind("Intercept", "Group", rcovs[i],
                             paste("Group *", rcovs[i]))
      slac$coef   <- cbind(tab4.col2, round(slac$coef, 3))
      tab4        <- rbind(tab4, slac$coef)
  } 
}
tmp <- plot.new()
tmp <- legend("topleft", legend=c("HCC", "BPD"), pch=c(1, 19),
              lty=c(2, 1), bty="n", lwd=2)
## write the coefficients to disk
write.csv(tab4, "../output/tables/fae-bpd_4_ancova.csv",
          row.names=FALSE, na="")
#-----------------------------------------------------------------------
# Covariates predicting betas 
#-----------------------------------------------------------------------
rois           <- c("beta_12_4248","beta12_2050",
                    "beta60_2830","beta_62_2614", "beta_8436")
roinames       <- c("L_Precuneus", "R_MFG", "R_TPJ", "L_STG", "L_dACC") 
par(mar = c(6, 1, 1, 6), oma = c(9, 6, 2, 6))
lf <- layout(matrix(c(1:6), 3, 2, byrow = FALSE))
#layout.show(lf)
for (i in 1:length(rois)) { 
  df.tab0.hcc    <- subset(df.tab0, group=="HCC")
  df.tab0.hcc$y  <- df.tab0.hcc[ , rois[i]]
  df.tab0.hcc$x  <- df.tab0.hcc[ , rcovs[1]]

  plot(x=df.tab0.hcc[ , rcovs[1]],
       y=df.tab0.hcc[ , rois[i]],
       xlim=c((min(df.tab0[, rcovs[1]], na.rm=TRUE) * 1.0),
              (max(df.tab0[, rcovs[1]], na.rm=TRUE) * 1.0)),
       ylim=c((min(df.tab0[, rois[i]] , na.rm=TRUE) * 1.1),
              (max(df.tab0[, rois[i]] , na.rm=TRUE) * 1.1)),
       col="white", xaxt="n", yaxt="n", xlab="", ylab="", bty="l")
  reg            <- plotreg(df.tab0.hcc, lwd=2, lty=3, pch=1)
  rhcc           <- ifelse(reg$coef[2]<0, -1, 1) *
                    round(sqrt(summary(reg)$r.squared), 2)

  df.tab0.bpd    <- subset(df.tab0, group=="BPD")
  df.tab0.bpd$y  <- df.tab0.bpd[ , rois[i]]
  df.tab0.bpd$x  <- df.tab0.bpd[ , rcovs[1]]
  reg            <- plotreg(df.tab0.bpd, lwd=2, lty=1, pch=19)
  rbpd           <- ifelse(reg$coef[2]<0, -1, 1) *
                    round(sqrt(summary(reg)$r.squared), 2)
  l1             <- bquote(paste(r[HCC], .(ifelse(rhcc>0, " =   ", " = ")),
                                 .(rhcc)))
  l2             <- bquote(paste(r[BPD], .(ifelse(rbpd>0, " =   ", " = ")),
                                 .(rbpd)))
  #tmp            <- legend("topright", legend=c(as.expression(l1),
  #                                              as.expression(l2)),
  #                         bty="n", cex=0.75, inset=c(0, -0.3), xpd=NA)
  axis(2, las=1)
  mtext(roinames[i], 3, 1, font=2, cex=0.7)
  mtext(expression(paste(bold(beta), " (a. u.)")), 2, 3, font=2)
  axis(1)
  if (i == 3) {
      mtext(rcovsl[1], 1, 3, font=2)
  }
}
plot.new()
mtext(rcovsl[1], 3, 3, font=2)
legend("topleft", legend=c("HCC", "BPD"), pch=c(1, 19),
       lty=c(2, 1), bty="n", lwd=2)

#-----------------------------------------------------------------------
# Dotplots for attributions and likability ratings
#-----------------------------------------------------------------------
tmp         <- par(mar = c(4, 0, 1, 1.5), oma = c(16, 6, 2, 22))
lf          <- layout(matrix(c(1:4), 2, 2, byrow = TRUE))
sdf         <- subset(df, qtypen==1 &
                  ((expversion==1 & run==2 & (trial==1 | trial==2)) |
                   (expversion==2 & run==1 & (trial==1 | trial==2)) |
                   (expversion==3 & run==2 & (trial==1 | trial==2)) |
                   (expversion==4 & run==1 & (trial==1 | trial==4))))

dfm         <- ddply(sdf, c("valence", "group"), summarise,
                     N=sum(!is.na(trait.attributions.val)),
                     m=mean(trait.attributions.val*100, na.rm=TRUE),
                     sd=sd(trait.attributions.val*100, na.rm=TRUE),
                     eb=sd/sqrt(N))

dfm$group   <- factor(dfm$group)
dfm$valence <- factor(dfm$valence)
tmp         <- plot(c(1:4), dfm$m, col="white", xaxt="n", xlab="",
                    bty="n", las=1, xlim=c(0.8, 4.2), ylim=c(0, 80))
tmp         <- plotdot(data.frame(x=c(1.25, 3.25), y=dfm$m[c(1, 3)],
                                  eb=dfm$eb[c(1.25, 3.25)]), pch=19,
                       lwd=1.1, l=TRUE, lty=1)
tmp         <- plotdot(data.frame(x=c(1.75, 3.75), y=dfm$m[c(2, 4)],
                                  eb=dfm$eb[c(1.75, 3.75)]), pch=1,
                       lwd=1.1, l=TRUE, lty=3)
tmp         <- axis(1, at=c(1:4), labels=c("", "", "", ""),
                    cex.axis=0.7, xpd=NA, tcl=0)
tmp         <- axis(1, at=c(1.5, 3.5), labels=c("Neg", "Pos"),
                    cex.axis=1, xpd=NA)
tmp         <- mtext("% Trait Attributions", 2, 3, font=2)
tmp         <- plot.new()
tmp         <- legend("bottomleft", legend=c("HCC", "BPD"),
                      pch=c(1, 19), lty=c(3, 1),
                      bty="n", cex=0.75, pt.cex=1)
dfm         <- ddply(subset(df, qtype=="like"),
                     c("valence", "id", "ttype", "group"), summarise,
                     N=sum(!is.na(resp)), mm=mean(resp, na.rm=TRUE))

dfm         <- ddply(subset(dfm, ttype=="dispositional"),
                     c("valence", "group"), summarise, N=sum(!is.na(m)),
                     m=mean(mm, na.rm=TRUE), sd=sd(mm, na.rm=TRUE),
                     eb=sd/sqrt(N)) 

dfm$group   <- factor(dfm$group)
dfm$valence <- factor(dfm$valence)
tmp         <- plot(c(1:4), dfm$m, col="white", xaxt="n", xlab="",
                    bty="n", las=1, xlim=c(0.8, 4.2), ylim=c(1, 6))
tmp         <- plotdot(data.frame(x=c(1.25, 3.25), y=dfm$m[c(1, 3)],
                                  eb=dfm$eb[c(1.25, 3.25)]), pch=19,
                       lwd=1.1, l=TRUE, lty=1)
tmp         <- plotdot(data.frame(x=c(1.75, 3.75), y=dfm$m[c(2, 4)],
                                  eb=dfm$eb[c(1.75, 3.75)]), pch=1,
                       lwd=1.1, l=TRUE, lty=3)
tmp         <- axis(1, at=c(1:4), labels=c("", "", "", ""),
                    cex.axis=0.7, xpd=NA, tcl=0)
tmp         <- axis(1, at=c(1.5, 3.5), labels=c("Neg", "Pos"),
                    cex.axis=1, xpd=NA)
tmp         <- mtext("Dispositional", 3, 0.5, cex=0.75, font=2)
tmp         <- mtext("Like Evaluations", 2, 3, font=2)
dfm         <- ddply(subset(df, qtype=="like"),
                     c("valence", "id", "ttype", "group"),
                     summarise, N=sum(!is.na(resp)),
                     mm=mean(resp, na.rm=TRUE))

dfm         <- ddply(subset(dfm, ttype=="situational"),
                     c("valence", "group"),
                     summarise, N=sum(!is.na(m)),
                     m=mean(mm, na.rm=TRUE),
                     sd=sd(mm, na.rm=TRUE),
                     eb=sd/sqrt(N))

dfm$group   <- factor(dfm$group)
dfm$valence <- factor(dfm$valence)
tmp         <- plot(c(1:4), dfm$m, col="white", axes=FALSE, xlab="",
                    bty="n", las=1, xlim=c(0.8, 4.2), ylim=c(1, 6))

tmp         <- plotdot(data.frame(x=c(1.25, 3.25), y=dfm$m[c(1, 3)],
                                  eb=dfm$eb[c(1.25, 3.25)]), pch=19,
                       lwd=1.1, lty=1, l=TRUE)

tmp         <- plotdot(data.frame(x=c(1.75, 3.75), y=dfm$m[c(2, 4)],
                                  eb=dfm$eb[c(1.75, 3.75)]), pch=1,
                       lwd=1.1, lty=3, l=TRUE)

tmp         <- axis(1, at=c(1:4), labels=c("", "", "", ""),
                    cex.axis=0.7, xpd=NA, tcl=0)

tmp         <- axis(1, at=c(1.5, 3.5), labels=c("Neg", "Pos"),
                    cex.axis=1, xpd=NA)

tmp         <- mtext("Situational", 3, 0.5, cex=0.75, font=2)

#-----------------------------------------------------------------------
# Dot plots (forest) for attributions and likability ratings
#-----------------------------------------------------------------------
tmp         <- par(mar = c(4, 0, 1, 1.5), oma = c(16, 6, 2, 22))
lf          <- layout(matrix(c(1:4), 2, 2, byrow = TRUE))
sdf         <- subset(df, qtypen==1 &
                  ((expversion==1 & run==2 & (trial==1 | trial==2)) |
                   (expversion==2 & run==1 & (trial==1 | trial==2)) |
                   (expversion==3 & run==2 & (trial==1 | trial==2)) |
                   (expversion==4 & run==1 & (trial==1 | trial==4))))

dfm         <- ddply(sdf, c("valence", "group"), summarise,
                     N=sum(!is.na(trait.attributions.val)),
                     m=mean(trait.attributions.val*100, na.rm=TRUE),
                     sd=sd(trait.attributions.val*100, na.rm=TRUE),
                     eb=sd/sqrt(N))

dfm$group   <- factor(dfm$group)

dev.off()
