#################################################################################
# R code for Visualizing, modeling and forecasting age distribution of fertility
#################################################################################

# packages <- c("Compositional", "psych", "ftsa", "meboot", "pracma", "reldist", "flexmix", "demography",
#               "MortalityLaws", "HMDHFDplus", "LaplacesDemon", "dplyr", "compositions", "easyCODA", "doMC",
#               "DescTools", "xtable", "DescTools", "fdadensity", "frechet", "reshape2", "ggplot2")

packages_windows <- c("Compositional", "psych", "ftsa", "meboot", "pracma", "reldist", "flexmix", "demography",
                      "MortalityLaws", "HMDHFDplus", "LaplacesDemon", "dplyr", "compositions", "easyCODA", "doParallel",
                      "DescTools", "xtable", "DescTools", "frechet", "fdadensity", "frechet", "reshape2", "ggplot2")


# Install packages not yet installed

# installed_packages <- packages %in% rownames(installed.packages())
# if (any(installed_packages == FALSE)) {
#   install.packages(packages[!installed_packages])
# }

installed_packages_windows <- packages_windows %in% rownames(installed.packages())

if (any(installed_packages_windows == FALSE)) {
  install.packages(packages_windows[!installed_packages_windows])
}

# Packages loading

# invisible(lapply(packages, library, character.only = TRUE))

invisible(lapply(packages_windows, library, character.only = TRUE))

# Function for selecting the number of components via EVR

select_K <- function(tau, eigenvalue)
{
  k_max = length(eigenvalue)
  k_all = rep(0, k_max-1)
  for(k in 1:(k_max-1))
  {
    k_all[k] = (eigenvalue[k+1]/eigenvalue[k])*ifelse(eigenvalue[k]/eigenvalue[1] > tau, 1, 0) + ifelse(eigenvalue[k]/eigenvalue[1] < tau, 1, 0)
  }
  K_hat = which.min(k_all)
  return(K_hat)
}

# Line-width constants (pt)
CURVE_LWD <- 1     # data curves
AXIS_LWD  <- 0.5   # axes / tick marks

# Apply before every figure block.  Relies on the output device being
# opened with pointsize = 12 so that cex = 1 gives 12 pt text.
# mar: bottom, left, top, right
# Top margin is larger because the horizontal y-axis label sits above the plot.
# Left margin is smaller because there is no rotated y-axis label.
std_par <- function(mar = c(4, 3, 4, 1)) {
  par(
    cex       = 1,
    cex.axis  = 1,
    cex.lab   = 1,
    cex.main  = 1,
    lwd       = AXIS_LWD,
    bty       = "n",
    mar       = mar
  )
}

# Add y-axis label as horizontal text sitting just above the top of the
# y-axis (left-aligned, same x as the axis).  Call AFTER the plot is drawn.
add_ylab <- function(label) {
  usr <- par("usr")
  text(x = usr[1], y = usr[4],
       labels = label, adj = c(0, 0),
       xpd = TRUE, cex = 1)
}

# Add x-axis label that terminates (right-justified) at the right end of
# the x-axis, sitting just below the axis tick labels.
# Call AFTER the plot is drawn.
add_xlab <- function(label, line = 2.8) {
  mtext(label, side = 1, las = 1, adj = 1,
        at = par("usr")[2], line = line, cex = 1)
}

# Draw a text label + arrow pointing to a curve (replaces legend entries).
# (x0, y0) : anchor of the text label
# (x1, y1) : arrow tip on the curve
# side      : "left" or "right" text justification relative to the arrow
curve_label <- function(label, x0, y0, x1, y1,
                        x_text = 0, y_text = 0,
                        adj = c(0.5, 0.5),
                        text_col = "black",
                        arrow_col = text_col, ...) {
  text(x0+x_text, y0+y_text, labels = label, cex = 1, 
       adj = adj, col = text_col, ...)
  arrows(x0 = x0, y0 = y0, x1 = x1, y1 = y1,
         length = 0.07, lwd = AXIS_LWD, col = arrow_col)
}

savefig = function (filename, height=10, width = (1 + sqrt(5))/2*height, type=c("eps","pdf","jpg","png"), pointsize = 10, family = "Helvetica", sublines = 0, toplines = 0, leftlines = 0, res=300) 
{
  type <- match.arg(type)
  filename <- paste(filename, ".", type, sep = "")
  if(type=="eps")
  {
    postscript(file = filename, horizontal = FALSE, 
               width = width/2.54, height = height/2.54, pointsize = pointsize, 
               family = family, onefile = FALSE, print.it = FALSE)
  }
  else if(type=="pdf")
  {
    pdf(file = filename, width=width/2.54, height=height/2.54, pointsize=pointsize,
        family=family, onefile=TRUE)
  }
  else if(type=="jpg")
  {
    jpeg(filename=filename, width=width, height=height, res=res,quality=100, units="cm")#, pointsize=pointsize*50)
  }
  else if(type=="png")
  {
    png(filename=filename, width=width, height=height, res=res, units="cm")#, pointsize=pointsize*50)
  }
  else
    stop("Unknown file type")
  par(mgp = c(2.2, 0.45, 0), tcl = -0.4, mar = c(3.2 + sublines + 0.25 * (sublines > 0), 
                                                 3.5 + leftlines, 1 + toplines, 1) + 0.1)
  par(pch = 1)
  invisible()
}

savepdf = function(...)
{
  savefig(...,type="pdf")
}


################################
# Import the ABS fertility data
################################

state = c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT")
ages = 15:49

# 2020:2024

AUS_fert_2020_2024 = read.csv("ABS_fert_2020_2024.csv", skip = 5, header = TRUE)
AUS_fert_collect_2020_2024 = matrix(NA, (8 * 35), 5)
AUS_national_fert_collect_2020_2024 = matrix(NA, 35, 5)
for(ik in 15:49)
{
  index = which(AUS_fert_2020_2024[,1]==paste("Age: ", ik, sep =""))
  print(ik); print(index)
  AUS_fert_collect_2020_2024[((ik-15)*8+1):((ik-14)*8),] = as.matrix(AUS_fert_2020_2024[(index+2):(index+9),3:7])
  AUS_national_fert_collect_2020_2024[(ik-14),] = as.matrix(AUS_fert_2020_2024[(index+1),3:7])
  rm(index)
}
colnames(AUS_fert_collect_2020_2024) = colnames(AUS_national_fert_collect_2020_2024) = 2020:2024
rownames(AUS_fert_collect_2020_2024) = rep(c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"), 35)
rownames(AUS_national_fert_collect_2020_2024) = ages

# 2017:2019

AUS_fert_2017_2019 = read.csv("ABS_fert_2017_2019.csv", skip = 5, header = TRUE)
AUS_fert_collect_2017_2019 = matrix(NA, (8*35), 3)
AUS_national_fert_collect_2017_2019 = matrix(NA, 35, 3)
for(ik in 15:49)
{
  index = which(AUS_fert_2017_2019[,1]==paste("Age: ", ik, sep =""))
  print(ik); print(index)
  AUS_fert_collect_2017_2019[((ik-15)*8+1):((ik-14)*8),] = as.matrix(AUS_fert_2017_2019[(index+2):(index+9),3:5])
  AUS_national_fert_collect_2017_2019[(ik-14),] = as.matrix(AUS_fert_2017_2019[(index+1),3:5])
  rm(index)
}
colnames(AUS_fert_collect_2017_2019) = colnames(AUS_national_fert_collect_2017_2019) = 2017:2019
rownames(AUS_fert_collect_2017_2019) = rep(c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"), 35)
rownames(AUS_national_fert_collect_2017_2019) = ages

# 2014:2016

AUS_fert_2014_2016 = read.csv("ABS_fert_2014_2016.csv", skip = 5, header = TRUE)
AUS_fert_collect_2014_2016 = matrix(NA, (8*35), 3)
AUS_national_fert_collect_2014_2016 = matrix(NA, 35, 3)
for(ik in 15:49)
{
  index = which(AUS_fert_2014_2016[,1]==paste("Age: ", ik, sep =""))
  print(ik); print(index)
  AUS_fert_collect_2014_2016[((ik-15)*8+1):((ik-14)*8),] = as.matrix(AUS_fert_2014_2016[(index+2):(index+9),3:5])
  AUS_national_fert_collect_2014_2016[(ik-14),] = as.matrix(AUS_fert_2014_2016[(index+1),3:5])
  rm(index)
}
colnames(AUS_fert_collect_2014_2016) = colnames(AUS_national_fert_collect_2014_2016) = 2014:2016
rownames(AUS_fert_collect_2014_2016) = rep(c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"), 35)
rownames(AUS_national_fert_collect_2014_2016) = ages

# 2000:2013

AUS_fert_2000_2013 = read.csv("ABS_fert_2000_2013.csv", skip = 5, header = TRUE)
AUS_fert_collect_2000_2013 = matrix(NA, (8*35), 14)
AUS_national_fert_collect_2000_2013 = matrix(NA, 35, 14)
for(ik in 15:49)
{
  index = which(AUS_fert_2000_2013[,1]==paste("Age: ", ik, sep =""))
  print(ik); print(index)
  AUS_fert_collect_2000_2013[((ik-15)*8+1):((ik-14)*8),] = as.matrix(AUS_fert_2000_2013[(index+2):(index+9),3:16])
  AUS_national_fert_collect_2000_2013[(ik-14),] = as.matrix(AUS_fert_2000_2013[(index+1),3:16])
  rm(index)
}
colnames(AUS_fert_collect_2000_2013) = colnames(AUS_national_fert_collect_2000_2013) = 2000:2013
rownames(AUS_fert_collect_2000_2013) = rep(c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"), 35)
rownames(AUS_national_fert_collect_2000_2013) = ages

# 1990:1999

AUS_fert_1990_1999 = read.csv("ABS_fert_1990_1999.csv", skip = 5, header = TRUE)
AUS_fert_collect_1990_1999 = matrix(NA, (8*35), 10)
AUS_national_fert_collect_1990_1999 = matrix(NA, 35, 10)
for(ik in 15:49)
{
  index = which(AUS_fert_1990_1999[,1]==paste("Age: ", ik, sep =""))
  print(ik); print(index)
  AUS_fert_collect_1990_1999[((ik-15)*8+1):((ik-14)*8),] = as.matrix(AUS_fert_1990_1999[(index+2):(index+9),3:12])
  AUS_national_fert_collect_1990_1999[(ik-14),] = as.matrix(AUS_fert_1990_1999[(index+1),3:12])
  rm(index)
}
colnames(AUS_fert_collect_1990_1999) = colnames(AUS_national_fert_collect_1990_1999) = 1990:1999
rownames(AUS_fert_collect_1990_1999) = rep(c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"), 35)
rownames(AUS_national_fert_collect_1990_1999) = ages

# 1975:1989

AUS_fert_1975_1989 = read.csv("ABS_fert_1975_1989.csv", skip = 5, header = TRUE)
AUS_fert_collect_1975_1989 = matrix(NA, (8*35), 15)
AUS_national_fert_collect_1975_1989 = matrix(NA, 35, 15)
for(ik in 15:49)
{
  index = which(AUS_fert_1975_1989[,1]==paste("Age: ", ik, sep =""))
  print(ik); print(index)
  AUS_fert_collect_1975_1989[((ik-15)*8+1):((ik-14)*8),] = as.matrix(AUS_fert_1975_1989[(index+2):(index+9),3:17])
  AUS_national_fert_collect_1975_1989[(ik-14),] = as.matrix(AUS_fert_1975_1989[(index+1),3:17])
  rm(index)
}
colnames(AUS_fert_collect_1975_1989) = colnames(AUS_national_fert_collect_1975_1989) = 1975:1989
rownames(AUS_fert_collect_1975_1989) = rep(c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT"), 35)
rownames(AUS_national_fert_collect_1975_1989) = ages

# together

AUS_national_fert_all = cbind(AUS_national_fert_collect_1975_1989, AUS_national_fert_collect_1990_1999, AUS_national_fert_collect_2000_2013,
                              AUS_national_fert_collect_2014_2016, AUS_national_fert_collect_2017_2019, AUS_national_fert_collect_2020_2024)

AUS_fert_all = cbind(AUS_fert_collect_1975_1989, AUS_fert_collect_1990_1999, AUS_fert_collect_2000_2013, 
                     AUS_fert_collect_2014_2016, AUS_fert_collect_2017_2019, AUS_fert_collect_2020_2024)
colnames(AUS_fert_all) = colnames(AUS_national_fert_all) = 1975:2024
age = 15:49
n_age = length(age)

year <- years <- 1975:2024
n_year <- length(year)  # 50

AUS_fert_NSW = AUS_fert_all[seq(1, 280, by = 8), ]
AUS_fert_VIC = AUS_fert_all[seq(2, 280, by = 8), ]
AUS_fert_QLD = AUS_fert_all[seq(3, 280, by = 8), ]
AUS_fert_SA  = AUS_fert_all[seq(4, 280, by = 8), ]
AUS_fert_WA  = AUS_fert_all[seq(5, 280, by = 8), ]
AUS_fert_TAS = AUS_fert_all[seq(6, 280, by = 8), ]
AUS_fert_NT  = AUS_fert_all[seq(7, 280, by = 8), ]
AUS_fert_ACT = AUS_fert_all[seq(8, 280, by = 8), ]

colnames(AUS_fert_NSW) = colnames(AUS_fert_VIC) = colnames(AUS_fert_QLD) = colnames(AUS_fert_SA) = 
  colnames(AUS_fert_WA) = colnames(AUS_fert_TAS) = colnames(AUS_fert_NT) = colnames(AUS_fert_ACT) = year

rownames(AUS_fert_NSW) = rownames(AUS_fert_VIC) = rownames(AUS_fert_QLD) = rownames(AUS_fert_SA) = 
  rownames(AUS_fert_WA) = rownames(AUS_fert_TAS) = rownames(AUS_fert_NT) = rownames(AUS_fert_ACT) = age

# Skewness of the age distribution: third standardised moment of f_{x,t}
# over age. NOT skewness(b), which treats the 35 counts as 35 observations.
skew_age <- function(b) {
  f <- b / sum(b); mu <- sum(age * f)
  s2 <- sum((age - mu)^2 * f)
  sum((age - mu)^3 * f) / s2^1.5
}

NSW_fert_skewness = apply(AUS_fert_NSW, 2, skew_age)
VIC_fert_skewness = apply(AUS_fert_VIC, 2, skew_age)
QLD_fert_skewness = apply(AUS_fert_QLD, 2, skew_age)
SA_fert_skewness  = apply(AUS_fert_SA,  2, skew_age)
WA_fert_skewness  = apply(AUS_fert_WA,  2, skew_age)
TAS_fert_skewness = apply(AUS_fert_TAS, 2, skew_age)
NT_fert_skewness  = apply(AUS_fert_NT,  2, skew_age)
ACT_fert_skewness = apply(AUS_fert_ACT, 2, skew_age)

# fts objects

AUS_fert_national_fts = fts(age, AUS_national_fert_all)
AUS_fert_NSW_fts = fts(age, AUS_fert_NSW)
AUS_fert_VIC_fts = fts(age, AUS_fert_VIC)
AUS_fert_QLD_fts = fts(age, AUS_fert_QLD)
AUS_fert_SA_fts  = fts(age, AUS_fert_SA)
AUS_fert_WA_fts  = fts(age, AUS_fert_WA)
AUS_fert_TAS_fts = fts(age, AUS_fert_TAS)
AUS_fert_NT_fts  = fts(age, AUS_fert_NT)
AUS_fert_ACT_fts = fts(age, AUS_fert_ACT)

######################
# normalise densities
######################

# using fdadensity package

AUS_den_normalise = t(normaliseDensities(dmatrix = t(AUS_national_fert_all), dSup = age))
NSW_den_normalise = t(normaliseDensities(dmatrix = t(AUS_fert_NSW), dSup = age))
VIC_den_normalise = t(normaliseDensities(dmatrix = t(AUS_fert_VIC), dSup = age))
QLD_den_normalise = t(normaliseDensities(dmatrix = t(AUS_fert_QLD), dSup = age))
SA_den_normalise  = t(normaliseDensities(dmatrix = t(AUS_fert_SA),  dSup = age))
WA_den_normalise  = t(normaliseDensities(dmatrix = t(AUS_fert_WA),  dSup = age))
TAS_den_normalise = t(normaliseDensities(dmatrix = t(AUS_fert_TAS), dSup = age))
NT_den_normalise  = t(normaliseDensities(dmatrix = t(AUS_fert_NT),  dSup = age))
ACT_den_normalise = t(normaliseDensities(dmatrix = t(AUS_fert_ACT), dSup = age))

# self normalize

AUS_den = sweep(AUS_national_fert_all, 2, colSums(AUS_national_fert_all), "/")
NSW_den = sweep(AUS_fert_NSW, 2, colSums(AUS_fert_NSW), "/")
VIC_den = sweep(AUS_fert_VIC, 2, colSums(AUS_fert_VIC), "/")
QLD_den = sweep(AUS_fert_QLD, 2, colSums(AUS_fert_QLD), "/")
SA_den  = sweep(AUS_fert_SA,  2, colSums(AUS_fert_SA),  "/")
WA_den  = sweep(AUS_fert_WA,  2, colSums(AUS_fert_WA),  "/")
TAS_den = sweep(AUS_fert_TAS, 2, colSums(AUS_fert_TAS), "/")
NT_den  = sweep(AUS_fert_NT,  2, colSums(AUS_fert_NT),  "/")
ACT_den = sweep(AUS_fert_ACT, 2, colSums(AUS_fert_ACT), "/")


#####################################################
# Gini coefficient for measuring spread distribution
#####################################################


Gini_NSW = Gini_VIC = Gini_QLD = Gini_SA = 
  Gini_WA = Gini_TAS = Gini_NT = Gini_ACT = 
  sd_NSW = sd_VIC = sd_QLD = sd_SA = 
  sd_WA = sd_TAS = sd_NT = sd_ACT =   vector("numeric", 50)
for(ik in 1:50)
{
  # Gini
  
  Gini_NSW[ik] = Gini(AUS_fert_NSW[,ik])
  Gini_VIC[ik] = Gini(AUS_fert_VIC[,ik])
  Gini_QLD[ik] = Gini(AUS_fert_QLD[,ik])
  Gini_SA[ik]  = Gini(AUS_fert_SA[,ik])
  Gini_WA[ik]  = Gini(AUS_fert_WA[,ik])
  Gini_TAS[ik] = Gini(AUS_fert_TAS[,ik])
  Gini_NT[ik]  = Gini(AUS_fert_NT[,ik])
  Gini_ACT[ik] = Gini(AUS_fert_ACT[,ik])
  
  # sd
  
  sd_NSW[ik] = sd(AUS_fert_NSW[,ik])
  sd_VIC[ik] = sd(AUS_fert_VIC[,ik])
  sd_QLD[ik] = sd(AUS_fert_QLD[,ik])
  sd_SA[ik]  = sd(AUS_fert_SA[,ik])
  sd_WA[ik]  = sd(AUS_fert_WA[,ik])
  sd_TAS[ik] = sd(AUS_fert_TAS[,ik])
  sd_NT[ik]  = sd(AUS_fert_NT[,ik])
  sd_ACT[ik] = sd(AUS_fert_ACT[,ik])
  print(ik); rm(ik)
}

#######
# mode
#######

NSW_den_mode = age[apply(NSW_den, 2, which.max)] # 32
VIC_den_mode = age[apply(VIC_den, 2, which.max)] # 32
QLD_den_mode = age[apply(QLD_den, 2, which.max)] # 31
SA_den_mode  = age[apply(SA_den,  2, which.max)] # 32
WA_den_mode  = age[apply(WA_den,  2, which.max)] # 32
TAS_den_mode = age[apply(TAS_den, 2, which.max)] # 31
NT_den_mode  = age[apply(NT_den,  2, which.max)] # 29
ACT_den_mode = age[apply(ACT_den, 2, which.max)] # 33

state_den_mode = rbind(NSW_den_mode, VIC_den_mode, QLD_den_mode, SA_den_mode, WA_den_mode, 
                       TAS_den_mode, NT_den_mode, ACT_den_mode)
rownames(state_den_mode) = state
colnames(state_den_mode) = year

t(state_den_mode[,c("1975", "1980", "1985", "1990", "1995", "2000",
                    "2005", "2010", "2015", "2020")])


# number of times where one state has the largest mode age out of 50 years

state_den_mode_count = state_den_least_count = list()
for(ik in 1:n_year)
{
  state_den_mode_count[[ik]]  = which(state_den_mode[,ik] == max(state_den_mode[,ik]))
  state_den_least_count[[ik]] = which(state_den_mode[,ik] == min(state_den_mode[,ik]))
  print(ik); rm(ik)    
}

table(unlist(state_den_mode_count))  # 10 23  5 15 11  4  6 38 
table(unlist(state_den_least_count)) #  9  2 12  7 10 25 34  3 

#######################
# Wasserstein distance
#######################
###############################
# compute Wasserstein distance
###############################

# NSW

NSW_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  NSW_den_list = NSW_den_list_2 = list()
  NSW_den_list$y = NSW_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    NSW_den_list_2$y = NSW_den_normalise[,ij]
    NSW_den_list$x = NSW_den_list_2$x = age
    NSW_W2_dist[ik,ij] = dist4den(d1 = NSW_den_list, d2 = NSW_den_list_2)
    rm(ij)
  }
  rm(NSW_den_list); rm(NSW_den_list_2); rm(ik)
}

NSW_den_median = which.min(colSums(NSW_W2_dist, na.rm=TRUE)) # 26
year[which.min(colSums(NSW_W2_dist, na.rm=TRUE))]            # 2000

savefig("Fig_3", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, NSW_den), xlab = "", ylab = "", main = "NSW", lwd = CURVE_LWD)
lines(age, NSW_den[,NSW_den_median], col = "black", lty = 2, lwd = CURVE_LWD)
add_ylab("Density")
curve_label("Median year", x0 = 45, y0 = 0.065,
            x1 = age[which.max(NSW_den[,NSW_den_median])],
            y1 = max(NSW_den[,NSW_den_median]))
dev.off()

# VIC

VIC_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  VIC_den_list = VIC_den_list_2 = list()
  VIC_den_list$y = VIC_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    VIC_den_list_2$y = VIC_den_normalise[,ij]
    VIC_den_list$x = VIC_den_list_2$x = age
    VIC_W2_dist[ik,ij] = dist4den(VIC_den_list, VIC_den_list_2)
    rm(ij)
  }
  rm(VIC_den_list); rm(VIC_den_list_2); rm(ik)
}

VIC_den_median = which.min(colSums(VIC_W2_dist, na.rm=TRUE)) # 25
year[which.min(colSums(VIC_W2_dist, na.rm=TRUE))]            # 1999

savefig("Fig_3_VIC", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, VIC_den), xlab = "", ylab = "", main = "VIC", lwd = CURVE_LWD)
lines(age, VIC_den[,VIC_den_median], col = "black", lty = 2, lwd = CURVE_LWD)
add_ylab("Density")
curve_label("Median year", x0 = 45, y0 = 0.065,
            x1 = age[which.max(VIC_den[,VIC_den_median])],
            y1 = max(VIC_den[,VIC_den_median]))
dev.off()


# QLD

QLD_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  QLD_den_list = QLD_den_list_2 = list()
  QLD_den_list$y = QLD_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    QLD_den_list_2$y = QLD_den_normalise[,ij]
    QLD_den_list$x = QLD_den_list_2$x = age
    QLD_W2_dist[ik,ij] = dist4den(QLD_den_list, QLD_den_list_2)
    rm(ij)
  }
  rm(QLD_den_list); rm(QLD_den_list_2); rm(ik)
}

QLD_den_median = which.min(colSums(QLD_W2_dist, na.rm=TRUE)) # 25
year[which.min(colSums(QLD_W2_dist, na.rm=TRUE))]            # 1999

savefig("Fig_3_QLD", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, QLD_den), xlab = "", ylab = "", main = "QLD", lwd = CURVE_LWD)
lines(age, QLD_den[,QLD_den_median], col = "black", lty = 2, lwd = CURVE_LWD)
add_ylab("Density")
add_xlab("Age")
curve_label("Median year", x0 = 45, y0 = 0.065,
            x1 = age[which.max(QLD_den[,QLD_den_median])],
            y1 = max(QLD_den[,QLD_den_median]))
dev.off()

# SA

SA_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  SA_den_list = SA_den_list_2 = list()
  SA_den_list$y = SA_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    SA_den_list_2$y = SA_den_normalise[,ij]
    SA_den_list$x = SA_den_list_2$x = age
    SA_W2_dist[ik,ij] = dist4den(SA_den_list, SA_den_list_2)
    rm(ij)
  }
  rm(SA_den_list); rm(SA_den_list_2); rm(ik)
}

SA_den_median = which.min(colSums(SA_W2_dist, na.rm=TRUE)) # 24
year[which.min(colSums(SA_W2_dist, na.rm=TRUE))]           # 1998

# WA

WA_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  WA_den_list = WA_den_list_2 = list()
  WA_den_list$y = WA_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    WA_den_list_2$y = WA_den_normalise[,ij]
    WA_den_list$x = WA_den_list_2$x = age
    WA_W2_dist[ik,ij] = dist4den(WA_den_list, WA_den_list_2)
    rm(ij)
  }
  rm(WA_den_list); rm(WA_den_list_2); rm(ik)
}

WA_den_median = which.min(colSums(WA_W2_dist, na.rm=TRUE)) # 25
year[which.min(colSums(WA_W2_dist, na.rm=TRUE))]           # 1999

# TAS

TAS_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  TAS_den_list = TAS_den_list_2 = list()
  TAS_den_list$y = TAS_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    TAS_den_list_2$y = TAS_den_normalise[,ij]
    TAS_den_list$x = TAS_den_list_2$x = age
    TAS_W2_dist[ik,ij] = dist4den(TAS_den_list, TAS_den_list_2)
    rm(ij)
  }
  rm(TAS_den_list); rm(TAS_den_list_2); rm(ik)
}

TAS_den_median = which.min(colSums(TAS_W2_dist, na.rm=TRUE)) # 25
year[which.min(colSums(TAS_W2_dist, na.rm=TRUE))]            # 1999

# NT

NT_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  NT_den_list = NT_den_list_2 = list()
  NT_den_list$y = NT_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    NT_den_list_2$y = NT_den_normalise[,ij]
    NT_den_list$x = NT_den_list_2$x = age
    NT_W2_dist[ik,ij] = dist4den(NT_den_list, NT_den_list_2)
    rm(ij)
  }
  rm(NT_den_list); rm(NT_den_list_2); rm(ik)
}

NT_den_median = which.min(colSums(NT_W2_dist, na.rm=TRUE)) # 24
year[which.min(colSums(NT_W2_dist, na.rm=TRUE))]           # 1998

# ACT

ACT_W2_dist = matrix(NA, n_year, n_year)
for(ik in 1:n_year)
{
  ACT_den_list = ACT_den_list_2 = list()
  ACT_den_list$y = ACT_den_normalise[,ik]
  index = (1:n_year)[-ik]
  for(ij in index)
  {
    ACT_den_list_2$y = ACT_den_normalise[,ij]
    ACT_den_list$x = ACT_den_list_2$x = age
    ACT_W2_dist[ik,ij] = dist4den(ACT_den_list, ACT_den_list_2)
    rm(ij)
  }
  rm(ACT_den_list); rm(ACT_den_list_2); rm(ik)
}

ACT_den_median = which.min(colSums(ACT_W2_dist, na.rm=TRUE)) # 27
year[which.min(colSums(ACT_W2_dist, na.rm=TRUE))]            # 2001

##############
# median year
##############

year[c(NSW_den_median, VIC_den_median, QLD_den_median, SA_den_median,
       WA_den_median, TAS_den_median, NT_den_median, ACT_den_median)]

# (original)                                 1998 1998 1999 1997 1997 1999 1998 2001
# (updated after including 2023 & 2024 data) 2000 1999 1999 1998 1999 1999 1998 2001


##################
# Continuous Mode
##################

# Kannisto's (2001) quadratic interpolation, as used for the modal age at
# death. Let x be the age with the most births and
#
#   A = d(x) - d(x-1),   Bq = d(x) - d(x+1),   both > 0 at an interior maximum.
#
# Fitting a quadratic through ages x-1, x, x+1 and taking its vertex gives
#
#   M = x + A/(A + Bq),
#
# written below in the equivalent symmetric second-difference form. ABS births
# by single year of age are counts for women aged x last birthday, i.e. the
# interval [x, x+1), so the count is placed at the midpoint of that interval;
# this keeps M within [x, x+1], so it never falls below the integer mode of
# Table 1. Vectorised over the columns of a 35 x k matrix of birth counts.


yrc <- as.character(year)
fert_array  <- array(NA_real_, c(8, n_age, n_year), dimnames = list(state, age, yrc))
for (s in state)
{
  fert_array[s, , ] <- get(paste0("AUS_fert_", s))
}

pmode <- function(Y) {
  m   <- max.col(t(Y), ties.method = "first")
  out <- ages[m] + 0.5                     
  int <- m > 1 & m < length(ages)          
  if (any(int)) {
    w  <- which(int)
    bp <- Y[cbind(m[int] - 1, w)]
    bm <- Y[cbind(m[int],     w)]
    bn <- Y[cbind(m[int] + 1, w)]
    d  <- bp - 2 * bm + bn                
    del <- ifelse(d < 0, (bp - bn) / (2 * d), 0)
    out[int] <- ages[m[int]] + 0.5 + pmax(pmin(del, 0.5), -0.5)
  }
  out
}

M_integer <- apply(fert_array, c(1, 3), function(b){age[which.max(b)]})    # integer, as Table 1
M_continuous <- apply(fert_array, c(1, 3), function(b){pmode(matrix(b, ncol = 1))})
N <- apply(fert_array, c(1, 3), sum)

# Check the refined estimate must lie inside the modal age interval [x, x+1]
stopifnot(all(M_continuous >= M_integer - 1e-9),
          all(M_continuous <= M_integer + 1 + 1e-9))


cat("\n===== TABLE 1 AS PUBLISHED (integer mode) =====\n")
year_tb1 <- as.character(seq(1975, 2020, by = 5))
print(t(M_integer[, year_tb1]))
raw <- setNames(rep(0, 8), state)
for (t in yrc) {
  top <- M_integer[, t] == max(M_integer[, t])
  raw[top] <- raw[top] + 1
}
cat("\n'# time of being mode' row (all ties counted):\n"); print(raw)
cat("sum =", sum(raw), "->", round(sum(raw) / n_year, 2),
    "states share the maximum in an average year\n")



#### Multinomial bootstrap for uncertainty computation ####
# For state $s$, year $t$, resample $n{s,t}$ births from Multinomial(n{s,t}, f_{x,t}), recompute continuous M_continuous, repeat B times.

NB <- 2000          # bootstrap replicates
set.seed(20260820)

boot_array <- array(NA_real_, c(8, n_year, NB), dimnames = list(state, yrc, NULL))
for (s in 1:8) for (t in 1:n_year) 
{
  b <- fert_array[s, , t]
  boot_array[s, t, ] <- pmode(rmultinom(NB, sum(b), b / sum(b)))
}


CI <- apply(boot_array, c(1, 2), quantile, c(0.025, 0.975))
W  <- CI[2, , ] - CI[1, , ]
dimnames(W) <- dimnames(M_continuous)

cat("\n===== PRECISION OF THE MODAL AGE =====\n")
prec <- data.frame(state = state, mean_births = round(rowMeans(N)),
                   CI_halfwidth = round(rowMeans(W) / 2, 2))
print(prec[order(prec$CI_halfwidth), ], row.names = FALSE)

cat("\n===== MODAL AGE WITH 95% BOOTSTRAP CIs =====\n")
for (t in year_tb1) {
  cat("\n", t, "\n", sep = "")
  print(data.frame(state = state, mode_int = M_integer[, t],
                   mode_cont = sprintf("%.2f", M_continuous[, t]),
                   CI95 = sprintf("[%.2f, %.2f]", CI[1, , t], CI[2, , t])),
        row.names = FALSE)
}

###  Pairwise Contrasts ###

# States are independent samples, so the difference is resampled directly.
sig <- matrix("", n_year, 7, dimnames = list(yrc, state[1:7]))
for (t in seq_len(n_year)) for (s in 1:7) {
  q <- quantile(boot_array[8, t, ] - boot_array[s, t, ], c(0.025, 0.975), names = FALSE)
  sig[t, s] <- if (q[1] > 0) "+" else if (q[2] < 0) "-" else "."
}
cat("\n===== ACT vs EACH STATE, PER YEAR =====\n")
cat("'+' ACT higher, '-' ACT lower, '.' indistinguishable at 5%\n\n")
print(sig[year_tb1, ])
cat("\nYears (of", n_year, ") in which ACT is significantly higher:\n")
print(colSums(sig == "+"))
cat("\nYears (of", n_year, ") indistinguishable:\n")
print(colSums(sig == "."))

cat("\n----- the NT 1995 cell -----\n")
for (s in c("NSW", "VIC", "SA", "ACT")) {
  d <- boot_array[s, "1995", ] - boot_array["NT", "1995", ]
  q <- quantile(d, c(0.025, 0.975), names = FALSE)
  cat(sprintf("  %-4s - NT: %+5.2f  95%% CI [%+.2f, %+.2f]  %s\n", s, mean(d), q[1], q[2],
              if (q[1] > 0) "distinguishable" else "NOT distinguishable"))
}
cat(sprintf("  NT 1995 mode %.2f  CI [%.2f, %.2f]  on %d births\n",
            M_continuous["NT", "1995"], CI[1, "NT", "1995"], CI[2, "NT", "1995"],
            sum(fert_array["NT", , "1995"])))


### expected years holding the highest mode ###
# Replaces the tie-inflated count; sums to n_year by construction.
winner  <- apply(boot_array, c(2, 3), which.max)
cnt_rep <- t(apply(winner, 2, function(w) tabulate(w, 8)))     # NB x 8
exp_cnt <- setNames(colMeans(cnt_rep), state)
ci_cnt  <- apply(cnt_rep, 2, quantile, c(0.025, 0.975))

cat("\n===== EXPECTED NUMBER OF YEARS WITH THE HIGHEST MODAL AGE =====\n")
tab <- data.frame(state = state, 
                  paper_count_ties = raw,
                  expected_years = round(exp_cnt, 1),
                  CI95 = sprintf("[%d, %d]", ci_cnt[1, ], ci_cnt[2, ]))
print(tab[order(-tab$expected_years), ], row.names = FALSE)
cat("expected total =", round(sum(exp_cnt), 1), "\n")

lead <- cnt_rep[, 8] - apply(cnt_rep[, -8], 1, max)
cat(sprintf("ACT margin over the next-best state: %.1f  95%% CI [%d, %d]\n",
            mean(lead), quantile(lead, .025, names = FALSE),
            quantile(lead, .975, names = FALSE)))
cat(sprintf("Pr(ACT holds the most years) = %.3f\n", mean(lead > 0)))

### Temporal shift ###
cat("\n===== MODAL AGE SHIFT 1975 -> 2024 =====\n")
for (s in state) {
  d <- boot_array[s, "2024", ] - boot_array[s, "1975", ]
  q <- quantile(d, c(0.025, 0.975), names = FALSE)
  cat(sprintf("  %-4s %.2f -> %.2f   %+5.2f yr  95%% CI [%+.2f, %+.2f]\n",
              s, M_continuous[s, "1975"], M_continuous[s, "2024"], mean(d), q[1], q[2]))
}
rng <- apply(M_integer, 2, function(x) diff(range(x)))
cat("\nWithin-year range across states: mean", round(mean(rng), 2),
    "| excluding NT:", round(mean(apply(M_integer[state != "NT", ], 2,
                                        function(x) diff(range(x)))), 2), "\n")
#############
# Gini Index
#############

dens <- function(b){b / sum(b)}
gini <- function(b) {
  f <- sort(dens(b)); n <- length(f)
  sum((2 * seq_len(n) - n - 1) * f) / (n^2 * mean(f))
}

# standard deviation of the age distribution, in years of age
sd_age <- function(b) {
  f <- dens(b); mu <- sum(age * f); sqrt(sum((age - mu)^2 * f))
}

# interquartile range, interpolated from the cumulative distribution function
iqr_age <- function(b) {
  f <- dens(b); cdf <- cumsum(f)
  q <- approx(c(0, cdf), c(age[1], age + 1), xout = c(0.25, 0.75),
              ties = "ordered")$y
  q[2] - q[1]
}

# Skewness of the age distribution: the third standardised moment of f_{x,t}
# over age. 
skew_age <- function(b) {
  f <- dens(b); mu <- sum(age * f)
  s2 <- sum((age - mu)^2 * f)
  sum((age - mu)^3 * f) / s2^1.5
}

apply_measure <- function(fun) apply(fert_array, c(1, 3), fun)
Gini <- apply_measure(gini)
SDv  <- apply_measure(sd_age)
IQRv <- apply_measure(iqr_age)
Skew <- apply_measure(skew_age)



################
# clr functions
################


# fdata: year x age
# ncomp_selection: EVR of K = 6
# fh: forecast horizon
# fore_method: forecasting method
# object_interest: point or interval forecasts
# B: number of bootstrap replication
# alpha: level of significance

clr_fun <- function(fdata, ncomp_selection, fh, fore_method, object_interest, B = 399, alpha)
{
  n_age = ncol(fdata)
  n_year = nrow(fdata)
  
  if(any(fdata == 0))
  {
    fdata = replace(fdata, which(fdata == 0), 10^-15)
  }
  h_x_t = CLR(fdata)$LR
  
  SVD_decomp = svd(h_x_t)
  if(ncomp_selection == "EVR")
  {
    ncomp = select_K(tau = 0.001, eigenvalue = SVD_decomp$d^2)
  }
  else if(ncomp_selection == "fixed")
  {
    ncomp = 6
  }
  else
  {
    warning("The number of retained component must be chosen by EVR or fixed at 6.")
  }
  basis = SVD_decomp$v[,1:ncomp]
  score = t(basis) %*% t(h_x_t)
  recon = basis %*% score
  resi = t(h_x_t) - recon
  
  # reconstruction (model in-sample fitting)
  
  recon = invCLR(t(recon))
  
  if(object_interest == "point")
  {
    # forecasts of principal component scores
    
    score_fore = matrix(NA, ncomp, fh)
    for(ik in 1:ncomp)
    {
      if(fore_method == "RWF_no_drift")
      {
        score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = FALSE)$mean
      }
      else if(fore_method == "RWF_drift")
      {
        score_fore[ik,] = rwf(as.numeric(score[ik,]), h = fh, drift = TRUE)$mean
      }
      else if(fore_method == "ETS")
      {
        score_fore[ik,] = forecast(ets(as.numeric(score[ik,])), h = fh)$mean
      }
      else if(fore_method == "ARIMA")
      {
        score_fore[ik,] = forecast(auto.arima(as.numeric(score[ik,])), h = fh)$mean
      }
      else
      {
        warning("Univariate time series forecasting method is not on the list.")
      }
    }
    
    # obtain forecasts in real-valued space
    
    fore_val = basis %*% score_fore
    fore_count = invCLR(t(fore_val))
    return(list(ncomp = ncomp, fore_count = fore_count))
  }
  else if(object_interest == "interval")
  {
    # determine in-sample forecast error for principal component scores
    
    olivia = matrix(NA, ncomp, fh)
    if(fore_method == "ets")
    {
      for(ij in 1:ncomp)
      {
        olivia[ij,] = forecast(ets(score[ij,]), h = fh)$mean
      }
    }
    else if(fore_method == "arima")
    {
      for(ij in 1:ncomp)
      {
        olivia[ij,] = forecast(auto.arima(score[ij,]), h = fh)$mean
      }
    }
    else if(fore_method == "rwf")
    {
      for(ij in 1:ncomp)
      {
        olivia[ij,] = rwf(score[ij,], h = fh, drift = TRUE)$mean
      }
    }
    else if(fore_method == "rw")
    {
      for(ij in 1:ncomp)
      {
        olivia[ij,] = rwf(score[ij,], h = fh, drift = FALSE)$mean
      }
    }
    else
    {
      warning("Forecasting method must be ets, arima, rwf or rw.")
    }
    forerr = matrix(NA, (n_year - ncomp - fh + 1), ncomp)
    for(i in fh:(n_year - ncomp))
    {
      k = i + (ncomp - fh)
      fore = matrix(NA, 1, ncomp)
      if(fore_method == "ets")
      {
        for(j in 1:ncomp)
        {
          fore[,j] = forecast(ets(score[j,1:k]), h = fh)$mean[fh]
        }
      }
      else if(fore_method == "arima")
      {
        for(j in 1:ncomp)
        {
          fore[,j] = forecast(auto.arima(score[j,1:k]), h = fh)$mean[fh]
        }
      }
      else if(fore_method == "rwf")
      {
        if(k <= 2)
        {
          for(j in 1:ncomp)
          {
            fore[,j] = score[j,k]
          }
        }
        if(k > 2)
        {
          for(j in 1:ncomp)
          {
            fore[,j] = rwf(score[j,1:k], h = fh, drift = TRUE)$mean[fh]
          }
        }
      }
      else if(fore_method == "rw")
      {
        if(k == 1)
        {
          for(j in 1:ncomp)
          {
            fore[,j] = score[j,1]
          }
        }
        if(k > 1)
        {
          for(j in 1:ncomp)
          {
            fore[,j] = rwf(score[j,1:k], h = fh, drift = FALSE)$mean[fh]
          }
        }
      }
      forerr[i - fh + 1,] = score[, k + fh] - fore
    }
    # bootstrapping residuals
    K = 1
    
    q = array(NA, dim = c(n_age, B, K, fh))
    for(j in 1:fh)
    {
      for(i in 1:n_age)
      {
        for(k in 1:K)
        {
          q[i,,k,j] = sample(resi[i,], size = B, replace = TRUE)
        }
      }
    }
    rm(i); rm(j); rm(k)
    # bootstrapping PC score errors
    ny = array(NA, dim = c(ncomp, B, fh))
    for(j in 1:fh)
    {
      for(i in 1:ncomp)
      {
        ny[i,,j] = sample(forerr[,i], size = B, replace = TRUE)
      }
    }
    rm(i); rm(j)
    # adding the PC score error to the predicted score
    oli = array(rep(olivia, B * fh), dim = c(ncomp, B, fh))
    fo = array(NA, dim = c(ncomp, B, fh))
    for(j in 1:fh)
    {
      for(i in 1:B)
      {
        fo[,i,j] = oli[,i,j] + ny[,i,j]
      }
    }
    rm(i); rm(j)
    # construct bootstrapped samples
    pred = array(NA, dim = c(n_age, B, K, fh))
    for(j in 1:fh)
    {
      for(i in 1:B)
      {
        for(k in 1:K)
        {
          pred[,i,k,j] = basis %*% fo[,i,j] + q[,i,k,j]
        }
      }
    }
    rm(i); rm(j); rm(k)
    
    pred_resize = array(NA, dim = c(n_age, B * K, fh))
    for(j in 1:fh)
    {
      for(i in 1:B)
      {
        pred_resize[, (((i-1)*K+1):(i*K)), ] = pred[,i,,j]
      }
    }
    rm(i); rm(j)
    
    # transform back
    
    d_x_t_star_fore = array(NA, dim = c(n_age, B * K, fh))
    for(iw in 1:fh)
    {
      for(ij in 1:(B * K))
      {
        d_x_t_star_fore[,ij,iw] = invCLR(t(pred_resize[,ij,iw]))
      }
    }
    rm(iw); rm(ij)
    return(list(PI = apply(d_x_t_star_fore, c(1, 3), quantile, c(alpha/2, 1-alpha/2)), ncomp = ncomp))
  }
  else
  {
    warning("object_interest must either be point or interval.")
  }
}


# fdata: n by p data matrix
# ncomp_selection: method for selecting the number of retained components
# fh: forecast horizon
# fore_method: forecasting method
# object_interest: point or interval forecasts
# B = 399: number of bootstrap replications
# alpha: level of significance

##################
# Point forecasts
##################

AUS_den_fore = clr_fun(fdata = t(AUS_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                       object_interest = "point", B = 399, alpha = 0.2)$fore_count

NSW_den_fore = clr_fun(fdata = t(NSW_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                       object_interest = "point", B = 399, alpha = 0.2)$fore_count

VIC_den_fore = clr_fun(fdata = t(VIC_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                       object_interest = "point", B = 399, alpha = 0.2)$fore_count

QLD_den_fore = clr_fun(fdata = t(QLD_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                       object_interest = "point", B = 399, alpha = 0.2)$fore_count

SA_den_fore = clr_fun(fdata = t(SA_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                      object_interest = "point", B = 399, alpha = 0.2)$fore_count

WA_den_fore = clr_fun(fdata = t(WA_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                      object_interest = "point", B = 399, alpha = 0.2)$fore_count

TAS_den_fore = clr_fun(fdata = t(TAS_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                       object_interest = "point", B = 399, alpha = 0.2)$fore_count

NT_den_fore = clr_fun(fdata = t(NT_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                      object_interest = "point", B = 399, alpha = 0.2)$fore_count

ACT_den_fore = clr_fun(fdata = t(ACT_den), ncomp_selection = "fixed", fh = 10, fore_method = "ARIMA", 
                       object_interest = "point", B = 399, alpha = 0.2)$fore_count

rownames(NSW_den_fore) = rownames(VIC_den_fore) = rownames(QLD_den_fore) = rownames(SA_den_fore) = 
  rownames(WA_den_fore) = rownames(TAS_den_fore) = rownames(NT_den_fore) = rownames(ACT_den_fore) = 1:10

colnames(NSW_den_fore) = colnames(VIC_den_fore) = colnames(QLD_den_fore) = colnames(SA_den_fore) = 
  colnames(WA_den_fore) = colnames(TAS_den_fore) = colnames(NT_den_fore) = colnames(ACT_den_fore) = age


overall_range = range(range(NSW_den), range(VIC_den), range(QLD_den), range(SA_den),
                      range(WA_den), range(TAS_den), range(NT_den), range(ACT_den))

overall_range_fore = range(range(NSW_den_fore), range(VIC_den_fore), range(QLD_den_fore), range(SA_den_fore),
                           range(WA_den_fore), range(TAS_den_fore), range(NT_den_fore), range(ACT_den_fore))


# data: n by p data matrix
# fh: forecast horizon
# fmethod: univariate time-series forecasting method
# ncomp_selection: EVR or K=6

compositions_clr_transformation <- function(data, fh, fmethod, ncomp_selection)
{
  n_age = ncol(data)
  clr_data = matrix(clr(data), nrow(data), ncol(data))
  
  SVD_decomp = svd(clr_data)
  if(ncomp_selection == "EVR")
  {
    ncomp = select_K(0.001, eigenvalue = SVD_decomp$d^2)
  }
  else if(ncomp_selection == "provide")
  {
    ncomp = 6
  }
  basis = SVD_decomp$v[,1:ncomp]
  score = t(basis) %*% t(clr_data)
  
  score_fore = matrix(NA, ncomp, 1)
  for(ik in 1:ncomp)
  {
    score_fore[ik,] = forecast(ets(as.numeric(score[ik,])), h = fh)$mean[fh]
  }
  fore_val = basis %*% score_fore
  
  clr_inv_data = as.numeric(clrInv(as.numeric(fore_val))) * 10^5
  return(clr_inv_data)
}

easyCoDa_clr_transformation <- function(data, fh, fmethod, ncomp_selection)
{
  n_age = ncol(data)
  if(any(data == 0))
  {
    data_set = replace(data, which(data == 0), 10^-5)
  }
  else
  {
    data_set = data
  }
  h_x_t = CLR(data_set)$LR
  
  SVD_decomp = svd(h_x_t)
  if(ncomp_selection == "EVR")
  {
    ncomp = select_K(tau = 0.001, eigenvalue = SVD_decomp$d^2)
  }
  else if(ncomp_selection == "provide")
  {
    ncomp = 6
  }
  basis = SVD_decomp$v[,1:ncomp]
  score = t(basis) %*% t(h_x_t)
  
  score_fore = matrix(NA, ncomp, 1)
  for(ik in 1:ncomp)
  {
    score_fore[ik,] = forecast(ets(as.numeric(score[ik,])), h = fh)$mean[fh]
  }
  fore_val = basis %*% score_fore
  clr_inv_data = invCLR(t(fore_val)) * 10^5
  return(clr_inv_data)
}

######
# CDF
######

# data: n by p data matrix
# fh: forecast horizon
# fmethod: univariate time series method
# ncomp_method: way for selecting the number of components

cdf_transformation <- function(data_mat, fh, fmethod, ncomp_selection)
{
  if(any(rowSums(data_mat) > 1))
  {
    data = data_mat/10^5
  }
  else
  {
    data = data_mat
  }
  n_age = ncol(data)
  data_cumsum_dum = matrix(NA, nrow(data), ncol(data))
  for(ij in 1:nrow(data))
  {
    data_cumsum_dum[ij,] = cumsum(data[ij,])
    rm(ij)
  }
  
  # check if any cumsum values equal to 0
  
  if(any(data_cumsum_dum == 0))
  {
    data_cumsum = replace(data_cumsum_dum, which(data_cumsum_dum == 0), 10^-5)
  }
  else
  {
    data_cumsum = data_cumsum_dum
  }
  rm(data_cumsum_dum)
  
  # logit transformation
  
  data_cumsum_logit = matrix(NA, nrow(data), (ncol(data) - 1))
  for(ij in 1:nrow(data))
  {
    data_cumsum_logit[ij,] = logit(data_cumsum[ij, 1:(ncol(data) - 1)])
    rm(ij)
  }
  rm(data_cumsum)
  
  # fitting a functional time series forecasting method
  
  if(ncomp_selection == "EVR")
  {
    ncomp = select_K(tau = 10^-3, eigenvalue = (svd(data_cumsum_logit)$d)^2)
  }
  else if(ncomp_selection == "provide")
  {
    ncomp = 6
  }
  data_cumsum_logit_fore = forecast(ftsm(fts(1:(n_age - 1), t(data_cumsum_logit)), order = ncomp), 
                                    h = fh, method = fmethod)
  
  data_cumsum_logit_fore_add = c(invlogit(data_cumsum_logit_fore$mean$y[,fh]), 1)
  data_cumsum_logit_fore_add_diff = c(data_cumsum_logit_fore_add[1], diff(data_cumsum_logit_fore_add))
  return(data_cumsum_logit_fore_add_diff * 10^5)
}












































