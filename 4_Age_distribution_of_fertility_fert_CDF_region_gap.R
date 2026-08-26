#####################
# clr transformation
#####################

# data: 
# ncomp_method: EVR or K = 6
# fh: forecast horizon
# fmethod: forecasting method
# object_interest: point or interval forecasts
# no_boot: number of bootstrap replications
# alpha: level of significance

AUS_den_sum_one = sweep(AUS_den, 2, colSums(AUS_den), "/")
NSW_den_sum_one = sweep(NSW_den, 2, colSums(NSW_den), "/")
VIC_den_sum_one = sweep(VIC_den, 2, colSums(VIC_den), "/")
QLD_den_sum_one = sweep(QLD_den, 2, colSums(QLD_den), "/")
SA_den_sum_one  = sweep(SA_den,  2, colSums(SA_den), "/")
WA_den_sum_one  = sweep(WA_den,  2, colSums(WA_den), "/")
TAS_den_sum_one = sweep(TAS_den, 2, colSums(TAS_den), "/")
NT_den_sum_one  = sweep(NT_den,  2, colSums(NT_den), "/")
ACT_den_sum_one = sweep(ACT_den, 2, colSums(ACT_den), "/")

# CDF

AUS_CDF = apply(AUS_den_sum_one, 2, cumsum)
NSW_CDF = apply(NSW_den_sum_one, 2, cumsum)
VIC_CDF = apply(VIC_den_sum_one, 2, cumsum)
QLD_CDF = apply(QLD_den_sum_one, 2, cumsum)
SA_CDF  = apply(SA_den_sum_one,  2, cumsum)
WA_CDF  = apply(WA_den_sum_one,  2, cumsum)
TAS_CDF = apply(TAS_den_sum_one, 2, cumsum)
NT_CDF  = apply(NT_den_sum_one, 2, cumsum)
ACT_CDF = apply(ACT_den_sum_one, 2, cumsum)

# gap between NSW and VIC or NSW and ACT

G_F_NSW_VIC = NSW_CDF - VIC_CDF
G_F_NSW_ACT = NSW_CDF - ACT_CDF

plot(fts(age, G_F_NSW_VIC), xlab = "Age", ylab = "Cumulative relative births")

# integral measure

S_F_NSW_VIC = apply(G_F_NSW_VIC, 2, function(diff) trapz(15:49, diff))
S_F_NSW_ACT = apply(G_F_NSW_ACT, 2, function(diff) trapz(15:49, diff))

#######################
# visualization ggplot
#######################

df_G_F_NSW_VIC <- melt(G_F_NSW_VIC)
colnames(df_G_F_NSW_VIC) <- c("Year", "Age", "Difference")

savefig("Fig_4a", type = "png", width = 15, height = 15, toplines = 0.8, pointsize = 12)
ggplot(df_G_F_NSW_VIC, aes(x = Age, y = Year, fill = Difference)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "CDF (NSW) − CDF (VIC)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5)) +
  labs(title = "Time Stochastic Dominance Map: NSW vs. VIC",
       x = "Year", y = "Age") + 
  theme(
    # Move X-axis title to the far right end
    axis.title.x = element_text(hjust = 1),
    # Move Y-axis title to the top end
    axis.title.y = element_text(angle = 0, hjust = 1, vjust = 1.02),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

##############################
# time-series difference plot
##############################

savefig("Fig_4b", type = "png", width = 15, height = 15, toplines = 0.8, pointsize = 12)
std_par()
plot(years, S_F_NSW_VIC, type = "l", col = "blue", lwd = CURVE_LWD,
     xlab = "", ylab = "",
     main = "Evolution of NSW vs. VIC Dominance", ylim = c(0, 0.8),
     xlim = c(1970, 2030))
abline(h = 0, lty = 2, lwd = AXIS_LWD)
add_ylab(expression("Dominance Measure " *italic(S)[italic(t)]))
add_xlab("Year")
dev.off()


fore_national_cdf <- function(data, ncomp_method, fh, fmethod, object_interest, 
                              no_boot = 1000, alpha = 0.2)
{
  data_mod = sweep(data, 1, rowSums(data), "/")
  data_cumsum_dum = matrix(NA, nrow(data), ncol(data))
  for(ij in 1:nrow(data))
  {
    data_cumsum_dum[ij,] = cumsum(data_mod[ij,])
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
  
  # log transformation
  
  data_cumsum_log = matrix(NA, nrow(data), (ncol(data) - 1))
  for(ij in 1:nrow(data))
  {
    data_cumsum_log[ij,] = log(data_cumsum[ij, 1:(ncol(data) - 1)])
    rm(ij)
  }
  rm(data_cumsum)
  rownames(data_cumsum_log) = years[1:nrow(data)]
  
  # fitting a functional time series forecasting method
  
  if(ncomp_method == "EVR")
  {
    ncomp = select_K(tau = 10^-3, eigenvalue = (svd(data_cumsum_logit)$d)^2)
  }
  else if(ncomp_method == "provide")
  {
    ncomp = 6
  }
  else
  {
    warning("The number of components is required.")
  }
  data_cumsum_log_fore = forecast(ftsm(fts(ages[1:(length(ages)-1)], t(data_cumsum_log)), order = ncomp), h = fh,
                                  method = fmethod, pimethod = "nonparametric", level = (1 - alpha) * 100, B = no_boot)
  
  if(object_interest == "point")
  {
    # h-step-ahead forecast
    
    data_cumsum_log_fore_add = rbind(exp(data_cumsum_log_fore$mean$y), rep(1, fh))
    data_cumsum_log_fore_add_diff = rbind(data_cumsum_log_fore_add[1,], apply(data_cumsum_log_fore_add, 2, diff))
    return(data_cumsum_log_fore_add_diff)
  }
  else if(object_interest == "interval")
  {
    # bootstrap forecasts
    
    data_cumsum_logit_fore_add_boot_diff = matrix(NA, ncol(data), no_boot)
    for(ik in 1:no_boot)
    {
      data_cumsum_logit_fore_add_boot = c(invlogit(data_cumsum_logit_fore$bootsamp[,ik,fh]), 1)
      data_cumsum_logit_fore_add_boot_diff[,ik] = c(data_cumsum_logit_fore_add_boot[1],
                                                    diff(data_cumsum_logit_fore_add_boot))
      rm(ik); rm(data_cumsum_logit_fore_add_boot)
    }
    return(t(apply(data_cumsum_logit_fore_add_boot_diff, 1, quantile, c(alpha/2, 1 - alpha/2))))
  }
  else
  {
    warning("object_interest must either be point or interval.")
  }
}

AUS_den_fore_CDF = fore_national_cdf(data = t(AUS_den), ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
NSW_den_fore_CDF = fore_national_cdf(data = t(NSW_den), ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
VIC_den_fore_CDF = fore_national_cdf(data = t(VIC_den), ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
QLD_den_fore_CDF = fore_national_cdf(data = t(QLD_den), ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
SA_den_fore_CDF  = fore_national_cdf(data = t(SA_den),  ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
WA_den_fore_CDF  = fore_national_cdf(data = t(WA_den),  ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
TAS_den_fore_CDF = fore_national_cdf(data = t(TAS_den), ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
NT_den_fore_CDF  = fore_national_cdf(data = t(NT_den),  ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")
ACT_den_fore_CDF = fore_national_cdf(data = t(ACT_den), ncomp_method = "provide", fh = 10, fmethod = "arima", object_interest = "point")

##############
# plot graphs
##############

savefig("Fig_5a", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, AUS_den), col = "gray", main = "AUS", ylab = "", xlab = "",
     ylim = c(0, 0.1), lwd = CURVE_LWD)
lines(fts(ages, AUS_den_fore_CDF), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
curve_label("Historical\ndata", y_text = 0.003, x0 = 20, y0 = 0.082,
            x1 = 25, y1 = 0.055, text_col = "gray")
curve_label("10-steps-ahead\nforecasts", x0 = 42, y0 = 0.082, y_text = 0.003,
            x_text = 0.5, x1 = 36, y1 = AUS_den_fore_CDF[22,5], text_col = "red")
add_xlab("Age")
dev.off()

savefig("Fig_5b", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, NSW_den), col = "gray", main = "NSW", ylab = "", xlab = "",
     ylim = c(0, 0.1), lwd = CURVE_LWD)
lines(fts(ages, NSW_den_fore_CDF), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
# curve_label("Historical\ndata",  x0 = 20, y0 = 0.082, x1 = 25, y1 = 0.055, col = "gray")
# curve_label("10-step\nforecasts",  x0 = 42, y0 = 0.082, x1 = 32, y1 = 0.060, col = "red")
add_xlab("Age")
dev.off()

savefig("Fig_5c", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, VIC_den), col = "gray", main = "VIC", ylab = "", xlab = "",
     ylim = c(0, 0.1), lwd = CURVE_LWD)
lines(fts(ages, VIC_den_fore_CDF), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
add_xlab("Age")
# curve_label("Historical\ndata",    x0 = 20, y0 = 0.082, x1 = 25, y1 = 0.055, col = "gray")
# curve_label("10-step\nforecasts",  x0 = 42, y0 = 0.082, x1 = 32, y1 = 0.060, col = "red")
add_xlab("Age")
dev.off()

savefig("Fig_5d", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, QLD_den), col = "gray", main = "QLD", ylab = "", xlab = "",
     ylim = c(0, 0.1), lwd = CURVE_LWD)
lines(fts(ages, QLD_den_fore_CDF), lwd = CURVE_LWD)
add_ylab("Age distribution\nof birth counts")
add_xlab("Age")
# curve_label("Historical\ndata",    x0 = 20, y0 = 0.082, x1 = 25, y1 = 0.055, col = "gray")
# curve_label("10-step\nforecasts",  x0 = 42, y0 = 0.082, x1 = 32, y1 = 0.060, col = "red")
dev.off()




plot(fts(ages, SA_den), col = "gray", ylab = "Age distribution of birth counts", xlab = "Age", ylim = c(0, 0.1))
lines(fts(ages, SA_den_fore_CDF))

plot(fts(ages, WA_den), col = "gray", ylab = "Age distribution of birth counts", xlab = "Age", ylim = c(0, 0.1))
lines(fts(ages, WA_den_fore_CDF))

plot(fts(ages, TAS_den), col = "gray", ylab = "Age distribution of birth counts", xlab = "Age", ylim = c(0, 0.1))
lines(fts(ages, TAS_den_fore_CDF))

plot(fts(ages, NT_den), col = "gray", ylab = "Age distribution of birth counts", xlab = "Age", ylim = c(0, 0.1))
lines(fts(ages, NT_den_fore_CDF))

plot(fts(ages, ACT_den), col = "gray", ylab = "Age distribution of birth counts", xlab = "Age", ylim = c(0, 0.1))
lines(fts(ages, ACT_den_fore_CDF))



#############
# region gap
#############

######################
#### STATE vs AUS ####
######################

G_F_NSW_AUS = NSW_CDF - AUS_CDF
G_F_VIC_AUS = VIC_CDF - AUS_CDF
G_F_QLD_AUS = QLD_CDF - AUS_CDF
G_F_SA_AUS  = SA_CDF - AUS_CDF
G_F_WA_AUS  = WA_CDF - AUS_CDF
G_F_TAS_AUS = TAS_CDF - AUS_CDF
G_F_NT_AUS  = NT_CDF - AUS_CDF
G_F_ACT_AUS = ACT_CDF - AUS_CDF

#############################
#### forecast region gap ####
#############################

# NSW

NSW_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_NSW_AUS))), h = 10)$mean$y)

savefig("Fig_6a", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_NSW_AUS), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "CDF difference: NSW vs. AUS",
     ylim = c(-0.025, 0.005))
add_ylab("Region gap")
add_xlab("Age")
dev.off()

savefig("Fig_6b", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, NSW_AUS_region_gap_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "Forecast from 2025 to 2034",
     ylim = c(-0.025, 0.005))
add_ylab("Region gap")
add_xlab("Age")
dev.off()

# VIC

VIC_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_VIC_AUS))), method = "ets", h = 10)$mean$y)

VIC_AUS_region_gap_fore_etsna <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_VIC_AUS))), method = "ets.na", h = 10)$mean$y)

savefig("Fig_6c", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_VIC_AUS), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "CDF difference: VIC vs. AUS",
     ylim = c(-0.08, 0.001))
add_ylab("Region gap")
add_xlab("Age")
dev.off()

savefig("Fig_6d", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, VIC_AUS_region_gap_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "Forecast from 2025 to 2034",
     ylim = c(-0.08, 0.001))
add_ylab("Region gap")
add_xlab("Age")
dev.off()


savefig("Fig_6d_etsna", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, VIC_AUS_region_gap_fore_etsna), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "Forecast from 2025 to 2034",
     ylim = c(-0.08, 0.001))
add_ylab("Region gap")
add_xlab("Age")
dev.off()

# QLD

QLD_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_QLD_AUS))), h = 10)$mean$y)

savefig("Fig_6e", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_QLD_AUS), xlab = "", ylab = "",  yaxt="none", lwd = CURVE_LWD,
     main = "CDF difference: QLD vs. AUS",
     ylim = c(0, 0.09))
axis(2, seq(0, 0.09, 0.015))
add_ylab("Region gap")
add_xlab("Age")
dev.off()

savefig("Fig_6f", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, QLD_AUS_region_gap_fore), xlab = "", ylab = "",  yaxt="none", lwd = CURVE_LWD,
     main = "Forecast from 2025 to 2034",
     ylim = c(0, 0.09))
axis(2, seq(0, 0.09, 0.015))
add_ylab("Region gap")
add_xlab("Age")
dev.off()

# SA

SA_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_SA_AUS))), h = 10)$mean$y)
## The forecasts are identical for various h

savefig("Fig_8a_SA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_SA_AUS), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "CDF difference: SA vs. AUS",
     ylim = range(range(G_F_SA_AUS), range(SA_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

savefig("Fig_8b_SA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, SA_AUS_region_gap_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "SA vs. AUS forecast (2025 to 2034)",
     ylim = range(range(G_F_SA_AUS), range(SA_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

# WA

WA_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_WA_AUS))), h = 10)$mean$y)
## The forecasts are identical for various h

savefig("Fig_8a_WA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_WA_AUS), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "CDF difference: WA vs. AUS",
     ylim = range(range(G_F_WA_AUS), range(WA_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

savefig("Fig_8b_WA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, WA_AUS_region_gap_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "WA vs. AUS forecast (2025 to 2034)",
     ylim = range(range(G_F_WA_AUS), range(WA_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

# TAS

TAS_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_TAS_AUS))), h = 10)$mean$y)

savefig("Fig_8a_TAS", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_TAS_AUS), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "CDF difference: TAS vs. AUS",
     ylim = range(range(G_F_TAS_AUS), range(TAS_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

savefig("Fig_8b_TAS", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, TAS_AUS_region_gap_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "TAS vs. AUS forecast (2025 to 2034)",
     ylim = range(range(G_F_TAS_AUS), range(TAS_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

# NT

NT_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_NT_AUS))), h = 10)$mean$y)

savefig("Fig_8a_NT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_NT_AUS), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "CDF difference: NT vs. AUS",
     ylim = range(range(G_F_NT_AUS), range(NT_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

savefig("Fig_8b_NT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, NT_AUS_region_gap_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "NT vs. AUS forecast (2025 to 2034)",
     ylim = range(range(G_F_NT_AUS), range(NT_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

# ACT

ACT_AUS_region_gap_fore <- FisherZInv(forecast(ftsm(fts(age, FisherZ(G_F_ACT_AUS))), h = 10)$mean$y)
## The forecasts are identical for various h

savefig("Fig_8a_ACT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, G_F_ACT_AUS), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "CDF difference: ACT vs. AUS",
     ylim = range(range(G_F_ACT_AUS), range(ACT_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

savefig("Fig_8b_ACT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(age, ACT_AUS_region_gap_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "ACT vs. AUS forecast (2025 to 2034)",
     ylim = range(range(G_F_ACT_AUS), range(ACT_AUS_region_gap_fore)))
add_ylab("Region gap, in CDF units")
add_xlab("Age")
dev.off()

################################
#### forecast national data ####
################################

# AUS

AUS_CDF_fore <- rbind(invlogit(forecast(ftsm(fts(age[1:34], logit(AUS_CDF[1:34,]))), h = 10)$mean$y), rep(1, 10))
NSW_CDF_fore <- NSW_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_7a", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, AUS_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "AUS CDF forecasts (2025 to 2034)")
add_ylab("Probability")
add_xlab("Age")
dev.off()

# NSW

savefig("Fig_7b", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, NSW_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "NSW CDF forecasts (2025 to 2034)")
add_ylab("Probability")
add_xlab("Age")
dev.off()

# VIC

VIC_CDF_fore <- VIC_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_7c", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, VIC_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "VIC CDF forecasts (2025 to 2034)")
add_ylab("Probability")
add_xlab("Age")
dev.off()

# QLD

QLD_CDF_fore <- QLD_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_7d", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, QLD_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "QLD CDF forecasts (2025 to 2034)")
add_ylab("Probability")
add_xlab("Age")
dev.off()

# SA

SA_CDF_fore <- SA_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_8d_SA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, SA_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "SA CDF (2025 to 2034)")
add_ylab("Cumulative distribution function")
add_xlab("Age")
dev.off()

# WA

WA_CDF_fore <- WA_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_8d_WA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, WA_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "WA CDF (2025 to 2034)")
add_ylab("Cumulative distribution function")
add_xlab("Age")
dev.off()

# TAS

TAS_CDF_fore <- TAS_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_8d_TAS", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, TAS_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "TAS CDF (2025 to 2034)")
add_ylab("Cumulative distribution function")
add_xlab("Age")
dev.off()

# NT

NT_CDF_fore <- NT_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_8d_NT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, NT_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "NT CDF (2025 to 2034)")
add_ylab("Cumulative distribution function")
add_xlab("Age")
dev.off()

# ACT

ACT_CDF_fore <- ACT_AUS_region_gap_fore + AUS_CDF_fore

savefig("Fig_8d_ACT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, ACT_CDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD, main = "ACT CDF (2025 to 2034)")
add_ylab("Cumulative distribution function")
add_xlab("Age")
dev.off()

#########################
#### from CDF to PDF ####
#########################

NSW_PDF_fore = VIC_PDF_fore = 
  QLD_PDF_fore = SA_PDF_fore = 
  WA_PDF_fore = TAS_PDF_fore = 
  NT_PDF_fore = ACT_PDF_fore = 
  AUS_PDF_fore = matrix(NA, 35, 10)

for(ik in 1:10)
{
  # AUS
  AUS_PDF_fore[,ik] = c(AUS_CDF_fore[1,ik], diff(AUS_CDF_fore[,ik]))
  
  # NSW
  NSW_PDF_fore[,ik] = c(NSW_CDF_fore[1,ik], diff(NSW_CDF_fore[,ik]))
    
  # VIC
  VIC_PDF_fore[,ik] = c(VIC_CDF_fore[1,ik], diff(VIC_CDF_fore[,ik]))
  
  # QLD
  QLD_PDF_fore[,ik] = c(QLD_CDF_fore[1,ik], diff(QLD_CDF_fore[,ik]))
  
  # SA
  SA_PDF_fore[,ik] = c(SA_CDF_fore[1,ik], diff(SA_CDF_fore[,ik]))
  
  # WA
  WA_PDF_fore[,ik] = c(WA_CDF_fore[1,ik], diff(WA_CDF_fore[,ik]))
  
  # TAS
  TAS_PDF_fore[,ik] = c(TAS_CDF_fore[1,ik], diff(TAS_CDF_fore[,ik]))
  
  # NT
  NT_PDF_fore[,ik] = c(NT_CDF_fore[1,ik], diff(NT_CDF_fore[,ik]))
  
  # ACT
  ACT_PDF_fore[,ik] = c(ACT_CDF_fore[1,ik], diff(ACT_CDF_fore[,ik]))
  
  rm(ik)
}

rownames(NSW_PDF_fore) = rownames(VIC_PDF_fore) = 
  rownames(QLD_PDF_fore) = rownames(SA_PDF_fore) = 
  rownames(WA_PDF_fore) = rownames(TAS_PDF_fore) = 
  rownames(NT_PDF_fore) = rownames(ACT_PDF_fore) = 
  rownames(AUS_PDF_fore) = ages

colnames(NSW_PDF_fore) = colnames(VIC_PDF_fore) = 
  colnames(QLD_PDF_fore) = colnames(SA_PDF_fore) = 
  colnames(WA_PDF_fore) = colnames(TAS_PDF_fore) =
  colnames(NT_PDF_fore) = colnames(ACT_PDF_fore) =
  colnames(AUS_PDF_fore) = 1:10

# AUS

savefig("Fig_8a", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, AUS_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "AUS PDF forecasts (2025 to 2034)",
     ylim = c(0,0.09), yaxt="none")
axis(2, seq(0, 0.09, 0.015))
add_ylab("Density")
add_xlab("Age")
dev.off()

# NSW

savefig("Fig_8b", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, NSW_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "NSW PDF forecasts (2025 to 2034)",
     ylim = c(0,0.09), yaxt="none")
axis(2, seq(0, 0.09, 0.015))
add_ylab("Density")
add_xlab("Age")
dev.off()

# VIC

savefig("Fig_8c", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, VIC_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "VIC PDF forecasts (2025 to 2034)",
     ylim = c(0,0.09), yaxt="none")
axis(2, seq(0, 0.09, 0.015))
add_ylab("Density")
add_xlab("Age")
dev.off()

# QLD

savefig("Fig_8d", width = 12, height = 12, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, QLD_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "QLD PDF forecasts (2025 to 2034)",
     ylim = c(0,0.09), yaxt="none")
axis(2, seq(0, 0.09, 0.015))
add_ylab("Density")
add_xlab("Age")
dev.off()

# SA

savefig("Fig_8f_SA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, SA_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "SA PDF (2025 to 2034)",
     ylim = range(range(AUS_PDF_fore), range(SA_PDF_fore)))
add_ylab("Density")
add_xlab("Age")
dev.off()

# WA

savefig("Fig_8f_WA", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, WA_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "WA PDF (2025 to 2034)",
     ylim = range(range(AUS_PDF_fore), range(WA_PDF_fore)))
add_ylab("Density")
add_xlab("Age")
dev.off()

# TAS

savefig("Fig_8f_TAS", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, TAS_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "TAS PDF (2025 to 2034)",
     ylim = range(range(AUS_PDF_fore), range(TAS_PDF_fore)))
add_ylab("Density")
add_xlab("Age")
dev.off()

# NT

savefig("Fig_8f_NT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, NT_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "NT PDF (2025 to 2034)",
     ylim = range(range(AUS_PDF_fore), range(NT_PDF_fore)))
add_ylab("Density")
add_xlab("Age")
dev.off()

# ACT

savefig("Fig_8f_ACT", width = 12, height = 10, toplines = 0.8, type = "png", pointsize = 12)
std_par()
plot(fts(ages, ACT_PDF_fore), xlab = "", ylab = "", lwd = CURVE_LWD,
     main = "ACT PDF (2025 to 2034)",
     ylim = range(range(AUS_PDF_fore), range(ACT_PDF_fore)))
add_ylab("Density")
add_xlab("Age")
dev.off()


