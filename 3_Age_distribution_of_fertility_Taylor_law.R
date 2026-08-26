#######################################################
# Taylor's law for original birth counts (not density)
# since mean and var can not be directly applied
#######################################################

# NSW

NSW_birth_mean = NSW_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    NSW_birth_var[iw] = var(AUS_fert_NSW[,iw])
    NSW_birth_mean[iw] = mean(AUS_fert_NSW[,iw])
    rm(iw)
}

NSW_birth_log_var = log(NSW_birth_var, base = 10)
NSW_birth_log_mean = log(NSW_birth_mean, base = 10)

NSW_birth_TL = round(lm(NSW_birth_log_var ~ NSW_birth_log_mean + 1)$coefficients, 4) # 1.1192
NSW_birth_R2 = round(summary(lm(NSW_birth_log_var ~ NSW_birth_log_mean + 1))$adj.r.squared, 4)

# VIC

VIC_birth_mean = VIC_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    VIC_birth_var[iw] = var(AUS_fert_VIC[,iw])
    VIC_birth_mean[iw] = mean(AUS_fert_VIC[,iw])
    rm(iw)
}

VIC_birth_log_var = log(VIC_birth_var, base = 10)
VIC_birth_log_mean = log(VIC_birth_mean, base = 10)

VIC_birth_TL = round(lm(VIC_birth_log_var ~ VIC_birth_log_mean + 1)$coefficients, 4) # 1.6525
VIC_birth_R2 = round(summary(lm(VIC_birth_log_var ~ VIC_birth_log_mean + 1))$adj.r.squared, 4)

# QLD

QLD_birth_mean = QLD_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    QLD_birth_var[iw] = var(AUS_fert_QLD[,iw])
    QLD_birth_mean[iw] = mean(AUS_fert_QLD[,iw])
    rm(iw)
}

QLD_birth_log_var = log(QLD_birth_var, base = 10)
QLD_birth_log_mean = log(QLD_birth_mean, base = 10)

QLD_birth_TL = round(lm(QLD_birth_log_var ~ QLD_birth_log_mean + 1)$coefficients, 4) # 1.3781
QLD_birth_R2 = round(summary(lm(QLD_birth_log_var ~ QLD_birth_log_mean + 1))$adj.r.squared, 4)


# SA

SA_birth_mean = SA_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    SA_birth_var[iw] = var(AUS_fert_SA[,iw])
    SA_birth_mean[iw] = mean(AUS_fert_SA[,iw])
    rm(iw)
}

SA_birth_log_var = log(SA_birth_var, base = 10)
SA_birth_log_mean = log(SA_birth_mean, base = 10)

SA_birth_TL = round(lm(SA_birth_log_var ~ SA_birth_log_mean + 1)$coefficients, 4) # 2.4441
SA_birth_R2 = round(summary(lm(SA_birth_log_var ~ SA_birth_log_mean + 1))$adj.r.squared, 4)

# WA

WA_birth_mean = WA_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    WA_birth_var[iw] = var(AUS_fert_WA[,iw])
    WA_birth_mean[iw] = mean(AUS_fert_WA[,iw])
    rm(iw)
}

WA_birth_log_var = log(WA_birth_var, base = 10)
WA_birth_log_mean = log(WA_birth_mean, base = 10)

WA_birth_TL = round(lm(WA_birth_log_var ~ WA_birth_log_mean + 1)$coefficients, 4) # 1.4510
WA_birth_R2 = round(summary(lm(WA_birth_log_var ~ WA_birth_log_mean + 1))$adj.r.squared, 4)

# TAS

TAS_birth_mean = TAS_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    TAS_birth_var[iw] = var(AUS_fert_TAS[,iw])
    TAS_birth_mean[iw] = mean(AUS_fert_TAS[,iw])
    rm(iw)
}

TAS_birth_log_var = log(TAS_birth_var, base = 10)
TAS_birth_log_mean = log(TAS_birth_mean, base = 10)

TAS_birth_TL = round(lm(TAS_birth_log_var ~ TAS_birth_log_mean + 1)$coefficients, 4) # 3.4595
TAS_birth_R2 = round(summary(lm(TAS_birth_log_var ~ TAS_birth_log_mean + 1))$adj.r.squared, 4)


# NT

NT_birth_mean = NT_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    NT_birth_var[iw] = var(AUS_fert_NT[,iw])
    NT_birth_mean[iw] = mean(AUS_fert_NT[,iw])
    rm(iw)
}

NT_birth_log_var = log(NT_birth_var, base = 10)
NT_birth_log_mean = log(NT_birth_mean, base = 10)

NT_birth_TL = round(lm(NT_birth_log_var ~ NT_birth_log_mean + 1)$coefficients, 4) # 1.2085
NT_birth_R2 = round(summary(lm(NT_birth_log_var ~ NT_birth_log_mean + 1))$adj.r.squared, 4)

# ACT

ACT_birth_mean = ACT_birth_var = vector("numeric", n_year)
for(iw in 1:n_year)
{
    ACT_birth_var[iw] = var(AUS_fert_ACT[,iw])
    ACT_birth_mean[iw] = mean(AUS_fert_ACT[,iw])
    rm(iw)
}

ACT_birth_log_var = log(ACT_birth_var, base = 10)
ACT_birth_log_mean = log(ACT_birth_mean, base = 10)

ACT_birth_TL = round(lm(ACT_birth_log_var ~ ACT_birth_log_mean + 1)$coefficients, 4) # 2.1489
ACT_birth_R2 = round(summary(lm(ACT_birth_log_var ~ ACT_birth_log_mean + 1))$adj.r.squared, 4)

# collect results

AUS_birth_TL = cbind(NSW_birth_TL, VIC_birth_TL, QLD_birth_TL, SA_birth_TL,
                     WA_birth_TL,  TAS_birth_TL, NT_birth_TL,  ACT_birth_TL)
colnames(AUS_birth_TL) = state
rownames(AUS_birth_TL) = c("Intercept", "Slope")

AUS_birth_TL_with_R2 = rbind(AUS_birth_TL, c(NSW_birth_R2, VIC_birth_R2, QLD_birth_R2,
                             SA_birth_R2, WA_birth_R2, TAS_birth_R2, NT_birth_R2, ACT_birth_R2))
rownames(AUS_birth_TL_with_R2)[3] = "Adjusted R2"

# export results

require(xtable)
xtable(AUS_birth_TL_with_R2, digits = 4)

