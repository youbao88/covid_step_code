###############################################################################
# Code to analyse time series of steps based on spectral analysis
library(mgcv)
library(lubridate)

dados = read.table("steps_covid_median_steps.csv", 
                   header = T, sep = ",", dec = ".")
dados$Date = as.Date(dados$Date, "%Y-%m-%d")
summary(dados$Date)
summary(dados$median)

plot(dados$Date, dados$median)

dados$Year = year(dados$Date)
table(dados$Year)

# data preparation (ny = number of years of data)
ny = 5                           # number of years of data
ndy = 365
nly = 1                          # number of leap years in the period
n = ny*ndy + nly                 # number of points in time series

# Choose number of sinusoidals and respective periods
m = 5     # Number of sinusoidals
m = 6     # Number of sinusoidals
m = 7     # Number of sinusoidals

# Periods of the sinusoidals (in days)
if (m == 5) {
  p = c(ndy/6, ndy/4, ndy/3, ndy/2, ndy) 
} else if (m == 6) {
  p = c(ndy/8, ndy/6, ndy/4, ndy/3, ndy/2, ndy)
} else if (m == 7) {
  p = c(ndy/12, ndy/8, ndy/6, ndy/4, ndy/3, ndy/2, ndy)
} else {
  p = c(ndy)                     # To be filled in if needed
}

# Reserves space for sinusoidals (sines and cosines) as well as  
# for index time series
t = seq(1, n, by = 1)
c = matrix(c(rep(0, n*m)), nrow = n, ncol = m, byrow = TRUE)
s = matrix(c(rep(0, n*m)), nrow = n, ncol = m, byrow = TRUE)

# generates sinusoidals (sines and cosines)
for (i in 1:n)
{
  for (j in 1:m)
  {
    c[i,j] = cos(((2*pi)/p[j])*i)
    s[i,j] = sin(((2*pi)/p[j])*i)
  }
}

# Sets up the design matrix
mat.design = as.data.frame(cbind(c, s))

# Name variables according to number of sinusoidals
if (m == 5) {
  names(mat.design) = c("cos1", "cos2", "cos3", "cos4", "cos5", "sin1", 
                        "sin2", "sin3", "sin4", "sin5")
} else if (m == 6) {
  names(mat.design) = c("cos1", "cos2", "cos3", "cos4", "cos5", "cos6", "sin1", 
                        "sin2", "sin3", "sin4", "sin5", "sin6")
} else if (m == 7) {
  names(mat.design) = c("cos1", "cos2", "cos3", "cos4", "cos5", "cos6", "cos7", 
                        "sin1", "sin2", "sin3", "sin4", "sin5", "sin6", "sin7")
} else {
  names(mat.design) = c("Complete")     # To be filled in if needed
}

# First date: 2017-01-01; Last date: 2021-12-31
as.Date(18992, origin = "1970-01-01")
sta = 18992 - n
time = seq(1:n) + sta
dat = as.Date(time, origin = "1970-01-01")
summary(dat)

mat.design = cbind(mat.design, dat, time)
mat.design$dat = as.Date(mat.design$dat, "%Y-%m-%d")
summary(mat.design$dat)

# merges steps data with sinusoidals (sines and cosines)
dados = merge(dados, mat.design, by.x = "Date", by.y = "dat")

# First model
if (m == 5) {
  spec = median ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + 
    sin4 + cos4 + sin5 + cos5
} else if (m == 6) {
  spec = median ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + 
    sin4 + cos4 + sin5 + cos5 + sin6 + cos6
}  else if (m == 7) {
  spec = median ~ sin1 + cos1 + sin2 + cos2 + sin3 + cos3 + 
    sin4 + cos4 + sin5 + cos5 + sin6 + cos6 + sin7 + cos7
} else {
  spec = median ~ 1       # To be filled in when needed
}

# Choose model according to number of sinusoidals
m1 = lm(spec, data = dados)
summary(m1)

dados$pred1 = fitted(m1)

plot(dados$Date, dados$median, col = "blue")
lines(dados$Date, dados$pred1, col = "red", lwd = 2)

dados$res1 = residuals(m1)

plot(dados$Date, dados$res1, col = "blue")
abline(0, 0, lty = 2, lwd = 2, col = "red")

a = mean(dados$median)

plot(dados$Date, dados$res1 + coef(m1)[1], col = "blue")
abline(coef(m1)[1], 0, lty = 2, lwd = 2, col = "red")

# Looks like df = 3 or 4 are good choices. df = 3 smooths the time series pretty well. 
ndf = ceiling(3*length(dados$time)/365)
ndf = ceiling(4*length(dados$time)/365)

m2 = gam(res1 ~ s(time, k = ndf), data = dados)
summary(m2)
dados$pred2 = fitted(m2)

plot(dados$Date, dados$res1 + coef(m1)[1], col = "blue")
lines(dados$Date, dados$pred2 + coef(m1)[1], lwd = 3, col = "red")
abline(coef(m1)[1], 0, lty = 2, lwd = 2, col = "red")

# preparing to save results
dd = subset(dados, select = c(1:2))

dd$y1 = dados$res1 + coef(m1)[1]
dd$f2 = dados$pred2 + coef(m1)[1]

# testing plot
plot(dd$Date, dd$y1, col = "blue")
lines(dd$Date, dd$f2, lwd = 3, col = "red")
abline(coef(m1)[1], 0, lty = 2, lwd = 2, col = "red")

# writing results
write.table(dd, "deseason_smoothed.csv", row.names = F, col.names = T, sep =";", dec = ",")
