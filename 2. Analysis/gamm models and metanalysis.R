
#Import packages
library(lubridate)
library(data.table)
library(dplyr)
library(broom)
library(mgcv)
library(metafor)
library(rio)
library(foreach)
library(doSNOW)
library(pbmcapply)


#Setting directories
setwd("G:/CES/CES Gemensamt/Folkhälsoappen/steps_and_covid/data")
#Data preparation
steps_data <-
    import(
      "G:/CES/CES Gemensamt/Folkhälsoappen/steps_and_covid/data/step_data_daily_aggregated_221114.csv ",
      header = T,
      sep = ","
    )%>% 
  mutate(Date = as.Date(as.character(Date), format = "%Y%m%d"))

#filter data to study period
steps_data <- setDT(steps_data)[Date >= "2019-01-01" & Date < "2022-01-01",]

#further cleaning and merging data sets
steps_data$one <- 1
count_obs <-aggregate(steps_data$one, by = list(steps_data$ParticipantUuid), sum)
colnames(count_obs) <- c("ParticipantUuid", "nobs")
count_obs$ID <- c(1:dim(count_obs)[1])

steps_data = merge(steps_data, count_obs, by = "ParticipantUuid")
steps_data = subset(steps_data,
                    select = c(ParticipantUuid, ID, Date, is_SE_holiday, nobs, Steps))
rm(count_obs)

# Add  waves and day of week (dow)
steps_data <- steps_data %>% 
  mutate(wave1 = as.factor(ifelse(Date >= "2020-03-11", 1, 0)),
         wave2 = as.factor(ifelse(Date >= "2020-09-01", 1, 0)),
         wave3 = as.factor(ifelse(Date >= "2021-02-01", 1, 0)),
         wave4 = as.factor(ifelse(Date >= "2021-07-01", 1, 0)),
         dow = relevel(as.factor(weekdays(Date)), 
                       ref = "måndag")
)

#create function that subsets data for each unique ID
subset_data <- function (X, i) {
  xx <- subset(X, select = -c(ParticipantUuid))[order(Date),]
  xx <- subset(xx, ID == i)
  
  date1 <- xx[1, 2]
  date1 <- mean(date1$Date)
  xx$time <- as.numeric(difftime(xx$Date,
                                 as.Date(unlist(date1)),
                                 units = c("days"))) + 1
  xx$is_SE_holiday <- as.factor(xx$is_SE_holiday)
  return(xx)
}

#Create a function to fit a GAMM model
fitGamm <- function(k) {
  Steps_data <- subset_data(steps_data, i = k)
  
  Holiday <- nlevels(Steps_data$is_SE_holiday)
  Dig_F <- 4 * ceiling(length(Steps_data$time) / 365)
  
  model <- if (Holiday == 2L) {
    gamm(
      Steps ~ dow + is_SE_holiday  + wave1 + wave2 + wave3 + wave4 +
        s(time, k = Dig_F),
      data = Steps_data ,
      method = "REML",
      correlation = corAR1(form = ~ time)
    )
  } else {
    gamm(
      Steps ~ dow  + wave1 + wave2 + wave3 + wave4 +
        s(time, k = Dig_F),
      data = Steps_data ,
      method = "REML",
      correlation = corAR1(form = ~ time)
    )
  }
  
  model_coef <-
    list(
      tidy (model$gam, parametric = T, conf.int = T) %>%
        filter(
          term == "wave11" | term == "wave21" | term == "wave31" |term == "wave41"
        )
    )
  
  return(model_coef)
}

## loop through individual models ##
#Setup parallel backend
n.cores <- detectCores() - 1
my.cluster <- makeCluster(n.cores,)
print(my.cluster)

registerDoSNOW(cl = my.cluster)
cat('If forEach is registered with parallel backend: ', foreach::getDoParRegistered())

# set data.table threads to 1
clusterEvalQ(cl = my.cluster, {
  library(data.table) # make sure to load the package on the cluster
  setDTthreads(1)     
})

#Setting progress bar
iterations <- length(unique(steps_data$ID))
pb <- progressBar(max = iterations, style = "ETA")
progress <- function(n)
  setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#Running the parallel model for fitting
gammModelResults <- foreach(
  i = unique(steps_data$ID),
  .packages = c('lubridate', 'mgcv', 'dplyr', 'broom'),
  .combine = "c",
  .options.snow = opts
) %dopar% {
  fitGamm(i)
  
}

#stop the cluster
stopCluster(my.cluster)
saveRDS(gammModelResults, file = "gammModelResults.RData")

#Close the progress bar
close(pb)
cat(
  'Sucessfully fitted all 560 individuals, saved the model as gammModelResults.RData under working directory'
)

#extract coefcients and se for further analysis
gamm_coeficients <- data.frame()

for (i in 1:560) {
  output = c(
    i,
    gammModelResults[[i]][1, ]$estimate,
    gammModelResults[[i]][1, ]$std.error,
    gammModelResults[[i]][2, ]$estimate,
    gammModelResults[[i]][2, ]$std.error,
    gammModelResults[[i]][3, ]$estimate,
    gammModelResults[[i]][3, ]$std.error,
    gammModelResults[[i]][4, ]$estimate,
    gammModelResults[[i]][4, ]$std.error
  )
  
  gamm_coeficients <- rbind(gamm_coeficients, output)
  
}

colnames(gamm_coeficients) <-
  c(
    "ID",
    "coef_wave1",
    "se_wave1",
    "coef_wave2",
    "se_wave2",
    "coef_wave3",
    "se_wave3",
    "coef_wave4",
    "se_wave4"
  )

##metanalysis and metaregression
#random effects models are fitted

#import and merge socio-demographic data

ses_data <-
  import(
    "G:/CES/CES Gemensamt/Folkhälsoappen/steps_and_covid/data/steps_scb_v2_221117.csv",
    header = T,
    sep = ","
  )
colnames(ses_data)[1] <- "ParticipantUuid"

#merge
ids <-steps_data[!duplicated(steps_data$ParticipantUuid),] #to get the ID
ses_data <-
  merge(ses_data, ids, by = "ParticipantUuid", all = T) %>%
  subset(select = c(
    ParticipantUuid,
    ID,
    alder,
    kon,
    yrke,
    dispinkke04,
    utb2kat,
    dispinkke04_q
  )) %>%
  merge(gamm_coeficients, by = "ID") %>%
  subset(!is.na(utb2kat)) %>%
  mutate(
    age_group = cut(
      alder,
      breaks = c(15, 30, 70, 100),
      right = F,
      labels = c("16-29", "30-69", "70-90")
    ),
    Office = relevel(
      recode_factor(
        factor(yrke),
        "På jobbet < 50%" = "< 50%",
        "På jobbet > 50%" = "> 50%",
        "Övriga" = "Other"
      ),
      ref = "< 50%"
    ),
    Income = relevel(factor(
      dispinkke04_q, labels = c("Q1", "Q2", "Q3", "Q4", "Q5")
    ), ref = "Q5"),
    Sex = factor(if_else(kon == 2, "Female", "Male")),
    Education = recode_factor(
      factor(utb2kat),
      "Hög utbildning" = "High",
      "Låg utbildning" = "Low"
    ),
  ) %>%
  subset(select = -c(kon, alder, utb2kat, yrke, dispinkke04, dispinkke04_q))


#loop thorugh all waves and vary covariates and extract results

#first reshape data to long
step_data_long <- ses_data %>% 
  pivot_longer(
    -c(ID, ParticipantUuid, age_group, Sex, Income, Education, Office),
    names_to = c(".value", "Wave"),
    names_sep = "_"
  )

step_data_long$Wave <-
  recode_factor(
    step_data_long$Wave,
    wave1 = 1,
    wave2 = 2,
    wave3 = 3,
    wave4 = 4
  )

##run and extract results

#loop through all and save all model contents  library(rlist)

metareg_models <- list() #storage

for (i in seq_along(levels(factor(step_data_long$Wave)))) {
  data = setDT(step_data_long)[Wave == i,] #extract data for each wave
  
  m1 = rma(yi = coef, sei = se, method = "ML", data = data, test = "knha")
  m2 = update(m1, mods = ~ age_group)
  m3 = update(m1, mods = ~ Sex)
  m4 = update(m1, mods = ~ Income)
  m5 = update(m1, mods = ~ Education)
  m6 = update(m1, mods = ~ Office)
  m7 = update(m1, mods = ~ age_group + Sex + Income + Education + Office) #joint model
  
  metareg_models [[paste0("Wave", i)]] = list(
    model1 = m1,
    model2 = m2,
    model3 = m3,
    model4 = m4,
    model5 = m5,
    model6 = m6,
    model7 = m7
  )
  
}

#saves results for further visualization 
saveRDS(metareg_models,
        file = "List of meta reg models per wave.RData") 


