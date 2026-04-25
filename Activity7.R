#Loading in libraries
library(ggplot2)
library(dplyr)
library(olsrr)
library(PerformanceAnalytics)
library(lubridate)
library(forecast)

## Multiple Regression Tutorial ----
#reading in greenhouse gas from reservoirs df
ghg <- read.csv("/cloud/project/activity07/Deemer_GHG_Data.csv")

# log transform methane fluxes
ghg$log.ch4 <- log(ghg$ch4+1)

#other log transformations
ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)

#identifying the different Region categories 
unique(ghg$Region)

# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)

# binary variable for alpine region
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)

# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)

# multiple regression
# creates a model object
mod.full <- lm(log.ch4 ~ airTemp +
                 log.age + mean.depth +
                 log.DIP +
                 log.precip + BorealV, data=ghg) #uses the data argument to specify dataframe

summary(mod.full)

#Checking assumptions
res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)

# shapiro-wilks test
shapiro.test(res.full)

#Residuals
plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:
reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)

# run stepwise
full.step <- ols_step_forward_aic(mod.full)
# view table
full.step 

# check full model
full.step$model

# plot AIC over time
plot(full.step )

#Predictions
# prediction with interval for predicting a point
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="prediction")

# look at prediction with 95% confidence interval of the mean
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="confidence")

#in class work ----
# You are asked to compile preliminary results from the Deemer data for climate policy 
# makers interested in understanding whether increased reservoir creation for 
# hydroelectric power would be expected to affect methane release. Develop a linear 
# model that improves on the model results from the tutorial. Ensure your model meets all 
# regression assumptions. 
# creates a model object
mod_new.full <- lm(log.ch4 ~ airTemp + log.age + mean.depth +
                 log.DIP +
                 log.precip + BorealV + HydroV + TropicalV + surface.area, data=ghg) #uses the data argument to specify dataframe

#Checking for multicolliniarity between precip and runoff
corr_precip_runoff <- lm(runoff ~ log.precip, data = ghg)
summary(corr_precip_runoff)

#checking for missing values
sum(is.na(ghg$volume)) #126
sum(is.na(ghg$chlorophyll.a)) #186
sum(is.na(ghg$HydroV)) #0
sum(is.na(ghg$surface.area)) #28

nrow(ghg$System)

summary(mod.full)
summary(mod_new.full)

#Checking assumptions
res_new.full <- rstandard(mod_new.full)
fit_new.full <- fitted.values(mod_new.full)

# qq plot
qqnorm(res_new.full, pch=19, col="grey50")
qqline(res_new.full)

# shapiro-wilks test
shapiro.test(res_new.full)

#Residuals
plot(fit_new.full,res_new.full, pch=19, col="grey50")
abline(h=0)

#to export csv file
write.csv(mod_new.full$coefficients, "/cloud/project/InClassRegression.csv", row.names = TRUE)


#Time Series ----
#tutorial
ETdat <- read.csv("/cloud/project/activity07/ETdata.csv")

#Displaying the different crop options
unique(ETdat$crop)

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

almondTrend <- almond_dec$trend
almondSeason <- almond_dec$seasonal

acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(almond_ts))

almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

model4 <- arima(almond_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4

# calculate fit
AR_fit1 <- almond_y - residuals(model1) 
AR_fit4 <- almond_y - residuals(model4)
#plot data
plot(almond_y)
# plot fit
points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newAlmond <- forecast(model4)
newAlmond

#make dataframe for plotting
newAlmondF <- data.frame(newAlmond)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newAlmondF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = almond, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(almond$date[1]),newAlmondF$dateF[24])+  # Plotting original data
  geom_line(data = newAlmondF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newAlmondF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")

# homework ----
#Question 1 FIXXXX: ----

#transforms co2 data using formula given
ghg$co2_transformed <- 1/(ghg$co2 + 1000)

#log transformation of chlorophyll.a
plot(ghg$chlorophyll.a)

ghg$log.chlp <- log(ghg$chlorophyll.a + 1)

plot(ghg$log.chlp)

#creates linear model to predict transformed co2 values
mod_co2.full <- lm(co2_transformed ~ airTemp + 
                     log.age +
                     mean.depth +
                     surface.area +
                     BorealV +
                     log.DIP +
                     log.precip +
                     HydroV +
                     TropicalV, data=ghg) 

summary(mod_co2.full)
nrow(ghg)
#checking for missing values
sum(!is.na(ghg$volume)) #126
sum(is.na(ghg$chlorophyll.a)) #186
sum(is.na(ghg$mean.depth)) #147


summary(mod_co2.full)
#number of observations in mod_co2.full
length(mod_co2.full$residuals)


#Checking assumptions
res_co2.full <- rstandard(mod_co2.full)
fit_co2.full <- fitted.values(mod_co2.full)

# qq plot
qqnorm(res_co2.full, pch=19, col="grey50")
qqline(res_co2.full)

# shapiro-wilks test
shapiro.test(res_co2.full)

#Residuals
plot(fit_co2.full,res_co2.full, pch=19, col="grey50")
abline(h=0)

#to export csv file
write.csv(mod_co2.full$coefficients, "/cloud/project/HomeWorkRegression.csv", row.names = TRUE)

#Question 3: ----

#adding almonds evapotranspiration data to the allcrops df
allcrops <- almond 
allcrops <- allcrops %>% rename(Almonds.ET.in = ET.in)

#list of the crop column names of interest
crop_names <- c("Pistachios", "Fallow/Idle Cropland", "Corn", "Grapes (Table/Raisin)")

#Plotting the decomposed plots of almonds created in tutorial
plot(almond_dec)

#looping through the crop_names list to plot all of the crops
for (i in 1:4) {
  current_crop <- ETdat %>% 
    filter(crop == crop_names[i]) %>% 
    group_by(date) %>%
    summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE))
  
  # Add ET data of crop to df to store data 
  allcrops <- merge(allcrops, current_crop, by = "date", all = TRUE)
  allcrops <- allcrops %>% rename(!!paste0(crop_names[i], ".ET.in") := ET.in)
  
  # ET time series
  current_crop_ts <- ts(current_crop$ET.in,
                  start = c(2016,1),
                  frequency= 12)
  
  # decompose almond ET time series
  current_crop_dec <- decompose(current_crop_ts)
  # plot decomposition
  plot(current_crop_dec)
}

#Question 4 ---- 
#Design an autoregressive model for pistachios and fallow/idle fields. 
#Forecast future evapotranspiration for each field so that water managers 
#can include estimates in their planning. Make a plot that includes 
#historical and forecasted evapotranspiration for the crops to present to 
#the water manager. Include a brief explanation of your autoregressive models.

unique(ETdat$crop)
# average fields for each month for pistachios
pistachio <- ETdat %>% # ET data
  filter(crop == "Pistachios") %>% # only use pistachio fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(pistachio, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# pistachio ET time series
pistachio_ts <- ts(pistachio$ET.in, 
                start = c(2016,1), 
                frequency= 12) 

# decompose almond ET time series
pistachio_dec <- decompose(pistachio_ts)
# plot decomposition
plot(pistachio_dec)

pistachioTrend <- pistachio_dec$trend
pistachioSeason <- pistachio_dec$seasonal

acf(na.omit(pistachio_ts), 
    lag.max = 24) 

pacf.plot <- pacf(na.omit(pistachio_ts))

pistachio_y <- na.omit(pistachio_ts)
model1 <- arima(pistachio_y ,  
                order = c(1,0,0)) 
model1

model2 <- arima(pistachio_y ,  
                order = c(2,0,0)) 

model3 <- arima(pistachio_y , # data 
                order = c(3,0,0)) # first number is AR order all other numbers get a 0 to keep AR format

#comparing different results
model1 #aic of 229.29
model2 #aic of 130.84
model3 #aic of 121.8
model4 #aic of 117.66


# calculate fit
AR_fit4_pistachio <- pistachio_y - residuals(model4) 

#plot data
plot(pistachio_y)

# plot fit
points(AR_fit4_pistachio, type = "l", col = "tomato3", lty = 2, lwd=2)

legend("topleft", c("data","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3"),
       bty="n")

#predict future ET using model 2
newPistachio <- forecast(model4)
newPistachio

#make dataframe for plotting
newPistachioF <- data.frame(newPistachio)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newPistachioF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = pistachio, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(pistachio$date[1]),newPistachioF$dateF[24])+  
  geom_line(data = newPistachioF, aes(x = dateF, y = Point.Forecast),
            col="red") + 
  geom_ribbon(data=newPistachioF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")

#repeat for fallow/idle fields

# average fields for each month for fallow/idle fields
fallowIdle <- ETdat %>% 
  filter(crop == "Fallow/Idle Cropland") %>% # only use Fallow/Idle fields
  group_by(date) %>% 
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(fallowIdle, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
fallowIdle_ts <- ts(fallowIdle$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose fallow/idle fields ET time series
fallowIdle_dec <- decompose(fallowIdle_ts)
# plot decomposition
plot(fallowIdle_dec)

fallowIdleTrend <- fallowIdle_dec$trend
fallowIdleSeason <- fallowIdle_dec$seasonal

acf(na.omit(fallowIdle_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(fallowIdle_ts))

fallowIdle_y <- na.omit(fallowIdle_ts)

#create different fallow/idle field models with varying orders up to 4
model1_fallow <- arima(fallowIdle_y , 
                       order = c(1,0,0))
model2_fallow <- arima(fallowIdle_y , 
                       order = c(2,0,0))
model3_fallow <- arima(fallowIdle_y , 
                order = c(3,0,0)) 
model4_fallow <- arima(fallowIdle_y , 
                       order = c(4,0,0))

#comparing different model results
model1_fallow #aic of 125.73
model2_fallow #aic of 93.65
model3_fallow #aic of 88.57
model4_fallow #aic of 87.74

# calculate fit
AR_fit3_fallow <- fallowIdle_y - residuals(model3_fallow)

#plot data
plot(fallowIdle_y)
# plot fit
points(AR_fit3_fallow, type = "l", col = "tomato3", lty = 2, lwd=2)
legend("topleft", c("data","AR3"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3"),
       bty="n")

newFallowIdle <- forecast(model3_fallow)
newFallowIdle

#make dataframe for plotting
newFallowIdleF <- data.frame(newFallowIdle)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newFallowIdleF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = fallowIdle, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(fallowIdle$date[1]),newFallowIdleF$dateF[24])+  # Plotting original data
  geom_line(data = newFallowIdleF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newFallowIdleF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")
