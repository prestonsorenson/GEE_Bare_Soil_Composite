#https://medium.com/applied-data-science/new-r-package-the-xgboost-explainer-51dd7d1aa211
library(ranger)
library(prospectr)
library(brnn)
library(earth)
library(Cubist)
library(xgboost)
library(Metrics)
library(ggplot2)
library(raster)
library(svMisc)


normalize <- function(x) {
  return ((x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x))))
}

setwd('/home/preston/OneDrive/Graduate Studies/Post_Doc/Shared/Data/Sask_Soil_Data')

#dat=read.csv("bare_soil_training_25Jan2021_PS.csv")
#dat_test_focal3=read.csv("bare_soil_training_focal3_26Jan2021_PS.csv")
dat=read.csv("bare_soil_training_focal10_26Jan2021_PS.csv")
dat=dat[!is.na(dat$b1),]
dat=dat[dat$b1!=0,]
dat=dat[dat$b1<5000,]
dat=dat[!is.na(dat$EXCH_CA),]

dat=dat[dat$EXCH_CA<80,]

#train-test split
train_index=read.csv('/home/preston/OneDrive/Papers/Historic_Carbon_Map/Analysis/train_index.csv')
test_index=read.csv('/home/preston/OneDrive/Papers/Historic_Carbon_Map/Analysis/test_index.csv')

dat_train=dat[dat$X %in% train_index$x,]
dat_test=dat[dat$X %in% test_index$x,]

#feature exploration
dat_sub=dat_train[,c(5, 22:28)]
cor(dat_sub)

plot(dat_sub$EXCH_CA~dat_sub$b1)
plot(dat_sub$EXCH_CA~dat_sub$b2)
plot(dat_sub$EXCH_CA~dat_sub$b3)
plot(dat_sub$EXCH_CA~dat_sub$b4)
plot(dat_sub$EXCH_CA~dat_sub$b5)
plot(dat_sub$EXCH_CA~dat_sub$b6)
plot(dat_sub$EXCH_CA~dat_sub$b7)

#model building
#ranger
model_val=ranger(EXCH_CA~b1+b2+b3+b4+b5+b6+b7, dat=dat_train, importance='impurity', quantreg = TRUE)
pred=predict(model_val, data=dat_test, type='quantiles', quantiles=c(0.25, 0.5, 0.75))
pred=predict(model_val, data=dat_test)
pred=pred$predictions
sort(importance(model_val))

boxplot(pred)
boxplot(dat_test$EXCH_CA)

plot(pred~dat_test$EXCH_CA, asp=1)
abline(0,1)

cor.test(pred, dat_test$EXCH_CA)
summary(lm(pred~dat_test$EXCH_CA))
rmse(dat_test$EXCH_CA, pred)
sd(dat_test$EXCH_CA)/rmse(dat_test$EXCH_CA, pred)

#create plot
plot_data=data.frame(dat_test$EXCH_CA, pred)
colnames(plot_data)=c("actual", 'predicted')
EXCH_CA_plot=ggplot(plot_data, (aes(x=actual, y=predicted))) + geom_point() + xlim(0, 80) + ylim(0, 80) + geom_abline() + xlab("Actual Exchangeable Calcium") + ylab("Predicted Exchangeable Calcium")
EXCH_CA_plot = EXCH_CA_plot + annotate("text", x=0, y = 80, label=expression(paste("R"^"2 ", "= 0.48"))) + annotate("text", x=0, y = 80, label=expression(paste("RMSE = 8.6 MEQ 100g"^""-1")))

EXCH_CA_plot

#predict results
model=ranger(EXCH_CA~b1+b2+b3+b4+b5+b6+b7, dat=dat, importance='impurity', quantreg = TRUE)
sort(importance(model))

raster_stack=stack('/media/preston/My Book/Saskatchewan/sk_bare_soil/sk_ls5_bare_soil_ndvi3_ndsi0_nbr1_focal10_filt/sk_ls5_bare_soil_ndsi_ndvi3_ndsi0_nbr1_focal10.tif')
raster_stack=raster_stack[[1:7]]

#raster_stack=aggregate(raster_stack, fact=3)

tiles=shapefile("/home/preston/OneDrive/Papers/Historic_Carbon_Map/Analysis/analysis_tiles.shp")

pred_results=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
#analysis loop
for (i in 1:length(tiles)){
  progress(i, max.value=length(tiles))
  tile_sub=tiles[tiles$id==i,]
  raster_sub=crop(raster_stack, tile_sub)
  raster_sub=rasterToPoints(raster_sub)
  xy=raster_sub[,1:2]
  raster_sub=raster_sub[,-c(1:2)]
  colnames(raster_sub)=c("b1", 'b2', 'b3','b4','b5', 'b6', 'b7')
  raster_sub=data.frame(raster_sub)
  pred_carb=predict(model, data=raster_sub, type='quantiles', quantiles=c(0.25, 0.5, 0.75))
  pred_carb=pred_carb$predictions
  pred_carb=data.frame(pred_carb)
  colnames(pred_carb)=c("pred_25", "pred_50", "pred_75")
  pred_carb$iqr=pred_carb$pred_75-pred_carb$pred_25
  pred_carb=data.frame(pred_carb, xy)
  pred_results=rbind(pred_results, pred_carb)
}

pred_results_25=pred_results[,c(1,5:6)]
pred_results_50=pred_results[,c(2,5:6)]
pred_results_75=pred_results[,c(3,5:6)]
pred_results_iqr=pred_results[,c(4,5:6)]
pred_results_iqr_ratio=pred_results[,c(4,5:6)]
pred_results_iqr_ratio$iqr=normalize(pred_results_iqr_ratio$iqr)

pred_results

coordinates(pred_results_25)=~x+y
coordinates(pred_results_50)=~x+y
coordinates(pred_results_75)=~x+y
coordinates(pred_results_iqr)=~x+y
coordinates(pred_results_iqr_ratio)=~x+y

pred_results_25=rasterFromXYZ(pred_results_25)
pred_results_50=rasterFromXYZ(pred_results_50)
pred_results_75=rasterFromXYZ(pred_results_75)
pred_results_iqr=rasterFromXYZ(pred_results_iqr)
pred_results_iqr_ratio=rasterFromXYZ(pred_results_iqr_ratio)

mask_1=raster_stack[[1]]
mask_1=crop(mask_1, pred_results_25)

pred_results_25=mask(pred_results_25, mask_1, maskvalue=0)
pred_results_50=mask(pred_results_50, mask_1, maskvalue=0)
pred_results_75=mask(pred_results_75, mask_1, maskvalue=0)
pred_results_iqr=mask(pred_results_iqr, mask_1, maskvalue=0)
pred_results_iqr_ratio=mask(pred_results_iqr_ratio, mask_1, maskvalue=0)

pred_stack=stack(pred_results_25, pred_results_50, pred_results_75, pred_results_iqr, pred_results_iqr_ratio)

names(pred_stack)=c("pred_results_25","pred_results_50","pred_results_75","pred_results_iqr","pred_results_iqr_ratio")

writeRaster(pred_stack, '/home/preston/OneDrive/Papers/Historic_Carbon_Map/Analysis/sk_predicted_clay.tif', format="GTiff", overwrite=TRUE)


