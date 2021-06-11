library(ranger)
library(prospectr)
library(Metrics)
library(ggplot2)
library(raster)
library(svMisc)
library(DescTools)

normalize <- function(x) {
  return ((x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x))))
}

setwd('/home/preston/OneDrive/Graduate Studies/Post_Doc/Shared/Data/Sask_Soil_Data')

dat=read.csv("bare_soil_training_focal10_26Jan2021_PS.csv")
dat=dat[!is.na(dat$b1),]
dat=dat[dat$b1!=0,]
dat=dat[dat$b1<5000,]
dat=dat[dat$CARB_ORG<10,]

boxplot(dat$CARB_ORG)

#train-test split
#original splits
#splits=kenStone(dat[,22:28], k=round(nrow(dat)*0.75), metric = 'mahal')

#dat_train=dat[splits$model,]
#dat_test=dat[splits$test,]
#load original splitting to keep same for all variables
train_index=read.csv('/home/preston/OneDrive/Papers/Historic_Carbon_Map/Analysis/train_index.csv')
test_index=read.csv('/home/preston/OneDrive/Papers/Historic_Carbon_Map/Analysis/test_index.csv')

`%notin%` <- Negate(`%in%`)

dat_train=dat[dat$X %notin% test_index$x,]
dat_test=dat[dat$X %in% test_index$x,]

#feature exploration
dat_sub=dat_train[,c(3, 22:28)]
cor(dat_sub)

dat_sub=dat_test[,c(3, 22:28)]
cor(dat_sub)

plot(dat_sub$CARB_ORG~dat_sub$b1)
plot(dat_sub$CARB_ORG~dat_sub$b2)
plot(dat_sub$CARB_ORG~dat_sub$b3)
plot(dat_sub$CARB_ORG~dat_sub$b4)
plot(dat_sub$CARB_ORG~dat_sub$b5)
plot(dat_sub$CARB_ORG~dat_sub$b6)
plot(dat_sub$CARB_ORG~dat_sub$b7)

#model building
#ranger
model=ranger(CARB_ORG~b1+b2+b3+b4+b5+b7, dat=dat_train, importance='impurity', quantreg = TRUE)
pred=predict(model, data=dat_test, type='quantiles', quantiles=c(0.25, 0.5, 0.75))
pred=predict(model, data=dat_test)
pred=pred$predictions
sort(importance(model))

boxplot(pred)
boxplot(dat_test$CARB_ORG)

plot(pred~dat_test$CARB_ORG, asp=1)
abline(0,1)

cor.test(pred, dat_test$CARB_ORG)
summary(lm(pred~dat_test$CARB_ORG))
rmse(dat_test$CARB_ORG, pred)
sd(dat_test$CARB_ORG)/rmse(dat_test$CARB_ORG, pred)
CCC(dat_test$CARB_ORG, pred)
bias(dat_test$CARB_ORG, pred)


#create plot
plot_data=data.frame(dat_test$CARB_ORG, pred)
colnames(plot_data)=c("actual", 'predicted')
carbon_plot=ggplot(plot_data, (aes(x=actual, y=predicted))) + geom_point() + xlim(0.5, 6) + ylim(0.5, 6) + geom_abline() + xlab("Actual Soil Organic Carbon (%)") + ylab("Predicted Soil Organic Carbon (%)")
carbon_plot = carbon_plot + annotate("text", x=1, y = 5.5, label=expression(paste("R"^"2 ", "= 0.55"))) + annotate("text", x=1, y=5.2, label="RMSE = 0.67%") + annotate("text", x=1, y=4.9, label=expression(paste(rho['c'], "= 0.71"))) + annotate("text", x=1, y=4.6, label='Bias = 0.04')
carbon_plot

#predict results
model=ranger(CARB_ORG~b1+b2+b3+b4+b5+b7, dat=dat, importance='impurity', quantreg = TRUE)
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

writeRaster(pred_stack, '/home/preston/OneDrive/Papers/Historic_Carbon_Map/Analysis/sk_predicted_carbon.tif', format="GTiff", overwrite=TRUE)


