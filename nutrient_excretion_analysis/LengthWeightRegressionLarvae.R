
library(ggplot2)
library(ggpubr)
library(car)
library(AICcmodavg)

caddis<-read.csv("Nutrient_excretion_ED\\LengthWeightRegressionLarvae.csv")
str(caddis)

ggplot(caddis, aes(x=sqrtWeight, y=Length))+
  geom_point()

cor(caddis_new$sqrtWeight, caddis_new$Length, method = c("pearson")) #0.975
cor.test(caddis_new$sqrtWeight, caddis_new$Length, method = c("pearson")) #p < 0.001

model1<-lm(data=caddis, Length~Weight*Site)
model1.1<-lm(data=caddis, Length~Weight + Site)
model1.2<-lm(data=caddis, Length~Weight)

models<-list(model1, model1.1, model1.2)
names(models) <- c("1.0", "1.1", "1.2")
aictab(cand.set = models) #model 1.2 is the best. Just include length~weight

summary(model1.2)
plot(model1.2)

sqrtWeight<-sqrt(caddis$Weight)
cbind(sqrtWeight, caddis)
summary(model2<-lm(data=caddis, Length~sqrt(caddis$Weight)))
plot(model2)

squaredWeight<-caddis$Weight^2
cbind(squaredWeight, caddis)
summary(model3<-lm(data=caddis, Length~Weight + squaredWeight))
plot(model3)

logWeight<-caddis$Weight
cbind(logWeight, caddis)
summary(model4<-lm(data=caddis, Length~logWeight))
plot(model4)

summary(model5<-lm(Length~poly(Weight,2), data=caddis))
plot(model5)


caddis_new<-cbind(caddis, sqrtWeight, squaredWeight, logWeight)

