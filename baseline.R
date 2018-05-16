library(dplyr)
library(reshape2)
library(chron)
library(ggplot2)
library(data.table)
theme_set(theme_bw())
#molecular weigh of docosanol
molecular_weigh=228.42
# number of C-H functional groups of docosanol
num_FG=15*2+2-1
sample_number=14
MyData <- matrix(0, ncol = sample_number+1, nrow = 1865)
baselined_spectra<- matrix(0, ncol = sample_number+1, nrow = 1865)
background_spectra<- matrix(0, ncol = sample_number+1, nrow = 1865)
MyData <- data.frame(MyData)
baselined_spectra <- data.frame(baselined_spectra)
background_spectra <- data.frame(background_spectra)
names(MyData)[1]="Wavenumber"
names(baselined_spectra)[1]="Wavenumber"
names(background_spectra)[1]="Wavenumber"
for (i in 1:sample_number){
names(MyData)[i+1]=paste("Sample",i) 
names(baselined_spectra)[i+1]=paste("Sample",i) 
names(background_spectra)[i+1]=paste("Sample",i) 
sample<-read.csv(file=paste("180503_pentadecanol_s",i+1,".0.DPT",sep=""), header=TRUE, sep=",")
nam <- paste("sample",i,sep="")
assign(nam, sample)
names(sample)<-c("Wavenumber","Absorption")
MyData[,"Wavenumber"]=sample[,"Wavenumber"]
baselined_spectra[,"Wavenumber"]=sample[,"Wavenumber"]
background_spectra[,"Wavenumber"]=sample[,"Wavenumber"]
MyData[,paste("Sample",i)]=sample[,"Absorption"]
}
#MyData <- read.csv(file="docosanol2.csv", header=TRUE, sep=",")
#baselined_spectra= data.table(Sample1=rep(0,1866),Sample2=rep(0,1866),Sample3=rep(0,1866),Sample4=rep(0,1866),Sample5=rep(0,1866),Sample6=rep(0,1866))
MyData_weighing <- read.csv(file="weighing_pentadecanol2.csv", header=TRUE, sep=",")
peak_area=rep(0,sample_number)
for (num in 1:sample_number){
sample=paste("Sample",num)
S=MyData[,sample]
Wn=MyData[,"Wavenumber"]
# deleting useless parts of yhe spectrum
nan_count=0
sample= rep(0, length(Wn))
w= rep(0, length(Wn))
for (i in 1:length(Wn)){
  if((Wn[i]<1500) ||(Wn[i]<2950 && Wn[i]>2700) || (Wn[i]<2400 && Wn[i]>2200) || (Wn[i]<2500)){
    sample[i]=S[i]
  w[i]=0
  nan_count=nan_count
  }
  else {
    sample[i]=S[i]
    w[i]=1
  }
}
i=1
j=1
sample_corrected= rep(0, length(Wn)-nan_count)
Wn_corrected= rep(0, length(Wn)-nan_count)
while (i<length(Wn)){
  i=i+1
  if (is.na(sample[i])){
 
  }
  else {
    sample_corrected[j]=sample[i]
    Wn_corrected[j]=Wn[i]
    j=j+1
  }
}
#plot(Wn, S,type='l')
#lines(spline(Wn, S, n = 201), col = 2)
spectra.spl <- smooth.spline(Wn_corrected, sample_corrected,w=w, df=6)
pp <- predict(spectra.spl, Wn, deriv = 0)
#lines(pp, col = 2)

a=pp[[2]]
background_spectra[[paste('Sample',num)]]=a
baseline=S-a
baselined_spectra[[paste('Sample',num)]]=baseline


#integration
integ=0
#Wn=rev(Wn)
#baseline=rev(baseline)

#plot(Wn, baseline,xlim=c(2800,3200),ylim=c(0,0.05),type='l')
#lines(c(0,4000),c(0,0), col = 2)

for (i in 1:length(Wn)){

if ((Wn[i]>2700) && (Wn[i]<3000)){

  integ=integ+(baseline[i]+baseline[i+1])/2*(-Wn[i+1]+Wn[i])
}
}
peak_area[num]=integ
print(integ)
}
weigh_mean=rep(0,sample_number)
for (num in 1:4){
sample=paste("Sample.",num,sep = "")
weigh_mean[num]=(mean(MyData_weighing[[sample]])-mean(MyData_weighing[['Empty']]))*1000
}
for (num in 5:11){
  sample=paste("Sample.",num,sep = "")
  weigh_mean[num]=(mean(MyData_weighing[[sample]])-mean(MyData_weighing[['Empty.2']]))*1000
}
for (num in 12:14){
  sample=paste("Sample.",num,sep = "")
  weigh_mean[num]=(mean(MyData_weighing[[sample]])-mean(MyData_weighing[['Empty.3']]))*1000
}

FG_abundance=weigh_mean*num_FG/molecular_weigh

x <-FG_abundance
y <- peak_area


mydata= data.frame(x, y)

linear = function(k) {
  z <- list(xx = format(coef(k)[1], digits = 2),
            yy = format(abs(coef(k)[2]), digits = 2),
            r2 = format(summary(k)$r.squared, digits = 3));
  if (coef(k)[2] >= 0)  {
    eq <- substitute(italic(hat(y)) == xx + yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)
  } else {
    eq <- substitute(italic(hat(y)) == xx - yy %.% italic(x)*","~~italic(r)^2~"="~r2,z)   
  }
  as.character(as.expression(eq));               
}

fo = y ~ x
linplot <- ggplot(data = mydata, aes(x = x, y = y)) + geom_smooth(method = "lm", se=FALSE, color="blue", formula = fo) +  geom_point()+ggtitle("1-Pentadecanol") +
  theme(plot.title = element_text(hjust = 0.5))
linplot1 = linplot + annotate("text", x =10, y = 30, label = linear(lm(fo, mydata)), colour="black", size = 3, parse=TRUE)
linplot1
linplot1 +  labs(y =expression('Absorbance(cm'^-1*')'))+ 
  labs(x = expression(paste("C-H Abundance(", mu, "mol)"))) +
  lims(x=c(0,NA),y=c(0,NA) )


lf <- melt(MyData, id.vars=c("Wavenumber"))
ggplot(data = lf %>% filter(variable=="Sample 3"), aes(x=Wavenumber, y=value)) +
  geom_line(aes(colour=variable))+
   lims(x=c(4000, 0))



lf2 <- melt(background_spectra, id.vars=c("Wavenumber"))
ggplot(data = lf2, aes(x=Wavenumber, y=value)) +
  geom_line(aes(colour=variable))+
  lims(x=c(4000, 0))

lf3 <- melt(baselined_spectra, id.vars=c("Wavenumber"))
ggplot(data=lf3, aes(x=Wavenumber, y=value)) +
  geom_line(aes(colour=variable))+
  lims(x=c(3400, 2500),y=c(-0.1, 0.8))



#######for residual
a=lm(fo, mydata)
aa=resid(a)
plot(mydata$x, aa,  ylab="Residuals", xlab="Waiting Time",    main="Old Faithful Eruptions") 
abline(0,0)
anova(a)
