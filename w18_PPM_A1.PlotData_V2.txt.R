#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[1] Figure 2
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.1] Hatzis_GSE25066 data  (classification)   --> data
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Hatzis_GSE25066/Clinical_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Breast_neoadjuvant_Symbol_oncotypeDX_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig2_data_oncotypeDX_discovery.rda"


myList = list(NULL)


data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info <- read.table(myinf1, header=T, sep="\t",  quote="")
source = c(rep("discovery", 310), rep("validation", 198))
info$source = source
info = info[!is.na(info$pathologic_response_pcr_rd), ]

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)
rawxx.data = data
se = c("mytf","age_years", "er_status_ihc","clinical_t_stage", "grade", "pam50_class", "drfs_even_time_years", "drfs_1_event_0_censored","pathologic_response_pcr_rd",  "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[!is.na(data$er_status_ihc), ]

#------------------------------
pos = row.names(data)[data$pathologic_response_pcr_rd=="pCR"]
neg = row.names(data)[data$pathologic_response_pcr_rd=="RD"]
pos = pos[!is.na(pos)]
neg = neg[!is.na(neg)]

xx = data
tmp = row.names(xx)
xx = xx$mytf
names(xx) = tmp
xx = sort(xx)
xx= names(xx)%in%pos
fp = 1-xx 
tp = xx
for(j in length(xx):2)
{
	fp[j-1]= fp[j]+fp[j-1]
	tp[j-1]= tp[j]+tp[j-1]
}
fp = fp/length(neg)
tp = tp/length(pos)	
xx = c(1, fp, 0)
yy = c(1, tp, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:length(tmp1))
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
## 0.7311121
myList[[1]] = myauc
names(myList)[1] = "AUC.OncotypeDX.Uni"
myList[[2]] = cbind(xx, yy)
names(myList)[2] = "ROC.OncotypeDX.Uni"
save(myList, file=myoutf1)

#-----------------------
library(randomForest)
xx = table(data$pathologic_response_pcr_rd)
sam.siz = min(xx)
myfit <- randomForest(pathologic_response_pcr_rd~mytf+age_years+er_status_ihc+clinical_t_stage, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
myfit

pdat = data[data$pathologic_response_pcr_rd=="pCR", ]
ndat = data[data$pathologic_response_pcr_rd=="RD", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = rep(0, 10)
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(pathologic_response_pcr_rd~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"pathologic_response_pcr_rd"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="pCR")
		bb = sum(myres[,"pCR"]<thr[i] & myres[,1]=="pCR" )
		cc = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="RD")
		dd = sum(myres[,"pCR"]<thr[i] & myres[,1]=="RD")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc[p] = sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}
myauc1 = mean(myauc)
myauc1		## 0.7459288
myx1 = myx/10
myy1 = myy/10

load(file= myoutf1)
myList[[3]] = myauc1
names(myList)[3] = "AUC.CPM.RF"
myList[[4]] = cbind(myx1, myy1)
names(myList)[4] = "ROC.CPM.RF"
save(myList, file=myoutf1)

#----------------------------------------------
tmp = predict(myfit, newdata = data, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)

mod1 = randomForest(class~ age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)
## 0.7868852

## 10-fold CV


pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = rep(0, 10)
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(class~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,1], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc[p] =  sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc2 = mean(myauc)
myx2 = myx/10
myy2 = myy/10
##  AUC = 0.8457247

load(file= myoutf1)
myList[[5]] = myauc2
names(myList)[5] = "AUC.PPM.RF"
myList[[6]] = cbind(myx2, myy2)
names(myList)[6] = "ROC.PPM.RF"
save(myList, file=myoutf1)


load(file= myoutf1)
myList[[7]] = mod1[["importance"]]
names(myList)[7] = "PPM.RF.RelativeImportance"
save(myList, file=myoutf1)


#----------------------------------------------

wilcox.test(data$age_years[data$class=="Y"], data$age_years[data$class=="N"] )
## p-value = 0.002288
> mean(data$age_years[data$class=="Y"])
[1] 51.03109
> mean(data$age_years[data$class=="N"])
[1] 46.94627

aa = sum(data$er_status_ihc=="P" & data$class=="Y", na.rm=T)
bb = sum(data$er_status_ihc=="P" & data$class=="N", na.rm=T)
cc = sum(data$er_status_ihc=="N" & data$class=="Y", na.rm=T)
dd = sum(data$er_status_ihc=="N" & data$class=="N", na.rm=T)
xx = matrix(c(aa, bb, cc, dd), 2,2)
chisq.test(xx)$p.value		##1.266434e-12

load(file= myoutf1)
tmp = rawxx.data[names(class),]
tmp = cbind(class, tmp)

myList[[8]] = tmp
names(myList)[8] = "Haztis.Class.Predictors"
save(myList, file=myoutf1)

#----------------------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[!is.na(data$er_status_ihc),]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(myfit, newdata= data, type="prob")
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(row.names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx,]

#------------------------------
pred = ifelse(risk[, "pCR"]>0.5, "pCR", "RD")
response = as.character(data$pathologic_response_pcr_rd)
res = data.frame(response, pred, predictability)
myorder = order(res[, "Y"], decreasing=T)
res = res[myorder,]

load(file= myoutf1)
myList[[9]] = res
names(myList)[9] = "dis2val.accVScscore"
save(myList, file=myoutf1)


xx = ifelse(res[,1]==res[,2],1,0)
for(k in 2:length(xx))
{
	xx[k] = xx[k-1]+xx[k]
}
acc = xx/(1:length(xx)) 
plot(acc, type="b", pch=20)

w.siz = 20
step = 1
cnum = ceiling((length(acc)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(acc))
	yy[k] = mean(acc[sta:end])
}
plot(yy, pch=20, type="b")
yy



#---------------
pred = risk[, "pCR"]
response = as.character(data$pathologic_response_pcr_rd)
res = data.frame(response, pred, predictability)
myorder = order(res[, "Y"], decreasing=T)
res = res[myorder,]
pos = row.names(res)[res$response=="pCR"]
neg = row.names(res)[res$response=="RD"]
pos = pos[!is.na(pos)]
neg = neg[!is.na(neg)]


w.siz = 50
step = 10
cnum = ceiling((nrow(res)-w.siz)/step)
auc = rep(0, cnum)
for(k in 1:cnum)
{
	sta = 1
	end = min(w.siz+(k-1)*step, nrow(res))	
	xx = res[sta:end,]
	tmp = row.names(xx)
	xx = xx$pred
	names(xx) = tmp
	xx = sort(xx)
	xx= names(xx)%in%pos
	fp = 1-xx 
	tp = xx
	for(j in length(xx):2)
	{
		fp[j-1]= fp[j]+fp[j-1]
		tp[j-1]= tp[j]+tp[j-1]
	}
	fp = fp/length(neg)
	tp = tp/length(pos)	
	xx = c(1, fp, 0)
	yy = c(1, tp, 0)
	tmp1 = tmp2 = rep(0,length(xx)-1)
	for(i in 1:length(tmp1))
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	auc[k] = sum(tmp1*tmp2)
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[1.2] Hatzis_GSE25066 data  (classification)   --> validation  data
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Hatzis_GSE25066/Clinical_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Breast_neoadjuvant_Symbol_oncotypeDX_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig2_data_oncotypeDX_validation.rda"


myList = list(NULL)
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info <- read.table(myinf1, header=T, sep="\t",  quote="")
source = c(rep("discovery", 310), rep("validation", 198))
info$source = source
info = info[!is.na(info$pathologic_response_pcr_rd), ]

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)
se = c("mytf","age_years", "er_status_ihc","clinical_t_stage", "grade", "pam50_class", "drfs_even_time_years", "drfs_1_event_0_censored","pathologic_response_pcr_rd",  "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[!is.na(data$er_status_ihc), ]

#------------------------------
pos = row.names(data)[data$pathologic_response_pcr_rd=="pCR"]
neg = row.names(data)[data$pathologic_response_pcr_rd=="RD"]
pos = pos[!is.na(pos)]
neg = neg[!is.na(neg)]

xx = data
tmp = row.names(xx)
xx = xx$mytf
names(xx) = tmp
xx = sort(xx)
xx= names(xx)%in%pos
fp = 1-xx 
tp = xx
for(j in length(xx):2)
{
	fp[j-1]= fp[j]+fp[j-1]
	tp[j-1]= tp[j]+tp[j-1]
}
fp = fp/length(neg)
tp = tp/length(pos)	
xx = c(1, fp, 0)
yy = c(1, tp, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:length(tmp1))
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
## 0.6284687
myList[[1]] = myauc
names(myList)[1] = "AUC.OncotypeDX.Uni"
myList[[2]] = cbind(xx, yy)
names(myList)[2] = "ROC.OncotypeDX.Uni"
save(myList, file=myoutf1)


library(randomForest)
xx = table(data$pathologic_response_pcr_rd)
sam.siz = min(xx)
myfit <- randomForest(pathologic_response_pcr_rd~mytf+age_years+er_status_ihc+clinical_t_stage, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
myfit

pdat = data[data$pathologic_response_pcr_rd=="pCR", ]
ndat = data[data$pathologic_response_pcr_rd=="RD", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = rep(0, 10)
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(pathologic_response_pcr_rd~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"pathologic_response_pcr_rd"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="pCR")
		bb = sum(myres[,"pCR"]<thr[i] & myres[,1]=="pCR" )
		cc = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="RD")
		dd = sum(myres[,"pCR"]<thr[i] & myres[,1]=="RD")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc[p] = sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}
myauc1 = mean(myauc)
myauc1		## 0.7099321
myx1 = myx/10
myy1 = myy/10

load(file= myoutf1)
myList[[3]] = myauc1
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx1, myy1)
names(myList)[4] = "ROC.PPM.RF"
save(myList, file=myoutf1)

#----------------------------------------------
tmp = predict(myfit, newdata = data, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)

mod1 = randomForest(class~ age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)


## 10-fold CV


pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = rep(0, 10)
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(class~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,1], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc[p] =  sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc2 = mean(myauc)
myx2 = myx/10
myy2 = myy/10
##  AUC = 0.6564373

load(file= myoutf1)
myList[[5]] = myauc2
names(myList)[5] = "AUC.PPM.RF"
myList[[6]] = cbind(myx2, myy2)
names(myList)[6] = "ROC.PPM.RF"
save(myList, file=myoutf1)


load(file= myoutf1)
myList[[7]] = myfit[["importance"]]
names(myList)[7] = "PPM.RF.RelativeImportance"
save(myList, file=myoutf1)


#----------------------------------------------

wilcox.test(data$age_years[data$class=="Y"], data$age_years[data$class=="N"] )
## p-value = 0.002288
> mean(data$age_years[data$class=="Y"])
[1] 51.03109
> mean(data$age_years[data$class=="N"])
[1] 46.94627

aa = sum(data$er_status_ihc=="P" & data$class=="Y", na.rm=T)
bb = sum(data$er_status_ihc=="P" & data$class=="N", na.rm=T)
cc = sum(data$er_status_ihc=="N" & data$class=="Y", na.rm=T)
dd = sum(data$er_status_ihc=="N" & data$class=="N", na.rm=T)
xx = matrix(c(aa, bb, cc, dd), 2,2)
chisq.test(xx)$p.value		##0.002990014

load(file= myoutf1)
myList[[8]] = data
names(myList)[8] = "Haztis.Class.Predictors"
save(myList, file=myoutf1)

#----------------------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[!is.na(data$er_status_ihc),]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(myfit, newdata= data, type="prob")
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(row.names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx,]

#------------------------------
pred = ifelse(risk[, "pCR"]>0.5, "pCR", "RD")
response = as.character(data$pathologic_response_pcr_rd)
res = data.frame(response, pred, predictability)
myorder = order(res[, "Y"], decreasing=T)
res = res[myorder,]

load(file= myoutf1)
myList[[9]] = res
names(myList)[9] = "val2dis.accVScscore"
save(myList, file=myoutf1)


xx = ifelse(res[,1]==res[,2],1,0)
for(k in 2:length(xx))
{
	xx[k] = xx[k-1]+xx[k]
}
acc = xx/(1:length(xx)) 
plot(acc, type="b", pch=20)

w.siz = 20
step = 1
cnum = ceiling((length(acc)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(acc))
	yy[k] = mean(acc[sta:end])
}
plot(yy, pch=20, type="b")
yy



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[2] Figure 3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.1] Hatzis_GSE25066 data  (classification)   --> data
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Hatzis_GSE25066/Clinical_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Breast_neoadjuvant_Symbol_E2F4_iRAS.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig3_data_E2F4_discovery.rda"


                 MeanDecreaseGini
mytf                    22.721297
age_years               15.273004
er_status_ihc           10.547034
clinical_t_stage         5.716373
> 

                 MeanDecreaseGini
age_years               13.411603
er_status_ihc            9.583210
clinical_t_stage         5.051283
pam50_class             11.251320


myList = list(NULL)
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info <- read.table(myinf1, header=T, sep="\t",  quote="")
source = c(rep("discovery", 310), rep("validation", 198))
info$source = source
info = info[!is.na(info$pathologic_response_pcr_rd), ]

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "all.intersection.ES"]
data = cbind(mytf, info)
se = c("mytf","age_years", "er_status_ihc","clinical_t_stage", "grade", "pam50_class", "drfs_even_time_years", "drfs_1_event_0_censored","pathologic_response_pcr_rd",  "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[!is.na(data$er_status_ihc), ]

#------------------------------
pos = row.names(data)[data$pathologic_response_pcr_rd=="pCR"]
neg = row.names(data)[data$pathologic_response_pcr_rd=="RD"]
pos = pos[!is.na(pos)]
neg = neg[!is.na(neg)]

xx = data
tmp = row.names(xx)
xx = xx$mytf
names(xx) = tmp
xx = sort(xx)
xx= names(xx)%in%pos
fp = 1-xx 
tp = xx
for(j in length(xx):2)
{
	fp[j-1]= fp[j]+fp[j-1]
	tp[j-1]= tp[j]+tp[j-1]
}
fp = fp/length(neg)
tp = tp/length(pos)	
xx = c(1, fp, 0)
yy = c(1, tp, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:length(tmp1))
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
## 0.7311121
myList[[1]] = myauc
names(myList)[1] = "AUC.E2F4.Uni"
myList[[2]] = cbind(xx, yy)
names(myList)[2] = "ROC.E2F4.Uni"
save(myList, file=myoutf1)


library(randomForest)
xx = table(data$pathologic_response_pcr_rd)
sam.siz = min(xx)
myfit <- randomForest(pathologic_response_pcr_rd~mytf+age_years+er_status_ihc+clinical_t_stage, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
myfit
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)

pdat = data[data$pathologic_response_pcr_rd=="pCR", ]
ndat = data[data$pathologic_response_pcr_rd=="RD", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(pathologic_response_pcr_rd~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"pathologic_response_pcr_rd"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="pCR")
		bb = sum(myres[,"pCR"]<thr[i] & myres[,1]=="pCR" )
		cc = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="RD")
		dd = sum(myres[,"pCR"]<thr[i] & myres[,1]=="RD")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc1 = myauc/10
myauc1		## 
myx1 = myx/10
myy1 = myy/10

load(file= myoutf1)
myList[[3]] = myauc1
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx1, myy1)
names(myList)[4] = "ROC.PPM.RF"
save(myList, file=myoutf1)

#----------------------------------------------
tmp = predict(myfit, newdata = data, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)

mod1 = randomForest(class~ age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(class~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T,data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,1], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc2 = myauc/10
myx2 = myx/10
myy2 = myy/10

load(file= myoutf1)
myList[[5]] = myauc2
names(myList)[5] = "AUC.PPM.RF"
myList[[6]] = cbind(myx2, myy2)
names(myList)[6] = "ROC.PPM.RF"
save(myList, file=myoutf1)


load(file= myoutf1)
myList[[7]] = myfit[["importance"]]
names(myList)[7] = "PPM.RF.RelativeImportance"
save(myList, file=myoutf1)


#----------------------------------------------

wilcox.test(data$age_years[data$class=="Y"], data$age_years[data$class=="N"] )
mean(data$age_years[data$class=="Y"])
mean(data$age_years[data$class=="N"])

aa = sum(data$er_status_ihc=="P" & data$class=="Y", na.rm=T)
bb = sum(data$er_status_ihc=="P" & data$class=="N", na.rm=T)
cc = sum(data$er_status_ihc=="N" & data$class=="Y", na.rm=T)
dd = sum(data$er_status_ihc=="N" & data$class=="N", na.rm=T)
xx = matrix(c(aa, bb, cc, dd), 2,2)
chisq.test(xx)$p.value		##

load(file= myoutf1)
myList[[8]] = data
names(myList)[8] = "Haztis.Class.Predictors"
save(myList, file=myoutf1)

#----------------------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[!is.na(data$er_status_ihc),]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(myfit, newdata= data, type="prob")
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(row.names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx,]

#------------------------------
pred = ifelse(risk[, "pCR"]>0.5, "pCR", "RD")
response = as.character(data$pathologic_response_pcr_rd)
res = data.frame(response, pred, predictability)
myorder = order(res[, "Y"], decreasing=T)
res = res[myorder,]

load(file= myoutf1)
myList[[9]] = res
names(myList)[9] = "dis2val.accVScscore"
save(myList, file=myoutf1)


xx = ifelse(res[,1]==res[,2],1,0)
for(k in 2:length(xx))
{
	xx[k] = xx[k-1]+xx[k]
}
acc = xx/(1:length(xx)) 
plot(acc, type="b", pch=20)

w.siz = 20
step = 1
cnum = ceiling((length(acc)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(acc))
	yy[k] = mean(acc[sta:end])
}
plot(yy, pch=20, type="b")
yy



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.2] Hatzis_GSE25066 data  (classification)   --> validation  data
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Hatzis_GSE25066/Clinical_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Breast_neoadjuvant_Symbol_E2F4_iRAS.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig3_data_E2F4_validation.rda"

                 MeanDecreaseGini
mytf                    20.122117
age_years               13.954295
er_status_ihc            2.854112
clinical_t_stage         3.523094

                 MeanDecreaseGini
age_years                6.921703
er_status_ihc            1.386954
clinical_t_stage         1.790621
pam50_class              3.329662


myList = list(NULL)
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info <- read.table(myinf1, header=T, sep="\t",  quote="")
source = c(rep("discovery", 310), rep("validation", 198))
info$source = source
info = info[!is.na(info$pathologic_response_pcr_rd), ]

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "all.intersection.ES"]
data = cbind(mytf, info)
se = c("mytf","age_years", "er_status_ihc","clinical_t_stage", "grade", "pam50_class", "drfs_even_time_years", "drfs_1_event_0_censored","pathologic_response_pcr_rd",  "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[!is.na(data$er_status_ihc), ]

#------------------------------
pos = row.names(data)[data$pathologic_response_pcr_rd=="pCR"]
neg = row.names(data)[data$pathologic_response_pcr_rd=="RD"]
pos = pos[!is.na(pos)]
neg = neg[!is.na(neg)]

xx = data
tmp = row.names(xx)
xx = xx$mytf
names(xx) = tmp
xx = sort(xx)
xx= names(xx)%in%pos
fp = 1-xx 
tp = xx
for(j in length(xx):2)
{
	fp[j-1]= fp[j]+fp[j-1]
	tp[j-1]= tp[j]+tp[j-1]
}
fp = fp/length(neg)
tp = tp/length(pos)	
xx = c(1, fp, 0)
yy = c(1, tp, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:length(tmp1))
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)

myList[[1]] = myauc
names(myList)[1] = "AUC.E2F4.Uni"
myList[[2]] = cbind(xx, yy)
names(myList)[2] = "ROC.E2F4.Uni"
save(myList, file=myoutf1)


library(randomForest)
xx = table(data$pathologic_response_pcr_rd)
sam.siz = min(xx)
myfit <- randomForest(pathologic_response_pcr_rd~mytf+age_years+er_status_ihc+clinical_t_stage, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
myfit
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)

pdat = data[data$pathologic_response_pcr_rd=="pCR", ]
ndat = data[data$pathologic_response_pcr_rd=="RD", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(pathologic_response_pcr_rd~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"pathologic_response_pcr_rd"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="pCR")
		bb = sum(myres[,"pCR"]<thr[i] & myres[,1]=="pCR" )
		cc = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="RD")
		dd = sum(myres[,"pCR"]<thr[i] & myres[,1]=="RD")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc1 = myauc/10
myauc1		## 
myx1 = myx/10
myy1 = myy/10

load(file= myoutf1)
myList[[3]] = myauc1
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx1, myy1)
names(myList)[4] = "ROC.PPM.RF"
save(myList, file=myoutf1)

#----------------------------------------------
tmp = predict(myfit, newdata = data, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)

mod1 = randomForest(class~ age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:4)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))		
		fit <- randomForest(class~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,1], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc2 = myauc/4
myx2 = myx/4
myy2 = myy/4


load(file= myoutf1)
myList[[5]] = myauc2
names(myList)[5] = "AUC.PPM.RF"
myList[[6]] = cbind(myx2, myy2)
names(myList)[6] = "ROC.PPM.RF"
save(myList, file=myoutf1)


load(file= myoutf1)
myList[[7]] = myfit[["importance"]]
names(myList)[7] = "PPM.RF.RelativeImportance"
save(myList, file=myoutf1)

                MeanDecreaseGini
age_years               26.226418
er_status_ihc            3.817344
clinical_t_stage         7.239049
pam50_class             10.108766

                 MeanDecreaseGini
mytf                    31.937517
age_years               20.430926
er_status_ihc            3.512371
clinical_t_stage         4.747124


#----------------------------------------------

wilcox.test(data$age_years[data$class=="Y"], data$age_years[data$class=="N"] )

mean(data$age_years[data$class=="Y"])
mean(data$age_years[data$class=="N"])

aa = sum(data$er_status_ihc=="P" & data$class=="Y", na.rm=T)
bb = sum(data$er_status_ihc=="P" & data$class=="N", na.rm=T)
cc = sum(data$er_status_ihc=="N" & data$class=="Y", na.rm=T)
dd = sum(data$er_status_ihc=="N" & data$class=="N", na.rm=T)
xx = matrix(c(aa, bb, cc, dd), 2,2)
chisq.test(xx)$p.value		##

load(file= myoutf1)
myList[[8]] = data
names(myList)[8] = "Haztis.Class.Predictors"
save(myList, file=myoutf1)

#----------------------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[!is.na(data$er_status_ihc),]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(myfit, newdata= data, type="prob")
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(row.names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx,]

#------------------------------
pred = ifelse(risk[, "pCR"]>0.5, "pCR", "RD")
response = as.character(data$pathologic_response_pcr_rd)
res = data.frame(response, pred, predictability)
myorder = order(res[, "Y"], decreasing=T)
res = res[myorder,]

load(file= myoutf1)
myList[[9]] = res
names(myList)[9] = "val2dis.accVScscore"
save(myList, file=myoutf1)


xx = ifelse(res[,1]==res[,2],1,0)
for(k in 2:length(xx))
{
	xx[k] = xx[k-1]+xx[k]
}
acc = xx/(1:length(xx)) 
plot(acc, type="b", pch=20)

w.siz = 20
step = 1
cnum = ceiling((length(acc)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(acc))
	yy[k] = mean(acc[sta:end])
}
plot(yy, pch=20, type="b")
yy



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.3] Hatzis_GSE25066 data  (classification)   --> gene70 (MammaPrint)
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Hatzis_GSE25066/Clinical_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Symbol_gene70_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig3_data_gene70_discovery.rda"


                 MeanDecreaseGini
mytf                    27.822950
age_years               14.827424
er_status_ihc            6.649719
clinical_t_stage         5.390300


                 MeanDecreaseGini
age_years               14.621984
er_status_ihc            8.995240
clinical_t_stage         3.790989
pam50_class             15.259727


myList = list(NULL)
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info <- read.table(myinf1, header=T, sep="\t",  quote="")
source = c(rep("discovery", 310), rep("validation", 198))
info$source = source
info = info[!is.na(info$pathologic_response_pcr_rd), ]

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)
se = c("mytf","age_years", "er_status_ihc","clinical_t_stage", "grade", "pam50_class", "drfs_even_time_years", "drfs_1_event_0_censored","pathologic_response_pcr_rd",  "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[!is.na(data$er_status_ihc), ]

#------------------------------
pos = row.names(data)[data$pathologic_response_pcr_rd=="pCR"]
neg = row.names(data)[data$pathologic_response_pcr_rd=="RD"]
pos = pos[!is.na(pos)]
neg = neg[!is.na(neg)]

xx = data
tmp = row.names(xx)
xx = xx$mytf
names(xx) = tmp
xx = sort(xx)
xx= names(xx)%in%pos
fp = 1-xx 
tp = xx
for(j in length(xx):2)
{
	fp[j-1]= fp[j]+fp[j-1]
	tp[j-1]= tp[j]+tp[j-1]
}
fp = fp/length(neg)
tp = tp/length(pos)	
xx = c(1, fp, 0)
yy = c(1, tp, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:length(tmp1))
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)
## 0.7925156
myList[[1]] = myauc
names(myList)[1] = "AUC.gene70.Uni"
myList[[2]] = cbind(xx, yy)
names(myList)[2] = "ROC.gene70.Uni"
save(myList, file=myoutf1)


library(randomForest)
xx = table(data$pathologic_response_pcr_rd)
sam.siz = min(xx)
myfit <- randomForest(pathologic_response_pcr_rd~mytf+age_years+er_status_ihc+clinical_t_stage, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
myfit
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)

pdat = data[data$pathologic_response_pcr_rd=="pCR", ]
ndat = data[data$pathologic_response_pcr_rd=="RD", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(pathologic_response_pcr_rd~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"pathologic_response_pcr_rd"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="pCR")
		bb = sum(myres[,"pCR"]<thr[i] & myres[,1]=="pCR" )
		cc = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="RD")
		dd = sum(myres[,"pCR"]<thr[i] & myres[,1]=="RD")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc1 = myauc/10
myauc1		## 
myx1 = myx/10
myy1 = myy/10

load(file= myoutf1)
myList[[3]] = myauc1
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx1, myy1)
names(myList)[4] = "ROC.PPM.RF"
save(myList, file=myoutf1)

#----------------------------------------------
tmp = predict(myfit, newdata = data, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)

mod1 = randomForest(class~ age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(class~age_years+er_status_ihc+clinical_t_stage+pam50_class,  sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,1], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc2 = myauc/10
myx2 = myx/10
myy2 = myy/10

load(file= myoutf1)
myList[[5]] = myauc2
names(myList)[5] = "AUC.PPM.RF"
myList[[6]] = cbind(myx2, myy2)
names(myList)[6] = "ROC.PPM.RF"
save(myList, file=myoutf1)


load(file= myoutf1)
myList[[7]] = myfit[["importance"]]
names(myList)[7] = "PPM.RF.RelativeImportance"
save(myList, file=myoutf1)


                 MeanDecreaseGini
age_years               29.515209
er_status_ihc            7.210708
clinical_t_stage         9.039587
pam50_class             15.085070


                 MeanDecreaseGini
mytf                    43.392418
age_years               27.049047
er_status_ihc            7.650928
clinical_t_stage         8.875785


#----------------------------------------------

wilcox.test(data$age_years[data$class=="Y"], data$age_years[data$class=="N"] )
mean(data$age_years[data$class=="Y"])
mean(data$age_years[data$class=="N"])

aa = sum(data$er_status_ihc=="P" & data$class=="Y", na.rm=T)
bb = sum(data$er_status_ihc=="P" & data$class=="N", na.rm=T)
cc = sum(data$er_status_ihc=="N" & data$class=="Y", na.rm=T)
dd = sum(data$er_status_ihc=="N" & data$class=="N", na.rm=T)
xx = matrix(c(aa, bb, cc, dd), 2,2)
chisq.test(xx)$p.value		##

load(file= myoutf1)
myList[[8]] = data
names(myList)[8] = "Haztis.Class.Predictors"
save(myList, file=myoutf1)

#----------------------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[!is.na(data$er_status_ihc),]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(myfit, newdata= data, type="prob")
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(row.names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx,]

#------------------------------
pred = ifelse(risk[, "pCR"]>0.5, "pCR", "RD")
response = as.character(data$pathologic_response_pcr_rd)
res = data.frame(response, pred, predictability)
myorder = order(res[, "Y"], decreasing=T)
res = res[myorder,]

load(file= myoutf1)
myList[[9]] = res
names(myList)[9] = "dis2val.accVScscore"
save(myList, file=myoutf1)


xx = ifelse(res[,1]==res[,2],1,0)
for(k in 2:length(xx))
{
	xx[k] = xx[k-1]+xx[k]
}
acc = xx/(1:length(xx)) 
plot(acc, type="b", pch=20)

w.siz = 20
step = 1
cnum = ceiling((length(acc)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(acc))
	yy[k] = mean(acc[sta:end])
}
plot(yy, pch=20, type="b")
yy




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[2.4] Hatzis_GSE25066 data  (classification)   --> validation  data (gene70)
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Hatzis_GSE25066/Clinical_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Symbol_gene70_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig3_data_gene70_validation.rda"


myList = list(NULL)
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info <- read.table(myinf1, header=T, sep="\t",  quote="")
source = c(rep("discovery", 310), rep("validation", 198))
info$source = source
info = info[!is.na(info$pathologic_response_pcr_rd), ]

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)
se = c("mytf","age_years", "er_status_ihc","clinical_t_stage", "grade", "pam50_class", "drfs_even_time_years", "drfs_1_event_0_censored","pathologic_response_pcr_rd",  "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[!is.na(data$er_status_ihc), ]

#------------------------------
pos = row.names(data)[data$pathologic_response_pcr_rd=="pCR"]
neg = row.names(data)[data$pathologic_response_pcr_rd=="RD"]
pos = pos[!is.na(pos)]
neg = neg[!is.na(neg)]

xx = data
tmp = row.names(xx)
xx = xx$mytf
names(xx) = tmp
xx = sort(xx)
xx= names(xx)%in%pos
fp = 1-xx 
tp = xx
for(j in length(xx):2)
{
	fp[j-1]= fp[j]+fp[j-1]
	tp[j-1]= tp[j]+tp[j-1]
}
fp = fp/length(neg)
tp = tp/length(pos)	
xx = c(1, fp, 0)
yy = c(1, tp, 0)
tmp1 = tmp2 = rep(0,length(xx)-1)
for(i in 1:length(tmp1))
{
	tmp1[i] = xx[i]-xx[i+1]
	tmp2[i] = (yy[i+1]+yy[i])/2	
}
myauc = sum(tmp1*tmp2)

myList[[1]] = myauc
names(myList)[1] = "AUC.gene70.Uni"
myList[[2]] = cbind(xx, yy)
names(myList)[2] = "ROC.gene70.Uni"
save(myList, file=myoutf1)


library(randomForest)
xx = table(data$pathologic_response_pcr_rd)
sam.siz = min(xx)
myfit <- randomForest(pathologic_response_pcr_rd~mytf+age_years+er_status_ihc+clinical_t_stage, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
myfit
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)

pdat = data[data$pathologic_response_pcr_rd=="pCR", ]
ndat = data[data$pathologic_response_pcr_rd=="RD", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(pathologic_response_pcr_rd~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"pathologic_response_pcr_rd"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="pCR")
		bb = sum(myres[,"pCR"]<thr[i] & myres[,1]=="pCR" )
		cc = sum(myres[,"pCR"]>=thr[i] & myres[,1]=="RD")
		dd = sum(myres[,"pCR"]<thr[i] & myres[,1]=="RD")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc1 = myauc/10
myauc1		## 
myx1 = myx/10
myy1 = myy/10

load(file= myoutf1)
myList[[3]] = myauc1
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx1, myy1)
names(myList)[4] = "ROC.PPM.RF"
save(myList, file=myoutf1)

#----------------------------------------------
tmp = predict(myfit, newdata = data, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$pathologic_response_pcr_rd, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)

mod1 = randomForest(class~ age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)
myauc = 0
myx = myy = rep(0, 101)

for(p in 1:10)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
		sam.siz = min(nrow(ptr), nrow(ntr))
		fit <- randomForest(class~age_years+er_status_ihc+clinical_t_stage+pam50_class, sampsize=c(sam.siz, sam.siz), replace=T, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,1], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc2 = myauc/10
myx2 = myx/10
myy2 = myy/10


load(file= myoutf1)
myList[[5]] = myauc2
names(myList)[5] = "AUC.PPM.RF"
myList[[6]] = cbind(myx2, myy2)
names(myList)[6] = "ROC.PPM.RF"
save(myList, file=myoutf1)


load(file= myoutf1)
myList[[7]] = myfit[["importance"]]
names(myList)[7] = "PPM.RF.RelativeImportance"
save(myList, file=myoutf1)

                 MeanDecreaseGini
age_years               30.586575
er_status_ihc            4.686668
clinical_t_stage         6.749625
pam50_class             13.981314


                 MeanDecreaseGini
mytf                    29.021895
age_years               22.498751
er_status_ihc            2.875848
clinical_t_stage         5.671855

#----------------------------------------------

wilcox.test(data$age_years[data$class=="Y"], data$age_years[data$class=="N"] )

mean(data$age_years[data$class=="Y"])
mean(data$age_years[data$class=="N"])

aa = sum(data$er_status_ihc=="P" & data$class=="Y", na.rm=T)
bb = sum(data$er_status_ihc=="P" & data$class=="N", na.rm=T)
cc = sum(data$er_status_ihc=="N" & data$class=="Y", na.rm=T)
dd = sum(data$er_status_ihc=="N" & data$class=="N", na.rm=T)
xx = matrix(c(aa, bb, cc, dd), 2,2)
chisq.test(xx)$p.value		##

load(file= myoutf1)
myList[[8]] = data
names(myList)[8] = "Haztis.Class.Predictors"
save(myList, file=myoutf1)

#----------------------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[!is.na(data$er_status_ihc),]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(myfit, newdata= data, type="prob")
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(row.names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx,]

#------------------------------
pred = ifelse(risk[, "pCR"]>0.5, "pCR", "RD")
response = as.character(data$pathologic_response_pcr_rd)
res = data.frame(response, pred, predictability)
myorder = order(res[, "Y"], decreasing=T)
res = res[myorder,]

load(file= myoutf1)
myList[[9]] = res
names(myList)[9] = "val2dis.accVScscore"
save(myList, file=myoutf1)


xx = ifelse(res[,1]==res[,2],1,0)
for(k in 2:length(xx))
{
	xx[k] = xx[k-1]+xx[k]
}
acc = xx/(1:length(xx)) 
plot(acc, type="b", pch=20)

w.siz = 20
step = 1
cnum = ceiling((length(acc)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(acc))
	yy[k] = mean(acc[sta:end])
}
plot(yy, pch=20, type="b")
yy














#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[3] Figure 4
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.1] collect all datasets
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Hatzis_GSE25066/Clinical_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Breast_neoadjuvant_Symbol_oncotypeDX_result.txt"
myinf3 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Breast_neoadjuvant_Symbol_E2F4_iRAS.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Hatzis_GSE25066_Symbol_gene70_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig4_Haztis_models_and_4otherdatasets.rda"

myList = list(NULL)
dat1 <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
dat2 <- read.table(myinf3, header=T, sep="\t", row.names=1, quote="")
dat3 <- read.table(myinf4, header=T, sep="\t", row.names=1, quote="")
data = cbind(dat1[, "score"], dat3[, "score"], dat2[, "all.intersection.ES"])
row.names(data) = row.names(dat1)
colnames(data) = c("oncotype", "gene70", "E2F4")

info <- read.table(myinf1, header=T, sep="\t",  quote="")
source = c(rep("discovery", 310), rep("validation", 198))
info$source = source
info = info[!is.na(info$pathologic_response_pcr_rd), ]
xx = as.character(info$grade)
xx[grep("Indeterminate", xx)]= "4"
info$grade = as.numeric(xx)

comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
se = c("age_years", "er_status_ihc","clinical_t_stage", "grade", "pam50_class", "drfs_even_time_years", "drfs_1_event_0_censored","pathologic_response_pcr_rd",  "source")
info = info[,se]
data = cbind(data, info)

myList[[1]] = data
names(myList)[1] = "Haztis.data"

#-------------------------------------
myinf1 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/USO_GSE23988_E2F4_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Iwamoto6dat/USO_GSE23988/Clinical_info.txt"
myinf3 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/USO_GSE23988_oncotypeDX_result.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/USO_GSE23988_Gene70_result.txt"
myinf5 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/USO_GSE23988_PAM50_result.txt"

onco = read.table(myinf3, sep="\t", header=T, row.names=1)
gen70 = read.table(myinf4, sep="\t", header=T, row.names=1)
data = read.table(myinf1, sep="\t", header=T, row.names=1)
info = read.table(myinf2, sep="\t", header=T, row.names=1)
PAM50 = read.table(myinf5, sep="\t", header=T, row.names=1)

comGen = intersect(row.names(data), row.names(info))
data = data[comGen,]
info = info[comGen,]
onco = onco[comGen,1]
gen70 = gen70[comGen,1]
PAM50 = PAM50[comGen,1]

mytf = data[,"all.intersection.ES"]
data =  cbind(PAM50, onco, gen70, mytf, info)
colnames(data)[1:4] = c("PAM50", "oncotype", "gene70", "E2F4")

myList[[2]] = data
names(myList)[2] = "USO.GSE23988.data"

#-------------------------------------
myinf1 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Desmedt_GSE16446_E2F4_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Desmedt_GSE16446/Clinical_info.txt"
myinf3 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Desmedt_GSE16446_oncotypeDX_result.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Desmedt_GSE16446_gene70_result.txt"
myinf5 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Desmedt_GSE16446_PAM50_result.txt"

onco = read.table(myinf3, sep="\t", header=T, row.names=1)
gen70 = read.table(myinf4, sep="\t", header=T, row.names=1)
data = read.table(myinf1, sep="\t", header=T, row.names=1)
info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
PAM50 = read.table(myinf5, sep="\t", header=T, row.names=1)

comGen = intersect(row.names(data), row.names(info))
data = data[comGen,]
info = info[comGen,]
onco = onco[comGen,1]
gen70 = gen70[comGen,1]
PAM50 = PAM50[comGen,1]

mytf = data[,"all.intersection.ES"]
data =  cbind(PAM50, onco, gen70, mytf, info)
colnames(data)[1:4] = c("PAM50", "oncotype", "gene70", "E2F4")

myList[[3]] = data
names(myList)[3] = "Desmedt.GSE16446.data"

#-------------------------------------
myinf1 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE20194_E2F4_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/MAQC_GSE20194/Clinical_info.txt"
myinf3 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE20194_oncotypeDX_result.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE20194_gene70_result.txt"
myinf5 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE20194_PAM50_result.txt"

onco = read.table(myinf3, sep="\t", header=T, row.names=1)
gen70 = read.table(myinf4, sep="\t", header=T, row.names=1)
data = read.table(myinf1, sep="\t", header=T, row.names=1)
info = read.table(myinf2, sep="\t", header=T, row.names=2, quote="")
PAM50 = read.table(myinf5, sep="\t", header=T, row.names=1)

mytf = data[,"all.intersection.ES"]
data =  cbind(PAM50, onco, gen70, mytf, info)
colnames(data)[1:4] = c("PAM50", "oncotype", "gene70", "E2F4")


myList[[4]] = data
names(myList)[4] = "MAQC_GSE20194.data"


#-------------------------------------
myinf1 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE22093_E2F4_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Iwamoto6dat/MAQC_GSE22093/Clinical_info.txt"
myinf3 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE22093_oncotypeDX_result.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE22093_Gene70_result.txt"
myinf5 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/MAQC_GSE22093_PAM50_result.txt"

onco = read.table(myinf3, sep="\t", header=T, row.names=1)
gen70 = read.table(myinf4, sep="\t", header=T, row.names=1)
data = read.table(myinf1, sep="\t", header=T, row.names=1)
info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
PAM50 = read.table(myinf5, sep="\t", header=T, row.names=1)

comGen = intersect(row.names(data), row.names(info))
data = data[comGen,]
info = info[comGen,]
onco = onco[comGen,1]
gen70 = gen70[comGen,1]
PAM50 = PAM50[comGen,1]

mytf = data[,"all.intersection.ES"]
data =  cbind(PAM50, onco, gen70, mytf, info)
colnames(data)[1:4] = c("PAM50", "oncotype", "gene70", "E2F4")

myList[[5]] = data
names(myList)[5] = "MAQC.GSE22093.data"

#-------------------------------------
myinf1 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Miyake_GSE32646_E2F4_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Miyake_GSE32646/Clinical_info.txt"
myinf3 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Miyake_GSE32646_oncotypeDX_result.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Miyake_GSE32646_gene70_result.txt"
myinf5 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Miyake_GSE32646_PAM50_result.txt"

onco = read.table(myinf3, sep="\t", header=T, row.names=1)
gen70 = read.table(myinf4, sep="\t", header=T, row.names=1)
data = read.table(myinf1, sep="\t", header=T, row.names=1)
row.names(data) = row.names(onco)
info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
PAM50 = read.table(myinf5, sep="\t", header=T, row.names=1)
row.names(PAM50) = row.names(info)

comGen = intersect(row.names(data), row.names(info))
data = data[comGen,]
info = info[comGen,]
onco = onco[comGen,1]
gen70 = gen70[comGen,1]
PAM50 = PAM50[comGen,1]

mytf = data[,"all.intersection.ES"]
data =  cbind(PAM50, onco, gen70, mytf, info)
colnames(data)[1:4] = c("PAM50", "oncotype", "gene70", "E2F4")

myList[[6]] = data
names(myList)[6] = "Miyake_GSE32646.data"


#-------------------------------------
myinf1 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Tabchy_GSE20271_E2F4_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Tabchy_GSE20271/Clinical_info.txt"
myinf3 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Tabchy_GSE20271_oncotypeDX_result.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Tabchy_GSE20271_Gene70_result.txt"
myinf5 = "/lorax/chenglab/cc59/WorSpa/m1_cancer/Breast/data/Tabchy_GSE20271_PAM50_result.txt"

onco = read.table(myinf3, sep="\t", header=T, row.names=1)
gen70 = read.table(myinf4, sep="\t", header=T, row.names=1)
data = read.table(myinf1, sep="\t", header=T, row.names=1)
info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
PAM50 = read.table(myinf5, sep="\t", header=T, row.names=1)

comGen = intersect(row.names(data), row.names(info))
data = data[comGen,]
info = info[comGen,]
onco = onco[comGen,1]
gen70 = gen70[comGen,1]
PAM50 = PAM50[comGen,1]

mytf = data[,"all.intersection.ES"]
data =  cbind(PAM50, onco, gen70, mytf, info)
colnames(data)[1:4] = c("PAM50", "oncotype", "gene70", "E2F4")

myList[[7]] = data
names(myList)[7] = "Tabchy.GSE20271.data"

names(myList)
save(myList, file=myoutf1)



colnames(myList[[1]])
colnames(myList[[2]])
colnames(myList[[3]])
colnames(myList[[4]])
colnames(myList[[5]])
colnames(myList[[6]])
colnames(myList[[7]])


> colnames(myList[[1]])
 [1] "oncotype"                   "gene70"                    
 [3] "E2F4"                       "age_years"                 
 [5] "er_status_ihc"              "clinical_t_stage"          
 [7] "grade"                      "pam50_class"               
 [9] "drfs_even_time_years"       "drfs_1_event_0_censored"   
[11] "pathologic_response_pcr_rd" "source"                    
> colnames(myList[[2]])
 [1] "oncotype"              "gene70"                "E2F4"                 
 [4] "Title"                 "pcr.v.rd"              "array.qc"             
 [7] "er.status"             "age"                   "grade"                
[10] "prechemo.tumor.size"   "prechemo.t.stage"      "prechemo.nodal.status"
> colnames(myList[[3]])
 [1] "oncotype"       "gene70"         "E2F4"           "Title"         
 [5] "agebin"         "t"              "n"              "grade"         
 [9] "her2fish"       "her2fishbin"    "top2atri"       "topoihc"       
[13] "esr1bimod"      "erbb2bimod"     "final_analysis" "pcr"           
[17] "dmfs_event"     "dmfs_time"      "os_event"       "os_time"       
> colnames(myList[[4]])
 [1] "oncotype"               "gene70"                 "E2F4"                  
 [4] "title"                  "CEL.file"               "source.name"           
 [7] "organism"               "age"                    "race"                  
[10] "ER_status"              "pCR_vs_RD"              "PR_status"             
[13] "molecule"               "label"                  "description"           
[16] "platform"               "Additional.information" "Tbefore"               
[19] "Nbefore"                "BMNgrd"                 "ER"                    
[22] "HER2.Status"            "Her2.IHC"               "Her2.FISH"             
[25] "Histology"              "Treatment.Code"         "Treatments.Comments"   
> colnames(myList[[5]])
 [1] "oncotype"                "gene70"                 
 [3] "E2F4"                    "Title"                  
 [5] "tissue"                  "pcr.v.rd"               
 [7] "array.qc"                "p53.status"             
 [9] "er.expression"           "er.immunohistochemistry"
[11] "age"                     "prechemo.t"             
[13] "prechemo.tumor.size"     "prechemo.n"             
[15] "bmn.grade"              
> colnames(myList[[6]])
 [1] "oncotype"            "gene70"              "E2F4"               
 [4] "Title"               "tissue"              "age"                
 [7] "clinical.t.stage"    "lymph.node.status"   "clinical.stage"     
[10] "histological.grade"  "er.status.ihc"       "pr.status.ihc"      
[13] "her2.status.fish"    "pathologic.response"
> colnames(myList[[7]])
 [1] "oncotype"                "gene70"                 
 [3] "E2F4"                    "Title"                  
 [5] "pcr.or.rd"               "array.qc"               
 [7] "dlda30.score"            "dlda30.pred.pcr1"       
 [9] "biopsy.date.mmddyy"      "age"                    
[11] "race"                    "histology"              
[13] "prechemo.t"              "prechemo.n"             
[15] "bmn.grade"               "er.positive"            
[17] "er.status"               "pr.positive"            
[19] "pr.status"               "her2.status"            
[21] "her2.ihc"                "her2.fish"              
[23] "preoperative.treatment"  "treatment.received.fac1"
[25] "randomized.fac1"         "surgery.type"           
[27] "surgery.date"            "post.chemo.size.cm"     
[29] "post.chemo..ln.total"   



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.2] created matched datasets (Haztis and 3 others)

rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig4_Haztis_models_and_4otherdatasets.rda"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig4_Haztis_models_and_3datasets_matched.rda"

load(file=myinf1)
names(myList)
##[1] "Haztis.data"           "USO.GSE23988.data"     "Desmedt.GSE16446.data"
## [4] "MAQC_GSE20194.data"    "MAQC.GSE22093.data"    "Miyake_GSE32646.data"   
## [7] "Tabchy.GSE20271.data" 

#------------------------------
data = myList[[1]]
se = c("oncotype", "gene70", "E2F4", "age_years", "er_status_ihc", "clinical_t_stage", "grade", "pathologic_response_pcr_rd", "source", "pam50_class")
data = data[,se]
colnames(data) = c("oncotype", "gene70", "E2F4", "age", "ER", "stage", "grade", "PCR.RD", "source", "PAM50")
myList[[1]] = data

data = myList[[2]]
se = c("oncotype", "gene70", "E2F4", "age", "er.status", "prechemo.t.stage", "grade", "pcr.v.rd", "PAM50")
data = data[,se]
colnames(data) = c("oncotype", "gene70", "E2F4", "age", "ER", "stage", "grade", "PCR.RD", "PAM50")
data$ER = ifelse(data$ER=="ERpos", "P", "N")
data$stage = paste("T", data$stage, sep="")
myList[[2]] = data


data = myList[[4]]
se = c("oncotype", "gene70", "E2F4", "age", "ER_status", "Tbefore", "BMNgrd", "pCR_vs_RD", "PAM50")
data = data[,se]
colnames(data) = c("oncotype", "gene70", "E2F4", "age", "ER", "stage", "grade", "PCR.RD", "PAM50")
data$stage = ifelse(!is.na(data$stage), paste("T", data$stage, sep=""), NA)
myList[[4]] = data



data = myList[[5]]
se = c("oncotype", "gene70", "E2F4", "age", "er.immunohistochemistry", "prechemo.t", "bmn.grade", "pcr.v.rd", "PAM50")
data = data[,se]
colnames(data) = c("oncotype", "gene70", "E2F4", "age", "ER", "stage", "grade", "PCR.RD", "PAM50")
data$ER = ifelse(data$ER=="ERpos", "P", "N")
data$stage = ifelse(!is.na(data$stage), paste("T", data$stage, sep=""), NA)
myList[[5]] = data


data = myList[[6]]
se = c("oncotype", "gene70", "E2F4", "age", "er.status.ihc", "clinical.t.stage", "histological.grade", "pathologic.response", "PAM50")
data = data[,se]
colnames(data) = c("oncotype", "gene70", "E2F4", "age", "ER", "stage", "grade", "PCR.RD", "PAM50")
data$ER = ifelse(data$ER=="positive", "P", "N")
xx = as.integer(as.character(data$stage))
data$stage = ifelse(!is.na(xx), paste("T", xx, sep=""), NA)
data$PCR.RD = ifelse(data$PCR.RD=="pCR", "pCR", "RD")
myList[[6]] = data

myList[[7]] =NULL
myList[[3]] =NULL


names(myList)

[1] "Haztis.data"          "USO.GSE23988.data"    "MAQC_GSE20194.data"  
[4] "MAQC.GSE22093.data"   "Miyake_GSE32646.data"

save(myList, file=myoutf1)



#---------------------------------------------------


library(survival)
library(randomForest)
#------------------------------
data = myList[[1]]
se = which(colnames(data)=="grade")
data = data[,-se]
data = data[!is.na(data$ER), ]
raw.data = data

data = raw.data
xx = table(data$PCR.RD)
sam.siz = min(xx)
myfit <- randomForest(PCR.RD~oncotype+age+ER+stage,  sampsize=c(sam.siz, sam.siz), replace=T,  data=data, ntree=10000)
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$PCR.RD, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)
mod.all = randomForest(class~ age+ER+stage,  sampsize=c(sam.siz, sam.siz), replace=T,  data=data, ntree=10000)
tmp = predict(mod.all, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)


data = raw.data
data = data[data$source=="discovery",]
xx = table(data$PCR.RD)
sam.siz = min(xx)
myfit <- randomForest(PCR.RD~oncotype+age+ER+stage,  sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$PCR.RD, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)
mod.dis = randomForest(class~ age+ER+stage,  sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod.dis, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)
## this is the best model

data = raw.data
data = data[data$source=="validation",]
xx = table(data$PCR.RD)
sam.siz = min(xx)
myfit <- randomForest(PCR.RD~oncotype+age+ER+stage,  sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$PCR.RD, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)
mod.val = randomForest(class~ age+ER+stage,  sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod.val, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[3.3] analysis
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig4_Haztis_models_and_3datasets_matched.rda"

load(file=myinf1)

library(survival)
library(randomForest)
#------------------------------
data = myList[[1]]
se = which(colnames(data)=="grade")
data = data[,-se]
data = data[!is.na(data$ER), ]
raw.data = data


data = raw.data
data = data[data$source=="discovery",]
xx = table(data$PCR.RD)
sam.siz = min(xx)
myfit <- randomForest(PCR.RD~oncotype+age+ER+stage, sampsize=c(sam.siz, sam.siz), replace=T,  data=data, ntree=10000)
tmp = predict(myfit, type="prob")
xx = ifelse(tmp[,"pCR"]>0.5, "pCR", "RD")
class = as.factor(ifelse(xx==data$PCR.RD, "Y", "N"))
data = cbind(class, data)
xx = table(data$class)
sam.siz = min(xx)
mod1 = randomForest(class~ age+ER+stage+PAM50,  sampsize=c(sam.siz, sam.siz), replace=T, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)
## this is the best model


#----------------------------
data = myList[[5]]
se = c("oncotype" , "age", "ER", "stage", "PCR.RD", "PAM50")
data = data[,se]
data$ER = as.factor(data$ER)
levels(data$ER) = levels(myList[[1]]$ER)
data$stage = as.factor(data$stage)
levels(data$stage) = levels(myList[[1]]$stage)
data = data[!is.na(data$PCR.RD),]


predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(myfit, newdata= data, type="prob")
risk = risk[!is.na(risk[,"pCR"]),]
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(row.names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx,]
risk = risk[comxx,]

#------------------------------
pred = ifelse(risk[, "pCR"]>0.5, "pCR", "RD")
response = as.character(data$PCR.RD)
res = cbind(response, pred)
myorder = order(predictability[, "Y"], decreasing=T)
res = res[myorder,]

xx = ifelse(res[,1]==res[,2],1,0)
for(k in 2:length(xx))
{
	xx[k] = xx[k-1]+xx[k]
}
acc = xx/(1:length(xx)) 

w.siz = 20
step = 1
cnum = ceiling((length(acc)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(acc))
	yy[k] = mean(acc[sta:end])
}
plot(yy, pch=20, type="b")
yy




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4] Figure 5     -- cox model (Curtis disovery <-> validation)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.1] Curtis discovery data  (Cox model) --> data
rm(list=ls())
myinf1A = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/discovery_Patient_tumor_info.txt"
myinf1B = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/validation_Patient_tumor_info.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/Curtis_oncotypeDX_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig5_Curtis_oncotypeDX_discovery_Cox.rda"

myList = list(NULL)

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)

xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="discovery",]

mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance

xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd
mycox$concordance

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)			#  this figure indicates the accuracy vary dramatically between samples
abline(v=50, col=2)

myList[[1]] = summary(mycox)$concordance
names(myList)[1] = "Cox.CI"
myList[[2]] = res
names(myList)[2] = "Sample.Specific.Accuracy"
save(myList, file= myoutf1)


#------
se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]

library(randomForest)
#----------------------------------------------
data = raw.data
data = data[data$source=="discovery",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1]
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)


#----------------------------------------------
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)

myauc = 0
myx = myy = rep(0, 101)

for(p in 1:5)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
	fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"class"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc = myauc/5
myx = myx/5
myy = myy/5

load(file= myoutf1)
myList[[3]] = myauc
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx, myy)
names(myList)[4] = "ROC.PPM.RF"
myList[[5]] = mod1[["importance"]]
names(myList)[5] = "PPM.RF.RelativeImportance"
myList[[6]] = data
names(myList)[6] = "Class.and.Predictors"
save(myList, file=myoutf1)


#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx, ]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## positive correlation 

load(file= myoutf1)
myList[[7]] = cv.res
names(myList)[7] = "Predicted.Validation.Result"
save(myList, file=myoutf1)


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]

w.siz = 200
step = 20
cnum = ceiling((nrow(tmp)-w.siz)/step)
yy1 = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(yy1))
	tmpxx = tmp[sta:end,]
	yy1[k] = (sum(tmpxx[, "concordance"]) + 0.5*sum(tmpxx[, "tie"]))/sum(tmpxx[, "count"])
}
plot(yy1, pch=20, type="b")


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]			## check here
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")
abline(h=cum.con[length(cum.con)], col=2)		## smoothed curves


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.2] Curtis validation data  (Cox model) --> data oncotype DX
rm(list=ls())
myinf1A = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/discovery_Patient_tumor_info.txt"
myinf1B = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/validation_Patient_tumor_info.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/Curtis_oncotypeDX_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig5_Curtis_oncotypeDX_validation_Cox.rda"

myList = list(NULL)

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)

xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="validation",]

mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance

xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd
mycox$concordance

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)			#  this figure indicates the accuracy vary dramatically between samples
abline(v=50, col=2)

myList[[1]] = summary(mycox)$concordance
names(myList)[1] = "Cox.CI"
myList[[2]] = res
names(myList)[2] = "Sample.Specific.Accuracy"
save(myList, file= myoutf1)


#------
se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]

library(randomForest)
#----------------------------------------------
data = raw.data
data = data[data$source=="validation",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1]
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)


#----------------------------------------------
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)

myauc = 0
myx = myy = rep(0, 101)

for(p in 1:5)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
	fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"class"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc = myauc/5
myx = myx/5
myy = myy/5
##  AUC = 0.8058482

load(file= myoutf1)
myList[[3]] = myauc
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx, myy)
names(myList)[4] = "ROC.PPM.RF"
myList[[5]] = mod1[["importance"]]
names(myList)[5] = "PPM.RF.RelativeImportance"
myList[[6]] = data
names(myList)[6] = "Class.and.Predictors"
save(myList, file=myoutf1)


#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx, ]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## positive correlation 

load(file= myoutf1)
myList[[7]] = cv.res
names(myList)[7] = "Predicted.Discovery.Result"
save(myList, file=myoutf1)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.3] Curtis discovery data  (Cox model) --> data  (gene70)
rm(list=ls())
myinf1A = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/discovery_Patient_tumor_info.txt"
myinf1B = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/validation_Patient_tumor_info.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/Curtis_gene70_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig5_Curtis_gene70_discovery_Cox.rda"

myList = list(NULL)

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)

xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="discovery",]

mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance

xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd
mycox$concordance

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)			#  this figure indicates the accuracy vary dramatically between samples
abline(v=50, col=2)

myList[[1]] = summary(mycox)$concordance
names(myList)[1] = "Cox.CI"
myList[[2]] = res
names(myList)[2] = "Sample.Specific.Accuracy"
save(myList, file= myoutf1)


#------
se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]

library(randomForest)
#----------------------------------------------
data = raw.data
data = data[data$source=="discovery",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1]
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)


#----------------------------------------------
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)

myauc = 0
myx = myy = rep(0, 101)

for(p in 1:5)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
	fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"class"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc = myauc/5
myx = myx/5
myy = myy/5

load(file= myoutf1)
myList[[3]] = myauc
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx, myy)
names(myList)[4] = "ROC.PPM.RF"
myList[[5]] = mod1[["importance"]]
names(myList)[5] = "PPM.RF.RelativeImportance"
myList[[6]] = data
names(myList)[6] = "Class.and.Predictors"
save(myList, file=myoutf1)


#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx, ]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## positive correlation 

load(file= myoutf1)
myList[[7]] = cv.res
names(myList)[7] = "Predicted.Validation.Result"
save(myList, file=myoutf1)


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]

w.siz = 200
step = 20
cnum = ceiling((nrow(tmp)-w.siz)/step)
yy1 = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(yy1))
	tmpxx = tmp[sta:end,]
	yy1[k] = (sum(tmpxx[, "concordance"]) + 0.5*sum(tmpxx[, "tie"]))/sum(tmpxx[, "count"])
}
plot(yy1, pch=20, type="b")


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]			## check here
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")
abline(h=cum.con[length(cum.con)], col=2)		## smoothed curves


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.4] Curtis validation data  (Cox model) --> data gene70
rm(list=ls())
myinf1A = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/discovery_Patient_tumor_info.txt"
myinf1B = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/validation_Patient_tumor_info.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/Curtis_gene70_result.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig5_Curtis_gene70_validation_Cox.rda"

myList = list(NULL)

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "score"]
data = cbind(mytf, info)

xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="validation",]

mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance

xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd
mycox$concordance

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)			#  this figure indicates the accuracy vary dramatically between samples
abline(v=50, col=2)

myList[[1]] = summary(mycox)$concordance
names(myList)[1] = "Cox.CI"
myList[[2]] = res
names(myList)[2] = "Sample.Specific.Accuracy"
save(myList, file= myoutf1)


#------
se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]

library(randomForest)
#----------------------------------------------
data = raw.data
data = data[data$source=="validation",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1]
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)


#----------------------------------------------
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)

myauc = 0
myx = myy = rep(0, 101)

for(p in 1:5)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
	fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"class"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc = myauc/5
myx = myx/5
myy = myy/5
##  AUC = 0.8058482

load(file= myoutf1)
myList[[3]] = myauc
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx, myy)
names(myList)[4] = "ROC.PPM.RF"
myList[[5]] = mod1[["importance"]]
names(myList)[5] = "PPM.RF.RelativeImportance"
myList[[6]] = data
names(myList)[6] = "Class.and.Predictors"
save(myList, file=myoutf1)


#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx, ]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## positive correlation 

load(file= myoutf1)
myList[[7]] = cv.res
names(myList)[7] = "Predicted.Discovery.Result"
save(myList, file=myoutf1)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.5] Curtis discovery data  (Cox model) --> data  (E2F4)
rm(list=ls())
myinf1A = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/discovery_Patient_tumor_info.txt"
myinf1B = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/validation_Patient_tumor_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/w4_surv/result/Breast/Breast_Curtis_Symbol_iRAS.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig5_Curtis_E2F4_discovery_Cox.rda"

myList = list(NULL)

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "all.intersection.ES"]
data = cbind(mytf, info)

xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="discovery",]

mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance

xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd
mycox$concordance

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)			#  this figure indicates the accuracy vary dramatically between samples
abline(v=50, col=2)

myList[[1]] = summary(mycox)$concordance
names(myList)[1] = "Cox.CI"
myList[[2]] = res
names(myList)[2] = "Sample.Specific.Accuracy"
save(myList, file= myoutf1)


#------
se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]

library(randomForest)
#----------------------------------------------
data = raw.data
data = data[data$source=="discovery",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1]
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)


#----------------------------------------------
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)

myauc = 0
myx = myy = rep(0, 101)

for(p in 1:5)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
	fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"class"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc = myauc/5
myx = myx/5
myy = myy/5

load(file= myoutf1)
myList[[3]] = myauc
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx, myy)
names(myList)[4] = "ROC.PPM.RF"
myList[[5]] = mod1[["importance"]]
names(myList)[5] = "PPM.RF.RelativeImportance"
myList[[6]] = data
names(myList)[6] = "Class.and.Predictors"
save(myList, file=myoutf1)


#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="validation",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx, ]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## positive correlation 

load(file= myoutf1)
myList[[7]] = cv.res
names(myList)[7] = "Predicted.Validation.Result"
save(myList, file=myoutf1)


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]

w.siz = 200
step = 20
cnum = ceiling((nrow(tmp)-w.siz)/step)
yy1 = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(yy1))
	tmpxx = tmp[sta:end,]
	yy1[k] = (sum(tmpxx[, "concordance"]) + 0.5*sum(tmpxx[, "tie"]))/sum(tmpxx[, "count"])
}
plot(yy1, pch=20, type="b")


#--------------------------------------------------------------
xx = cv.res[cv.res[, "count"]>nrow(data)*0.2,]			## check here
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")
abline(h=cum.con[length(cum.con)], col=2)		## smoothed curves


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[4.6] Curtis validation data  (Cox model) --> data E2F4
rm(list=ls())
myinf1A = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/discovery_Patient_tumor_info.txt"
myinf1B = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Curtis/validation_Patient_tumor_info.txt"
myinf2 = "/lorax/chenglab/cc59/WorSpa/w4_surv/result/Breast/Breast_Curtis_Symbol_iRAS.txt"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig5_Curtis_E2F4_validation_Cox.rda"

myList = list(NULL)

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
#------------------------------
info.A <- read.table(myinf1A, header=T, sep="\t",  quote="")
info.B <- read.table(myinf1B, header=T, sep="\t",  quote="")
source = c(rep("discovery", nrow(info.A)), rep("validation", nrow(info.B)))
info = rbind(info.A, info.B)
info = cbind(info, source)
info = info[!is.na(info$last_follow_up_status),]
info = info[!is.na(info$T),]
info[,1] = as.character(info[,1])
xx = info[,1]
xx = gsub("-", ".", xx)
info[,1] = xx
row.names(info) = xx
xx = rep(0, nrow(info))
xx[info[, "last_follow_up_status"]=="d-d.s."] = 1
e.rfs = xx
t.rfs = as.numeric(info[, "T"])
info = cbind(t.rfs, e.rfs, info)

#------------------------------
comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
mytf = data[, "all.intersection.ES"]
data = cbind(mytf, info)

xx = data$lymph_nodes_positive
xx = ifelse(xx>0, 1, 0)
data$lymph_nodes_positive = xx
xx = apply(is.na(data), 2, sum)
data = data[, xx<500]
se = grep("CT", data$Treatment)
CT = rep("no", nrow(data))
CT[se] = "yes"
se = grep("RT", data$Treatment)
RT = rep("no", nrow(data))
RT[se] = "yes"
se = grep("HT", data$Treatment)
HT = rep("no", nrow(data))
HT[se] = "yes"
data = cbind(CT, RT, HT, data)
se = c("CT", "RT","HT", "mytf","t.rfs", "e.rfs", "age_at_diagnosis", "menopausal_status_inferred", "grade","size", "stage","lymph_nodes_positive", "NPI", "histological_type","ER_IHC_status", "cellularity", "Pam50Subtype", "ER.Expr", "Her2.Expr","PR.Expr", "source")
data = data[,se]
raw.data = data


library(survival)
#------------------------------
data = raw.data
data = data[data$source=="validation",]

mycox = coxph(Surv(t.rfs, e.rfs)~mytf, data) 
summary(mycox)$concordance

xx = cbind(data$t.rfs, mycox$linear.predictors, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd
mycox$concordance

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

plot(count, accuracy, pch=20)			#  this figure indicates the accuracy vary dramatically between samples
abline(v=50, col=2)

myList[[1]] = summary(mycox)$concordance
names(myList)[1] = "Cox.CI"
myList[[2]] = res
names(myList)[2] = "Sample.Specific.Accuracy"
save(myList, file= myoutf1)


#------
se = which(res[, "count"]>nrow(data)*0.2)
res = res[se,]

library(randomForest)
#----------------------------------------------
data = raw.data
data = data[data$source=="validation",]
comxx = intersect(row.names(data), row.names(res))
acc = res[comxx, "accuracy"]
thr = summary(mycox)$concordance[1]
plot(density(acc))
abline(v=thr)
sum(acc>=thr)
class = ifelse(acc>=thr, "Y", "N")
data = cbind(class, data[comxx,])
data=rfImpute(class~., data=data)
data = cbind(acc, data)


#----------------------------------------------
mod1 = randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=data, ntree=10000)
tmp = predict(mod1, type="prob")
xx = ifelse(tmp[,"Y"]>0.5, "Y", "N")
sum(xx==data$class)
sum(xx==data$class)/nrow(data)

## 10-fold CV

pdat = data[data$class=="Y", ]
ndat = data[data$class=="N", ]
kk = 10
psiz = floor(nrow(pdat)/kk)
nsiz = floor(nrow(ndat)/kk)

myauc = 0
myx = myy = rep(0, 101)

for(p in 1:5)
{
	myres = NULL
	pdat = pdat[sample(1:nrow(pdat)),]
	ndat = ndat[sample(1:nrow(ndat)),]
	for(k in 1:kk)
	{
		cat("\r\r\r", p, "-->", k)
		if(k==kk)
		{
			se = ((k-1)*psiz+1):nrow(pdat)
		}else
		{
			se = ((k-1)*psiz+1):(k*psiz)	
		}
		ptr = pdat[-se,]	
		pte = pdat[se,]
		if(k==kk)
		{
			se = ((k-1)*nsiz+1):nrow(ndat)
		}else
		{
			se = ((k-1)*nsiz+1):(k*nsiz)	
		}
		ntr = ndat[-se,]	
		nte = ndat[se,]
		
		tr = rbind(ptr, ntr)
		te = rbind(pte, nte)
	fit <- randomForest(class~ age_at_diagnosis + grade + size + stage + lymph_nodes_positive + Pam50Subtype + ER.Expr + Her2.Expr + CT + RT + HT, data=tr, ntree=10000)
		tmp = predict(fit, te[,-1], type="prob")
		tmp = data.frame(te[,"class"], tmp)
		myres = rbind(myres, tmp)
	}


	thr = (1:99)*0.01
	yy =  xx =  rep(0, length(thr))
	fdr = rep(0,99)
	for(i in 1:length(thr))
	{
		aa = sum(myres[,"Y"]>=thr[i] & myres[,1]=="Y")
		bb = sum(myres[,"Y"]<thr[i] & myres[,1]=="Y" )
		cc = sum(myres[,"Y"]>=thr[i] & myres[,1]=="N")
		dd = sum(myres[,"Y"]<thr[i] & myres[,1]=="N")
		fdr[i] = aa/sum(myres[,2]>=thr[i])
		yy[i] = aa/(aa+bb)
		xx[i] = cc/(cc+dd)
	}
	xx = c(1, xx, 0)
	yy = c(1, yy, 0)
	tmp1 = tmp2 = rep(0,100)
	for(i in 1:100)
	{
		tmp1[i] = xx[i]-xx[i+1]
		tmp2[i] = (yy[i+1]+yy[i])/2	
	}
	myauc = myauc + sum(tmp1*tmp2)
	myx = myx+xx
	myy = myy+yy
}

myauc = myauc/5
myx = myx/5
myy = myy/5
##  AUC = 0.8058482

load(file= myoutf1)
myList[[3]] = myauc
names(myList)[3] = "AUC.PPM.RF"
myList[[4]] = cbind(myx, myy)
names(myList)[4] = "ROC.PPM.RF"
myList[[5]] = mod1[["importance"]]
names(myList)[5] = "PPM.RF.RelativeImportance"
myList[[6]] = data
names(myList)[6] = "Class.and.Predictors"
save(myList, file=myoutf1)


#----------------------------------------------
#------------------------------
data = raw.data
data = data[data$source=="discovery",]
data = data[data$Pam50Subtype!="NC", ]
predictability = predict(mod1, newdata=data, type="prob")		## predictability score
risk = predict(mycox, newdata= data, na.action=na.exclude)
comxx = intersect(row.names(predictability), row.names(data))
comxx = intersect(names(risk), comxx)
data = data[comxx,]
predictability = predictability[comxx, ]
risk = risk[comxx]

##
xx = cbind(data$t.rfs, risk, data$e.rfs)
row.names(xx) = row.names(data)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(data)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
cv.res = cbind(predictability, f.mat, count, accuracy)
cv.res = cv.res[order(cv.res[,"Y"], decreasing=T), ]
plot(cv.res[,"Y"], cv.res[, "accuracy"])
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs")
cor(cv.res[,"Y"], cv.res[, "accuracy"], use="pairwise.complete.obs", method="s")
## positive correlation 

load(file= myoutf1)
myList[[7]] = cv.res
names(myList)[7] = "Predicted.Discovery.Result"
save(myList, file=myoutf1)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[5] Figure 6
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.1] Curtis model --> Rehman data  oncotype DX
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/WorSpa/w4_surv/result/Breast/Breast_UrRehman_GSE47561_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Ur-Rehman/UrRehman_intrinsic_subtype_PAM50_result.txt"
myinf3 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Ur-Rehman/Info/Breast_UrRehman_all_info.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/data/Curtis_predictability_OncotypeDX_models.rda"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig6_oncotypeDX_Curtis2Rehman.rda"

myList = list(NULL)

data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
pam50 = read.table(myinf2, header=T, sep="\t", quote="", row.names=1)
info <- read.table(myinf3, header=T, sep="\t", quote="", row.names=2)
comxx = intersect(row.names(data), row.names(info))
comxx = intersect(row.names(pam50), comxx)
E2F4 = data[comxx,8]
Pam50 = pam50[comxx,1]
info = info[comxx,]
data = cbind(E2F4, Pam50, info)
colnames(data)[1:8] = c("E2F4", "Pam50", "Dataset", "age", "ER", "grade", "size", "LN")
xx = as.numeric(as.character(data$age))
data$age = xx
xx = as.character(data$ER)
xx[xx=="Positive"] = 1
xx[xx=="Negative"] = 0
xx = as.integer(xx)
data$ER = xx
xx = as.character(data$grade)
xx = as.integer(xx)
data$grade = xx
xx = as.numeric(as.character(data$size))
data$size = xx
xx = as.character(data$LN)
xx[xx=="Positive"] = 1
xx[xx=="Negative"] = 0
xx = as.integer(xx)
data$LN = xx

#---------------------------------------------
load(file = myinf4)
library(randomForest)

score1 = predict(discovery.ppm.model, newdata= data, type="prob")[, "Y"]
score2 = predict(validation.ppm.model, newdata= data, type="prob")[, "Y"]
attr(validation.ppm.model$terms,"dataClasses")

data = cbind(score1, score2, data)

#---------------------------------------------
xx = data[order(data$score1, decreasing=T), ]
t.surv = as.numeric(as.character(xx$RFS_time))
e.surv = as.integer(as.character(xx$relapse_free_survival))
xx = cbind(t.surv, e.surv, xx)
xx = xx[!is.na(xx$t.surv), ]
xx = xx[!is.na(xx$e.surv), ]
myxx = xx

##
xx = cbind(myxx$t.surv, myxx$E2F4, myxx$e.surv)
row.names(xx) = row.names(myxx)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(myxx)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

myList[[1]] = res
names(myList)[[1]] = "discovery.model.result"


#----------------------------------------
xx = res[res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")


#---------------------------------------------
## based on score2 
#---------------------------------------------
xx = data[order(data$score2, decreasing=T), ]
t.surv = as.numeric(as.character(xx$RFS_time))
e.surv = as.integer(as.character(xx$relapse_free_survival))
xx = cbind(t.surv, e.surv, xx)
xx = xx[!is.na(xx$t.surv), ]
xx = xx[!is.na(xx$e.surv), ]
myxx = xx

##
xx = cbind(myxx$t.surv, myxx$E2F4, myxx$e.surv)
row.names(xx) = row.names(myxx)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(myxx)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

myList[[2]] = res
names(myList)[[2]] = "validation.model.result"
save(myList, file= myoutf1)

#----------------------------------------
xx = res[res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.2] Curtis model --> Rehman data  (gene70)
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/WorSpa/w4_surv/result/Breast/Breast_UrRehman_GSE47561_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Ur-Rehman/UrRehman_intrinsic_subtype_PAM50_result.txt"
myinf3 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Ur-Rehman/Info/Breast_UrRehman_all_info.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/data/Curtis_predictability_gene70_models.rda"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig6_gene70_Curtis2Rehman.rda"

myList = list(NULL)

data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
pam50 = read.table(myinf2, header=T, sep="\t", quote="", row.names=1)
info <- read.table(myinf3, header=T, sep="\t", quote="", row.names=2)
comxx = intersect(row.names(data), row.names(info))
comxx = intersect(row.names(pam50), comxx)
E2F4 = data[comxx,8]
Pam50 = pam50[comxx,1]
info = info[comxx,]
data = cbind(E2F4, Pam50, info)
colnames(data)[1:8] = c("E2F4", "Pam50", "Dataset", "age", "ER", "grade", "size", "LN")
xx = as.numeric(as.character(data$age))
data$age = xx
xx = as.character(data$ER)
xx[xx=="Positive"] = 1
xx[xx=="Negative"] = 0
xx = as.integer(xx)
data$ER = xx
xx = as.character(data$grade)
xx = as.integer(xx)
data$grade = xx
xx = as.numeric(as.character(data$size))
data$size = xx
xx = as.character(data$LN)
xx[xx=="Positive"] = 1
xx[xx=="Negative"] = 0
xx = as.integer(xx)
data$LN = xx

#---------------------------------------------
load(file = myinf4)
library(randomForest)

score1 = predict(discovery.ppm.model, newdata= data, type="prob")[, "Y"]
score2 = predict(validation.ppm.model, newdata= data, type="prob")[, "Y"]
attr(validation.ppm.model$terms,"dataClasses")

data = cbind(score1, score2, data)

#---------------------------------------------
xx = data[order(data$score1, decreasing=T), ]
t.surv = as.numeric(as.character(xx$RFS_time))
e.surv = as.integer(as.character(xx$relapse_free_survival))
xx = cbind(t.surv, e.surv, xx)
xx = xx[!is.na(xx$t.surv), ]
xx = xx[!is.na(xx$e.surv), ]
myxx = xx

##
xx = cbind(myxx$t.surv, myxx$E2F4, myxx$e.surv)
row.names(xx) = row.names(myxx)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(myxx)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

myList[[1]] = res
names(myList)[[1]] = "discovery.model.result"


#----------------------------------------
xx = res[res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")


#---------------------------------------------
## based on score2 
#---------------------------------------------
xx = data[order(data$score2, decreasing=T), ]
t.surv = as.numeric(as.character(xx$RFS_time))
e.surv = as.integer(as.character(xx$relapse_free_survival))
xx = cbind(t.surv, e.surv, xx)
xx = xx[!is.na(xx$t.surv), ]
xx = xx[!is.na(xx$e.surv), ]
myxx = xx

##
xx = cbind(myxx$t.surv, myxx$E2F4, myxx$e.surv)
row.names(xx) = row.names(myxx)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(myxx)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

myList[[2]] = res
names(myList)[[2]] = "validation.model.result"
save(myList, file= myoutf1)

#----------------------------------------
xx = res[res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[5.3] Curtis model --> Rehman data  (E2F4)
rm(list=ls())
myinf1 = "/lorax/chenglab/cc59/WorSpa/w4_surv/result/Breast/Breast_UrRehman_GSE47561_Symbol_iRAS.txt"
myinf2 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Ur-Rehman/UrRehman_intrinsic_subtype_PAM50_result.txt"
myinf3 = "/lorax/chenglab/cc59/PubDat/Cancer/BreastCancer/Ur-Rehman/Info/Breast_UrRehman_all_info.txt"
myinf4 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/data/Curtis_predictability_E2F4_models.rda"

myoutf1 = "/lorax/chenglab/cc59/WorSpa/w18_PPM/Figure/Fig6_E2F4_Curtis2Rehman.rda"

myList = list(NULL)

data = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
pam50 = read.table(myinf2, header=T, sep="\t", quote="", row.names=1)
info <- read.table(myinf3, header=T, sep="\t", quote="", row.names=2)
comxx = intersect(row.names(data), row.names(info))
comxx = intersect(row.names(pam50), comxx)
E2F4 = data[comxx,8]
Pam50 = pam50[comxx,1]
info = info[comxx,]
data = cbind(E2F4, Pam50, info)
colnames(data)[1:8] = c("E2F4", "Pam50", "Dataset", "age", "ER", "grade", "size", "LN")
xx = as.numeric(as.character(data$age))
data$age = xx
xx = as.character(data$ER)
xx[xx=="Positive"] = 1
xx[xx=="Negative"] = 0
xx = as.integer(xx)
data$ER = xx
xx = as.character(data$grade)
xx = as.integer(xx)
data$grade = xx
xx = as.numeric(as.character(data$size))
data$size = xx
xx = as.character(data$LN)
xx[xx=="Positive"] = 1
xx[xx=="Negative"] = 0
xx = as.integer(xx)
data$LN = xx

#---------------------------------------------
load(file = myinf4)
library(randomForest)

score1 = predict(discovery.ppm.model, newdata= data, type="prob")[, "Y"]
score2 = predict(validation.ppm.model, newdata= data, type="prob")[, "Y"]
attr(validation.ppm.model$terms,"dataClasses")

data = cbind(score1, score2, data)

#---------------------------------------------
xx = data[order(data$score1, decreasing=T), ]
t.surv = as.numeric(as.character(xx$RFS_time))
e.surv = as.integer(as.character(xx$relapse_free_survival))
xx = cbind(t.surv, e.surv, xx)
xx = xx[!is.na(xx$t.surv), ]
xx = xx[!is.na(xx$e.surv), ]
myxx = xx

##
xx = cbind(myxx$t.surv, myxx$E2F4, myxx$e.surv)
row.names(xx) = row.names(myxx)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(myxx)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

myList[[1]] = res
names(myList)[[1]] = "discovery.model.result"


#----------------------------------------
xx = res[res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")


#---------------------------------------------
## based on score2 
#---------------------------------------------
xx = data[order(data$score2, decreasing=T), ]
t.surv = as.numeric(as.character(xx$RFS_time))
e.surv = as.integer(as.character(xx$relapse_free_survival))
xx = cbind(t.surv, e.surv, xx)
xx = xx[!is.na(xx$t.surv), ]
xx = xx[!is.na(xx$e.surv), ]
myxx = xx

##
xx = cbind(myxx$t.surv, myxx$E2F4, myxx$e.surv)
row.names(xx) = row.names(myxx)
colnames(xx) = c("time", "score", "event")
xx[, "score"] = -xx[, "score"]

cc = 0
dd = 0
tt = 0
cnum = nrow(xx)
f.mat= matrix(0, cnum, 3)
colnames(f.mat) = c("concordance", "discordance", "tie")
row.names(f.mat) = row.names(myxx)
f.mat= as.data.frame(f.mat)
for(i in 1:(cnum-1))
  for(j in (i+1):cnum)
  {
  	if(xx[i,"event"]==0 & xx[j,"event"]==0) next
  	x1 = xx[i,"time"]-xx[j,"time"]
  	x2 = xx[i,"score"]-xx[j,"score"]
  	
  	if(xx[i,"event"]==1 & xx[j,"event"]==1)
  	{
  		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
  		if(x1*x2>0)  
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
  		if(x1*x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
  		next
  	}
   	if(xx[i,"event"]==1 & xx[j,"event"]==0)
 	{
 		if(x1>0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2<0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2>0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 		next
 	}
   	if(xx[i,"event"]==0 & xx[j,"event"]==1)
 	{
 		if(x1<0)  next
 		if(x2==0) 
  		{	
  			tt=tt+1
  			f.mat[i,"tie"] = f.mat[i,"tie"] +1
  			f.mat[j,"tie"] = f.mat[j,"tie"] +1
  		}
 		if(x2>0) 
  		{
  			cc = cc+1
   			f.mat[i,"concordance"] = f.mat[i,"concordance"] +1
  			f.mat[j,"concordance"] = f.mat[j,"concordance"] +1 			
  		}
 		if(x2<0) 
  		{
  			dd = dd+1
    		f.mat[i,"discordance"] = f.mat[i,"discordance"] +1
  			f.mat[j,"discordance"] = f.mat[j,"discordance"] +1 						
  		}
 	}
  }
  
#----------------------------------------------
cc
dd

res = f.mat
count = apply(f.mat, 1, sum)
accuracy = (f.mat[, "concordance"] + 0.5 * f.mat[, "tie"])/count
res = cbind(f.mat, count, accuracy)

myList[[2]] = res
names(myList)[[2]] = "validation.model.result"
save(myList, file= myoutf1)

#----------------------------------------
xx = res[res[, "count"]>nrow(data)*0.2,]
tmp = xx[, c("concordance", "discordance", "tie", "count")]
for(k in 2:nrow(tmp))
{
	tmp[k,] = tmp[k,] + tmp[k-1,]
}

cum.con = (tmp[,"concordance"] + tmp[, "tie"]*0.5)/tmp[, "count"]
plot(cum.con, pch=20, type="b")

w.siz = 200
step = 10
cnum = ceiling((length(cum.con)-w.siz)/step)
yy = rep(0, cnum)
for(k in 1:cnum)
{
	sta = (k-1)*step+1
	end = min(sta+w.siz-1, length(cum.con))
	yy[k] = mean(cum.con[sta:end])
}
plot(yy, pch=20, type="b")




















