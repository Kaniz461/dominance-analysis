# dominance-analysis

### R code for poisson regression model##

rm(list=ls())
library(haven)
k<-read_sav('C:/the/bdhs 2017-18/main1.sav')
attach(k)
names(k)
library(MASS)
library(corpcor)
varlist<-c("ANC", "DIVISION_NEW","Wealth_index","place","Birth_order",
"media_exposure","education","problem","wanted","age")
library(yhat)
k$ANC=as.numeric(k$ANC)
k$DIVISION_NEW=as.factor(k$DIVISION_NEW)
k$Wealth_index=as.factor(k$Wealth_index)
k$Birth_order=as.factor(k$Birth_order)
k$place=as.factor(k$place)
k$media_exposure=as.factor(k$media_exposure)
k$education=as.factor(k$education)
k$age=as.factor(k$age)
k$problem=as.factor(k$problem)
k$wanted=as.factor(k$wanted)
lm.out<-glm(ANC~DIVISION_NEW+Wealth_index+place+Birth_order+media_exposure
+education+problem+wanted+age,data=k,family=poisson())
a<-summary(lm.out)
a
b<-a$coefficients
c<-b[,1]
c
exp(c)
library(modEvA)
Dsquared(lm.out)
dv<-"ANC"
ivlist<-c("DIVISION_NEW","Wealth_index","place","Birth_order","media_exposure",
"education","problem","wanted","age")
ilist <- unlist(ivlist)
cols <- length(ilist)
dv<-as.factor(dv)
prmtn <- function(n) { 
	##function permutations from package e1071 
	if (n == 1) 
	return(matrix(1)) 
	else if (n < 2) 
	stop("n must be a positive integer") 
	z <- matrix(1) 
	for (i in 2:n) { 
		x <- cbind(z, i) 
		a <- c(1:i, 1:(i - 1)) 
		z <- matrix(0, ncol = ncol(x), nrow = i * nrow(x)) 
		z[1:nrow(x), ] <- x 
		for (j in 2:i - 1) { 
			z[j * nrow(x) + 1:nrow(x), ] <- x[, a[1:i + j]] 
		} 
	} 
	dimnames(z) <- NULL 
	z 
}

#prmtn<-prmtn(cols)[i,]
aps_formula <- function(dv, ilist, prm)
{
	fd <- paste(dv, "~", sep="")
	p <- sort(prm)
	fi <- paste(ilist[p], collapse="+")
	formula <- paste(fd, fi, sep="")
}




lmgm <- function(dataMatrix, dv, ivlist)
{
	
	dom <- matrix(, nrow = dim(prmtn(cols))[1], ncol=cols)
	order_v <- vector(mode = "character", length = dim(prmtn(cols))[1])
	r_sq <- vector(mode="numeric")
	tmp_r <- vector(mode = "numeric", length = cols)
	b <- vector(mode = "numeric", length = cols)
	
	
	count <- 0
	for (i in 1:dim(prmtn(cols))[1])
	{
		prmtn<-prmtn(cols)[i,]
		order_v[i] <-  paste(ivlist[prmtn], collapse="")
		for(j in 1:cols)
		{
			
			f <- aps_formula(dv, ivlist, prmtn[1:j])
			
			if(j == 1) a <- 0
			else a <- b[j-1]
			if(is.na(r_sq[f]))
			{  
				r_sq[f]<-Dsquared(glm(f,dataMatrix,family=poisson()))
				count <- count+1	
			}
			b[j] <- r_sq[f]
			tmp_r[j] = b[j] - a
		}
		dom[i, ] <- tmp_r[order(prmtn)]
	}
	rownames(dom) <- order_v
	CR=apply(dom,2,mean)
	return(list(DA = dom, CR = CR, C = count))
}


g<-lmgm(k,dv,ivlist)
g
A<-sort(g$CR)
apply(g$DA,1,sum)
poisson<-(g$CR/sum(g$CR))*100

## R-code for finding dispersion statistics value##


#From model fitting
rss<-sum(residuals(lm.out,type="pearson")^2)
disp=rss/lm.out$df.residual
round(disp,3)

#Using formula
mu<-predict.glm(lm.out,type="response")
pearson<-sum((ANC-mu)^2/mu)
dispersion<-pearson/lm.out$df.residual
round(dispersion,3)

#R-code for calculating Z-score test

case(1): Dean and lawless (1993)
z_dean=((ANC-mu)^2-ANC)/(sqrt(2)*mu)
z_score=lm(z_dean~1)
summary(z_score)


case(2):Winkelmann (2008)
z_win=((ANC-mu)^2-ANC)/(2*mu)
z_win=lm(z_win~1)
summary(z_win)


case(3):Cameron (2008)
z_cam=((ANC-mu)^2-ANC)/mu
z_cam=lm(z_cam~1)
summary(z_cam)


##R-code for calculating Lagrange multiplier test 
n=nrow(k)
lag_multi=((sum(mu^2)-n*mean(mu))^2)/(2*sum(mu^2))
lag_multi
pchisq(lag_multi,1,lower.tail=F)


##R-code for Negative binomial model
rm(list=ls())
library(haven)
k<-read_sav('C:/thesis/bdhs 2017-18/main1.sav')
k
attach(k)
names(k)

library(MASS)
library(corpcor)

varlist<-c("ANC", "DIVISION_NEW","Wealth_index","Birth_order","media_exposure",
"education","age","problem","wanted","place")

library(yhat)
k$ANC=as.numeric(k$ANC)
k$DIVISION_NEW=as.factor(k$DIVISION_NEW)
k$Wealth_index=as.factor(k$Wealth_index)

k$Birth_order=as.factor(k$Birth_order)
k$place=as.factor(k$place)
k$media_exposure=as.factor(k$media_exposure)
k$education=as.factor(k$education)
k$age=as.factor(k$age)
k$problem=as.factor(k$problem)
k$wanted=as.factor(k$wanted)


lm.out<-glm.nb(ANC~DIVISION_NEW+Wealth_index+place+Birth_order+media_exposure
+education+problem+wanted+age,data=k)
a<-summary(lm.out)
a
b<-a$coefficients
c<-b[,1]
c
exp(c)
library(modEvA)
Dsquared(lm.out)


dv<-"ANC"
ivlist<-c("DIVISION_NEW","Wealth_index","place","Birth_order","media_exposure",
"education","age","problem","wanted")
ilist <- unlist(ivlist)
cols <- length(ilist)
dv<-as.factor(dv)
prmtn <- function(n) { 
	##function permutations from package e1071 
	if (n == 1) 
	return(matrix(1)) 
	else if (n < 2) 
	stop("n must be a positive integer") 
	z <- matrix(1) 
	for (i in 2:n) { 
		x <- cbind(z, i) 
		a <- c(1:i, 1:(i - 1)) 
		z <- matrix(0, ncol = ncol(x), nrow = i * nrow(x)) 
		z[1:nrow(x), ] <- x 
		for (j in 2:i - 1) { 
			z[j * nrow(x) + 1:nrow(x), ] <- x[, a[1:i + j]] 
		} 
	} 
	dimnames(z) <- NULL 
	z 
}

#prmtn<-prmtn(cols)[i,]
aps_formula <- function(dv, ilist, prm)
{
	fd <- paste(dv, "~", sep="")
	p <- sort(prm)
	fi <- paste(ilist[p], collapse="+")
	formula <- paste(fd, fi, sep="")
}




lmgm <- function(dataMatrix, dv, ivlist)
{
	
	dom <- matrix(, nrow = dim(prmtn(cols))[1], ncol=cols)
	order_v <- vector(mode = "character", length = dim(prmtn(cols))[1])
	r_sq <- vector(mode="numeric")
	tmp_r <- vector(mode = "numeric", length = cols)
	b <- vector(mode = "numeric", length = cols)
	
	
	count <- 0
	for (i in 1:dim(prmtn(cols))[1])
	{
		prmtn<-prmtn(cols)[i,]
		order_v[i] <-  paste(ivlist[prmtn], collapse="")
		for(j in 1:cols)
		{
			
			f <- aps_formula(dv, ivlist, prmtn[1:j])
			
			if(j == 1) a <- 0
			else a <- b[j-1]
			if(is.na(r_sq[f]))
			{  
				r_sq[f]<-Dsquared(glm.nb(f,dataMatrix))
				count <- count+1	
			}
			b[j] <- r_sq[f]
			tmp_r[j] = b[j] - a
		}
		dom[i, ] <- tmp_r[order(prmtn)]
	}
	rownames(dom) <- order_v
	CR=apply(dom,2,mean)
	return(list(DA = dom, CR = CR, C = count))
}


n<-lmgm(k,dv,ivlist)
n
neg<-(n$CR/sum(n$CR))*100
neg




#From Dominance analysis
#Poisson distribution
percentage<-c(9.67,4.27,28.78,19.90,9.35,3.12,20.39,1.78,2.69)
variables<-c("Place","Region","educational","Wealth index","Birth order","Age",
"Media exposure","Wanted","Distance")
a<-data.frame(variables,percentage)
barplot(a$percentage,names=a$variables,col = "cyan3")

#Negative binomial distribution
percentages<-c(9.09,4.76,28.80,19.72,9.55,1.67,20.63,2.73,3.02)
b<-data.frame(variables,percentages)
barplot(b$percentages,names=b$variables,col = "coral2")
