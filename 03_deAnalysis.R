# File: 03_deAnalysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Modelling and selecting DE genes
# Date: 07/09/2018

## load the data
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 30) AND (MetaFile.comment like "%count%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n)

## load the metadata i.e. covariates
q = paste0('select Sample.* from Sample where Sample.idData = 30')
dfSample = dbGetQuery(db, q)
dim(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)


## make count matrix
colnames(dfCounts)
mCounts = as.matrix(dfCounts[,-c(1:8)])
rownames(mCounts) = dfCounts$Accession
mData = mCounts
dim(mData)

## remove NAs by adding a jitter
table(is.na(mData))
mData[is.na(mData)] = runif(5384, 0.01)
dim(na.omit(mData))
dim(mData)

summary(mData)
# drop the samples where average across rows is less than 2
i = rowMeans(mData)
table( i < 2)
mData = mData[!(i< 2),]
dim(mData)

## perform DE analysis
## delete sample section after testing
hist(mData)
plot(density(mData))
mData.norm = mData
# set.seed(123)
# i = sample(1:nrow(mData.norm), 10, replace = F)
# dfData = data.frame(t(log(mData.norm[i,]+0.5)))

dfData = data.frame(t(log(mData.norm+0.5)))
dfData = stack(dfData)
dim(dfData)
f = factor(dfSample$group1)
levels(f)
dfData$fBatch = factor(dfSample$group1)
dfData$fAdjust1 = factor(dfSample$group2)
dfData$fAdjust2 = factor(dfSample$group3)
dfData$Coef = factor(dfData$fBatch:dfData$ind)
dfData$Coef.adj1 = factor(dfData$fAdjust1:dfData$ind)
dfData$Coef.adj2 = factor(dfData$fAdjust2:dfData$ind)
dim(dfData)
dfData = droplevels.data.frame(dfData)
dfData = dfData[order(dfData$Coef, dfData$Coef.adj1, dfData$Coef.adj2), ]
str(dfData)

# # ## setup the model
library(lme4)
fit.lme1 = lmer(values ~ 1  + (1 | Coef) + (1 | Coef.adj1) + (1 | Coef.adj2), data=dfData)
summary(fit.lme1)
fit.lme2 = lmer(values ~ 1 + (1 | Coef) + (1 | Coef.adj1), data=dfData)
summary(fit.lme2)
fit.lme3 = lmer(values ~ 1 + (1 | Coef) + (1 | Coef.adj2), data=dfData)
summary(fit.lme3)
fit.lme4 = lmer(values ~ 1 + (1 | Coef), data=dfData)
summary(fit.lme4)

anova(fit.lme1, fit.lme2, fit.lme3, fit.lme4)

ran = ranef(fit.lme1, condVar=F)

plot((fitted(fit.lme1)), resid(fit.lme1), pch=20, cex=0.7)
lines(lowess((fitted(fit.lme1)), resid(fit.lme1)), col=2)
plot(density(dfData$values))
lines(density(fitted(fit.lme1)))
lines(density(fitted(fit.lme2)))
lines(density(fitted(fit.lme3)))
lines(density(fitted(fit.lme4)))

plot(resid(fit.lme1), dfData$fBatch)
plot(resid(fit.lme1), dfData$fAdjust1)
plot(resid(fit.lme1), dfData$fAdjust2)
plot(resid(fit.lme1), dfData$ind)

## setup the stan model
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='tResponse3RandomEffectsMultipleScales.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$values), 2*sd(dfData$values))
# # ## set initial values
# ran = ranef(fit.lme1)
# r1 = ran$Coef
# r2 = ran$Coef.adj1
# r3 = ran$Coef.adj2
# 
# r1 = rep(0, times=nlevels(dfData$Coef))
# r2 = rep(0, times=nlevels(dfData$Coef.adj1))
# #r3 = rep(0, times=nlevels(dfData$ind))
# initf = function(chain_id = 1) {
#   list(sigmaRan1 = 0.1, sigmaRan2=2, #rGroupsJitter1=r1, rGroupsJitter2=r2,
#        iSize=8)
# }

### try a t model without mixture
## subset the data to get the second level of nested parameters
## this is done to avoid loops in the stan script to map the scale parameters
## of each ind/gene to the respective set of coefficients for jitters
d = dfData[!duplicated(dfData$Coef), ]

lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj1),
                 Nclusters3=nlevels(dfData$Coef.adj2),
                 NScaleBatches1 = nlevels(dfData$ind), # to add a separate scale term for each gene
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj1),
                 NgroupMap3=as.numeric(dfData$Coef.adj2),
                 NBatchMap1=as.numeric(d$ind), # this is where we use the second level mapping
                 y=dfData$values, 
                 gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=3,
                    pars=c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaRan3',
                           'nu', 'sigmaPop', #'mu',
                           'rGroupsJitter1', 'rGroupsJitter2', 'rGroupsJitter3'),
                    cores=3, control=list(adapt_delta=0.99, max_treedepth = 12))#, init=initf, 
#save(fit.stan, file='temp/fit.stan.mul.scale.10.sep.rds')

print(fit.stan, c('betas', 'sigmaRan2', 'sigmaRan3', 'sigmaPop', 'nu'), digits=3)
traceplot(fit.stan, 'betas')
traceplot(fit.stan, c('sigmaRan2', 'sigmaRan3'))

# m = extract(fit.stan)$mu
# plot(density(dfData$values))
# lines(density(colMeans(m)))
# dim(m)
# apply(m[sample(1:nrow(m), 100),], 1, function(x) lines(density(x)))
## get the coefficient of interest - Modules in our case from the random coefficients section
mCoef = extract(fit.stan)$rGroupsJitter1
dim(mCoef)
#dim(extract(fit.stan)$sigmaRan1)
## get the intercept at population level
iIntercept = as.numeric(extract(fit.stan)$betas)
# add the intercept to each random effect variable, to get the full coefficient
mCoef = sweep(mCoef, 1, iIntercept, '+')

## function to calculate statistics for differences between coefficients
getDifference = function(ivData, ivBaseline){
  stopifnot(length(ivData) == length(ivBaseline))
  # get the difference vector
  d = ivData - ivBaseline
  # get the z value
  z = mean(d)/sd(d)
  # get 2 sided p-value
  p = pnorm(-abs(mean(d)/sd(d)))*2
  return(list(z=z, p=p))
}

## split the data into the comparisons required
d = data.frame(cols=1:ncol(mCoef), mods=levels(dfData$Coef))
## split this factor into sub factors
f = strsplit(as.character(d$mods), ':')
d = cbind(d, do.call(rbind, f))
head(d)
colnames(d) = c(colnames(d)[1:2], c('fBatch', 'ind'))
d$split = factor(d$ind)

levels(d$fBatch)
## repeat this for each comparison

## get a p-value for each comparison
l = tapply(d$cols, d$split, FUN = function(x, base='Pre-Control', deflection='Pre-Recovered') {
  c = x
  names(c) = as.character(d$fBatch[c])
  dif = getDifference(ivData = mCoef[,c[deflection]], ivBaseline = mCoef[,c[base]])
  r = data.frame(ind= as.character(d$ind[c[base]]), coef.base=mean(mCoef[,c[base]]), 
                 coef.deflection=mean(mCoef[,c[deflection]]), zscore=dif$z, pvalue=dif$p)
  r$difference = r$coef.deflection - r$coef.base
  #return(format(r, digi=3))
  return(r)
})

dfResults = do.call(rbind, l)
dfResults$adj.P.Val = p.adjust(dfResults$pvalue, method='BH')

### plot the results
dfResults$logFC = dfResults$difference
dfResults$P.Value = dfResults$pvalue
head(rownames(dfResults))
library(org.Hs.eg.db)

df = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(dfResults$ind), columns = 'SYMBOL', keytype = 'UNIPROT')
i = match(as.character(dfResults$ind), as.character(df$UNIPROT))
df = df[i,]
dfResults$SYMBOL = df$SYMBOL
identical(as.character(dfResults$ind), as.character(df$UNIPROT))
## produce the plots 
range(dfResults$logFC)
f_plotVolcano(dfResults, 'Pre-Recovered vs Pre-Control', fc.lim=range(dfResults$logFC))

m = tapply(dfData$values, dfData$ind, mean)
i = match(rownames(dfResults), names(m))
m = m[i]
identical(names(m), rownames(dfResults))
plotMeanFC((m), dfResults, 0.1, '')
table(dfResults$pvalue < 0.05)
table(dfResults$adj.P.Val < 0.05)
## save the results 
write.csv(dfResults, file='results/DEAnalysisPre-RecoveredVsPre-Control.xls')

s = extract(fit.stan)$sigmaRan1
dim(s)
hist(colMeans(s))
s = colMeans(s)
names(s) = levels(dfData$ind)
