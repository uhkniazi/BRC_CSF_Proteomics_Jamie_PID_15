# Name: scratch.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 30/01/2019
# Desc: Collection of odd tasks

### making the heatmap
if (!require(NMF)) stop('R package NMF needs to be installed.')
mData = read.csv('dataExternal/expressionJamie.csv', header=T, row.names=1)
str(mData)
cn = gsub('\\.\\d', '', colnames(mData))
colnames(mData) = cn

mData = as.matrix(mData)

table(is.na(mData))
mData[is.na(mData)] = runif(52, 0.01)
dim(na.omit(mData))
dim(mData)
mData = log(mData+0.5)

# standardize the variables
s = apply(mData, 1, sd)
## remove any variables with sd 0
f = s <= 0
table(f)
s = s[!f]
mData = mData[!f,]
mData = t(scale(t(mData)))
# cluster the samples
hc = hclust(dist(t(mData)))
# cluster the variables
hcV = hclust(dist(mData))
# sanity check
range(mData)
ivScale = c(-3, 3)
# threshhold the values
mData[mData < ivScale[1]] = ivScale[1]
mData[mData > ivScale[2]] = ivScale[2]
# draw the heatmap  color='-RdBu:50'
col=c('blue', 'black', 'red')
aheatmap(mData, color=col, breaks=0, scale='none', Rowv = hcV, annColors=NA, Colv=NA)


##############################################################################################

## utility function for plotting
getms = function(f){
  m = mean(f)
  se = sd(f)
  m.up = m+1.96*se
  m.down = m-1.96*se
  ret= c(m, m.up, m.down)
  names(ret) = c('m', 'm.up', 'm.down')
  return(ret)
}


mTreatment = mCoef[,aqp4]
colnames(mTreatment) = as.character(d$fBatch[aqp4])
df = apply(mTreatment, 2, getms)
x = 1:ncol(mTreatment)

#par(p.old)
plot(x, df['m',], ylim=c(min(df), max(df)), pch=20, xlab='', main='Model Estimated Average Expression for AQP4',
     ylab='log Average', xaxt='n')
axis(1, at = x, labels = colnames(mTreatment), las=2, cex.axis=0.7)
for(l in 1:ncol(df)){
  lines(x=c(x[l], x[l]), y=df[c(2,3),l], lwd=0.5)
}
abline(h = 0, col='grey')


