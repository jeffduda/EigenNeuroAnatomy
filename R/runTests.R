library(ANTsR)
library(ggplot2)
#library(randomForest)
#library(spls)
#library(BGLR)

#source("~/pkg/EigenNeuroAnatomy/R/explainedVariance.R")

filelist = "~/pkg/EigenNeuroAnatomy/data/grossman_info.csv"
info = read.csv(filelist)
filenames = as.character(info$file)
nSubjects = length(filenames)

# use pre-defined mask, or data-specific
#mask = antsImageRead( "/data/grossman/pipedream2018/templates/OASIS/priors/priors2.nii.gz")
meanImg = antsImageRead(filenames[1])
for ( i in 2:nSubjects ) {
  meanImg = meanImg + antsImageRead(filenames[i])
}
meanImg = meanImg / nSubjects
mask = meanImg*0
mask[ meanImg > 0.5 ] = 1

ref = resampleImage(mask, c(4,4,4))

tx = createAntsrTransform( dimension=3, type="Euler3DTransform" )

#mask = applyAntsrTransform( tx, mask, reference=ref, interpolation="Linear" )
#mask[mask<0.5] = 0
#mask[mask>0] = 1


# optional sigma value for smoothing - faster than smoothing in call to sparseDecom?
# could we "cheat" here: down-sample data, get regions, then upsample regions to full resolution
#  this would save time, help with smaller data sets possibly?
#mat = imageListToMatrix( filenames, mask )
mat = matrix( 0, nrow=nSubjects, ncol=sum(mask) )
for ( i in 1:nSubjects ) {
  print(filenames[i])
  img = antsImageRead(filenames[i])

  # If downsampling
  # img = applyAntsrTransformToImage( tx, img, ref, interpolation="linear")

  # Intra-subject, relative measure
  # img = img / mean( img[mask>0] ) # intra-subject normalization

  # Smoothing
  # img = smoothImage( img, sigma=2 )

  mat[i,] = img[mask>0]

}

train = sample(nSubjects)[1:80]
testIdx = rep(1,99)
testIdx[train] = 0
testIdx = which(testIdx==1)

# group-based relative measure or nuisance regressors
# mat = resid( lm( mat ~ rowMeans(mat)) )

# scale data only ( for alternate decomp methods )
# mat = scale(mat, center=FALSE)

# center & scale data
origMat = mat
mat = scale(mat)

# ICA-whiten
mat = icawhiten( mat, nSubjects )

nComp = 60
eanatStruct = sparseDecom( inmatrix=mat[train,],
 inmask = mask,
 sparseness=1.0/nComp,
 nvecs = nComp,
 cthresh=1000,
 its = 2,
 mycoption=1,
 maxBased = F,
 verbose=T )

eanat = eanatStruct$eigenanatomyimages

e = explainedVariance( mat[train,], eanat, by.component=T)
reor = sort(e, decreasing=T, index.return=T)
eanat = eanat[reor$ix,]

trainDat = data.frame(Age=info$AgeatMRI[train], Predictor=mat[train,] %*% t(eanat) )
form = "Age ~ Predictor.1"
for ( i in 2:nComp ) {
  form = paste0(form, " + ", "Predictor.", i)
}
form = as.formula(form)
ageMod = lm( form, trainDat )
agePredTrain = predict(ageMod)
trainDat$PredictedAge = agePredTrain
trainDat$Group = rep("Train", dim(trainDat)[1])

testDat = data.frame(Age=info$AgeatMRI[testIdx], Predictor = mat[testIdx,] %*% t(eanat))
agePredTest = predict( ageMod, newdata=testDat )
testDat$PredictedAge = agePredTest
testDat$Group = rep("Test", dim(testDat)[1])

g.dat = rbind(trainDat, testDat)

m1 = min( c(g.dat$Age, g.dat$PredictedAge))
m2 = max( c(g.dat$Age, g.dat$PredictedAge))
g = ggplot(g.dat, aes(x=Age,y=PredictedAge,group=Group,colour=Group)) + geom_point() + xlim(m1,m2) + ylim(m1,m2)
print(g)
