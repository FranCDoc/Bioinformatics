library(gplots)
library(ggplot2)
library(mRMRe)
library(LiblineaR)

load("/Users/marcelochiesa/Desktop/EXAMS\ LLN/LGBIO2010/Glioblastoma.Rdata")

#-----------------------------POINT-2-----------------------------#

# Reading data
dmatrix<- data.matrix(data)
columns<-(ncol(dmatrix))-1
rows<-nrow(dmatrix)

# Mean and variance
v<-rep(0L,columns)
m<-rep(0L,columns)
for (i in 1:columns){
  v[i]<-var(dmatrix[1:rows,i])
  m[i]<-mean(dmatrix[1:rows,i])
}

# Rank variances
vrank<-order(v) #index
v.ranked<-v[vrank]
plot(1:54613,v.ranked[1:54613], xlab= "probesets", ylab = "variance value")

# Rank normalized matrix by variance
data.ranked<-dmatrix[,vrank]

# Taking the 25% with higher variance
limit<-round(columns*0.75)
print(columns-limit) # Q2.2
print(v.ranked[limit]) # Q2.2

# Top five by variance
five<-columns-4
top5<-data.ranked[1:rows,five:columns] 
print(top5) # Q2.3
print(v.ranked[five:columns]) # Q2.3

# Keep the higher 25% by variance
top25<-data.ranked[1:rows,limit:columns]

#-----------------------------POINT-3-----------------------------#

# Feature normalization
top25.scale<-scale(top25,center=TRUE,scale=TRUE)
top25.matrix<-cbind(top25.scale,data[1:rows,columns+1])
top25.columns<-(ncol(top25.matrix))-1

v.norm<-rep(0L,top25.columns)
m.norm<-rep(0L,top25.columns)

for (i in 1:top25.columns){
 v.norm[i]<-var(top25.scale[1:rows,i])
 m.norm[i]<-mean(top25.scale[1:rows,i])
}

label1<-rep(0L,rows)
label2<-rep(0L,rows)
p.value<-rep(0L,top25.columns)

# Find distribution of genes among classes
for (j in 1:top25.columns){
  cont_label1<-1
  cont_label2<-1
  for (i in 1:rows){
    if (top25.matrix[i,top25.columns+1]==1){
      label1[cont_label1]<-top25.matrix[i,j]
      cont_label1<-cont_label1+1
    }
    if (top25.matrix[i,top25.columns+1]==2){
      label2[cont_label2]<-top25.matrix[i,j]
      cont_label2<-cont_label2+1
    }
  }
  # Find p-value in t-test
  test<-t.test(label1[1:cont_label1-1],label2[1:cont_label2-1])
  p.value[j]<-test[["p.value"]]
}

# Rank p-values without correction
p.value.rank<- order(p.value) #index
p.value.ranked<-p.value[p.value.rank]

p.value.5PC<- sum(p.value.ranked < 0.05) #Q3.2
print(p.value.5PC)
plot(1:top25.columns,p.value.ranked)

# Rank top 25% according to p-value
top25.matrix.ranked.p<-top25.matrix[,p.value.rank]

# Bonferroni correction (BC)
BC.threshold<- 0.05/top25.columns # New threshold

# Rank p-values with BC correction
# First method
p.value.BC.lower<- sum(p.value.ranked < BC.threshold) #Q3.2
print(p.value.BC.lower)
# Second method
BC<-p.adjust(p.value.ranked, method = "bonferroni")
BC.lower<- sum(BC < 0.05) #Q3.2
print(BC.lower)

# Top ten with BC correction
top10.BC<-rbind(top25.matrix.ranked.p[0,1:10],BC[1:10])
print(top10.BC) #Q3.3

# False Discovery Rate correction (FDR)
# First method
for (i in 1:top25.columns){
  if ((p.value.ranked[i]*top25.columns/i) < 0.05){
    index<- i
  }
}
print(index) #Q3.4
# Second method
FDR<-p.adjust(p.value.ranked, method = "fdr")
FDR.lower<- sum(FDR < 0.05) #Q3.4
print(FDR.lower)

# Top ten FDR
top10.FDR<-rbind(top25.matrix.ranked.p[0,1:10],FDR[1:10])
print(top10.FDR) #Q3.4

#-----------------------------POINT-4-----------------------------#
# Heatmap
rc<-rainbow(nrow(top25.matrix.ranked.p[1:rows,1:50]),start=0,end=0.2)
cc<-rainbow(ncol(top25.matrix.ranked.p[1:rows,1:50]),start=0,end=0.2)
heatmap.2(top25.matrix.ranked.p[1:rows,1:50],ColSideColors=cc,RowSideColors=rc) #Q4.1,Q4.2

# Colors
label.colors<-rep(0L,rows)
for (i in 1:rows){
  if (top25.matrix[i,top25.columns+1]==1){
    label.colors[i] <- "red"
  }
  else {
    label.colors[i] <- "blue"
  }
}

# 2D plot of 2 most significant genes
print(top25.matrix.ranked.p[0,1:2]) #names
plot(top25.matrix.ranked.p[1:rows,1],top25.matrix.ranked.p[1:rows,2],xlab="X_22349",ylab="X_27963",col=label.colors) #Q4.3

# 2D plot of 2 least significant genes
print(top25.matrix.ranked.p[0,7953:7954]) #names
plot(top25.matrix.ranked.p[1:rows,7953],top25.matrix.ranked.p[1:rows,7954],xlab="X_53516",ylab="X_48207",col=label.colors) #Q4.3

#-----------------------------POINT-5-----------------------------# 

# Train and test sets
train <-top25.matrix[1:82,1:13654]
test.svm <-top25.matrix[83:103,1:13654]
class <- top25.matrix[1:82,13655] 
vali.class<-top25.matrix[83:103,13655] 

# Support Vector Machine

# We used the optimization configuration suggested in this link: 
# https://www.rdocumentation.org/packages/LiblineaR/versions/2.10-8/topics/LiblineaR

# Find the best model with the best cost parameter via 10-fold cross-validations
tryTypes=c(0:7)
tryCosts=c(1000,1,0.001)
bestCost=NA
bestAcc=0
bestType=NA

for(ty in tryTypes){
  for(co in tryCosts){
    acc=LiblineaR(data=train,target=class,type=ty,cost=co,bias=1,cross=10,verbose=FALSE)
    cat("Results for C=",co," : ",acc," accuracy.\n",sep="")
    if(acc>bestAcc){
      bestCost=co
      bestAcc=acc
      bestType=ty
    }
  }
}

cat("Best model type is:",bestType,"\n")
cat("Best cost is:",bestCost,"\n")
cat("Best accuracy is:",bestAcc,"\n")

# Re-train best model with best cost value.
SVM=LiblineaR(data=train,target=class,type=bestType,cost=bestCost,bias=1,verbose=FALSE)


# Rank by weight
w<-abs(SVM$W[1,1:13654])
w.rank<-order(w) #index
w.ranked<-w[w.rank]

# Top ten by weight
plot(1:13654,w.ranked[1:13654])
print(w.ranked[(13654-9):13654])

# Respective p-value ranking
weight.names<-names(w.ranked[(13654-9):13654])
print(weight.names)
p.value.names<-colnames(top25.matrix.ranked.p)

position<-rep(0L,10)
for (j in 1:10){
  for(i in 1:top25.columns){
    if (weight.names[j]==p.value.names[i]){
      position[j]<-i #index in p-value rank
    }
  }
}

print(w.ranked[(13654-9):13654])
print(position)


# Prediction based on SVM
prediction<-predict(SVM, test.svm) # Validation set
matrix<-table(prediction$predictions,vali.class)

# Confusion matrix on SVM
# TN|FN
# FP|TP
print(matrix) 

# Classification accuracy on SVM
BCR = mean(c(matrix[1,1]/sum(matrix[,1]),matrix[2,2]/sum(matrix[,2])))
print(BCR)

#-----------------------------POINT-6-----------------------------#

feature_data<- mRMR.data(data = data.frame(top25.matrix))

# Mutual information of genes with the label
# V13655 is the name of labels
mutual_info<- mim(feature_data)
MI.1<- mutual_info[,13655]

MI.order.index<-order(MI.1) #index
print(MI.order.index[13645:13655])
MI.ranked<-MI.1[MI.order.index]
print(MI.ranked[13645:13655])

# Select genes according to the mRMR algorithm
data.mrmr<- mRMR.data(data = data.frame(top25.matrix))

fs<- mRMR.classic(data.mrmr, target_indices=13655, feature_count=10)

top25.matrix[1,t(fs@filters[["13655"]])]
print(fs@scores[["13655"]])
