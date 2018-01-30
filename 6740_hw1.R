myData = read.table(file = "semeion.csv", sep = ",")
myX = data.matrix(myData[,1:256])

####kmeans 
k=10
km <- kmeans(myX, k, nstart = 20)
#### initial estimates of gamaIK
gamaIK = matrix(0,1593,k)
for(i in 1:1593){
  clusI = km$cluster[i]
  gamaIK[i,clusI] = 1
}
#### initial estimates of miuK
NK <- colSums(gamaIK)
sumGamaX = matrix(0,10,256)
miuK = matrix(0,10,256)
for(k in 1:10){
  for(i in 1:1593){
    sumGamaX[k,] = sumGamaX[k,] + gamaIK[i,k] * myX[i,]
  }
  miuK[k,] = 1/NK[k] * sumGamaX[k,]
}
paiK = 1/1593 * NK
#### initial estimates of sigma
MNK <- list(matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
            ,matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
            ,matrix(0,256,256),matrix(0,256,256))
sumMIUX <- list(matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
                ,matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
                ,matrix(0,256,256),matrix(0,256,256))
for(k in 1:10){
  for(i in 1:1593){
    sumMIUX[[k]] = sumMIUX[[k]] + gamaIK[i,k] * as.matrix((myX[i,]-miuK[k,])) %*% t(as.matrix((myX[i,]-miuK[k,])))
  }
  MNK[[k]] = 1/NK[k] * sumMIUX[[k]]
}
d=256
q=6
qplus = q+1
VQK <- list()
for(k in 1:10){VQK[[k]] = eigen(MNK[[k]])$vectors[,1:q]}
squarE <- list()
for(k in 1:10){
  squarE[[k]] = 1/(d-q)*sum(eigen(MNK[[k]])$values[qplus:d])
}
WQK = list()
for(k in 1:10){
  WQK[[k]] = matrix(0,q,q)
  for( i in 1:q){
    WQK[[k]][i,i] = sqrt(eigen(MNK[[k]])$values[i]-squarE[[k]])
  }
  WQK[[k]] = VQK[[k]] %*% WQK[[k]]
}
Ksigma <- list()
for(k in 1:10){
  Ksigma[[k]] = WQK[[k]] %*% t(WQK[[k]]) + squarE[[k]] * diag(d)
}
###########################
loglik = matrix(0,25,1)
library(mvtnorm)
####iteration begins
for (nstep in 1:25){
  #### E step
  prob <- matrix(0,1593,10)
  for(k in 1:10){
    prob[,k] = dmvnorm(myX, mean = miuK[k,], sigma = Ksigma[[k]], log = FALSE)
  }
  sumPaiProb = matrix(0,1593,1)
  for(k in 1:10){
    sumPaiProb[,1] = sumPaiProb[,1] + paiK[k] * prob[,k]
  }
  gamaIK = matrix(0,1593,10)
  for(k in 1:10){
    for(i in 1:1593){
      gamaIK[i,k] = paiK[k] * prob[i,k]/sumPaiProb[i,]
    }
  }

  logprob = 0
  for(i in 1:1593){
    logprob = logprob + log(sumPaiProb[i,])
  }
  loglik[nstep,1] = logprob
  #### M step
  NK <- colSums(gamaIK)
  sumGamaX = matrix(0,10,256)
  miuK = matrix(0,10,256)
  for(k in 1:10){
    for(i in 1:1593){
      sumGamaX[k,] = sumGamaX[k,] + gamaIK[i,k] * myX[i,]
    }
    miuK[k,] = 1/NK[k] * sumGamaX[k,]
  }
  paiK = 1/1593 * NK

  MNK <- list(matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
              ,matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
              ,matrix(0,256,256),matrix(0,256,256))
  sumMIUX <- list(matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
                  ,matrix(0,256,256),matrix(0,256,256),matrix(0,256,256),matrix(0,256,256)
                  ,matrix(0,256,256),matrix(0,256,256))
  for(k in 1:10){
    for(i in 1:1593){
      sumMIUX[[k]] = sumMIUX[[k]] + gamaIK[i,k] * as.matrix((myX[i,]-miuK[k,])) %*% t(as.matrix((myX[i,]-miuK[k,])))
    }
    MNK[[k]] = 1/NK[k] * sumMIUX[[k]]
  }
  VQK <- list()
  for(k in 1:10){VQK[[k]] = eigen(MNK[[k]])$vectors[,1:q]}
  squarE <- list()
  for(k in 1:10){
    squarE[[k]] = 1/(d-q)*sum(eigen(MNK[[k]])$values[qplus:d])
  }
  WQK = list()
  for(k in 1:10){
    WQK[[k]] = matrix(0,q,q)
    for( i in 1:q){
      WQK[[k]][i,i] = sqrt(eigen(MNK[[k]])$values[i]-squarE[[k]])
    }
    WQK[[k]] = VQK[[k]] %*% WQK[[k]]
  }
  Ksigma <- list()
  for(k in 1:10){
    Ksigma[[k]] = WQK[[k]] %*% t(WQK[[k]]) + squarE[[k]] * diag(d)
  }
}

# plot the iteration no. vs loglikelihood graph
piclog = loglik[1:25,]
iteration = length(piclog)
plot(x=1:iteration, y=piclog, type="o",xlab="iteration no.",ylab = "loglikelihood", main = "p=0")
AIC = -2*loglik[25,] + 2*(d*q+1-q*(q-1)/2)
####################################################
### COUNT THE MIS-CATEGORIZATION RATE
myLabel=apply(myData[,257:266],1,function(xx){
  return(which(xx=="1")-1)
})
myNumber = data.matrix(myData[,257:266])
col<-colSums(myNumber)
show(col)
colno<-vector(mode="numeric",length=0)
# k column in gama
for(i in 1:10){
  a<-colSums(gamaIK[myLabel==i-1,])
  colno[i]<-which(a==max(a),arr.ind=T)
}
colcluster = matrix(0,1593,1)
for(i in 1:1593){
  colcluster[i,] = which.max(gamaIK[i,])
}
newLabel = matrix(0,1593,10)
for(i in 1:1593){
  newLabel[i,colcluster[i,]] = 1
}
#### count the remaining 1(which is the wrong numbers)
wrongNO = matrix(0,10,1)
for(i in 1:10){
  b = matrix(myNumber[,i])
  bplus = colno[i]
  c = matrix(newLabel[,bplus])
  d = b-c
  wrongNO[i,] = length(which(d[,1]==1))
}
#### wrong rate
wrongrate = matrix(0,10,1)
for(i in 1:10){
  wrongrate[i,1] = wrongNO[i,1]/col[i]
}
### PLOT IMAGE
listcluster <- list()
for(k in 1:10){
  kplus = colno[k]
  listcluster[[k]] = rbind(myX[which(colcluster[,1]==kplus),])
}
par(mai = c(0.03,0.03,0.03,0.03),mfrow = c(10,6))
for(kcol in 1:10){
  kplus = colno[kcol]
  image(t(matrix(miuK[kplus,], byrow=TRUE,16,16)[16:1,]),col=gray(0:100/100),axes=FALSE)
  indices = sample(dim(listcluster[[kcol]])[1], 5)
  for(i in indices){
    image(t(matrix(listcluster[[kcol]][i,], byrow=TRUE,16,16)[16:1,]),col=gray(0:1/1),axes=FALSE)
  }
}
