### Iterator - simulates network expression values over 20 iterations based on an given 
### interaction network, and a randomly generated initial gene expression status. Returns
### a list that contains the interaction network, network stability, and expression values
load("/home/shea/Quick.RData")
iterator <- function(interactionNetwork,nIters,expression){
  nGenes<-nrow(interactionNetwork)
  newExpression <- expression
  for(i in 1:nIters){
##       multiplication <- matrix(rep.int(newExpression, times=10),10,10)
##      result <- interactionNetwork * multiplication
##      x <- colSums(result) 
### this direction of multiplication means the interactionNetwork is the influence of rows on columns
      x <- t(t(newExpression)%*%interactionNetwork)
      newExpression <- ((tanh(2*x)+1)/2)
      expression <- newExpression 
  }
  
  return(list(interactionNetwork = interactionNetwork, expression=expression))
}

resume<-function(){
  while(i<expSpan){
    fitnesses<- life(population, expressionOne, expressionTwo, popSize, nIterations, para1, para2)
    aveFit<-rbind(aveFit,mean(fitnesses))
    bestFit<-rbind(bestFit,max(fitnesses,na.rm = FALSE))
    offspring=sample(1:popSize,popSize,replace=T,prob=fitnesses)
    trackOffspring<-cbind(trackOffspring,offspring)
    population<-population[,,offspring]
    lineage<-lineage[offspring]
    allLineages<- cbind(allLineages,lineage)
    i<-i+1
    newGen<-evolution(population, popSize, i)
    population<-newGen$newPop
    mutation<-rbind(mutation,newGen$value)
    if(i==5){
      test<-population
    }
  }
}

life<-function(population, expressionOne, expressionTwo, popSize, nIterations, para1, para2){
  ends<-iterator(population[,,1], nIterations, expressionOne)$expression
  for(i in 2:popSize){
    end<-iterator(population[,,i], nIterations, expressionOne)$expression
    ends<-cbind(ends,end)
  }
  endsTwo<-iterator(population[,,1], nIterations, expressionTwo)$expression
  for(i in 2:popSize){
    endTwo<-iterator(population[,,i], nIterations, expressionTwo)$expression
    endsTwo<-cbind(endsTwo,endTwo)
  }
  fitnesses<-fitness(ends[,1],endsTwo[,1],para1,para2)
  for(i in 2:popSize){
    fitnesses<- c(fitnesses,fitness(ends[,i],endsTwo[,i], para1, para2))
  }
  return(fitnesses)
} 

change<- function(expressionOne){
  nGenes<-nrow(expressionOne)
  expressionTwo <-matrix(c(expressionOne),nrow=nGenes, ncol=1)
  changeGene <-10
  print(changeGene)
  expressionTwo[changeGene,]<-0.85
  if(expressionTwo[changeGene,]==expressionOne[changeGene,]){
    expressionTwo[changeGene,]<-0.35
  }
  return (expressionTwo)
}
evolution <- function (arg1, popSize, i)
{
  rBoolean<- array(runif(nGenes*nGenes*popSize,0,1), dim=c(nGenes,nGenes,popSize))<0.0005
  index<-which(rBoolean)
  rValue<-rBoolean*rnorm(nGenes*nGenes*popSize,0,sd=0.05)
  value<-rValue[rBoolean]
  gen<- rep(i,length(value))
  value<-cbind(value,index)
  value<-cbind(value,gen)
  newPop<-arg1+(rValue)
  return(list(newPop=newPop, value=value))
}

fitness <- function(X,Y,para1,para2){

  distance<-0
  distance<-sqrt(((sum(abs(X-para1)))^2)+((sum(abs(Y-para2)))^2))
  return(1/(distance+1))
}

nGenes=20
nIterations=100
expSpan= 100000
popSize=3000
para1<-c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0)
print(para1)
para2<-c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)
print(para2)
population<-originalPopulation
expressionTwo <-change(expressionOne)
lineage <- c(1:popSize)
allLineages<-lineage
test<-originalPopulation
mutation<-matrix(c(0,0,0),nrow=1,ncol=3)
i=0
trackOffspring<-lineage
fitnesses<- life(population, expressionOne, expressionTwo, popSize, nIterations, para1, para2)
aveFit<-c(mean(fitnesses))
bestFit<-c(max(fitnesses, na.rm = FALSE))
while(i<expSpan){
  fitnesses<- life(population, expressionOne, expressionTwo, popSize, nIterations, para1, para2)
  aveFit<-rbind(aveFit,mean(fitnesses))
  bestFit<-rbind(bestFit,max(fitnesses,na.rm = FALSE))
  offspring=sample(1:popSize,popSize,replace=T,prob=fitnesses)
  trackOffspring<-cbind(trackOffspring,offspring)
  population<-population[,,offspring]
  lineage<-lineage[offspring]
  allLineages<- cbind(allLineages,lineage)
  i<-i+1
  newGen<-evolution(population, popSize, i)
  population<-newGen$newPop
  mutation<-rbind(mutation,newGen$value)
  if(i%%3000==0){
    save(allLineages,aveFit,bestFit, expSpan, popSize,expressionOne,expressionTwo,i,nGenes,fitnesses,originalPopulation,population,trackOffspring, mutation, nIterations, para1, para2, file="/home/shea/ProgressionControl.RData")
  }
}

incFit<-unlist(as.matrix(sort.int(fitnesses, index.return=TRUE))[2])
finalMatrix<- c(population[,,incFit[1]])
for(i in 2:popSize){
  finalMatrix<-rbind(finalMatrix,c(population[,,incFit[i]]))
}
repopulation<-function(number)
{
  newPop<-originalPopulation
  j<-2
  for(i in 1:number){
    mutes<-rep(0, nGenes*nGenes*popSize)
    while(i== mutation[j,3]){
      mutes[mutation[j,2]]<-mutation[j,1]
      j<-j+1
    }
    newPop<-newPop+mutes
    offspring<-trackOffspring[,i+1]
    newPop<-newPop[,,offspring]
  }
  return (newPop)
}
development<-function(population, expressionOne, expressionTwo, popsize, nIterations){
  ends<-sum(iterator(population[,,1], nIterations, expressionOne)$expression[,nIterations+1])
  for(i in 2:popSize){
    end<-sum(iterator(population[,,i], nIterations, expressionOne)$expression[,nIterations+1])
    ends<-cbind(ends,end)
  }
  endsTwo<-iterator(population[,,1], nIterations, expressionTwo)$expression[,nIterations+1]
  for(i in 2:popSize){
    endTwo<-sum(iterator(population[,,i], nIterations, expressionTwo)$expression[,nIterations+1])
    endsTwo<-cbind(endsTwo,endTwo)
  }
  plot(ends,type="l")
  plot(endsTwo,type="l")
}
save(allLineages,aveFit,bestFit, expSpan, popSize,expressionOne,expressionTwo,i,nGenes,fitnesses,originalPopulation,population,trackOffspring, mutation, nIterations, para1, para2, file="/home/shea/finalCon.RData")

nGenes<-20
from<-c(1:(nGenes*nGenes))
to<-c(1:(nGenes*nGenes))
weight<-c(1:(nGenes*nGenes))
for(i in 1:nGenes){
  for(j in 1:nGenes){
    from[((i-1)*nGenes)+j]<-i
    to[((i-1)*nGenes)+j]<-j
    weight[((i-1)*nGenes)+j]<-originalPopulation[i,j,1]
    j<-j+1
  }
  i<-i+1
}
dfEdges= data.frame(from, to, weight)
g1<-graph_from_data_frame(dfEdges, directed = TRUE, vertices = NULL)
circle<-layout_in_circle(g1, order=1:20)
g1 <- delete.edges(g1, E(g1)[ abs(weight) < 0.3 ])
plot.igraph(g1, vertex.size=20, edge.width=abs(E(g1)$weight)*8, 
            edge.color=ifelse(weight > 0, "blue","red"), layout=circle)


for(i in 1:nGenes){
  for(j in 1:nGenes){
    from[((i-1)*nGenes)+j]<-i
    to[((i-1)*nGenes)+j]<-j
    weight[((i-1)*nGenes)+j]<-population[i,j,1]
    j<-j+1
  }
  i<-i+1
}
dfEdges= data.frame(from, to, weight)
g1<-graph_from_data_frame(dfEdges, directed = TRUE, vertices = NULL)
circle<-layout_in_circle(g1, order=1:20)
g1 <- delete.edges(g1, E(g1)[ abs(weight) < 0.3 ])
plot.igraph(g1, vertex.size=20, edge.width=abs(E(g1)$weight)*8, 
            edge.color=ifelse(weight > 0, "blue","red"), layout=circle)


plot(bestFit,type="l")
lines(aveFit,col="red")
km=kmeans(originalPopulation[,,1],2)
print(km$betweenss/km$totss)
km=kmeans(population[,,1],2)
print(km$betweenss/km$totss)

