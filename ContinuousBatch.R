#You can use code you wrote for the correlation exercise here.
source("ContinuousFunctions")
setwd("~/Desktop/UTK/Spring_2016/PhyloMeth/ContinuousTraits")

tree <- read.tree("Eurycea_Tree")

# Here's the old, bad data. 
#discrete.data <- round(runif(length(tree$tip.label)))
#names(discrete.data)<-tree$tip.label

#continuous.data <- rnorm(nTips(tree),mean=25,sd=8)
#names(continuous.data)<-tree$tip.label

# Here I'm simulating a discrete trait.
q <- list(rbind(c(-12,12), c(12,-12)))
sim.discrete.1 <- sim.char(tree, q, model='discrete', n=1)
sim.discrete.2 <- sim.char(tree, q, model='discrete', n=1)
# Here I'm reformatting those trait data so that they work in downstream analyses. The format they were in had a weird string of text at the top that was causing problems.
new.discrete.1 <- as.vector(sim.discrete.1)
new.discrete.2 <- as.vector(sim.discrete.2)
names(new.discrete.1) <- row.names(sim.discrete.1)
names(new.discrete.2) <- row.names(sim.discrete.2)

full.discrete <- data.frame(new.discrete.1,new.discrete.2)

# Here I'm simulating a continuous trait.
sim.continuous.1 <- sim.char(tree, 1, n=1)
sim.continuous.2 <- sim.char(tree, 10000, root=50, n=1)
# Here I'm reformatting those trait data so that they work in downstream analyses. The format they were in had a weird string of text at the top that was causing problems.
new.continuous.1 <- as.vector(sim.continuous.1)
new.continuous.2 <- as.vector(sim.continuous.2)
names(new.continuous.1) <- row.names(sim.continuous.1)
names(new.continuous.2) <- row.names(sim.continuous.2)

full.continuous <- data.frame(new.continuous.1,new.continuous.2)

cleaned.discrete <- CleanData(tree, full.discrete)
cleaned.continuous <- CleanData(tree, full.continuous)

VisualizeData(cleaned.continuous)
VisualizeData(cleaned.discrete)

#First, start basic. What is the rate of evolution of your trait on the tree? 

BM1 <- fitContinuous(cleaned.continuous$phy, cleaned.continuous$data[,2], model="BM")
print(paste("The rate of evolution is", BM1[[4]]$sigsq, "in units of", "mm per x fraction of a nucleotide substitution"))
#Important: What are the rates of evolution? In what units?

OU1 <- fitContinuous(tree, cleaned.continuous$data[,2], model="OU")

par(mfcol=c(1,2))
plot(tree, show.tip.label=FALSE)
ou.tree <- rescale(tree, model="OU", alpha=OU1[[4]]$alpha)
plot(ou.tree, show.tip.label=FALSE)
#How are the trees different?
# They're very similar, but some branch lengths differ slightly.

#Compare trees
AIC.BM1 <- BM1[[4]]$aic
AIC.OU1 <- OU1[[4]]$aic
delta.AIC.BM1 <- AIC.BM1-min(c(AIC.BM1,AIC.OU1))
delta.AIC.OU1 <- AIC.OU1-min(c(AIC.BM1,AIC.OU1))

#OUwie runs:
#This takes longer than you may be used to. 
#We're a bit obsessive about doing multiple starts and in general
#performing a thorough numerical search. It took you 3+ years
#to get the data, may as well take an extra five minutes to 
#get an accurate answer.

#First, we need to assign regimes. The way we do this is with ancestral state estimation of a discrete trait.
#We can do this using ace() in ape, or similar functions in corHMM or diversitree. Use only one discrete char
one.discrete.char <- new.discrete.1
reconstruction.info <- ace(one.discrete.char, tree, type="discrete", method="ML", CI=TRUE)
best.states <- apply(reconstruction.info$lik.anc, 1, which.max)

#NOW ADD THESE AS NODE LABELS TO YOUR TREE

labeled.tree <- tree
labeled.tree$node.label<-best.states

tips<-rownames(cleaned.continuous$data)
newer.continuous <- data.frame(tips,cleaned.discrete$data[,1],cleaned.continuous$data[,1])
colnames(newer.continuous) <- c("tips","regime","data")

nodeBased.OUMV <- OUwie(labeled.tree, newer.continuous,model="OUMV", simmap.tree=FALSE, diagn=FALSE)
print(nodeBased.OUMV)
#What do the numbers mean?

#Now run all OUwie models:
models <- c("BM1","BMS","OU1","OUMV","OUM","OUMA","OUMVA")
results <- lapply(models, RunSingleOUwieModel, phy=labeled.tree, data=newer.continuous)
AICc.values<-sapply(results, "[[", "AICc")
names(AICc.values)<-models
AICc.values<-AICc.values-min(AICc.values)

print(AICc.values) #The best model is the one with smallest AICc score

best<-results[[which.min(AICc.values)]] #store for later

print(best) #prints info on best model


#We get SE for the optima (see nodeBased.OUMV$theta) but not for the other parameters. Let's see how hard they are to estimate. 
#First, look at ?OUwie.fixed to see how to calculate likelihood at a single point.
?OUwie.fixed

#Next, keep all parameters but alpha at their maximum likelihood estimates (better would be to fix just alpha and let the others optimize given this constraint, but this is harder to program for this class). Try a range of alpha values and plot the likelihood against this.
alpha.values<-seq(from= 1 , to= 100 , length.out=50)

#keep it simple (and slow) and do a for loop:
likelihood.values <- rep(NA, length(alpha.values))
for (iteration in sequence(length(alpha.values))) {
	likelihood.values[iteration] <- OUwie.fixed(labeled.tree, newer.continuous, model="OUMA", alpha=rep(alpha.values[iteration],2), sigma.sq=best$solution[2,], theta=best$theta[,1])$loglik
}

plot(x= alpha.values , y= likelihood.values, xlab="Alpha Values", ylab="Log.Likelihood", type="l", bty="n")

points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
text(x=best$solution[1,1], y=best$loglik, "unconstrained best", pos=4, col="red")

#a rule of thumb for confidence for likelihood is all points two log likelihood units worse than the best value. Draw a dotted line on the plot to show this
abline(h=best$loglik-2, lty="dotted") #Two log-likelihood 


#Now, let's try looking at both theta parameters at once, keeping the other parameters at their MLEs
require("akima")

nreps<-400
theta1.points<-c(best$theta[1,1], rnorm(nreps-1, best$theta[1,1], 5*best$theta[1,2])) #center on optimal value, have extra variance
theta2.points<-c(best$theta[2,1], rnorm(nreps-1, best$theta[2,1], 5*best$theta[2,2])) #center on optimal value, have extra variance
likelihood.values<-rep(NA,nreps)

for (iteration in sequence(nreps)) {
	likelihood.values[iteration] <- OUwie.fixed(labeled.tree, newer.continuous, model="OUMA", alpha=best$solution[1,], sigma.sq=best$solution[2,], theta=c(theta1.points[iteration], theta2.points[iteration]))$loglik
}
#think of how long that took to do 400 iterations. Now remember how long the search took (longer).

likelihood.differences<-(-(likelihood.values-max(likelihood.values)))

#We are interpolating here: contour wants a nice grid. But by centering our simulations on the MLE values, we made sure to sample most thoroughly there
interpolated.points<-interp(x=theta1.points, y=1000000000000000*theta2.points, z= likelihood.differences, linear=FALSE, extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), yo=seq(min(1000000000000000*theta2.points), max(1000000000000000*theta2.points), length = 400))
	
contour(interpolated.points, xlim=range(c(theta1.points, 100000000000000*theta2.points)),ylim=range(c(theta1.points, 100000000000000*theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)

points(x=newer.continuous$X[which(new.continuous$Reg==1)],y=rep(min(c(theta1.points, theta2.points)), length(which(new.continuous$Reg==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
points(y=newer.continuous$X[which(new.continuous$Reg==2)],x=rep(min(c(theta1.points, theta2.points)), length(which(new.continuous$Reg==2))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis


#The below only works if the discrete trait rate is low, so you have a good chance of estimating where the state is.
#If it evolves quickly, hard to estimate where the regimes are, so some in regime 1 are incorrectly mapped in
#regime 2 vice versa. This makes the models more similar than they should be.
#See Revell 2013, DOI:10.1093/sysbio/sys084 for an exploration of this effect.
library(phytools)
trait.ordered<-data.frame(trait[,2], trait[,2],row.names=trait[,1])
trait.ordered<- trait.ordered[tree$tip.label,]
z<-trait.ordered[,1]
names(z)<-rownames(trait.ordered)
tree.mapped<-make.simmap(tree,z,model="ER",nsim=1)
leg<-c("black","red")
names(leg)<-c(1,2)
plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)

simmapBased<-OUwie(tree.mapped,trait,model="OUMV", simmap.tree=TRUE, diagn=FALSE)
print(simmapBased)
#How does this compare to our best model from above? Should they be directly comparable?
print(best)