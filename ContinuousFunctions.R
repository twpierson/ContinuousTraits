library(ape) #utility fns
library(geiger) #utilty fns
library(OUwie)

#You can use code you wrote for the correlation exercise here.


VisualizeData <- function(data) {
	pdf("Visualize.pdf")
	par(mfrow=c(1,2))
	plot(data[[1]])
	hist(data[[2]])
	dev.off()
}

CleanData <- function(phy, data) {
	treedata(phy,data)
}

RunSingleOUwieModel<-function(model, phy, data) {
	print(paste("Now starting model",model))
	return(OUwie(phy, data, model, simmap.tree=FALSE, diagn=FALSE))	
}
