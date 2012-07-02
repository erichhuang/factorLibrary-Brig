library(lattice)
load("betaCatenin/unsupervisedFits.Rda")
bcat <- fits
load("egfr/unsupervisedFits.Rda")
egfr <- fits
load("pi3k/unsupervisedFits.Rda")
pi3k <- fits
load("src/unsupervisedFits.Rda")
src <- fits
load("e2f3/unsupervisedFits.Rda")
e2f3 <- fits
load("myc/unsupervisedFits.Rda")
myc <- fits
load("ras/unsupervisedFits.Rda")
ras <- fits

all <- list(bcat=bcat, egfr=egfr, pi3k=pi3k, src=src, e2f3=e2f3, myc=myc, ras=ras)
gics <- sapply(all, function(x){ x[[1]]$statistics$singular.values[,1]})
nms <- sapply(gics, function(x){ names(x)})
in.common <- names(which(table(unlist(nms)) == 7))

mat <- matrix(0, nr=length(in.common), ncol=length(all))
rownames(mat) <- in.common
colnames(mat) <- names(all)

for(i in 1:length(gics)){
	mat[in.common,i] <- gics[[i]][in.common]
}

obj <- 
data.frame(fic=as.numeric(mat),
		oncogene = rep(colnames(mat),each=nrow(mat)))
png(file="allFIC.png")
densityplot(obj$fic, 
		groups=obj$oncogene, 
		plot.points=FALSE,
		auto.key=list(columns=4))
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"Factor Library\" \"text/plain\" ", 'allFIC.png'); system(a)
