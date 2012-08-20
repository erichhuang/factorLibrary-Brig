library(mGenomics)
library(synapseClient)
library(snm)

ent <- loadEntity('syn138518')
#fits <- runWorkflow(ent$cacheDir, workflow="snm")
#save(fits,file="unsupervisedFits.Rda")

newEnt <- Data(list(name="ectopicMYC", 
				parentId = "syn275939", 
				platform=names(fits),						
				tissueType="breast",
				numSamples=ncol(exprs(fits[[1]][[1]])),
				species="Homo sapiens"))
annotValue(newEnt,'processingAlgorithm') <- 'snmUnsupervised'
annotValue(newEnt,'cellLine') <- 'MCF'
annotValue(newEnt,'perturbation') <- 'MYC'
newEnt <- addObject(newEnt, as.list(fits[[1]]), unlist=TRUE)
newEnt <- createEntity(newEnt)
newEnt <- storeEntity(newEnt)

ent <- loadEntity("syn319533")
dat <- exprs(ent$objects$eset)

# Perform a differential expression test.
tmp <- ent$objects$eset@protocolData@data$ScanDate
dates <- sapply(strsplit(tmp," "), function(x){ x[1]})
times <- sapply(strsplit(tmp," "), function(x){ x[2]})
perturbation <- ifelse(grepl('M', colnames(dat)), "MYC","GFP")
obj <- data.frame(perturbation=perturbation,
		scanDates=dates,
		scanTimes=times)
rownames(obj) <- colnames(dat)
write.table(obj,file="ectopicMYC_metadata.txt",sep="\t",quote=FALSE)

mdEnt <- Data(list(name="ectopicMYC_metadata", 
				parentId = propertyValue(newEnt,'parentId'), 
				platform=propertyValue(newEnt, 'platform'),
				tissueType=propertyValue(newEnt, 'tissueType'),
				disease=propertyValue(newEnt, 'disease'),
				numSamples=propertyValue(newEnt, 'numSamples'),
				species=propertyValue(newEnt, 'species')))
annotations(mdEnt) <- annotations(newEnt)
addObject(mdEnt, obj, "ectopicMYC_metadata")
addFile(mdEnt, "ectopicMYC_metadata.txt")
mdEnt <- createEntity(mdEnt)
mdEnt <- storeEntity(mdEnt)

mdEnt <- loadEntity('syn319619')
table(mdEnt$objects$ectopicMYC_metadata$perturbation, mdEnt$objects$ectopicMYC_metadata$scanDates)
X <- model.matrix(~ factor(perturbation))
sig <- calcSig(dat, X)

myc <- names(which(sig$cfs[,2] > 0.5))
library(mg.hgu133plus2.db); 
hsym <- as.character(mg.hgu133plus2SYMBOL)
hsym[1:5]
load("~/Desktop/studies.Rda")
load("~/Documents/Randoms/gic_explore.Rda")
myc.sig <- calc.signature.stats("",gpl570.gic, sort(as.character(hsym[myc])))
studies[names(which(myc.sig[2,] > 0.5)),1:2]

rownames(dat) <- hsym[rownames(dat)]
th <- intersect(rownames(gpl570.gic), rownames(dat))
dat2 <- dat[th,]
gpl2 <- gpl570.gic[th,]

sig <- calcSig(dat2, X)
Z <- model.matrix(~ myc.sig[2,])
sig2 <- calcSig(gpl2, Z)





qry <- synapseQuery('select id, name from entity where entity.parentId == "syn301181" and entity.name=="GSE10070_processed"')
newEnt <- loadEntity(qry$entity.id)
gics  <- newEnt$objects$hgu133plus2$statistics$singular.values[,1]
lab <- rep("Else",length(gics)); 
names(lab) <- names(gics); 
lab[myc] <- "In MYC Signature"
densityplot(gics,groups=lab,
		auto.key=list(columns=2))

new.exprs <- exprs(newEnt$objects$hgu133plus2$eset)
u <- fs(new.exprs[myc,])
par(mfrow=c(1,3))
plot(u$d)
plot(u$v[,1])
plot(u$u[,1])
cors <- cor(t(new.exprs[myc,]))


#GSE10739
qry <- synapseQuery('select id, name from entity where entity.parentId == "syn301181" and entity.name=="GSE10739_processed"')
newEnt <- loadEntity(qry$entity.id)
gics  <- newEnt$objects$hgu133plus2$statistics$singular.values[,1]
lab <- rep("Else",length(gics)); 
names(lab) <- names(gics); 
lab[myc] <- "In MYC Signature"
densityplot(gics,groups=lab,
		auto.key=list(columns=2))

new.exprs <- exprs(newEnt$objects$hgu133plus2$eset)
u <- fs(new.exprs[myc,])
par(mfrow=c(1,3))
plot(u$d)
plot(u$v[,1])
plot(u$u[,1])
cors <- cor(t(new.exprs[myc,]))

#GSE10847
qry <- synapseQuery('select id, name from entity where entity.parentId == "syn301181" and entity.name=="GSE10847_processed"')
newEnt <- loadEntity(qry$entity.id)
gics  <- newEnt$objects$hgu133plus2$statistics$singular.values[,1]
lab <- rep("Else",length(gics)); 
names(lab) <- names(gics); 
lab[myc] <- "In MYC Signature"
densityplot(gics,groups=lab,
		auto.key=list(columns=2))

new.exprs <- exprs(newEnt$objects$hgu133plus2$eset)
u <- fs(new.exprs[myc,])
par(mfrow=c(1,3))
plot(u$d)
plot(u$v[,1])
plot(u$u[,1])
cors <- cor(t(new.exprs[myc,]))

gses <- as.character(unlist(studies[names(which(myc.sig[2,] > 0.5)),1]))

lapply(gses[7:length(gses)],function(gse){
			cat(gse,"\t")
			gse <- paste(gse,'_processed',sep="")
			qry <- synapseQuery(paste('select id, name from entity where entity.parentId == "syn301181" and entity.name=="',gse,'"',sep=""))
			if(is.null(qry)){
				cat("not found\n")
				return(list());
			}
			cat("found\n")
			newEnt <- loadEntity(qry$entity.id)
			gics  <- newEnt$objects$hgu133plus$statistics$singular.values[,1]
			lab <- rep("Else",length(gics)); 
			names(lab) <- names(gics); 
			lab[myc] <- "In MYC Signature"
			png(file=paste(gse,"_densityplot.png",sep=""),width=960)
			print(densityplot(gics,groups=lab,
					auto.key=list(columns=2)))
	dev.off()
			
			new.exprs <- exprs(newEnt$objects$hgu133plus2$eset)
			u <- fs(new.exprs[myc,])
			png(file=paste(gse,"_svd.png",sep=""),width=960)
			par(mfrow=c(1,3))
			plot(u$d)
			plot(u$v[,1])
			plot(u$u[,1])
			dev.off()
			u
		}) -> us 


# Take an SVD of the data.  Look at eigenweights and first couple of eigengenes
u.full <- fs(dat)

png(file='pvalues.png')
hist(sig$pval, xlab="P values", ylab="Ras Effect on Gene Expression Variation")
dev.off()

png(file="u_unsup_d.png")
barplot(u.full$d, ylab="Prop Variance Explained by Each Eigengene")
dev.off()

png(file="eg1.png")
plot(u.full$v[,1], xlab="Samples",ylab="Eigengene 1")
dev.off()

png(file="eg2.png")
plot(u.full$v[,2], xlab="Samples",ylab="Eigengene 2")
dev.off()

a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'pvalues.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'u_unsup_d.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'eg1.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'eg2.png'); system(a)

# Explore workflow output. Calculate the total sums of squares of the data.
dat.m <- rowMeans(dat)
dat.c <- sweep(dat,1,dat.m)
tSSQ_dat <- sum(dat.c^2)

# Now we want to remove the effects of the biological treatment
# on the data.  The reason is we are interested in identifying 
# any latent structure. In order to do so we calcualte the
# residual sums of squares and take a singular value decomposition
# of the data.  Note the singular values are weighted to sum to 1.
res <- dat - t(X %*% solve(t(X) %*% X) %*% t(X) %*% t(dat))
rSSQ_dat <- sum(res^2)
u <- fs(res)

png(file="pVE_rEGs.png")
barplot(round(rSSQ_dat * u$d,3) / tSSQ_dat, ylab="Prop tSSQ Explained by Each Eigengene")
dev.off()

a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'pVE_rEGs.png'); system(a)

sva.fit <- sva(dat, bio.var=X, n.sv=1, num.iter=30, diagnose=FALSE)
# Now, we'll take a look at the estimated basis vectors
png(file="svaEGs.png")
#par(mfrow=c(2,1))
plot(sva.fit$svd[[30]]$v[,1], xlab="Samples", ylab="First Basis Vector")
#plot(sva.fit$svd[[30]]$v[,2], xlab="Samples", ylab="Second Basis Vector")
dev.off()

a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'svaEGs.png'); system(a)

# Get the null probes
Z <- model.matrix(~ sva.fit$svd[[30]]$v[,1])
dat2 <- snm(dat, bio.var=X, adj.var=Z, rm.adj=TRUE)$norm.dat
sig2 <- calcSig(dat2, X, Z)

png(file="supPvalues.png")
hist(sig2$pval)
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'supPvalues.png'); system(a)

u2 <- fs(dat2)
png(file="supD.png")
barplot(u2$d, ylab="Prop Variance Explained by Each Eigengene")
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'supD.png'); system(a)

png(file="eg1_sup.png")
plot(u2$v[,1], xlab="Samples",ylab="Eigengene 1")
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"MYC\" \"text/plain\" ", 'eg1_sup.png'); system(a)

# Make the synapse object
new.obj <- buildmGenomicsObject(dat2, 'hthgu133a', fits$hthgu133a$statistics)
newEnt <- Data(list(name = paste(propertyValue(ent, "name"),"_processed",sep=""),
				parentId = 'syn275939'))
annotations(newEnt) <- annotations(ent)

eset <- new.obj$eset
statistics <- new.obj$statistics
annotValue(newEnt, "processingAlgorithm") <- "snmSupervised"
# Attach eset and statistics
newEnt <- addObject(newEnt, eset)
newEnt <- addObject(newEnt, statistics)
# Create new entity in synapse
newEnt <- createEntity(newEnt)
newEnt <- storeEntity(newEnt)






