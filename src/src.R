library(mGenomics)
library(synapseClient)
library(snm)

ent <- loadEntity('syn138525')
#fits <- runWorkflow(ent$cacheDir, workflow="snm")

newEnt <- Data(list(name="ectopicSRC", 
				parentId = "syn275939", 
				platform=names(fits),						
				tissueType="breast",
				numSamples=ncol(exprs(fits[[1]][[1]])),
				species="Homo sapiens"))
annotValue(newEnt,'processingAlgorithm') <- 'snmUnsupervised'
annotValue(newEnt,'cellLine') <- 'MCF'
annotValue(newEnt,'perturbation') <- 'SRC'
newEnt <- addObject(newEnt, as.list(fits[[1]]), unlist=TRUE)
newEnt <- createEntity(newEnt)
newEnt <- storeEntity(newEnt)

ent <- loadEntity("syn317803")
dat <- exprs(ent$objects$eset)

# Perform a differential expression test.
tmp <- ent$objects$eset@protocolData@data$ScanDate
dates <- sapply(strsplit(tmp," "), function(x){ x[1]})
times <- sapply(strsplit(tmp," "), function(x){ x[2]})
perturbation <- ifelse(grepl('GFP', colnames(dat)), "GFP","SRC")
obj <- data.frame(perturbation=perturbation,
		scanDates=dates,
		scanTimes=times)
rownames(obj) <- colnames(dat)
write.table(obj,file="ectopicSRC_metadata.txt",sep="\t",quote=FALSE)

mdEnt <- Data(list(name="ectopicSRC_metadata", 
				parentId = propertyValue(newEnt,'parentId'), 
				platform=propertyValue(newEnt, 'platform'),
				tissueType=propertyValue(newEnt, 'tissueType'),
				disease=propertyValue(newEnt, 'disease'),
				numSamples=propertyValue(newEnt, 'numSamples'),
				species=propertyValue(newEnt, 'species')))
annotations(mdEnt) <- annotations(newEnt)
addObject(mdEnt, obj, "ectopicSRC_metadata")
addFile(mdEnt, "ectopicSRC_metadata.txt")
mdEnt <- createEntity(mdEnt)
mdEnt <- storeEntity(mdEnt)

mdEnt <- loadEntity('syn317806')
table(mdEnt$objects$ectopicSRC_metadata$perturbation, mdEnt$objects$ectopicSRC_metadata$scanDates)


fits2 <- runWorkflow(ent$cacheDir, workflow="snm", bio.var=X, adj.var=model.matrix(~ u$v[,1]), rm.adj=TRUE)


# Get the data
# Perform a differential expression test.
X <- model.matrix(~ factor(perturbation))
sig <- calcSig(dat, X)
# Take an SVD of the data.  Look at eigenweights and first couple of eigengenes
u.full <- fs(dat)

png(file='pvalues.png')
hist(sig$pval, xlab="P values", ylab="SRC on Gene Expression Variation")
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

a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'pvalues.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'u_unsup_d.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'eg1.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'eg2.png'); system(a)

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
a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'pVE_rEGs.png'); system(a)

sva.fit <- sva(dat, bio.var=X, n.sv=2, num.iter=30, diagnose=TRUE)
# Now, we'll take a look at the estimated basis vectors
png(file="svaEGs.png")
plot(sva.fit$svd[[30]]$v[,1])
plot(sva.fit$svd[[30]]$v[,2])
plot(sva.fit$svd[[30]]$v[,3])
plot(sva.fit$svd[[30]]$v[,4])
dev.off()

a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'svaEGs.png'); system(a)



ids <- which(sig$pval > 0.6)
u2 <- fs(dat[ids,])
snm.fit2 <- snm(dat, 
		bio.var=X, 
		adj.var=model.matrix(~u2$v[,1:2]), 
		rm.adj=TRUE)
sig2 <- calcSig(snm.fit2$norm.dat,X)
u3 <- fs(snm.fit2$norm.dat)

png(file="prePost.png")
par(mfrow=c(2,2))
plot(u$v[,1], main="Unsupervised First Eigengene")
plot(u3$v[,1], main="Supervised First Eigengene")
plot(sig$cfs[,2], sig2$cfs[,2], xlab="Unsupervised SRC Effect", ylab="Supervised SRC Effect")
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'prePost.png'); system(a)

load("~/Documents/Randoms/gic_explore.Rda")
sign <- hsym[names(sort(-sig$cfs[,2]))[1:300]]
src.sig <- calc.signature.stats("",gpl570.gic, sign)

png(file="ks.png")
plot(src.sig[2,], xlab="Studies", ylab="Signature Enrichment")
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'ks.png'); system(a)

rownames(dat) <- hsym[rownames(dat)]
th <- intersect(rownames(gpl570.gic), rownames(dat))
dat2 <- dat[th,]
gpl2 <- gpl570.gic[th,]

sig <- calcSig(dat2, X)
Z <- model.matrix(~ src.sig[2,])
sig2 <- calcSig(gpl2, Z)


ent <- loadEntity('syn138525')

library(affy)
abatch <- ReadAffy(filenames=list.files(ent$cacheDir,full.names=TRUE))
int.var <- bio.var <- adj.var <- gene2row <- NULL
rm.adj <- TRUE
int <- exprs(abatch)
# Gets the gene2row object that helps us define probe sets.
gene2row.cdf <- getGene2Row(gene2row, annotation(abatch))
# Loads the pm data and calls snm
pms <- int[unlist(gene2row.cdf), ]
pms[pms <= 1] <- 1
data <- log2(pms)
if(is.null(int.var)) {
	int.var <- data.frame(array = factor(1:ncol(data)))
}
# Calls snm
if(verbose) {cat("Normalizing Data\n")}
sendLogMessage("Normalizing Data",logFile)
snm.fit <- snm(data, bio.var, adj.var, int.var=int.var, diagnose=FALSE, rm.adj=rm.adj,verbose=FALSE)
colnames(snm.fit$norm.dat) <- sampleNames(abatch)


int.var <- data.frame(scanDate=factor(ifelse(dates=='03/11/04',1,2)))
snm.fit.2 <- snm(data, bio.var, adj.var, int.var=int.var, diagnose=FALSE, rm.adj=rm.adj,verbose=FALSE)

gene2row.tmp <- split(1:nrow(pms), rep(names(gene2row.cdf), sapply(gene2row.cdf,length)))

# Calls EPSA
if(verbose) {cat("Summarizing Data\n")}
sendLogMessage("Summarizing Data",logFile)	
fits <- fit.pset(snm.fit$norm.dat, gene2row.tmp)











