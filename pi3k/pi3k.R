library(mGenomics)
library(synapseClient)
library(snm)
synapseLogin('brig.mecham@sagebase.org','letmein')

ent <- loadEntity('syn138520')
fits <- runWorkflow(ent$cacheDir, workflow="snm")

newEnt <- Data(list(name="ectopicPI3K", 
				parentId = "syn275939", 
				platform=names(fits),						
				tissueType="breast",
				numSamples=ncol(exprs(fits[[1]][[1]])),
				species="Homo sapiens"))
annotValue(newEnt,'processingAlgorithm') <- 'snmUnsupervised'
annotValue(newEnt,'cellLine') <- 'MCF'
annotValue(newEnt,'perturbation') <- 'PI3K'
newEnt <- addObject(newEnt, as.list(fits[[1]]), unlist=TRUE)
newEnt <- createEntity(newEnt)
newEnt <- storeEntity(newEnt)

ent <- loadEntity("syn317781")


# Get the data
dat <- exprs(ent$objects$eset)
# Perform a differential expression test.
tmp <- ent$objects$eset@protocolData@data$ScanDate
dates <- sapply(strsplit(tmp," "), function(x){ x[1]})
times <- sapply(strsplit(tmp," "), function(x){ x[2]})
perturbation <- ifelse(grepl('P11', colnames(dat)), "PI3K","GFP")
obj <- data.frame(perturbation=perturbation,
		scanDates=dates,
		scanTimes=times)
rownames(obj) <- colnames(dat)
write.table(obj,file="ectopicPi3k_metadata.txt",sep="\t",quote=FALSE)

mdEnt <- Data(list(name="ectopicPi3k_metadata", 
				parentId = propertyValue(newEnt,'parentId'), 
				platform=propertyValue(newEnt, 'platform'),
				tissueType=propertyValue(newEnt, 'tissueType'),
				disease=propertyValue(newEnt, 'disease'),
				numSamples=propertyValue(newEnt, 'numSamples'),
				species=propertyValue(newEnt, 'species')))
annotations(mdEnt) <- annotations(newEnt)
addObject(mdEnt, obj, "ectopicPi3k_metadata")
addFile(mdEnt, "ectopicPi3k_metadata.txt")
mdEnt <- createEntity(mdEnt)
mdEnt <- storeEntity(mdEnt)

mdEnt <- loadEntity('syn317785')
table(mdEnt$objects$ectopicPi3k_metadata$perturbation, mdEnt$objects$ectopicPi3k_metadata$scanDates)





X <- model.matrix(~ factor(mdEnt$objects$ectopicPi3k_metadata$perturbation))
sig <- calcSig(dat, X)
# Take an SVD of the data.  Look at eigenweights and first couple of eigengenes
u.full <- fs(dat)

png(file='pvalues.png')
hist(sig$pval, xlab="P values", ylab="Beta-Catenin Effect on Gene Expression Variation")
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

png(file="cfs.png")
plot(sig$cfs, xlab="Average in GFP",ylab="Effect of PI3K")
dev.off()

a <- paste("./add.attachments.py METAGENOMICS \"PI3K\" \"text/plain\" ", 'pvalues.png'); system(a)
#a <- paste("./add.attachments.py METAGENOMICS \"Factor Library\" \"text/plain\" ", 'cfs.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"PI3K\" \"text/plain\" ", 'u_unsup_d.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"PI3K\" \"text/plain\" ", 'eg1.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"PI3K\" \"text/plain\" ", 'eg2.png'); system(a)

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

a <- paste("./add.attachments.py METAGENOMICS \"PI3K\" \"text/plain\" ", 'pVE_rEGs.png'); system(a)

sva.fit <- sva(dat, bio.var=X, n.sv=4, num.iter=30, diagnose=FALSE)
# Now, we'll take a look at the estimated basis vectors
png(file="svaEGs.png")
plot(sva.fit$svd[[30]]$v[,1])
plot(sva.fit$svd[[30]]$v[,2])
plot(sva.fit$svd[[30]]$v[,3])
plot(sva.fit$svd[[30]]$v[,4])
dev.off()

a <- paste("./add.attachments.py METAGENOMICS \"PI3K\" \"text/plain\" ", 'svaEGs.png'); system(a)

load("~/Desktop/studies.Rda")
load("~/Documents/Randoms/gic_explore.Rda")
sign <- hsym[names(sort(-sig$cfs[,2]))[1:300]]
pi3k.sig <- calc.signature.stats("",gpl570.gic, sign)
studies[names(which(pi3k.sig[2,] > 0.5)),1:2]

png(file="ks.png")
plot(pi3k.sig[2,], xlab="Studies", ylab="Signature Enrichment")
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"SRC\" \"text/plain\" ", 'ks.png'); system(a)

rownames(dat) <- hsym[rownames(dat)]
th <- intersect(rownames(gpl570.gic), rownames(dat))
dat2 <- dat[th,]
gpl2 <- gpl570.gic[th,]

sig <- calcSig(dat2, X)
Z <- model.matrix(~ pi3k.sig[2,])
sig2 <- calcSig(gpl2, Z)

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















