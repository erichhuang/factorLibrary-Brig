library(mGenomics)
library(synapseClient)
library(snm)
synapseLogin('brig.mecham@sagebase.org','letmein')

ent <- loadEntity('syn138511')
fits <- runWorkflow(ent$cacheDir, workflow="snm")

newEnt <- Data(list(name="ectopicEGFR", 
				parentId = "syn275939", 
				platform=names(fits),						
				tissueType="breast",
				numSamples=ncol(exprs(fits[[1]][[1]])),
				species="Homo sapiens"))
annotValue(newEnt,'processingAlgorithm') <- 'snmUnsupervised'
annotValue(newEnt,'cellLine') <- 'MCF'
annotValue(newEnt,'perturbation') <- 'EGFR'
newEnt <- addObject(newEnt, as.list(fits[[1]]), unlist=TRUE)
newEnt <- createEntity(newEnt)
newEnt <- storeEntity(newEnt)

ent <- loadEntity("syn317799")

# Get the data
dat <- exprs(ent$objects$eset)
# Perform a differential expression test.
tmp <- ent$objects$eset@protocolData@data$ScanDate
dates <- sapply(strsplit(tmp," "), function(x){ x[1]})
times <- sapply(strsplit(tmp," "), function(x){ x[2]})
perturbation <- ifelse(grepl('GFP', colnames(dat)), "GFP","EGFR")
obj <- data.frame(perturbation=perturbation,
		scanDates=dates,
		scanTimes=times)
rownames(obj) <- colnames(dat)
write.table(obj,file="ectopicEGFR_metadata.txt",sep="\t",quote=FALSE)

mdEnt <- Data(list(name="ectopicEGFR_metadata", 
				parentId = propertyValue(newEnt,'parentId'), 
				platform=propertyValue(newEnt, 'platform'),
				tissueType=propertyValue(newEnt, 'tissueType'),
				disease=propertyValue(newEnt, 'disease'),
				numSamples=propertyValue(newEnt, 'numSamples'),
				species=propertyValue(newEnt, 'species')))
annotations(mdEnt) <- annotations(newEnt)
addObject(mdEnt, obj, "ectopicEGFR_metadata")
addFile(mdEnt, "ectopicEGFR_metadata.txt")
mdEnt <- createEntity(mdEnt)
mdEnt <- storeEntity(mdEnt)

mdEnt <- loadEntity('syn317801')
table(mdEnt$objects$ectopicEGFR_metadata$perturbation, mdEnt$objects$ectopicEGFR_metadata$scanDates)

X <- model.matrix(~ factor(perturbation))
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

a <- paste("./add.attachments.py METAGENOMICS \"Beta-Catenin\" \"text/plain\" ", 'pvalues.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"Beta-Catenin\" \"text/plain\" ", 'u_unsup_d.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"Beta-Catenin\" \"text/plain\" ", 'eg1.png'); system(a)
a <- paste("./add.attachments.py METAGENOMICS \"Beta-Catenin\" \"text/plain\" ", 'eg2.png'); system(a)

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

a <- paste("./add.attachments.py METAGENOMICS \"Beta-Catenin\" \"text/plain\" ", 'pVE_rEGs.png'); system(a)

sva.fit <- sva(dat, bio.var=X, n.sv=4, num.iter=30, diagnose=FALSE)
# Now, we'll take a look at the estimated basis vectors
png(file="svaEGs.png")
plot(sva.fit$svd[[30]]$v[,1])
plot(sva.fit$svd[[30]]$v[,2])
plot(sva.fit$svd[[30]]$v[,3])
plot(sva.fit$svd[[30]]$v[,4])
dev.off()

a <- paste("./add.attachments.py METAGENOMICS \"Beta-Catenin\" \"text/plain\" ", 'svaEGs.png'); system(a)
