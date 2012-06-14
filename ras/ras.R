library(mGenomics)
library(synapseClient)
library(snm)
synapseLogin('brig.mecham@sagebase.org','letmein')

ent <- loadEntity('syn138522')
fits <- runWorkflow(ent$cacheDir, workflow="snm")

# Get the data
dat <- exprs(fits$hthgu133a[[1]])
# Perform a differential expression test.
treatment <- ifelse(grepl('M', list.files(ent$cacheDir)), "Myc", "GFP")
X <- model.matrix(~ factor(treatment))
sig <- calcSig(dat, X)
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
	





