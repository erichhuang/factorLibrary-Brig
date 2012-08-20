library(mGenomics)
library(synapseClient)
library(snm)

ent <- loadEntity('syn138507')
fits <- runWorkflow(ent$cacheDir, workflow="snm")

# Get the data
dat <- exprs(fits$hgu133plus2[[1]])
# Perform a differential expression test.
treatment <- ifelse(grepl('E2F3', list.files(ent$cacheDir)), "GFP", "E2F3")
X <- model.matrix(~ factor(treatment))
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
par(mfrow=c(2,2))
plot(sva.fit$svd[[30]]$v[,1])
plot(sva.fit$svd[[30]]$v[,2])
plot(sva.fit$svd[[30]]$v[,3])
plot(sva.fit$svd[[30]]$v[,4])
dev.off()
a <- paste("./add.attachments.py METAGENOMICS \"Beta-Catenin\" \"text/plain\" ", 'svaEGs.png'); system(a)
