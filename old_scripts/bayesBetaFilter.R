table <- read.table(commandArgs()[5], sep="\t",header=F);

betas <- matrix(0, length(table[,1]), 3);

refb <- table[, 7] - table[, 8];
misb <- table[, 8];

for (i in 1:length(table[,1])) {

    # extract two freqent alleles

    betas[i, 1] <- qbeta(0.1, misb[i] + 1, refb[i] + 1);
    betas[i, 2] <- (misb[i] + 1) / (refb[i] + misb[i] + 2);
    betas[i, 3] <- qbeta(0.9, misb[i] + 1, refb[i] + 1);

}


#write.table(cbind( table[betas[,1] > 0.05,], betas[betas[,1] > 0.05,]), commandArgs()[6], sep="\t", row.names=F, col.names=F);
# write.table(cbind( table[misb > 1 & betas[,1] > 0.05,,drop=FALSE], betas[misb > 1 & betas[,1] > 0.05,,drop=FALSE]), commandArgs()[6], sep="\t", row.names=F, col.names=F);
write.table(cbind( table[,,drop=FALSE], betas[,,drop=FALSE]), commandArgs()[6], sep="\t", row.names=F, col.names=F);





