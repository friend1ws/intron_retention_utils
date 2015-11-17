table <- read.table(commandArgs()[5], sep="\t",header=F);

pvalues <- rep(0, length(table[,1]));

# uniq tumor:9, normal:25
for (i in 1:length(table[,1])) {

	# extract two freqent alleles
    tum_alt <- table[i, 9];
    nor_alt <- table[i, 10];
    tum_ref <- table[i, 7] - table[i, 9];
    nor_ref <- table[i, 8] - table[i, 10];

    # print(as.integer(c(tumor[ind1], normal[ind1], tumor[ind2], normal[ind2])));

	Fmatrix <- matrix(as.integer(c(tum_alt, nor_alt, tum_ref, nor_ref)), 2, 2);
	pvalues[i] <- fisher.test(Fmatrix)$p.value;

}


my_trans <- function(x = 0) {
	return( max(0, -log10(x)));
}


kekka <- apply(as.matrix(pvalues), c(1,2), my_trans);
write.table(cbind(table[kekka > -log10(0.1),], kekka[kekka > -log10(0.1)]), commandArgs()[6], sep="\t", row.names=F, col.names=F);






