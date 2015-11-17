library(VGAM)

lambda <- 0.5;

marginalLikelihood <- function(params, data) {
		
	alpha <- params[1];
	beta <- params[2];
	vec <-data;
	dlen <- length(vec);
	
	As <- vec[seq(1, dlen, 2)];
	Bs <- vec[seq(2, dlen, 2)];
	
	# print(As);
	# print(Bs);
	
	ML <- 0;
	for (i in 1:length(As)) {
		ML <- ML + ( lgamma(As[i] + Bs[i] + 1) - lgamma(As[i] + 1) - lgamma(Bs[i] + 1));
		ML <- ML - ( lgamma(alpha + beta + As[i] + Bs[i]) - lgamma(alpha + As[i]) - lgamma(beta + Bs[i]) );
		ML <- ML + ( lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) );
	}
	
	ML <- ML - lambda * log(alpha + beta);

	return(-ML);	
	
}


marginalLikelihood_PG <- function(params, data) {
                
        r <- params[1];
        p <- params[2];
     
        vec <-data;
        dlen <- length(vec);
        
        ks <- vec[seq(1, dlen, 2)];
        weight <- vec[seq(2, dlen, 2)];
        
        ML <- 0;
        for (i in 1:length(ks)) {
                ML <- ML + ks[i] * log(1 - p) + r * log(p) / weight[i];
                ML <- ML + lgamma(ks[i] + r / weight[i]) - lgamma(r / weight[i]) - lgamma(ks[i] + 1);
        }
        return(-ML);    
        
}



my_trans <- function(x = 0) {
    return( max(0, -log10(x)));
}


bdata <- read.table(commandArgs()[5], sep="\t",header=F);

tmapB_L <- read.table(commandArgs()[6], sep="\t",header=F);
mapB_L <- as.numeric(tmapB_L[,2] / 100000000);

tmapB_t <- read.table(commandArgs()[7], sep="\t",header=F);
mapB_t <- as.numeric(tmapB_t[1]) / 100000000;

print(mapB_L);
print(mapB_t);


infoData <- bdata[,1:6];
tumData <- bdata[,7];
norData <- bdata[,8];
refData <- bdata[,9:ncol(bdata)];

positiveList <- rep(0, nrow(bdata));
resultList <- matrix(0, nrow(bdata), 13);

for (i in 1:nrow(infoData)) {

    print(infoData[i,]);
	
	tumBases <- as.integer(unlist(strsplit(as.character(tumData[i]), ",")));
	norBases <- as.integer(unlist(strsplit(as.character(norData[i]), ",")));
	
	if (pbinom(norBases[2], norBases[1], 0.5) < 1) {
		
		derror <- rep(0, 2 * ncol(refData));
		
		for (j in 1:ncol(refData)) {
			
			refBases <- as.integer(unlist(strsplit(as.character(refData[i, j]), ",")));
			
			derror[2 * j - 1] <- refBases[2];
			derror[2 * j] <- refBases[1] - refBases[2];
			
		}

        ref_tum <- tumBases[1] - tumBases[2];
        alt_tum <- tumBases[2];
        ref_nor <- norBases[1] - norBases[2];
        alt_nor <- norBases[2];

        Fmatrix <- matrix(as.integer(c(alt_tum, alt_nor, ref_tum, ref_nor)), 2, 2);
        fisPV <- max(0, -log10(fisher.test(Fmatrix)$p.value));

        print(fisPV);
        if (fisPV < 1) {
            next;
        }       
   


        # print(derror_p);
        # print(derror_n);

        alpha <- 0.1;
        beta <- 1;
        if (sum(derror) > 0) {		
		    res <- constrOptim(c(20, 20), marginalLikelihood, grad=NULL, ui=matrix(c(1, 0, 1, 0, 1, 1), 3, 2), ci=c(0.1, 1, 1), data=derror);
		    alpha <- res$par[1];
            beta <- res$par[2];
        }

		lpv <- 0;
		lpv <- 0;
		
        if (tumBases[2] > 0) {
            pv <- pbetabinom.ab(tumBases[2] - 1, tumBases[1], alpha, beta);
            if (1 - pv < 1e-60) {
                lpv <- 60;
            } else {
                lpv <- -log10(1 - pv);
            }
        }

		misRatio_tum <- (tumBases[2]) / (tumBases[1]);
		misRatio_nor <- (norBases[2]) / (norBases[1]);

        print(lpv);
		if (lpv > -log10(0.01)) {
			positiveList[i] <- 1
		} else {
            next;
        }

        derror_PG <- derror;
        derror_PG[seq(2, length(derror_PG), 2)] <- mapB_L;

        r_PG <- 10;
        p_PG <- 0.001;
        if (sum(derror) > 0) {
            res <- constrOptim(c(10, 0.1), marginalLikelihood_PG, grad=NULL, ui=matrix(c(1, 0, 0, 0, 1, -1), 3, 2), ci=c(1, 0.01, -0.99), data=derror_PG);
            r_PG <- res$par[1];
            p_PG <- res$par[2];
        }

        
        GP_PV <- max(0, -log10(1 - pnbinom(tumBases[2], r_PG / mapB_t, p_PG)));
        print(GP_PV);
        if (GP_PV == Inf) {
            GP_PV <- 60;
        }
        # print(c(tumBases[2], mapB_t, derror_PG));
        # print(c(GP_PV));
        # print(c(r_PG, p_PG));
	
        resultList[i,] <- c(misRatio_tum, ref_tum + alt_tum, alt_tum, misRatio_nor, ref_nor + alt_nor, alt_nor, lpv, fisPV, GP_PV, alpha, beta, r_PG, p_PG);
        print(resultList[i,]);
	}
	
}



write.table(cbind(infoData[positiveList == 1,,drop=FALSE], tumData[positiveList == 1,drop=FALSE], norData[positiveList == 1,drop=FALSE], resultList[positiveList == 1,,drop=FALSE]), commandArgs()[8], sep="\t", row.names=F, col.names=F);






