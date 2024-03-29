# the messy one with all the output and print statements. For messing about with.
suppressMessages(library(glmnet))
suppressMessages(library(methods))
"%&%" <- function(a,b) paste(a, b, sep = "")

TW_CV_model <- function(expression_RDS, geno_file, gene_annot_RDS, snp_annot_RDS, n_k_folds, alpha, out_dir, tis, chrom, snpset, window, seed = NA) {
  expression <- readRDS(expression_RDS)
  allout <- tis %&% "_genesall_" %&% chrom %&% ".txt"
  # for (i in 1:ncol(expression)) {
  #   gene = colnames(expression)[i]
  #   write(gene,file=allout,append=TRUE) #21024
  # }
  class(expression) <- 'numeric'
  genotype <- read.table(geno_file, header = TRUE, row.names = 1, stringsAsFactors = FALSE,check.names=FALSE)
  # Transpose genotype for glmnet
  genotype <- t(genotype)
  gene_annot <- readRDS(gene_annot_RDS)
  gene_annot <- subset(gene_annot, gene_annot$chr == chrom)
  snp_annot <- readRDS(snp_annot_RDS)
  print("snp_annot:")
  print(head(snp_annot))

  rownames(gene_annot) <- gene_annot$gene_id
  # Subset expression data to only include genes with gene_info
  expression <- expression[, intersect(colnames(expression), rownames(gene_annot))]
  # subset expression to only those samples with genotypes
  expression <- expression[rownames(expression) %in% rownames(genotype),]

  interout <- tis %&% "_genesinter_" %&% chrom %&% ".txt"
  # for (i in 1:ncol(expression)) {
  #   gene = colnames(expression)[i]
  #   write(gene,file=interout,append=TRUE) #15978
  # }
  exp_samples <- rownames(expression)
  exp_genes <- colnames(expression)
  n_samples <- length(exp_samples)
  n_genes <- length(exp_genes)
  seed <- ifelse(is.na(seed), sample(1:2016, 1), seed)
  log_df <- data.frame(chrom, n_genes, seed, alpha)
  colnames(log_df) <- c('chr', 'n_genes', 'seed_for_cv', 'alpha')
  write.table(log_df, file = out_dir %&% tis %&% '_chr' %&% chrom %&% '_elasticNet_model_log.txt', quote = FALSE, row.names = FALSE, sep = "\t")
  set.seed(seed)
  groupid <- sample(1:n_k_folds, length(exp_samples), replace = TRUE)
  print(groupid)
  print(length(groupid))

  resultsarray <- array(0, c(length(exp_genes), 9))
  dimnames(resultsarray)[[1]] <- exp_genes
  resultscol <- c("gene", "alpha", "cvm", "lambda.iteration", "lambda.min", "n.snps", "R2", "pval", "genename")
  dimnames(resultsarray)[[2]] <- resultscol
  workingbest <- out_dir %&% "working_TW_" %&% tis %&% "_exp_" %&% n_k_folds %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_chr" %&% chrom %&% ".txt"
  write(resultscol, file = workingbest, ncolumns = 9, sep = "\t")

  weightcol <- c("gene","rsid","ref","alt","beta","alpha")
  workingweight <- out_dir %&% "TW_" %&% tis %&% "_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% "_weights_chr" %&% chrom %&% ".txt"
  write(weightcol, file = workingweight, ncol = 6, sep = "\t")
  expgenesout <- tis %&% "_genesexpgenes_" %&% chrom %&% ".txt"
  cisgenosout <- tis %&% "_genescisgenos_" %&% chrom %&% ".txt"
  covariance_out <- out_dir %&% tis %&% '_chr' %&% chrom %&% '_snpset_' %&% snpset %&% '_alpha_' %&% alpha %&% "_covariances.txt"

  for (i in 1:length(exp_genes)) {
    cat(i, "/", length(exp_genes), "\n")
    gene <- exp_genes[i]
    # Reduce genotype data to only include SNPs within specified window of gene.
    geneinfo <- gene_annot[gene,]
    start <- geneinfo$start - window
    end <- geneinfo$end + window
    # Pull cis-SNP info
    cissnps <- subset(snp_annot, snp_annot$pos >= start & snp_annot$pos <= end)
    print(head(cissnps))

    # Pull cis-SNP genotypes
    # cisgenos <- data.frame(genotype[,intersect(colnames(genotype), cissnps$varID), drop = FALSE]) #join by varID
    cisgenos <- data.frame(genotype[,intersect(colnames(genotype), cissnps$rsid), drop = FALSE]) #join by rsID
      # print(nrow(cisgenos))
      # print(ncol(cisgenos))
    # Reduce cisgenos to only include SNPs with at least 1 minor allele in dataset
    cm <- colMeans(cisgenos, na.rm = TRUE)
    minorsnps <- subset(colMeans(cisgenos), cm > 0 & cm < 2)
    minorsnps <- names(minorsnps)
    cisgenos <- data.frame(cisgenos[,minorsnps, drop = FALSE])
      # print(typeof(cisgenos))
    missingout <- tis %&% "_genesmissing_" %&% chrom %&% ".txt"
    if (ncol(cisgenos) < 2) { #only do this if there are _no_ SNPs nearby.
      # Need 2 or more cis-snps for glmnet.
      write(gene,file=missingout,append=TRUE) # 264 here
      bestbetas <- data.frame()
    } else {
      # Pull expression data for gene
      exppheno <- expression[,gene]
      # Scale for fastLmPure to work properly
      exppheno <- scale(exppheno, center = TRUE, scale = TRUE)
      exppheno[is.na(exppheno)] <- 0
      rownames(exppheno) <- rownames(expression)
      # reduce cisgenos to only include individuals we have exp data for
      #      print(head(cisgenos))
      #      print(head(exppheno))
      cisgenos <- cisgenos[rownames(cisgenos) %in% rownames(exppheno),]
      # print(nrow(cisgenos))
      # print(ncol(cisgenos))
      # print(typeof(cisgenos))
      # print("before (genes/snps)")
      # print(nrow(exppheno))
      # print(ncol(exppheno))
      # print(typeof(exppheno))

      #-----debug-----
      #reduce exppheno to only include individuals we have geno data
      # print(head(exppheno[,1]))
      exppheno <- as.matrix(exppheno[rownames(exppheno) %in% rownames(cisgenos),])
      #match the sampleid
      cisgenos<-cisgenos[sapply(rownames(exppheno), function(x) which(x==rownames(cisgenos))),]
      #-----debug-----
      print("after (genes/snps)")
      print(nrow(exppheno))
      # print(sapply(cisgenos,typeof))
      print(nrow(cisgenos))

      # Run Cross-Validation
      # parallel = TRUE is slower on tarbell
      write(gene,file=cisgenosout,append=TRUE) #
      bestbetas <- tryCatch(
        {
        #set seed-debug
        set.seed(as.numeric(sub('^....','',gene)))
        print(ncol(cisgenos))
        print(dim(cisgenos))
        print(dim(exppheno))
        print(dim(as.matrix(cisgenos)))
        print(groupid)
        fit <- cv.glmnet(as.matrix(cisgenos),
                           as.vector(exppheno),
                           nfolds = n_k_folds,
                           alpha = alpha,
                           keep = TRUE,
                           foldid = groupid,
                           parallel = FALSE)
        # Pull info from fit to find the best lambda
        fit.df <- data.frame(fit$cvm, fit$lambda, 1:length(fit$cvm))
        # Needs to be min or max depending on cv measure (MSE min, AUC max, ...)
        best.lam <- fit.df[which.min(fit.df[,1]),]
        cvm.best <- best.lam[,1]
        lambda.best <- best.lam[,2]
        # Position of best lambda in cv.glmnet output
        nrow.best <- best.lam[,3]
        # Get the betas from the best lambda value
        ret <- as.data.frame(fit$glmnet.fit$beta[,nrow.best])
        print("ret")
        print(head(ret))
        for (row in 1:length(ret[ret == 0.0])) {write(gene,file=expgenesout,append=TRUE)}
        ret[ret == 0.0] <- NA
        # Pull the non-zero betas from model
        as.vector(ret[which(!is.na(ret)),])
        },
        error = function(cond) { #not getting fired
          # Should fire only when all predictors have 0 variance.
          message('Error with gene ' %&% gene %&% ', index ' %&% i)
          message(geterrmessage())
          return(data.frame())
        }
      )
    }
    print(bestbetas)
    if (length(bestbetas) > 0) {
      names(bestbetas) <- rownames(ret)[which(!is.na(ret))]
      # Pull out the predictions at the best lambda value.
      pred.mat <- fit$fit.preval[,nrow.best]
      res <- summary(lm(exppheno~pred.mat))
      genename <- as.character(gene_annot[gene, 3])
      rsq <- res$r.squared
      pval <- res$coef[2,4]
      resultsarray[gene,] <- c(gene, alpha, cvm.best, nrow.best, lambda.best, length(bestbetas), rsq, pval, genename)
      # Output best shrunken betas for PrediXcan
      bestbetalist <- names(bestbetas)
      print("bestbetalist:")
      print(bestbetalist)
      print(head(snp_annot))
      bestbetainfo <- snp_annot[bestbetalist,]
      betatable <- as.matrix(cbind(bestbetainfo,bestbetas))
      print(betatable)
      write_covariance(gene, cisgenos, betatable[,"rsid"], betatable[,"varID"], covariance_out)
      # Output "gene", "rsid", "refAllele", "effectAllele", "beta"
      # For future: To change rsid to the chr_pos_ref_alt_build label, change "rsid" below to "varID".
      betafile <- cbind(gene,betatable[,"rsid"],betatable[,"refAllele"],betatable[,"effectAllele"],betatable[,"bestbetas"], alpha)
      # Transposing betafile necessary for correct output from write() function
      write(t(betafile), file = workingweight, ncolumns = 6, append = TRUE, sep = "\t")
      write(resultsarray[gene,], file = workingbest, ncolumns = 9, append = TRUE, sep = "\t")
    } else {
      genename <- as.character(gene_annot[gene,3])
      resultsarray[gene,1] <- gene
      resultsarray[gene,2:8] <- c(alpha,NA,NA,NA,0,NA,NA)
      resultsarray[gene,9] <- genename
    }
  }
  write.table(resultsarray,file=out_dir %&% "TW_" %&% tis %&% "_chr" %&% chrom %&% "_exp_" %&% n_k_folds %&% "-foldCV_elasticNet_alpha" %&% alpha %&% "_" %&% snpset %&% ".txt",quote=F,row.names=F,sep="\t")
}

write_covariance <- function(gene, cisgenos, model_rsids, model_varIDs, covariance_out) {
  print(model_varIDs)
  print("rsids")
  print(model_rsids)
  # model_geno <- cisgenos[,model_varIDs, drop=FALSE]
  model_geno <- cisgenos[,model_rsids, drop=FALSE]
  geno_cov <- cov(model_geno)
  cov_df <- data.frame(gene=character(),rsid1=character(),rsid2=character(), covariance=double())
  for (i in 1:length(model_rsids)) {
    for (j in i:length(model_rsids)) {
      cov_df <- rbind(cov_df, data.frame(gene=gene,rsid1=model_rsids[i], rsid2=model_rsids[j], covariance=geno_cov[i,j]))
    }
  }
  write.table(cov_df, file = covariance_out, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}
