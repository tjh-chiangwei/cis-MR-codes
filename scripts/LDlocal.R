library(TwoSampleMR)

#exposure随机id
random_string <- function(n=1, len=6)
{
  randomString <- c(1:n)
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    len, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}

#本地ld clump
clump_data2 <- function(dat, clump_kb=250, clump_r2=0.01, clump_p1=1)
{
  # .Deprecated("ieugwasr::ld_clump()")
  
  pval_column <- "pval.exposure"
  
  if(!is.data.frame(dat))
  {
    stop("Expecting data frame returned from format_data")
  }
  
  if("pval.exposure" %in% names(dat) & "pval.outcome" %in% names(dat))
  {
    message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
  } else if(!"pval.exposure" %in% names(dat) & "pval.outcome" %in% names(dat))
  {
    message("pval.exposure column not present, using pval.outcome column for clumping.")
    pval_column <- "pval.outcome"
  } else if(! "pval.exposure" %in% names(dat))
  {
    message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
    dat$pval.exposure <- 0.99
  } else {
    pval_column <- "pval.exposure"
  }
  
  if(! "id.exposure" %in% names(dat))
  {
    dat$id.exposure <- random_string(1)
  }
  
  d <- data.frame(rsid=dat$SNP, pval=dat[[pval_column]], id=dat$id.exposure)
  out <- ieugwasr::ld_clump(d, 
                            clump_kb=clump_kb, 
                            clump_r2=clump_r2, 
                            clump_p=clump_p1,
                            plink_bin = "C:/Users/chian/plink/plink.exe",
                            bfile = "C:/Users/chian/plink/EUR"
                            )
  
  keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, out$id)
  return(dat[keep, ])
}


#本地获取matrix
ld_matrix2 <- function(snps)
{
  # .Deprecated("ieugwasr::ld_matrix()")
  ieugwasr::ld_matrix(variants=snps,
                      with_alleles = TRUE,
                      pop = "EUR",
                      plink_bin = "C:/Users/chian/plink/plink.exe",
                      bfile = "C:/Users/chian/plink/EUR"
  )
}

#转换数据，源代码有bug。
dat_to_MRInput2 <- function(dat, get_correlations=FALSE)
{
  out <- plyr::dlply(dat, c("exposure", "outcome"), function(x)
  {
    x <- plyr::mutate(x)
    message("Converting:")
    message(" - exposure: ", x$exposure[1])
    message(" - outcome: ", x$outcome[1])
    if(get_correlations)
    {
      message(" - obtaining LD matrix")
      ld <- ld_matrix2(unique(x$SNP))
      
      #排序，bug所在。
      rns=unlist(data.frame(strsplit(rownames(ld), split="_"))[1,])
      ord <- match(rns,x$SNP)
      ld <- ld[order(ord), order(ord)]
      rownames(ld)=gsub("TRUE","T",rownames(ld))
      colnames(ld)=gsub("TRUE","T",colnames(ld))
      
      out <- harmonise_ld_dat(x, ld) 
      if(is.null(out))
      {
        return(NULL)
      }
      x <- out$x
      ld <- out$ld
      
      MendelianRandomization::mr_input(
        bx = x$beta.exposure,
        bxse = x$se.exposure,
        by = x$beta.outcome,
        byse = x$se.outcome,
        exposure = x$exposure[1],
        outcome = x$outcome[1],
        snps = x$SNP,
        effect_allele=x$effect_allele.exposure,
        other_allele=x$other_allele.exposure,
        eaf = x$eaf.exposure,
        correlation = ld
      )
      
    } else {
      MendelianRandomization::mr_input(
        bx = x$beta.exposure,
        bxse = x$se.exposure,
        by = x$beta.outcome,
        byse = x$se.outcome,
        exposure = x$exposure[1],
        outcome = x$outcome[1],
        snps = x$SNP,
        effect_allele=x$effect_allele.exposure,
        other_allele=x$other_allele.exposure,
        eaf = x$eaf.exposure
      )
    }
  })
  return(out)
}

ToMR2 <- function(mydata){
  if (nrow(mydata) == 1) {
    mydata2 = dat_to_MRInput2(mydata)
  } else
    mydata2 = dat_to_MRInput2(mydata,get_correlation=TRUE) 
} #转化为MR文件

