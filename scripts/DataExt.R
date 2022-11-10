library(TwoSampleMR) 

exp <-extract_instruments(
  outcomes='prot-*', 
  p1=1e-5,
  clump=F,
  access_token = NULL
) 
#OR local summary data
exp <- format_data(
  exp,
  type = 'exposure',#±©Â¶
  snps = exp$variant_id,
  head = Ture,
  phenotype_col = "phenotype",
  snp_col = "variant_id",#SNP
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col ="effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  ncase_col = "N_CASES",
  ncontrol_col = "N_CONTROLS",
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)

mygene=GWAS_ID[*,] #dataset for gene locations and GWAS_IDs

Exposure=subset(exp,exp$chr.exposure == mygene$Chro[1] &
                  exp$pos.exposure>=(mygene$Start[1]-200000)&
                  exp$pos.exposure<=(mygene$End[1]+200000)
) #cis-SNPs

Exposure=clump_data(Exposure,clump_r2=0.4,clump_kb=250) 
Exposure$R2=get_r_from_pn(Exposure$pval.exposure, Exposure$samplesize.exposure)
Exposure$'F' <- (N-k-1)/k*Exposure$R2/(1-Exposure$R2) #N and k first
Exposure=Exposure[Exposure$F>=10,] 

Outcome <- format_data(
  gData,
  type = 'outcome', 
  snps = Exposure$SNP,
  head = Ture,
  phenotype_col = "phenotype", 
  snp_col = "variant_id",
  beta_col = "beta", 
  se_col = "standard_error", 
  effect_allele_col ="effect_allele", 
  other_allele_col = "other_allele", 
  eaf_col = "effect_allele_frequency", 
  pval_col = "p_value", 
  ncase_col = "N_CASES", 
  ncontrol_col = "N_CONTROLS", 
  chr_col = "chromosome", 
  pos_col = "base_pair_location" 
) 
head(Outcome) 

mydata <- harmonise_data(
  exposure_dat=Exposure,
  outcome_dat=Outcome,
  action= 2
) 
