#==============================================================================#
# Diego Londono-Correa Using Genomic SEM to explore genetic associations among mental disorders, cognitive abilities, education, and personality
#==============================================================================#

#install.packages("devtools")
#library(devtools)
#install_github("GenomicSEM/GenomicSEM")
library(GenomicSEM)
library(data.table)

#==============================================================================#
# Step -1: Clean data
#==============================================================================#
setwd("~/Library/CloudStorage/Box-Box/SEM Project/Genomic_SEM/new_GSEM")


### LANGUAGE FACTOR EISING 2022   
nread_sum <- fread("nread_sum", data.table=FALSE)
wr_sum <- fread("wr_sum", data.table=FALSE)
sp_sum <- fread("sp_sum", data.table=FALSE)
pa_sum <- fread("pa_sum", data.table=FALSE)

#### RUNNING GENOMIC SEM   

#==============================================================================#
# Step 0: Set input variables and arguments
#==============================================================================#

#set paths to reference data.
hm3 = "w_hm3.snplist"
ld = "eur_w_ld_chr/"
wld = "eur_w_ld_chr/"
ref = "reference.1000G.maf.0.005.txt"
#set SNP filters for analysis.
info.filter = 0.90
maf.filter = 0.01

files = c("nread_sum",
          "wr_sum",
          "sp_sum",
          "pa_sum")

#set names for input traits.
trait.names = c("NW_READING",
                "W_READING",
                "SPELLING",
                "PHONEME_AWARENESS")



#set names of munged traits
traits = c("NW_READING.sumstats.gz",
           "W_READING.sumstats.gz",
           "SPELLING.sumstats.gz",
           "PHONEME_AWARENESS.sumstats.gz")



#==============================================================================#
# Step 2: Run multivariable LDSC
#==============================================================================#

#set sample size for input traits.
sample.prev <-     c(NA,NA,NA,NA)
population.prev <- c(NA,NA,NA,NA)

#To have enough buffer size for this analysis
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 6) 

##run multivariable LDSC
LDSCoutput <- ldsc(traits = traits,
                   sample.prev = sample.prev,
                   population.prev = population.prev,
                   ld = ld,
                   wld = wld,
                   trait.names = trait.names)

##optional command to save the ldsc output in case you want to use it in a later R session.

saveRDS(LDSCoutput,file="Model_SNP_effects/LANG/Language/LDSCoutput.rds")
LDSCoutput  <- readRDS("Model_SNP_effects/LANG/Language/LDSCoutput.rds")

#Save S matrix
write.csv (LDSCoutput$S, 'Model_SNP_effects/LANG/Language/LDSC_Matrix_S.csv' , row.names = F)

#Save Intercepts
write.csv (LDSCoutput$I, 'Model_SNP_effects/LANG/Language/LDSC_Matrix_I.csv' , row.names = F)

#SEs
k<-nrow(LDSCoutput$S)
SE<-matrix(0, k, k)
SE[lower.tri(SE,diag=TRUE)] <-sqrt(diag(LDSCoutput$V))
SE[upper.tri(SE)] <- t(SE)[upper.tri(SE)]
rownames(SE) <- trait.names 
colnames(SE) <- trait.names
write.csv (SE, 'Model_SNP_effects/LANG/Language/SE.csv' , row.names = T)


#==============================================================================#
# Step 3b: User-specified  Model
#==============================================================================#


#Unit Variance identification
model_uvar<-'#Unit Variance identification
#Factors
 LANG=~ NA*NW_READING + W_READING + start(0.4)*SPELLING + PHONEME_AWARENESS 
         
         #Variances: 
         LANG~~1*LANG'


output_uvar <- usermodel(LDSCoutput,estimation="DWLS",model=model_uvar)
output_uvar$modelfit
#     chisq df   p_chisq      AIC CFI      SRMR
#df 1.147588  2 0.5633838 17.14759   1 0.0141329

write.csv(output_uvar$results, 'Language/output_uvar.csv', row.names = F)

#> output_uvar$results
#               lhs op               rhs Unstand_Est         Unstand_SE STD_Genotype    STD_Genotype_SE     STD_All      p_value
#1              LANG =~        NW_READING 0.533700785 0.0298257533413965  0.996328873 0.0556796241547889 0.996328866 1.314317e-71
#4              LANG =~         W_READING 0.440120107 0.0244767346563447  0.955172304 0.0531207239972446 0.955172309 2.737008e-72
#3              LANG =~          SPELLING 0.508966087 0.0275455879307562  0.993109604 0.0537477640953504 0.993109610 3.149569e-76
#2              LANG =~ PHONEME_AWARENESS 0.480065535 0.0358690797816269  0.951833655 0.0711182045885845 0.951833739 7.517548e-41
#6        NW_READING ~~        NW_READING 0.002102912 0.0195002769510296  0.007328791 0.0679595550733353 0.007328791 9.141225e-01
#9         W_READING ~~         W_READING 0.018608468 0.0145928493117919  0.087645859 0.0687323372504766 0.087645860 2.022466e-01
#8          SPELLING ~~          SPELLING 0.003607074 0.0231288611769643  0.013733303 0.0880584314523581 0.013733303 8.760681e-01
#7 PHONEME_AWARENESS ~~ PHONEME_AWARENESS 0.023914701 0.0304889778641993  0.094012516  0.119857150836175 0.094012533 4.328219e-01
#5              LANG ~~              LANG 1.000000000                     1.000000000                    1.000000000           NA


#Unit Loading identification

#Selected nonword reading as the canonical indicator because of higher loading and h2 0.273


model_uload <-'#Unit Loading identification, Cholesky model without external traits Y vector 
#Factors
 LANG=~ 1*NW_READING + W_READING + SPELLING +   PHONEME_AWARENESS'


output_uload <- usermodel(LDSCoutput,estimation="DWLS",model=model_uload)
output_uload$modelfit
#     chisq df   p_chisq      AIC CFI      SRMR
#df 1.147589  2 0.5633837 17.14759   1 0.0141329
write.csv(output_uload$results, 'Language/output_uload.csv', row.names = F)

## Loadings on Lang

#                lhs op               rhs Unstand_Est         Unstand_SE STD_Genotype    STD_Genotype_SE     STD_All      p_value
#               LANG =~         W_READING 0.824657055 0.0539974626809894  0.958691797 0.0627738818997956 0.955172301 1.172527e-52
#3              LANG =~          SPELLING 0.953654348 0.0697747824093407  0.996768898 0.0729292883904083 0.993109621 1.585217e-42
#2              LANG =~ PHONEME_AWARENESS 0.899503137 0.0674671891946002  0.955340864 0.0716552953061139 0.951833664 1.498867e-40
#6        NW_READING ~~        NW_READING 0.002102927  0.019500276967468  0.007328813 0.0679595539774196 0.007328813 9.141219e-01
#9         W_READING ~~         W_READING 0.018608463 0.0145928496820409  0.087645875 0.0687323372665093 0.087645875 2.022467e-01
#8          SPELLING ~~          SPELLING 0.003607098 0.0231288606844847  0.013733281 0.0880584323383176 0.013733281 8.760672e-01
#7 PHONEME_AWARENESS ~~ PHONEME_AWARENESS 0.023914714  0.030488977813523  0.094012676  0.119857151408135 0.094012675 4.328216e-01
#5              LANG ~~              LANG 0.284836525 0.0318360559512781  0.992671190  0.110950432457296 0.992671190 3.653437e-19
#1              LANG =~        NW_READING 1.000000000                     1.000000000                    0.996328855           NA


#==============================================================================#
# Step 4: Multivariate GWAS with SNP associations --- Language
#==============================================================================#

## NEW LDSCS, without external traits Y vector. For GWAS, model with SNP effects
#set paths to reference data.
hm3 = "w_hm3.snplist"
ld = "eur_w_ld_chr/"
wld = "eur_w_ld_chr/"
ref = "reference.1000G.maf.0.005.txt"
#set SNP filters for analysis.
info.filter = 0.90
maf.filter = 0.01

## prepare sumstats ##


files = c("nread_sum",
          "wr_sum",
          "sp_sum",
          "pa_sum")

#set names for input traits.
trait.names = c("NW_READING",
                "W_READING",
                "SPELLING",
                "PHONEME_AWARENESS")



#set names of munged traits
traits = c("NW_READING.sumstats.gz",
           "W_READING.sumstats.gz",
           "SPELLING.sumstats.gz",
           "PHONEME_AWARENESS.sumstats.gz")

##run multivariable LDSC

#set sample size for input traits.
sample.prev <-     c(NA,NA,NA,NA)
population.prev <- c(NA,NA,NA,NA)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 6) 
LDSCoutput2 <- ldsc(traits = traits,
                    sample.prev = sample.prev,
                    population.prev = population.prev,
                    ld = ld,
                    wld = wld,
                    trait.names = trait.names,
                    stand = TRUE)

##optional command to save the ldsc output in case you want to use it in a later R session.
saveRDS(LDSCoutput2,file="Language/LDSCoutput2.rds")
LDSCoutput  <- readRDS("Language/LDSCoutput2.rds")

save(LDSCoutput2, file="Language/LDSCoutput2.RData")
LDSCoutput2 <- load("Language/LDSCoutput2.RData")

N=c(NA,NA,NA,NA)

## Prepare Summary statistics for GWAS

se.logit = c(F,F,F,F)
OLS = c(T,T,T,T)
linprob = c(F,F,F,F)


model_sumstats <- sumstats(files = files,
                           ref = ref,
                           trait.names = trait.names,
                           se.logit = se.logit,
                           OLS = OLS,
                           linprob = linprob,
                           N = N,
                           info.filter = info.filter,
                           maf.filter = maf.filter,
                           keep.indel = FALSE,
                           parallel = FALSE,
                           cores = NULL)

saveRDS(model_sumstats,file="Language/model_sumstats.rds")
model_sumstats <- readRDS("Language/model_sumstats.rds")

write.table(x = model_sumstats, file = "Language/model_sumstats.txt",
            quote = F, sep = " ",row.names = F, col.names = T)

#write in the model

model_snps <-'#Model with SNP effects
#Factors
 LANG=~ 1*NW_READING  + 0.959*W_READING + 0.997*SPELLING + 0.955*PHONEME_AWARENESS 
         
#SNP effects    
LANG ~ SNP'


#specify what parts of the model you want to save
sub<-c('LANG ~ SNP')

#example. Real was processed on TACC
GWAS_Run <-userGWAS(covstruc=LDSCoutput2,SNPs=model_sumstats,model=model_snps,sub=sub,  smooth_check = T)


