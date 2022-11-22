# LDPRED2 MEDIANTE SUMSTATS Y DATA TARGET PROPIOS, LIMITADO A SNPS DE HAPMAP3
# GUILLERMO SÁEZ GARCIA - JUNIO 2022

library(remotes)
library(data.table) 
remotes::install_github("https://github.com/privefl/bigsnpr.git")

library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

setwd("/media/jjanssen/DATA/projects/genetics/Metodo LDpred2/Outcome")

# Los autores de LDpred2 recomiendan restringir el análisis a los SNPs de HapMap3
info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"))


#### A) CARGA Y FILTRADO DE DATOS

# summary statistics (SCZ)
setwd("/media/jjanssen/DATA/projects/genetics/ARCHIVOS NECESARIOS")
dataDiscovery <- read.table('PGC3_leave_out_cibersam', header=TRUE, stringsAsFactors = FALSE)  ## 7583883 variantes

# Cohorte CIBERSAM (target sample)
setwd("./PGC_imput_QC_filtered")
dataTarget <- read.table('PGC_imput_QC_filtered.bim', header = FALSE, stringsAsFactors = FALSE) ## 9340561 variantes


### FILTRAMOS DISCOVERY ###

setwd("/media/jjanssen/DATA/projects/genetics/Metodo LDpred2/Outcome")
head(dataDiscovery)

# Filtro de variantes imputadas a mala calidad
temp_discovery1<-subset(dataDiscovery, INFO >= 0.8)  ## 6938323 variantes

# Eliminamos cvariantes con strand ambiguity (las TA y CG). En este paso eliminamos tb indels
temp_discovery1$POS<-paste(temp_discovery1$CHR, temp_discovery1$BP, sep=":")
temp_discovery1$C<-paste(temp_discovery1$A1, temp_discovery1$A2, sep="")
temp_discovery2<-subset(temp_discovery1, !(C %in% c("AT","TA", "CG", "GC"))) ## 5884966 variantes

# Creo un archivo con frecuencia de cada variante por si aparece repetida, y eliminamos repetidas
FRECUENCIAS_discovery<-as.data.frame(table(temp_discovery2$POS))
temp_discovery2_UNIQ<-temp_discovery2[temp_discovery2$POS %in% FRECUENCIAS_discovery$Var1[FRECUENCIAS_discovery$Freq < 2], ]
temp_discovery2_UNIQ$POS<-sub("chr", "", temp_discovery2_UNIQ$POS)
dataDiscovery <- temp_discovery2_UNIQ
head(dataDiscovery)              ## 5883000 variantes

# Filtro columnas que no necesito. Me quedo con la informacion relevante, renombro columnas y cambiamos a BETA
temp_discovery3<-dataDiscovery[,c("SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P","Nca","Nco","POS")]
temp_discovery3$OR = log(temp_discovery3$OR)
temp_discovery3$N = rowSums(temp_discovery3[,c(9,10)])
temp_discovery3<-temp_discovery3[,c("SNP", "CHR", "BP", "A1", "A2", "OR", "SE","N","P","POS")]
colnames(temp_discovery3)<-c("rsid", "chr", "pos", "a0", "a1", "beta","beta_se","N","p","POS")
dataDiscovery <- temp_discovery3
rm(list=ls(pattern="^temp"))  ## Elimino objetos temporales


### FILTRAMOS TARGET ###

head(dataTarget)

# Cambio nombre de columnas
colnames(dataTarget) <- c("V1", "varID", "V3", "pos", "A1", "A2")

# Creo un archivo con frecuencia de cada variante por si aparece repetida, y eliminamos repetidas
FRECUENCIAS_target<-as.data.frame(table(dataTarget$varID))
temp_target_UNIQ<-dataTarget[dataTarget$varID %in% FRECUENCIAS_target$Var1[FRECUENCIAS_target$Freq < 2], ]
dataTarget <- temp_target_UNIQ  ## 9340561 variantes
rm(list=ls(pattern="^temp"))         ## Elimino objetis temporales
rm(list=ls(pattern="^FRECUENCIAS"))  ## Elimino objetis temporales



### HACEMOS LA INTERSECCI?N ENTRE LOS DOS ARCHIVOS PARA QUEDARNOS CON LAS VARIANTES PRESENTES EN AMBOS ARCHIVOS

# Me quedo con las variantes de un archivo presentes en el otro, y les cambio el nombre para distinguirlos
dataTarget_filtered <- dataTarget[which(dataTarget$varID %in% dataDiscovery$POS),]                 ## 5692221 variantes
dataDiscovery_filtered <- dataDiscovery[which(dataDiscovery$POS %in% dataTarget_filtered$varID),]  ## 5692221 variantes


# Exporto el archivo
setwd("/media/jjanssen/DATA/projects/genetics/Metodo LDpred2/Outcome")
dataDiscovery_filtered<-dataDiscovery_filtered[,c("rsid", "chr", "pos", "a0", "a1", "beta","beta_se","N","p")]
write.table(dataDiscovery_filtered, "dataDiscovery_filtered.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Exportamos las variantes filtradas en la target
varIDs <- data.frame(dataTarget_filtered$varID)
write.table(varIDs, 'variableIDs_forFiltering.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)


# Creamos nuevos archivos de la target sample en formato PLINK .bed .bim and .fam solo con las variantes antes filtradas
# Eliminamos tambi?n las variantes dentro de la regi?n del MHC (chr6:26Mb - 33Mb)
setwd("/media/jjanssen/DATA/projects/genetics/ARCHIVOS NECESARIOS")
system2("./plink", args=c("--bfile ./PGC_imput_QC_filtered/PGC_imput_QC_filtered", '--extract "../Metodo LDpred2/Outcome/variableIDs_forFiltering.txt"', "--exclude range ./MHC_filter.txt", "--make-bed", '--out "../Metodo LDpred2/Outcome/PGC_imput_QC_new_filtered"'))

head(dataDiscovery_filtered)
colnames(dataDiscovery_filtered) <- c("rsid","chr","pos","a0","a1","beta","beta_se","n_eff","p")
sumstats <- dataDiscovery_filtered

# Nos quedamos con los SNPs de HapMap3
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]

setwd("/media/jjanssen/DATA/projects/genetics/Metodo LDpred2/Outcome")
write.table(sumstats, "sumstats.txt", sep="\t", quote=FALSE, row.names=FALSE)



#### B) CALCULAMOS LA LD MATRIX

setwd("/media/jjanssen/DATA/projects/genetics/Metodo LDpred2/Outcome")
sumstats  <- read.table('sumstats.txt', header=TRUE, stringsAsFactors = FALSE) 

# Conseguimos el número máximo de cores
NCORES <- nb_cores()
# Abrimos un archivo temporal
tmp <- tempfile(tmpdir = "tmp-data")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Inicializamos variables para almacenar el LD score y LD matrix
corr <- NULL
ld <- NULL
# Inicializamos otra variable para tener el orden de las muestras en el archivo bed
fam.order <- NULL
# Preprocesamos el archivo .bed
snp_readBed("/media/jjanssen/DATA/projects/genetics/Metodo LDpred2/Outcome/PGC_imput_QC_new_filtered.bed")
# now attach the genotype object
obj.bigSNP <- snp_attach("PGC_imput_QC_new_filtered.rds")

# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a0", "a1")
# perform SNP matching
info_snp <- snp_match(sumstats, map)
# Asignamos genotype a una variable 
genotype <- obj.bigSNP$genotypes
# Renombramos la estructura de map
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")

### CALCULAMOS LD
for (chr in 1:22) {
  # Extraemos los SNPs que están incluidos en cada chr
  print(chr)
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$'_NUM_ID_'[ind.chr]
  
  # Calculamos una puntuación de LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = NCORES,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}


# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))


### Perform LD score regression

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]



# Ajustamos df_beta utilizando el modo automático (hay varios modos, utilizamos auto)
multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
  ncores = NCORES
)
beta_auto <- sapply(multi_auto, function(auto)
  auto$beta_est)


if(is.null(obj.bigSNP)){
  obj.bigSNP <- snp_attach("PGC_imput_QC_new_filtered.rds")
}
genotype <- obj.bigSNP$genotypes



#### C) JUNTAR VCOVARIABLES DE AN?LISIS Y PRS A DISTINTOS UMBRALES DE P VALOR. Primera parte

# Calculamos mediante PCA las componentes de ancestralidad.
setwd("/media/jjanssen/DATA/projects/genetics/Metodo LDpred2/Outcome/")

system2('/media/jjanssen/DATA/projects/genetics/ARCHIVOS NECESARIOS/plink', args=c("--bfile PGC_imput_QC_new_filtered", "--indep-pairwise 500 1 0.1", "--out PCA_PRUNED"))
system2('/media/jjanssen/DATA/projects/genetics/ARCHIVOS NECESARIOS/plink', args=c("--bfile PGC_imput_QC_new_filtered", "--extract PCA_PRUNED.prune.in", "--make-bed", "--out PGC_imput_QC_new_filtered_PRUNED"))
system2('/media/jjanssen/DATA/projects/genetics/ARCHIVOS NECESARIOS/plink', args=c("--bfile PGC_imput_QC_new_filtered_PRUNED", "--pca", "--out PCA"))

# Seleccionamos el archivo de eigenvectors (valores de PCA para las primeras 20 PC en los sujetos de la target sample) y exportamos las 10 primeras
PCA <- read.table('PCA.eigenvec', header=FALSE, stringsAsFactors = FALSE)  
colnames(PCA) = c("FID", "IID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20")
PCA_final <-PCA[ ,c("FID", "IID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10")]
write.table(PCA_final, "PCA.txt", sep="\t", quote=FALSE, row.names=FALSE)
# Usamos la informaci?n fenot?pica en el archivo FAM. Y creamos un archivo con las covariables (Edad, sexo, componentes principales)
familyInfo <- read.table("PGC_imput_QC_new_filtered.fam", header=FALSE)
familyInfo <- data.frame(FID= familyInfo$V1, IID = familyInfo$V2, scode= familyInfo$V5, dcode = familyInfo$V6)
# scode = 1 (Males), scode = 2 (Females)
# dcode = 1 (Case),  dcode = 2 (Control), dcode = -9 (No Diagnosis Info)
head(familyInfo)
PCA       <- read.table("PCA.txt", header=T, stringsAsFactors = FALSE)
covars    <- merge(PCA, familyInfo, by = c("IID", "FID"))
age       <- read.table('/media/jjanssen/DATA/projects/genetics/ARCHIVOS NECESARIOS/EDADES.txt', header=T, stringsAsFactors = FALSE)
covars    <- merge(covars, age, by = c("IID", "FID"))

write.table(covars, "covariables_SCZ.txt", row.names = FALSE, col.names = TRUE, dec = ".")
covars <- read.table("covariables_SCZ.txt", header = TRUE, stringsAsFactors = FALSE)



#### D) CALCULAMOS LAS TODAS LAS PRS

ind.test <- 1:nrow(genotype)
pred_auto <-
  big_prodMat(genotype,
              beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`,
              ncores = NCORES)


# scale the PRS generated from AUTO
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-
  rowMeans(beta_auto[,
                     abs(pred_scaled -
                           median(pred_scaled)) <
                       3 * mad(pred_scaled)])

# Obtenemos las PRS finales en pred_auto
pred_auto <-  
  big_prodVec(genotype,
              final_beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)


# Guardamos las PRS como un archivo .RData 
score <- pred_auto
saveRDS(score, file="score.RData")



#### E) CONTROL DE CALIDAD Y ORDENAR INDIVIDUOS

covars <- read.table("covariables_SCZ.txt", header = TRUE, stringsAsFactors = FALSE)
covars_filter <- covars[covars$dcode %in% c(1, 2),] # Los pheno -9 se quitan.
variables_SCZ <- covars_filter
score <- readRDS("score.RData")

# hemos hecho un control de calidad en variables_SCZ (quitando los individuos con -9 en dcode)
# por lo que debemos eliminar esos individuos también en fam.order y score

fam.order_filtered <- fam.order[fam.order$FID %in% variables_SCZ$FID]

# ordenamos los individuos en variables_SCZ según el orden de fam_order (que es el mismo que el de score)
SCZ = variables_SCZ[match(fam.order_filtered$FID, variables_SCZ$FID),]

index_delete <- which(obj.bigSNP$fam$affection == -9)  #31 individuos

# creamos el archivo final con todas las covariables y scores de cada individuo
SCZ$PRS <- score[-c(index_delete)]
colnames(SCZ)[16] <- "SCORE"

write.table(SCZ, "SCZ.txt", dec= ".", row.names = FALSE)



### F) REGRESI?N LOG?STICA PARA EVALULAR PRS 

library(rms)
library(pROC)
variables_SCZ <- read.table('SCZ.txt', header=TRUE, stringsAsFactors = FALSE)  

# Definimos variables de regresi?n log?stica. Metemos las 10 primeras PCs, sexo y edad, adem?s de PRS.
pheno <- as.numeric(variables_SCZ$dcode)
PC1  <- variables_SCZ$C1
PC2  <- variables_SCZ$C2
PC3  <- variables_SCZ$C3
PC4  <- variables_SCZ$C4
PC5  <- variables_SCZ$C5
PC6  <- variables_SCZ$C6
PC7  <- variables_SCZ$C7
PC8  <- variables_SCZ$C8
PC9  <- variables_SCZ$C9
PC10 <- variables_SCZ$C10
SEX  <- variables_SCZ$scode
AGE  <- variables_SCZ$AGE


# Comparamos un modelo de regresi?n solo incluyendo covariables frente a otro que incluye adem?s el PRS. De esa comparaci?n extraemos la varianza explicada (pseudoR2) por la varible PRS
# El P valor de la varible PRS la extraemos del modelo log?stico (glm)
# Hacemos un bucle para la comparaci?n en cada uno de los umbrales de significaci?n seleccionados. habr? un valor de R2 y P para cada umbral. Y adem?s, en cada permutaci?n
# Sacamos el valor de AUC tb, pero en ese caso usamos un modelo sin covaribles

H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)
# H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX , scale = TRUE) #si no quisiera controlar por edad
print(H0)


dataRanges <- c() 
score <- variables_SCZ[,"SCORE"]  
H1 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, scale = TRUE)
DELTA_R2 <- H1$stats[10] - H0$stats[10]
R2 <- DELTA_R2 * 100  

phenoN <- pheno - 1 #  cmabiamos de 2/1 a 1/0 para el modelo de regresi?n
model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, family = binomial())
P <- summary(model)

AUC = roc(pheno, score)

# Guardamos los rsultados de la regresi?n en el dataframe
data_SCZ_lrm <- data.frame(R2 = R2, P = P$coefficients[14,4], AUC = AUC$auc)



### EXPORTAMOS EL GR?FICO DE LA CURVA ROC PARA EL PRS DEL UMBRAL P<0.1 COMO MEJOR PREDICTOR
## m?s info: https://rpubs.com/Wangzf/pROC ; https://github.com/xrobin/pROC

rocobj <- plot.roc(variables_SCZ$dcode, variables_SCZ$SCORE,
                   main = "Confidence intervals", 
                   percent=TRUE,
                   ci = TRUE,                  # compute AUC (of AUC by default)
                   print.auc = TRUE)           # print the AUC (will contain the CI)
ciobj <- ci.se(rocobj,                         # CI of sensitivity
               specificities = seq(0, 100, 5)) # over a select set of specificities
plot(ciobj, type = "shape", col = "#1c61b6AA")     # plot as a blue shape
plot(ci(rocobj, of = "thresholds", thresholds = "best")) # add one threshold

dev.copy(png,'ROC_plot.png')
dev.off()


#### G) CAMBIAMOS LA R2 A LIABILITY SCALE 
## Esto lo hacemos porque estamos ustilizando una pseudoR2, por ser un rasgo dicot?mico (caso/control), y hay que ajustar este valor de R2 por prevalencia
## m?s info: http://www.nealelab.is/blog/2017/9/13/heritability-201-types-of-heritability-and-how-we-estimate-it


## Usamos una prevalencia estimada para SCZ del 1% en la poblaci?n general, y una proporci?n de casos en la muestra target de un 55%
k = 0.01
p = 0.5525

h2l_R2N <- function(k, r2n, p)  # Creamos la funci?n para cambiar la pseudoR2 a R2 en liability scale
{
  # k baseline disease risk
  # r2n Nagelkerke's attributable to genomic profile risk score
  # proportion of sample that are cases
  # calculates proportion of variance explained on the liability scale
  # from ABC at http://www.complextraitgenomics.com/software/
  # Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  h2l_R2N <- cc * e * r2n / (1 + cc * e * theta * r2n)
}

data_SCZ_lrm$R2_liab = (h2l_R2N(k, ((data_SCZ_lrm$R2)*0.01), p))*100
write.table(data_SCZ_lrm, "data_SCZ_lrm.txt", dec= ",", row.names = FALSE)



### H) REGRESI?N LOG?STICA PARA EVALULAR PRS  -  PERO HACEMOS PERMUTACIONES SOBRE UN RANDOM DE PSEUDOCASOS Y PSEUDOCONTROLES. bootstrap distribution (with replacement) sobre el mismo n?mero de casos y controles. 
## En algunas ocasiones puede ser necesario exportar no valores concretos sino intervalos de confianza. De este modo, conseguimos intervalos permutando sujetos de la propia muestra

library("dplyr")

## Leemos archivo
table <- read.table("SCZ.txt", header=T, stringsAsFactors = FALSE)

## Establecemos el n?mero de simulaciones que queremos y hacemos el loop
n.simul <- 1000                

dataRanges <- c() 
# Seleccionamos los 1927 casos y 1561 controles pero al azar con remplazamiento. Por tanto, se repetir?n algunos y faltar?n otros en cada permutaci?n. Sirve como an?lisis de sensibilidad
for(i in 1:n.simul) 
{
  cases = subset(table, dcode == 2)
  controls = subset(table, dcode == 1)
  pseudocases = cases[sample(nrow(cases),size = 1927, replace = TRUE),]
  pseudocontrols = controls[sample(nrow(cases),size = 1561, replace = TRUE),]
  table2 = rbind(pseudocases,pseudocontrols)
  
  #Variables_SCZ for logistic regression (check names)
  pheno <- as.numeric(table2$dcode)
  PC1  <- table2$C1
  PC2  <- table2$C2
  PC3  <- table2$C3
  PC4  <- table2$C4
  PC5  <- table2$C5
  PC6  <- table2$C6
  PC7  <- table2$C7
  PC8  <- table2$C8
  PC9  <- table2$C9
  PC10 <- table2$C10
  SEX  <- table2$scode
  AGE  <- table2$AGE
  
  # Comparamos un modelo de regresi?n solo incluyendo covariables frente a otro que incluye adem?s el PRS. De esa comparaci?n extraemos la varianza explicada (pseudoR2) por la varible PRS
  # El P valor de la varible PRS la extraemos del modelo log?stico (glm)
  # Hacemos un bucle para la comparaci?n en cada uno de los umbrales de significaci?n seleccionados. habr? un valor de R2 y P para cada umbral. Y adem?s, en cada permutaci?n
  
  H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)
  
  
  score <- table2[,"SCORE"]  # tenemos que ponerle esta sintaxis para que coja score1 y no score10,11,12 a la vez
  H1 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, scale = TRUE)
  DELTA_R2 <- H1$stats[10] - H0$stats[10]
  R2 <- DELTA_R2 * 100  
  
  AUC = roc(pheno, score)
  
  # Guardamos los resultados de la regresi?n en el dataframe
  df <- data.frame(R2 = R2, AUC = AUC$auc, perm = paste0( "perm_", i))
  dataRanges <- rbind(dataRanges, df)
}

data_SCZ_lrm_perm = dataRanges

# A?adimos el valor en liability scale
data_SCZ_lrm_perm$R2_liab = (h2l_R2N(k, ((data_SCZ_lrm_perm$R2)*0.01), p))*100
write.table(data_SCZ_lrm_perm, "data_SCZ_lrm_perm.txt", dec= ",", row.names = FALSE)



### Creamos una tabla donde aparecen los resultados en promedio, con los intervalos de confianza

table <- read.table("data_SCZ_lrm_perm.txt", header=T, stringsAsFactors = FALSE, dec = ",")
table$AUC = as.numeric(table$AUC)

datos_nuevos <- c() 

MEAN_liab = mean(table$R2_liab)
se = sd(table$R2_liab)
CI_inf_liab = MEAN_liab - (se*1.96)
CI_sup_liab = MEAN_liab + (se*1.96)

MEAN_r2N = mean(table$R2)
se = sd(table$R2)
CI_inf_r2N = MEAN_r2N - (se*1.96)
CI_sup_r2N = MEAN_r2N + (se*1.96)

MEAN_AUC = mean(table$AUC)
se = sd(table$AUC)
CI_inf_AUC = MEAN_AUC - (se*1.96)
CI_sup_AUC = MEAN_AUC + (se*1.96)

datos_nuevos <- data.frame(R2 = MEAN_r2N, CI_inf_R2 = CI_inf_r2N, CI_sup_R2 = CI_sup_r2N, R2_liab = MEAN_liab, CI_inf_R2_liab = CI_inf_liab, CI_sup_R2_liab = CI_sup_liab, AUC = MEAN_AUC, CI_inf_AUC = CI_inf_AUC, CI_sup_AUC = CI_sup_AUC )

write.table(datos_nuevos, "FINAL_data_SCZ_lrm_perm_CIs.txt", dec= ",", row.names = FALSE)



### I) RESULTADOS DE LA COMPARACI?N POR DECILES.
## Para ello, asignamos un decil a cada sujeto dentro del umbral de significaci?n de mayor R2 (en este caso P < v0.1).
## Despu?s hacemos una comparaci?n de la probabilidad de ser un caso con SCZ dentro de cada decilr respecto al medio (5?) o respecto al decil m?s bajo (1?)

library(base)
library(rms)
library(dplyr)

## Cargamos la tabla con todos los scores y covariables, y creamos la varible decil asignando a cada sujeto un valor del 1..10 en funci?n de su score en P<0.1
SCZ <- read.table('SCZ.txt', header=TRUE, stringsAsFactors = FALSE)  
SCZ$decile <- ntile(SCZ$SCORE, 10)


## COMPARACI?N FRRENTE AL DECIL 5 (PROMEDIO)
## Hacemos una subselecci?n de sujetos en base a su decil, inlcuyendo siempre el 5? para hacer la comparaci?n
y <- c(1,2,3,4,6,7,8,9,10)
for(i in y)
{
  variables_SCZ <- SCZ[ which(SCZ$decile == 5 | SCZ$decile == i), ]
  write.table(variables_SCZ, paste0("data",i), dec= ",", row.names = FALSE)
}

## Cargamos cada archivo y calculamos el OR de la comparaci?n de probabilidad de ser caso en el decil seeleccionado en comparado con el decil 5
dataRanges <- c() 
for(i in y)
{
  variables_SCZ <- read.table(paste0("data",i), header = TRUE, stringsAsFactors = FALSE, dec = ",")
  variables_SCZ$comparison <- ifelse(variables_SCZ$decile == 5, 0, 1)
  
  # Definimos variables de regresi?n log?stica. Metemos las 10 primeras PCs, sexo y edad. Metemos la variable PRS del siguiente modo:
  # decil X == 1
  # decil 5 == 0
  
  pheno <- as.numeric(variables_SCZ$dcode)
  PC1  <- variables_SCZ$C1
  PC2  <- variables_SCZ$C2
  PC3  <- variables_SCZ$C3
  PC4  <- variables_SCZ$C4
  PC5  <- variables_SCZ$C5
  PC6  <- variables_SCZ$C6
  PC7  <- variables_SCZ$C7
  PC8  <- variables_SCZ$C8
  PC9  <- variables_SCZ$C9
  PC10 <- variables_SCZ$C10
  SEX   <- variables_SCZ$scode
  AGE   <- variables_SCZ$AGE
  SCORE <- variables_SCZ$comparison
  
  phenoN <- pheno - 1 #Coded for glm
  model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + SCORE, family = binomial())
  P <- summary(model)
  
  coef <- summary(model)$coeff
  OR = exp(coef[14,1])
  CI_inf = exp(coef[14,1] - 1.96*coef[14,2])
  CI_sup = exp(coef[14,1] + 1.96*coef[14,2])
  
  # Guardamos en archivo y exportamos
  df <- data.frame(OR = OR, CI_inf = CI_inf, CI_sup = CI_sup, P = P$coefficients[14,4], Decile = paste0( "Decile_", i))
  dataRanges <- rbind(dataRanges, df)
}
data_decils_lgm = dataRanges
write.table(data_decils_lgm, "data_SCZ_lrm_DECIL_VS_5th.txt", dec= ",", row.names = FALSE)



## COMPARACI?N FRRENTE AL DECIL UNO (MENOR RIESGO POSIBLE)
## Hacemos una subselecci?n de sujetos en base a su decil, inlcuyendo siempre el 5? para hacer la comparaci?n
y <- c(2,3,4,5,6,7,8,9,10)
for(i in y)
{
  variables_SCZ <- SCZ[ which(SCZ$decile == 1 | SCZ$decile == i), ]
  write.table(variables_SCZ, paste0("data",i), dec= ",", row.names = FALSE)
}

## Cargamos cada archivo y calculamos el OR de la comparaci?n de porbabilidad de ser caso en el decil seeleccionado en comparado con el decil 1
dataRanges <- c() 
for(i in y)
{
  variables_SCZ <- read.table(paste0("data",i), header = TRUE, stringsAsFactors = FALSE, dec = ",")
  variables_SCZ$comparison <- ifelse(variables_SCZ$decile == 1, 0, 1)
  
  # Definimos variables de regresi?n log?stica. Metemos las 10 primeras PCs, sexo y edad. Metemos la variable PRS del siguiente modo:
  # decil X == 1
  # decil 1 == 0
  
  pheno <- as.numeric(variables_SCZ$dcode)
  PC1  <- variables_SCZ$C1
  PC2  <- variables_SCZ$C2
  PC3  <- variables_SCZ$C3
  PC4  <- variables_SCZ$C4
  PC5  <- variables_SCZ$C5
  PC6  <- variables_SCZ$C6
  PC7  <- variables_SCZ$C7
  PC8  <- variables_SCZ$C8
  PC9  <- variables_SCZ$C9
  PC10 <- variables_SCZ$C10
  SEX   <- variables_SCZ$scode
  AGE   <- variables_SCZ$AGE
  SCORE <- variables_SCZ$comparison
  
  phenoN <- pheno - 1 #Coded for glm
  model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + SCORE, family = binomial())
  P <- summary(model)
  
  coef <- summary(model)$coeff
  OR = exp(coef[14,1])
  CI_inf = exp(coef[14,1] - 1.96*coef[14,2])
  CI_sup = exp(coef[14,1] + 1.96*coef[14,2])
  
  # Guardamos en archivo y exportamos
  df <- data.frame(OR = OR, CI_inf = CI_inf, CI_sup = CI_sup, P = P$coefficients[14,4], Decile = paste0( "Decile_", i))
  dataRanges <- rbind(dataRanges, df)
}
data_decils_lgm = dataRanges
write.table(data_decils_lgm, "data_SCZ_lrm_DECIL_VS_1st.txt", dec= ",", row.names = FALSE)





### J) GR?FICOS DE LOS AN?LISIS DE VARIANZA Y COMPARCI?N POR DECILES

library(ggplot2)

####################################
## GR?FICOS DE VARIANZA EXPLICADA ##
####################################

# Resultados de varianza se encuentran en el archivo "FINAL_data_SCZ_lrm_perm_CIs"
# No se crean gr?ficos de varianza explicada al no haber ninguna comparaci?n

#########################################
## GR?FICOS DE COMPARACI?N POR DECILES ##
#########################################

###  COMPARACIONES FRENTE AL 5 DECIL

table2=read.table("data_SCZ_lrm_DECIL_VS_5th.txt", header = T, dec = ",", stringsAsFactors = FALSE)
table2$Decile <- sub("Decile_", "", table2$Decile)
table2$Decile <- factor(table2$Decile, levels = table2$Decile)


# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gr?fico de R2. Redondeamos a 3 d?gitos, menos cuando sea exponencial que sustitu?mos la e por x10
table2$print.p <- round(table2$P, digits = 3)
table2$print.p[!is.na(table2$print.p) & table2$print.p == 0] <- format(table2$P[!is.na(table2$print.p) & table2$print.p == 0], digits = 2)
table2$print.p <- sub("e", "*x*10^", table2$print.p)

# Cambiamos el valor de decil a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = table2, aes(x = factor(Decile), y = OR, fill = -log10(P), color = -log10(P))) +
  
  # Especificamos que queremos el valor de P-value sobre el l?mite superior de la errorbar, e inclinado
  geom_text(aes(label = paste(print.p), y = CI_sup), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T, color = "black") +
  
  # Especificamos el rango de la escala Y como 1.6 por el m?ximo de R2, para dejar suficiente espacio para barras y p valores (esto puede tener que retocarse).
  scale_y_continuous(limits = c(0, max(table2$OR) * 1.6)) +
  
  # Especificamos los labels de los ejes X e Y.
  xlab(expression(paste("PRS decile"))) +
  ylab(expression(paste("OR (comparison against 5th decile)"))) +
  
  # configuramos la errorbarplot especificando su posici?mn en eje X e Y, as? como formato.
  geom_errorbar( aes( x=factor(Decile), ymin=CI_inf, ymax=CI_sup, color = -log10(P)), width=0.4, alpha=0.9, size=1.3) + geom_point() +
  
  # A?adimos una linea horizontal indicando el valor neutro (OR = 1)
  geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.8) +
  
  # Especificamos los colores de las barras en funci?n del grado de asociaci?n (seg?n pusimos arriba en fill een barplot). Ponemos un midpoint seg?n el grado de asociaci?n que tengamos. Y la leyenda relativa
  scale_colour_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-3,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +
  
  # Le a?ado colores a los puntos tb del mismo modo que las barras, pero le quito la leyenda para que no aparezcan las dos (guide = "none")
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-3,
    guide = "none") +
  
  # Cambios est?ticos del gr?fico seg?n queramos
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gr?fico
ggsave("DECILES_COMP_5th.png", height = 10, width = 12)




###  COMPARACIONES FRENTE AL PRIMER DECIL

table2=read.table("data_SCZ_lrm_DECIL_VS_1st.txt", header = T, dec = ",", stringsAsFactors = FALSE)
table2$Decile <- sub("Decile_", "", table2$Decile)
table2$Decile <- factor(table2$Decile, levels = table2$Decile)


# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gr?fico de R2. Redondeamos a 3 d?gitos, menos cuando sea exponencial que sustitu?mos la e por x10
table2$print.p <- round(table2$P, digits = 3)
table2$print.p[!is.na(table2$print.p) & table2$print.p == 0] <- format(table2$P[!is.na(table2$print.p) & table2$print.p == 0], digits = 2)
table2$print.p <- sub("e", "*x*10^", table2$print.p)

# Cambiamos el valor de decil a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = table2, aes(x = factor(Decile), y = OR, fill = -log10(P), color = -log10(P))) +
  
  # Especificamos que queremos el valor de P-value sobre el l?mite superior de la errorbar, e inclinado
  geom_text(aes(label = paste(print.p), y = CI_sup), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T, color = "black") +
  
  # Especificamos el rango de la escala Y como 1.7 por el m?ximo de R2, para dejar suficiente espacio para barras y p valores (esto puede tener que retocarse).
  scale_y_continuous(limits = c(0, max(table2$OR) * 1.7)) +
  
  # Especificamos los labels de los ejes X e Y.
  xlab(expression(paste("PRS decile"))) +
  ylab(expression(paste("OR (comparison against 1st decile)"))) +
  
  # configuramos la errorbarplot especificando su posici?mn en eje X e Y, as? como formato.
  geom_errorbar( aes( x=factor(Decile), ymin=CI_inf, ymax=CI_sup, color = -log10(P)), width=0.4, alpha=0.9, size=1.3) + geom_point() +
  
  # A?adimos una linea horizontal indicando el valor neutro (OR = 1)
  geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.8) +
  
  # Especificamos los colores de las barras en funci?n del grado de asociaci?n (seg?n pusimos arriba en fill een barplot). Ponemos un midpoint seg?n el grado de asociaci?n que tengamos. Y la leyenda relativa
  scale_colour_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-3,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +
  
  # Le a?ado colores a los puntos tb del mismo modo que las barras, pero le quito la leyenda para que no aparezcan las dos (guide = "none")
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-3,
    guide = "none") +
  
  # Cambios est?ticos del gr?fico seg?n queramos
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gr?fico
ggsave("DECILES_COMP_1st.png", height = 10, width = 12)
datos_nuevos
