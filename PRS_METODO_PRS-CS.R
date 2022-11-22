
# Script para cálculo de PRS empleando PRS-CS
# JUNIO 2022 - Javi Peñas.
# ==============================================
# Description -
#   El script utiliza os datos de genotipado de una muestra objetivo (cohorte de casos - controles CIBERSAM) que se quiere analizar (dataTarget)
#		y los summary statistics de GWAS de esquizofrenia (dataDiscovery)
#
# ==============================================




# R
library(base)

### A) LEEMOS LOS ARCHIVOS QUE VAMOS A UTILIZAR
# summary statistics (SCZ)
setwd("/media/jjanssen/DATA/projects/genetics/ARCHIVOS_NECESARIOS")
dataDiscovery <- read.table('PGC3_leave_out_cibersam', header=TRUE, stringsAsFactors = FALSE)  ## 7583883 variantes

# Cohorte CIBERSAM (target sample)
setwd("./PGC_imput_QC_filtered")
dataTarget <- read.table('PGC_imput_QC_filtered.bim', header = FALSE, stringsAsFactors = FALSE) ## 9340561 variantes





### B) FILTRAMOS DISCOVERY ###
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

# Filtro columnas que no necesito. Me quedo con la informacion relevante, y renombro columnas
temp_discovery3<-dataDiscovery[,c("CHR", "SNP", "BP", "A1", "A2", "FRQ_A_65397", "FRQ_U_92509", "OR", "SE", "P", "ngt", "POS")]
colnames(temp_discovery3)<-c("CHR", "SNP_ID", "BP", "A1", "A2", "FRQ_A_65397", "FRQ_U_92509", "OR", "SE", "P", "ngt", "SNP")
dataDiscovery <- temp_discovery3
rm(list=ls(pattern="^temp"))     ## Elimino objetis temporales




### C) FILTRAMOS TARGET ##
head(dataTarget)

# Cambio nombre de columnas
colnames(dataTarget) <- c("V1", "varID", "V3", "pos", "A1", "A2")

# Creo un archivo con frecuencia de cada variante por si aparece repetida, y eliminamos repetidas
FRECUENCIAS_target<-as.data.frame(table(dataTarget$varID))
temp_target_UNIQ<-dataTarget[dataTarget$varID %in% FRECUENCIAS_target$Var1[FRECUENCIAS_target$Freq < 2], ]
dataTarget <- temp_target_UNIQ  ## 9340561 variantes
rm(list=ls(pattern="^temp"))         ## Elimino objetis temporales
rm(list=ls(pattern="^FRECUENCIAS"))  ## Elimino objetis temporales





### D) HACEMOS LA INTERSECCIÓN ENTRE LOS DOS ARCHIVOS PARA QUEDARNOS CON LAS VARIANTES PRESENTES EN AMBOS ARCHIVOS, y exportamos los archivos para la predicción

# Me quedo con las variantes de un archivo presentes en el otro, y les cambio el nombre para distinguirlos
dataTarget_filtered <- dataTarget[which(dataTarget$varID %in% dataDiscovery$SNP),]                 ## 5692221 variantes
dataDiscovery_filtered <- dataDiscovery[which(dataDiscovery$SNP %in% dataTarget_filtered$varID),]  ## 5692221 variantes

## como tenemos los archivos de referncia de 1000G en formato código de rs, usamos esos valores e los archivos. Para ello:
data_ALL = merge(dataTarget_filtered,dataDiscovery_filtered, by.x ="varID", by.y="SNP")

# Exportamos el archivo de discovery antes de hacer clumping
setwd("/media/jjanssen/DATA/projects/genetics/Metodo PRS-CS/Outcome")
dataDiscovery_filtered_FINAL<-data_ALL[,c("SNP_ID","A1.y", "A2.y", "OR", "P")]
colnames(dataDiscovery_filtered_FINAL)<-c("SNP","A1", "A2", "OR", "P")
write.table(dataDiscovery_filtered_FINAL, "SUMSTATS.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Exportamos las variantes filtradas en la target
varIDs <- data.frame(data_ALL$varID)

SNP_to_POS <- data.frame(data_ALL$varID, data_ALL$SNP_ID)
# write.table(dataDiscovery_filtered, 'dataDiscovery_filtered.txt', row.names = FALSE, col.names = TRUE)
write.table(varIDs, 'variableIDs_forFiltering.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(SNP_to_POS, 'vSNP_to_POS_forFiltering.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)





#### E) PREPARAMOS ARCHIVOS PARA PRS-CS

# Creamos nuevos archivos de la target sample en formato PLINK .bed .bim and .fam solo con las variantes antes filtradas
# Eliminamos también las variantes dentro de la región del MHC (chr6:26Mb - 33Mb)
setwd("/media/jjanssen/DATA/projects/genetics/ARCHIVOS_NECESARIOS")
system2("./plink", args=c("--bfile ./PGC_imput_QC_filtered/PGC_imput_QC_filtered", '--extract "../Metodo PRS-CS/Outcome/variableIDs_forFiltering.txt"', "--exclude range ./MHC_filter.txt", "--make-bed", '--out "../Metodo PRS-CS/Outcome/PGC_imput_QC_new_filtered"'))


## cambiamos el nuevo BIM por los códigos SNP en lugar de CHR:POS
setwd("/media/jjanssen/DATA/projects/genetics/Metodo PRS-CS/Outcome")
dataTarget <- read.table('PGC_imput_QC_new_filtered.bim', header = FALSE, stringsAsFactors = FALSE) 
data_BIM = merge(dataTarget,SNP_to_POS, by.x ="V2", by.y="data_ALL.varID")

## ordenamos por CHR y POS para que el orden del nuevo BIM coincida con lo esperado en el BED file 
data_BIM <- data_BIM[order(as.numeric(as.character(data_BIM$V1)),as.numeric(as.character(data_BIM$V4))),]
data_BIM <- data_BIM[,c("V1", "data_ALL.SNP_ID", "V3", "V4", "V5", "V6")]
write.table(data_BIM, 'PGC_imput_QC_new_filtered.bim', row.names = FALSE, col.names = FALSE, quote = FALSE) # sustituímos este archivo por el de la target que tiene códigos en chr:pos



#### F) INSTALAMOS PRS-CS, CORREMNOS PRS-CS y CANLCULAMOS PRS POR SUJETO

#### INSTALAMOS PROGRAMA Y CARGAMOS PAQUETES
## En python

git clone https://github.com/getian107/PRScs.git

## PRScs requires Python packages scipy (https://www.scipy.org/) and h5py (https://www.h5py.org/) installed.

conda activate myenv2
conda install h5py
conda install scipy


## ref_dir: diretorio del archivo de referencia. DSescargamos el archivo de 1000G de la población en cuestión para ese programa (https://github.com/getian107/PRScs)
## bim_prefix: directorio al Archivo BIM de la target
## sst_file:sumstats en formato SNP; A1; A2; OR(o BETA); P
## n_gwas: N
## dir output
## ressto de parámetros opcionales

./PRScs/PRScs.py \
--ref_dir=ldblk_1kg_eur \
--bim_prefix=Outcome/PGC_imput_QC_new_filtered \
--sst_file=Outcome/SUMSTATS.txt \
--n_gwas=157906 \
--out_dir=results 

cat results_pst*.txt >> scores.txt ## PRSCS nos da valores de score por SNP, por lo que tenemos que concatenar los archivos de todos los cromosomas y formatearlo para usarlo en PLINK en la función --score
 
## R

## cargamos el archivo de scores por SNP y lo formateamos para plink

setwd("/media/jjanssen/DATA/projects/genetics/Metodo PRS-CS/Outcome")
scores <- read.table('../scores.txt', header = FALSE, stringsAsFactors = FALSE) 
scores_FINAL_rs <-scores[,c("V2", "V4", "V6")]

SNP_to_POS <- read.table("vSNP_to_POS_forFiltering.txt", header=F, stringsAsFactors = FALSE)
scores_FINAL = merge(scores_FINAL_rs,SNP_to_POS, by.x ="V2", by.y="V2")
scores_FINAL <- scores_FINAL[,c("V1","V4","V6")]
write.table(scores_FINAL, 'scores_FINAL.txt', row.names = FALSE, col.names = T, quote = FALSE)


# Calculamos PRS
# Empleamos el archivo SCORE
# Usamos PLINK como en el método clásico pero esta vez sin umbrales.
setwd("/media/jjanssen/DATA/projects/genetics/ARCHIVOS_NECESARIOS")
system2("./plink", args=c("--bfile ./PGC_imput_QC_filtered/PGC_imput_QC_filtered", "--score '../Metodo PRS-CS/Outcome/scores_FINAL.txt'", "--out '../Metodo PRS-CS/Outcome/SCZ_PRS_SCORE'"))





### G) JUNTAR COVARIABLES DE ANÁLISIS Y PRS 

setwd("/media/jjanssen/DATA/projects/genetics/Metodo PRS-CS/Outcome")


system2('/media/jjanssen/DATA/projects/genetics/ARCHIVOS_NECESARIOS/plink', args=c("--bfile PGC_imput_QC_new_filtered", "--indep-pairwise 500 1 0.1", "--out PCA_PRUNED"))
system2('/media/jjanssen/DATA/projects/genetics/ARCHIVOS_NECESARIOS/plink', args=c("--bfile PGC_imput_QC_new_filtered", "--extract PCA_PRUNED.prune.in", "--make-bed", "--out PGC_imput_QC_new_filtered_PRUNED"))
system2('/media/jjanssen/DATA/projects/genetics/ARCHIVOS_NECESARIOS/plink', args=c("--bfile PGC_imput_QC_new_filtered_PRUNED", "--pca", "--out PCA"))
# Seleccionamos el archivo de eigenvectors (valores de PCA para las primeras 20 PC en los sujetos de la target sample) y exportamos las 10 primeras
PCA <- read.table('PCA.eigenvec', header=FALSE, stringsAsFactors = FALSE)  
colnames(PCA) = c("FID", "IID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20")
PCA_final <-PCA[ ,c("FID", "IID", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10")]
write.table(PCA_final, "PCA.txt", sep="\t", quote=FALSE, row.names=FALSE)


# Usamos la información fenotípica en el archivo FAM. Y creamos un archivo con las covariables (Edad, sexo, componentes principales)
familyInfo <- read.table("PGC_imput_QC_new_filtered.fam", header=FALSE)
familyInfo <- data.frame(FID= familyInfo$V1, IID = familyInfo$V2, scode= familyInfo$V5, dcode = familyInfo$V6)

# scode = 1 (Males), scode = 2 (Females)
# dcode = 1 (Case),  dcode = 2 (Control), dcode = -9 (No Diagnosis Info)

head(familyInfo)

PCA       <- read.table("PCA.txt", header=T, stringsAsFactors = FALSE)
covars    <- merge(PCA, familyInfo, by = c("IID", "FID"))
age       <- read.table("/media/jjanssen/DATA/projects/genetics/ARCHIVOS_NECESARIOS/EDADES.txt", header=T, stringsAsFactors = FALSE)
covars    <- merge(covars, age, by = c("IID", "FID"))
write.table(covars, "covariables_SCZ.txt", row.names = FALSE, col.names = TRUE, dec = ".")

covars <- read.table("covariables_SCZ.txt", header = TRUE, stringsAsFactors = FALSE)
covars_filter <- covars[covars$dcode %in% c(1, 2),] # Los pheno -9 se quitan. 
variables_SCZ <- covars_filter

fin <- read.table("SCZ_PRS_SCORE.profile", header = TRUE, stringsAsFactors =  FALSE)
score <- data.frame(FID= fin$FID, IID= fin$IID, fin$SCORE)
colnames(score)[3] <- paste0("SCORE")
variables_SCZ <- merge(variables_SCZ, score, by = c("IID", "FID"))

write.table(variables_SCZ, "SCZ.txt", dec= ",", row.names = FALSE)





### H) REGRESIÓN LOGÍSTICA PARA EVALULAR PRS 

library(rms)
library(pROC)
variables_SCZ <- read.table('SCZ.txt', dec= ",", header=TRUE, stringsAsFactors = FALSE)  

# Definimos variables de regresión logística. Metemos las 10 primeras PCs, sexo y edad, además de PRS.
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
score  <- variables_SCZ$SCORE


# Comparamos un modelo de regresión solo incluyendo covariables frente a otro que incluye además el PRS. De esa comparación extraemos la varianza explicada (pseudoR2) por la varible PRS
# El P valor de la varible PRS la extraemos del modelo logístico (glm)
# Sacamos el valor de AUC tb, pero en ese caso usamos un modelo sin covaribles

H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)
# H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX , scale = TRUE) #si no quisiera controlar por edad
print(H0)


H1 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, scale = TRUE)
DELTA_R2 <- H1$stats[10] - H0$stats[10]
R2 <- DELTA_R2 * 100  
  
phenoN <- pheno - 1 #  cmabiamos de 2/1 a 1/0 para el modelo de regresión
model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, family = binomial())
P <- summary(model)

AUC = roc(pheno, score)
 	
# Guardamos los rsultados de la regresión en el dataframe
df <- data.frame(R2 = R2, P = P$coefficients[14,4], AUC = AUC$auc, LAMBDA = 0.002069138, S = 0.9)
data_SCZ_lrm = df #


### EXPORTAMOS EL GRÁFICO DE LA CURVA ROC PARA EL PRS DEL UMBRAL P<0.1 COMO MEJOR PREDICTOR
## más info: https://rpubs.com/Wangzf/pROC ; https://github.com/xrobin/pROC

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




#### I) CAMBIAMOS LA R2 A LIABILITY SCALE 
## Esto lo hacemos porque estamos ustilizando una pseudoR2, por ser un rasgo dicotómico (caso/control), y hay que ajustar este valor de R2 por prevalencia
## más info: http://www.nealelab.is/blog/2017/9/13/heritability-201-types-of-heritability-and-how-we-estimate-it

## Usamos una prevalencia estimada para SCZ del 1% en la población general, y una proporción de casos en la muestra target de un 55%
k = 0.01
p = 0.5525

h2l_R2N <- function(k, r2n, p)  # Creamos la función para cambiar la pseudoR2 a R2 en liability scale
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







### J) REGRESIÓN LOGÍSTICA PARA EVALULAR PRS  -  PERO HACEMOS PERMUTACIONES SOBRE UN RANDOM DE PSEUDOCASOS Y PSEUDOCONTROLES. bootstrap distribution (with replacement) sobre el mismo número de casos y controles. 
## En algunas ocasiones puede ser necesario exportar no valores concretos sino intervalos de confianza. De este modo, conseguimos intervalos permutando sujetos de la propia muestra

library("dplyr")

## Leemos archivo
table <- read.table("SCZ.txt", header=T, dec= ",", stringsAsFactors = FALSE)

## Establecemos el número de simulaciones que queremos y hacemos el loop
n.simul <- 1000               

dataRanges <- c() 
# Seleccionamos los 1927 casos y 1561 controles pero al azar con remplazamiento. Por tanto, se repetirán algunos y faltarán otros en cada permutación. Sirve como análisis de sensibilidad
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
	score  <- table2$SCORE

	H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)


	# Comparamos un modelo de regresión solo incluyendo covariables frente a otro que incluye además el PRS. De esa comparación extraemos la varianza explicada (pseudoR2) por la varible PRS
	# El P valor de la varible PRS la extraemos del modelo logístico (glm)
	# Hacemos un bucle para la comparación en cada uno de los umbrales de significación seleccionados. habrá un valor de R2 y P para cada permutación

	H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)
	# H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX , scale = TRUE) #si no quisiera controlar por edad
	print(H0)


	H1 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, scale = TRUE)
	DELTA_R2 <- H1$stats[10] - H0$stats[10]
	R2 <- DELTA_R2 * 100  
  
	phenoN <- pheno - 1 #  cmabiamos de 2/1 a 1/0 para el modelo de regresión
	model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, family = binomial())
	P <- summary(model)

	AUC = roc(pheno, score)
 	
	# Guardamos los rsultados de la regresión en el dataframe
	df <- data.frame(R2 = R2, P = P$coefficients[14,4], AUC = AUC$auc, LAMBDA = 0.002069138, S = 0.9, perm = paste0( "perm_", i))
	data_SCZ_lrm = df #

	H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)
  	dataRanges <- rbind(dataRanges, df)
	}

data_SCZ_lrm_perm = dataRanges



# Añadimos el valor en liability scale
data_SCZ_lrm_perm$R2_liab = (h2l_R2N(k, ((data_SCZ_lrm_perm$R2)*0.01), p))*100
write.table(data_SCZ_lrm_perm, "data_SCZ_lrm_perm.txt", dec= ",", row.names = FALSE)



### Creamos una tabla donde aparecen los resultados en promedio, con los intervalos de confianza

table <- read.table("data_SCZ_lrm_perm.txt", header=T, stringsAsFactors = FALSE, dec = ",")
table$AUC = as.numeric(table$AUC)

 	
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

# Le añadimos el valor de P de la comparación inicial mediante glm sin bootstrap.
table2 <- read.table("data_SCZ_lrm.txt", header=T, stringsAsFactors = FALSE, dec = ",")
P = table2$P

df <- data.frame(R2 = MEAN_r2N, CI_inf_R2 = CI_inf_r2N, CI_sup_R2 = CI_sup_r2N, P = P, R2_liab = MEAN_liab, CI_inf_R2_liab = CI_inf_liab, CI_sup_R2_liab = CI_sup_liab, AUC = MEAN_AUC, CI_inf_AUC = CI_inf_AUC, CI_sup_AUC = CI_sup_AUC )
write.table(df, "FINAL_data_SCZ_lrm_perm_CIs_P.txt", dec= ",", row.names = FALSE)



### K) RESULTADOS DE LA COMPARACIÓN POR DECILES.
## Para ello, asignamos un decil a cada sujeto 
## Después hacemos una comparación de la probabilidad de ser un caso con SCZ dentro de cada decilr respecto al medio (5º) o respecto al decil más bajo (1º)

library(base)
library(rms)
library(dplyr)

## Cargamos la tabla con todos los scores y covariables, y creamos la varible decil asignando a cada sujeto un valor del 1..10 en función de su score en P<0.1
SCZ <- read.table('SCZ.txt', header=TRUE, dec= ",", stringsAsFactors = FALSE)  
SCZ$decile <- ntile(SCZ$SCORE, 10)


## COMPARACIÓN FRRENTE AL DECIL 5 (PROMEDIO)
## Hacemos una subselección de sujetos en base a su decil, inlcuyendo siempre el 5º para hacer la comparación
y <- c(1,2,3,4,6,7,8,9,10)
for(i in y)
	{
  	variables_SCZ <- SCZ[ which(SCZ$decile == 5 | SCZ$decile == i), ]
 	write.table(variables_SCZ, paste0("data",i), dec= ",", row.names = FALSE)
	}

## Cargamos cada archivo y calculamos el OR de la comparación de porbabilidad de ser caso en el decil seeleccionado en comparado con el decil 5
dataRanges <- c() 
for(i in y)
	{
  	variables_SCZ <- read.table(paste0("data",i), header = TRUE, stringsAsFactors = FALSE, dec = ",")
  	variables_SCZ$comparison <- ifelse(variables_SCZ$decile == 5, 0, 1)

  	# Definimos variables de regresión logística. Metemos las 10 primeras PCs, sexo y edad. Metemos la variable PRS del siguiente modo:
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




## COMPARACIÓN FRRENTE AL DECIL UNO (MENOR RIESGO POSIBLE)
## Hacemos una subselección de sujetos en base a su decil, inlcuyendo siempre el 5º para hacer la comparación
y <- c(2,3,4,5,6,7,8,9,10)
for(i in y)
	{
  	variables_SCZ <- SCZ[ which(SCZ$decile == 1 | SCZ$decile == i), ]
 	write.table(variables_SCZ, paste0("data",i), dec= ",", row.names = FALSE)
	}

## Cargamos cada archivo y calculamos el OR de la comparación de porbabilidad de ser caso en el decil seeleccionado en comparado con el decil 1
dataRanges <- c() 
for(i in y)
	{
  	variables_SCZ <- read.table(paste0("data",i), header = TRUE, stringsAsFactors = FALSE, dec = ",")
  	variables_SCZ$comparison <- ifelse(variables_SCZ$decile == 1, 0, 1)

  	# Definimos variables de regresión logística. Metemos las 10 primeras PCs, sexo y edad. Metemos la variable PRS del siguiente modo:
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









### L) GRÁFICOS DE COMPARCIÓN POR DECILES

### gráficos de varianza explicada no tienen sentido porque no hay diferentes p thresholds

library(ggplot2)



#########################################
## GRÁFICOS DE COMPARACIÓN POR DECILES ##
#########################################

###  COMPARACIONES FRENTE AL 5 DECIL

table2=read.table("data_SCZ_lrm_DECIL_VS_5th.txt", header = T, dec = ",", stringsAsFactors = FALSE)
table2$Decile <- sub("Decile_", "", table2$Decile)
table2$Decile <- factor(table2$Decile, levels = table2$Decile)


# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gráfico de R2. Redondeamos a 3 dígitos, menos cuando sea exponencial que sustituímos la e por x10
table2$print.p <- round(table2$P, digits = 3)
table2$print.p[!is.na(table2$print.p) & table2$print.p == 0] <- format(table2$P[!is.na(table2$print.p) & table2$print.p == 0], digits = 2)
table2$print.p <- sub("e", "*x*10^", table2$print.p)

# Cambiamos el valor de decil a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = table2, aes(x = factor(Decile), y = OR, fill = -log10(P), color = -log10(P))) +

    # Especificamos que queremos el valor de P-value sobre el límite superior de la errorbar, e inclinado
    geom_text(aes(label = paste(print.p), y = CI_sup), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T, color = "black") +

	# Especificamos el rango de la escala Y como 1.6 por el máximo de R2, para dejar suficiente espacio para barras y p valores (esto puede tener que retocarse).
    scale_y_continuous(limits = c(0, max(table2$OR) * 1.7)) +

	# Especificamos los labels de los ejes X e Y.
    xlab(expression(paste("PRS decile"))) +
    ylab(expression(paste("OR (comparison against 5th decile)"))) +

	# configuramos la errorbarplot especificando su posiciómn en eje X e Y, así como formato.
	geom_errorbar( aes( x=factor(Decile), ymin=CI_inf, ymax=CI_sup, color = -log10(P)), width=0.4, alpha=0.9, size=1.3) + geom_point() +

	# Añadimos una linea horizontal indicando el valor neutro (OR = 1)
	geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.8) +

	# Especificamos los colores de las barras en funciòn del grado de asociación (según pusimos arriba en fill een barplot). Ponemos un midpoint según el grado de asociación que tengamos. Y la leyenda relativa
    scale_colour_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +

	# Le añado colores a los puntos tb del mismo modo que las barras, pero le quito la leyenda para que no aparezcan las dos (guide = "none")
        scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        guide = "none") +

    # Cambios estéticos del gráfico según queramos
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))


# Guardamos el gráfico
ggsave("DECILES_COMP_5th.png", height = 10, width = 12)




###  COMPARACIONES FRENTE AL PRIMER DECIL

table2=read.table("data_SCZ_lrm_DECIL_VS_1st.txt", header = T, dec = ",", stringsAsFactors = FALSE)
table2$Decile <- sub("Decile_", "", table2$Decile)
table2$Decile <- factor(table2$Decile, levels = table2$Decile)


# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gráfico de R2. Redondeamos a 3 dígitos, menos cuando sea exponencial que sustituímos la e por x10
table2$print.p <- round(table2$P, digits = 3)
table2$print.p[!is.na(table2$print.p) & table2$print.p == 0] <- format(table2$P[!is.na(table2$print.p) & table2$print.p == 0], digits = 2)
table2$print.p <- sub("e", "*x*10^", table2$print.p)

# Cambiamos el valor de decil a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = table2, aes(x = factor(Decile), y = OR, fill = -log10(P), color = -log10(P))) +

    # Especificamos que queremos el valor de P-value sobre el límite superior de la errorbar, e inclinado
    geom_text(aes(label = paste(print.p), y = CI_sup), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T, color = "black") +

	# Especificamos el rango de la escala Y como 1.7 por el máximo de R2, para dejar suficiente espacio para barras y p valores (esto puede tener que retocarse).
    scale_y_continuous(limits = c(0, max(table2$OR) * 1.9)) +

	# Especificamos los labels de los ejes X e Y.
    xlab(expression(paste("PRS decile"))) +
    ylab(expression(paste("OR (comparison against 1st decile)"))) +

	# configuramos la errorbarplot especificando su posiciómn en eje X e Y, así como formato.
	geom_errorbar( aes( x=factor(Decile), ymin=CI_inf, ymax=CI_sup, color = -log10(P)), width=0.4, alpha=0.9, size=1.3) + geom_point() +

	# Añadimos una linea horizontal indicando el valor neutro (OR = 1)
	geom_hline(yintercept=1, linetype="dashed", color = "grey", size=0.8) +

	# Especificamos los colores de las barras en funciòn del grado de asociación (según pusimos arriba en fill een barplot). Ponemos un midpoint según el grado de asociación que tengamos. Y la leyenda relativa
    scale_colour_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +

	# Le añado colores a los puntos tb del mismo modo que las barras, pero le quito la leyenda para que no aparezcan las dos (guide = "none")
        scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-3,
        guide = "none") +

    # Cambios estéticos del gráfico según queramos
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gráfico
ggsave("DECILES_COMP_1st.png", height = 10, width = 12)





