# Script para cálculo de PRS empleando el método clásico de podado informativo de variantes y umbrales de asociación
# Enero 2022 - Javi Peñas.
# ==============================================
# Description -
#   El script utiliza os datos de genotipado de una muestra objetivo (cohorte de casos - controles CIBERSAM) que se quiere analizar (dataTarget)
#		y los summary statistics de GWAS de esquizofrenia (dataDiscovery)
#
# ==============================================


library(base)

### A) LEEMOS LOS ARCHIVOS QUE VAMOS A UTILIZAR
# summary statistics (SCZ)
setwd("C:/Users/guill/OneDrive - Universidad Rey Juan Carlos (1)/Escritorio/Ingeniería Biomédica URJC/TFG/ARCHIVOS NECESARIOS")
dataDiscovery <- read.table('PGC3_leave_out_cibersam', header=TRUE, stringsAsFactors = FALSE)  ## 7583883 variantes

# Cohorte CIBERSAM (target sample)
setwd("C:/Users/guill/OneDrive - Universidad Rey Juan Carlos (1)/Escritorio/Ingeniería Biomédica URJC/TFG/ARCHIVOS NECESARIOS/PGC_imput_QC_filtered")
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


### D) HACEMOS LA INTERSECCIÓN ENTRE LOS DOS ARCHIVOS PARA QUEDARNOS CON LAS VARIANTES PRESENTES EN AMBOS ARCHIVOS

# Me quedo con las variantes de un archivo presentes en el otro, y les cambio el nombre para distinguirlos
dataTarget_filtered <- dataTarget[which(dataTarget$varID %in% dataDiscovery$SNP),]                 ## 5692221 variantes
dataDiscovery_filtered <- dataDiscovery[which(dataDiscovery$SNP %in% dataTarget_filtered$varID),]  ## 5692221 variantes



### E) CLUMPING - Eliminar variantes en desequilibrio de ligamiento (LD; es decir, no independientes)

# Exportamos el archivo de discovery antes de hacer clumping
setwd("C:/Users/guill/OneDrive - Universidad Rey Juan Carlos (1)/Escritorio/Ingeniería Biomédica URJC/TFG/ARCHIVOS NECESARIOS/Cálculos PRS")
write.table(dataDiscovery_filtered, "discovery_paraCLUMP.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Exportamos las variantes filtradas en la target
varIDs <- data.frame(dataTarget_filtered$varID)
# write.table(dataDiscovery_filtered, 'dataDiscovery_filtered.txt', row.names = FALSE, col.names = TRUE)
write.table(varIDs, 'variableIDs_forFiltering.txt', row.names = FALSE, col.names = FALSE, quote = FALSE)


# Creamos nuevos archivos de la target sample en formato PLINK .bed .bim and .fam solo con las variantes antes filtradas
# Eliminamos también las variantes dentro de la región del MHC (chr6:26Mb - 33Mb)
system2("../plink", args=c("--bfile ../PGC_imput_QC_filtered/PGC_imput_QC_filtered", "--extract variableIDs_forFiltering.txt", "--exclude range ../MHC_filter.txt", "--make-bed", "--out PGC_imput_QC_new_filtered"))


# Hacemos CLUMPING usando  plink con el archivo de la target ya filtrado, y usando como muestra discovery el exportado en los pasos anteriores
# Mediante el clumping vamos a quedarnos con las variantes independientes (no en LD - usando la propia target sample para ello) priorizando las más asociadas con SCZ, a partir de la información de la muestra discovery
# En un primer paso exportamos la lista de variantes resultante y posteriormente creamos un nuevo archivo de la target solo con esas variantes, que serán las únicas que se considerarán para el cálculo de PRS posterior
# Como parámetros de clumping usamos r2 < 0.1 y en una ventana de 500 Kb.
# IMPORTANTE: PLINK necesita reconargs=ocer el ID de la variante (columna SNP) y el valor de asociación de la muestra discovery (columna P). Aquí usamos como ID la columna chr:posición que creamos antes y es común a la discovery y a la target. Si no coinciden los IDs va a dar error
system2("../plink", c("--bfile PGC_imput_QC_new_filtered", "--clump discovery_paraCLUMP.txt", "--clump-p1 1", "--clump-p2 1", "--clump-kb 500", "--clump-r2 0.1", "--out PGC_imput_QC_filtered_CLUMPED"))
CLUMPED_PLINK<-read.table("PGC_imput_QC_filtered_CLUMPED.clumped", header=T)
write.table(CLUMPED_PLINK$SNP, "CLUMPED_PLINK.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
system2("../plink", args=c("--bfile PGC_imput_QC_new_filtered", "--extract CLUMPED_PLINK.txt", "--make-bed", "--out PGC_imput_QC_filtered_CLUMPED"))  ## 210691 variantes 


# NOTAS: BFILE LLAMA A LOS 3 ARCHIVOS QUE TENGAN LOS FORMATOS PLINK
# EXCLUIR VARIANTES EN UN RANGO DE MHC 
# MAKE BED HACER ARCHIVO EN EL MISMO FORMATO PLINK
# R2 te quedas con las que estén menos correlacionadas, más independientes /// p valor p1 lo cojo todo no filtro por p valor
# "--clump discovery_paraCLUMP.txt" de todas las variantes que yo tengo, prioriza la que tenga el p valor más bajo.



### F)  CALCULO DE PUNTUCIONES PRS  #####
head(dataDiscovery_filtered)
head(dataTarget_filtered)

# IMPORTANTE! Si hay OR como variable hay que cambiara a Beta, usando logOR para convertirla a Beta
dataDiscovery_filtered$BETA = log(dataDiscovery_filtered$OR)

# Creamos los archivos necesarios
# IMPORTANTE: Es necesario conocer cuál es el alelo sobre el que está calculado el efecto - si nos equivocamos veremos unos cálculos a la inversa. Consultar en el readme qué alelo es al que sde refiere en el OR
SCZ_SCORE <- dataDiscovery_filtered[,c("SNP", "A1", "BETA")]
SCZ_P <- dataDiscovery_filtered[,c("SNP", "P")]
write.table(SCZ_SCORE, "SCZ_SCORE.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(SCZ_P, "SCZ_P.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Calculamos PRS
# Empleamos el archivo SCORE, el P valor para seleccionar variantes y calcular puntuaciones
# Usamos un archivo de umbrales de P valor (Q_RANGES): En el definimos los umbrales que queremos emplear para distintos cálculos
system2("../plink", args=c("--bfile PGC_imput_QC_filtered_CLUMPED", "--score SCZ_SCORE.txt", "--q-score-file SCZ_P.txt", "--q-score-range ../Q_RANGES.txt", "--out SCZ_PRS_SCORE"))

# q score range tiene en cuenta el p valor (p value threslhold)
# primero calcular para cada umbral (tercera columna), cada vez se incluyen más variantes 
# cálculos de puntuaciones para cada umbral
# puntuaciones para cada individuo y cada umbral (S)
# nombres de indivuos en fam con su score 



### G) JUNTAR VCOVARIABLES DE ANÁLISIS Y PRS A DISTINTOS UMBRALES DE P VALOR.
# Calculamos mediante PCA las componentes de ancestralidad.
system2("../plink", args=c("--bfile PGC_imput_QC_new_filtered", "--indep-pairwise 500 1 0.1", "--out PCA_PRUNED"))
system2("../plink", args=c("--bfile PGC_imput_QC_new_filtered", "--extract PCA_PRUNED.prune.in", "--make-bed", "--out PGC_imput_QC_new_filtered_PRUNED"))
system2("../plink", args=c("--bfile PGC_imput_QC_new_filtered_PRUNED", "--pca", "--out PCA"))
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
age       <- read.table("../EDADES.txt", header=T, stringsAsFactors = FALSE)
covars    <- merge(covars, age, by = c("IID", "FID"))
write.table(covars, "covariables_SCZ.txt", row.names = FALSE, col.names = TRUE, dec = ".")
#G: edad en que fueron diagnosticados
covars <- read.table("covariables_SCZ.txt", header = TRUE, stringsAsFactors = FALSE)
covars_filter <- covars[covars$dcode %in% c(1, 2),] # Los pheno -9 se quitan.
variables_SCZ <- covars_filter
for (j in 1:12)  # Pongo 12 porque he generado 12 PRS para los 12 umbrales del archivo Q_RANGES
{
  fin <- read.table(paste0("SCZ_PRS_SCORE.S",j,".profile"), header = TRUE, stringsAsFactors =  FALSE)
  score <- data.frame(FID= fin$FID, IID= fin$IID, fin$SCORE)
  colnames(score)[3] <- paste0("SCORE_S",j)
  variables_SCZ <- merge(variables_SCZ, score, by = c("IID", "FID"))
}
write.table(variables_SCZ, "SCZ.txt", dec= ".", row.names = FALSE)
# G: generar 1 archivo con el score en cada individuo para los distintos umbrales



### H) REGRESIÓN LOGÍSTICA PARA EVALULAR PRS 

library(rms)
library(pROC)
variables_SCZ <- read.table('SCZ.txt', header=TRUE, stringsAsFactors = FALSE)  

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

#G: I es binario, logistica fucnión de probabilidad de sacar 0 o 1. 

# Comparamos un modelo de regresión solo incluyendo covariables frente a otro que incluye además el PRS. De esa comparación extraemos la varianza explicada (pseudoR2) por la varible PRS
# El P valor de la varible PRS la extraemos del modelo logístico (glm)
# Hacemos un bucle para la comparación en cada uno de los umbrales de significación seleccionados. habrá un valor de R2 y P para cada umbral. Y además, en cada permutación
# Sacamos el valor de AUC tb, pero en ese caso usamos un modelo sin covaribles

H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)
# H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX , scale = TRUE) #si no quisiera controlar por edad
print(H0)

# G: Esdtimar efecto de prs sobre fenotipo, para ello comparar modelo con y sin prs. así se extrae el efecto de prs
# y corregir por covar.; pred de todo menos pred de todo menos PRS.


dataRanges <- c() 
for (j in 1:12) 
{
  score <- variables_SCZ[,grep(paste0("^SCORE_S",j,"$"), colnames(variables_SCZ))]  # tenemos que ponerle esta sintaxis para que coja score1 y no score10,11,12 a la vez
  H1 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, scale = TRUE)
  DELTA_R2 <- H1$stats[10] - H0$stats[10]
  R2 <- DELTA_R2 * 100  
  
  phenoN <- pheno - 1 #  cmabiamos de 2/1 a 1/0 para el modelo de regresión
  model <- glm(phenoN ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, family = binomial())
  P <- summary(model)
  
  AUC = roc(pheno, score)
  
  # Guardamos los rsultados de la regresión en el dataframe
  df <- data.frame(R2 = R2, P = P$coefficients[14,4], AUC = AUC$auc, Range = paste0( "Q_RANGE_S", j))
  dataRanges <- rbind(dataRanges, df)
}

data_SCZ_lrm = dataRanges # Mejor umbral S9 P < 0.1


# G: delta R2 hace la R2 del 1 menos la del 0, resta del R2; al final porcentaje
# glm calcular p valor; al final p valor de la varible prs score , y resta del r2
# AUC de curvas roc en mod log, prob de que se clasifique bien; clasificado como caso o control de forma errónea? mezcla entre sensilbilidad y especificiadad
# cada línea es un umbral con 1 valor de r2 (explica prs), p valor, AUC
# si está por encima de diagonal, predice más que un modelo aleatorio, valor de 0.8 es predicción buena


### EXPORTAMOS EL GRÁFICO DE LA CURVA ROC PARA EL PRS DEL UMBRAL P<0.1 COMO MEJOR PREDICTOR
## más info: https://rpubs.com/Wangzf/pROC ; https://github.com/xrobin/pROC

rocobj <- plot.roc(variables_SCZ$dcode, variables_SCZ$SCORE_S9,
                   main = "Confidence intervals", 
                   percent=TRUE,
                   ci = TRUE,                  # compute AUC (of AUC by default)
                   print.auc = TRUE)           # print the AUC (will contain the CI)
ciobj <- ci.se(rocobj,                         # CI of sensitivity
               specificities = seq(0, 100, 5)) # over a select set of specificities
plot(ciobj, type = "shape", col = "#1c61b6AA")     # plot as a blue shape
plot(ci(rocobj, of = "thresholds", thresholds = "best")) # add one threshold




#### I) CAMBIAMOS LA R2 A LIABILITY SCALE 
## Esto lo hacemos porque estamos ustilizando una pseudoR2, por ser un rasgo dicotómico (caso/control), y hay que ajustar este valor de R2 por prevalencia
## más info: http://www.nealelab.is/blog/2017/9/13/heritability-201-types-of-heritability-and-how-we-estimate-it

# G: R2 era la resta de los R2, es una forma no ajustada matematic; al ser regreesiónm logistica; más casos que en prevalencia; si prev es mayor que proporción de casos
# por lo que no es representativa; R2 en liability scale
# k es prevalencia; p es proporción de casos en la muestra
# crea función para ajustar este R2; añade a la tabla un R2 liability scale

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



# G: bootstrap es análisis de permutaciones para intervalos de confianza teniendo en cuenta la avriabilidad entre individuos, selecciona un grupo de x elementos y hay reemplazamiento
# los 1927 casos son cogidos al azar y con reemplazamiento; 



### J) REGRESIÓN LOGÍSTICA PARA EVALULAR PRS  -  PERO HACEMOS PERMUTACIONES SOBRE UN RANDOM DE PSEUDOCASOS Y PSEUDOCONTROLES. bootstrap distribution (with replacement) sobre el mismo número de casos y controles. 
## En algunas ocasiones puede ser necesario exportar no valores concretos sino intervalos de confianza. De este modo, conseguimos intervalos permutando sujetos de la propia muestra

library("dplyr")

## Leemos archivo
table <- read.table("SCZ.txt", header=T, stringsAsFactors = FALSE)

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
  
  # Comparamos un modelo de regresión solo incluyendo covariables frente a otro que incluye además el PRS. De esa comparación extraemos la varianza explicada (pseudoR2) por la varible PRS
  # El P valor de la varible PRS la extraemos del modelo logístico (glm)
  # Hacemos un bucle para la comparación en cada uno de los umbrales de significación seleccionados. habrá un valor de R2 y P para cada umbral. Y además, en cada permutación
  
  H0 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE, scale = TRUE)
  
  for (j in 1:12) 
  {
    score <- table2[,grep(paste0("^SCORE_S",j,"$"), colnames(table2))]  # tenemos que ponerle esta sintaxis para que coja score1 y no score10,11,12 a la vez
    H1 <- lrm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + SEX + AGE + score, scale = TRUE)
    DELTA_R2 <- H1$stats[10] - H0$stats[10]
    R2 <- DELTA_R2 * 100  
    
    AUC = roc(pheno, score)
    
    # Guardamos los resultados de la regresión en el dataframe
    df <- data.frame(R2 = R2, AUC = AUC$auc, Range = paste0( "Q_RANGE_S", j), perm = paste0( "perm_", i))
    dataRanges <- rbind(dataRanges, df)
  }
}

data_SCZ_lrm_perm = dataRanges

#G: por cada umbral habrá 1000 cáluclos, por haber hecho 1000 permutaciones; y genero intervaos de confianza de todo al tener 1000

# Añadimos el valor en liability scale
data_SCZ_lrm_perm$R2_liab = (h2l_R2N(k, ((data_SCZ_lrm_perm$R2)*0.01), p))*100
write.table(data_SCZ_lrm_perm, "data_SCZ_lrm_perm.txt", dec= ",", row.names = FALSE)



### Creamos una tabla donde aparecen los resultados en promedio, con los intervalos de confianza

table <- read.table("data_SCZ_lrm_perm.txt", header=T, stringsAsFactors = FALSE, dec = ",")
table$AUC = as.numeric(table$AUC)

datos_nuevos <- c() 
for(i in 1:12)
{
  q_range <- table[table$Range == paste0("Q_RANGE_S",i),]
  
  MEAN_liab = mean(q_range$R2_liab)
  se = sd(q_range$R2_liab)
  CI_inf_liab = MEAN_liab - (se*1.96)
  CI_sup_liab = MEAN_liab + (se*1.96)
  
  MEAN_r2N = mean(q_range$R2)
  se = sd(q_range$R2)
  CI_inf_r2N = MEAN_r2N - (se*1.96)
  CI_sup_r2N = MEAN_r2N + (se*1.96)
  
  MEAN_AUC = mean(q_range$AUC)
  se = sd(q_range$AUC)
  CI_inf_AUC = MEAN_AUC - (se*1.96)
  CI_sup_AUC = MEAN_AUC + (se*1.96)
  
  df <- data.frame(Range = paste0( "Q_RANGE_S", i), R2 = MEAN_r2N, CI_inf_R2 = CI_inf_r2N, CI_sup_R2 = CI_sup_r2N, R2_liab = MEAN_liab, CI_inf_R2_liab = CI_inf_liab, CI_sup_R2_liab = CI_sup_liab, AUC = MEAN_AUC, CI_inf_AUC = CI_inf_AUC, CI_sup_AUC = CI_sup_AUC )
  datos_nuevos = rbind(datos_nuevos, df)
}

write.table(datos_nuevos, "data_SCZ_lrm_perm_CIs.txt", dec= ",", row.names = FALSE)


# Le añadimos el valor de P de la comparación inicial mediante glm si bootstrap.
table2 <- read.table("data_SCZ_lrm.txt", header=T, stringsAsFactors = FALSE, dec = ",")
threshold <- read.table("../Q_RANGES.txt", header=FALSE, stringsAsFactors = FALSE, dec = ".")


TABLA_FINAL = cbind(threshold[,3],datos_nuevos,table2[,2])
names(TABLA_FINAL)[names(TABLA_FINAL) == "table2[, 2]"] <- 'P_value' # renombreamos la columna que acabamos de añadir
names(TABLA_FINAL)[names(TABLA_FINAL) == "threshold[, 3]"] <- 'Threshold' # renombreamos la columna que acabamos de añadir



write.table(TABLA_FINAL, "FINAL_data_SCZ_lrm_perm_CIs_P.txt", dec= ",", row.names = FALSE)



# G: deciles represe si agrupo indiv por deciles (por su prs) y digo qué sujeto pertenece a cada decidles, subconj ordenados; cuál es la razón de diferencia de mayor decil y menos (mayor riesgo y menos)



### K) RESULTADOS DE LA COMPARACIÓN POR DECILES.
## Para ello, asignamos un decil a cada sujeto dentro del umbral de significación de mayor R2 (en este caso P < v0.1).
## Después hacemos una comparación de la probabilidad de ser un caso con SCZ dentro de cada decilr respecto al medio (5º) o respecto al decil más bajo (1º)

library(base)
library(rms)
library(dplyr)

## Cargamos la tabla con todos los scores y covariables, y creamos la varible decil asignando a cada sujeto un valor del 1..10 en función de su score en P<0.1
SCZ <- read.table('SCZ.txt', header=TRUE, stringsAsFactors = FALSE)  
SCZ$decile <- ntile(SCZ$SCORE_S9, 10)


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








### L) GRÁFICOS DE LOS ANÁLISIS DE VARIANZA Y COMPARCIÓN POR DECILES

library(ggplot2)

####################################
## GRÁFICOS DE VARIANZA EXPLICADA ##
####################################

# Cargamos el archivo final con todos los valores para exportar
prs.result <- read.table("FINAL_data_SCZ_lrm_perm_CIs_P.txt", header=T, stringsAsFactors = FALSE, dec = ",")

# Generamos un buen formato para lo P-valores del modeo, y ponerlos sobre el gráfico de R2. Redondeamos a 3 dígitos, menos cuando sea exponencial que sustituímos la e por x10
prs.result$print.p <- round(prs.result$P_value, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) & prs.result$print.p == 0] <- format(prs.result$P[!is.na(prs.result$print.p) & prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)


## GRÁFICO PARA R2 EN LIABILITY SCALE
# cargamos el gráfico en ggplot. 
# Cambiamos el valor de threshold a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2_liab)) +
  
  # Especificamos que queremos el valor de P-value sobre el límite superior de la errorbar, e inclinado
  geom_text(aes(label = paste(print.p), y = CI_sup_R2_liab), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T) +
  
  # Especificamos el rango de la escala Y como 1.25 por el máximo de R2_liab, para dejar suficiente espacio para barras y p valores.
  scale_y_continuous(limits = c(0, max(prs.result$R2_liab) * 1.25)) +
  
  # Especificamos los labels de los ejes X e Y en itálica.
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PRS model fit: liability ", R ^ 2))) +
  
  # configuramos la barplot que es el tipo de gráfico que queremos. Configuramos en fill para que el clor sea funció de su grado de asociación.
  geom_bar(aes(fill = -log10(P_value)), stat = "identity") +
  
  # configuramos la errorbarplot especificando su posiciómn en eje X e Y, así como formato.
  geom_errorbar( aes(x=factor(Threshold), ymin=CI_inf_R2_liab, ymax=CI_sup_R2_liab), width=0.4, colour="black", alpha=0.9, size=1.3) + 
  
  # Especificamos los colores de las barras en funciòn del grado de asociación (según pusimos arriba en fill een barplot). Ponemos un midpoint según el grado de asociación que tengamos. Y la leyenda relativa
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-10,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +
  
  # Cambios estéticos del gráfico según queramos
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gráfico
ggsave("R2_LIAB_PLOT.png", height = 16, width = 12)



## GRÁFICO PARA PSEUDO R2 DE NAGGELKERKE.
# cargamos el gráfico en ggplot. 
# Cambiamos el valor de threshold a FACTOR para que se pueda mantener el orden creciente de P valor.
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
  
  # Especificamos que queremos el valor de P-value sobre el límite superior de la errorbar, e inclinado
  geom_text(aes(label = paste(print.p), y = CI_sup_R2), vjust = -1.5, hjust = 0, angle = 45, cex = 4, parse = T) +
  
  # Especificamos el rango de la escala Y como 1.25 por el máximo de R2, para dejar suficiente espacio para barras y p valores.
  scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
  
  # Especificamos los labels de los ejes X e Y en itálica.
  xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
  ylab(expression(paste("PRS model fit: pseudo-", R ^ 2))) +
  
  # configuramos la barplot que es el tipo de gráfico que queremos. Configuramos en fill para que el color sea funció de su grado de asociación.
  geom_bar(aes(fill = -log10(P_value)), stat = "identity") +
  
  # configuramos la errorbarplot especificando su posiciómn en eje X e Y, así como formato.
  geom_errorbar( aes(x=factor(Threshold), ymin=CI_inf_R2, ymax=CI_sup_R2), width=0.4, colour="black", alpha=0.9, size=1.3) + 
  
  # Especificamos los colores de las barras en funciòn del grado de asociación (según pusimos arriba en fill een barplot). Ponemos un midpoint según el grado de asociación que tengamos. Y la leyenda relativa
  scale_fill_gradient2(
    low = "dodgerblue",
    high = "firebrick",
    mid = "dodgerblue",
    midpoint = 1e-10,
    name = bquote(atop(-log[10] ~ model, italic(P) - value),)) +
  
  # Cambios estéticos del gráfico según queramos
  theme_classic() + theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1))

# Guardamos el gráfico
ggsave("R2_NAG_PLOT.png", height = 16, width = 12)






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
  scale_y_continuous(limits = c(0, max(table2$OR) * 1.6)) +
  
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
  scale_y_continuous(limits = c(0, max(table2$OR) * 1.7)) +
  
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



