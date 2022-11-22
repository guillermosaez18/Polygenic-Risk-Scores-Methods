# Estad?sticas de comparaci?n de resultados entre m?todos de PRS

library(dplyr)
library(nortest)


##### ESTAD?STICAS CL?SICO

setwd("C:/Users/guillermo.saez/OneDrive - BABEL/Escritorio/Personal/TFG/Escrito/Stats/Clásico")

data_perm <- read.table('data_SCZ_lrm_perm.txt', header=TRUE, stringsAsFactors = FALSE)
data_perm$R2 <- scan(text=data_perm$R2, dec=",", sep=".")
data_perm$R2_liab <- scan(text=data_perm$R2_liab, dec=",", sep=".")

# Coger n = 84 elementos aleatorios de cada Q_Range para que el método clásico tenga el mismo número 
# de valores que los demás métodos (1000) pero obteniendo valores de todos los rangos
data_classic <- data_perm %>% group_by(Range) %>% sample_n(size = 84)

# De los 84*12=1008 elementos, cogemos 1000 aleatorios
data_classic <- data_classic[sample(nrow(data_classic),1000),]

# visualiza cuántos elementos hay por cada Q_Range
data_classic %>% count(Range)

data_auc <- as.data.frame(data_classic$AUC)
data_R2 <- as.data.frame(data_classic$R2)
data_R2liab <- as.data.frame(data_classic$R2_liab)


##### ESTAD?STICAS LDPRED2

setwd("C:/Users/guillermo.saez/OneDrive - BABEL/Escritorio/Personal/TFG/Escrito/Stats/LDpred2")

data_perm <- read.table('data_SCZ_lrm_perm.txt', header=TRUE, stringsAsFactors = FALSE)
data_perm$R2 <- scan(text=data_perm$R2, dec=",", sep=".")
data_perm$R2_liab <- scan(text=data_perm$R2_liab, dec=",", sep=".")

data_auc$LDpred2 <- as.data.frame(data_perm$AUC)
data_R2$LDpred2 <- as.data.frame(data_perm$R2)
data_R2liab$LDpred2 <- as.data.frame(data_perm$R2_liab)


##### ESTAD?STICAS LASSOSUM

setwd("C:/Users/guillermo.saez/OneDrive - BABEL/Escritorio/Personal/TFG/Escrito/Stats/LASSOSUM")

data_perm <- read.table('data_SCZ_lrm_perm.txt', header=TRUE, stringsAsFactors = FALSE)
data_perm$R2 <- scan(text=data_perm$R2, dec=",", sep=".")
data_perm$R2_liab <- scan(text=data_perm$R2_liab, dec=",", sep=".")

data_auc$Lassosum <- as.data.frame(data_perm$AUC)
data_R2$Lassosum <- as.data.frame(data_perm$R2)
data_R2liab$Lassosum <- as.data.frame(data_perm$R2_liab)


##### ESTAD?STICAS PRS-CS

setwd("C:/Users/guillermo.saez/OneDrive - BABEL/Escritorio/Personal/TFG/Escrito/Stats/PRS-CS")

data_perm <- read.table('data_SCZ_lrm_perm.txt', header=TRUE, stringsAsFactors = FALSE)
data_perm$R2 <- scan(text=data_perm$R2, dec=",", sep=".")
data_perm$R2_liab <- scan(text=data_perm$R2_liab, dec=",", sep=".")

data_auc$PRS_CS <- as.data.frame(data_perm$AUC)
data_R2$PRS_CS <- as.data.frame(data_perm$R2)
data_R2liab$PRS_CS <- as.data.frame(data_perm$R2_liab)


##### ESTAD?STICAS MEGAPRS

setwd("C:/Users/guillermo.saez/OneDrive - BABEL/Escritorio/Personal/TFG/Escrito/Stats/MEGAPRS")

data_perm <- read.table('data_SCZ_lrm_perm.txt', header=TRUE, stringsAsFactors = FALSE)
data_perm$R2 <- scan(text=data_perm$R2, dec=",", sep=".")
data_perm$R2_liab <- scan(text=data_perm$R2_liab, dec=",", sep=".")

data_auc$MEGAPRS <- as.data.frame(data_perm$AUC)
data_R2$MEGAPRS <- as.data.frame(data_perm$R2)
data_R2liab$MEGAPRS <- as.data.frame(data_perm$R2_liab)




data_auc$LDpred2 <- unlist(data_auc$LDpred2)
data_auc$MEGAPRS <- unlist(data_auc$MEGAPRS)
data_auc$PRS_CS <- unlist(data_auc$PRS_CS)
data_auc$Lassosum <- unlist(data_auc$Lassosum)

data_R2$LDpred2 <- unlist(data_R2$LDpred2)
data_R2$MEGAPRS <- unlist(data_R2$MEGAPRS)
data_R2$PRS_CS <- unlist(data_R2$PRS_CS)
data_R2$Lassosum <- unlist(data_R2$Lassosum)

data_R2liab$LDpred2 <- unlist(data_R2liab$LDpred2)
data_R2liab$MEGAPRS <- unlist(data_R2liab$MEGAPRS)
data_R2liab$PRS_CS <- unlist(data_R2liab$PRS_CS)
data_R2liab$Lassosum <- unlist(data_R2liab$Lassosum)


colnames(data_auc) <- c("Classic","LDpred2","Lassosum","PRS_CS","MEGAPRS")
colnames(data_R2) <- c("Classic","LDpred2","Lassosum","PRS_CS","MEGAPRS")
colnames(data_R2liab) <- c("Classic","LDpred2","Lassosum","PRS_CS","MEGAPRS")




# ANÁLISIS ESTADÍSTICO

###### AUC 

# An?lisis de normalidad (Normal si p-valor > 0.05)

lillie.test(data_auc$Classic)[2]  # 1.097323e-67
lillie.test(data_auc$LDpred2)[2]  # 0.5153529
lillie.test(data_auc$MEGAPRS)[2]  # 0.3653507
lillie.test(data_auc$PRS_CS)[2]   # 0.3210209
lillie.test(data_auc$Lassosum)[2] # 0.9547421


# Resultado: Las distribuciones de AUC de todos los m?todos excepto el Cl?sico son normales, 
# por lo que debemos usar el test de Kruskal-Wallis


# reformat
AUC <- data.frame(
  Method = c(rep("Classic", 1000), rep("LDpred2", 1000), rep("Lassosum", 1000), rep("PRS-CS", 1000), rep("MEGAPRS", 1000)),
  values = c(data_auc$Classic,data_auc$LDpred2,data_auc$Lassosum,data_auc$PRS_CS,data_auc$MEGAPRS))


kruskal.test(values ~ Method, data = AUC) # Kruskal-Wallis chi-squared = 3012.3, df = 4, p-value < 2.2e-16

# p-valor < 0.05 significa que existen diferencias significativas entre grupos
# hay que aplicar Mann-Whitney  



wilcox.test(data_auc$LDpred2,data_auc$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999968, p-value < 2.2e-16
# 95 percent confidence interval: 0.06258434 0.06633565
# difference in location 0.06438861



wilcox.test(data_auc$Lassosum,data_auc$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 998953, p-value < 2.2e-16
# 95 percent confidence interval: 0.05114255 0.05491598
# difference in location  0.0529517



wilcox.test(data_auc$PRS_CS,data_auc$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999813, p-value < 2.2e-16
# 95 percent confidence interval:0.05688552 0.06060394
# difference in location  0.05868113 



wilcox.test(data_auc$MEGAPRS,data_auc$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999968, p-value < 2.2e-16
# 95 percent confidence interval: 0.06290833 0.06663152
# difference in location  0.06471818 



###### R2

# An?lisis de normalidad (Normal si p-valor > 0.05)

lillie.test(data_R2$Classic)[2]  # 5.283723e-59
lillie.test(data_R2$LDpred2)[2]  # 0.02142314
lillie.test(data_R2$MEGAPRS)[2]  # 0.2892598
lillie.test(data_R2$PRS_CS)[2]   # 0.5555235
lillie.test(data_R2$Lassosum)[2] # 0.0427544

# Resultado: Las distribuciones de R2 son normales para MEGAPRS, PRS-CS, 
# pero no para el m?todo Cl?sico, LDpred2 y Lassosum, por lo que debemos usar el test de Kruskal-Wallis


# reformat
R2 <- data.frame(
  Method = c(rep("Classic", 1000), rep("LDpred2", 1000), rep("Lassosum", 1000), rep("PRS-CS", 1000), rep("MEGAPRS", 1000)),
  values = c(data_R2$Classic,data_R2$LDpred2,data_R2$Lassosum,data_R2$PRS_CS,data_R2$MEGAPRS))


kruskal.test(values ~ Method, data = R2) # Kruskal-Wallis chi-squared = 2980.2, df = 4, p-value < 2.2e-16


# p-valor < 0.05 significa que existen diferencias significativas entre grupos
# hay que aplicar Mann-Whitney  



wilcox.test(data_R2$LDpred2,data_R2$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999819, p-value < 2.2e-16
# 95 percent confidence interval: 8.841734 9.439391
# difference in location 9.130082



wilcox.test(data_R2$Lassosum,data_R2$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999930, p-value < 2.2e-16
# 95 percent confidence interval: 9.930251 10.529780
# difference in location  10.21891 



wilcox.test(data_R2$PRS_CS,data_R2$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 998773, p-value < 2.2e-16
# 95 percent confidence interval: 7.796518 8.374624
# difference in location  8.074971  



wilcox.test(data_R2$MEGAPRS,data_R2$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999915, p-value < 2.2e-16
# 95 percent confidence interval: 9.677237 10.268495
# difference in location  9.961408  




###### R2 liability

# An?lisis de normalidad (Normal si p-valor > 0.05)

lillie.test(data_R2liab$Classic)[2]  # 4.676183e-54
lillie.test(data_R2liab$LDpred2)[2]  # 0.005909314
lillie.test(data_R2liab$MEGAPRS)[2]  # 0.3477516
lillie.test(data_R2liab$PRS_CS)[2]   # 0.7856797
lillie.test(data_R2liab$Lassosum)[2] # 0.01302962 

# Resultado: Las distribuciones de R2 liability son normales para MEGAPRS, PRS-CS, 
# pero no para el m?todo Cl?sico, LDpred2 y Lassosum, por lo que debemos usar el test de Kruskal-Wallis


# reformat
R2liab <- data.frame(
  Method = c(rep("Classic", 1000), rep("LDpred2", 1000), rep("Lassosum", 1000), rep("PRS-CS", 1000), rep("MEGAPRS", 1000)),
  values = c(data_R2liab$Classic,data_R2liab$LDpred2,data_R2liab$Lassosum,data_R2liab$PRS_CS,data_R2liab$MEGAPRS))


kruskal.test(values ~ Method, data = R2liab) # Kruskal-Wallis chi-squared = 2980.2, df = 4, p-value < 2.2e-16


# p-valor < 0.05 significa que existen diferencias significativas entre grupos
# hay que aplicar Mann-Whitney  



wilcox.test(data_R2liab$LDpred2,data_R2liab$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999819, p-value < 2.2e-16
# 95 percent confidence interval: 4.501511 4.797526
# difference in location 4.644786 



wilcox.test(data_R2liab$Lassosum,data_R2liab$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999930, p-value < 2.2e-16
# 95 percent confidence interval: 5.089232 5.384689
# difference in location  5.232398  



wilcox.test(data_R2liab$PRS_CS,data_R2liab$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 998773, p-value < 2.2e-16
# 95 percent confidence interval: 3.945652 4.229042
# difference in location  4.082035   



wilcox.test(data_R2liab$MEGAPRS,data_R2liab$Classic, alternative = "two.sided",
            paired = FALSE, conf.int = 0.95)

# W = 999915, p-value < 2.2e-16
# 95 percent confidence interval: 4.951336 5.243177
# difference in location  5.092925   



##### Kruskal Wallis sólo nuevos métodos

AUC <- data.frame(
  Method = c(rep("LDpred2", 1000), rep("Lassosum", 1000), rep("PRS-CS", 1000), rep("MEGAPRS", 1000)),
  values = c(data_auc$LDpred2,data_auc$Lassosum,data_auc$PRS_CS,data_auc$MEGAPRS))


kruskal.test(values ~ Method, data = AUC)


R2 <- data.frame(
  Method = c(rep("LDpred2", 1000), rep("Lassosum", 1000), rep("PRS-CS", 1000), rep("MEGAPRS", 1000)),
  values = c(data_R2$LDpred2,data_R2$Lassosum,data_R2$PRS_CS,data_R2$MEGAPRS))


kruskal.test(values ~ Method, data = R2)


R2liab <- data.frame(
  Method = c(rep("LDpred2", 1000), rep("Lassosum", 1000), rep("PRS-CS", 1000), rep("MEGAPRS", 1000)),
  values = c(data_R2liab$LDpred2,data_R2liab$Lassosum,data_R2liab$PRS_CS,data_R2liab$MEGAPRS))


kruskal.test(values ~ Method, data = R2liab)



library(ggplot2)

ggplot(data = AUC, mapping = aes(x = Method, y = values, colour = Method)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Methods")+
  ylab("AUC")


ggplot(data = R2, mapping = aes(x = Method, y = values, colour = Method)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Methods")+
  ylab("R2")


ggplot(data = R2_liab, mapping = aes(x = Method, y = values, colour = Method)) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Methods")+
  ylab("R2  liability scale")

