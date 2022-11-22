# SCRIPT MEGAPRS 

######### Antes de usar: Descarga los 64 Annotations de:
# http://dougspeed.com/bldldak/


######### Antes de usar: Descarga el archivo highld.txt de:
# https://dougspeed.com/high-ld-regions/



############ PASO 1
### MEGAPRS: Calculate predictor-predictor correlations
# https://dougspeed.com/megaprs/

../ldak5.2.linux --calc-cors pred_pred_cors --bfile PGC_imput_QC_new_filtered --window-kb 3000 --max-threads 28 | tee ../Output_P1.txt

# For each predictor, there are on average 2067.62 other predictors with correlation squared at least 0.001884
# The correlations are saved in files with prefix pred_pred_cors



############ PASO 2
### To calculate taggings assuming the BLD-LDAK Model
# https://dougspeed.com/calculate-taggings/
# https://dougspeed.com/per-predictor-heritabilities/


../ldak5.2.linux --cut-weights sections --bfile PGC_imput_QC_new_filtered --max-threads 28 | tee ../Output_P2.txt

# Thinning complete: 2428322 predictors kept (saved in sections/thin.in), 3263765 lost (sections/thin.out), 42 trivial (sections/thin.trivial)
# If the subsequent cutting of weights fails, or if you wish to re-cut (with different parameters) you can avoid having to thin again by adding "--no-thin DONE"

# The 2428322 predictors have been split into 2438 sections (average length 1223 predictors, longest 1772); details saved in sections/section.details
# Ideally each section would have length less than the number of samples (3519), and this is the case

# Section details saved in sections/section.details and sections/section.number

../ldak5.2.linux --calc-weights-all sections --bfile PGC_imput_QC_new_filtered --max-threads 28 | tee -a ../Output_P2.txt

# Reading weights from 2438 sections
# Merged weights saved in sections/weights.all, with a condensed version in sections/weights.short

mv sections/weights.short bld65 | tee -a ../Output_P2.txt

../ldak5.2.linux --calc-tagging BLD-LDAK --bfile PGC_imput_QC_new_filtered --ignore-weights YES --power -.25  --window-kb 1000 --annotation-number 65 --annotation-prefix bld --save-matrix YES --max-threads 28 | tee -a ../Output_P2.txt

# Taggings saved in BLD-LDAK.tagging, with heritability matrix in BLD-LDAK.matrix



############ PASO 3
### Estimate heritabilities
# https://dougspeed.com/per-predictor-heritabilities/

../ldak5.2.linux --sum-hers bld.ldak --tagfile BLD-LDAK.tagging --summary dataDiscovery_filtered.txt --matrix BLD-LDAK.matrix --max-threads 28 | tee ../Output_P3.txt


# Warning, 14367 (0) of the 5692087 predictors have negative expected heritability (test statistic); this suggests an over-complicated heritability model
# Main results saved in bld.ldak.hers, bld.ldak.cats, bld.ldak.share, bld.ldak.enrich and bld.ldak.extra
# Calculating per-predictor heritabilities for 5692087 predictors
# Estimates of per-predictor heritabilities saved in bld.ldak.ind.hers and bld.ldak.ind.hers.positive



############ PASO 4
### Calculate High-LD Regions:
# https://dougspeed.com/high-ld-regions/
# https://dougspeed.com/gene-based-analysis/

../ldak5.2.linux --cut-genes highld --bfile PGC_imput_QC_new_filtered --genefile highld.txt | tee ../Output_P4.txt

# 24 of the 24 genes were found, spanning 193747 unique predictors and divided into 1 partition; the longest gene contains 37132 predictors
# Details saved in highld/genes.details and highld/genes.predictors.used, with distances of predictors from genes in highld/genes.distances




############ PASO 5
### MEGAPRS: Construct the prediction model
# https://dougspeed.com/megaprs/

../ldak5.2.linux --mega-prs pred_model --model bayesr --ind-hers bld.ldak.ind.hers --summary dataDiscovery_filtered.txt --cors pred_pred_cors --cv-proportion 0.1 --high-LD ./highld/genes.predictors.used --window-kb 1000 --max-threads 28 | tee ../Output_P5_B.txt

# The prediction model is saved in <outfile>.effects.
# This file will have five columns, providing the predictor name, its A1 and A2 alleles, the average number of A1 alleles, then its raw effect (relative to the A1 allele) 
# and is ready to be used for Calculating Scores



############ PASO 6
### Calculate Scores
# https://dougspeed.com/profile-scores/

../ldak5.2.linux --calc-scores SCZ_PRS_SCORE --scorefile pred_model.effects --bfile PGC_imput_QC_new_filtered --power 0 --max-threads 28 | tee ../Output_P6_B.txt


