# Abstract- Metabolic consequences of sex-reversal in two lizard species: a test of the like genotype and like phenotype hypotheses
Vertebrate sex is typically determined genetically, but in many ectotherms sex can be determined by genes (Genetic Sex Determination: GSD), temperature (Temperature-dependent Sex Determination: TSD), or interactions between genes and temperature during development. Temperature dependent sex determination may involve GSD systems with either male or female heterogamety (XX/XY or ZZ/ZW). In either case, temperature overrides chromosomal sex determination to cause a mismatch between genetic sex and phenotypic sex (sex-reversal). Evolutionary transitions in sex determination can occur rapidly when selection favours the reversed sex over their concordant phenotypic sex or when frequency-dependent selection eliminates heterogametic individuals from populations. To investigate the consequences of sex-reversal on offspring fitness, we measured two energy-driven traits often linked to fitness (metabolism and growth rate) and 6-month survival in two species of reptile with different patterns of temperature-induced sex-reversal. Male sex-reversal occurs in Bassiana duperreyi when chromosomal females (female XX) develop male phenotypes (maleSR XX), while female sex-reversal occurs in Pogona vitticeps when chromosomal males (male ZZ) develop female phenotypes (femaleSR ZZ). We show, for B. duperreyi, metabolism in maleSR XX was like that of male XY, that is, reflective of phenotypic sex and lower than genotypic sex. In contrast, for P. vitticeps, femaleSR ZZ metabolism was intermediate between male ZZ and female ZW metabolic rate. While we find weak evidence that sex-reversal conveys a fitness advantage in either species, energetic processes may still constrain the distribution of sex-reversal in nature.

# Metadata
Folder "final.analysis.data" you will find data used for all analysis and figures (metabolism, surivial & growth). This includes two csv files of O2 data for both species (Bassiana.finalO2.sexreversal.analysis.data.clean.csv & Pogona.finalO2.sexreversal.analysis.data.clean.csv). There are also x2 RDS files that are used for standard deviation prediction curves found in Figure 2 (B/D). Growth and data for both species is found csv file "growth.bassiana.pogona.csv"


# Scripts
The script "sexreversal.analysis.R" has all the code for the BRMS models for each species and the trait being measured (metabolism, growth), with their appropriate checks. Also there you'll find analysis for survivorship and growth in this script. BRMS models have been saved in the folder "models" which was then used in the .rmd file "Results.Rmd". All figures in scripts were saved in folder "Final.Figures". Headings will provide detail of code that follows

The script "Supp_SMR_02_Bassiana_Pog.R" has all the code for the BRMS models for SMR for each species , with their appropriate checks. BRMS models have been saved in the folder "SMR_models" within the "models" folder, which was then used in the .rmd file "Supplementary.Results.Rmd". All supp figures in scripts were saved in folder "Final.Figures". Headings will provide detail of code that follows


# Rmd file
rmd file "Results.Rmd" has results, and tables as appeared in the text. Code chunks detail analysis that is being conducted for each section. 

rmd file "Supplementary.Results.Rmd" has supplementary results that include nest data and SMR analysis, and tables. Code chunks detail analysis for this section. 

