###Script for Using Fishers Exact Test in PM Characterization Study###
##Author: Patrick Woods##
##Date: March 13, 2024##
#Note: the counts data file has four columns: Reference_Resistant, Alternate_Resistant, Reference_Susceptible, and Alternate_Susceptible.#
#Note: the counts data file has 23 rows containing allele count data, one row for each locus to be tested.#

library(readxl)
counts <- read_xlsx(file.choose()) #read in allele counts data
counts <- as.data.frame(counts) #convert allele counts data to data.frame object

loci <- 1:nrow(counts) #this will determine the number of tests

p_vals <- data.frame() #initialize an empty dataframe to store p_values in

#Set up a for loop that creates a 2x2 table for each locus to be tested#
for (i in loci) {
  table <- matrix(c(counts$Reference_Resistant[i],counts$Alternate_Resistant[i],counts$Reference_Susceptible[i],counts$Alternate_Susceptible[i]), nrow = 2)
  ft <- fisher.test(table)
  p_vals <- rbind(p_vals, ft$p.value)
}                

#Determine bonferroni correction threshold

bf_threshold <- 0.05/loci #0.002173913 is the adjusted 5% significance threshold

