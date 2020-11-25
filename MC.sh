#!/bin/csh
#$ -M yxu6@nd.edu
#$ -m abe
#$ -pe smp 4
#$ -q long ##*@@emichaellab
#$ -N STAR_World 

module load R

R CMD BATCH STAR_MC_Parallel.r
