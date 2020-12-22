#!/bin/csh
#$ -M yxu6@nd.edu
#$ -m abe
#$ -pe smp 24
#$ -q debug ##*@@emichaellab
#$ -N STAR_World

module load R

R CMD BATCH STAR_MC_Parallel_1220N.r
