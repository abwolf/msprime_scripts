#!/bin/bash
#SBATCH --get-user-env
#SBATCH --mem=5G
#SBATCH --qos=1hr
#SBATCH --time=1:00:00

source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

mdl=$( echo Sriram )
seed=$( cat ../Sriram.chr_list | awk 'BEGIN {OFS="\t"} NR==(2000 + '$SLURM_ARRAY_TASK_ID' ) {print $0}' )
#seed=$seed

n1=$( echo 0.0)
n2=$( echo 0.0)

m_AF_B=$( echo 0.0 )
m_B_AF=$( echo 0.0 )
m_AF_EU=$( echo 0.0 )
m_EU_AF=$( echo 0.0 )

eur=$( echo 1006 )   #1006
asn=$( echo 1008 )   #1008 ; 2040
len=$( awk 'BEGIN {print 1e6}' )
pop=$( echo modHum )

tag=$(echo "$mdl"_"$pop"_"$seed"_n1_"$n1"_mAfB_"$m_AF_B"_mBAf_"$m_B_AF"_mAfEu_"$m_AF_EU"_mEuAf_"$m_EU_AF")
dir=~/SimulatedDemographic/msprime/
sstardir=~/SimulatedDemographic/Sstar/


echo mdl: $mdl	seed: $seed n1: $n1 n2: $n2 eur: $eur asn: $asn		tag: $tag
echo **RUN MSPRIME SIMULATION AND OUTPUT VCF**
cmd=$( echo " python $dir/bin/msprime.Admixture_Simulate.py \n
			-p $pop \n
			-o $mdl \n
			-s $seed \n
			-i 2 \n
			-n $n1 \n
			-d $n2 \n
			-e $eur \n
			-a $asn \n
			-r 216 \n
			--migration_AF_B $m_AF_B \n
			--migration_B_AF $m_B_AF \n
			--migration_AF_AS 0 \n
			--migration_AS_AF 0 \n
			--migration_AF_EU $m_AF_EU \n
			--migration_EU_AF $m_EU_AF \n
			--migration_EU_AS 0 \n
			--migration_AS_EU 0 \n
			-c ILS \n
			-l $len  \n
			| gzip -c - > $tag.bed.gz")
echo -e $cmd

eval $( echo -e $cmd )

echo fin simulation
date
