#!/bin/bash

#SBATCH --mem=1G
#SBATCH --time 0-1
#SBATCH --error=%x_%A_%a.out
#SBATCH --output=%x_%A_%a.out

set -euo pipefail
date

module load anaconda3
PS1=""
conda activate msprime_scripts

readarray -t SWEEP_ARR <<<"$SWEEP_OPTIONS"

START=$(($SLURM_ARRAY_TASK_ID * $BATCH_SIZE))
END=$(($START + $BATCH_SIZE))
END=$(($END < ${#SWEEP_ARR[@]} ? $END : ${#SWEEP_ARR[@]}))

for (( ind=$START; ind<$END; ind++)); do

    echo $ind
    fname=$(printf "%s" "${SWEEP_ARR[$ind]}" | tr -d ' ')
    time python $MSPRIME_LOC \
        $BASE_OPTIONS \
        ${SWEEP_ARR[$ind]} \
        -s ${SEED:-$RANDOM} \
        --haplo "haplo$fname.gz"

done
date
