#!/bin/bash

#SBATCH --job-name=tRNAscan
#SBATCH -D /nfs3_ib/nfs-ip34/home/def-labolcf/programs/labolcf/my_phastest
#SBATCH -o /nfs3_ib/nfs-ip34/home/def-labolcf/programs/labolcf/my_phastest/JOBS/ERR017368assembly/tmp/tRNAscan/log/slurm-%A_%a.out
#SBATCH --account=def-labolcf
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G

echo "load env"
module load gcc/9.3.0 trnascan-se/2.0.12 fraggenescan/1.31 aragorn/1.2.41 barrnap/0.9 blast+/2.13.0 prodigal/2.6.3 mugqic/ucsc/v387
echo "get fasta"
export FA_IN=$(ls "/nfs3_ib/nfs-ip34/home/def-labolcf/programs/labolcf/my_phastest/JOBS/${JOB_ID}/tmp/tRNAscan/*.fa" | awk "NR==$SLURM_ARRAY_TASK_ID")
echo "fasta: $FA_IN"
b=$(basename $FA_IN)

cp $FA_IN $SLURM_TMPDIR/

export OUT_DIR=/nfs3_ib/nfs-ip34/home/def-labolcf/programs/labolcf/my_phastest/JOBS/${JOB_ID}/tmp/tRNAscan/out
echo "running tRNAscan-SE"
tRNAscan-SE -B -o $SLURM_TMPDIR/$b.out $SLURM_TMPDIR/$b --thread 12
echo "copy result back to $OUT_DIR"
cp $SLURM_TMPDIR/$b.out $OUT_DIR/
echo "done!"
