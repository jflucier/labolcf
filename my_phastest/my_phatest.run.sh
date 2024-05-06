#!/bin/bash

module load gcc/9.3.0 trnascan-se/2.0.12 fraggenescan/1.31 aragorn/1.2.41 barrnap/0.9 blast+/2.13.0 prodigal/2.6.3 mugqic/ucsc/v387

export FA_IN="/home/def-labolcf/programs/test_my_phastest/62_CNEO.1-Contigs.fasta"
export job_id="62_CNEO.1"

export PHASTEST_HOME=/home/def-labolcf/programs/labolcf/my_phastest
export jobs_dir="$PHASTEST_HOME/JOBS";
export scripts_dir="$PHASTEST_HOME/scripts";  # dir of executables.
export database_dir="$PHASTEST_HOME/DB";
export PATH=${PHASTEST_HOME}../bin:$PATH

## options
export single_diamond=1
export complete_annotation=1
# lite
export anno_mode="-l"
#fasta input
export flag="-s"

## Blast options
export use_split_db_vir=0; # do not split up the viral DB
export target_query_pieces_vir=208;
export threads=24;
export virus_database="prophage_virus.db"; # virus db name
export virus_header_database="prophage_virus_header_lines.db"
export virus_database_path="${database_dir}/${virus_database}"
export virus_header_database_path="${database_dir}/${virus_header_database}"
export bac_database="swissprot.db"
export bac_header_database="swissprot_header.db";
export bac_database_path="${database_dir}/${bac_database}"
export bac_header_database="${database_dir}/${bac_header_database}"

rm -fr ${PHASTEST_HOME}/JOBS/$job_id
mkdir -p $jobs_dir/$job_id
cp $FA_IN ${PHASTEST_HOME}/JOBS/$job_id/$job_id.fna

if [ -f "$jobs_dir/$job_id/$job_id.fasta" ]; then
    mv $jobs_dir/$job_id/$job_id.fasta $jobs_dir/$job_id/$job_id.fna
fi

if [ "$flag" = "-c" ]; then
    mv $jobs_dir/$job_id/$job_id.fna $jobs_dir/$job_id/${job_id}_original.fna
fi

cd $jobs_dir/$job_id/

echo "Handle fna file ...";
if [ "$flag" = "-c" ]; then
    perl $scripts_dir/make_contig_position.pl ${job_id}_original.fna
    echo 'Regenerate fna file....'
    perl $scripts_dir/find_common_element.pl $job_id.fna $job_id
    cp ${job_id}_original.fna ${job_id}_filtered.fna
    perl $scripts_dir/fix_fna_lines.pl ${job_id}_filtered.fna
    echo "Running FragGeneScan on ${job_id}.fna"
    fraggenescan_dir="$jobs_dir/$num/tmp/fraggenescan";
    #bash $scripts_dir/call_fraggenescan_parallel.sh $fraggenescan_dir ${job_id}
    perl $scripts_dir/call_fraggenescan.pl ${job_id}.fna
fi

if [ "$flag" = "-s" ]; then
    mv ${job_id}.fna ${job_id}_original.fna
    perl -ne '
    chomp($_);
    if($_ =~ /\>/){
      my $h = substr($_,1);
      print ">gi|00000000|ref|NC_000000| $h\n";
    }
    else{
      print $_ . "\n";
    }
    ' ${job_id}_original.fna > ${job_id}.fna
    perl $scripts_dir/fix_fna_lines.pl ${job_id}.fna

    echo "Running Prodigal on ${job_id}.fna"
    prodigal -i ${job_id}.fna -o ${job_id}.gff -f gff
    perl $scripts_dir/format_prodigal.pl ${job_id}
fi

echo "Generating ptt file ..."
perl $scripts_dir/change_to_ptt_format.pl ${job_id}.fna ${job_id}.predict ${job_id}.ptt

echo "Generating faa file..."
perl $scripts_dir/change_to_protein_seq.pl ${job_id}.fna ${job_id}.predict ${job_id}.faa

echo "Running phage search ..."
echo "Performing BLAST search on virus db."
export pepfile="$job_id.faa"
export infofile="$job_id.ptt"

export blast_v_dir="$PHASTEST_HOME/JOBS/$job_id/tmp/blast_v"
mkdir -p $blast_v_dir
cp $PWD/$pepfile $blast_v_dir

#cd $blast_v_dir
export blast_out_fmt="6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore"; # tabular
blastp \
-db $virus_database_path \
-outfmt "${blast_out_fmt}" \
-evalue 0.0001 \
-query $blast_v_dir/$pepfile \
-out $blast_v_dir/${pepfile}_out \
-num_threads $threads \
-seg no

cat $blast_v_dir/${pepfile}_out > $blast_v_dir/${pepfile}_blast_out
cp $blast_v_dir/${pepfile}_blast_out $PWD/ncbi.out

### not sure this is needed
# split tRNAscan input for faster processing
export tRNAscan_dir="$PHASTEST_HOME/JOBS/$job_id/tmp/tRNAscan"
mkdir -p $tRNAscan_dir/out
mkdir -p $PHASTEST_HOME/JOBS/log
faSplit sequence $job_id.fna 200 "${tRNAscan_dir}/$job_id"

# submit and wait
#export tmp_jobid="62_CNEO"
rm $PHASTEST_HOME/JOBS/log/*
sbatch --array=1-196 --export=JOB_ID="$job_id" /nfs3_ib/nfs-ip34/home/def-labolcf/programs/labolcf/my_phastest/tRNAscan_task.sh

#grep -Le 'done!' tmp/tRNAscan/log/*
#wc -l tmp/tRNAscan/ERR017368assembly121.fa
#wc -l tmp/tRNAscan/ERR017368assembly142.fa
#wc -l tmp/tRNAscan/ERR017368assembly148.fa
#wc -l tmp/tRNAscan/ERR017368assembly186.fa
#
#grep -e '>' tmp/tRNAscan/ERR017368assembly121.fa | wc -l
#grep -e '>' tmp/tRNAscan/ERR017368assembly142.fa | wc -l
#grep -e '>' tmp/tRNAscan/ERR017368assembly148.fa | wc -l
#grep -e '>' tmp/tRNAscan/ERR017368assembly186.fa | wc -l

### not sure this is needed
#echo "find tRNA sequences using tRNAscan..."
#tRNAscan-SE -B -o tRNAscan.out $job_id.fna --thread $threads

echo "find tRNA sequences using aragorn..."
aragorn -v -m -o tmRNA_aragorn.out $job_id.fna
echo "find rRNA sequences using barrnap..."
barrnap --threads $threads --outseq rRNA_barrnap.out $job_id.fna
echo "extract rRNA from results..."
perl $scripts_dir/extract_RNA.pl $job_id extract_RNA_result.txt.tmp $flag
echo "running make_RNA_png_input"
perl $scripts_dir/make_RNA_png_input.pl extract_RNA_result.txt.tmp RNA_output.out

export NC="NC_000000";
export gi='N/A';
mkdir -p $jobs_dir/$job_id/${NC}_dir
cd $jobs_dir/$job_id/${NC}_dir

echo "Scanning for phage regions ..."
if [[ -s "$jobs_dir/$job_id/${job_id}_contig_positions.txt" ]]
then
  perl $scripts_dir/scan.pl \
  -n $jobs_dir/$job_id/$job_id.fna \
  -a $jobs_dir/$job_id/$job_id.faa  \
  -t $jobs_dir/$job_id/RNA_output.out \
  -m $jobs_dir/$job_id/tmRNA_aragorn.out \
  -b $jobs_dir/$job_id/ncbi.out \
  -p $jobs_dir/$job_id/$job_id.ptt \
  -use 5 -c $jobs_dir/$job_id/${job_id}_contig_positions.txt >${NC}_phmedio.txt
else
  perl $scripts_dir/scan.pl \
  -n $jobs_dir/$job_id/$job_id.fna \
  -a $jobs_dir/$job_id/$job_id.faa  \
  -t $jobs_dir/$job_id/RNA_output.out \
  -m $jobs_dir/$job_id/tmRNA_aragorn.out \
  -b $jobs_dir/$job_id/ncbi.out \
  -p $jobs_dir/$job_id/$job_id.ptt \
  -use 5 >${NC}_phmedio.txt
fi

# 0 if prophage region only, 1 for annotating all proteins in genome.
export anno_flag=0
perl $scripts_dir/non_hit_region_pro_to_faa.pl \
$jobs_dir/$job_id/$job_id.faa \
${NC}_phmedio.txt \
$jobs_dir/$job_id/$job_id.faa.non_hit_pro_region \
$anno_flag

# skip bacterial blast: see line 774 - 842 if required in labolcf/my_phastest/phastest.jfl.pl
echo "Annotating proteins in regions ..."
perl $scripts_dir/annotation.pl \
$NC $job_id \
$virus_header_database_path $bac_database_path \
$jobs_dir/$job_id/ncbi.out.non_hit_pro_region \
${NC}_phmedio.txt $jobs_dir/$job_id/RNA_output.out $jobs_dir/$job_id/$job_id.ptt $flag $anno_flag

echo "Extracting proteins in regions ..."
perl $scripts_dir/extract_protein.pl $job_id ${NC}_phmedio.txt $job_id.ptt extract_result.txt $anno_flag

echo "Get true regions ..."
perl $scripts_dir/get_true_region.pl ${NC}_phmedio.txt extract_result.txt true_defective_prophage.txt
cp true_defective_prophage.txt  $jobs_dir/$job_id/
cp true_defective_prophage.txt  $jobs_dir/$job_id/summary.txt

cp extract_result.txt $jobs_dir/$job_id/detail.txt

perl $scripts_dir/make_json.pl extract_result.txt  true_defective_prophage.txt $job_id $flag ../$job_id.fna
cp json_input_regions $jobs_dir/$job_id/predicted_phage_regions.json
cp json_input $jobs_dir/$job_id/predicted_genes.json

cd $jobs_dir/$job_id
perl $scripts_dir/make_region_DNA.pl $jobs_dir $job_id

cd ${jobs_dir}/${job_id}
rm -rf \
${NC}_dir \
extract_RNA_result.tmp \
extract_RNA_result.txt.tmp \
ncbi.out \
tmRNA_aragorn.out \
tRNAscan* \
true_defective_prophage.txt \
$job_id.faa \
$job_id.faa.non_hit_pro_region \
${NC}_phmedio.txt_bk \
tmp \
$job_id.predict \
$job_id.ptt \
$job_id.fna.fai \
ncbi.out \
rRNA_barrnap.out


#ls *predict
#ls *process
#ls *ptt
#ls *gbk
#ls *.txt.old
#ls image.png
#ls *f