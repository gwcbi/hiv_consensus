#! /bin/bash

SN="haphpipe_consensus.sh"
read -r -d '' USAGE <<EOF
$SN [sample_dir] [reference_fasta] <adapters_fasta>

Consensus calling from fastq files

EOF

#--- Read command line args
[[ -n "$1" ]] && [[ "$1" == '-h' ]] && echo "$USAGE" && exit 0

[[ -n "$1" ]] && samp="$1"
[[ -n "$2" ]] && ref="$2"
[[ -n "$3" ]] && adapters="$3"

[[ -z ${samp+x} ]] && echo "FAILED: sample_dir is not set" && echo "$USAGE" && exit 0
[[ -z ${ref+x} ]] && echo "FAILED: reference_fasta is not set" && echo "$USAGE" && exit 0

module unload python
module load miniconda3
source activate haphpipe

[[ -e /scratch ]] && export TMPDIR=/scratch
[[ -e /scratch ]] && MAXPROC=$(nproc) || MAXPROC=8

echo "[---$SN---] ($(date)) Starting $SN"

#--- Check that fastq files exist
[[ ! -e "$samp/00_raw/original_1.fastq" && ! -e "$samp/00_raw/original_1.fastq.gz" ]] &&\
 echo "[---$SN---] ($(date)) FAILED: file $samp/00_raw/original_1.fastq does not exist" &&\
 exit 1

gz1=$([[ -e "$samp/00_raw/original_1.fastq.gz" ]] && echo ".gz" || echo "")

[[ ! -e "$samp/00_raw/original_2.fastq" && ! -e "$samp/00_raw/original_2.fastq.gz" ]] &&\
 echo "[---$SN---] ($(date)) FAILED: file $samp/00_raw/original_2.fastq does not exist" &&\
 exit 1

echo "[---$SN---] ($(date)) Sample:    $samp"

gz2=$([[ -e "$samp/00_raw/original_2.fastq.gz" ]] && echo ".gz" || echo "")

#--- Check that reference exists
[[ ! -e "$ref" ]] && echo "[---$SN---] ($(date)) FAILED: reference index $ref does not exist" && exit 1
echo "[---$SN---] ($(date)) Reference: $ref"

#--- Check that GTF exists
# [[ ! -e "$refgtf" ]] && echo "[---$SN---] ($(date)) FAILED: reference gtf $refgtf does not exist" && exit 1
# echo "[---$SN---] ($(date)) Ref GTF:   $refgtf"

#--- Check adapters if provided
if [[ -n "$adapters" ]]; then
    [[ ! -e "$adapters" ]] && echo "[---$SN---] ($(date)) FAILED: adapters $adapters does not exist" && exit 1
    echo "[---$SN---] ($(date)) Adapters:  $adapters"
    aparam="--adapter_file $adapters"
else
    echo "[---$SN---] ($(date)) Adapters:  $adapters"
    aparam=""
fi

echo "[---$SN---] ($(date)) Ad param:  $aparam"

# Get the current directory
prev=$(pwd)

#--- Start the timer
t1=$(date +"%s")


##########################################################################################
# Step 1: Trim reads.
##########################################################################################
stage="trim_reads"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/00_trim

if [[ -e $samp/00_trim/read_1.fq && -e $samp/00_trim/read_2.fq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage read_1.fq,read_2.fq"
else
    substage=""
    read -r -d '' cmd <<EOF
hp_assemble trim_reads --ncpu $MAXPROC\
 $aparam\
 --fq1 $samp/00_raw/original_1.fastq$gz1\
 --fq2 $samp/00_raw/original_2.fastq$gz2\
 --outdir $samp/00_trim
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi
    
    # Symlink to rename files
    ln -fs trimmed_1.fastq $samp/00_trim/read_1.fq
    ln -fs trimmed_2.fastq $samp/00_trim/read_2.fq
fi

##########################################################################################
# Step 2: Align to HXB2 and call consensus
##########################################################################################
stage="align"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/02_align

if [[ -e $samp/02_align/consensus.fasta ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage consensus.fasta"
else
    # Subsample reads - 20 percent, seed=1
    # samtools view -bs 1.2 $samp/00_bless/reads.bam > $samp/07_refine1/sub.bam
    # picard SamToFastq I=$samp/07_refine1/sub.bam F=$samp/07_refine1/sub_1.fq F2=$samp/07_refine1/sub_2.fq
    
    substage="align_reads"
    read -r -d '' cmd <<EOF 
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $ref\
 --fq1 $samp/00_trim/read_1.fq\
 --fq2 $samp/00_trim/read_2.fq\
 --rgid $samp\
 --bt2_preset fast-local\
 --outdir $samp/02_align
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi
    
    substage="call_variants"
    read -r -d '' cmd <<EOF 
hp_assemble call_variants --ncpu $MAXPROC\
 --ref_fa $ref\
 --aln_bam $samp/02_align/aligned.bam\
 --emit_all\
 --outdir $samp/02_align
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi

    module load viral-ngs

    substage="vcf_to_fasta"
    read -r -d '' cmd <<EOF
assembly.py vcf_to_fasta\
 --min_coverage 1 \
  $samp/02_align/variants.vcf.gz\
 $samp/02_align/consensus.fasta
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi

    module unload viral-ngs
fi

# 
# rm -f $samp/07_refine1/sub.bam
# rm -f $samp/07_refine1/sub_1.fq
# rm -f $samp/07_refine1/sub_2.fq
# rm -f $samp/07_refine1/aligned.bam
# rm -f $samp/07_refine1/aligned.bam.bai
# 



##########################################################################################
# Step 3: Align to previous consensus and call refined
##########################################################################################
stage="refine"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/03_refine

if [[ -e $samp/03_refine/consensus.fasta ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage consensus.fasta"
else
    # Subsample reads - 20 percent, seed=1
    # samtools view -bs 1.2 $samp/00_bless/reads.bam > $samp/07_refine1/sub.bam
    # picard SamToFastq I=$samp/07_refine1/sub.bam F=$samp/07_refine1/sub_1.fq F2=$samp/07_refine1/sub_2.fq
    
    substage="align_reads"
    read -r -d '' cmd <<EOF 
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $samp/02_align/consensus.fasta\
 --fq1 $samp/00_trim/read_1.fq\
 --fq2 $samp/00_trim/read_2.fq\
 --rgid $samp\
 --bt2_preset fast-local\
 --outdir $samp/03_refine
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi
    
    substage="call_variants"
    read -r -d '' cmd <<EOF 
hp_assemble call_variants --ncpu $MAXPROC\
 --ref_fa $samp/02_align/consensus.fasta\
 --aln_bam $samp/03_refine/aligned.bam\
 --emit_all\
 --outdir $samp/03_refine
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi

    module load viral-ngs

    substage="vcf_to_fasta"
    read -r -d '' cmd <<EOF
assembly.py vcf_to_fasta\
 --min_coverage 1 \
  $samp/03_refine/variants.vcf.gz\
  $samp/03_refine/consensus.fasta
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi

    module unload viral-ngs

    sed -i "s/^>/>${samp}_consensus /" $samp/03_refine/consensus.fasta
fi

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

