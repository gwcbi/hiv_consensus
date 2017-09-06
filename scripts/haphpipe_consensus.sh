#! /bin/bash

SN="haphpipe_consensus.sh"
read -r -d '' USAGE <<EOF
$SN [sample_dir] [reference_fasta] <adapters_fasta>

Consensus calling from fastq files

EOF

##########################################################################################
# Initialize: Command line args, input files, environment variables
##########################################################################################
#--- Read command line args
[[ -n "$1" ]] && [[ "$1" == '-h' ]] && echo "$USAGE" && exit 0

[[ -n "$1" ]] && samp="$1"
[[ -n "$2" ]] && ref="$2"
[[ -n "$3" ]] && adapters="$3"

[[ -z ${samp+x} ]] && echo "FAILED: samp is not set" && echo "$USAGE" && exit 0
[[ -z ${ref+x} ]] && echo "FAILED: reference_fasta is not set" && echo "$USAGE" && exit 0

module unload python
module load miniconda3
source activate haphpipe

# Set environment variables for java and TMPDIR
export _JAVA_OPTIONS="-Xms1g -Xmx32g"
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
mkdir -p $samp/01_trim

if [[ -e $samp/01_trim/read_1.fq && -e $samp/01_trim/read_2.fq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage read_1.fq,read_2.fq"
else
    substage=""
    read -r -d '' cmd <<EOF
hp_assemble trim_reads --ncpu $MAXPROC\
 $aparam\
 --fq1 $samp/00_raw/original_1.fastq$gz1\
 --fq2 $samp/00_raw/original_2.fastq$gz2\
 --outdir $samp/01_trim
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
    ln -fs trimmed_1.fastq $samp/01_trim/read_1.fq
    ln -fs trimmed_2.fastq $samp/01_trim/read_2.fq
    
fi


##########################################################################################
# Step 1b: Convert FASTQ to unaligned BAM file
##########################################################################################
if [[ -e $samp/01_trim/reads.bam ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage reads.bam"
else
    substage="fq_to_bam"
    read -r -d '' cmd <<EOF
picard FastqToSam SM=$samp RG=$samp\
 F1=$samp/01_trim/read_1.fq F2=$samp/01_trim/read_2.fq\
 O=$samp/01_trim/reads.bam
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi
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
    samtools view -bs 1.2 $samp/01_trim/reads.bam > $samp/02_align/sub.bam
    picard SamToFastq I=$samp/02_align/sub.bam F=$samp/02_align/sub_1.fq F2=$samp/02_align/sub_2.fq
    
    substage="align_reads"
    read -r -d '' cmd <<EOF 
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $ref\
 --fq1 $samp/02_align/sub_1.fq\
 --fq2 $samp/02_align/sub_2.fq\
 --rgid $samp\
 --bt2_preset fast-local\
 --no_realign\
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

rm -f $samp/02_align/sub.bam
rm -f $samp/02_align/sub_1.fq
rm -f $samp/02_align/sub_2.fq
# rm -f $samp/02_align/aligned.bam
# rm -f $samp/02_align/aligned.bam.bai


##########################################################################################
# Step 3: Align to previous consensus and call refined
##########################################################################################
stage="refine"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/03_refine

if [[ -e $samp/03_refine/consensus.fasta ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage consensus.fasta"
else    
    substage="align_reads"
    read -r -d '' cmd <<EOF 
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $samp/02_align/consensus.fasta\
 --fq1 $samp/01_trim/read_1.fq\
 --fq2 $samp/01_trim/read_2.fq\
 --rgid $samp\
 --bt2_preset sensitive-local\
 --no_realign\
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

    sed -i "s/^>/>${samp}./" $samp/03_refine/consensus.fasta
fi

##########################################################################################
# Step 4a: Fix consensus - align_reads
##########################################################################################
stage="fix_consensus"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/04_fixed
cp $samp/03_refine/consensus.fasta $samp/04_fixed/consensus.fa

if [[ -e $samp/04_fixed/final.bam ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage final.bam"
else
    substage="align_reads"
    read -r -d '' cmd <<EOF
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $samp/04_fixed/consensus.fa\
 --fq1 $samp/01_trim/read_1.fq\
 --fq2 $samp/01_trim/read_2.fq\
 --rgid $samp\
 --bt2_preset very-sensitive-local\
 --markdup\
 --outdir $samp/04_fixed
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi

    # Create filtered alignment including DUPLICATES
    # Keep 2 (PROPER_PAIR)
    # Remove 2828 (UNMAP,MUNMAP,SECONDARY,QCFAIL,SUPPLEMENTARY)
    samtools view -b -f 2 -F 2828 $samp/04_fixed/aligned.bam > $samp/04_fixed/final.bam
    samtools index $samp/04_fixed/final.bam

    # Create final alignment excluding DUPLICATES
    # Keep 2 (PROPER_PAIR)
    # Remove 3852 (UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY)
    samtools view -b -f 2 -F 3852 $samp/04_fixed/final.bam > $samp/04_fixed/final_nodup.bam
    samtools index $samp/04_fixed/final_nodup.bam

fi

##########################################################################################
# Step 4b: Fix consensus - call_variants
##########################################################################################
if [[ -e $samp/04_fixed/final.vcf.gz ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage final.vcf.gz"
else
    # Call variants
    substage="call_variants"
    read -r -d '' cmd <<EOF 
hp_assemble call_variants --ncpu $MAXPROC\
 --ref_fa $samp/04_fixed/consensus.fa\
 --aln_bam $samp/04_fixed/final.bam\
 --outdir $samp/04_fixed
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi
    
    mv $samp/04_fixed/variants.vcf.gz $samp/04_fixed/final.vcf.gz
    mv $samp/04_fixed/variants.vcf.gz.tbi $samp/04_fixed/final.vcf.gz.tbi

    substage="call_variants_nodup"
    read -r -d '' cmd <<EOF 
hp_assemble call_variants --ncpu $MAXPROC\
 --ref_fa $samp/04_fixed/consensus.fa\
 --aln_bam $samp/04_fixed/final_nodup.bam\
 --outdir $samp/04_fixed
EOF

    echo -e "[---$SN---] ($(date)) $stage $substage command:\n\n$cmd\n"
    eval $cmd

    if [[ $? -eq 0 ]]; then
        echo "[---$SN---] ($(date)) COMPLETED: $stage $substage"
    else
        echo "[---$SN---] ($(date)) FAILED: $stage $substage"
        exit 1
    fi
    
    mv $samp/04_fixed/variants.vcf.gz $samp/04_fixed/final_nodup.vcf.gz
    mv $samp/04_fixed/variants.vcf.gz.tbi $samp/04_fixed/final_nodup.vcf.gz.tbi
    
fi


##########################################################################################
# Step 4c: Generate summary
##########################################################################################
if [[ -e $samp/04_fixed/summary.txt ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage summary.txt"
else
    # Create the header
    echo -ne "SAMP\tRAW\tCLEAN\tALNRATE\tPROPER" > $samp/04_fixed/summary.txt
    echo -ne "\tPRRT_LEN\tPRRT_RC\tPRRT_COV" >> $samp/04_fixed/summary.txt
    echo -ne "\tINT_LEN\tINT_RC\tINT_COV" >> $samp/04_fixed/summary.txt
    echo -e "\tENV_LEN\tENV_RC\tENV_COV" >> $samp/04_fixed/summary.txt

    #--- Number of read pairs before QC
    if [[ -z ${numraw+x} ]]; then
        if [[ -e $samp/00_raw/original_1.fastq ]]; then
            numraw=$(( $(wc -l < $samp/00_raw/original_1.fastq) / 4 ))
        else
            numraw=$(( $(gunzip -c $samp/00_raw/original_1.fastq.gz | wc -l ) / 4 ))
        fi
    fi
    
    #--- Number of read pairs after QC
    numclean=$(head -n1 PL-94/04_fixed/bowtie2.out  | cut -d' ' -f1)
    
    #--- Alignment rate (calculated by bowtie2)
    alnrate=$(grep 'overall alignment rate' $samp/04_fixed/bowtie2.out | sed 's/% overall alignment rate//')
    
    #--- Number of aligned proper pairs
    # Keep 66 (PROPER_PAIR,READ1)
    # Remove 268 (UNMAP,MUNMAP,SECONDARY)
    numpro=$(samtools view -f 66 -F 268 $samp/04_fixed/final.bam | wc -l)
    
    #--- Print to summary.txt
    echo -ne "$samp\t$numraw\t$numclean\t$alnrate\t$numpro" >> $samp/04_fixed/summary.txt
    
    #--- Number of reads aligned to each chromosome
    samtools idxstats $samp/04_fixed/final.bam > $samp/04_fixed/final.idxstat.txt

    for x in "PRRT" "INT" "ENV"; do
        n=$(grep -P "^>.+$x" $samp/04_fixed/consensus.fa | sed 's/^>//')
        if [[ -n $n ]]; then
            # l1=$(samtools view -H $samp/04_fixed/final.bam | grep -P "^@SQ.+$n" | sed 's/.\+LN://')
            l1=$(grep "$n" $samp/04_fixed/final.idxstat.txt | cut -f2)
            rc1=$(grep "$n" $samp/04_fixed/final.idxstat.txt | cut -f3)            
            c1=$(samtools depth -r $n  $samp/04_fixed/final.bam | grep -vc "\s0$")
        else
            l1="0"
            rc1="0"
            c1="0"
        fi
        echo -ne "\t$l1\t$rc1\t$c1" >> $samp/04_fixed/summary.txt
    done
    echo '' >> $samp/04_fixed/summary.txt
fi

echo "[---$SN---] ($(date)) Consensus assembly summary:"
column -t $samp/04_fixed/summary.txt

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

