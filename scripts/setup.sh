#! /bin/bash

### Symlink the data #####################################################################

cd /lustre/groups/cbi/Projects/HIV_prospective_consensus/
mkdir -p Data && cd Data

for f in /import/cbi/Archive/HIV_prospective/DCC_pilot/complete_set_dec02/PL*.fastq.gz; do
    echo "$f"
    ln -s "$f"
done

for f in /import/cbi/Archive/HIV_prospective/DCC\ 4/DC\ Cohort\ 4_121216/Fastq\ files/PL*.fastq.gz; do
    echo "$f"
    ln -s "$f"
done

for f in /import/cbi/Archive/HIV_prospective/DCC\ 5/DCC_5\ Fastq/PL*.fastq.gz; do
    echo "$f"
    ln -s "$f"
done

for f in /import/cbi/Archive/HIV_prospective/DCC_6_7/DCC_6_7/PL*.fastq.gz; do
    echo "$f"
    ln -s "$f"
done

for f in /import/cbi/Archive/HIV_prospective/DCC_8/PL*.fastq.gz; do
    echo "$f"
    ln -s "$f"
done

ls *_R1_*.fastq.gz | cut -f1 -d'_' > samples.txt

cd ..

### Setup reference directories ##########################################################
mkdir -p ref && cd ref
cp /c1/apps/trimmomatic/Trimmomatic-0.33/adapters/NexteraPE-PE.fa adapters.fa
echo '' >> adapters.fa

# Also, download the FASTA for HXB2 K03455 from NCBI
# Save as HIV_B.K03455.HXB2.fasta 

# Create references for the sequenced regions

##########################################################################################
# This is python:
python << EOF
from Bio import SeqIO
s = SeqIO.read('HIV_B.K03455.HXB2.fasta', 'fasta')
with open('HXB2.regions.fasta', 'w') as outh:
    print >>outh, '>HIV_B.PRRT'
    ss = str(s[2084:3869].seq)
    print >>outh, '\n'.join(ss[i:i+60] for i in range(0, len(ss), 60))
    print >>outh, '>HIV_B.INT'
    ss = str(s[4229:5096].seq)
    print >>outh, '\n'.join(ss[i:i+60] for i in range(0, len(ss), 60))
    print >>outh, '>HIV_B.ENV'
    ss = str(s[6224:8795].seq)
    print >>outh, '\n'.join(ss[i:i+60] for i in range(0, len(ss), 60))

EOF
##########################################################################################

cd ..


### Create analysis directories ##########################################################
mkdir -p Analysis && cd Analysis

cp ../Data/samples.txt .

cat samples.txt | while read s; do
    echo "Copying data for $s"
    mkdir -p $s/00_raw
    cp ../Data/${s}_*R1*.fastq.gz $s/00_raw/original_1.fastq.gz
    cp ../Data/${s}_*R2*.fastq.gz $s/00_raw/original_2.fastq.gz    
done

# Symlink scripts here to make paths shorter
ln -s ../scripts/haphpipe_consensus.sh
ln -s ../scripts/run_array.sh

cd ..
