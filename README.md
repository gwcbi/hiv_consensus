# HIV Consensus pipeline

Created this pipeline for reference-guided assembly of HIV prospective data.
The reads are first trimmed using Trimmomatic, then aligned to HXB2.
Genotype calls (VCF) and consensus (FASTA) is generated from these alignments.
A second (optional) refinement step re-aligns the reads to
the consensus generated in the first pass and calls a refined consensus
sequence. This (may) result in increased proportion of mapped reads.

The following scripts are used to generate the data:

+ **`setup.sh`** Describes how the project directory was set up.
+ **`haphpipe_consensus.sh`** HAPHPIPE pipeline for generating consensus
+ **`run_array.sh`** Script to run all samples in a SLURM array job

After setting up the directory according to `setup.sh`, all the analysis can be run
using the following array job (from within the `Analysis` directory):

```bash
sbatch \
 -N 1 -p short,defq -t 240 \
 -a 1-$(wc -l < samples.txt) \
 --export idlist=samples.txt \
 run_array.sh

```
