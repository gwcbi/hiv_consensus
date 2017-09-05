#! /bin/bash

if [[ -n ${SLURM_ARRAY_TASK_ID+x} && -n ${idlist+x} ]]; then
    samp=$(sed -n "$SLURM_ARRAY_TASK_ID"p  $idlist)
else
    exit 1
fi

. haphpipe_consensus.sh $samp ../ref/HIV_B.K03455.HXB2.fasta ../ref/adapters.fa &> $samp/hp.log

