#! /bin/bash

if [[ -n ${SLURM_ARRAY_TASK_ID+x} && -n ${idlist+x} ]]; then
    samp=$(sed -n "$SLURM_ARRAY_TASK_ID"p  $idlist)
else
    exit 1
fi

cmd=". haphpipe_consensus.sh $samp ../ref/HXB2.regions.fasta ../ref/adapters.fa &> $samp/hp.log"
echo -e "Command line:\n$cmd"
eval $cmd

exit 0
