# Haphpipe Consensus

First follow `scripts/setup.sh` to setup the analysis directory.

Launch the array jobs like this within the `Analysis` directory

```bash
sbatch -N 1 -a 1-$(wc -l < samples.txt) -p short,defq -t 240 run_array.sh
```
