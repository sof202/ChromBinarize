# ONT methylation processing

This is a selection of scripts that will process bed files that are produced
from ONT's pipeline. Files have the form:

|Chromosome|Start|End|methylation-type|coverage|strand|percent-methylation
|----------|-----|---|----------------|--------|------|-------------------|

## Running

In order to run these scripts you will need to first fill out the config file
provided. Ideally you would then put this config file next to your data.
Then, you would call scripts sequentially using SLURM workload manager with:
```bash
sbatch path/to/script path/to/config/file
```

## Script purposes

The binomial script will try to find an estimate for the probability that an
erroneous read can occur on each of the extremes (probability of incorrect 
unmethylated read for methylated sites and probability of incorrect methylated 
read for unmethylated sites). It then uses these values with a binomial 
distribution to determine the sites that are:

1) Methylated
2) Unmethylated
3) Neither (not necessarily hemi-methylated)

From here the other script will attempt to binarise the methylation data so
that it can be used suitably with ChromHMM. It employs the same approach that
ChromHMM does (using a poisson distribution to tell apart peaks from the
background). However, unlike ChromHMM, it is only considering bins/windows
where reads do in fact exist to calculate the 'background' signal. It's not
really a background signal, instead this process is more trying to discern
methylated CpG rich areas against 'randomly' situated methylated sites.

## Supplementary scripts

These scripts net you some extra statistics if you want them.

- a1_erroneous_rate_plot.sh will show you how the value of p (binomial 
distribution parameter) changes with different reference datasets.
- b1_CpG_robustness.sh will show you how correlated the methylation status
of nearby CpGs are. Ideally, close CpGs are highly correlated.

