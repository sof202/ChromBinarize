<p align="center">
  <img src="https://github.com/sof202/ChromBinarize/assets/147140110/811a4728-d701-4d99-a8a6-d3cbd6f145b1" />
</p>

</p>
<p align="center">
    <a href="https://img.shields.io/codefactor/grade/github/sof202/ChromBinarize" alt="CodeFactor">
        <img src="https://img.shields.io/codefactor/grade/github/sof202/ChromBinarize" /></a>
    <a href="https://img.shields.io/github/commit-activity/m/sof202/ChromBinarize" alt="Commit activity">
        <img src="https://img.shields.io/github/commit-activity/m/sof202/ChromBinarize" /></a>
</p>

# ChromHMM Binarization Tools

This is a selection of scripts that will convert various bed files for ONT, 
oxBS and WGBS datasets into a format compliant with ChromHMM. 

ChromHMM is great at binarizing at a simple level, but struggles for datasets
that are not traditionally peak called. In addition to this, 'better' peak
calling algorithms (like MACS) exist for ChIP-Seq and ATAC-Seq datasets. As
such, a separate suite of scripts that binarize these data sets is proposed
here.

## Running

In order to run these scripts you will need to first fill out the config file
(template provided). Ideally you would then put this config file next to your 
data (though, realistically you can put this anywhere you wish).
Then, you would call scripts sequentially using SLURM workload manager with:
```bash
sbatch path/to/script path/to/config/file
```

## Included scripts

The scripts in this repository are split into the following categories:

- [ONT](#ont)
- [BS-Seq](#bs-seq)
- [ChIP-Seq](#chip-seq)
- [supplementary scripts](#supplementary)

## ONT

This pipeline expects input bed files of the following format (standard output
of ONT's [modkit](https://github.com/nanoporetech/modkit)):

|Chromosome|Start|End|methylation-type|coverage|strand|percent-methylation|
|----------|-----|---|----------------|--------|------|-------------------|

These scripts have been created in an attempt to binarize the methylation and
hydroxymethylation calls that come out of ONT data. The process is split into
two steps:

### Purification 

In this step, 'poor sites' are removed from the dataset. Here 'poor sites' is 
defined as: 

> Sites that have a low read depth or are unlikely to be true 
(hydroxy)methylation signal.

This is done as there exists two sources of error in the 'percent methylation'
metric that is given by ONT's pipeline (modified basecaller like `dorado` 
alongside bed file generator like `modkit`). The first is errors in the 
original signal track generated from the sequencing experiment. The second is
errors from the base caller itself. We don't know what the true methylation
signal looks like, so we can only remove sites that are 'unlikely' to be true
methylation signal. 

Sites are determined to be true methylation signal if they have a high enough
read depth (user determined) and they pass a statistical test using a binomial
distribution. The test attempts to answer this question:


> How probable is it that unmethylated cells are erroroneous or 
unrepresentative cell types (of the total population)?

This question is answered using a binomial distribution using probabilities 
gathered from a 'good' reference set. By default, this 'good' reference set is
a subset of the original data that only has sites with >95% methylation and
very high read depth (>500). It is assumed that all reads in this dataset that
are unmethylated are erroneous (or unrepresentative of the total population).
Thus this subsetted dataset can be used to estimate a probability to use with
the binomial distribution that can be applied to the whole dataset.

### Binarization

This step converts the now 'purified' ONT reads into the binarized data format
that ChromHMM expects. That is, a single vector of 0s and 1s for each genomic
bin/window/region.

It is hard to binarize methylation data. Identifying regions that are ~200bp in
size that have a single CpG with methylation signal as 'regions characterised
by methylation' feels wrong. However, restricting genomic bins that are 
'characterised by methylation' to those with lots of methylated sites is very
strict (resulting in almost a 0 vector). To combat this problem, this pipeline
gets the best of both worlds. The binarization pulls out 'dense' and 'sparse'
regions of (hydroxy)methylation.

Dense regions of (hydroxy)methylation are defined as:

> Bins that have abnormally high numbers of methylated CpGs in comparison to
the background

'Abnormally high' is determined using the Poisson distribution with a threshold
of 0.0001.

Sparse regions of (hydroxy)methylation are defined as the remaining bins that
have at least one site that is methylated.

It is up to the user to decide whether to include each of these binary files
in their subsequent analysis (you could even combine them if you wish).

## BS-Seq

For the most part, the BS-Seq scripts follow the same process as the ONT
pipeline. The main difference is the expected format of the input file:

|Chromosome|Start|End|number of methylated reads|total reads|
|----------|-----|---|--------------------------|-----------|

If you have oxidative bisulphite sequencing data (alongside regular bisulphite 
sequencing data), you can use them in tandem to extract sites that are exactly
5mC. Oxidative bisulphite sequencing allows you to discern which sites are
hydroxymethylated, and as such this can be used to disentangle your bisulphite
sequencing data into 5mC and 5hmC.

## ChIP-Seq

It is currently assumed that your ChIP-Seq data has already been peak called
by an external program such as [MACS](https://github.com/macs3-project/MACS).

This is a very simple script that converts the peaks called by such a program
into the binary format that is expected by ChromHMM. Make sure your input bed
file for this script is of the format specified 
[here](https://macs3-project.github.io/MACS/docs/callpeak.html#output-files).

All that is really required is the chromosome, start and end fields. Just
ensure that these fields are indeed corresponding with the called peaks and
not something else.

## Supplementary

There are 4 different scripts included in this folder. They are mainly here to
empower the user to further inspect their data.

### Erroneous rate plot

This script is in place so that the user can inspect how the arbitrary decision
of choosing a 'good' reference set changes the parameter used with the binomial
distribution in the purificaiton step (for ONT and BS-Seq).

In an ideal world, this arbitrary decision should not matter, *i.e.* peturbing
the values will not change the output parameter too heavily. This isn't always
the case, and so it is wise to check if this bothers you.

### CpG robustness

This script was born out of distrust in a certain base caller. If a basecaller's
output is to be trusted, you at least expect it to be consistent with itself.
This script plots how similar the (hydroxy)methylation calls are between
nearby CpGs. Ideally, close CpGs should generally be of the same methylation
status. It would indeed be strange if CpGs with no methylation are usually
situated next to fully methylated CpGs. This script allows you to check this
hypothesis for your data.

### WGBS comparison

This script allows you to compare two methylation datasets. One being from ONT
data, and the other being from WGBS. Ideally, if the two datasets are of the
same cell type and from the same sample, they should match up pretty well. 

This script outputs 3 histograms. 

1) The absolute change in read depth between the same sites (sense checking)
2) The absolute change in methylation percent at the same site
3) The methylation percent in ONT vs the methylation percent in WGBS 
(2d histogram)

### oxBS comparison

This is the same as WGBS comparison, but now expects oxBS data instead. This
allows the user to specifically look at how the hydroxymethylation calls 
compare between ONT and oxidative bisulphite sequencing. 

