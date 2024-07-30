<p align="center">
  <img src="https://github.com/sof202/ChromBinarize/assets/147140110/811a4728-d701-4d99-a8a6-d3cbd6f145b1" />
</p>

</p>
<p align="center">
    <a href="https://www.codefactor.io/repository/github/sof202/chrombinarize" alt="CodeFactor">
        <img src="https://img.shields.io/codefactor/grade/github/sof202/ChromBinarize?style=for-the-badge&color=red" /></a>
    <a href="https://github.com/sof202/ChromBinarize/commits/main/" alt="Commit activity">
        <img src="https://img.shields.io/github/commit-activity/m/sof202/ChromBinarize?style=for-the-badge&color=red" /></a>
    <a href="https://github.com/sof202/ChromBinarize/blob/main/LICENSE" alt="License">
        <img src="https://img.shields.io/github/license/sof202/ChromBinarize?style=for-the-badge&color=red" /></a>
</p>

# ChromHMM Binarization Tools

This is a selection of scripts that will convert various bed files for ONT, 
oxBS and WGBS datasets into a format compliant with ChromHMM. 

ChromHMM is great at binarizing at a simple level, but struggles for datasets
that are not traditionally peak called. In addition to this, 'better' peak
calling algorithms (like MACS) exist for ChIP-Seq and ATAC-Seq datasets. As
such, a separate suite of scripts that binarize these datasets (into a format
recognised by ChromHMM) is proposed here.

>[!NOTE]
>In the following README (and greater repository), the word 'methylation' means
>*any* type of DNA methylation. As such, when more precise language is required,
>you will see instead '5mC' or '5hmC' (*etc.*). If at any point the wording 
>feels ambiguous when it shouldn't be, please raise an issue.

## Setup
In order to run these scripts you will need to first fill out the config file
(template provided in `./config-setup.txt`). It is recommended that you put 
this config file near your data (note: this is not a requirement, you can 
actually put this file anywhere you wish).

Next run the setup script with:

```bash
./setup
```

This setup script requires user input for removing SLURM directives and also
when setting up the conda environment. This was a conscious decision as you may
want to check what is being installed by conda first. Also, this setup script
will take quite some time due to the dependency tree (~49 packages) for R. 

You will see the following message on success:

```bash
[1] "success"
```

## Usage
After completing setup, run scripts sequentially using SLURM workload manager:

```bash
sbatch path/to/script path/to/config/file
```

> [!NOTE]
> If you want to get a quick summary of what a script does, run the script
> without any positional parameters (you can just run it like a normal bash
> script in this case, `sbatch` is not required).

## Software Requirements 

This pipeline requires a unix-flavoured OS and requires the following software
to be installed. Versions are those that were used during testing, lower minor
version numbers are likely to still work.

- [bash](https://www.gnu.org/software/bash/) (>=4.2.46(2))
- [SLURM Workload Manager](https://slurm.schedmd.com/overview.html) (>=20.02.3)
- [GNU awk](https://www.gnu.org/software/gawk/) (>=4.0.2)
- [GNU gzip](https://www.gnu.org/software/gzip/) (>=1.5)

The following software and R packages are installed for you in the `setup`
script:

- [Bedtools](https://github.com/arq5x/bedtools2) (>=v2.29.2)
- [R](https://www.r-project.org) (>=4.2.1)
    - [dplyr](https://dplyr.tidyverse.org)
    - [data.table](https://github.com/Rdatatable/data.table)
    - [fitdistrplus](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
    - Supplementary scripts only:
        - [grid](https://github.com/cran/grid)
        - [gridExtra](https://github.com/baptiste/gridExtra)
        - [cowplot](https://github.com/wilkelab/cowplot)

## Included scripts

The scripts in this repository are split into the following categories:

- [ONT](#ont)
- [BS-Seq](#bs-seq)
- [ChIP-Seq and ATAC-Seq](#chip-seq-and-atac-seq)
- [Supplementary scripts](#supplementary)

## ONT

This pipeline expects input bed files of the following format (standard output
of ONT's [modkit](https://github.com/nanoporetech/modkit)):

|Chromosome|Start|End|Methylation type|Coverage|Strand|Percent methylation|
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

> [!NOTE]
> Optionally, you can pass the `-c` flag into the purification step to further
> filter the CpGs in the reference set to only be those that exist within
> known CpG islands (CGIs). CpG islands are usually unmethylated, making them
> a more reliable set (for CpGs that are truly unmethylated).
>
> However, if you are convinced that modified basecallers (such as dorado)
> perform worse when it comes to calling the methylation status of isolated
> CpGs: Do not use this option. Isolated CpGs have less context around them
> and so modified basecallers may have less accuracy.
>
> Two cpg islands files are provided out of the box in the 'references'
> directory of this repository (hg19 and hg38). These were obtained from
> UCSC's [data integrator](https://genome.ucsc.edu/cgi-bin/hgIntegrator).
> If your data is not in these assemblies, you might be able to find the
> reference CGI file for your dataset here aswell.

Once a suitable reference set is created, sites are determined to be true 
methylation signal if:
- They have a high enough read depth (user determined) 
- They pass a statistical test using a binomial distribution. The test 
attempts to answer this question:

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

It is hard to binarize methylation data. For example, identifying regions that are 200bp in
size that have a single CpG with methylation signal as 'regions characterised
by methylation' feels wrong. However, restricting genomic bins that are 
'characterised by methylation' to those with lots of methylated sites is very
strict (resulting in almost a 0 vector). To combat this problem, this pipeline
gets the best of both worlds. The binarization pulls out 'dense' and 'sparse'
regions of (hydroxy)methylation.

Dense regions of (hydroxy)methylation are defined as:

> Bins that have significantly high numbers of methylated CpGs 
in comparison to the background

'Significantly high' is determined using the beta distribution
(with a default threshold of 0.001).

Sparse regions of (hydroxy)methylation are defined as the remaining bins that
have at least one site that is methylated.

It is up to the user to decide whether to include each of these binary files
in their subsequent analysis (you could even combine them if you wish).

## BS-Seq

For the most part, the BS-Seq scripts follow the same process as the ONT
pipeline. The main difference is the expected format of the input file. The
assumed format is that of the output from 
[wgbs_tools](https://github.com/nloyfer/wgbs_tools). Specifically the `beta2bed`
command.

|Chromosome|Start|End|Number of methylated reads|Total reads|
|----------|-----|---|--------------------------|-----------|


### Extracting 5hmC signal

If you have oxidative bisulphite sequencing data (alongside regular bisulphite 
sequencing data), you can use them in tandem to extract sites that are exactly
5hmC. Oxidative bisulphite sequencing allows you to discern which sites are
exactly 5mC. As such, this can be used to disentangle your bisulphite
sequencing data into 5mC and 5hmC.

In order to extract sites that have significant 5hmC signal, we unfortunately
can't just remove the oxBS signal from the WGBS signal and use the same
'purification' method from the normal pipeline. This is because it is (at least
currently) very unlikely that a site will be registered to have close to 100%
signal in WGBS but close to 0% signal in oxBS. It's definitely plausible that
this scenario *can happen* at certain sites in the genome, but (at least in our
data) this is not seen anywhere even in neuronal cells. 

We also can't run the pipeline for both WGBS and oxBS then just remove the
registered 5mC signal from the full methylation signal given by the WGBS data 
(*i.e* perform a vector subtraction between the binary files). This is because
a bin could feasibly be described by both 5mC and 5hmC, the marks are not
mutually exclusive within regions of the genome.

As a result, we take a different approach for extracting 5hmC. What we really 
care about is the answer to the following question:

> Is there significant evidence to suggest that a site has a stronger signal
in the WGBS data versus the oxBS data?

If the answer to this question is yes, then we can infer that the contribution
to the signal that 5hmC makes is significant. As such we can say that the site
in question is significantly hydroxymethylated.

## Peak called data (ChIP, ATAC *etc.*)

It is currently assumed that your peak called data has already been peak 
called by an external program such as 
[MACS](https://github.com/macs3-project/MACS).

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

### Changing bin size

This is a script that allows the user to change the bin size used for a
binary file. This can be very useful when working with multiple modalities.
For example, you may want to use a bin size of 150bp for ATAC data, but use a
bin size of 300 for your BS-Seq data. This script allows you to do both of
these things, but still keep the binary files compatible within ChromHMM. You
can now convert the 300bp binary file into a 150bp binary file so that both
binary files have the same number of lines.

This works in the simplest way possible. Each bin in the new binary file
(that lies within a bin of the original binary file) will take on the same 
value as in the original binary file.

> [!WARNING]
> It is not reccomended that you convert a binary file to have a larger bin
> size using this script. This is because you can't make the assumption that
> you would still have a peak in a wider area. Imagine you combined 3 bins
> together where only one bin had a peak, you probably would agree that the
> combined bin shouldn't have a peak. But what about the case where you
> combine two bins, one with a peak and the other without. What should you do
> then? It is better to just avoid this entirely and re-call the binarization
> script for a larger bin size.

To learn how the script works (in terms of positional arguments), run the
script without any arguments.
