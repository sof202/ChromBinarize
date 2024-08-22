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



