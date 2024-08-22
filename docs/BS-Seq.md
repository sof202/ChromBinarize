## BS-Seq

For the most part, the BS-Seq scripts follow the same process as the [ONT
pipeline](./ONT.md). The main difference is the expected format of the input
file. The assumed format is that of the output from
[wgbs_tools](https://github.com/nloyfer/wgbs_tools). Specifically the
`beta2bed` command.

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

