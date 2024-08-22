## Supplementary Scripts

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

