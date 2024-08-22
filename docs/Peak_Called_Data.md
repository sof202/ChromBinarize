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

