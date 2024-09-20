<p align="center">
  <img src="https://github.com/sof202/ChromBinarize/assets/147140110/811a4728-d701-4d99-a8a6-d3cbd6f145b1" />
</p>

</p>
<p align="center">
    <a href="https://github.com/sof202/ChromBinarize/actions/workflows/test.yml">
      <img src="https://img.shields.io/github/actions/workflow/status/sof202/ChromBinarize/test.yml?style=for-the-badge&color=red" />
    </a>
    <a href="https://github.com/sof202/ChromBinarize/commits/main/" alt="Commit activity">
        <img src="https://img.shields.io/github/commit-activity/m/sof202/ChromBinarize?style=for-the-badge&color=red" />
    </a>
    <a href="https://github.com/sof202/ChromBinarize/blob/main/LICENSE" alt="License">
        <img src="https://img.shields.io/github/license/sof202/ChromBinarize?style=for-the-badge&color=red" />
    </a>

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
when setting up conda environments. This was a conscious decision as you may
want to check what is being installed by conda first. Also, this setup script
will take quite some time due to the dependency tree (~49 packages) for R. 

You will see the following message on success:

```bash
[1] "success"
```

If you do not see this success message, please open up an
[issue](https://github.com/sof202/ChromBinarize/issues/new?assignees=&labels=bug&projects=&template=bug-report.yaml&title=%5BBug%5D%3A+).

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
- [Conda](https://conda.io/projects/conda/en/latest/index.html)
    - Any installation will do, this has worked on Miniconda 4.5.2 (from 2020)
    - Make sure conda can be found on your `PATH` (check with `which conda`)
- [GNU awk](https://www.gnu.org/software/gawk/) (>=4.0.2)
- [GNU gzip](https://www.gnu.org/software/gzip/) (>=1.5)

The following software and R packages are installed for you in the `setup`
script:

- [Bedtools](https://github.com/arq5x/bedtools2) (v2.29.2)
- [R](https://www.r-project.org) (4.4.1)
    - [dplyr](https://dplyr.tidyverse.org)
    - [data.table](https://github.com/Rdatatable/data.table)
    - [fitdistrplus](https://cran.r-project.org/web/packages/fitdistrplus/index.html)
    - Supplementary scripts only:
        - [grid](https://github.com/cran/grid)
        - [gridExtra](https://github.com/baptiste/gridExtra)
        - [cowplot](https://github.com/wilkelab/cowplot)

## Further documentation

Please consult the [wiki](https://sof202.github.io/ChromBinarize) for further
documentation on specific scripts.
