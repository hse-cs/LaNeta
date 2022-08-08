# LaNeta
Statistical method for estimating parameters (e. g. the timing) of two pulse model based on three locus Linkage Disequilibrium and Local Ancestry.


# Installation

For performance reasons, we use cython to speed up calculations, so you need
to compile `.pyx` by yourself. For this, you need a working toolchain for building C
code (gcc and clang are known to work).

First, install the dependencies

```
$ python3 -m pip install numpy>=1.19.5 cython
$ python3 -m pip install msprime
```

To compile .pyx you should use

```
$ python3 setup.py build_ext -i
```

## Quick start
We use .vcf files for admixed and two source populations and simple txt .map and .pop files with morgan units:

To specify all necessary use these parameters:
`-vcf vcf/dir.vcf.gz`
`-m morgans.map`
`-p populations.pop`

Also you need to specify which populations from `.pop` are admixed and source:

`-p0 ADM` - admixed population

`-p1 SRC1` - admixed 1 time

`-p2 SRC2` - admixed 2 times
---
morgans .map format:
```
  CHR_NAME VAR_ID POS_PB POS_MORGANS
```
Make sure that there are no duplicates in this file!

populations .txt format:
```
  SAMPLE_ID POPULATION
```

---

Example for real data:
```
python3 laneta.py -b 0.01 -vcf mer.vcf.gz -p populations.pop -m morgans.map -p0 CLM -p1 YRI -mt 0.94 -jk -nmt
```
If you specify only one source population, admixed population is separated into two equal-sized groups. These groups are used as the admixed and the missing source population.


## All settings
### flags
`-ht` if you use haplotype data instead of genotype.

`-jk` for calculation of confidence intervals using jackknife by leaving out each chromosome.

`-nmt` if you want `mt` to be used only for estimating frequencies from missing source population and not for estimating parameters.

### bracket parameters
`-b` sets the distance between centres of brackets for FFT (default is 0.01).

`-r` sets the radius of brackets for FFT (default is a half of `b`).
### files
`-vcf` specifies .vcf file that contains data for all populations.

`-p` .pop file that contains samples with indicated population.

`-m` .map file with chromosome name(1-22), var id, var position(bp), var position(cm).

`-p0` name of the admixed population in population .txt file.

`-p1` name of the first source population in population .txt file.

`-p2` name of the second source in popultion .txt file.

### two pulse parameters
`-m1`, `-m2`, `-mt` used for setting adm. proportions for times estimation. If you know the total ancestry proportion mt you can specify it for better estimation of other parameters.
`-t1`, `-t2` used for setting adm. times.
### cm parameters
`-min` and `-max` specifies min and max genetic distance for estimations (in cantimorgans).

## Output
Output is a tab delimited list of parameters:
T1 - time between two admixture events, T2 - time to the most recent admixture event,
M1 - admixture proportion of the first admixture event,
M2 - admixture proportion of the most recent admixture event.


If jackknife used, output have additional lines:

jk bias

jk var

jk 95% conf. interval

## Data preparation for analysis

Files which are required for data preparation are places in utilites folder. You need your vcf.gz file and its index file, also check if `bcftools`, `vcftools` and `plink` are installed.

1. Execute *get_interpolation_files.sh* to get recombination maps for future interpolation. It will be in folder map.
2. To execute *preparation.sh* you will have to enter the name of your file without .vcf. As a result, you will have .txt file with 4 columns: `CHROM` `ID` `POS` `GEN_POS`. The forth column refers to morgan units.
3. Repeat with all populations you have
4. Merge all your population file by running
  ```
    bcftools merge pop1.vcf pop2.vcf -Oz -o merged.vcf.gz
    plink --vcf merged.vcf.gz --chr 1-22 --snps-only --recode vcf --out filtered_merged
    bcftools view -m2 -M2 -c20 -v snps -Oz -o filtered_merged.vcf.gz filtered_merged
  ```
  filtered_merged.vcf.gz is the file for your further analysis.
