# LaNeta
Statistical method for estimating the timing of multiple admixture events based on three locus Linkage Disequilibrium and Local Ancestry.


# Installation

For performance reasons, we use cython to speed up calculations, so you need
to compile `.pyx` by yourself. For this, you need a working toolchain for building C
code (gcc and clang are known to work).

First, install the dependencies

```
$ python -m pip install numpy>=1.19.5 cython
$ python -m pip install msprime
```

To compile .pyx you should use

```
$ python3 setup.py build_ext -i
```

# Settings
## Simulations


You can use the method with simulated data for two pulse model. We use msprime for simulations.
To simulate genetic data you should use `-ms` flag and `simulator_setup.txt`. Use `-seed` parameter to set random seed.
Example for simulations:
```
python3 laneta.py -ms -e 0.01 -mt 0.16
```

## Real data
We use .vcf files for admixed and two source populations and .txt file with morgan units:

To specify .vcf  and .txt use these parameters:
`-vcf vcf/dir.vcf.gz`
`-m morgans.txt`
`-pf pop.txt`

Also you need to specify which populations from `pop.txt` are admixed and source:
`-p0 ADM`
`-p1 SRC1`
`-p2 SRC2`

morgans .txt format:
```
  CHR_NUMBER VAR_ID POS_PB POS_MORGANS
```

populations .txt format:
```
  SAMPLE POPULATION
```

Example for real data:
```
python3 laneta.py -e 0.01 -vcf data_yri_clm/mer.vcf.gz -pf data_yri_clm/pop.txt -m data_yri_clm/map.txt -p0 CLM -p1 YRI -mt 0.95
```
If you specify only one source population, admixed population is separated into two equal-sized groups. These groups are used as the admixed and the missing source population.


## All settings
### flags
`-gt` if you use genotype data instead of haplotype.

`-af` is needed if you want to estimate and subtract affine term that due to population substructure.

`-jk` for calculation of confidence intervals using jackknife by leaving out each chromosome.

`-ms` for running on simulations.
### bracket parameters
`-e` sets the distance between centres of brackets for FFT (default is 0.01).

`-r` sets the radius of brackets for FFT (default is half of `e`).
### files
`vcf` specifies .vcf file that contains data for all populations.

`pf` .txt file that contains samples with indicated population.

`-m` .txt file with chromosome name(1-22), var id, var position(bp), var position(cm).

`-p0` name of the admixed population in population .txt file.

`-p1` name of the first source population in population .txt file.

`-p2` name of the second source in popultion .txt file.

### proportions
`-m1`, `-m2`, `-mt` used for setting adm. proportions for times estimation. You can specify m1 and m2 or the total ancestry proportion mt.
### cm parameters
`-min` and `-max` specifies min and max genetic distance for estimations (in cantimorgans).
### other
`-seed` specifies random seed.



## Data preparation for analysis

Files which are required for data preparation are places in utilites folder. You need your vcf.gz file and its index file, also check if `bcftools`, `vcftools` and `plink` are installed.

1. Execute *get_interpolation_files.sh* to get recombination maps for future interpolation. It will be in folder map.
2. To execute *preparation.sh* you will have to enter the name of your file without .vcf. At the end you will have .txt file with 4 columns: `CHROM` `ID` `POS` `GEN_POS`. The forth column refers to morgan units.
3. Repeat with all populations you have
4.
```
  comm -12 pop1.gen.pos.txt pop2.gen.pos.txt | awk '{print("%s\n", $2)}' > common_snps
  vcftools --vcf ../pop1.pop1.vcf --snps common_snps --recode --recode-INFO-all | bcftools view -Oz -o filtered.pop1.vcf
  vcftools --vcf ../pop2.pop2.vcf --snps common_snps --recode --recode-INFO-all | bcftools view -Oz -o filtered.pop2.vcf
 ```
