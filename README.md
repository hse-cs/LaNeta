# LaNeta
Statistical method for estimating the timing of multiple admixture events based on three locus Linkage Disequilibrium.


## Installation

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
$ python setup.py build_ext -i
```

## Settings
### Simulations


You can use the method with simulated data for two pulse model. We use msprime for simulations.
To simulate genetic data you should use `-ms` flag and `simulator_setup.txt`. Use `-seed` parameter to set random seed.
Example for simulations:
```
python thld.py -ms -e 0.01 -mt 0.16
```

### Real data
We use .vcf files for admixed and two source populations and .txt file with morgan units:

To specify .vcf  and .txt use these parameters:
`-a admixed/population/dir.vcf.gz`
`-s1 first/source/pop/dir.vcf.gz`
`-s2 second/source/pop/dir.vcf.gz`
`-m morgans.txt`

.txt format:
```
  CHR_NUMBER VAR_ID POS_PB POS_MORGANS
```
Example for real data:
```
python thld.py -a admixed/population/dir.vcf.gz -s1 first/source/pop/dir.vcf.gz -s2 second/source/pop/dir.vcf.gz -m yri_clm/mapfile.txt -e 0.01 -cm 30 -m2 0.9925 -m1 0.4175 -af
```
If you specify only one source population, admixed population is separated into two equal-sized groups. These groups are used as the admixed and the missing source population.


### Another settings

`-e` sets the length of brackets for FFT (default is 0.01).

`-af` flag is needed if you want to estimate and subtract affine term that due to population substructure.

`-cm` parameter specifies max genetic distance for estimations (in cantimorgans).

`-m1`, `-m2`, `-mt` parameters can be used for setting adm. proportions for times estimation. You can specify m1 and m2 or the total ancestry proportion mt.    
