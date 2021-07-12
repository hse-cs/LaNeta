# ThLd
Statistical method for estimating the timing of multiple admixture events based on three locus Linkage Disequilibrium.


Installation
--------------------

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

### Settings
Now it is possible to use the method only with simulated data for two pulse model. We use msprime for simulations. 
To simulate genetic data you should use `-ms` flag and `tap.txt`
file:
```#tap.txt
10  10
0.2 0.2
```
In the first row you should put times `T1` and `T2`,
and in the second row you should put proportions `m1` and `m2`. 

Use `-cn` or `--chrn` to set the number of chromosomes (default is 20).

Use `-p` to set the number of bins for FFT (default is 300).

Example:
```
python setup.py build_ext --inplace
python thld.py -ms -p 1000 -cn 23

```
