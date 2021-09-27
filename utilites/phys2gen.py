import pandas as pd
import matplotlib.pyplot as plt
import scipy
import os, sys


# gen_map is hapmap for hg19

population_name = ' '.join(sys.argv[1:])

for i in range(1, 23):
        filename = "map/chr" + str(i) + ".interpolated_genetic_map"
        gen_map = pd.read_csv(filename, delim_whitespace=True, names=['ID','POS','GEN_POS'])
        filename = "gen." + population_name + ".chr" + str(i) + ".txt"
        snp = pd.read_csv(filename, sep='\t', names = ['CHROM', 'POS', 'ID'])
        map_can = snp.merge(gen_map, how='left', on=['ID', 'POS'])
        if map_can['GEN_POS'].isnull().values.any() == True and map_can['GEN_POS'].isnull().values.all() == False:
                c = map_can.interpolate().dropna()
                c = c[['CHROM', 'ID', 'POS', 'GEN_POS' ]]
                filename = population_name + ".chr" + str(i) +".txt"
                c.to_csv(filename, header=None, index=None, sep=' ')
        else:
                print("There is nothing to evaluate, check your files")
