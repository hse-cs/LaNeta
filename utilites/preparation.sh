#!/bin/sh


read -p "Inser name of ONE population to be transformed to morgans " answer

for i in (echo ${answer} | tr " " "\n")
  do
    plink --vcf ../${i}.vcf.gz --chr 1-22 --snps-only --make-bed --out ${i}
    plink --bfile ${i} --recode vcf --out ${i}.${i}
    bcftools query -f '%CHROM\t%POS\t%ID\n' ${i}.${i}.vcf >> gen.${i}.txt
    for j in {1..22}; do grep -w ${j} gen.${i}.txt >> gen.${i}.chr${j}.txt; done
    srun python phys2gen.py ${i}
    srun cat ${i}.chr{1..22}.txt >> ${i}.gen.pos.txt
    rm ${i}.{bim,bed,fam}
    rm ${i}.chr{1..22}.txt
    rm gen.${i}.chr{1..22}.txt
    rm gen.${i}.txt
  done
 
  
