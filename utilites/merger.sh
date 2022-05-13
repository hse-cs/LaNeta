#!/bin/sh

sample=$1

find $sample -name '*.vcf' \
| xargs -tI{} bgzip -f {} > {}.gz

find $sample -name '*.gz' \
| xargs -tI{} tabix -C -p vcf {}

vsf_files=`find $sample -name '*.gz' | sort -k 5 -t _ -n`
echo "concating: \n $vsf_files \n to $sample/sample_$sample.vcf.gz"
vcf-concat $vsf_files | bgzip -c > $sample/sample_$sample.vcf.gz

tabix -C -p vcf $sample/sample_$sample.vcf.gz

exit 0
