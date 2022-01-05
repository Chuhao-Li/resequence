#!/bin/bash
if [ ! -d report/data ]; then mkdir report/data; fi

for i in $(cat samples.list); do 
    cp -l BWA/${i}/${i}.final.vcf report/data/${i}.snv.vcf
    cp -l cnv/${i}/${i}.del.xls report/data/${i}.low.xls
    cp -l cnv/${i}/${i}.high.xls report/data/${i}.high.xls
    cp -l delly/${i}/${i}.filtered.vcf report/data/${i}.sv.vcf
done

cp -l resequence.html report
cp -lr snapshot report

cd report
cp SV/clustered.SV.txt data
cp SNV/clustered.SNV.txt data
rm -r SV SNV
tar -zcvf data.tar.gz data
rm -r data
cd ../

# generate_html.py readme.md
# project_manage.py publish this
