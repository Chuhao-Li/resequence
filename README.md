# Introduction
Bacterial variants(structural variant and single nucleotide variant) calling pipeline by short read NGS data. 

# INSTALL

The following tools should be in the environment: 
- samtools v1.10
- bcftools v1.10.2
- bwa v0.7.17-r1188
- java openjdk v11.0.13
- snakemake v5.4.5
- fastqc v0.11.9
- multiqc v1.7
- python v3.7.3

The following tools has pre-installed in this package: 
- vcfutils.pl
- picard
- annovar v2018-04-16
- delly v0.8.5

The following python package should have beed install: 
- pandas

to install this package: 
``` bash
git clone https://github.com/Chuhao-Li/resequence.git
```

# How to run

1. Edit conf.json

2. Run by the pipeline using example files as input: 
``` bash
cd /path/to/package_directory

snakemake -s resequence.snk.py --config BASE_DIR=./ --configfile conf.json
```

The output should be in the `test_resequence_output` directory. 

# Input files

- reference genomic sequence file(fasta format)
- reference annotation file(gff3 format)
- **cleaned** paired-end short read sequence files(fastq format)

# Result files
After the pipeline finished, you would like to first read the following two files: 
- `report/SNV/clustered.SNV.txt`
- `report/SNV/clustered.SNV.html`
- `report/SV/clustered.SV.txt`
- `report/SV/clustered.SV.html`

Each kind of variants corresponds to two output files with the same content but different format(html and txt). To make the result more readable, variants are clustered by start position. 

To further understand the SV details, you may want to read the following files: 
- `delly/sample1/sample1.filtered.vcf`
- `cnv/sample1/sample1.high.xls`
- `cnv/sample1/sample1.del.xls`

To further understand the SNV details, you may want to read the following files:
- `BWA/sample1/sample1.final.vcf`

To view specific variants in genome browser like igv, you may need the following files: 
- `BWA/sample1/sample1.bam`
- `reference/reference.fasta`
- `reference/reference.gff`

# Pipeline detail
1. Align short reads to reference sequence. 
2. Call SNP by samtools + bcftools, and then annotate SNP by annovar. 
3. Call SV by delly + samtools depth. Why use samtools here? SV results in improper-paired reads or splited-reads at the end of SVs. Delly use these reads to identify SVs. However, some kinds of SV can not be identified by delly: repeat-mediated SV. As we knowed, under the recombination system, repeat can result in SVs like deletion and inversion. These SVs are surrounded by repeat, where reads are marked as low quality by aligner for multi-map, thus can not be called by delly. In a compromise, we can identify deletion by read depth, this is why samtools is included in this pipeline. 
