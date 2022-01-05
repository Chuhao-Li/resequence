########
#config#
########
import os

base_dir = os.path.abspath(config['BASE_DIR'])
vcfutils = base_dir + "/vcfutils.pl"
PICARD = base_dir + '/picard.jar'
delly = base_dir + '/delly'
python_path = 'python3'
script_path = base_dir + '/script'
annovar_path = base_dir + '/annovar'
os.environ['PATH'] = script_path + ':' + os.environ['PATH']


# 样品列表转化成字典，方便取值。
sample_list = config['samples']
with open('samples.list', 'w') as f:
    for s in sample_list:
        print(s, file=f)

sample_dict = {}
for name, read1, read2 in sample_list:
    read1 = os.path.abspath(read1)
    read2 = os.path.abspath(read2)
    sample_dict[name] = [read1, read2]
config['samples'] = sample_dict
config["reference_gff"] = os.path.abspath(config["reference_gff"])
config["reference_fasta"] = os.path.abspath(config["reference_fasta"])
config["out_dir"] = os.path.abspath(config["out_dir"])

reference_dir = "{}/reference".format(config["out_dir"])
if not os.path.exists(reference_dir):
    os.system("mkdir -p {}/reference".format(reference_dir))
    os.system("cp {} {}/reference.gff".format(config["reference_gff"], reference_dir))
    os.system("cp {} {}/reference.fasta".format(config["reference_fasta"], reference_dir))

# print(config)
goals = ['qc/multiqc_report.html']

if config['SV'] == "true":
    # for i in config['samples']:
    #     goals.append("delly/{sample}/{sample}.filtered.filtered.vcf".format(sample=i))
    goals.append("report/SV/clustered.SV.txt")

if config['SNV'] == "true":
    # for i in config['samples']:
    #     goals.append("snv_anno/{sample}/{sample}_snp.exonic_variant_function".format(sample=i))
    goals.append("report/SNV/clustered.SNV.txt")
 
workdir: config['out_dir']

########
# run  #
########

rule all:
    input:
        goals


def qc_in(x):
    return {'fq1': config['samples'][x.sample][0], 
            'fq2': config['samples'][x.sample][1]}

rule quality_control:
    # 自动检测感觉还是不靠谱。
    input:
        [j for i in config['samples'].values() for j in i]
    output:
        qc = 'qc/multiqc_report.html'
    run:
        shell('if [ ! -d qc ]; then mkdir qc; fi')
        threads = len(config['samples']) * 2
        shell('fastqc -t {t} -o qc {fq}'.format(t=str(threads), fq=' '.join([j for i in config['samples'].values() for j in i])))
        shell('if [ -d qc/multiqc_data ]; then rm -r qc/multiqc_data; fi')
        shell('if [ -f qc/multiqc_report.html]; then rm qc/multiqc_report.html; fi')
        shell('multiqc -o qc qc')

# snakemake里面的语法错误不容易debug……运行之后不像普通的python脚本那样告诉你哪里有问题。
def bwa_in(x):
    return {'fq1': config['samples'][x.sample][0], 
            'fq2': config['samples'][x.sample][1],
            'fa' : config['reference_fasta']}
            # 'qc' : 'qc/{}.qc.ok'.format(x.sample)}

rule bwa:
    input:
        unpack(bwa_in)
    output:
        bam="BWA/{sample}/{sample}.bam",
        stats="BWA/{sample}/{sample}.stats"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index {input.fa} && "
        "if [ -f {output.bam} ]; then rm {output.bam}.* ;fi && "
        "bwa mem -t 30 {input.fa} {input.fq1} {input.fq2} |samtools sort -o {output.bam} && "
        "samtools stats -@ 20 {output.bam} >{output.stats} && "
        "samtools index {output.bam}"
# view map rate: for i in `ls BWA/`; do awk -F '\t' -v sample="${i}" '$2=="raw total sequences:"{total=$3} $2=="reads mapped:"{mapped=$3} END{print sample, mapped/total}' OFS='\t' BWA/${i}/*.stats; done

rule deletions:
    input:
        'BWA/{s}/{s}.bam'.format(s=i) for i in config['samples']
    output:
        ['cnv/{sample}/{sample}.del.xls'.format(sample=i) for i in config['samples']]
    run:
        temp = []
        for s in config['samples']:
            shell('samtools depth -a BWA/{s}/{s}.bam >cnv/{s}/{s}.depth'.format(s=s))
            shell('{python_path} {script_path}/copy_number_variant.py cnv/{s}/{s}.depth cnv/{s}/{s}.del.xls cnv/{s}/{s}.high.xls > cnv/{s}/{s}.stat'.format(python_path=python_path, script_path=script_path, s=s))
            # temp.append('cnv/{s}/{s}.del.xls'.format(s=s))
        # shell('cluster_SV.py ' + ' '.join(temp) + ' >cnv/cluster.del.xls')

rule delly:
    #find sv between 00_lasI and 01_lasI with wildtype
    #require: samtools, java, delly, picard
    input:
        bam=rules.bwa.output.bam,
        fa=config['reference_fasta']
    output:
        filtered_bam=temp("delly/{sample}/{sample}.sorted.filtered.bam"),
        markDup=temp("delly/{sample}/{sample}.markDuplicates.bam"),
        mdm=temp("delly/{sample}/marked_dup_metrics.txt"),
        bai=temp("delly/{sample}/{sample}.markDuplicates.bam.bai"),
        bcf="delly/{sample}/{sample}.bcf",
        filtered_vcf="delly/{sample}/{sample}.filtered.vcf"
    shell:
        "samtools view -bF 4 {input.bam} > {output.filtered_bam} && "
        "java -jar {PICARD} MarkDuplicates I={output.filtered_bam} O={output.markDup} M={output.mdm} && "
        "samtools index {output.markDup} && "
        "{delly} call -o {output.bcf} -g {input.fa} {output.markDup} && "
        "bcftools filter -i 'PE > 10' {output.bcf} -o {output.filtered_vcf}"

rule call_snp:
    # snp between delta_lasI(00 and 01) and wildtype
    #export PATH=/bioinfo_tools/src/bcftools/misc/:$PATH
    #require: samtools, bcftools, vcfutils.pl
    input:
        bam=rules.bwa.output.bam,
        fa=config['reference_fasta']
    output:
        raw_bcf=temp("BWA/{sample}/{sample}.raw.bcf"),
        bcf=temp("BWA/{sample}/{sample}.bcf"),
        final_vcf="BWA/{sample}/{sample}.final.vcf"
    conda:
        "envs/call_snp.yaml",
    shell:
        "samtools faidx {input.fa} && "
        "samtools mpileup -g -f {input.fa} {input.bam} > {output.raw_bcf} && "
        "bcftools call -c --ploidy 1 {output.raw_bcf} > {output.bcf} && "
        "bcftools view {output.bcf} | {vcfutils} varFilter - > {output.final_vcf}"

rule snv_anno:
    #export PATH=/mnt/sdb1/home/lch/.local/lib/annovar/:$PATH
    #require: gffread gtfToGenePred annovar 
    input:
        vcf_snp=rules.call_snp.output.final_vcf,
        fa = config['reference_fasta'], 
        gff = config['reference_gff']
    output:
        genepred=temp("snv_anno/{sample}/db/{sample}_refGene.txt"),
        refGeneMrna=temp("snv_anno/{sample}/db/{sample}_refGeneMrna.fa"),
        avinput_snp=temp("snv_anno/{sample}/{sample}_snp.avinput"),
        func_anno_snp="snv_anno/{sample}/{sample}_snp.exonic_variant_function",
    shell:
        "sed -i 's/EC_number/ec_number/g' {input.gff} && "
        "{annovar_path}/gff3ToGenePred {input.gff} {output.genepred} && "
        "{annovar_path}/retrieve_seq_from_fasta.pl --format refGene --seqfile {input.fa} {output.genepred} --out {output.refGeneMrna} && "
        "{annovar_path}/convert2annovar.pl --format vcf4 {input.vcf_snp} >{output.avinput_snp} && "
        "{annovar_path}/annotate_variation.pl -buildver {wildcards.sample} -out snv_anno/{wildcards.sample}/{wildcards.sample}_snp {output.avinput_snp} snv_anno/{wildcards.sample}/db/"

#
rule sv_filter:
    input:
        rules.delly.output.filtered_vcf
    output:
        "delly/{sample}/{sample}.filtered.filtered.vcf"
    run:
        shell('/home/lch/Project/2020-07-09_test_resequence_pipeline/script/delly_filter.py {input} >{output}')

rule cluster_SV:
    input:
        ['delly/{s}/{s}.filtered.vcf'.format(s=s) for s in config['samples']] + \
                ['cnv/{s}/{s}.del.xls'.format(s=s) for s in config['samples']] + \
                ['cnv/{s}/{s}.high.xls'.format(s=s) for s in config['samples']]
    output:
        'report/SV/clustered.SV.txt', 
        'report/SV/clustered.SV.html'
    run:
        shell('{python_path} {script_path}/special_SV.py {x}'.format(x=' '.join([i for i in input]), python_path=python_path, script_path=script_path))
        shell('{python_path} {script_path}/create_report.py SV'.format(python_path=python_path, script_path=script_path))

rule cluster_SNV:
    input:
        ['BWA/{s}/{s}.final.vcf'.format(s=s) for s in config['samples']] + \
                ['snv_anno/{s}/{s}_snp.exonic_variant_function'.format(s=s) for s in config['samples']]
    output:
        'report/SNV/clustered.SNV.txt',
        'report/SNV/clustered.SNV.html'
    run:
        shell('{python_path} {script_path}/special_SNV.py'.format(python_path=python_path, script_path=script_path))
        shell('{python_path} {script_path}/create_report.py SNV'.format(python_path=python_path, script_path=script_path))

rule sv_anno:
    input:
        vcf_sv=rules.delly.output.filtered_vcf,
        fa = config['reference_fasta'], 
        gff = config['reference_gff']
    output:
        genepred=temp("sv_anno/{sample}/db/{sample}_refGene.txt"),
        refGeneMrna=temp("sv_anno/{sample}/db/{sample}_refGeneMrna.fa"),
        avinput_sv=temp("sv_anno/{sample}/{sample}_sv.avinput"),
        func_anno_sv="sv_anno/{sample}/{sample}_sv.exonic_variant_function"
    shell:
        "sed -i 's/EC_number/ec_number/g' {input.gff} && "
        "{annovar_path}/gff3ToGenePred {input.gff} {output.genepred} && "
        "{annovar_path}/retrieve_seq_from_fasta.pl --format refGene --seqfile {input.fa} {output.genepred} --out {output.refGeneMrna} && "
        "{annovar_path}/convert2annovar.pl --format vcf4 {input.vcf_sv} >{output.avinput_sv} && "
        "{annovar_path}/annotate_variation.pl -buildver {wildcards.sample} -out sv_anno/{wildcards.sample}/{wildcards.sample}_sv {output.avinput_sv} sv_anno/{wildcards.sample}/db/"
