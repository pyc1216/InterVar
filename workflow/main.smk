def get_webfrom(wildcards):
    '''
    Some database have different webfrom
    Return the webfrom parameter based on dbname
    '''
    return 'ucsc' if wildcards.dbname == "rmsk" else 'annovar'

def get_real_output(wildcards):
    '''Map dbname to the actual downloaded file name'''
    special_cases = {
        "1000g2015aug": "hg19_ALL.sites.2015_08.txt",
        "ALL.sites.2015_08": "hg19_ALL.sites.2015_08.txt"  # 兼容旧命名
    }
    return f"humandb/{special_cases.get(wildcards.dbname, f'hg19_{wildcards.dbname}.txt')}"

rule annovar_download:
    '''
    Download annovar humandb
    Refer to the following url: https://gist.github.com/fo40225/f135b50b3e47d0997098264c3d28e590
    Unified download rules: Handle special file names through params and output virtual tag files.
    '''
    output:
        touch("humandb/{dbname}.download.done")  # 虚拟文件标记任务完成
    params:
        webfrom=get_webfrom,
        real_output=get_real_output,  # 实际文件名
        dbdir=lambda w: f"humandb/{w.dbname}"  # 兼容解压到目录的情况（可选）
    log:
        "logs/annovar_download.{dbname}.log"
    singularity:
        config['singularity']['annovar']
    shell:
        # 下载到实际路径（自动创建目录）
        "annotate_variation.pl --downdb --webfrom {params.webfrom} "
        "--buildver hg19 {wildcards.dbname} humandb 1>{log} 2>&1; "
        # 检查文件是否存在（可选验证）
        "[ -f {params.real_output} ] && touch {output} || (echo 'Error: File not found'; exit 1)"

rule download_all:
    '''
    Rely on the virtual tag files of all databases to trigger the download.
    '''
    input:
        expand(
            "humandb/{dbname}.download.done",
            dbname=[
                "esp6500siv2_all", "1000g2015aug", "avsnp147", 
                "dbnsfp42a", "dbscsnv11", "gnomad_genome", 
                "clinvar_20210501", "refGene", "ensGene", 
                "knownGene", "rmsk"
            ]
        )

# 规则1：处理ex1（AVinput格式）
rule run_example_1:
    input:
        "example/ex1.avinput"
    output:
        "example_out/ex1.hg19_multianno.txt.intervar"
    params:
        input_type="AVinput",
        #output_prefix=lambda w, output: str(output[0]).replace('.hg19_multianno.txt.intervar', ''),
        output_prefix='example_out/ex1'
    log:
        "logs/run_example_1.log"
    singularity:
        config['singularity']['annovar']
    shell:
        "python Intervar.py -c config.ini -i {input} --input_type {params.input_type} "
        "-o {params.output_prefix} 1>{log} 2>&1"

# 规则2：处理ex2（多样本VCF）
rule run_example_2:
    input:
        "example/ex2.vcf"
    output:
        "example_out/ex2.TUMOR.hg19_multianno.txt.intervar"
    params:
        input_type="VCF_m",
        output_prefix='example_out/ex2'
    log:
        "logs/run_example_2.log"
    singularity:
        config['singularity']['annovar']
    shell:
        "python Intervar.py -c config.ini -i {input} --input_type {params.input_type} "
        "-o {params.output_prefix} 1>{log} 2>&1"

# 规则3：处理ex3（单样本VCF）
rule run_example_3:
    input:
        "example/ex3.vcf"
    output:
        "example_out/ex3.hg19_multianno.txt.intervar"
    params:
        input_type="VCF",
        output_prefix='example_out/ex3'
    log:
        "logs/run_example_3.log"
    singularity:
        config['singularity']['annovar']
    shell:
        "python Intervar.py -c config.ini -i {input} --input_type {params.input_type} "
        "-o {params.output_prefix} 1>{log} 2>&1"

# 聚合所有结果
rule run_example_all:
    input:
        "example_out/ex1.hg19_multianno.txt.intervar",
        "example_out/ex2.TUMOR.hg19_multianno.txt.intervar",
        "example_out/ex3.hg19_multianno.txt.intervar"