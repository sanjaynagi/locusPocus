


def getBamURL(sampleID):

    ena_df = pd.read_csv("results/alignments/ena_alignments.csv")
    ERZ = ena_df.query("sample_id == @sampleID")['ena_analysis'].to_list()[0]
    bamPath = f"ftp://ftp.sra.ebi.ac.uk/vol1/{ERZ[:6]}/{ERZ}/{sampleID}.bam"

    return(bamPath)






rule downloadBams:
    output:
        temp("results/alignments/{sampleID}.bam")
    log:
        dl="logs/downloadBams/{sampleID}.log",
    params:
        bamPath = getBamURL,
    shell:
        """
        wget {params.bamPath} -O {output} 2> {log.dl}
        """

rule indexBams:
     input:
        "results/alignments/{sampleID}.bam"
     output:
        "results/alignments/{sampleID}.bam.bai"
     log:
        "logs/index_bams/{sampleID}_index.log"
     shell:
        "samtools index {input} {output} 2> {log}"


rule subsetBams:
    input:
        "results/alignments/{sampleID}.bam",
        "results/alignments/{sampleID}.bam.bai"
    output:
        "results/alignments_{name}/{sampleID}.bam",
        "results/alignments_{name}/{sampleID}.bam.bai"
    log:
        "logs/subset_bams/{sampleID}.{name}.log"
    params:
        region = config['Locus']['region']
    shell:
        """
        samtools view {input} {params.region} -bS -o {output} 2> {log}
        samtools index {output} 2>> {log}
        """


rule haplotype_caller_gvcf:
    input:
        bams = expand("results/alignments_{name}/{sampleID}.bam", name=config['Locus']['name'], sampleID = metadata['sample_id']),
        bais = expand("results/alignments_{name}/{sampleID}.bam.bai" , name=config['Locus']['name'], sampleID = metadata['sample_id']),
        ref = config['reference']
    output:
        gvcf = "results/indels/gvcfs/{sampleID}.g.vcf"
    log:
        "logs/gatk/haplotypeCaller/{sampleID}.log"
    params:
        extra="",  # optional
        java_opts="",  # optional
    wrapper:
        "v1.7.0/bio/gatk/haplotypecaller"



rule genomics_db_import:
    input:
        gvcfs = expand("results/indels/gvcfs/{sampleID}.g.vcf", sampleID = metadata['sample_id'])
    output:
        db=directory("results/indels/db"),
    log:
        "logs/gatk/genomicsdbimport.log",
    params:
        intervals="ref",
        db_action="create",  # optional
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.7.0/bio/gatk/genomicsdbimport"


rule genotype_gvcfs:
    input:
        genomicsdb="results/indels/db",  # combined gvcf over multiple samples
        ref=config['reference']
    output:
        vcf="results/indels/calls.{name}.vcf",
    log:
        "logs/gatk/genotypegvcfs.{name}.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    resources:
        mem_mb=1024
    wrapper:
        "v1.7.0/bio/gatk/genotypegvcfs"


rule gatk_select:
    input:
        vcf="results/indels/calls.{name}.vcf",
        ref=config['reference']
    output:
        vcf="results/indels/indels.{name}.vcf",
    log:
        "logs/gatk/select/indels_{name}.log",
    params:
        extra="--select-type-to-include INDEL",  # optional filter arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=1024,
    wrapper:
        "v1.7.0/bio/gatk/selectvariants"



# rule prepMultiAlleles_indels_cnvs:
#     input:
#     output:
#     log:
#     shell: