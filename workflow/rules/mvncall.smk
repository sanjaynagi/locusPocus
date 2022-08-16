



rule write_haps_for_mvncall:
    input:
        nb = "workflow/notebooks/haps_for_mvncall.ipynb",
    output:
        nb = "results/notebooks/haps_for_mvncall_{name}.ipynb",
        html = "results/notebooks/haps_for_mvncall_{name}.html",
        vcf = "results/phasing/{name}.vcf",
        scaffold = "results/phasing/{name}.haps",
        sample_file = "results/phasing/{name}.sample"
    log:
        "logs/haps_for_mvncall/{name}.log"
    params:
        dataset = config['LocusOfInterest']['name'],
        contig = config['LocusOfInterest']['contig'],
        ag3_sample_sets = ag3_sample_sets,
        start = config['LocusOfInterest']['start'],
        end = config['LocusOfInterest']['end']
    shell:
        """
        papermill {input.nb} {output.nb} -k locusPocus -p dataset {params.dataset} -p cohorts {params.ag3_sample_sets}
        -p contig {params.contig} -p start {params.start} -p end {params.end}

        python -m nbconvert {output.nb} --to html --stdout --no-input \
             --ExecutePreprocessor.kernel_name=locusPocus > {output.html} 
        """

rule mvncall:
    input:
        vcf = "results/phasing/{name}.vcf",
        scaffold = "results/phasing/{name}.haps",
        sample_file = "results/phasing/haps.sample",
    output:
        phasedVCF = "results/{name}_phasedMulti.vcf",
    log:
        "logs/mvncall/{name}.log"
    threads: 16
    params:
        start = config['LocusOfInterest']['start'],
        end = config['LocusOfInterest']['end'],
    shell:
        """
        /home/sanj/mvncall_v1.0_x86_64_dynamic/mvncall  --sample-file {input.sample_file} \
        --glfs {input.vcf} --scaffold-file {input.scaffold} --o {output.phasedVCF} --int {params.start} {params.end} 2> {log}
        """
