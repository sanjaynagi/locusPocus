



rule write_haps_for_mvncall:
    input:
        "resources/.kernel.set",
    output:
        nb = "results/notebooks/haps_for_mvncall_{name}.ipynb",
        html = "results/notebooks/haps_for_mvncall_{name}.html",
        vcf = "results/phasing/{name}.genotypes.vcf",
        scaffold = "results/phasing/{name}.haps",
        sample_file = "results/phasing/{name}.sample"
    log:
        "logs/haps_for_mvncall/{name}.log"
    params:
        nb = "workflow/notebooks/haps_for_mvncall.ipynb",
        dataset = config['Locus']['name'],
        contig = config['Locus']['contig'],
        ag3_sample_sets = lambda set: f"'{{cohorts:{ag3_sample_sets}}}'",
        start = config['Locus']['start'],
        end = config['Locus']['end']
    shell:
        """
        papermill {params.nb} {output.nb} -k locusPocus -p dataset {params.dataset} -y {params.ag3_sample_sets} \
        -p contig {params.contig} -p start {params.start} -p end {params.end}

        python -m nbconvert {output.nb} --to html --stdout --no-input \
             --ExecutePreprocessor.kernel_name=locusPocus > {output.html}
        """











rule mvncall:
    input:
        vcf = "results/phasing/{name}.genotypes.vcf",
        scaffold = "results/phasing/{name}.haps",
        sample_file = "results/phasing/{name}.sample",
    output:
        phasedVCF = "results/phasing/{name}.phasedMulti.vcf",
    log:
        "logs/mvncall/{name}.log"
    threads: 16
    params:
        start = config['Locus']['start'],
        end = config['Locus']['end'],
    shell:
        """
        /home/sanj/apps/mvncall_v1.0_x86_64_dynamic/mvncall --sample-file {input.sample_file} \
        --glfs {input.vcf} --scaffold-file {input.scaffold} --o {output.phasedVCF} --int {params.start} {params.end} --numsnps 80 --lambda 0.1 2> {log}
        """



### split multialleles
#bcftools norm results/phasing/coeae1f.genotypes.vcf -m - snps > results/phasing/coeae1f.genotypes.split.vcf