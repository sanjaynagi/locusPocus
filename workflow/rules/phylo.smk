rule write_haplotype_fastas:
    input:
        "resources/.kernel.set",
        nb = "workflow/notebooks/write_phylo_hap_fastas.ipynb",
        outgroups = config['phylo']['outgroup']['alleles'],
        phase1_pos = lambda w : config['phylo']['outgroup']['pos'].format(contig=config['Locus']['contig']),
        regions_yaml = "resources/phylo_regions.yml"
    output:
        nb = "results/notebooks/write_phylo_hap_fastas.ipynb",
        html = "results/notebooks/write_phylo_hap_fastas.html",
        fasta = expand("results/phylo/{locus_name}_{region}.fasta", locus_name = config['Locus']['name'], region = config['phylo']['regions'].keys()),
        fasta_metadata = expand("results/phylo/{locus_name}_{region}.metadata.tsv",locus_name = config['Locus']['name'], region = config['phylo']['regions'].keys())
    log:
        "logs/write_phylo_hap_fastas.log"
    params:
        dataset = config['Locus']['name'],
        contig = config['Locus']['contig'],
        ag3_sample_sets = lambda set: f"'cohorts':{ag3_sample_sets}",
        regions = lambda r: f"'regions': {list(config['phylo']['regions'].values())}",
        region_names = lambda r: f"'region_names': {list(config['phylo']['regions'].keys())}",
        flanking = config['phylo']['flanking'],
        remove_2la_hets = config['phylo']['remove_2la_hets']
    shell:
        """
        papermill {input.nb} {output.nb} -k locusPocus -p locus_name {params.dataset} -y '{{{params.ag3_sample_sets}, {params.regions}, {params.region_names}}}' \
        -p contig {params.contig} -p outgroup_h5_path {input.outgroups} -p phase1_positions_array_path {input.phase1_pos} -p flanking {params.flanking} 2> {log}

        python -m nbconvert {output.nb} --to html --stdout --no-input --ExecutePreprocessor.kernel_name=locusPocus > {output.html}
        """

rule iqtree:
    input:
        fasta = "results/phylo/{name}_{region}.fasta",
    output:
        phy = "results/phylo/{name}_{region}.fasta.uniqueseq.phy",
        model = "results/phylo/{name}_{region}.fasta.model.gz",
        treefile = "results/phylo/{name}_{region}.fasta.treefile",
        iqtree = "results/phylo/{name}_{region}.fasta.iqtree",
        bionj = "results/phylo/{name}_{region}.fasta.bionj"
    log:
        "logs/iqtree/{name}_{region}.log"
    params:
        bootstraps = 10000
    threads: 1
    shell:
        """
        iqtree -s {input.fasta} -B {params.bootstraps} -T {threads} --redo 2> {log}
        """


# rule plot_phylo:
#     input:
#        fasta = "results/phylo/{name}_{region}.fasta",
#        fasta_metadata = "results/phylo/{locus_name}_{region}.metadata.tsv"
#     output:
#         fastas = 
#     params:
#         regions = 
#     shell:
#     """

#     """