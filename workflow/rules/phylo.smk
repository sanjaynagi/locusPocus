

# rule write_haplotype_fastas:
#     input:
#         cohorts = ag3_sample_sets
#     output:
#         fastas = "{name}_{region}.fasta"
#         metadata = "{name}_{region}.fasta.metadata.tsv"
#     params:
#         regions = 
#     shell:
#         """

#         """

rule iqtree:
    input:
        fasta = "results/phylo/{name}_{region}.fasta"
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
        iqtree -s {input.fasta} -B {params.bootstraps} -T {threads} 2> {log}
        """


# rule plot_phylo:
#     input:
#         cohorts
#     output:
#         fastas = 
#     params:
#         regions = 
#     shell:
#     """

#     """