

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
        fastas = "{name}_{region}.fasta"
    output:
        phy = "{name}_{region}.fasta.uniqueseq.phy",
        model = "{name}_{region}.fasta.model.gz",
        treefile = "{name}_{region}.fasta.treefile",
        iqtree = "{name}_{region}.fasta.iqtree",
        bionj = "{name}_{region}.fasta.bionj"
    log:
        "logs/iqtree/{name}_{region}.log"
    params:
        bootstraps = 10000
    threads: 16
    shell:
        """
        iqtree -s {input.fasta} -B {bootstraps} -T {threads} 2> {log}
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