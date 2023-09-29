# [Single-cell transcriptomic profiling of human pancreatic islets reveals genes responsive to glucose exposure over 24 hours](https://www.biorxiv.org/content/10.1101/2023.06.06.543931.abstract)

Code used in [Single-cell transcriptomic profiling of human pancreatic islets reveals genes responsive to glucose exposure over 24 hours](https://www.biorxiv.org/content/10.1101/2023.06.06.543931.abstract). 

## Preprocessing and quality control procedures

1. [CellRanger v3.1.0](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger#workflows)
2. [Decontex workflow](https://github.com/CollinsLabBioComp/snakemake-decontx/tree/2d59ca7960c4712091db1dae6440ef54ec3d05c6)
3. [Quality control procedures and cell type labelling workflow](https://github.com/wtsi-hgi/nf_scrna_qc/tree/main)

## Analysis

4. Time interpolation: `scripts/time_interpolation.py`
5. [Differential gene expression workflow](https://github.com/CollinsLabBioComp/nextflow-sc_dge/tree/802a91345513d06acc25e1f010968f9e6229e187)
6. [PoPs analysis workflow](https://github.com/CollinsLabBioComp/snakemake-polygenic_priority_score/tree/3603383bd57abd30f8fdd0cfeeffcbae37ef9729)


