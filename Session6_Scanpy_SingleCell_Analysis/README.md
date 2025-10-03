# Session 6: Single-Cell Analysis with Python/Scanpy

Materials for Session 6 for single-cell analysis with Scanpy, led by Dr. Dominik Otto, Postdoctoral Research Fellow in the Setty Lab at Fred Hutch Cancer Center.

## Overview

This session introduces single-cell RNA-seq analysis using Python's **scverse** ecosystem, with a focus on **Scanpy** and **AnnData**. Students will learn the foundational concepts of the Python single-cell framework and complete a hands-on PBMC analysis workflow.

## Session Structure

### Pre-Session Setup (~2 minutes)
Students should complete the setup before the session, or can catch up quickly at the start:
- Follow [`Setup.md`](Setup.md) to install `uv` and create the Python environment
- Install required packages: `scanpy`, `palantir`, `matplotlib`

### Part 1: Framework Introduction (Notebook 01)
**Instructor-led demonstration** covering:
- **AnnData** data structure and design philosophy
- **Scanpy** ecosystem and core functionality
- **scverse** community and tool landscape (scvi-tools, muon, cellrank, palantir, etc.)
- **Standardization efforts**: best practices, data formats, metadata ontologies
- **Python vs R**: similarities, differences, and interoperability (rpy2/anndata2ri examples)

**Notebook**: [`01_Introduction_to_Python_scverse.ipynb`](01_Introduction_to_Python_scverse.ipynb)

### Part 2: Hands-On Analysis (Notebook 02)
**Interactive practical session** with a complete PBMC workflow:
1. **Quality Control** - Filter cells by UMI counts and mitochondrial content
2. **Normalization** - Library size normalization and log-transformation
3. **Feature Selection** - Highly variable gene identification
4. **Dimensionality Reduction** - PCA, UMAP, diffusion maps
5. **Clustering** - Leiden algorithm for cell grouping
6. **Doublet Detection** - Scrublet on raw counts (with layer swapping)
7. **Imputation** - MAGIC for visualization
8. **Cell Type Annotation** - Marker-based annotation

**Notebook**: [`02_PBMC_RNA_Analysis.ipynb`](02_PBMC_RNA_Analysis.ipynb)

## Key Learning Objectives

By the end of this session, students will:
1. Understand the **AnnData** structure and why it's central to Python single-cell analysis
2. Be familiar with the **scverse** ecosystem and available tools
3. Know how to perform a **complete scRNA-seq analysis** from raw data to cell type annotation
4. Understand **best practices** and community standards
5. Be able to **interoperate between Python and R** when needed

## Resources

### Official Documentation
- [scverse.org](https://scverse.org/) - Central hub for all scverse projects
- [Scanpy documentation](https://scanpy.readthedocs.io/)
- [AnnData documentation](https://anndata.readthedocs.io/)
- [Palantir documentation](https://palantir.readthedocs.io/)

### Best Practices
- [Single-cell best practices book](https://www.sc-best-practices.org/) - Community-driven guide

### Additional Tutorials
- [Scanpy tutorials](https://scanpy-tutorials.readthedocs.io/)
- [Single-cell primers](https://github.com/settylab/single-cell-primers) - More example workflows

### Getting Help
- [scverse Discourse](https://discourse.scverse.org/) - Community forum
- GitHub issues for individual tools
- Stack Overflow with tags: `scanpy`, `single-cell`

## Notes

- The tutorial uses **PBMC data** (10X Genomics) which downloads automatically
- The workflow emphasizes **best practices** and links to resources throughout
