# DGE Tool Choice Benchmark for MPS Transcriptomics

This repository provides reproducible analysis workflows for comparing two widely used RNA-seq differential gene expression tools, **edgeR** and **DESeq2**, in the context of transcriptomic model evaluation, cross-study robustness, and microphysiological systems-oriented open science.

The repository accompanies the preprint:

**Rezapour, Mostafa. _Tool Choice Matters: Evaluating edgeR vs. DESeq2 for Sensitivity, Robustness, and Cross-Study Performance._ arXiv preprint arXiv:2601.04122, 2026.**

## Why this repository matters

Microphysiological systems, organ-on-chip models, engineered tissues, and related human-relevant platforms increasingly rely on transcriptomic profiling to evaluate biological fidelity, disease modeling, injury response, drug response, and cross-platform reproducibility.

In these settings, computational analysis is not a neutral final step. The choice of differential gene expression tool can influence which genes are identified, which pathways appear enriched, and which biological conclusions are carried forward.

This repository provides a transparent framework for evaluating how edgeR and DESeq2 differ across sensitivity, robustness, pathway-level interpretation, and cross-study performance.

## Connection to PhysioVerse-OSE

PhysioVerse-OSE aims to support open, reproducible, and standardized evaluation of microphysiological and tissue-engineered systems. This repository contributes to that mission by providing a benchmark workflow for transcriptomic analysis decisions.

In the OSE context, this tool helps users ask:

- Do edgeR and DESeq2 identify the same biological signals?
- Which genes are shared across tools, and which are tool-specific?
- Are tool-specific gene sets useful for downstream classification?
- Do results generalize across independent studies?
- How should transcriptomic evidence be interpreted when benchmarking MPS fidelity?

## Repository contents

The repository is organized around the manuscript figures and reproducibility workflows.

- R scripts for differential gene expression analysis
- Python scripts for data summarization, comparison, and visualization
- Notebooks used for figure generation
- Final figure outputs in PNG and PDF format
- Documentation for reproducibility and reuse

## What is included

This repository includes code and figure-generation workflows.

## What is not included

Raw RNA-seq count matrices, large processed datasets, and sensitive sample-level metadata are not included.

Users should obtain public datasets from their original sources and format them according to the instructions provided in the documentation.

## Study scope

The accompanying analysis compares edgeR and DESeq2 across publicly available bulk RNA-seq datasets covering infectious and fibrotic disease contexts. The workflows examine:

1. Sensitivity to sample size
2. Robustness to outlier perturbation
3. Concordance and divergence of differentially expressed genes
4. Pathway-level interpretation
5. Within-dataset classification using tool-specific gene sets
6. Cross-study generalization across independent SARS-CoV-2 datasets

## Intended users

This repository is intended for:

- MPS and organ-on-chip investigators using RNA-seq for model evaluation
- Computational biologists comparing DGE workflows
- Open science teams developing reproducible benchmark pipelines
- Reviewers and collaborators evaluating transcriptomic robustness
- PhysioVerse-OSE users interested in standardized biological model assessment

## Citation

If you use this repository, please cite:

Rezapour, Mostafa. _Tool Choice Matters: Evaluating edgeR vs. DESeq2 for Sensitivity, Robustness, and Cross-Study Performance._ arXiv preprint arXiv:2601.04122, 2026.

## Status

This repository is intended as a reproducible companion repository for the manuscript and as an OSE-aligned transcriptomics benchmark tool. The current release focuses on figure reproducibility and transparent documentation. Future versions may include command-line wrappers, standardized input validators, and automated benchmark reports for MPS transcriptomic datasets.

## License

Please see the LICENSE file for reuse terms.
