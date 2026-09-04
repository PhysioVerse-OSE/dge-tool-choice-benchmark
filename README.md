<div align="center">

<a href="https://physioverse.org/">
  <img src="https://physioverse.org/branding/physioverse-logo.png" alt="PhysioVerse" width="620" />
</a>

# DGE Tool Choice Benchmark

### Evaluating edgeR and DESeq2 across sensitivity, robustness, and cross-study performance

Reproducible analysis workflows supporting transparent comparison of differential gene expression methods across multiple RNA-seq datasets.

<a href="https://physioverse.org/tools"><strong>PhysioVerse Tools</strong></a>
&nbsp;&nbsp;•&nbsp;&nbsp;
<a href="README_reproducibility.md"><strong>Reproducibility Guide</strong></a>
&nbsp;&nbsp;•&nbsp;&nbsp;
<a href="CITATION.cff"><strong>Citation</strong></a>

</div>

---

## Overview

This repository provides the reproducible analysis workflows associated with:

**Rezapour, Mostafa. _Tool Choice Matters: Evaluating edgeR vs. DESeq2 for Sensitivity, Robustness, and Cross-Study Performance._ arXiv preprint arXiv:2601.04122, 2026.**

The benchmark compares **edgeR** and **DESeq2** across RNA-seq datasets representing viral infection, bacterial infection, and fibrotic lung disease.

The workflows evaluate how differential gene expression tool choice can affect detected genes, downstream biological interpretation, robustness, and cross-study predictive performance.

---

## Benchmark Scope

<table>
<tr>
<td width="50%" valign="top">

### Differential Expression Performance

- sample-size sensitivity
- outlier robustness
- DEG concordance and divergence
- tool-specific gene sets

</td>
<td width="50%" valign="top">

### Downstream Biological Performance

- pathway-level concordance
- within-dataset classification
- training-gene selection
- cross-study SARS-CoV-2 prediction

</td>
</tr>
</table>

---

## Why Tool Choice Matters

RNA-seq differential expression analysis is widely used to characterize biological responses and compare experimental systems.

Different analytical methods can identify overlapping but non-identical sets of differentially expressed genes. These differences can propagate into pathway analysis, biomarker selection, predictive modeling, and biological interpretation.

This repository provides a transparent and reproducible framework for examining those differences across multiple datasets and analytical settings.

---

## Repository Contents

The repository is organized around manuscript figures and reproducibility workflows.

```text
figure 1/          Sample-size sensitivity and outlier robustness
figure 2 and 3/    DEG concordance and within-dataset classification
figure 4 and 5/    Pathway-level concordance
figure 6/          Cross-study training-gene selection
figure 7/          Cross-study held-out classification
data/              Data preparation guidance
docs/              Reproduction and contextual documentation
pathways/          Pathway-resource documentation
```

The repository includes R scripts for differential gene expression analysis, Python scripts for summarization and visualization, notebooks used in figure-generation workflows, final figure outputs, reproducibility documentation, and citation metadata.

---

## Data Availability

Raw RNA-seq count matrices and large processed datasets are not bundled with this repository.

The analyses use public datasets from their original repositories. Users should download the required datasets and format them according to the reproducibility documentation.

| Accession / Source | Biological Context |
|---|---|
| `GSE196134` | RSVB stimulation |
| `GSE234118` | Mpox infection |
| `GSE115785` | Ebola virus infection |
| `PMC8202013` | SARS-CoV-2 and influenza |
| `GSE161731` | Bacterial pneumonia, influenza, SARS-CoV-2, and healthy controls |
| `GSE231693` | Idiopathic pulmonary fibrosis and healthy lung |
| `GSE152418` | SARS-CoV-2 PBMC |
| `GSE171110` | SARS-CoV-2 whole blood |

For complete dataset requirements and figure-specific instructions, see [`README_reproducibility.md`](README_reproducibility.md).

---

## Reproduce the Analysis

### [Open the Reproducibility Guide →](README_reproducibility.md)

Additional documentation:

| Resource | Purpose |
|---|---|
| [`docs/figure-reproduction.md`](docs/figure-reproduction.md) | Figure reproduction guidance |
| [`docs/mps-and-ose-context.md`](docs/mps-and-ose-context.md) | Context for MPS and open-science use |
| [`data/README.md`](data/README.md) | Data preparation guidance |
| [`pathways/README.md`](pathways/README.md) | Pathway resources |

---

## Software

The workflows use both **R** and **Python**. Primary R packages include edgeR and DESeq2. Primary Python packages and execution details are documented in [`README_reproducibility.md`](README_reproducibility.md).

---

## Connection to PhysioVerse

This benchmark contributes to the PhysioVerse goal of supporting transparent and reproducible analysis of data generated from advanced in vitro and human-relevant model systems.

Within the broader ecosystem, the repository provides a concrete example of how analytical choices can be evaluated rather than treated as interchangeable.

<p align="center">
  <a href="https://physioverse.org/"><strong>Visit PhysioVerse</strong></a>
  &nbsp;&nbsp;•&nbsp;&nbsp;
  <a href="https://physioverse.org/tools"><strong>Explore Tools</strong></a>
  &nbsp;&nbsp;•&nbsp;&nbsp;
  <a href="https://github.com/PhysioVerse-OSE"><strong>PhysioVerse GitHub</strong></a>
</p>

---

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

If you use this repository or its associated workflows, please cite the associated work described in that file.

---

## Repository Status

This repository serves as the reproducible companion resource for the DGE tool-choice benchmark and provides the code, documentation, and figure-generation workflows needed to inspect and reproduce the analyses.
