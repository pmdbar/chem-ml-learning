# ChemML Learning Project

This repository contains tools and prototypes developed as part of my cheminformatics and AI/ML engineering training. The focus is on small-molecule drug discovery workflows, including property calculation, molecular descriptors, and RDKit-based utilities.

---

## **Features**
- RDKit-based small-molecule property computation
- Command-line tool (CLI) supporting:
  - **Basic Properties** (Wildman-Crippen logP, MW, Canonical SMILES) 
  - **Extended descriptors** (TPSA,  HBS,  HBA, rotatable bonds)
- CSV → CSV processing pipeline
- Modular Python code in a standard `src/` package layout
- Editable install ('pip install -e .') so ''chemml' can be imported anywhere
- Automated tests (pytest)
- Fully version-controlled with Git + GitHub

---

## **Installation**
```bash
conda create -n chemml python=3.10
conda activate chemml
conda install numpy pandas rdkit -c conda-forge
```
---

## **Usage**
### Basic Properties only
```bash
python -m chemml.compute_logp \
    --input molecules.csv \
    --output molecules_basic.csv \
    --feature-set basic
```
### Extended Properties
```bash
python -m chemml.compute_logp \
    --input molecules.csv \
    --output molecules_extended.csv \
    --feature-set extended
```
### Output CSV Columns
For feature-set basic:
- SMILES  
- name  
- logP  
- mw  
- canonical_SMILES

For feature-set extended:
- SMILES  
- name  
- logP  
- mw  
- canonical_SMILES
- tpsa
- hbd
- hba
- rotatable_bonds
---

## **Project Purpose**
This project serves as an evolving sandbox for:
- Python engineering fundamentals
- Package structure & CLI design
- Testing and reproducibility
- Cheminformatics workflows
- Property calculation pipelines

As the AI/ML curriculum progresses, this repository will expand with new features, modules, and experiments.

---
