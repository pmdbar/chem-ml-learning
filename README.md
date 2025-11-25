# ChemML Learning Project

This repository contains tools and prototypes developed as part of my cheminformatics and AI/ML engineering training. The focus is on small-molecule drug discovery workflows, including property calculation, molecular descriptors, and RDKit-based utilities.

---

## **Features**
- RDKit-based small-molecule property computation
- Command-line tool (CLI) supporting:
  - **logP** (Wildman–Crippen)
  - **Molecular Weight (MW)**
- CSV → CSV processing pipeline
- Modular Python code in a standard `src/` package layout
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
```bash
python -m chemml.compute_logp     --input molecules.csv     --output molecules_with_props.csv
```

### Output CSV Columns
- smiles  
- name  
- logP  
- mw  

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
