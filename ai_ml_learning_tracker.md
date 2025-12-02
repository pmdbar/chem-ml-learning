# AI/ML Learning Tracker

This document organizes your 24-week curriculum, checklists, reflections, terminal cheat sheets, and resources. Use it alongside our ongoing chat for Q&A, debugging, and weekly guidance.

---

## **1. Curriculum Overview (24 Weeks)**

A high-level roadmap of your progression from Python foundations → chem-ML → deep learning → capstone project.

* **Weeks 1–4:** Python, Git, testing, CLI tools, data wrangling  
* **Weeks 5–8:** ML basics, chemical featurization, baseline models  
* **Weeks 9–16:** Chemistry-focused ML, pipelines, tuning, interpretability  
* **Weeks 17–24:** Deep learning (PyTorch/GNNs), generative models, capstone project

---

## **2. Weekly Checklists**

### **Week 1 – Environment + Python Reboot**

*

*Notes: Week 1 completed. Environment, tooling, CLI, tests, and GitHub all set up.*

---

### **Week 2 – CLI Enhancements, Refactor, and Project Structure**

* [x] Extend the CLI tool with at least one new feature (added MW)  
* [x] Refactor code into a simple package structure (`src/chemml/` with modules)  
* [x] Update tests to work with the new structure  
* [x] Run pytest successfully after refactor  
* [x] Make at least one new commit + push to GitHub documenting the Week 2 changes

---

## **3. Project & Portfolio Tracker**

### **Project 1: Baseline Molecular Property Predictor**

* Status:  
* Repo link:  
* Next steps:  

### **Project 2: Mini Hackathon Project**

* Status:  
* Repo link:  
* Next steps:  

### **Capstone Project (Weeks 22–24)**

* Idea:  
* Problem Statement:  
* Users:  
* Data Sources:  
* Architecture Plan:  

---

## **4. Terminal & Command Line Cheat Sheet**

### **Conda Basics**

* `conda create -n chemml python=3.10`  
* `conda activate chemml`  
* `conda install rdkit -c conda-forge`

### **Git Basics**

* `git clone <url>`  
* `git status`  
* `git add .`  
* `git commit -m "message"`  
* `git push`  
* `git pull`

### **Python Execution**

* `python script.py`  
* `python -m module.submodule`

---

## **5. Python Patterns & Snippets**

### **Reading a CSV**
```python
import pandas as pd
df = pd.read_csv("file.csv")
```

### **RDKit Molecule from SMILES**
```python
from rdkit import Chem
mol = Chem.MolFromSmiles("c1ccccc1")
```

### **Function Template**
```python
def compute_property(smiles: str) -> float:
    # TODO: implement
    pass
```

---

## **6. Reflection Section**

### **Week 1 Reflection**

* What I learned: I learned about the importance of GIT commits and the relationship with GitHub  
* What was challenging: Staying focused on the bigger picture while dealing with lots of python code.  
* What I want to improve next week: I would like to learn more about python code structure  

---

### **Week 2 Reflection**

**What I learned**  
This week I learned the importance of clean project structure, modularity, and—most importantly—testing as a safety net. Separating chemistry logic into a dedicated module (`chemml.properties`) and isolating CLI responsibilities made the codebase more maintainable. When I ran pytest after the refactor, it immediately revealed issues I wouldn’t have caught easily on my own, including typos, name mismatches, and import problems. Seeing tests fail in specific, helpful ways made the value of testing “click.”

**A vignette from the week**  
After reorganizing the package, pytest reported one passing test and two failing ones. Instead of derailing progress, those failures precisely guided me to the broken spots. It turned debugging into a straightforward, almost mechanical process. For the first time, I experienced tests not as overhead, but as engineering tools—like guardrails that made refactoring safe.

**What was challenging**  
Adapting to a true package structure and managing imports across modules took deliberate effort. Transitioning from Nano to VS Code also involved a learning curve, but ultimately made coding smoother.

**What I want to improve next week**  
Next week I want to deepen my understanding of Python module design, write cleaner ML pipelines, and continue building the testing habits that support scalable engineering.

---

## **7. Resources Library**

### Python & Software Engineering
*

### Machine Learning
*

### Cheminformatics
*

### Deep Learning / GNNs
*

### Deployment / Tools
*

---

