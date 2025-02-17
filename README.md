# Prion Formation and Polygon-Based Simulations

## 📌 Example Implementation

All these processes are demonstrated in **`main.ipynb`**, providing a structured example of their execution.

---

## 🧬 Abstract

This repository contains the code for a computational framework that models **prion formation and transformation** using polygon-based structures. Prions are misfolded proteins that propagate their structural arrangement to neighboring proteins, leading to both functional and pathological outcomes. Understanding prion dynamics is essential for biological and synthetic applications, yet computational models for **experimental design, hypothesis testing, and control** remain limited.

This project identifies **key prionic properties** and implements a biologically inspired model using **simple mechanical structures capable of complex conformational changes**. The framework includes tools for generating, analyzing, and validating prion-like behavior through computational simulations. A **prototypical mechanical prion** is designed and experimentally validated, demonstrating the utility of this approach. This repository provides a foundation for **studying and manipulating prionic behavior** in both natural and engineered systems.

---

## 📂 Project Structure

### 🔹 0. Import

Handles all necessary **imports and dependencies** for running the scripts.

### 🔹 1. External Polygons

- Sets up the **database** for polygon structures.
- Initializes **simulation values** and generates random polygons for analysis.

### 🔹 2. Pair Selection for Prion Formation

- Establishes a **structured database** for systematic searching.
- Identifies **healthy-prion pairs** and tracks their structural transformations.
- Uses **`polygon_matching_main.py`** for polygon matching and database operations.

### 🔹 3. Internal Nodes

- **`internal_node_main.py`** manages **generation, comparison, and validation** of internal node configurations.
- Integrates **MATLAB** for **optimized internal node placement** while ensuring connectivity constraints.

### 🔹 4. NEB Computation

#### 🏗 **Nudged Elastic Band (NEB) Simulation for Polygon Transitions**

This module performs **NEB simulations** to study the transition of a **healthy polygon** into a **prion-like structure** under physical constraints. The goal is to compute:

✔ **Reaction pathways**\
✔ **Binding energies**\
✔ **Morphological transformations** between polygonal structures.

### 🔹 5. Dynamical Simulation

#### 🔬 **Simulation and Protein Binding **

This script models **protein binding and molecular interactions** using **LAMMPS**. It includes:

✔ **Molecule insertion**\
✔ **Energy minimization**\
✔ **Dynamic simulations**\
✔ **Mathematical utilities for evaluating transformations**

### 🔹 6. Markov Analysis

- Implements **Markov analysis** for statistical modeling of **polygon transitions**.
- Enables **state change predictions** within the system.

---

📢 **For inquiries and contributions, feel free to reach me at ouellet@seas.upenn.edu** 🚀

