# MA-Scripts

This repository contains configuration files and scripts for a scenario-based
resilience analysis of the Danish energy system using PyPSA-Eur.
The study analyses disruptive events for the years 2025 and 2035.

---

## Prerequisites

- A working installation of **PyPSA-Eur**
- Python environment and dependencies as defined by PyPSA-Eur
- `snakemake`

This project is intended to be run within or alongside a PyPSA-Eur setup.

---

## Usage

### 1. Generate base networks (2025 and 2035)

Use the provided PyPSA-Eur configuration files to generate the base networks:

- One configuration for the 2025 system
- One configuration for the 2035 sector-coupled system

Make sure paths, filenames, and scenario names in the config files match your
local setup and the location where the base networks are stored.

---

### 2. Run scenario simulations

Execute the relevant Snakemake files to run all disruptive event scenarios based
on the generated base networks.

Each Snakemake workflow will:

- Apply the defined disruptions
- Solve the network
- Generate and save plots and intermediate results in the corresponding
  scenario folders

---

### 3. Post-processing and metrics

After all scenarios have been executed, aggregated resilience metrics and summary
results are automatically computed and stored as CSV files for further analysis.

---

**Note:**  
Paths, config filenames, and network locations may need to be updated depending
on how PyPSA-Eur is installed and where outputs are saved.
