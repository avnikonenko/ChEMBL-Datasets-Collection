# ChEMBL Datasets Collection

![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
![Python: >=3.6](https://img.shields.io/badge/Python-≥3.6-green.svg)

A Python tool to automate retrieval and preprocessing of bioactivity datasets from [ChEMBL](https://www.ebi.ac.uk/chembl/), ready for QSAR/ML modeling.


---

## Features

- Fetch bioactivity data via the **ChEMBL REST API** or a **local SQLite database**
- Classify compounds by activity type (agonist, antagonist, inhibitor, etc.) using assay description keyword matching
- Build **classification** (active/inactive labels) or **regression** (mean pActivity values) datasets
- Convert units automatically (uM, mM → nM) and calculate plog values
- Remove controversial data (e.g., compounds active as both agonist and antagonist)
- Batch-process multiple targets in parallel with literature-recommended thresholds

---

## Installation

```bash
conda create -n dataset-collect python=3.10
conda activate dataset-collect
conda install -c conda-forge rdkit pandas requests numpy
```

---

## Quick Start

**Single target — regression dataset:**
```bash
python Dataset_collection.py -i CHEMBL301 --t reg --activity Ki IC50 -v plog_value -a inhibitor
```

**Multiple targets from a file:**
```bash
cat example/chembl_ids | xargs -I {} python Dataset_collection.py -i {} --t reg --activity EC50 IC50 -v plog_value -a inhibitor
```

**Batch processing with recommended thresholds (parallel):**
```bash
python dataset_script.py -i example/chembl29_targets_test.csv -n 4
```

See [example/README.md](example/README.md) for a description of all example input files and [example/run_examples.sh](example/run_examples.sh) for ready-to-run commands.

---

## Scripts

| Script | Purpose |
|---|---|
| `Dataset_collection.py` | Single target: fetch, classify, preprocess, save |
| `dataset_script.py` | Batch targets from CSV with literature-recommended thresholds |

---

## Pipeline

### 1. Fetch activity data
Retrieves from ChEMBL API or local SQLite database:
`ChEMBL ID`, `SMILES`, `activity value`, `units`, `bioactivity type` (Ki, IC50, EC50…), `assay type`, `assay description`, `pChEMBL value`, `operator`, `target`, `reference`

### 2. Classify activity type
Assigns `act` (broad) and `type_act` (specific) columns by keyword matching on `assay_description`.

Supported types:

| `act` | `type_act` examples |
|---|---|
| agonist | full agonist, partial agonist |
| antagonist | competitive antagonist, non-competitive antagonist, inverse agonist |
| inhibitor | competitive inhibitor, non-competitive inhibitor, uncompetitive inhibitor |
| activator | activator |
| allosteric | allosteric |
| blockator | blockator |
| binding | binding |
| other | (no classification) |

> **Order matters:** list more specific types before general ones.
> Use `-a "partial agonist" "agonist"`, not `-a "agonist" "partial agonist"`.

To inspect the full keyword tree:
```bash
python Dataset_collection.py --show_tree_act
```

### 3. Preprocess
- Convert units to nM; calculate plog values
- Add standardized SMILES if provided via `-s`
- Calculate molecular weight (RDKit)

### 4. Label activity (classification) / aggregate (regression)

**Classification:** Labels each experiment as `active`, `inactive`, or `undefined` based on value thresholds and operators.

**Regression:** Groups by (SMILES, bioactivity type, units, act, type_act) and computes mean plog value.

### 5. Remove controversial data
- Compounds active in opposite-type assays (e.g., agonist AND antagonist)
- Classification: compounds with both `active` and `inactive` labels in similar experiments
- Regression: compounds with plog range ≥ 0.5, or non-`=` operators

### 6. Save output
- `{prefix}_{act}_{type}.csv` — dataset per activity type
- `{prefix}_{type}_log.csv` — full log of all data
- `{prefix}_stat.csv` — statistics (active/inactive/undefined counts)

---

## Usage: Dataset_collection.py

```
python Dataset_collection.py -i ChEMBLID --t class/reg --activity Ki IC50 [options]
```

| Argument | Description | Default |
|---|---|---|
| `-i, --input` | Target ChEMBL ID | required |
| `--t` | Dataset type: `class` or `reg` | required |
| `--activity` | Bioactivity types to include (e.g. `Ki IC50 EC50`) | required |
| `-o, --output` | Output path prefix | `{ChEMBLID}/{ChEMBLID}` |
| `-a, --act` | Activity types to extract | `other` |
| `-c, --contr` | Pairs of contradictory activity types to filter | — |
| `-v, --value` | Use `value` (nM) or `plog_value` | `value` |
| `-d, --db` | Path to local ChEMBL SQLite database | API |
| `-s, --smi_std` | Path to standardized SMILES file | — |
| `--active_value` | Threshold value for active label | `6` |
| `--active_op` | Operator for active threshold (`>=`, `<=`, `>`, `<`) | `>=` |
| `--inactive_value` | Threshold value for inactive label | `6` |
| `--inactive_op` | Operator for inactive threshold | `<` |
| `--sep` | CSV separator | `\t` |
| `--show_tree_act` | Print activity keyword tree and exit | — |

---

## Usage: dataset_script.py

```
python dataset_script.py -i targets.csv [options]
```

| Argument | Description | Default |
|---|---|---|
| `-i, --input` | CSV file with target metadata | required |
| `-t, --type_dataset` | Dataset type: `class` or `reg` | `class` |
| `-d, --db` | Local ChEMBL SQLite database | API |
| `-s, --smi_std` | Standardized SMILES file | — |
| `-n, --ncpu` | Number of CPU cores | `1` |
| `--sep` | CSV separator | `\t` |

**Input CSV columns:**

| Column | Required | Description |
|---|---|---|
| `chembl_id` | yes | Target ChEMBL ID |
| `target_class` | no | Protein class for recommended thresholds |
| `cutoff_act` | no | Custom active threshold (overrides `target_class`) |
| `cutoff_inact` | no | Custom inactive threshold |
| `activity` | no | Semicolon-separated bioactivity types, e.g. `Ki;IC50` |
| `category` | no | Semicolon-separated activity types, e.g. `agonist;antagonist` |

**Recommended thresholds by protein class** (Bosc et al., J. Cheminf. 2019, 11, 4):

| Protein class | Active threshold (pKi or pIC50) |
|---|---|
| enzyme | ≥ 7.5 |
| membrane_protein, membrane_receptor | ≥ 7 |
| epigenetic_factor, epigenetic regulator, transporter | ≥ 6 |
| ion_channel | ≥ 5 |

---


## Dependencies

- [RDKit](https://www.rdkit.org/)
- [Pandas](https://pandas.pydata.org/)
- [NumPy](https://numpy.org/)
- [Requests](https://requests.readthedocs.io/)

---

## License

MIT — see [LICENSE](LICENSE)

## Authors
Aleksandra Ivanova, Pavel Polishchuk (2021)
