# ChEMBL Datasets Collection — Technical Documentation

## Table of Contents

1. [Architecture](#architecture)
2. [Dataset_collection.py](#dataset_collectionpy)
   - [BaseDataSet](#basedataset)
   - [ClassDataSet](#classdataset)
   - [ReggDataSet](#reggdataset)
   - [start()](#start-function)
3. [dataset_script.py](#dataset_scriptpy)
4. [Activity Classification System](#activity-classification-system)
5. [Data Processing Details](#data-processing-details)
6. [Output File Formats](#output-file-formats)
7. [Local Database Schema](#local-database-schema)

---

## Architecture

```
Dataset_collection.py
├── tree_act                  # Activity keyword classification dictionary
├── BaseDataSet               # Base class: data loading, preprocessing, saving
│   ├── ClassDataSet          # Adds: status labels, classification filtering
│   └── ReggDataSet           # Adds: mean aggregation, regression filtering
└── start()                   # CLI entry point

dataset_script.py
├── cuttof_target_type()      # Maps protein class → thresholds
├── pool_map_iter()           # Worker for multiprocessing pool
└── get_from_csv()            # Reads CSV, dispatches parallel jobs
```

---

## Dataset_collection.py

### BaseDataSet

Parent class implementing data loading, activity classification, unit conversion, and file saving.

#### `__init__(self, fname, act_types, sep='\t')`

| Parameter | Type | Description |
|---|---|---|
| `fname` | str | ChEMBL target ID or path to local file |
| `act_types` | list[str] | Activity types to extract (e.g. `['inhibitor', 'agonist']`) |
| `sep` | str | CSV separator for file input |

#### `load_data(db=None)`

Loads bioactivity data from ChEMBL.

- If `db` is provided: queries local SQLite database
- Otherwise: fetches from ChEMBL REST API (paginated, 1000 records per request)

**API endpoint:**
```
https://www.ebi.ac.uk/chembl/api/data/activity?format=json&limit=1000&offset={n}&target_chembl_id={id}
```

**Columns loaded:**

| Column | Source |
|---|---|
| `cmp` | `molecule_chembl_id` |
| `molecule_pref_name` | molecule name |
| `operator` | relation (`=`, `>`, `<`, `>=`, `<=`) |
| `value` | standard value |
| `pchembl_value` | ChEMBL pActivity value |
| `units` | standard units |
| `bioactivity_type` | standard type (Ki, IC50, EC50…) |
| `assay_type` | B (binding), F (functional) |
| `target_chembl_id` | target ID |
| `target_pref_name` | target name |
| `assay_description` | full assay text |
| `smiles` | canonical SMILES |
| `reference` | publication reference |

#### `add_standardized_smiles(smi_file)`

Merges standardized SMILES from an external `.smi` file (tab-separated: SMILES, ChEMBL ID).
Adds column `smiles_std`. If no file is provided, `smiles_std = smiles` (canonical SMILES).

#### `calc_MW()`

Calculates molecular weight using `RDKit.Chem.Descriptors.MolWt` from the standardized SMILES column.
Adds column `MW`.

#### `save_data(output_prefix, dataset_type)`

Saves results split by `act` value. For each activity type in `act_types`:
- Saves `{prefix}_{act}_{type}.csv` — filtered dataset
- Saves `{prefix}_{type}_log.csv` — full data log
- Saves `{prefix}_stat.csv` — statistics table

---

### ClassDataSet

Extends `BaseDataSet` for binary classification tasks.

#### `add_status(active_value, active_op, inactive_value, inactive_op)`

Labels each row as `active`, `inactive`, or `undefined`.

Active condition: `value {active_op} active_value`
Inactive condition: `value {inactive_op} inactive_value`
Otherwise: `undefined`

Valid operators: `<=`, `>=`, `!=`, `==`, `>`, `<`

Unit filter: only rows with `units` in `['nM', 'uM', 'mM', 'pM', '%']` are considered.

#### `filt_data(act_types, contr_types)`

Runs the full classification filtering pipeline:
1. `_filt_experiment_type()` — marks experiments as `norm` or `undefined`
2. `_filt_contr_status()` — identifies compounds with conflicting active/inactive labels
3. `filt_contr_data(contr_types)` — removes compounds active in opposite-type assays

#### `_filt_contr_status()`

Marks a compound as `controversial` if it has both `active` and `inactive` status within the same experiment group (grouped by: SMILES, bioactivity_type, units, act, type_act).

#### `filt_contr_data(contr_types)`

For each pair in `contr_types` (e.g., `[('agonist', 'antagonist')]`), removes any compound that appears as active in both types.

---

### ReggDataSet

Extends `BaseDataSet` for regression tasks.

#### `filt_data(act_types, contr_types)`

Runs the full regression filtering pipeline:
1. `_filt_experiment_type()` — marks experiments as `norm` or `undefined`
2. `_filt_diff_results()` — removes inconsistent measurements
3. `filt_contr_data(contr_types)` — removes compounds with opposite activity

#### `_filt_diff_results()`

Within each experiment group (SMILES, bioactivity_type, units, act, type_act):
- Marks as `del` if `max(plog) - min(plog) >= 0.5`
- Marks as `del` if any operator is not `=`

#### `mean_duplicate()`

Aggregates duplicate experiments by computing the mean `plog_value` per group.

---

### `start()` function

CLI entry point. Orchestrates the full pipeline:

```python
start(
    fname,          # ChEMBL ID or input file
    dataset_type,   # 'class' or 'reg'
    activity,       # list of bioactivity types
    act_types,      # list of activity categories
    contr_types,    # list of contradictory type pairs
    value_type,     # 'value' or 'plog_value'
    db,             # path to local SQLite DB or None
    smi_std,        # path to standardized SMILES or None
    active_value,   # threshold for active label
    active_op,      # operator for active threshold
    inactive_value, # threshold for inactive label
    inactive_op,    # operator for inactive threshold
    output,         # output path prefix
    sep             # CSV separator
)
```

**Pipeline steps:**
1. Instantiate `ClassDataSet` or `ReggDataSet`
2. `load_data(db)` — fetch from API or local DB
3. `add_standardized_smiles(smi_std)` — merge SMILES if provided
4. Unit conversion to nM, plog calculation
5. `add_status(...)` — classification only
6. `filt_data(act_types, contr_types)` — filter
7. `mean_duplicate()` — regression only
8. `calc_MW()` — if SMILES available
9. `save_data(output, dataset_type)` — write output files

---

## dataset_script.py

### `cuttof_target_type(target_class)`

Returns `(active_threshold, inactive_threshold)` based on protein class.

| `target_class` | Active (pActivity ≥) | Inactive (pActivity <) |
|---|---|---|
| `enzyme` | 7.5 | 7.5 |
| `membrane_protein` | 7.0 | 7.0 |
| `membrane_receptor` | 7.0 | 7.0 |
| `epigenetic_factor` | 6.0 | 6.0 |
| `transporter` | 6.0 | 6.0 |
| `ion_channel` | 5.0 | 5.0 |
| default | 6.0 | 6.0 |

Reference: Bosc et al., *J. Cheminf.* 2019, 11, 4.

### `pool_map_iter(args)`

Worker function for `multiprocessing.Pool`. Unpacks arguments and calls `start()`.

### `get_from_csv(csv_file, dataset_type, db, smi_std, ncpu, sep)`

Reads the input CSV, builds argument tuples for each target, and dispatches them using `multiprocessing.Pool(ncpu)`.

---

## Activity Classification System

### `tree_act` Dictionary Structure

```python
tree_act = {
    'act_class': {
        'type_act': [
            [any_keywords],      # at least one must be present
            [forbidden_keywords], # none of these may be present
            [all_keywords]        # all of these must be present
        ]
    }
}
```

### How matching works (`__check_act`)

For a given assay description and a set of keyword rules:
1. Check that at least one word from `any_keywords` appears (case-insensitive)
2. Check that no word from `forbidden_keywords` appears
3. Check that all words from `all_keywords` appear

If all conditions pass → the compound is assigned that `type_act`.

### Assignment order

`__add_act_descr()` iterates over `act_types` in the order provided by the user.
**Once a type is assigned, it is not overwritten.**
This is why more specific types must come first: `"partial agonist"` before `"agonist"`.

### Example: agonist rules

```python
'agonist': {
    'partial agonist': [
        ['partial agonist'],
        ['antagonist', 'inhibit'],
        []
    ],
    'full agonist': [
        ['full agonist'],
        ['antagonist'],
        []
    ],
    'agonist': [
        ['agonist', 'agonism', 'agonistic'],
        ['antagonist', 'partial', 'inverse'],
        []
    ]
}
```

---

## Data Processing Details

### Unit Conversion

`__convert_to_nmol()` normalizes all values to nM:

| Original unit | Conversion |
|---|---|
| `uM` | × 1000 |
| `mM` | × 1,000,000 |
| `M` | × 1,000,000,000 |
| `pM` | × 0.001 |
| `nM` | × 1 (no change) |

### plog Value Calculation

```
plog_value = -log10(value_in_nM / 1e9)
           = 9 - log10(value_in_nM)
```

If ChEMBL provides `pchembl_value` (already computed), it is used directly; otherwise `plog_value` is calculated from the raw value.

### Molecular Weight

Calculated from `smiles_std` using RDKit:
```python
from rdkit import Chem
from rdkit.Chem import Descriptors
mol = Chem.MolFromSmiles(smiles)
mw = Descriptors.MolWt(mol)
```

---

## Output File Formats

### Dataset file: `{prefix}_{act}_{type}.csv`

Tab-separated (default). Contains only rows with `exp_type == 'norm'` and the specified `act`.

Key columns:

| Column | Type | Description |
|---|---|---|
| `cmp` | str | ChEMBL compound ID |
| `smiles_std` | str | Standardized or canonical SMILES |
| `act` | str | Broad activity type |
| `type_act` | str | Specific activity type |
| `value` | float | Activity value in nM |
| `plog_value` | float | -log10(value/1e9) |
| `bioactivity_type` | str | Ki, IC50, EC50… |
| `units` | str | nM (after conversion) |
| `status` | str | active / inactive / undefined (classification only) |
| `MW` | float | Molecular weight |

### Log file: `{prefix}_{type}_log.csv`

Full data including all rows (norm, undefined, del, controversial). Used for auditing.

### Statistics file: `{prefix}_stat.csv`

Counts per activity type:

| Column | Description |
|---|---|
| `act` | Activity type |
| `active` | Number of active compounds |
| `inactive` | Number of inactive compounds |
| `undefined` | Number of undefined compounds |

---

## Local Database Schema

When using `-d` (local SQLite ChEMBL database), the following tables are queried:

| Table | Key columns used |
|---|---|
| `target_dictionary` | `chembl_id`, `tid` |
| `assays` | `tid`, `assay_id`, `assay_type`, `description` |
| `activities` | `assay_id`, `molregno`, `standard_value`, `standard_units`, `standard_type`, `standard_relation`, `pchembl_value`, `doc_id` |
| `molecule_dictionary` | `molregno`, `chembl_id`, `pref_name` |
| `compound_structures` | `molregno`, `canonical_smiles` |
| `docs` | `doc_id`, `authors`, `journal`, `year`, `volume`, `first_page` |

The SQL query performs a multi-table JOIN across all these tables filtered by `target_chembl_id`.
