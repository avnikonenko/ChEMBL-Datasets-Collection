# Example files

## Standardized SMILES

| File | Description |
|---|---|
| `chembl_33_step4.smi.tar.gz` | Full standardized SMILES from ChEMBL 33. Tab-separated, columns: SMILES, ChEMBL compound ID. Pass with `-s` / `--smi_std`; `.tar.gz` is read directly without manual extraction. |

## Target list for `dataset_script.py`

`dataset_script.py` takes a tab-separated CSV with a required `chembl_id` column and optional columns that control thresholds and activity filtering.

| File | Columns used | Description |
|---|---|---|
| `chembl29_test.csv` | `chembl_id`, `target_class` | Full ChEMBL export; `target_class` drives class-specific thresholds (Bosc et al., J. Cheminf. 2019). |
| `targets_custom_cutoffs.csv` | `chembl_id`, `cutoff_act`, `cutoff_inact`, `activity`, `experiment_type` | Per-target numeric cutoffs; overrides any `target_class` thresholds. `activity` and `experiment_type` are `;`-separated lists of bioactivity types (e.g. `Ki;IC50`). |
| `targets_category.csv` | `chembl_id`, `target_class`, `category` | `category` restricts which activity mechanism types are retained (`;`-separated, e.g. `inhibitor;allosteric`). Allowed values: `allosteric`, `antagonist`, `agonist`, `activator`, `blockator`, `inhibitor`, `binding`, `other`. |

### Supported `target_class` values

| Value | Active threshold (plog) | Notes |
|---|---|---|
| `enzyme` | ≥ 7.5 | |
| `epigenetic_factor` / `epigenetic regulator` | ≥ 6 | Both spellings accepted |
| `transporter` | ≥ 6 | |
| `ion_channel` | ≥ 5 | |
| `membrane_protein` / `membrane_receptor` | ≥ 7 | |
| *(any other)* | ≥ 6 | Default thresholds |

## Other

| File | Description |
|---|---|
| `chembl_ids` | Plain list of ChEMBL target IDs, one per line. Not used by the scripts directly; provided as a reference ID set. |
| `run_examples.sh` | Ready-to-run shell commands demonstrating common invocations of both `Dataset_collection.py` and `dataset_script.py`. |
