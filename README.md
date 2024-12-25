# ChEMBL datasets preparation
**Dataset collection and preprocessing script** 
a Python-based tool to automate the retrieval and preprocessing of datasets from the ChEMBL database, tailored to specific activity types (e.g., agonists, antagonists, inverse agonist
## Pipeline
**1. Collect activity data for preset target by ChEMBLID from ChEMBL(from api version or local db)**  
Received data:
- ChEMBLID of compound
- molecule_pref_name
- operator: <,>,=,>=,<=
- activity value
- pchembl_value: log of the activity value 
- units: nM, %, mM etc 
- bioactivity_type: Ki, Kd, IC50, EC50 etc.
- assay_type: B, F etc.
- target_chembl_id
- target_pref_name
- assay_description
- smiles
- reference  

**2. Assign type of activity**  
_Columns "act" and "type_act"._  
_"act"_: Search combination of key words in the assay_description column to assign act of compound among chosen by user.  
Also we check some combination of key words which are not allowed to present in the assay description.

_For example._ Type of act is "activator".  
 Search all words specified in the list in the assay description:   
 ['activation', 'activity at', 'activity towards', 'activity against', 'activity of', 'induction',
                 'activation'].  
Then check if there are not any "anti key" words:  ['inhibition', 'inhibitory','inhibitors', 'inhibitor', 'inhibit', 'antagonist', 'inactivation']. 
```
'activator': {'activator': [
                ['activation', 'activity at', 'activity towards', 'activity against', 'activity of', 'induction',
                 'activation'],
                ['inhibition', 'inhibitory','inhibitors', 'inhibitor', 'inhibit', 'antagonist', 'inactivation'], []]}, 

```
Allowed types:
- allosteric
- antagonist
- agonist
- inhibitor
- activator
- blockator
- binding
- other (use this option if you want just collect data without any specification of act)

_"type act"_
More detail specifictaion of act.
Allowed types:
- allosteric
- inverse agonist
- competitive antagonist
- non-competitive antagonist
- antagonist
- partial agonist
- full agonist
- agonist
- activator
- blockator
- competitive inhibitor
- uncompetitive inhibitor
- non-competitive inhibitor
- inhibitor
- binding
- other

You can find more details using --show_tree_act argument  
**! Warning:** The order of parameters in act/act_type is important. Set it in order from most specific type_act to least. Because once a type is assigned, it cannot be overwritten, and the order of assignment is equal to the order of the given arguments. Example: -a "partial agonist" "agonist" but not -a "agonist" "partial agonist".    
```
python Dataset_collection.py --show_tree_act
```
Tree_act dictionary structure: {act_class: and its key words [[any],[not],[all]]}

**3. Preprocessing**
- Remove compounds with undefined value of activity
- Add standardized smiles column if special file was provided (see example)
- Calculate MW  

_Classification dataset_
 - Add column _Status_  
 Set it as "Active" if:
  1) experiment has bioactivity_type value listed by user in activity argument. _Exp: Ki Kd IC50_
  2) Value of activity in range of given threshold (value <= given threshold)
  3) Units of activity in nM, uM, mM (internal conversation to nM. Original data is saved in Log_Note column of log files)
  
_Regression dataset_
- Get mean value for compound between similar experiments  
(Group by (standardized or canonical) smiles, bioactivity type, units, act and type of act)

**4. Remove controversial activity**
- Remove compound if it is active in experiments of the opposite type  
_Example_: remove compound which is active as agonist and antagonist in different experiments towards the target

_Classification dataset_
- Remove compound if it has more then 1 status of activity (active and inactive) in the similar experiments

_Regression dataset_
 - Remove compound if difference between min and max log of value equals or greater than 0.5 in the similar experiments.
 - Remove compound if at least one of the group of similar experiments operator not equals '=' 

**5. Save processed data**
 - Remove duplicates of compounds in the similar experiments
 - Split and save compounds by given by user act values (Exp: inhibitor, agonist etc.)
 - Save log file containing all original data
 - Calculate and save statistics data (Number of active, inactive, undefined compounds in each act group)

## Dependency
* **Python (>=3.6)**
* **RDKit**
* **Numpy**
* **Pandas**
* **Requests**

## Installation
```
conda create -n dataset-collect
conda activate dataset-collect
```
**Install dependencies from conda**
```
conda install -c conda-forge rdkit pandas requests
```

## Usage
```
source activate dataset-collect
```
**1. Set all values via command line**
```
cat example/chembl_ids | xargs -I {1} python Dataset_collection.py -i {1} --t reg --activity EC50 IC50 -v plog_value -a inhibitor
```
**2. Set all values via recommended thresholds**   
Using thresholds recommended for different families of protein targets in the paper of Bosc et al. _[N. Bosc, F. Atkinson, E. Felix, A. Gaulton, A. Hersey, A. R. Leach,
J. Cheminf. 2019, 11, 4.]_.   
Compounds are labeled active if their pKi or pIC50
is >= 6 for enzyme targets and >= 7.5 for membrane proteins, and inactive otherwise.  
**Using protein classes:** _epigenetic_factor, transporter, ion_channel, enzyme, membrane_protein, membrane_receptor_
```
python dataset_script.py -i  example/chembl29_test.csv -n 4 
```
**3. Collect data from ChEMBL without specification of act**
Dataset_collection.py -i your_ChEMBLID --t reg --activity EC50 IC50 -v plog_value -a other


## Command-line
```
$ python Dataset_collection.py -h
usage: Dataset_collection.py [-h] -i ChEMBLID --t class/reg --activity
                             [Ki IC50 EC50 [Ki IC50 EC50 ...]] [-o output]
                             [-a [binding inhibitor activator [binding inhibitor activator ...]]]
                             [-c [inhibitor activator [inhibitor activator ...]]]
                             [-v value/plog_value] [-d ChEMBL29.db]
                             [-s smi_std.smi] [--active_value number]
                             [--active_op operator] [--inactive_value number]
                             [--inactive_op operator] [--sep separator]
                             [--show_tree_act]

optional arguments:
  -h, --help            show this help message and exit
  -i ChEMBLID, --input ChEMBLID
  --t class/reg
  --activity [Ki IC50 EC50 [Ki IC50 EC50 ...]]
  -o output, --output output
  -a [binding inhibitor activator [binding inhibitor activator ...]], --act [binding inhibitor activator [binding inhibitor activator ...]]
  -c [inhibitor activator [inhibitor activator ...]], --contr [inhibitor activator [inhibitor activator ...]]
  -v value/plog_value, --value value/plog_value
  -d ChEMBL29.db, --db ChEMBL29.db
                        If value is not set api ]version of ChEMBL will be used (default: None)
  -s smi_std.smi, --smi_std smi_std.smi
                        standardized smiles (default: None)
  --active_value number
  --active_op operator
  --inactive_value number
  --inactive_op operator
  --sep separator
  --show_tree_act
```

## License
BSD-3
