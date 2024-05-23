# Project of Practical Programming in Chemistry

## Description
This Jupyter notebook demonstrates the use of custom functions to filter and remove items from a database. It uses Pandas to load data from a CSV file and applies several filtering functions to prepare the data for further analysis.

## Prerequisites
To run this notebook, make sure you have installed Python and the following packages:
- Pandas
- rdkit
- Chem (imported from rdkit)
- PandasTools (imported from rdkit)


You can install the necessary dependencies by running the following command:
```bash
pip install pandas

conda create -c conda-forge -n my-rdkit-env rdkit
(do not forget to activate environnement using :conda activate my-rdkit-env)
```

## Installation:
```bash
git clone https://github.com/victoirelang/Project.git
cd Project
pip install -e .
```


## Data
explain here that we have our own database we created from sratch

## Notebooks
here explain we have other README files which go more into details for example for our functions 

## src
contains scripts + functions 
### [Functions Folder:](https://github.com/victoirelang/Project/tree/main/functions)
Here, we'll explain what the different functions of our project do.

#### filtering.ipnyb

These functions filter entries in a pandas DataFrame database according to the presence of a specific keyword in a designated column. 
They are particularly useful for isolating specific chemical datas, such as the presence of "Isopropyl" or “Sodium Hydroxide” in the 'Molecule' column.

#### newdatabase.ipnyb

These functions allow you to exclude rows from a pandas DataFrame database based on the presence of a specific value in a given column. In our project, they are useful for eliminating creams containing specific bad molecules, and so help to rank and select the best creams.

#### ranking_for_pregnancy.ipnyb

By creating dictionaries, these functions have made it possible to rank the various creams in ascending order of their harmfulness to pregnant women. This ranking is based on the number of molecules hazardous to pregnancy.


More details concerning each function can be found in a seperate README file in the "function module".

#### visualizing.ipnyb

expliquer ce que fait la fonction, qu';on a fait pour avoir une facon visuelle de filtrer, qu'on a utiliser RDKIT (va voir tt les trucs que j'importe), qu'on a eu des pb parce qu'on melangaient les structures smarts et smiles et qu'on a tout fait avec smiles. Dire donc ce que ca rend etc regarde la fonction et explique puis les conclusions qu'on tire


## Tests

