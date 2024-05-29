# Project of Practical Programming in Chemistry 🤩🤩🤩

## Description
This Jupyter notebook demonstrates the use of custom functions to filter and remove items from a database. It uses Pandas to load data from a CSV file and applies several filtering functions to prepare the data for further analysis.

## Prerequisites ‼️
To run this notebook, make sure you have installed Python and the following packages:
- Pandas
- rdkit
- Chem (imported from rdkit)
- PandasTools (imported from rdkit)
- Draw (imported from rdkit.Chem)
- display (imported from IPython.display)


You can install the necessary dependencies by running the following command:
```bash
pip install pandas

conda create -c conda-forge -n my-rdkit-env rdkit
(do not forget to activate environnement using :conda activate my-rdkit-env)
```

## Installation 🤭
```bash
git clone https://github.com/victoirelang/Project.git
cd Project
pip install -e .
```


## [Data Folder:](https://github.com/victoirelang/Project/tree/main/data)
To carry out our comparative analysis of hydrating creams, it was necessary to have a database compiling detailed information on the compounds present in each of the 10 selected creams. Since no such database existed, we took the initiative to create our own. For each cream, we researched and documented the active compounds, filling in 22 columns of information for each one. These columns include data on the chemical formula, physical and chemical properties, functionality, origin of the compounds, as well as details on safety, toxicity, and suitability for different age groups. 
This customized database allowed us to then perform precise comparisons between products. It was developed with the help of ChatGPT and various websites, ensuring the completeness and accuracy of the information collected.

## [Notebooks Folder:](https://github.com/victoirelang/Project/tree/main/notebooks)
More details concerning each function can be found in a seperate README file in the "function module".

## [src Folder:](https://github.com/victoirelang/Project/tree/main/src)

This folder contains scripts and functions.
## [Functions Folder:](https://github.com/victoirelang/Project/tree/main/functions)
Here, we'll explain what the different functions of our project do.

#### filtering.ipnyb

These functions filter entries in a pandas DataFrame database according to the presence of a specific keyword in a designated column. 
They are particularly useful for isolating specific chemical datas, such as the presence of "Isopropyl" or “Sodium Hydroxide” in the 'Molecule' column.

#### newdatabase.ipnyb

These functions allow you to exclude rows from a pandas DataFrame database based on the presence of a specific value in a given column. In our project, they are useful for eliminating creams containing specific bad molecules, and so help to rank and select the best creams.

#### ranking_for_pregnancy.ipnyb

By creating dictionaries, these functions have made it possible to rank the various creams in ascending order of their harmfulness to pregnant women. This ranking is based on the number of molecules hazardous to pregnancy.

#### visualizing.ipnyb

Using RDKit, this function enables the user to visualize the molecules present in a given cream, and colors molecules deemed harmful. Some molecules contain harmful groups, and are only partially colored, since only these groups are colored. However, in this case, the molecules are not really harmful, since the harmful groups are contained in the molecule and not present on their own: their harmful effect is therefore cancelled out, or at least greatly reduced.
Finally, this function posed a problem when it was first created, as Smiles and smarts structures were initially mixed. We therefore ended up using only Smiles structures.

## [Tests Folder:](https://github.com/victoirelang/Project/tree/main/tests)

For each function, a test was executed. All tests were successful.

## Limits of our project 🫣

It is important to note that our project has certain limitations. Firstly, our “homemade” database may contain errors, as it was manually assembled and verified. Additionally, we analyzed only ten creams which means our sampling is limited and may not reflect the full variety of products available on the market. This is due to time constraints since filling in the data base for every compound was very long. 
Another point to consider is that our sorting method only takes into account the presence of certain molecules, without evaluating their concentration, which can play a crucial role in the harmfulness of the products. Furthermore, we opted to sort based on the negative aspects of compounds. However, there are many other ways to analyze the data, for example we could have highlighted the benefits of some molecules for the skin by filtrating differently. 
These limitations should be taken into account when interpreting the results of our analysis.
