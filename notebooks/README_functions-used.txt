NOTE: We decided to create 2 function files to seperate the filtering functions from the functions creating the new said 'better' databases

# Use of the function `filter_Isopropyl`

The `filter_Isopropyl` function filters molecules in a pandas DataFrame database according to the presence of a specific keyword in the 'Molecule' column. 

## Example of use

Suppose you have a DataFrame named `db` which contains a 'Molecule' column. To filter out all entries containing the keyword 'Isopropyl', use :

```python
from nom_du_package.filtering import filter_molecules

filtered_db = filter_molecules(db)
print(filtered_db)

## The function used on our project

We can see that the 'La Rosée' cream contains Isopropyl, so we're going to implement a new function that selectively removes these creams.

# Use of the function `new_database_no_Isopropyl`

The `new_database_no_Isopropyl` function allows you to exclude rows from a pandas DataFrame database based on the presence of a specific value in a given column. This function is useful for eliminating unwanted data before carrying out more detailed analyses.

## Example of use

Suppose you have a DataFrame named `db` which contains a 'Cream' column. To exclude all entries where the 'Cream' column contains the value 'a', use :

```python
import pandas as pd
from nom_du_package.filtering import new_database_no_Isopropyl

## Example of creating a DataFrame

db = pd.DataFrame({
    'Cream': ['a', 'd', 'b', 'g'],
    'Molecule': ['Isopropyl alcohol', 'Ethanol', 'Isopropyl acetate', 'Water']
})

## Apply function to exclude 'a'

filtered_db = new_database_no_Isopropyl(db, column_name='Cream', excluded_value='a')
print(filtered_db)

## The function used on our project

A new data base without the “La Rosee” cream, which turned out to be bad, is obtained.

## External use

Here, our function selectively filters La Rosee. To selectively filter creams containing Isopropyl, simply replace “La Rosée” in the function with the name of the cream. If several creams are to be filtered, modify 'excluded_value=“La Rosée” ' and write :excluded_values=[“a”, “b”,(etc)]).

# Use of the function `filter_Isopropyl`

The `filter_NaOH` function is designed to filter entries in a pandas DataFrame according to the presence of a specific keyword in a designated column. It is particularly useful for isolating specific chemical data, such as the presence of “Sodium Hydroxide” in the 'Molecule' column.

## Example of use

Suppose you have a DataFrame named `db` which contains a 'Molecule' column. To filter out all entries containing the keyword 'Sodium Hydroxide', use the `filter_NaOH` function as follows:

```python
import pandas as pd
from nom_du_package.filtering import filter_NaOH

## Example of DataFrame

db = pd.DataFrame({
    'Molecule': ['Sodium Hydroxide', 'Water', 'Sodium Chloride', 'Sodium Hydroxide solution']
})

### Application of the 'Sodium Hydroxide' filter function

filtered_db = filter_NaOH(db, column_name='Molecule', keyword='Sodium Hydroxide')
print(filtered_db)

### Expected results

Molecule
0           Sodium Hydroxide
3  Sodium Hydroxide solution

## The function used on our project

The result is a list of all creams containing NaOH.

## External use

Here, our function is applied to our filtered database (without Isopropyl). To filter NaOH from any database, simply apply it to the database you're using.

# Use of the function `new_database_no_NaOH`

The `new_database_no_NaOH` function allows you to exclude rows from a pandas DataFrame database based on the presence of a specific value in a given column. This function is useful for eliminating creams containing NaOH.
This function has the same utility as the `new_database_no_Isopropyl` function.

## Example of use

Suppose you have a DataFrame named `db` which contains a 'Cream' column. To exclude all entries where the 'Cream' column contains the value 'd', use :

```python
import pandas as pd
from nom_du_package.filtering import new_database_no_Isopropyl

## Example of creating a DataFrame

db = pd.DataFrame({
    'Cream': ['c', 'd', 'e', 'f'],
    'Molecule': ['Isopropyl alcohol', 'Sodium Hydroxyde', 'Isopropyl acetate', 'Water']
})

## Apply function to exclude 'd'

database_filtered_2 = new_database_no_Isopropyl(db, column_name='Cream', excluded_value='d')
print(database_filtered_2)

## The function used on our project

A new data base is obtained, without the creams 'La Roche-Posay', 'Nivea', 'Avene' and 'La Prairie' which proved to be bad.

## External use

Here, as before, our function filters our cream of interest. To selectively filter NaOH-containing creams, simply replace “Nivea” in the function with the name of your cream. If several creams are to be filtered, modify 'excluded_value=“Nivea” ' and write :excluded_values=[“d”, “e”, (etc)]).

# Use of the function `filter_allother`

The `filter_allother` function is designed to filter pandas DataFrame entries based on the presence of several specific keywords in a designated column. It is particularly useful for isolating specific chemical data such as the presence of Linalool|Citronellol|Limonene|Benzyl Alcohol|Benzyl Salicylate, which are potentially irritating agents in the 'Molecule' column.
Note that these active ingredients are not really harmful to the skin, but it's better to have a cream without them.

## Example of use

Suppose you have a DataFrame named `db` which contains a 'Molecule' column. To filter out all entries containing the keywords 'Linalool|Citronellol|Limonene|Benzyl Alcohol|Benzyl Salicylate', use the `filter_allother` function as follows:

```python
import pandas as pd
from nom_du_package.filtering import filter_NaOH

### Example of a DataFrame

db = pd.DataFrame({
    'Molecule': ['Sodium Hydroxide', 'Linalool', 'Sodium Chloride', 'Limonene','Water','Citronellol']
})

### Application of the filter function

filtered_db = filter_allother(db, column_name='Molecule', keyword='Linalool|Citronellol|Limonene|Benzyl Alcohol|Benzyl Salicylate')
print(filtered_db)

### Expected results

             Molecule
1           Linalool
3           Limonene
5           Citronellol

## The function used on our project

The result is a list of all creams containing these active ingredients (note that some creams contain more than one).

## External use

Here, our function is applied to our filtered database (without Isopropyl or NaOH). To filter these assets from any database, simply apply the program to the database in question.

# Use of the function `new_database_no_other`

The `new_database_no_other` function allows you to exclude rows from a pandas DataFrame database based on the presence of a specific value in a given column. This function is useful for eliminating creams containing several of the specific active ingredients mentioned above.
This function has the same utility as the `new_database_no_Isopropyl` and `new_database_no_NaOH` functions.

### Example of use

data = {
    'Cream': ['La Rosée', 'Nivea', 'La Roche-Posay', 'Dove', 'Avene', 'La Prairie', 'Other']
}
df = pd.DataFrame(data)

### Using the function to exclude multiple brands

excluded_brands = ["La Rosée", "La Roche-Posay", "Nivea", "Avene", "La Prairie"]
filtered_df = filter_exclusions_simple(df, 'Cream', excluded_brands)
print(filtered_df)

### Expected results

             Cream
0           La Rosée
1           Nivea
2           La Roche Posay
4           Avene
5           La Prairie

## The function used on our project

We obtain a new data base without cremes 'Nuxe', 'Diadermine', 'L'Oréal', 'La Neige'.

## External use

Here, as before, our function filters our creams of interest. To selectively filter creams containing these active ingredients, simply replace the creams mentioned in the function with the name of your cream. If several creams are to be filtered, modify 'excluded_value=“[Nuxe”, “Diadermine”...] ' and write :excluded_values=[“d”, “e”, (etc)]).

# Use of the function `filter_pregnancy_risk'

This function filters molecules in a pandas DataFrame database according to the presence of a risk for pregnant women in the 'Pregnancy Risk' column.

## Example of use

Let's suppose you have a DataFrame named "data" which contains a 'Pregnancy Risk' column. To filter out all entries containing a risk for pregnant women, use :

from nom_du_package.filtering import filter_pregnancy_risk

filtered_data = filter_pregnancy_risk(data)
print(filtered_data)

## The function used on our project

We can see that some molecules present a risk for pregnancy. Thus, we're going to implement a new function that ranks creams according to how dangerous they are for pregnant women. 

# Use of the function 'dangerous_pregnancy'

This function filters a database to identify molecules hazardous to pregnant women, and calculates the percentage of such molecules for each cream.

## Example of use

Suppose you have a DataFrame named database which contains the columns 'Cream', 'Molecule', and 'Pregnancy Risk'. To calculate the percentage of dangerous molecules for each cream, use :

from nom_du_package.filtering import dangerous_pregnancy

pourcentages = dangerous_pregnancy(database)
print(pourcentages)

## The function used on our project

This function allows us to calculate and display the percentage of molecules hazardous to pregnancy per cream, and therefore to know which ones to prioritize for pregnant women.

# Use of the function 'visualize_molecules_for_cream'

This function displays the molecules of a certain category of creams with specific atoms highlighted in a pandas DataFrame.

## Prerequesites

Before running this function, you must use first:

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG, display

## Example of use

Suppose you have a DataFrame named df containing 'Cream' and 'Smiles' columns. To view the molecules for a specific cream, use :

from nom_du_package.visualization import visualize_molecules_for_cream

visualize_molecules_for_cream(df, 'Nom_de_la_Crème')

## The function used on our project

The function allows us to visualize the groups and molecules deemed harmful in each cream. It should be noted that when only part of a molecule is colored, the group deemed harmful is then anhilated, and the molecule is in fact no longer irritating, or at least much less so for humans.
Moreover, this function posed a problem when it was first created, as Smiles and smarts structures were initially mixed. We therefore ended up using only Smiles structures.