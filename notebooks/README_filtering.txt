NOTE: We decided to create 2 function files to seperate the filtering functions from the functions creating the new, 'better' databases
###### Utilisation de la fonction filter_Isopropyl

La fonction `filter_Isopropyl` permet de filtrer les molécules d'une base de données pandas DataFrame en fonction de la présence d'un mot-clé spécifique dans la colonne 'Molecule'. 

### Exemple d'utilisation

Supposons que vous avez un DataFrame nommé `db` qui contient une colonne 'Molecule'. Pour filtrer toutes les entrées contenant le mot-clé 'Isopropyl', utilisez :

```python
from nom_du_package.filtering import filter_molecules

filtered_db = filter_molecules(db)
print(filtered_db)

### La fonction utilisée sur notre projet: 

On peut voir que la creme 'La Rosée' contient de l'Isopropyl, on va donc implementer une nouvelle fonction qui retire selectivement ces crèmes.









###### Utilisation de la fonction `new_database_no_Isopropyl`

La fonction `new_database_no_Isopropyl` permet d'exclure des lignes d'une base de données pandas DataFrame en fonction de la présence d'une valeur spécifique dans une colonne donnée. Cette fonction est utile pour éliminer les données non désirées avant d'effectuer des analyses plus détaillées.

### Exemple d'utilisation

Supposons que vous avez un DataFrame nommé `db` qui contient une colonne 'Cream'. Pour exclure toutes les entrées où la colonne 'Cream' contient la valeur 'a', utilisez :

```python
import pandas as pd
from nom_du_package.filtering import new_database_no_Isopropyl

### Création d'un DataFrame exemple
db = pd.DataFrame({
    'Cream': ['a', 'd', 'b', 'g'],
    'Molecule': ['Isopropyl alcohol', 'Ethanol', 'Isopropyl acetate', 'Water']
})

### Application de la fonction pour exclure 'a'
filtered_db = new_database_no_Isopropyl(db, column_name='Cream', excluded_value='a')
print(filtered_db)

### La fonction utilisée sur notre projet: 
On obtient une nouvelle data base sans la creme "La Rosee" qui s'est revelee mauvaise

### Utilisation externe:
ici, notre fonction filtre selectivement La Rosee. Pour filtrer selectivement les cremes contenant de l'Isopropyl d'une autre molecule, il faut simplement remplacer "La Rosée" dans la fonction par le nom de sa creme. Si plusieurs cremes doivent etre filtrer, modifier au niveau de 'excluded_value="La Rosée" ' et ecrire :excluded_values=["a", "b",(etc)]).









###### Utilisation de la fonction `filter_NaOH`

La fonction `filter_NaOH` est conçue pour filtrer les entrées d'un DataFrame pandas en fonction de la présence d'un mot-clé spécifique dans une colonne désignée. Elle est particulièrement utile pour isoler les données chimiques spécifiques telles que la présence de "Sodium Hydroxide" dans la colonne 'Molecule'.

### Exemple d'utilisation

Supposons que vous disposez d'un DataFrame nommé `db` qui contient une colonne 'Molecule'. Pour filtrer toutes les entrées contenant le mot-clé 'Sodium Hydroxide', utilisez la fonction `filter_NaOH` comme suit :

```python
import pandas as pd
from nom_du_package.filtering import filter_NaOH

# Exemple de DataFrame
db = pd.DataFrame({
    'Molecule': ['Sodium Hydroxide', 'Water', 'Sodium Chloride', 'Sodium Hydroxide solution']
})

# Application de la fonction pour filtrer par 'Sodium Hydroxide'
filtered_db = filter_NaOH(db, column_name='Molecule', keyword='Sodium Hydroxide')
print(filtered_db)

#Resultats attendus:
                     Molecule
0           Sodium Hydroxide
3  Sodium Hydroxide solution

### La fonction utilisée sur notre projet: 
On obtient une liste de toutes les cremes qui contiennent du NaOH.

### Utilisation externe:
ici, notre fonction est apliquee sur notre base de donnée filtrée (sans Isopropyl). Pour filtrer le NaOH de n'importe quelle base de donnée, simplement appliquée dans l'execution du programme sur la base de donnée utilisée.









###### Utilisation de la fonction `new_database_no_NaOH`

La fonction `new_database_no_Isopropyl` permet d'exclure des lignes d'une base de données pandas DataFrame en fonction de la présence d'une valeur spécifique dans une colonne donnée. Cette fonction est utile pour éliminer les cremes contenant du NaOH.
Cette fonction a la meme utilite que la fonction `new_database_no_Isopropyl`

### Exemple d'utilisation

Supposons que vous avez un DataFrame nommé `db` qui contient une colonne 'Cream'. Pour exclure toutes les entrées où la colonne 'Cream' contient la valeur 'd', utilisez :

```python
import pandas as pd
from nom_du_package.filtering import new_database_no_Isopropyl

### Création d'un DataFrame exemple
db = pd.DataFrame({
    'Cream': ['c', 'd', 'e', 'f'],
    'Molecule': ['Isopropyl alcohol', 'Sodium Hydroxyde', 'Isopropyl acetate', 'Water']
})

### Application de la fonction pour exclure 'd'
database_filtered_2 = new_database_no_Isopropyl(db, column_name='Cream', excluded_value='d')
print(database_filtered_2)

### La fonction utilisée sur notre projet: 
On obtient une nouvelle data base sans les cremes 'La Roche-Posay','Nivea','Avene' et 'La Prairie' qui se sont revelees mauvaises

### Utilisation externe:
ici, comme precedemment, notre fonction filtre notre creme d'interet. Pour filtrer selectivement les cremes contenant du NaOH, il faut simplement remplacer "Nivea" dans la fonction par le nom de sa creme. Si plusieurs cremes doivent etre filtrer, modifier au niveau de 'excluded_value="Nivea" ' et ecrire :excluded_values=["d", "e", (etc)]).










###### Utilisation de la fonction `filter_allother`
La fonction `filter_allother` est conçue pour filtrer les entrées d'un DataFrame pandas en fonction de la présence de plusieurs mots-clés spécifiques dans une colonne désignée. Elle est particulièrement utile pour isoler les données chimiques spécifiques telles que la présence de Linalool|Citronellol|Limonene|Benzyl Alcohol|Benzyl Salicylate qui sont des agents pouvant eventuellement etre irritant dans la colonne 'Molecule'.
À noter que ces actifs ne sont pas vraiment dangeureux pour la peau, mais il est mieux d'avoir une creme sans.

### Exemple d'utilisation

Supposons que vous disposez d'un DataFrame nommé `db` qui contient une colonne 'Molecule'. Pour filtrer toutes les entrées contenant les mots-clés 'Linalool|Citronellol|Limonene|Benzyl Alcohol|Benzyl Salicylate', utilisez la fonction `filter_allother` comme suit :

```python
import pandas as pd
from nom_du_package.filtering import filter_NaOH

# Exemple de DataFrame
db = pd.DataFrame({
    'Molecule': ['Sodium Hydroxide', 'Linalool', 'Sodium Chloride', 'Limonene','Water','Citronellol']
})

# Application de la fonction pour filtrer par 'Sodium Hydroxide'
filtered_db = filter_NaOH(db, column_name='Molecule', keyword='Sodium Hydroxide')
print(filtered_db)

#Resultats attendus:
                     Molecule
1           Linalool
3           Limonene
5           Citronellol

### La fonction utilisée sur notre projet: 
On obtient une liste de toutes les cremes qui contiennent ces actifs (à noter que certaines crèmes en contiennent plusieurs).

### Utilisation externe:
ici, notre fonction est apliquee sur notre base de donnée filtrée (sans Isopropyl ni NaOH). Pour filtrer ces actifs de n'importe quelle base de donnée, simplement appliquée dans l'execution du programme sur la base de donnée utilisée.












###### Utilisation de la fonction `new_database_no_other`

La fonction `new_database_no_other` permet d'exclure des lignes d'une base de données pandas DataFrame en fonction de la présence d'une valeur spécifique dans une colonne donnée. Cette fonction est utile pour éliminer les cremes contenant plusieursactifs specifiques mentionnés au dessus.
Cette fonction a la meme utilite que les fonctions `new_database_no_Isopropyl` et `new_database_no_NaOH`


# Example usage
data = {
    'Cream': ['La Rosée', 'Nivea', 'La Roche-Posay', 'Dove', 'Avene', 'La Prairie', 'Other']
}
df = pd.DataFrame(data)

# Using the function to exclude multiple brands
excluded_brands = ["La Rosée", "La Roche-Posay", "Nivea", "Avene", "La Prairie"]
filtered_df = filter_exclusions_simple(df, 'Cream', excluded_brands)
print(filtered_df)


#Resultats attendus:
                    Cream
0           La Rosée
1           Nivea
2           La Roche Posay
4           Avene
5           La Prairie



### La fonction utilisée sur notre projet: 
On obtient une nouvelle data base sans les cremes 'Nuxe', 'Diadermine', 'L'Oréal', 'La Neige'
### Utilisation externe:
ici, comme precedemment, notre fonction filtre nos cremes d'interet. Pour filtrer selectivement les cremes contenant ces actifs, il faut simplement remplacer les cremes mentionnées dans la fonction par le nom de sa creme. Si plusieurs cremes doivent etre filtrer, modifier au niveau de 'excluded_value="[Nuxe","Diadermine"...] ' et ecrire :excluded_values=["d", "e", (etc)]).



































