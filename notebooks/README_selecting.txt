This document is to test the functions created on our database.

It helps us insure they were programmed properly.

It is used after each new function implemented.

######Prerequisite:
Pour exécuter ce notebook, assurez-vous d'avoir installé Python et les packages suivants :
- Pandas

Vous pouvez installer les dépendances nécessaires en exécutant la commande suivante :
```bash
pip install pandas


#####Usage:
first, use the following command:
    import pandas as pd

make sure you have imported your database properly, you can do so using the following command right at the start:
   url="https://raw.githubusercontent.com/path/to/your/database"
   database= pd.read_csv(url, sep=';') ###reads your database
   database ###shows you your database


since the scrpit and functions aren't in the same folder in this package, make sure to use the following line:
import sys
sys.path.append('../')  # Ajoute le dossier parent à sys.path

also, make sure before using each function, you import them the following way:
   from functions.nameoffunctionfile import name_of_function

