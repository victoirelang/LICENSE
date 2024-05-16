This document is to test the functions created on our database.

It helps us insure they were programmed properly.

It is used after each new function implemented.

# Use:

First, use the following command:

    import pandas as pd

Make sure you have imported your database properly, you can do so using the following command right at the start:

   url="https://raw.githubusercontent.com/path/to/your/database"
   database= pd.read_csv(url, sep=';') ###reads your database
   database ###shows you your database


Since the scrpit and functions aren't in the same folder in this package, make sure to use the following line:

import sys
sys.path.append('../')  # Ajoute le dossier parent à sys.path

Also, make sure before using each function, you import them the following way:

   from functions.nameoffunctionfile import name_of_function