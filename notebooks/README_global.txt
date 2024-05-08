To use this package with it's functions, first, write the following comands:

      import pandas as pd
      url="path/to/file"
###for us, the path is the following: "https://raw.githubusercontent.com/victoirelang/Project/main/data/database.csv"
      database= pd.read_csv(url, sep=';')
###lets show the database to ensure everything is correct
      database 

      