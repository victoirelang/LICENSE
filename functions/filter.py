 import pandas as pd

def filter_Isopropyl(database, column_name="Molecule", keyword="Isopropyl"):
    database_contains_Isopropyl= database[database[column_name].str.contains(keyword, na=False)]
    return database_contains_Isopropyl



