import pandas as pd

def filter_Isopropyl(database, column_name="Molecule", keyword="Isopropyl"):
    database_contains_Isopropyl= database[database[column_name].str.contains(keyword, na=False)]
    return database_contains_Isopropyl

def filter_NaOH(database_filtered_1, column_name="Molecule", keyword="Sodium Hydroxide"):
    database_contains_NaOH= database_filtered_1[database_filtered_1[column_name].str.contains(keyword, na=False)]
    return database_contains_NaOH

def filter_allother(database, column_name="Molecule", keywords="Linalool|Citronellol|Limonene|Benzyl Alcohol|Benzyl Salicylate"):
    database_filtered_2=database[database[column_name].str.contains(keywords, na=False)]
    return database_filtered_2

