import pandas as pd

def new_database_no_Isopropyl(database, column_name="Cream", excluded_value="La Rosée"):
    database_filtered_1= database[database[column_name] != excluded_value]
    return database_filtered_1


def new_database_no_NaOH(database_filtered_1, column_name, excluded_values=["La Roche-Posay","Nivea","Avene","La Prairie"]:
    database_filtered_2= database_filtered_1[database_filtered_1[column_name].isin(excluded_values)]
    return database_filtered_2
    