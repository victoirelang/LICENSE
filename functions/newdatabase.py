import pandas as pd

def new_database_no_Isopropyl(database, column_name="Cream", excluded_value="La Rosée"):
    database_filtered_1= database[database[column_name] != excluded_value]
    return database_filtered_1


def new_database_no_NaOH(database, column_name="Cream", excluded_values=["La Roche-Posay","Nivea","Avene","La Prairie"]):
    database_filtered_2= database[~database[column_name].isin(excluded_values)]
    return database_filtered_2

def new_database_no_other(database, column_name="Cream", excluded_values=["Nuxe","Diadermine","L'Oréal","La Neige"]):
    return pd.DataFrame([row for idx, row in database.iterrows() if row[column_name] not in excluded_values])
    