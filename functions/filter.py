import pandas as pd

def a(database, column_name="Molecule", keyword="Isopropyl"):
    """
    Filtre le DataFrame pour inclure uniquement les lignes où la colonne spécifiée contient le mot-clé donné.

    Parameters:
    - database (pd.DataFrame): Le DataFrame à filtrer.
    - column_name (str): Nom de la colonne sur laquelle appliquer le filtre. Par défaut, 'Molecule'.
    - keyword (str): Le mot-clé à rechercher dans la colonne spécifiée. Par défaut, 'Isopropyl'.

    Returns:
    - pd.DataFrame: Un nouveau DataFrame contenant uniquement les lignes où la colonne spécifiée contient le mot-clé.

    Example:
    >>> db = pd.DataFrame({'Molecule': ['Isopropyl alcohol', 'Ethanol', 'Isopropyl acetate', 'Water']})
    >>> filtered_db = filter_for_keyword(db)
    >>> print(filtered_db)
               Molecule
    0  Isopropyl alcohol
    2  Isopropyl acetate
    """
    database_contains_Isopropyl= database[database[column_name].str.contains(keyword, na=False)]
    return database_contains_Isopropyl



def filter_NaOH(database_filtered_1, column_name="Molecule", keyword="Sodium Hydroxide"):
    """
    Filtre le DataFrame pour inclure uniquement les lignes où la colonne spécifiée contient la chaîne spécifiée.

    Parameters:
    - database (pd.DataFrame): Le DataFrame à filtrer.
    - column_name (str): Nom de la colonne sur laquelle appliquer le critère de recherche. Par défaut, 'Molecule'.
    - keyword (str): La chaîne de caractères à rechercher dans la colonne spécifiée. Par défaut, 'Sodium Hydroxide'.

    Returns:
    - pd.DataFrame: Un nouveau DataFrame contenant uniquement les lignes où la colonne spécifiée contient la chaîne de recherche.

    Example:
    >>> db = pd.DataFrame({'Molecule': ['Sodium Hydroxide', 'Water', 'Sodium Chloride', 'Sodium Hydroxide solution']})
    >>> filtered_db = filter_contains(db)
    >>> print(filtered_db)
                         Molecule
    0           Sodium Hydroxide
    3  Sodium Hydroxide solution
    """
    database_contains_NaOH= database_filtered_1[database_filtered_1[column_name].str.contains(keyword, na=False)]
    return database_contains_NaOH