import pandas as pd
def filter_pregnancy_risk(data):
    dataPR = data[data["Pregnancy Risk"] != "No"]
    print("Data base with molecules dangerous for pregnant women")
    return dataPR
    
import pandas as pd

def dangerous_pregnancy(database):
    dataPR = database[database["Pregnancy Risk"] != "No"]
    nombre_molecules_par_creme = database.groupby('Cream')['Molecule'].nunique()
    dictionnaire_molecules = nombre_molecules_par_creme.to_dict()
    nb_risky_molecules = dataPR.groupby('Cream')['Molecule'].nunique()
    dictionnaire_riskymolecules = nb_risky_molecules.to_dict()
    dictionnaire_complet = {cream: 0 for cream in database['Cream'].unique()}
    dictionnaire_complet.update(dictionnaire_riskymolecules)
    pourcentage_molecules_risquees = {}
    for creme, total in dictionnaire_molecules.items():
        if creme in dictionnaire_riskymolecules:
            risque = dictionnaire_riskymolecules[creme]
            pourcentage = (risque / total) * 100
            pourcentage_arrondi = round(pourcentage, 2)
        else:
            pourcentage_arrondi = 0
        pourcentage_molecules_risquees[creme] = f"{pourcentage_arrondi}%"
    print("Pourcentage de molécules risquées par crème:")
    for creme, pourcentage in pourcentage_molecules_risquees.items():
        print(f"{creme}: {pourcentage}")
    return pourcentage_molecules_risquees

url = "https://raw.githubusercontent.com/blanchebillarant/ppchem/main/project/data.csv"
database = pd.read_csv(url, sep=';')
resultat = adangerous_pregnancy(database)
