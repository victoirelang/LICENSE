import pandas as pd
def filter_pregnancy_risk(data):
    dataPR = data[data["Pregnancy Risk"] != "No"]
    print("Data base with molecules dangerous for pregnant women")
    return dataPR
    

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
        pourcentage_molecules_risquees[creme] = pourcentage_arrondi

    sorted_cremes = sorted(pourcentage_molecules_risquees.items(), key=lambda x: x[1])
    result = {creme: f"{pourcentage}%" for creme, pourcentage in sorted_cremes}

    return result

    print("Voici le classement des crèmes pour le critère de pregnancy:")
    print(result)
    





