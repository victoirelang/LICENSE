import pandas as pd
def filter_pregnancy_risk(data):
    dataPR = data[data["Pregnancy Risk"] != "No"]
    print("Data base with molecules dangerous for pregnant women")
    return dataPR
    
def found_percentage_molecules_danger(dictionnaire_molecules, dictionnaire_riskymolecules):
    percentage_molecules_danger = {}
    for creme, total in dictionnaire_molecules.items():
        if creme in dictionnaire_riskymolecules:
            risque = dictionnaire_riskymolecules[creme]
            pourcentage = (risque / total) * 100
            pourcentage_arrondi = round(pourcentage, 2)  
        else:
            pourcentage_arrondi = 0  
        
        pourcentage_molecules_risquees[creme] = f"{pourcentage_arrondi}%"

    return pourcentage_molecules_risquees
