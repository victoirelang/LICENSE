import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display

def visualize_molecules_for_cream(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    # Filtrer les molécules pour la catégorie de crème spécifiée
    filtered_df = df[df['Cream'] == cream_name]
    # Supprimer les lignes avec des valeurs manquantes dans 'Smiles'
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    # Supprimer les lignes avec des valeurs non-string dans 'Smiles'
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    # Définir les couleurs pour les atomes spécifiques
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (1, 0.5, 0), # Orange 1
        'Citronellol': (1, 0.6, 0), # Orange 2
        'Limonene': (1, 0.7, 0),   # Orange 3
        'Benzyl Alcohol': (1, 0.8, 0), # Orange 4
        'Benzyl Salicylate': (1, 0.9, 0) # Orange 5
    }
    
    # Définir les motifs SMARTS pour les composés spécifiques
    smarts_patterns = {
        'Isopropyl': 'CC(C)O',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(C)CC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)(C)CO',
        'Limonene': 'CC1=CC=CCC1(C)C',
        'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',
        'Benzyl Salicylate': 'C1=CC=C(C=C1)OC(=O)C2=CC=CC=C2'
    }
    
    mols = []
    legends = []
    
    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        # Ignorer les SMILES invalides
        if mol is None:
            continue
        
        highlight_dict = {}
        # Identifier les atomes à mettre en évidence
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        # Ajouter la molécule et le dictionnaire de surbrillance à la liste
        legends.append(row['Smiles'])  # Ajouter le SMILES en tant que légende
        mols.append((mol, highlight_dict))
    
    # Générer la grille d'images
    mols_to_draw = [m[0] for m in mols]  # Extraire les molécules
    highlight_dicts = [m[1] for m in mols]  # Extraire les dictionnaires de surbrillance
    # Créer une image de grille avec les molécules et les surbrillances
    img = Draw.MolsToGridImage(mols_to_draw, molsPerRow=4, subImgSize=(200, 200),
                               legends=legends, highlightAtomLists=[list(h.keys()) for h in highlight_dicts],
                               highlightAtomColors=[h for h in highlight_dicts])
    
    # Afficher l'image de la grille
    display(img)
import pandas as pd
from rdkit import Chem

def check_smiles_validity(df, smiles_column='Smiles'):
    """
    Vérifie la validité des chaînes SMILES dans un DataFrame.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    smiles_column (str): Nom de la colonne contenant les chaînes SMILES.
    
    Returns:
    dict: Dictionnaire avec les résultats de la vérification. Les clés sont les indices des lignes du DataFrame,
          et les valeurs sont des booléens indiquant si la chaîne SMILES est valide (True) ou non (False).
    """
    validity_dict = {}
    
    for index, row in df.iterrows():
        smi = row[smiles_column]
        if not isinstance(smi, str):
            validity_dict[index] = False
            continue
        
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            validity_dict[index] = False
        else:
            validity_dict[index] = True
    
    return validity_dict


