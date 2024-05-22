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
        
        # Vérifier si la molécule est valide
        if mol is None:
            print(f"Erreur de parsing SMILES pour: {smi}")
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
        if mol:
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

# Exemple d'utilisation
import os
import pandas as pd

# Charger la base de données
current_dir = os.path.dirname(os.path.abspath(__file__))
database_path = os.path.join(current_dir, '../../data/database.csv')
df = pd.read_csv(database_path, sep=';')

# Définir le nom de la crème
cream_name = 'Nuxe'  # Remplacez par le nom de la crème que vous souhaitez visualiser

# Appeler la fonction de visualisation
visualize_molecules_for_cream(df, cream_name)
