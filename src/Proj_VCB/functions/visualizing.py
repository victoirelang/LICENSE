import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
from rdkit.Chem.Draw import rdMolDraw2D

def visualize_molecules_for_cream1(df, cream_name):
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
        'Linalool': (0.6, 0, 0.6),   # Violet 1
        'Citronellol': (0.7, 0, 0.7), # Violet 2
        'Limonene': (0.8, 0, 0.8),   # Violet 3
        'Benzyl Alcohol': (0.9, 0, 0.9), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
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
    atom_colors = []
    
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
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    # Générer la grille d'images sans légendes et avec des images plus grandes
    img = Draw.MolsToGridImage(mols, molsPerRow=1, subImgSize=(1000, 1000),
                               legends=None, highlightAtomLists=[list(color.keys()) for color in atom_colors],
                               highlightAtomColors=atom_colors)
    
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



import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG


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
        'Linalool': (0.6, 0, 0.6),   # Violet 1
        'Citronellol': (0.7, 0, 0.7), # Violet 2
        'Limonene': (0.8, 0, 0.8),   # Violet 3
        'Benzyl Alcohol': (0.9, 0, 0.9), # Violet 4
        'Benzyl Salicylate': (1, 0, 1)   # Violet 5
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
    atom_colors = []

    # Parcourir les molécules filtrées et préparer les images
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                raise ValueError("Mol is None")
        except Exception as e:
            print(f"Erreur de parsing SMILES pour: {smi}, erreur: {e}")
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
        mols.append(mol)
        atom_colors.append(highlight_dict)

    # Définir un style de dessin personnalisé pour s'assurer que tous les atomes non surlignés sont en noir
    options = Draw.MolDrawOptions()
    options.atomPalette = {6: (0, 0, 0), 7: (0, 0, 0), 8: (0, 0, 0), 1: (0, 0, 0)}  # Carbone, Azote, Oxygène, Hydrogène en noir
    
    # Utiliser rdMolDraw2D pour dessiner les molécules avec des options personnalisées
    drawer = rdMolDraw2D.MolDraw2DSVG(500 * 3, 500 * ((len(mols) + 2) // 3), 500, 500)
    drawer.SetDrawOptions(options)
    drawer.DrawMolecules(mols, highlightAtomLists=[list(color.keys()) for color in atom_colors],
                         highlightAtomColors=atom_colors)
    drawer.FinishDrawing()
    
    # Afficher l'image de la grille
    svg = drawer.GetDrawingText().replace('svg:', '')
    display(SVG(svg))