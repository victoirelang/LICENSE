import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
def highlight_atoms(mol, highlight_dict):
    drawer = Draw.MolDraw2DCairo(300, 300)
    atom_colors = {idx: color for idx, color in highlight_dict.items()}
    atom_indices = list(highlight_dict.keys())
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=atom_indices, highlightAtomColors=atom_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
def visualize_molecules_for_cream(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    filtered_df = df[df['Cream'] == cream_name]
    filtered_df = filtered_df.dropna(subset=['SMILES'])  
    filtered_df = filtered_df[filtered_df['SMILES'].apply(lambda x: isinstance(x, str))]  
    
    
    color_map = {
        'Isopropyl': (0, 1, 0),  # Vert
        'NaOH': (0, 0, 1),       # Bleu
        'Linalool': (1, 0.5, 0), # Orange 1
        'Citronellol': (1, 0.6, 0), # Orange 2
        'Limonene': (1, 0.7, 0),   # Orange 3
        'Benzyl Alcohol': (1, 0.8, 0), # Orange 4
        'Benzyl Salicylate': (1, 0.9, 0) # Orange 5
    }
    
    smarts_patterns = {
        'Isopropyl': 'CC(C)O',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(C)CC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)(C)CO',
        'Limonene': 'CC1=CC=CCC1(C)C',
        'Benzyl Alcohol': 'C1=CC=C(C=C1)CO',
        'Benzyl Salicylate': 'C1=CC=C(C=C1)OC(=O)C2=CC=CC=C2'
    }
    
    for index, row in filtered_df.iterrows():
        smi = row['SMILES']
        mol = Chem.MolFromSmiles(smi)
        
        highlight_dict = {}
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        if mol:
            img_data = highlight_atoms(mol, highlight_dict)
            display(Draw.MolToImage(mol, highlightAtoms=highlight_dict.keys(), highlightAtomColors=highlight_dict))
            with open(f"molecule_{index}.png", "wb") as f:
                f.write(img_data)
        highlight_dict = {}
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        if mol:
            img_data = highlight_atoms(mol, highlight_dict)
            display(Draw.MolToImage(mol, highlightAtoms=highlight_dict.keys(), highlightAtomColors=highlight_dict))
            with open(f"molecule_{index}.png", "wb") as f:
                f.write(img_data)

