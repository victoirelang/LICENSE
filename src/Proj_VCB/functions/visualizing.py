import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from IPython.display import display
from PIL import Image

def highlight_atoms(mol, highlight_dict):
    drawer = Draw.MolDraw2DCairo(300, 300)
    atom_colors = {idx: color for idx, color in highlight_dict.items()}
    atom_indices = list(highlight_dict.keys())
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=atom_indices, highlightAtomColors=atom_colors)
    drawer.FinishDrawing()
    return Image.open(drawer.GetDrawingText())


def visualize_molecules_for_cream(df, cream_name):
    filtered_df = df[df['Cream'] == cream_name]
    filtered_df = filtered_df.dropna(subset=['Smiles']) 
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))] 
    
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
    
    mols = []
    legends = []
    
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        if mol is None:
            print(f"Erreur de parsing SMILES pour: {smi}")
            continue
        
        highlight_dict = {}
        for compound, smarts in smarts_patterns.items():
            pattern = Chem.MolFromSmarts(smarts)
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        if mol:
            legends.append(row['Smiles'])  # Ajouter le SMILES en tant que légende
            mols.append((mol, highlight_dict))
    
    mols_to_draw = [m[0] for m in mols]
    highlight_dicts = [m[1] for m in mols]
    img = Draw.MolsToGridImage(mols_to_draw, molsPerRow=4, subImgSize=(200, 200),
                               legends=legends, highlightAtomLists=[list(h.keys()) for h in highlight_dicts],
                               highlightAtomColors=[h for h in highlight_dicts])
    
    display(img)
    img.save("molecules_grid.png")

