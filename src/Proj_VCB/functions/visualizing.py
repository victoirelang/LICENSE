import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG, display

def visualize_molecules_for_cream(df, cream_name):
    """
    Affiche les molécules d'une certaine catégorie de crèmes avec des atomes spécifiques mis en évidence.
    
    Parameters:
    df (pd.DataFrame): DataFrame contenant les données.
    cream_name (str): Nom de la crème à visualiser.
    """
    
    filtered_df = df[df['Cream'] == cream_name]
    
    filtered_df = filtered_df.dropna(subset=['Smiles'])
    
    filtered_df = filtered_df[filtered_df['Smiles'].apply(lambda x: isinstance(x, str))]
    
    print(f"Total molécules pour la crème {cream_name}: {len(filtered_df)}")
    
    
    color_map = {
        'Isopropyl': (0, 1, 0),  
        'NaOH': (0, 0, 1),       
        'Linalool': (0.25, 0, 0.25),   
        'Citronellol': (0.4, 0, 0.4), 
        'Limonene': (0.6, 0, 0.6),   
        'Benzyl Alcohol': (0.8, 0, 0.8), 
        'Benzyl Salicylate': (1, 0, 1)   
    }
    
    
    smiles_patterns = {
        'Isopropyl': 'CC(C)O',
        'NaOH': '[Na+].[OH-]',
        'Linalool': 'CC(=C)CCC1=CC=C(C=C1)O',
        'Citronellol': 'CC(C)CCC(C)CO',
        'Limonene': 'CC1=CC[C@H](C=C1)C(C)C',
        'Benzyl Alcohol': 'c1ccccc1CO',
        'Benzyl Salicylate': 'COC1=CC=CC=C1C(=O)OCC2=CC=CC=C2'


    }
    mols = []
    atom_colors = []
    
    
    for index, row in filtered_df.iterrows():
        smi = row['Smiles']
        mol = Chem.MolFromSmiles(smi)
        
        # Ignorer les SMILES invalides
        if mol is None:
            continue
        
        highlight_dict = {}
        
        for compound, smiles in smiles_patterns.items():
            pattern = Chem.MolFromSmiles(smiles)  
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                print(f"Match trouvé pour {compound} dans la molécule {smi}")
            for match in matches:
                for idx in match:
                    highlight_dict[idx] = color_map[compound]
        
        
        mols.append(mol)
        atom_colors.append(highlight_dict)
    
    
    for mol, highlights in zip(mols, atom_colors):
        drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
        drawer.drawOptions().useBWAtomPalette()
        drawer.DrawMolecule(mol, highlightAtoms=list(highlights.keys()), highlightAtomColors=highlights)
        drawer.FinishDrawing()
        
        
        svg = drawer.GetDrawingText().replace('svg:', '')
        display(SVG(svg))





