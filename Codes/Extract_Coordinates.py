"""
Author: Farzan Soleymani
Date: Dec 01-2023
"""
import os
import pickle
import pdbreader
import pandas as pd 
from tqdm import tqdm

""" ########## Extract_Coordinates:
        This code takes the directory where the PDB files are saved and 
        extracts alpha-carbon Coordinates and AA sequences. The output is a dictionary
        with "AA_sq" and "Coordinate" keys and it will be saved as pickle file.
    ########## """
def Extract_Coordinates(directory):
    extensions = [".pdb"+str(i) for i in range(1,32)]
    os.chdir(directory)
    chain = {}
    protein = {}
    AA_Chain = {}
    Ca_Chain = {}
    for item in tqdm(os.listdir(directory)):
        for extension in extensions:
            if item.endswith(extension):
                # Parse and get basic information
                id_ = item.replace(extension, "")
                protein[id_] = pdbreader.read_pdb(item)
                if 'ATOM' not in protein[id_].keys():
                    continue
                else:
                    chain[id_] = set(list(protein[id_]['ATOM']['chain']))
                    ca = {}
                    c_a = {}
                    seq = {}
                    for ch in chain[id_]:
                        ca[ch] = protein[id_]['ATOM'].loc[
                            protein[id_]['ATOM']['name'] == 'CA',
                            ['chain','x','y','z','resname']].loc[
                            protein[id_]['ATOM'].loc[protein[id_]['ATOM']['name'] == 'CA',
                            ['chain','x','y','z','resname']]['chain']==ch].drop('chain', axis=1)
                        # Filter sequence with length bigger than 5 and less than 1024
                        if ca[ch].shape[0] <= 3:
                            continue
                        else:
                            seq[ch] = list(ca[ch]['resname'])
                            c_a[ch] = ca[ch].drop('resname', axis=1).to_numpy()
                    # ---------------------------------------------------------- #   
                        ID = id_+'_'+ch
                        Ca_Chain[ID] = c_a[ch]
                        AA_Chain[ID] = seq[ch]
    # -------------------------------------------------------------------------- #   
    Info = dict()
    Info['AA_sq'] = AA_Chain
    Info['Coordinate'] = Ca_Chain
    return Info

""" ########## Save the Dataset ########## """
if __name__ == "__main__":
    print('Saving the coordinate of alpha carbons!')
    directory = os.getcwd() + '/Dataset/PDB_alpha_C'
    Proteins = Extract_Coordinates(directory)
    with open("PDB_Backbone.pkl","wb") as file:
        pickle.dump(Proteins,file)
        file.close()
