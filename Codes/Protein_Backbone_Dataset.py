"""
Author: Farzan Soleymani
Date: Nov 20-2023
"""
""" ########## Processing Protein backbone coordinates and distance matrices ########## """
import os
import torch
import pickle
import numpy as np
import pandas as pd 
from Functions import Functions
from scipy.spatial import distance
# ========================================= #
padding = Functions("/Dataset/").padding
encode_CT = Functions("/Dataset/").encode_CT
standardize = Functions("/Dataset/").standardize
# ========================================= #
class Protein_Backbone():
    def __init__(self, directory):
        super().__init__()
        self.directory = directory
    # ========================================= #
    #                 Distance Matrix           #
    # ========================================= #
    def Distance_Matrix(self, Pad_Length, Directory, files):
        Coordinates = {}
        Protein_backbone = {}
        Dist_Protein_backbone = {}
        for file in files: 
            coordinate = pd.read_excel(Directory +'/'+file, 
                            names=['X_coordinate', 'Y_coordinate', 'Z_coordinate']).to_numpy()
                # -------------------------------------------- #
            if coordinate.shape[0] > Pad_Length:
                cut_coordinate = coordinate[:Pad_Length, :]
                # cut_coordinate = self.standardize(cut_coordinate)
                
                Coordinates[file.replace('.xlsx','')] = cut_coordinate
                # -------------------------------------------- #
                Protein_backbone[file.replace('.xlsx','')] = cut_coordinate
                Dist_Protein_backbone[file.replace('.xlsx','')] = distance.cdist(
                                padding(Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                                padding(Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                                'euclidean')
                # -------------------------------------------- #
            else:
                coordinates_ = padding(Pad_Length, coordinate)
                # coordinates_ = self.standardize(coordinates_)                
                Coordinates[file.replace('.xlsx','')] = coordinates_
                # -------------------------------------------- #
                Protein_backbone[file.replace('.xlsx','')] = Coordinates[file.replace('.xlsx','')][:,:]
                Dist_Protein_backbone[file.replace('.xlsx','')] = distance.cdist(
                                padding(Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                                padding(Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                                'euclidean')
                
        with open('Dataset/Distance_Protein_Backbone.pkl', 'wb') as file:
            pickle.dump(Dist_Protein_backbone, file)
        with open('Dataset/Protein_Backbone.pkl', 'wb') as file:
            pickle.dump(Protein_backbone, file)
    # ========================================= #
    #         Normalize Distance Matrix         #
    # ========================================= #
    def Dist_Norm(self, distance_matrix):
        return (distance_matrix - np.min(distance_matrix))/(np.max(distance_matrix) - np.min(distance_matrix))
    # ========================================= #
    #         Encoded Protein Backbone Seq      #
    # ========================================= #
    def Encoded_Backbone_Seq(self, Pad_Length, Directory):
        # Directory = 'Dataset/AA_Seq_main.csv'
        AA_dataset = pd.read_csv(Directory)
        Encoded_AA = encode_CT(Pad_Length, AA_dataset)
        with open('Dataset/Encoded_Backbone_Seq.pkl', 'wb') as file:
            pickle.dump(Encoded_AA, file)
    # ========================================= #
    #           Padding Protein backbone        #
    # ========================================= #
    def Padding_Protein_Backbone(self, Pad_Length, Encoded_AA):
        # Pad_Length = 128
        Backbone_Seq = {}
        for key, value in Encoded_AA.items():
            Backbone_Seq[key] = np.zeros((Pad_Length, value.shape[-1]))
            Backbone_Seq[key][:value[:Pad_Length].shape[0],
                           :value[:Pad_Length].shape[-1]] = value[:Pad_Length]
        # -----------------------------------------------------
        with open('Dataset/Padded_Backbone_Seq.pkl', 'wb') as file:
            pickle.dump(Backbone_Seq, file)

# ------------ calculate distance matrix ------------ #
Pad_Length = 128
data_dir = os.getcwd() +'/Dataset/PDB_alpha_C'
files = os.listdir(data_dir)
Load_Data = Protein_Backbone(data_dir).Distance_Matrix(Pad_Length, data_dir, files)

# ------------ Proteins ------------ #
with open('Dataset/Distance_Protein_Backbone.pkl', 'rb') as file:
    Proteins = pickle.load(file)
    file.close()
Proteins_PDB_ID = Proteins.keys()
Proteins_dist = []
for key, value in Proteins.items():
    value_scaled = Protein_Backbone(data_dir).Dist_Norm(value)
    Proteins_dist.append(torch.FloatTensor(value_scaled[None,:,:]))

torch.save(Proteins_dist, 'Dataset/Proteins_Tensor.pt')
torch.save(list(Proteins.keys()), 'Dataset/Proteins_PDB_ID.pt')

# ------------ Encoded Backbone ------------ #
Pad_Length = 64
Directory = 'Dataset/AA_Seq_main.csv'
Encoded_Backbone_Seq = Protein_Backbone(Directory).Encoded_Backbone_Seq(Pad_Length, Directory)

# ------------ Padding encoded backbone ------------ #
file = open('Dataset/Encoded_Backbone_Seq_64.pkl', 'rb')
Dataset = pickle.load(file)
file.close()
Padding_Protein_Backbone = Protein_Backbone(Directory).Padding_Protein_Backbone(Pad_Length, Dataset)
