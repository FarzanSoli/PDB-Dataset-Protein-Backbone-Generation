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
class C_a_Distance_Matrix():
    def __init__(self, Pad_Length, Directory):
        super().__init__()
        self.directory = Directory
        self.Pad_Length = Pad_Length
        self.padding = Functions("/Dataset/").padding
        self.files = os.listdir(os.getcwd() + self.directory)
    def Distance_Matrix(self):
        Motif = {}
        Scaffold = {}
        Dist_motif = {}
        Coordinates = {}
        Dist_scaffold = {}
        for file in self.files: 
            coordinate = pd.read_excel(self.directory+'/'+file, names=['X_coordinate', 'Y_coordinate', 'Z_coordinate']).to_numpy()
            # --------------------------------- #
            if coordinate.shape[0] > self.Pad_Length:
                cut_coordinate = coordinate[:self.Pad_Length, :]
                Coordinates[file.replace('.xlsx','')] = cut_coordinate
                Motif[file.replace('.xlsx','')] = cut_coordinate[:int(cut_coordinate.shape[0]*0.3),:]
                Dist_motif[file.replace('.xlsx','')] = distance.cdist(self.padding(int(Motif[file.replace('.xlsx','')].shape[0])+1,
                                                                                       Motif[file.replace('.xlsx','')]),
                                                                      self.padding(int(Motif[file.replace('.xlsx','')].shape[0])+1, 
                                                                                       Motif[file.replace('.xlsx','')]),
                                                                         'euclidean')
                Scaffold[file.replace('.xlsx','')] = cut_coordinate[Motif[file.replace('.xlsx','')].shape[0]:, :]
                Dist_scaffold[file.replace('.xlsx','')] = distance.cdist(self.padding(self.Pad_Length, Scaffold[file.replace('.xlsx','')]), 
                                                                         self.padding(self.Pad_Length, Scaffold[file.replace('.xlsx','')]), 
                                                                         'euclidean')
            # --------------------------------- #
            else:
                Coordinates[file.replace('.xlsx','')] = self.padding(self.Pad_Length, coordinate)
                Motif[file.replace('.xlsx','')] = Coordinates[file.replace('.xlsx','')][:int(Coordinates[file.replace('.xlsx','')].shape[0]*0.3),:]
                Dist_motif[file.replace('.xlsx','')] = distance.cdist(self.padding(int(Motif[file.replace('.xlsx','')].shape[0])+1,
                                                                                       Motif[file.replace('.xlsx','')]),
                                                                      self.padding(int(Motif[file.replace('.xlsx','')].shape[0])+1, 
                                                                                       Motif[file.replace('.xlsx','')]),
                                                                      'euclidean')
                Scaffold[file.replace('.xlsx','')] = Coordinates[file.replace('.xlsx','')][Motif[file.replace('.xlsx','')].shape[0]:, :]
                Dist_scaffold[file.replace('.xlsx','')] = distance.cdist(self.padding(self.Pad_Length, Scaffold[file.replace('.xlsx','')]), 
                                                                         self.padding(self.Pad_Length, Scaffold[file.replace('.xlsx','')]), 
                                                                         'euclidean')
        with open('Dataset/Distance_Scaffold.pkl', 'wb') as file:
            pickle.dump(Dist_scaffold, file)
        with open('Dataset/Distance_Motif.pkl', 'wb') as file:
                pickle.dump(Dist_motif, file)
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
Load_Data = C_a_Distance_Matrix(Pad_Length, data_dir).Distance_Matrix()

# ------------ Proteins ------------ #
with open('Dataset/Distance_C_a_Distance_Matrix.pkl', 'rb') as file:
    Proteins = pickle.load(file)
    file.close()
Proteins_PDB_ID = Proteins.keys()
Proteins_dist = []
for key, value in Proteins.items():
    value_scaled = C_a_Distance_Matrix(data_dir).Dist_Norm(value)
    Proteins_dist.append(torch.FloatTensor(value_scaled[None,:,:]))

torch.save(Proteins_dist, 'Dataset/Proteins_Tensor.pt')
torch.save(list(Proteins.keys()), 'Dataset/Proteins_PDB_ID.pt')

# ------------ Encoded Backbone ------------ #
Pad_Length = 64
Directory = 'Dataset/AA_Seq_main.csv'
Encoded_Backbone_Seq = C_a_Distance_Matrix(Directory).Encoded_Backbone_Seq(Pad_Length, Directory)

# ------------ Padding encoded backbone ------------ #
file = open('Dataset/Encoded_Backbone_Seq_64.pkl', 'rb')
Dataset = pickle.load(file)
file.close()
Padding_Protein_Backbone = C_a_Distance_Matrix(Directory).Padding_Protein_Backbone(Pad_Length, Dataset)
