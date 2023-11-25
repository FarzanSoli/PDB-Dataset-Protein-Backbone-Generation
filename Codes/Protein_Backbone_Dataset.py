"""
Author: Farzan Soleymani
Date: Dec 01-2023
"""
""" ########## Calculating Protein backbone distance matrices ########## """
import os
import torch
import pickle
import numpy as np
import pandas as pd 
from Functions import Functions
from scipy.spatial import distance
# ========================================= #
class C_a_Distance_Matrix():
    def __init__(self, Pad_Length, Directory, Pickle_file):
        super().__init__()
        self.directory = Directory
        self.Pad_Length = Pad_Length
        self.padding = Functions("/Dataset/").padding
        self.Pickle_file = Pickle_file
        self.encode_CT = Functions("/Dataset/").encode_CT
        self.files = os.listdir(os.getcwd() + self.directory)
        self.standardize = Functions("/Dataset/").standardize
    # ========================================= #
    #      Distance Matrix of CSV dataset       #
    # ========================================= #
        """ 
            Calculate distance matrix given a dataset saved as ".csv".
        """
    def Distance_Matrix_CSV(self):
        Coordinates = {}
        Protein_backbone = {}
        Dist_Protein_backbone = {}
        for file in self.files: 
            coordinate = pd.read_excel(self.directory +'/'+file, 
                            names=['X_coordinate', 'Y_coordinate', 'Z_coordinate']).to_numpy()
                # -------------------------------------------- #
            if coordinate.shape[0] > self.Pad_Length:
                cut_coordinate = coordinate[:self.Pad_Length, :]
                # cut_coordinate = self.standardize(cut_coordinate)
                
                Coordinates[file.replace('.xlsx','')] = cut_coordinate
                # -------------------------------------------- #
                Protein_backbone[file.replace('.xlsx','')] = cut_coordinate
                Dist_Protein_backbone[file.replace('.xlsx','')] = distance.cdist(
                            self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                            self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                            'euclidean')
                # -------------------------------------------- #
            else:
                coordinates_ = self.padding(self.Pad_Length, coordinate)
                # coordinates_ = self.standardize(coordinates_)                
                Coordinates[file.replace('.xlsx','')] = coordinates_
                # -------------------------------------------- #
                Protein_backbone[file.replace('.xlsx','')] = Coordinates[file.replace('.xlsx','')][:,:]
                Dist_Protein_backbone[file.replace('.xlsx','')] = distance.cdist(
                            self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                            self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                            'euclidean')
                
        with open('Dataset/Distance_Protein_Backbone.pkl', 'wb') as file:
            pickle.dump(Dist_Protein_backbone, file)
        with open('Dataset/Protein_Backbone.pkl', 'wb') as file:
            pickle.dump(Protein_backbone, file)
    # ========================================= #
    #     Distance Matrix of Pickle dataset     #
    # ========================================= #
        """ 
            Calculate distance matrix given a dataset saved as ".pkl".
        """
    def Distance_Matrix_PKL(self):
        Coordinates = {}
        Protein_backbone = {}
        Dist_Protein_backbone = {}
        with(self.Pickle_file,'rb') as file:
            coordinate = pickle.load(file)
            file.close()
        # -------------------------------------------- #
        if coordinate.shape[0] > self.Pad_Length:
            cut_coordinate = coordinate[:self.Pad_Length, :]
            # cut_coordinate = self.standardize(cut_coordinate)
            Coordinates[file.replace('.xlsx','')] = cut_coordinate
            # -------------------------------------------- #
            Protein_backbone[file.replace('.xlsx','')] = cut_coordinate
            Dist_Protein_backbone[file.replace('.xlsx','')] = distance.cdist(
                        self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                        self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                        'euclidean')
        # -------------------------------------------- #
        else:
            coordinates_ = self.padding(self.Pad_Length, coordinate)
            # coordinates_ = self.standardize(coordinates_)                
            Coordinates[file.replace('.xlsx','')] = coordinates_
            # -------------------------------------------- #
            Protein_backbone[file.replace('.xlsx','')] = Coordinates[file.replace('.xlsx','')][:,:]
            Dist_Protein_backbone[file.replace('.xlsx','')] = distance.cdist(
                        self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                        self.padding(self.Pad_Length, Protein_backbone[file.replace('.xlsx','')]), 
                        'euclidean')
            
        with open('Dataset/Distance_Protein_Backbone.pkl', 'wb') as file:
            pickle.dump(Dist_Protein_backbone, file)
        with open('Dataset/Protein_Backbone.pkl', 'wb') as file:
            pickle.dump(Protein_backbone, file)
    # ========================================= #
    #         Normalize Distance Matrix         #
    # ========================================= #
    def Dist_Norm(self, distance_matrix):
        return (distance_matrix - np.min(distance_matrix))/(
                    np.max(distance_matrix) - np.min(distance_matrix))
    # ========================================= #
    #         Encoded Protein Backbone Seq      #
    # ========================================= #
    def Encoded_Backbone_Seq(self, Pad_Length, Directory):
        # Directory = 'Dataset/AA_Seq_main.csv'
        AA_dataset = pd.read_csv(Directory)
        Encoded_AA = self.encode_CT(Pad_Length, AA_dataset)
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
if __name__ == "__main__":
    """ Save the distance matrix coordinate for a given backbone length. """
    Pad_Length = 128
    data_dir = os.getcwd() +'/Dataset/PDB_alpha_C'
    files = os.listdir(data_dir)
    # For CSV dataset
    Load_Data = C_a_Distance_Matrix(Pad_Length, data_dir).Distance_Matrix_CSV()
    # For Pickle dataset
    Load_Data = C_a_Distance_Matrix(Pad_Length, data_dir).Distance_Matrix_PKL()
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
    """ Save encoded protein backbone AA sequence. """
    Pad_Length = 64
    Directory = 'Dataset/AA_Seq_main.csv'
    Encoded_Backbone_Seq = C_a_Distance_Matrix(Directory).Encoded_Backbone_Seq(Pad_Length, Directory)
    
    # ------------ Padding encoded backbone ------------ #
    """ Padding the encoded protein backbone AA sequence. """
    file = open('Dataset/Encoded_Backbone_Seq_64.pkl', 'rb')
    Dataset = pickle.load(file)
    file.close()
    Padding_Protein_Backbone = C_a_Distance_Matrix(Directory).Padding_Protein_Backbone(Pad_Length, Dataset)
