"""
Author: Farzan Soleymani
Date: Dec 01-2023
"""
""" ########## Essential functions ########## """
import os
import copy
import torch
import zipfile
import numpy as np
import pandas as pd 
import gzip, shutil
from AA_features import features
# ========================================= #
class Functions():
    def __init__(self, directory):
        super().__init__()
        self.directory = directory
    # ========================================= #
    #                Fix AA sequenc             #
    # ========================================= #
    """ 
        Filling missing parts of the given sequence.
    """
    def missing_seq(self, Dataset):
        df = Dataset.loc[Dataset['Sequence'].str.contains('X') & 
                           Dataset['Sequence'].str.startswith('X') &
                           Dataset['Sequence'].str.endswith('X')]
        for i in ['A','C','D','E','F','G','H','I','K','L',
                  'M','N','P','Q','R','S','T','V','W','Y']:
            kn_seq = df.loc[Dataset['Sequence'].str.contains(i)]
        result = pd.concat([df, kn_seq]).drop_duplicates(keep=False)
        return result
    # ========================================= #
    #                 Read Dataset              #
    # ========================================= #
    def Dataset_Reader(self, dataset):
        Dataset = pd.concat([pd.read_csv('Dataset/'+dataset), 
                             self.missing_seq(pd.read_csv(
                            'Dataset/'+dataset))]).drop_duplicates(keep=False)
        return Dataset
    # ========================================= #
    #                   Padding                 #
    # ========================================= #
    """ 
        padding given sequence to a specific length.
    """
    
    def padding(self, pad_len, matrix):
        mat = np.zeros((pad_len, 3))
        mat[:matrix.shape[0], :matrix.shape[-1]] = matrix
        return mat
    # ========================================= #
    #                  Unzip files              #
    # ========================================= #
    """ 
        Unzipping files in a given directory.
    """
    def unzip(self, directory, file):
        # file = 'PDB alpha-C.zip'
        # directory = "/Dataset/"
        with zipfile.ZipFile(os.getcwd()+directory+file, 'r') as zip_ref:
            zip_ref.extractall(os.getcwd()+directory)
        
    # ========================================= #
    #                 Extract gzip              #
    # ========================================= #
    """ 
        Extract Gzip files.
    """
    def gz_extract(self):
        extension = ".gz"
        os.chdir(self.directory)
        for item in os.listdir(self.directory): # loop through items in dir
          if item.endswith(extension): # check for ".gz" extension
              gz_name = os.path.abspath(item) # get full path of files
              # get file name for file within -> removes '.cif'
              file_name = (os.path.basename(gz_name)).rsplit('.',1)[0] 
              with gzip.open(gz_name, "rb") as f_in, open(file_name, "wb") as f_out:
                  # Copy the contents of source file to destination file
                  shutil.copyfileobj(f_in, f_out)
              os.remove(gz_name) # delete zipped file
    # ========================================= #
    #               Encoding Amio acids         #
    # ========================================= #
    """
        Normalizing amino acid physicochemical features.
    """
    def Normalize_AA(self):
        props = []
        for t in range(len(features().AA_prop_keys)):
            prop = []
            for i in range(len(list(features().AA_dict))):
                prop.append(features().Amino_acids[
                    list(features().AA_dict)[i]][features().AA_prop_keys[t]])
            props.append(np.array(prop))
        # ------------------------------------- #
        Norm_props = []
        for t in range(len(features().AA_prop_keys)):
            Norm_props.append(
                (props[t]-np.min(props[t]))/(np.max(props[t])-np.min(props[t])))
        # ------------------------------------- #
        AA_props = []
        for t in range(len(list(features().AA_dict))):
            aa_props = []
            for i in range(len(features().AA_prop_keys)):
                aa_props.append(Norm_props[i][t])
            AA_props.append(aa_props)
        return Norm_props, AA_props
    # ========================================= #
    #               Amino acid Encoding         #
    # ========================================= #
    """
        Encode amino acids based on their physicochemical features.
    """
    def encode_CT(self, Pad_Length, dataframe):
        encoding = dict(zip(list(features().AA_dict), self.Normalize_AA()[1]))
        # ------------------------------------- #
        Encoded_AA = {}
        for index, row in dataframe.iterrows():
            Encoded_AA[row['PDB_ID']] = np.array([encoding[c.upper()] for c in row['Sequence']])
        # ------------------------------------- #
        encoded_AA = np.zeros((Pad_Length, len(features().AA_prop_keys)))
        Encoded_AA_padded = {}
        for key, value in Encoded_AA.items():
            if value.shape[0] > Pad_Length:
                encoded_AA[:,:] = value[:Pad_Length,:]
            else: 
                encoded_AA[:value.shape[0],:] = value
            Encoded_AA_padded[key] = encoded_AA
        return Encoded_AA_padded
    # ========================================= #
    #                 READ Fasta file           #
    # ========================================= #
    """
        Read ".fasta" file and extract dataset.
    """
    def read_fasta(self, fasta_file, comment="#"):
        # with gzip.open(fasta_file,"r") as f:
        with open(fasta_file, "r") as file:
            id_ = None
            seq_id = []
            sequence = []
            sequences = []
            # loop through each line in the file
            for line in file:
                # If the line starts with a ">" character, it is a new sequence identifier
                if line.startswith(">"):
                    # If this is not the first sequence, print the previous one
                    if id_ is not None:
                        seq_id.append(id_)
                        sequences.append(''.join(sequence))
                    # Get the new sequence identifier and reset the sequence variable
                    id_ = line.strip()[1:]
                    sequence = []
                # Otherwise, it is part of the sequence, so append it to the sequence variable
                else:
                    sequence.append(line.upper())
            if id_ is not None:
                seq_id.append(id_)
                sequences.append(''.join(sequence))
            return list(zip(seq_id, sequences))
    # ========================================= #
    #                 Normalize                 #
    # ========================================= #
    """
        Normalize given alpha-carbon coordinates. 
    """
    def Normalize(self, coordinates):
        x_norm = (coordinates[:,0] - np.min(coordinates[:,0]))/(np.max(coordinates[:,0]) - np.min(coordinates[:,0]))
        y_norm = (coordinates[:,1] - np.min(coordinates[:,1]))/(np.max(coordinates[:,1]) - np.min(coordinates[:,1]))
        z_norm = (coordinates[:,2] - np.min(coordinates[:,2]))/(np.max(coordinates[:,2]) - np.min(coordinates[:,2]))
        Normalized_coordinate = np.transpose(np.array([x_norm, y_norm, z_norm]))

        return Normalized_coordinate
    # ========================================= #
    #                 Standardize               #
    # ========================================= #
    """
        Standardize given alpha-carbon coordinates. 
    """
    def standardize(self, coordinates):
        x_standard = (coordinates[:,0] - coordinates[:,0].mean())/(np.std(coordinates[:,0]))
        y_standard = (coordinates[:,1] - coordinates[:,1].mean())/(np.std(coordinates[:,1]))
        z_standard = (coordinates[:,2] - coordinates[:,2].mean())/(np.std(coordinates[:,2]))
        standard_coordinate = np.transpose(np.array([x_standard, y_standard, z_standard]))
        return standard_coordinate


if __name__ == "__main__": 
    padding = Functions("/Dataset/").padding
    encode_CT = Functions("/Dataset/").encode_CT
    Normalize = Functions("/Dataset/").Normalize
    standardize = Functions("/Dataset/").standardize