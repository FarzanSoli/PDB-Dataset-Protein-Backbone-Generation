"""
Author: Farzan Soleymani
Date: Dec 01-2023
"""
import numpy as np
# =============================================================================
class features():
    def __init__(self):
        # ================================================================================= #
        #                                   Amino acid parameters                           #
        # ================================================================================= #
        # Helix prob = alpha | Sheet prob = beta |
        # Hydrophobicity = H_1 | Hydrophilicity = H_2 |
        # Polarity = P_1 | Polarizability = P_2 | isoelectric_PH = P_i |
        # side chain net charge number = NCN | solvent accessible surface area = SASA |
        self.ALA = {'P_1':8.1, 'P_2':0.046,'vol':1.00,'H_1':0.62, 'H_2':-0.5,'P_i': 6.11,'alpha':0.42,'beta':0.23,'NCN':0.007187,'SASA':1.181}
        self.CYS = {'P_1':5.50,'P_2':0.128,'vol':2.43,'H_1':0.29, 'H_2':-1.0,'P_i': 6.35,'alpha':0.17,'beta':0.41,'NCN':-0.03661,'SASA':1.461}
        self.ASP = {'P_1':13.0,'P_2':0.105,'vol':2.78,'H_1':-0.9, 'H_2':3.0, 'P_i': 2.95,'alpha':0.25,'beta':0.20,'NCN':-0.02382,'SASA':1.587}
        self.GLU = {'P_1':12.3,'P_2':0.151,'vol':3.78,'H_1':-0.74,'H_2':3.0, 'P_i': 3.09,'alpha':0.42,'beta':0.21,'NCN':0.006802,'SASA':1.862}
        self.PHE = {'P_1':5.20,'P_2':0.29, 'vol':5.89,'H_1':1.19, 'H_2':-2.5,'P_i': 5.67,'alpha':0.30,'beta':0.38,'NCN':0.037552,'SASA':2.228}
        self.GLY = {'P_1':9.0, 'P_2':0.00, 'vol':0.00,'H_1':0.48, 'H_2':0.0, 'P_i': 6.07,'alpha':0.13,'beta':0.15,'NCN':0.179052,'SASA':0.881}
        self.HIS = {'P_1':20.4,'P_2':0.23, 'vol':4.66,'H_1':-0.4, 'H_2':-0.5,'P_i': 7.69,'alpha':0.27,'beta':0.30,'NCN':-0.01069,'SASA':2.025}
        self.ILE = {'P_1':5.20,'P_2':0.186,'vol':4.00,'H_1':1.38, 'H_2':-1.8,'P_i': 6.04,'alpha':0.30,'beta':0.45,'NCN':0.021631,'SASA':1.810}
        self.LYS = {'P_1':11.3,'P_2':0.219,'vol':4.77,'H_1':-1.5, 'H_2':3.0, 'P_i': 9.99,'alpha':0.32,'beta':0.27,'NCN':0.017708,'SASA':2.258}
        self.LEU = {'P_1':4.90,'P_2':0.186,'vol':4.00,'H_1':1.06, 'H_2':-1.8,'P_i': 6.04,'alpha':0.39,'beta':0.31,'NCN':0.051672,'SASA':1.931}
        self.MET = {'P_1':5.70,'P_2':0.221,'vol':4.43,'H_1': 0.64,'H_2':-1.3,'P_i': 5.71,'alpha':0.38,'beta':0.32,'NCN':0.002683,'SASA':2.034}
        self.ASN = {'P_1':11.6,'P_2':0.134,'vol':2.95,'H_1':-0.78,'H_2':2.0, 'P_i': 6.52,'alpha':0.21,'beta':0.22,'NCN':0.005392,'SASA':1.655}
        self.PRO = {'P_1':8.0, 'P_2':0.131,'vol':2.72,'H_1': 0.12,'H_2':0.0, 'P_i': 6.80,'alpha':0.13,'beta':0.34,'NCN':0.239530,'SASA':1.468}
        self.GLN = {'P_1':10.5,'P_2':0.180,'vol':3.95,'H_1':-0.85,'H_2':0.2, 'P_i': 5.65,'alpha':0.36,'beta':0.25,'NCN':0.049211,'SASA':1.932}
        self.ARG = {'P_1':10.5,'P_2':0.291,'vol':6.13,'H_1':-2.53,'H_2':3.0, 'P_i':10.74,'alpha':0.36,'beta':0.25,'NCN':0.043587,'SASA':2.560}
        self.SER = {'P_1':9.20,'P_2':0.062,'vol':1.60,'H_1':-0.18,'H_2':0.3, 'P_i': 5.70,'alpha':0.20,'beta':0.28,'NCN':0.004627,'SASA':1.298}
        self.THR = {'P_1':8.60,'P_2':0.108,'vol':2.60,'H_1':-0.05,'H_2':-0.4,'P_i': 5.60,'alpha':0.21,'beta':0.36,'NCN':0.003352,'SASA':1.525}
        self.VAL = {'P_1':5.90,'P_2':0.140,'vol':3.00,'H_1':1.08, 'H_2':-1.5,'P_i': 6.02,'alpha':0.27,'beta':0.49,'NCN':0.057004,'SASA':1.645}
        self.TRP = {'P_1':5.40,'P_2':0.409,'vol':8.08,'H_1':0.81, 'H_2':-3.4,'P_i': 5.94,'alpha':0.32,'beta':0.42,'NCN':0.037977,'SASA':2.663}
        self.TYR = {'P_1':6.20,'P_2':0.298,'vol':6.47,'H_1':0.26, 'H_2':-2.3,'P_i': 5.66,'alpha':0.25,'beta':0.41,'NCN':0.023599,'SASA':2.368}
        # ================================================================================= #
        self.AA_prop_keys = ['P_1','P_2','vol','H_1','H_2','P_i','alpha','beta','NCN','SASA']
        # --------------------------------------------------------------------------------- #
        self.AAs = [self.ALA, self.CYS, self.ASP, self.GLU, self.PHE, self.GLY, self.HIS,
                    self.ILE, self.LYS, self.LEU, self.MET, self.ASN, self.PRO, self.GLN, 
                    self.ARG, self.SER, self.THR, self.VAL, self.TRP,self.TYR]
        # --------------------------------------------------------------------------------- #
        p_1, p_2, vol = [],[],[]
        H_1, H_2, pi = [],[],[]
        alpha, beta, ncn, sasa = [],[],[],[]
        for i in range(len(self.AAs)):
            pi.append(self.AAs[i]['P_i'])
            p_1.append(self.AAs[i]['P_1'])
            p_2.append(self.AAs[i]['P_2'])
            vol.append(self.AAs[i]['vol'])
            H_1.append(self.AAs[i]['H_1'])
            H_2.append(self.AAs[i]['H_2'])
            ncn.append(self.AAs[i]['NCN'])
            sasa.append(self.AAs[i]['SASA'])
            beta.append(self.AAs[i]['beta'])
            alpha.append(self.AAs[i]['alpha'])
        # --------------------------------------------------------------------------------- #
        """ The unknown amino acid features are calcualted as the median of standard amino acid 
            features. """
        self.UNK = {'P_1':np.median(p_1),'P_2':np.median(p_2), 
                    'vol':np.median(vol)
                    ,'H_1':np.median(H_1),'H_2':np.median(H_2),
                    'P_i':np.median(pi) 
                    ,'SASA':np.median(sasa)
                    ,'alpha':np.median(alpha),'beta':np.median(beta),
                    'NCN':np.median(ncn)}
        # ================================================================================= #
        self.AA_dict = {'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7, 'K':8,'L':9,'M':10,
                       'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19,'X':20}
        # --------------------------------------------------------------------------------- #
        self.Amino_acids = {'A':self.ALA,'C':self.CYS,'D':self.ASP,'E':self.GLU,'F':self.PHE,
                            'G':self.GLY,'H':self.HIS,'I':self.ILE,'K':self.LYS,'L':self.LEU,
                            'M':self.MET,'N':self.ASN,'P':self.PRO,'Q':self.GLN,'R':self.ARG,
                            'S':self.SER,'T':self.THR,'V':self.VAL,'W':self.TRP,
                            'Y':self.TYR,'X': self.UNK}
        # --------------------------------------------------------------------------------- #
        
        