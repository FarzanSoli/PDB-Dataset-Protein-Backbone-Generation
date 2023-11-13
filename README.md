# PDB Dataset-Protein Backbone Generation
 This repository contains the codes to retrieve 3-D coordinates of the protein backbone from [PDB database][1]. Additionally, You can find the pre-processed dataset containing the [alpha-carbon coordinates.][2]. You can find the amino acid sequences of the proteins in the AA_Seq_main.zip file. 
 
The "AA_features.py" code includes physicochemical features of the 20 standard amino acids. Additionally, the physicochemical features of the unknown amino acid (X) is computed using the median of the known features. Each amino acid can be encoded using their respective [physicochemical features][3]. 

In order to download the alpha-carbon coordinates, run the code "Fetch_Dataset.py". 

"Protein_Backbone_Dataset.py" and "Functions.py" codes contain the required code to extract alpha-carbon coordinates and compute their respective distance matrices.


[1]: https://files.wwpdb.org/pub/pdb/data/biounit/PDB/divided/

[2]: https://uottawa-my.sharepoint.com/personal/fsole078_uottawa_ca/_layouts/15/guestaccess.aspx?docid=074356a5d40844683a4c31fea17cfa56a&authkey=AZdy7nF7VP0wl_IGgEkcxc8&e=iOswpm

[3]: https://www.sciencedirect.com/science/article/pii/S2001037023000296
