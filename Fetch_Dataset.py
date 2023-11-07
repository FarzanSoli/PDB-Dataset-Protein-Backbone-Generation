""" ########## Processing C-alpha files ########## """
import wget
import pickle
import requests
import numpy as np
import pandas as pd 
from lxml import etree
from io import StringIO
import os, sys, gemmi, json
from Functions import Functions
from scipy.spatial import distance
# ========================================= #
import pdbreader
# ========================================= #

directory = os.getcwd() + '\Dataset\PDB_alpha_C'
if not os.path.exists(directory):
    os.mkdir(directory)
    print('The new directory is created')

""" ########## Download PDB files ########## """
class Download_PDB:
    def __init__(self, main_url):
        super().__init__()

        self.main_url = main_url
        self.Pop_list = ["/pub/pdb/data/biounit/PDB/divided/", 
                    "https://www.wwpdb.org/ftp/pdb-ftp-sites", "?C=M;O=A","?C=S;O=A","?C=N;O=D"]
        self.Pop_dict = self.Pop_list[1:]
    # ------------------------------------------ #
    def getLinks(self, url):
        print("Getting links from: " + url)
        session = requests.Session()
        page = session.get(url)
        html = page.content.decode("utf-8")
        tree = etree.parse(StringIO(html), parser=etree.HTMLParser())
        refs = tree.xpath("//a")    
        return list(set([link.get('href', '') for link in refs]))
    # ------------------------------------------ #
    def getsubLinks(self,url):
        print("Getting links from: " + url)
        session = requests.Session()
        page = session.get(url)
        html = page.content.decode("utf-8")
        tree = etree.parse(StringIO(html), parser=etree.HTMLParser())
        refs = tree.xpath("//a")    
        return list(set([link.get('href', '') for link in refs]))
    # ------------------------------------------ #
    def Download(self,directory):
        urls = self.getLinks(main_url)
        urls.remove('/pub/pdb/data/biounit/PDB/')
        remove_links = []
        for i in self.Pop_dict:
            remove_links.append(main_url+i)
        Links = {}
        for url in urls:
            Suburl = []
            for suburl in self.getsubLinks(main_url+url):
                Suburl.append(suburl)
            for item in self.Pop_list:
                if item in Suburl:
                    Suburl.remove(item)
            Links[url] = Suburl
            if url in self.Pop_dict:
                del Links[url]
        # ------------------------------------------
        for url, codes in Links.items():
            for code in codes:
                if main_url+url+code in remove_links:
                    continue
                else:
                    wget.download(main_url+url+code, directory)
    # ------------------------------------------
main_url = 'https://files.wwpdb.org/pub/pdb/data/biounit/PDB/divided/'
LINKS = Download_PDB(main_url).Download(directory)
Extract_dataset = Functions(directory).gz_extract()

directory = os.getcwd() + '/Dataset/PDB'
""" ########## Extract_Coordinates ########## """
def C_a_dist_mat(directory):
    extensions = [".pdb"+str(i) for i in range(1,10)]
    os.chdir(directory)
    C_a = {}
    chain = {}
    AA_sq = {}
    protein = {}
    Distance = {}
    for item in os.listdir(directory):
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
                    dist = {}
                    for ch in chain[id_]:
                        ca[ch] = protein[id_]['ATOM'].loc[
                            protein[id_]['ATOM']['name'] == 'CA',
                            ['chain','x','y','z','resname']].loc[
                            protein[id_]['ATOM'].loc[protein[id_]['ATOM']['name'] == 'CA',
                            ['chain','x','y','z','resname']]['chain']==ch].drop('chain', axis=1)
                        if ca[ch].shape[0] <= 5 or ca[ch].shape[0] >= 512:
                            continue
                        else:
                            seq[ch] = list(ca[ch]['resname'])
                            c_a[ch] = ca[ch].drop('resname', axis=1).to_numpy()
                            dist[ch] = distance.cdist(c_a[ch], c_a[ch],'euclidean')
                    # ---------------------------------------------------------- #   
                    if len(c_a) != 0:
                        C_a[id_] = c_a
                        AA_sq[id_] = seq
                        Distance[id_] = dist
    # -------------------------------------------------------------------------- #   
    Info = dict()
    # Info['Coordinate'] = C_a
    Info['AA_sq'] = AA_sq
    Info['Distance'] = Distance
    return Info

# ================================== #
Proteins = C_a_dist_mat(directory)
f = open("PDB.pkl","wb")
pickle.dump(Proteins,f)
f.close()

