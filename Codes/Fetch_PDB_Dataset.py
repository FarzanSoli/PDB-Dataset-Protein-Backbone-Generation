""" ########## Processing C-alpha files ########## """
import os
import wget
import requests
import pandas as pd 
from lxml import etree
from io import StringIO
from Functions import Functions
# ========================================= #

# ========================================= #

directory = os.getcwd() + '\Dataset\PDB_alpha_C'
if not os.path.exists(directory):
    os.mkdir(directory)
    print('The new directory is created')

""" ########## Download PDB files 
        A dataset is downloaded from PDB database which contains alpha-carbon coordinates and 
        AA seqeuences among other featues. 
    ########## """
class Download_PDB:
    def __init__(self, main_url):
        super().__init__()

        self.main_url = main_url
        self.Pop_list = ["/pub/pdb/data/biounit/PDB/divided/", 
                         "https://www.wwpdb.org/ftp/pdb-ftp-sites", 
                         "?C=M;O=A","?C=S;O=A","?C=N;O=D"]
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
if __name__ == "__main__":
    main_url = 'https://files.wwpdb.org/pub/pdb/data/biounit/PDB/divided/'
    LINKS = Download_PDB(main_url).Download(directory)
    Extract_dataset = Functions(directory).gz_extract()

