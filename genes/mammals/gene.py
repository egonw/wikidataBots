#!usr/bin/env python
# -*- coding: utf-8 -*-

'''
Authors: 
  Andra Waagmeester (andra' at ' micelio.be)

This file is part of ProteinBoxBot.

ProteinBoxBot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ProteinBoxBot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ProteinBoxBot.  If not, see <http://www.gnu.org/licenses/>.
'''

__author__ = 'Andra Waagmeester'
__license__ = 'GPL'

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../ProteinBoxBot_Core")
import PBB_Core
import PBB_Debug
import PBB_login
import PBB_settings
import PBB_Functions
import ProteinBoxBotKnowledge
import urllib
import urllib3
import certifi
import copy
import traceback
import sys
import mygene_info_settings

from raven import Client

try:
    import simplejson as json
except ImportError as e:
    import json
    
"""
This is the human-genome specific part of the ProteinBoxBot. Its purpose is to enrich Wikidata with
human gene and known external identifiers.
  
"""
    
class genome(object):
    def __init__(self, object):
        counter = 0
        self.genomeInfo = object
        print "Getting all {} genes in Wikidata".format(self.genomeInfo["name"])
        self.content = json.loads(self.download_genes(self.genomeInfo["name"]))
        self.gene_count = self.content["total"]
        self.genes = self.content["hits"]
        self.logincreds = PBB_login.WDLogin(PBB_settings.getWikiDataUser(), PBB_settings.getWikiDataPassword())

        entrezWikidataIds = dict()
        wdqQuery = "CLAIM[703:{}] AND CLAIM[351]".format(self.genomeInfo["wdid"].replace("Q", ""))
        InWikiData = PBB_Core.WDItemList(wdqQuery, "351")
        
        '''
        Below a mapping is created between entrez gene ids and wikidata identifiers.
        '''
        for geneItem in InWikiData.wditems["props"]["351"]:
            entrezWikidataIds[str(geneItem[2])] = geneItem[0]

        for gene in self.genes:
          try:         
                if str(gene["entrezgene"]) in entrezWikidataIds.keys():
                   gene["wdid"] = 'Q'+str(entrezWikidataIds[str(gene["entrezgene"])])
                else:
                   gene["wdid"] = None 
                gene["logincreds"] = self.logincreds
                gene["genomeInfo"] = self.genomeInfo
                geneClass = mammal_gene(gene)
                if str(geneClass.entrezgene) in entrezWikidataIds.keys():
                    geneClass.wdid = 'Q'+str(entrezWikidataIds[str(geneClass.entrezgene)])
                    print geneClass.wdid
                else: 
                    geneClass.wdid = None 
                counter = counter +1
                if counter == 100:
                    self.logincreds = PBB_login.WDLogin(PBB_settings.getWikiDataUser(), PBB_settings.getWikiDataPassword())

                
          except:
              print "There has been an except"
              print "Unexpected error:", sys.exc_info()[0]

              f = open('/tmp/exceptions_{}.txt'.format(self.genomeInfo["name"]), 'a')
              f.write(str(gene["entrezgene"])+"\n")
              traceback.print_exc(file=f)
              f.close()
              

    def download_genes(self, species):
        """
        Downloads the latest list of human genes from mygene.info through the URL specified in mygene_info_settings
        """
        urllib.urlretrieve (mygene_info_settings.getGenesUrl().format(species), species+"_genes.json")
        file = open(species+"_genes.json", 'r')
        return file.read()
        
class mammal_gene(object):
    def __init__(self, object):
        self.genomeInfo = object["genomeInfo"]
        self.content = object
        self.entrezgene = object["entrezgene"]
        self.name = object["name"]
        self.logincreds = object["logincreds"]
        gene_annotations = json.loads(self.annotate_gene())
        print "gene_annotations"
        PBB_Debug.prettyPrint(gene_annotations)
        print "object"
        PBB_Debug.prettyPrint(object)
        self.annotationstimestamp = gene_annotations["_timestamp"]
        self.wdid = object["wdid"]
        
        # symbol
        if isinstance(gene_annotations["symbol"], list):
             self.symbol = object["symbol"]
        else:
             self.symbol = [object["symbol"]]
        
        # HGNC
        if "HGNC" in gene_annotations:
            if isinstance(gene_annotations["HGNC"], list): 
                self.hgnc = gene_annotations["HGNC"]
            else:
                self.hgnc = [gene_annotations["HGNC"]]
        else:
            self.hgnc = None
            
        # Ensembl Gene & transcript
        if "ensembl" in gene_annotations:
            if "gene" in gene_annotations["ensembl"]:
                if isinstance(gene_annotations["ensembl"]["gene"], list): 
                    self.ensembl_gene = gene_annotations["ensembl"]["gene"]
                else:
                    self.ensembl_gene = [gene_annotations["ensembl"]["gene"]]
            else:
                self.ensembl_gene = None
            
            if "transcript" in gene_annotations["ensembl"]:
                if isinstance(gene_annotations["ensembl"]["transcript"], list): 
                    self.ensembl_transcript = gene_annotations["ensembl"]["transcript"]
                else:
                    self.ensembl_transcript = [gene_annotations["ensembl"]["transcript"]]
            else:
                self.ensembl_transcript = None
        # Homologene
        if "homologene" in gene_annotations:
            if isinstance(gene_annotations["homologene"]["id"], list): 
                self.homologene = [str(i) for i in gene_annotations["homologene"]["id"]]
            else:
                self.homologene = [str(gene_annotations["homologene"]["id"])]
        else:
            self.homologene = None
        # Refseq 
        if "refseq" in gene_annotations:
            if "rna" in gene_annotations["refseq"]:
                if isinstance(gene_annotations["refseq"]["rna"], list): 
                    self.refseq_rna = gene_annotations["refseq"]["rna"]
                else:
                    self.refseq_rna = [gene_annotations["refseq"]["rna"]]
            else :
                self.refseq_rna = None
        else :
            self.refseq_rna = None       
        if "genomic_pos" in gene_annotations:
            self.genomic_pos =[]
            if (isinstance(gene_annotations["genomic_pos"], list)):
                for i in range(len(gene_annotations["genomic_pos"])):
                    if gene_annotations["genomic_pos"][i]["chr"] in ProteinBoxBotKnowledge.chromosomes.keys():
                           self.genomic_pos.append(ProteinBoxBotKnowledge.chromosomes[gene_annotations["genomic_pos"][i]["chr"]])

            if isinstance(gene_annotations["genomic_pos"], list): 
                self.genomic_pos = gene_annotations["genomic_pos"]

            else:
                self.genomic_pos = [gene_annotations["genomic_pos"]]
        else:
            self.genomic_pos = None
        
        # type of Gene
        if "type_of_gene" in gene_annotations:
            f = open('/tmp/geneType.txt', 'a')
            self.type_of_gene = []
            f.write(gene_annotations["type_of_gene"]+"\n")
            f.close()
            if gene_annotations["type_of_gene"]=="ncRNA":
                self.type_of_gene.append("427087")
            if gene_annotations["type_of_gene"]=="snRNA":
                self.type_of_gene.append("284578") 
            if gene_annotations["type_of_gene"]=="snoRNA":
                self.type_of_gene.append("284416")
            if gene_annotations["type_of_gene"]=="rRNA":
                self.type_of_gene.append("215980")
            if gene_annotations["type_of_gene"]=="tRNA":
                self.type_of_gene.append("201448")
            if gene_annotations["type_of_gene"]=="pseudo":
                self.type_of_gene.append("277338") 
            if gene_annotations["type_of_gene"]=="protein-coding":
                self.type_of_gene.append("20747295")                          
        else:
            self.type_of_gene = None
        # Reference section           
        gene_reference = {
                    'ref_properties': [u'P248', u'P143', 'TIMESTAMP'],
                    'ref_values': [self.genomeInfo["release"], u'Q20641742' , 'TIMESTAMP']
                }
                   
        references = dict()
        
        data2add = dict()
        data2add["P279"] = ["7187"]
        references['P279'] = [copy.deepcopy(gene_reference)]
        data2add["P703"] = [self.genomeInfo["wdid"].replace("Q", "")]
        references['P703'] = [copy.deepcopy(gene_reference)]    
        data2add['P351'] = [str(self.entrezgene)]
        references['P351'] = [copy.deepcopy(gene_reference)]
        data2add['P353'] = self.symbol
        references['P353'] = [copy.deepcopy(gene_reference)]
        # references['P353'] = gene_reference
        
        if "type_of_gene" in vars(self):
            if self.type_of_gene != None:
                for i in range(len(self.type_of_gene)):
                    data2add["P279"].append(self.type_of_gene[i])
                    references['P279'].append(copy.deepcopy(gene_reference))
                    
        if "ensembl_gene" in vars(self):
            if self.ensembl_gene != None:
                data2add["P594"] = self.ensembl_gene
                references['P594'] = []
                for i in range(len(self.ensembl_gene)):
                    references['P594'].append(copy.deepcopy(gene_reference))
        if "ensembl_transcript" in vars(self):
            if self.ensembl_transcript != None:
                data2add['P704'] = self.ensembl_transcript
                references['P704'] = []
                for i in range(len(self.ensembl_transcript)):
                    references['P704'].append(copy.deepcopy(gene_reference))
        if "hgnc" in vars(self):
            if self.hgnc != None:
                data2add['P354'] = self.hgnc
                references['P354'] = []
                for i in range(len(self.hgnc)):
                    references['P354'].append(copy.deepcopy(gene_reference))
        if "homologene" in vars(self):
            if self.homologene != None:
                data2add['P593'] = self.homologene
                references['P593'] = []
                for i in range(len(self.homologene)):
                    references['P593'].append(copy.deepcopy(gene_reference))
        if "refseq_rna" in vars(self):
            if self.refseq_rna != None:
                data2add['P639'] = self.refseq_rna
                references['P639'] = []
                for i in range(len(self.refseq_rna)):
                    references['P639'].append(copy.deepcopy(gene_reference))
        if "genomic_pos" in object:
            if (isinstance(object["genomic_pos"], list)):
               chromosome = object["genomic_pos"][0]["chr"]
            else: chromosome = object["genomic_pos"]["chr"]
            data2add['P1057'] =  chromosomes[str(chromosome)]
            references['P1057'] = gene_reference    

        if "alias" in gene_annotations.keys():
            if isinstance(gene_annotations["alias"], list):
                self.synonyms = []
                for alias in gene_annotations["alias"]:
                    self.synonyms.append(alias)
            else:
               self.synonyms = [gene_annotations["alias"]]
            for syn in self.symbol:
               self.synonyms.append(syn)
            print self.synonyms
        else:
            self.synonyms = None
        if self.wdid != None:   
            wdPage = PBB_Core.WDItemEngine(self.wdid, self.name, data = data2add, server="www.wikidata.org", references=references, domain="genes")
            wdPage.set_description(description=self.genomeInfo['name']+' gene', lang='en')
            if self.synonyms != None:
                wdPage.set_aliases(aliases=self.synonyms, lang='en', append=True)
            print self.wdid
            self.wd_json_representation = wdPage.get_wd_json_representation() 
            wdPage.write(self.logincreds)
        else:
            for key in data2add.keys():
                if len(data2add[key]) == 0:
                    data2add.pop(key, None)
            for key in references.keys():
                if len(references[key]) == 0:
                    references.pop(key, None)
            wdPage = PBB_Core.WDItemEngine(item_name=self.name, data=data2add, server="www.wikidata.org", references=references, domain="genes")
            wdPage.set_description(description=self.genomeInfo['name']+' gene', lang='en')
            if self.synonyms != None:
                wdPage.set_aliases(aliases=self.synonyms, lang='en', append=True)
            self.wd_json_representation = wdPage.get_wd_json_representation() 

            wdPage.write(self.logincreds)
               
    def annotate_gene(self):
        "Get gene annotations from mygene.info"
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        request = http.request("GET", mygene_info_settings.getGeneAnnotationsURL()+str(self.entrezgene))
        return request.data
        
        
 
        
        
        
        
        