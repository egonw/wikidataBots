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


sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")
import ProteinBoxBot_Core.PBB_Core as PBB_Core
import ProteinBoxBot_Core.PBB_Debug as PBB_Debug
import ProteinBoxBot_Core.PBB_login as PBB_login
import ProteinBoxBot_Core.PBB_settings as PBB_settings
import genes.mammals.ProteinBoxBotKnowledge as ProteinBoxBotKnowledge
import requests
import copy
import traceback
from SPARQLWrapper import SPARQLWrapper, JSON




from time import gmtime, strftime
import time


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
        self.start = time.time()
        self.genomeInfo = object["speciesInfo"][object["species"]]
        self.speciesInfo = object["speciesInfo"]
        # print("Getting all {} genes in Wikidata".format(self.genomeInfo["name"]))
        self.content = self.download_genes(self.genomeInfo["taxid"])
        self.gene_count = self.content["total"]
        self.genes = self.content["hits"]
        self.logincreds = PBB_login.WDLogin(PBB_settings.getWikiDataUser(), os.environ['wikidataApi'])
        fast_run_base_filter = {'P351': '', 'P703': 'Q15978631'}
        fast_run = True
        entrezWikidataIds = dict()
        uniprotwikidataids = dict()


        wdqQuery = "CLAIM[703:{}] AND CLAIM[351]".format(self.genomeInfo["wdid"].replace("Q", ""))
        InWikiData = PBB_Core.WDItemList(wdqQuery, wdprop="351")

        sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
        sparql.setQuery("""
        PREFIX wd: <http://www.wikidata.org/entity/>
        PREFIX wdt: <http://www.wikidata.org/prop/direct/>
        select DISTINCT * where {
            {?item wdt:P351 ?ncbi_gene .}
        }
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

        for result in results["results"]["bindings"]:
            try:
                geneItem = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
                ncbi_geneValue = result["ncbi_gene"]["value"]
                entrezWikidataIds[str(ncbi_geneValue)] = ncbi_geneValue

            except Exception as e:
                print(traceback.format_exc())


        for gene in self.genes:
            try:
                if str(gene["entrezgene"]) in entrezWikidataIds.keys():
                    gene["wdid"] = 'Q' + str(entrezWikidataIds[str(gene["entrezgene"])])
                else:
                    gene["wdid"] = None
                gene["uniprotwikidataids"] = uniprotwikidataids
                gene["logincreds"] = self.logincreds
                gene["genomeInfo"] = self.genomeInfo
                gene["speciesInfo"] = self.speciesInfo
                gene["start"] = self.start
                gene["fast_run"] = fast_run
                gene["fast_run_filter"] = fast_run_base_filter

                geneClass = mammal_gene(gene)
                if str(geneClass.entrezgene) in entrezWikidataIds.keys():
                    geneClass.wdid = 'Q' + str(entrezWikidataIds[str(geneClass.entrezgene)])
                else:
                    geneClass.wdid = None
                counter = counter + 1
                if counter == 100:
                    self.logincreds = PBB_login.WDLogin(PBB_settings.getWikiDataUser(), os.environ['wikidataApi'])


            except Exception as e:
                print(traceback.format_exc())
                PBB_Core.WDItemEngine.log('ERROR', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
                        main_data_id=gene["entrezgene"],
                        exception_type=type(e),
                        message=e.__str__(),
                        wd_id='-',
                        duration=time.time() - self.start
                    ))

    def download_genes(self, species):
        """
        Downloads the latest list of human genes from mygene.info through the URL specified in mygene_info_settings
        """
        r = requests.get("http://mygene.info/v2/query?q=__all__&species={}&entrezonly=true&size=100000".format(species))
        return r.json()


class mammal_gene(object):
    def __init__(self, object):
        """
        :type self: object
        """
        self.start = object["start"]
        self.entrezgene = object["entrezgene"]
        self.uniprotwikidataids = object["uniprotwikidataids"]
        gene_annotations = self.annotate_gene()
        self.genomeInfo = object["speciesInfo"][str(gene_annotations['taxid'])]
        self.content = object
        self.name = gene_annotations["name"]
        self.logincreds = object["logincreds"]
        self.fast_run = object["fast_run"]
        self.fast_run_filter = object["fast_run_filter"]


        if "_timestamp" in gene_annotations.keys():
            self.annotationstimestamp = gene_annotations["_timestamp"]
        self.wdid = object["wdid"]

        # symbol:
        self.symbol = gene_annotations["symbol"]

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
            else:
                self.refseq_rna = None
        else:
            self.refseq_rna = None

            # MGI
        if "MGI" in gene_annotations:
            if isinstance(gene_annotations["MGI"], list):
                self.MGI = gene_annotations["MGI"]
            else:
                self.MGI = [gene_annotations["MGI"]]
        else:
            self.MGI = None

        self.chromosome = None
        self.startpost = None
        self.endpos = None
        if "genomic_pos" in gene_annotations:
            if isinstance(gene_annotations["genomic_pos"], list):
                self.chromosome = []
                self.startpos = []
                self.endpos = []
                self.strand = []
                for i in range(len(gene_annotations["genomic_pos"])):
                    if gene_annotations["genomic_pos"][i]["chr"] in ProteinBoxBotKnowledge.chromosomes[
                        self.genomeInfo["name"]].keys():
                        self.chromosome.append(ProteinBoxBotKnowledge.chromosomes[self.genomeInfo["name"]][
                                                   gene_annotations["genomic_pos"][i]["chr"]])
                        self.startpos.append(gene_annotations["genomic_pos"][i]["start"])
                        self.endpos.append(gene_annotations["genomic_pos"][i]["end"])
                        self.strand.append(gene_annotations["genomic_pos"][i]["strand"])
            else:
                self.chromosome = []
                self.startpos = []
                self.endpos = []
                self.strand = []
                if gene_annotations["genomic_pos"]["chr"] in ProteinBoxBotKnowledge.chromosomes[
                    self.genomeInfo["main_key"]].keys():
                    self.chromosome.append(ProteinBoxBotKnowledge.chromosomes[self.genomeInfo["main_key"]][
                                               gene_annotations["genomic_pos"]["chr"]])
                    self.startpos.append(gene_annotations["genomic_pos"]["start"])
                    self.endpos.append(gene_annotations["genomic_pos"]["end"])
                    self.strand.append(gene_annotations["genomic_pos"]["strand"])

        self.encodes = None


        self.chromosomeHg19 = None
        self.startposHg19 = None
        self.endposHg19 = None
        if "genomic_pos_hg19" in gene_annotations:
            if isinstance(gene_annotations["genomic_pos_hg19"], list):
                self.chromosomeHg19 = []
                self.startposHg19 = []
                self.endposHg19 = []
                for i in range(len(gene_annotations["genomic_pos_hg19"])):
                    if gene_annotations["genomic_pos_hg19"][i]["chr"] in ProteinBoxBotKnowledge.chromosomes[
                        self.genomeInfo["name"]].keys():
                        self.chromosomeHg19.append(ProteinBoxBotKnowledge.chromosomes[self.genomeInfo["name"]][
                                                       gene_annotations["genomic_pos_hg19"][i]["chr"]])
                        self.startposHg19.append(gene_annotations["genomic_pos_hg19"][i]["start"])
                        self.endposHg19.append(gene_annotations["genomic_pos_hg19"][i]["end"])
            else:
                self.chromosomeHg19 = []
                self.startposHg19 = []
                self.endposHg19 = []
                if gene_annotations["genomic_pos_hg19"]["chr"] in ProteinBoxBotKnowledge.chromosomes[
                    self.genomeInfo["name"]].keys():
                    self.chromosomeHg19.append(ProteinBoxBotKnowledge.chromosomes[self.genomeInfo["name"]][
                                                   gene_annotations["genomic_pos_hg19"]["chr"]])
                    self.startposHg19.append(gene_annotations["genomic_pos_hg19"]["start"])
                    self.endposHg19.append(gene_annotations["genomic_pos_hg19"]["end"])

        # type of Gene
        if "type_of_gene" in gene_annotations:
            self.type_of_gene = []
            if gene_annotations["type_of_gene"] == "ncRNA":
                self.type_of_gene.append("Q427087")
            if gene_annotations["type_of_gene"] == "snRNA":
                self.type_of_gene.append("Q284578")
            if gene_annotations["type_of_gene"] == "snoRNA":
                self.type_of_gene.append("Q284416")
            if gene_annotations["type_of_gene"] == "rRNA":
                self.type_of_gene.append("Q215980")
            if gene_annotations["type_of_gene"] == "tRNA":
                self.type_of_gene.append("Q201448")
            if gene_annotations["type_of_gene"] == "pseudo":
                self.type_of_gene.append("Q277338")
            if gene_annotations["type_of_gene"] == "protein-coding":
                self.type_of_gene.append("Q20747295")
        else:
            self.type_of_gene = None
        # Reference section  
        # Prepare references
        refStatedIn = PBB_Core.WDItemID(value=self.genomeInfo["release"], prop_nr='P248', is_reference=True)
        #refStatedIn.overwrite_references = True
        refImported = PBB_Core.WDItemID(value='Q20641742', prop_nr='P143', is_reference=True)
        #refImported.overwrite_references = True
        refEntrezId = PBB_Core.WDString(value=str(self.entrezgene), prop_nr='P351', is_reference=True)
        #refEntrezId.overwrite_references = True
        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
        refRetrieved = PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
        #refRetrieved.overwrite_references = True
        gene_reference = [refStatedIn, refEntrezId, refImported, refRetrieved]

        refMiriamStatedIn = PBB_Core.WDItemID(value="Q16335166", prop_nr='P248', is_reference=True)
        #refMiriamStatedIn.overwrite_references = True
        refMiriamUrl = PBB_Core.WDUrl(value="http://www.ebi.ac.uk/miriam/main/collections/MIR:00000069", prop_nr='P854', is_reference=True)
        #refMiriamUrl.overwrite_references = True
        miriam_reference = [refMiriamStatedIn, refMiriamUrl]

        refStatedInEnsembl = PBB_Core.WDItemID(value= 'Q21996330', prop_nr='P248', is_reference=True)
        #refStatedInEnsembl.overwrite_references = True
        refImportedEnsembl = PBB_Core.WDItemID(value='Q1344256', prop_nr='P143', is_reference=True)
        #refImportedEnsembl.overwrite_references = True

        ensembl_reference = [refStatedInEnsembl, refImportedEnsembl, refRetrieved]

        genomeBuildQualifier = PBB_Core.WDItemID(value=self.genomeInfo["genome_assembly"], prop_nr='P659',
                                                 is_qualifier=True)
        genomeBuildPreviousQualifier = PBB_Core.WDItemID(value=self.genomeInfo["genome_assembly_previous"],
                                                         prop_nr='P659', is_qualifier=True)

        prep = dict()
        prep['P2888'] = [PBB_Core.WDUrl(value="http://identifiers.org/ncbigene/"+str(self.entrezgene), prop_nr='P2888',
                                          references=[copy.deepcopy(miriam_reference)])]
        prep['P703'] = [PBB_Core.WDItemID(value=self.genomeInfo['wdid'], prop_nr='P703',
                                          references=[copy.deepcopy(gene_reference)])]
        if self.genomeInfo["name"] == "human":
            prep['P353'] = [
                PBB_Core.WDString(value=self.symbol, prop_nr='P353', references=[copy.deepcopy(gene_reference)])]
        prep['P351'] = [
            PBB_Core.WDExternalID(value=str(self.entrezgene), prop_nr='P351', references=[copy.deepcopy(gene_reference)])]

        prep['P279'] = [PBB_Core.WDItemID(value='Q7187', prop_nr='P279', references=[copy.deepcopy(gene_reference)])]
        if "type_of_gene" in vars(self):
            if self.type_of_gene != None:
                for i in range(len(self.type_of_gene)):
                    prep['P279'].append(PBB_Core.WDItemID(value=self.type_of_gene[i], prop_nr='P279',
                                                          references=[copy.deepcopy(gene_reference)]))

        if "ensembl_gene" in vars(self):
            prep['P594'] = []
            if self.ensembl_gene != None:
                for ensemblg in self.ensembl_gene:
                    prep['P594'].append(
                        PBB_Core.WDExternalID(value=ensemblg, prop_nr='P594', references=[copy.deepcopy(gene_reference)]))
            else:
                prep['P594'].append([PBB_Core.WDBaseDataType.delete_statement(prop_nr='P594')])



        if "ensembl_transcript" in vars(self):
            if self.ensembl_transcript != None:
                prep['P704'] = []
                for ensemblt in self.ensembl_transcript:
                    prep['P704'].append(
                        PBB_Core.WDExternalID(value=ensemblt, prop_nr='P704', references=[copy.deepcopy(gene_reference)]))


        if "hgnc" in vars(self):
            if self.hgnc != None:
                prep['P354'] = []
                for hugo in self.hgnc:
                    prep['P354'].append(
                        PBB_Core.WDExternalID(value=hugo, prop_nr='P354', references=[copy.deepcopy(gene_reference)]))

        if "homologene" in vars(self):
            if self.homologene != None:
                prep['P593'] = []
                for ortholog in self.homologene:
                    prep['P593'].append(
                        PBB_Core.WDExternalID(value=ortholog, prop_nr='P593', references=[copy.deepcopy(gene_reference)]))

        if "refseq_rna" in vars(self):
            if self.refseq_rna != None:
                prep['P639'] = []
                for refseq in self.refseq_rna:
                    prep['P639'].append(
                        PBB_Core.WDExternalID(value=refseq, prop_nr='P639', references=[copy.deepcopy(gene_reference)]))

        if "chromosome" in vars(self):
            prep['P1057'] = []
            if self.chromosome != None:
                for chrom in list(set(self.chromosome)):
                    prep['P1057'].append(
                        PBB_Core.WDItemID(value=chrom, prop_nr='P1057', references=[copy.deepcopy(gene_reference)]))

        if "startpos" in vars(self):
            if not 'P644' in prep.keys():
                prep['P644'] = []
            if self.startpos != None:
                for pos in self.startpos:
                    for strand in self.strand:
                        if strand == -1:
                            strandQualifier = PBB_Core.WDItemID(value='Q22809711', prop_nr='P2549',
                                                 is_qualifier=True)
                        if strand == 1:
                            strandQualifier = PBB_Core.WDItemID(value='Q22809680', prop_nr='P2549',
                                                 is_qualifier=True)

                    prep['P644'].append(
                        PBB_Core.WDString(value=str(pos), prop_nr='P644', references=[copy.deepcopy(ensembl_reference)],
                                          qualifiers=[copy.deepcopy(genomeBuildQualifier)]))
        if "endpos" in vars(self):
            if not 'P645' in prep.keys():
                prep['P645'] = []
            if self.endpos != None:
                for pos in self.endpos:
                    for strand in self.strand:
                        if strand == -1:
                            strandQualifier = PBB_Core.WDItemID(value='Q22809711', prop_nr='P2549',
                                                 is_qualifier=True)
                        if strand == 1:
                            strandQualifier = PBB_Core.WDItemID(value='Q22809680', prop_nr='P2549',
                                                 is_qualifier=True)

                    prep['P645'].append(
                        PBB_Core.WDString(value=str(pos), prop_nr='P645', references=[copy.deepcopy(ensembl_reference)],
                                          qualifiers=[copy.deepcopy(genomeBuildQualifier)]))

        if "strand" in vars(self):
            if not 'P2548' in prep.keys():
                prep['P2548'] = []
            for strand in self.strand:
                if strand == -1:
                    prep['P2548'].append(
                        PBB_Core.WDItemID(value='Q22809711', prop_nr='P2548', references=[copy.deepcopy(gene_reference)],
                                          qualifiers=[copy.deepcopy(genomeBuildQualifier)]))
                if strand == 1:
                    prep['P2548'].append(
                        PBB_Core.WDItemID(value='Q22809680', prop_nr='P2548', references=[copy.deepcopy(gene_reference)],
                                          qualifiers=[copy.deepcopy(genomeBuildQualifier)]))

        if "startposHg19" in vars(self):
            if not 'P644' in prep.keys():
                prep['P644'] = []
            if self.startposHg19 != None:
                for pos in self.startposHg19:
                    prep['P644'].append(
                        PBB_Core.WDString(value=str(pos), prop_nr='P644', references=[copy.deepcopy(ensembl_reference)],
                                          qualifiers=[copy.deepcopy(genomeBuildPreviousQualifier)]))
        if "endposHg19" in vars(self):
            if not 'P644' in prep.keys():
                prep['P645'] = []
            if self.endposHg19 != None:
                for pos in self.endposHg19:
                    prep['P645'].append(
                        PBB_Core.WDString(value=str(pos), prop_nr='P645', references=[copy.deepcopy(ensembl_reference)],
                                          qualifiers=[copy.deepcopy(genomeBuildPreviousQualifier)]))

        if "MGI" in vars(self):
            prep['P671'] = []
            if self.MGI != None:
                for mgi in self.MGI:
                    prep['P671'].append(PBB_Core.WDExternalID(value=mgi, prop_nr='P671',
                                        references=[copy.deepcopy(gene_reference)]))

        if "alias" in gene_annotations.keys():
            if isinstance(gene_annotations["alias"], list):
                self.synonyms = []
                for alias in gene_annotations["alias"]:
                    self.synonyms.append(alias)
            else:
                self.synonyms = [gene_annotations["alias"]]
            self.synonyms.append(self.name)
        else:
            self.synonyms = None

        data2add = []
        for key in prep.keys():
            for statement in prep[key]:
                data2add.append(statement)

        if self.wdid != None:
          # if self.encodes != None:
            wdPage = PBB_Core.WDItemEngine(self.wdid, data=data2add, server="www.wikidata.org",
                                           domain="genes", fast_run=self.fast_run, fast_run_base_filter=self.fast_run_filter)

            wdPage.set_label(self.symbol, lang='en')
            wdPage.set_label(self.symbol, lang='nl')
            wdPage.set_label(self.symbol, lang='fr')
            wdPage.set_label(self.symbol, lang='de')
            wdPage.set_label(self.symbol, lang='srn')
            wdPage.set_label(self.symbol, lang='es')
            wdPage.set_label(self.symbol, lang='pt')
            wdPage.set_label(self.symbol, lang='sv')

            if wdPage.get_description(lang='en') == "" or wdPage.get_description(lang='en') == "human gene":
                wdPage.set_description(description="gene of the species "+self.genomeInfo['name'], lang='en')
            if wdPage.get_description(lang='fr') == "" or wdPage.get_description(lang='fr') == "gène humain" or wdPage.get_description(lang='fr') == "gène humaine" or wdPage.get_description(lang='fr') == "gene humain":
                wdPage.set_description(description="gène de l'espèce " + self.genomeInfo['name'], lang='fr')
            if wdPage.get_description(lang='nl') == "" or wdPage.get_description(lang='nl') == "menselijk gen" or wdPage.get_description(lang='nl') == "menselijkgen" or wdPage.get_description(lang='nl') == "Een Menselijk gen":
                wdPage.set_description(description="gen van de soort "+self.genomeInfo['name'], lang='nl')
            if wdPage.get_description(lang='de') == "" or wdPage.get_description(lang='de') == "gen":
                wdPage.set_description(description="Gen der Spezies " + self.genomeInfo['name'], lang='de')
            if wdPage.get_description(lang='es') == "":
                wdPage.set_description(description="gen de la especie "+self.genomeInfo['name'], lang='es')
            if wdPage.get_description(lang='pt') == "":
                wdPage.set_description(description="gene da espécie "+self.genomeInfo['name'], lang='pt')
            if wdPage.get_description(lang='sv') == "":
                wdPage.set_description(description="Genen från arten "+self.genomeInfo['name'], lang='sv')
            if wdPage.get_description(lang='srn') == "":
                wdPage.set_description(description="gen fu a sortu "+self.genomeInfo['name'], lang='srn')

            if self.synonyms != None:
                wdPage.set_aliases(aliases=self.synonyms, lang='en', append=True)
            self.wd_json_representation = wdPage.get_wd_json_representation()
            PBB_Debug.prettyPrint(self.wd_json_representation)
            PBB_Debug.prettyPrint(data2add)
            # print(self.wd_json_representation)

            wdPage.write(self.logincreds)

        else:
            wdPage = PBB_Core.WDItemEngine(item_name=self.name, data=data2add, server="www.wikidata.org",
                                           domain="genes", fast_run=self.fast_run, fast_run_base_filter=self.fast_run_filter)
            wdPage.set_label(self.symbol, lang='en')
            wdPage.set_label(self.symbol, lang='nl')
            wdPage.set_label(self.symbol, lang='fr')
            wdPage.set_label(self.symbol, lang='de')
            wdPage.set_label(self.symbol, lang='srn')
            wdPage.set_label(self.symbol, lang='es')
            wdPage.set_label(self.symbol, lang='pt')
            wdPage.set_label(self.symbol, lang='sv')

            if wdPage.get_description(lang='en') == "" or wdPage.get_description(lang='en') == "human gene":
                wdPage.set_description(description="gene of the species "+self.genomeInfo['name'], lang='en')
            if wdPage.get_description(lang='fr') == "" or wdPage.get_description(lang='fr') == "gène humain" or wdPage.get_description(lang='fr') == "gène humaine":
                wdPage.set_description(description="gène de l'espèce " + self.genomeInfo['name'], lang='fr')
            if wdPage.get_description(lang='nl') == "" or wdPage.get_description(lang='nl') == "menselijk gen":
                wdPage.set_description(description="gen van de soort "+self.genomeInfo['name'], lang='nl')
            if wdPage.get_description(lang='de') == "" or wdPage.get_description(lang='de') == "gen":
                wdPage.set_description(description="Gen der Spezies " + self.genomeInfo['name'], lang='de')
            if wdPage.get_description(lang='es') == "":
                wdPage.set_description(description="gen de la especie "+self.genomeInfo['name'], lang='es')
            if wdPage.get_description(lang='pt') == "":
                wdPage.set_description(description="gene da espécie "+self.genomeInfo['name'], lang='pt')
            if wdPage.get_description(lang='sv') == "":
                wdPage.set_description(description="Genen från arten "+self.genomeInfo['name'], lang='sv')
            if wdPage.get_description(lang='srn') == "":
                wdPage.set_description(description="gen fu a sortu "+self.genomeInfo['name'], lang='srn')

            if self.synonyms != None:
                wdPage.set_aliases(aliases=self.synonyms, lang='en', append=True)
            self.wd_json_representation = wdPage.get_wd_json_representation()
            self.wdid = wdPage.write(self.logincreds)

            PBB_Core.WDItemEngine.log('INFO', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
                    main_data_id=str(self.entrezgene),
                    exception_type='',
                    message="",
                    wd_id=self.wdid,
                    duration=time.time()-self.start
                ))


    def annotate_gene(self):
        # "Get gene annotations from mygene.info"     
        r = requests.get("http://mygene.info/v2/gene/" + str(self.entrezgene))
        return r.json()
        # return request.data
