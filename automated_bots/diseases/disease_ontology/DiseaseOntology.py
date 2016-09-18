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

import time
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../ProteinBoxBot_Core")
import ProteinBoxBot_Core.PBB_login as PBB_login
import ProteinBoxBot_Core.PBB_settings as PBB_settings
import ProteinBoxBot_Core.PBB_Core as PBB_Core
import ProteinBoxBot_Core.PBB_Debug as PBB_Debug
import xml.etree.cElementTree as ET
import sys
import requests
import DiseaseOntology_settings
try:
    import simplejson as json
except ImportError as e:
    import json
import copy
import traceback
from time import gmtime, strftime
import os
import pprint
from SPARQLWrapper import SPARQLWrapper, JSON


class diseaseOntology():
    def __init__(self):
        self.start = time.time()
        self.content = ET.fromstring(self.download_disease_ontology())


        self.logincreds = PBB_login.WDLogin(PBB_settings.getWikiDataUser(), os.environ['wikidataApi'])

        # Get all WikiData entries that contain a WikiData ID
        print("Getting all terms with a Disease Ontology ID in WikiData")
        doWikiData_id = dict()

        #
        sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
        sparql.setQuery("""
        PREFIX wd: <http://www.wikidata.org/entity/>
        PREFIX wdt: <http://www.wikidata.org/prop/direct/>
        select DISTINCT * where {
            {?item wdt:P699 ?doid .}
        }
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        for result in results["results"]["bindings"]:
            try:
                diseaseItem = result["item"]["value"].replace("http://www.wikidata.org/entity/", "")
                doidValue = result["doid"]["value"]
                doWikiData_id[str(doidValue)]=diseaseItem # diseaseItem[2] = DO identifier, diseaseItem[0] = WD identifier


            except Exception as e:
                print(traceback.format_exc())



        print("Getting latest version of Disease Ontology from Github")
        r = requests.get("https://api.github.com/repos/DiseaseOntology/HumanDiseaseOntology/git/refs")
        test = r.json()
        sha = test[0]["object"]["sha"]
        githubReferenceUrl = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/"+sha+"/src/ontology/doid.owl"




        for doClass in self.content.findall('.//owl:Class', DiseaseOntology_settings.getDoNameSpaces()):
          try:
            res = doClass.findall('.//oboInOwl:id', DiseaseOntology_settings.getDoNameSpaces())
            if len(res)==0:
                continue

            disVars = []
            disVars.append(doClass)
            disVars.append(githubReferenceUrl)
            disVars.append(doWikiData_id)
            disVars.append(self.logincreds)
            disVars.append(self.start)
            disVars.append(res[0].text)

            
            diseaseClass = disease(disVars)          
            

          except Exception as e:
              PBB_Core.WDItemEngine.log('ERROR', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
                        main_data_id=diseaseClass.do_id,
                        exception_type=type(e),
                        message=e.__str__(),
                        wd_id='-',
                        duration=time.time() - self.start
                    ))
              f = open('/tmp/Diseaseexceptions.txt', 'a')
              # f.write("Unexpected error:", sys.exc_info()[0]+'\n')
              f.write(diseaseClass.do_id+"\n")
              #f.write(diseaseClass.wd_json_representation)
              traceback.print_exc(file=f)
              f.close()

    def download_disease_ontology(self):
        """
        Downloads the latest version of disease ontology from the URL specified in DiseaseOntology_settings
        """
        r = requests.get(DiseaseOntology_settings.getdoUrl())
        return r.text

class  disease(object):
    def __init__(self, object):
        """
        constructor
        :param wd_do_content: Wikidata item id
        :param do_id: Identifier of the disease in Disease Ontology
        :param label: Primary label of the disease in Disease Ontology
        :param synonyms: All synonyms for the disease captured in the Disease Ontology
        :param xrefs: a dictionary with all external references of the Disease captured in the Disease Ontology
        """
        # Reference section
        doVersionURL = object[1]
        doClass = object[0]         
        self.logincreds = object[3]
        self.wd_doMappings = object[2]
        self.start = object[4]
        
        self.wd_do_content = doClass
        self.do_id = object[5]

        self.name = self.getDoValue(self.wd_do_content, './/rdfs:label')[0].text
        classDescription = self.getDoValue(self.wd_do_content, './/oboInOwl:hasDefinition/oboInOwl:Definition/rdfs:label')
        if len(classDescription)>0:
            self.description = classDescription[0].text

        if self.do_id in object[2].keys():
            self.wdid = "Q"+str(object[2][self.do_id])
        else:
            self.wdid = None
        if len(self.getDoValue(self.wd_do_content, './/owl:deprecated')) > 0 and self.getDoValue(self.wd_do_content, './/owl:deprecated')[0].text == "true":
            self.rank = "deprecated"
        else:
            self.rank = "normal"

            
        self.synonyms = []
        for synonym in self.getDoValue(self.wd_do_content, './/oboInOwl:hasExactSynonym'):
            self.synonyms.append(synonym.text)
        
        self.subclasses = []
        for subclass in self.getDoValue(self.wd_do_content, './/rdfs:subClassOf'):
            parts = subclass.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource').split("DOID_")
            if len(parts)>1:
                self.subclasses.append("DOID:"+parts[1])
            if "DOID:4" in self.subclasses:
                self.subclasses.remove("DOID:4")
        
        self.xrefs = dict()
        for xref in self.getDoValue(self.wd_do_content, './/oboInOwl:hasDbXref'):
            if not xref.text.split(":")[0] in self.xrefs.keys():
                self.xrefs[xref.text.split(":")[0]] = []
            self.xrefs[xref.text.split(":")[0]].append(xref.text.split(":")[1])



        refStatedIn = PBB_Core.WDUrl(value=doVersionURL, prop_nr='P1065', is_reference=True)
        refStatedIn.overwrite_references = True
        refImported = PBB_Core.WDItemID(value=5282129, prop_nr='P248', is_reference=True)
        refImported.overwrite_references = True
        refDoId = PBB_Core.WDString(value=str(self.do_id), prop_nr='P699', is_reference=True)
        refDoId.overwrite_references = True
        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
        refRetrieved = PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
        refRetrieved.overwrite_references = True
        do_reference = [refImported, refDoId, refRetrieved, refStatedIn]
        prep = dict()
        if self.do_id != "DOID:4":
            prep["P279"] = [PBB_Core.WDItemID(value='Q12136', prop_nr='P279', references=[copy.deepcopy(do_reference)], rank=self.rank)]
        # Subclass of disease
            for subclass in self.subclasses:
                if subclass in self.wd_doMappings.keys():
                    prep["P279"].append(PBB_Core.WDItemID(value=self.wd_doMappings[subclass], prop_nr='P279', references=[copy.deepcopy(do_reference)], rank=self.rank))


        if "Orphanet" in self.xrefs.keys():
            prep["P1550"] = []
            if isinstance(self.xrefs["Orphanet"], list):
                for id in self.xrefs["Orphanet"]:
                    prep["P1550"].append(PBB_Core.WDString(value=self.xrefs["Orphanet"], prop_nr='P1550', references=[copy.deepcopy(do_reference)], rank=self.rank))
            else:
                prep["P1550"] = [PBB_Core.WDString(value=self.xrefs["Orphanet"], prop_nr='P1550', references=[copy.deepcopy(do_reference)], rank=self.rank)]

        #disease Ontology

        prep["P699"] = [PBB_Core.WDString(value=self.do_id, prop_nr='P699', references=[do_reference], rank=self.rank)]

        if "url" in self.xrefs.keys():
            if isinstance(self.xrefs["url"], list):
                for i in self.xrefs["url"]:
                    if "//en.wikipedia.org/wiki/" in i:
                        wikilink = self.i.replace("//en.wikipedia.org/wiki/", "").replace("_", "")
                    else:
                        wikilink = None
            else:
                if "//en.wikipedia.org/wiki/" in self.xrefs["url"]:
                    wikilink = self.xrefs["url"].replace("//en.wikipedia.org/wiki/", "").replace("_", "")
                else:
                    wikilink = None
        else:
            wikilink = None

        if "ICD10CM" in self.xrefs.keys():
            prep["P494"] = []
            if isinstance(self.xrefs["ICD10CM"], list):
                for id in self.xrefs["ICD10CM"]:
                    prep["P494"].append(PBB_Core.WDString(value=id, prop_nr='P494', references=[copy.deepcopy(do_reference)], rank=self.rank))
            else:
                prep["P494"] = [PBB_Core.WDString(value=self.xrefs["ICD10CM"], prop_nr='P494', references=[copy.deepcopy(do_reference)], rank=self.rank)]

        if "ICD9CM" in self.xrefs.keys():
            prep["P493"] = []
            if isinstance(self.xrefs["ICD9CM"], list):
                for id in self.xrefs["ICD9CM"]:
                    prep["P493"].append(PBB_Core.WDString(value=id, prop_nr='P493', references=[copy.deepcopy(do_reference)], rank=self.rank))
            else:
                prep["P493"] = [PBB_Core.WDString(value=self.xrefs["ICD9CM"], prop_nr='P493', references=[copy.deepcopy(do_reference)], rank=self.rank)]

        if "MSH" in self.xrefs.keys():
            prep["P486"] = []
            if isinstance(self.xrefs["MSH"], list):
                for id in self.xrefs["MSH"]:
                    prep["P486"].append(PBB_Core.WDString(value=id, prop_nr='P486', references=[copy.deepcopy(do_reference)], rank=self.rank))
            else:
                prep["P486"] = [PBB_Core.WDString(value=self.xrefs["MSH"], prop_nr='P486', references=[copy.deepcopy(do_reference)], rank=self.rank)]

        if "NCI" in self.xrefs.keys():
            prep["P1748"] = []
            if isinstance(self.xrefs["NCI"], list):
                for id in self.xrefs["NCI"]:
                    prep["P1748"].append(PBB_Core.WDString(value=id, prop_nr='P1748', references=[copy.deepcopy(do_reference)], rank=self.rank))
            else:
                prep["P1748"] = [PBB_Core.WDString(value=self.xrefs["NCI"], prop_nr='P1748', references=[copy.deepcopy(do_reference)], rank=self.rank)]

        if "OMIM" in self.xrefs.keys():
            prep["P492"] = []
            if isinstance(self.xrefs["OMIM"], list):
                for id in self.xrefs["OMIM"]:
                    prep["P492"].append(PBB_Core.WDString(value=id, prop_nr='P492', references=[copy.deepcopy(do_reference)], rank=self.rank))
            else:
                prep["P492"] = [PBB_Core.WDString(value=self.xrefs["OMIM"], prop_nr='P492', references=[copy.deepcopy(do_reference)], rank=self.rank)]

        # add exact Matc
        refMiriamStatedIn = PBB_Core.WDItemID(value="Q16335166", prop_nr='P248', is_reference=True)
        refMiriamStatedIn.overwrite_references = True
        refMiriamUrl = PBB_Core.WDUrl(value="http://www.ebi.ac.uk/miriam/main/collections/MIR:00000233", prop_nr='P854', is_reference=True)
        refMiriamUrl.overwrite_references = True
        miriam_reference = [refMiriamStatedIn, refMiriamUrl]
        prep["P2888"] = [PBB_Core.WDUrl(value="http://purl.obolibrary.org/obo/"+self.do_id.replace(":", "_"), prop_nr='P2888', references=[copy.deepcopy(do_reference)]),
                         PBB_Core.WDUrl(value="http://identifiers.org/doid/"+self.do_id, prop_nr='P2888', references=[copy.deepcopy(miriam_reference)])]



        refStatedIn = PBB_Core.WDUrl(value="http://www.ebi.ac.uk/miriam/main/collections/MIR:00000233", prop_nr='P854', is_reference=True)
        refStatedIn.overwrite_references = True

        #io_reference = [refStatedIn]
        #prep["P1709"].append(PBB_Core.WDUrl(value="http://identifiers.org/doid/"+self.do_id, prop_nr='P1709', references=[copy.deepcopy(io_reference)], rank=self.rank))

        data2add = []
        for key in prep.keys():
            for statement in prep[key]:
                data2add.append(statement)
                # print(statement.prop_nr, statement.value)

        #if self.wdid is not None:
        #    wdPage = PBB_Core.WDItemEngine(self.wdid, item_name=self.name, data=data2add, server="www.wikidata.org", domain="diseases",append_value=['P279'])
        #else:
        wdPage = PBB_Core.WDItemEngine(item_name=self.name, data=data2add, server="www.wikidata.org", domain="diseases", append_value=['P279'])

        # wdPage.set_description(description='Human disease', lang='en')
        if wikilink is not None:
            wdPage.set_sitelink(site="enwiki", title = wikilink)
        if self.synonyms is not None:
             wdPage.set_aliases(aliases=self.synonyms, lang='en', append=True)
        self.wd_json_representation = wdPage.get_wd_json_representation()
        # PBB_Debug.prettyPrint(self.wd_json_representation)

        if self.rank == "normal":
            wdPage.write(self.logincreds)
            PBB_Core.WDItemEngine.log('INFO', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
                        main_data_id=self.do_id,
                        exception_type='',
                        message="succesfully added",
                        wd_id=self.wdid,
                        duration=time.time()-self.start
                    ))


        
    def getDoValue(self, doClass, doNode):
        returnValue = doClass.findall(doNode, DiseaseOntology_settings.getDoNameSpaces())
        return returnValue
        
     
                