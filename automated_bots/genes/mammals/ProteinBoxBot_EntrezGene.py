#!usr/bin/env python
# -*- coding: utf-8 -*-

'''
Author:Andra Waagmeester (andra@waagmeester.net)

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
# Load the path to the PBB_Core library
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../../ProteinBoxBot_Core")
import PBB_Core
import PBB_Debug
import PBB_Functions
import PBB_login
import PBB_settings

# Resource specific 
import gene

import traceback
from datetime import date, datetime, timedelta


try:

    speciesInfo = dict()
    speciesInfo["homo_sapiens"] = dict()
    speciesInfo["9606"] = dict()
    speciesInfo["mus_musculus"] = dict()
    speciesInfo["10090"] = dict()
    speciesInfo["rattus_norvegicus"] = dict()
    speciesInfo["10116"] = dict()
    speciesInfo["macaca_nemestrina"] = dict()
    speciesInfo["9545"] = dict()
    
    speciesInfo["homo_sapiens"]["taxid"] = "9606"
    speciesInfo["homo_sapiens"]["wdid"] = "Q15978631"
    speciesInfo["homo_sapiens"]["main_key"] = "homo_sapiens"
    speciesInfo["homo_sapiens"]["name"] = "Homo sapiens"
    speciesInfo["homo_sapiens"]["release"] = "Q20950174"
    speciesInfo["homo_sapiens"]["genome_assembly"] = "Q20966585"
    speciesInfo["homo_sapiens"]["genome_assembly_previous"] = "Q21067546"

    speciesInfo["9606"]["taxid"] = "9606"
    speciesInfo["9606"]["wdid"] = "Q15978631"
    speciesInfo["9606"]["name"] = "Homo sapiens"
    speciesInfo["9606"]["main_key"] = "homo_sapiens"
    speciesInfo["9606"]["release"] = "Q20950174"
    speciesInfo["9606"]["genome_assembly"] = "Q20966585"
    speciesInfo["9606"]["genome_assembly_previous"] = "Q21067546"

    speciesInfo["mus_musculus"]["taxid"] = "10090"
    speciesInfo["mus_musculus"]["wdid"] = "Q83310"
    speciesInfo["mus_musculus"]["name"] = "Mus musculus"
    speciesInfo["mus_musculus"]["main_key"] = "mus_musculus"
    speciesInfo["mus_musculus"]["release"] = "Q20973051"
    speciesInfo["mus_musculus"]["genome_assembly"] = "Q20973075"
    speciesInfo["mus_musculus"]["genome_assembly_previous"] = "Q20973075"

    speciesInfo["10090"]["taxid"] = "10090"
    speciesInfo["10090"]["wdid"] = "Q83310"
    speciesInfo["10090"]["main_key"] = "mus_musculus"
    speciesInfo["10090"]["name"] = "Mus musculus"
    speciesInfo["10090"]["release"] = "Q20973051"
    speciesInfo["10090"]["genome_assembly"] = "Q20973075"
    speciesInfo["10090"]["genome_assembly_previous"] = "Q20973075"

    speciesInfo["rattus_norvegicus"]["taxid"] = "10116"
    speciesInfo["rattus_norvegicus"]["wdid"] = "Q184224"
    speciesInfo["rattus_norvegicus"]["main_key"] = "rattus_norvegicus"
    speciesInfo["rattus_norvegicus"]["name"] = "Rattus norvegicus"
    speciesInfo["rattus_norvegicus"]["release"] = "Q19296606"
    speciesInfo["rattus_norvegicus"]["genome_assembly"] = "Q21578759"
    speciesInfo["rattus_norvegicus"]["genome_assembly_previous"] = "None"

    speciesInfo["10116"]["taxid"] = "10116"
    speciesInfo["10116"]["wdid"] = "Q184224"
    speciesInfo["10116"]["main_key"] = "rattus_norvegicus"
    speciesInfo["10116"]["name"] = "Rattus norvegicus"
    speciesInfo["10116"]["release"] = "Q19296606"
    speciesInfo["10116"]["genome_assembly"] = "Q21578759"
    speciesInfo["10116"]["genome_assembly_previous"] = "None"

    speciesInfo["macaca_nemestrina"]["taxid"] = "9545"
    speciesInfo["macaca_nemestrina"]["wdid"] = "Q618026"
    speciesInfo["macaca_nemestrina"]["main_key"] = "macaca_nemestrina"
    speciesInfo["macaca_nemestrina"]["name"] = "Macaca nemestrina"
    speciesInfo["macaca_nemestrina"]["release"] = "Q24067511"
    speciesInfo["macaca_nemestrina"]["genome_assembly"] = "Q24067518"
    speciesInfo["macaca_nemestrina"]["genome_assembly_previous"] = "None"

    speciesInfo["9545"]["taxid"] = "9545"
    speciesInfo["9545"]["wdid"] = "Q618026"
    speciesInfo["9545"]["name"] = "Macaca nemestrina"
    speciesInfo["9545"]["main_key"] = "macaca_nemestrina"
    speciesInfo["9545"]["release"] = "Q24067511"
    speciesInfo["9545"]["genome_assembly"] = "Q24067518"
    speciesInfo["9545"]["genome_assembly_previous"] = "None"



    if len(sys.argv) == 1:
        print("Please provide one of the following species as argument: "+ str(speciesInfo.keys()))
        print("Example: python ProteinBoxBot_EntrezGene.py homo_sapiens")
        sys.exit()
    else:
        if not sys.argv[1] in speciesInfo.keys():
            print(sys.argv[1] + " is not (yet) supported.")
            sys.exit()

    tempvar = dict()
    tempvar["speciesInfo"] = speciesInfo
    tempvar["species"] = sys.argv[1]
    genome = gene.genome(tempvar)
    
      
except (Exception):
    print(traceback.format_exc())
    # client.captureException()  