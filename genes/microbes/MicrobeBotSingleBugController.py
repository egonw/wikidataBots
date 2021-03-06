import MicrobeBotResources as MBR
import MicrobeBotGenes as MBG
import MicrobeBotProteins as MBP
import MicrobeBotEncoder as MBE
import sys
import os
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")
import PBB_login
import datetime

__author__ = 'timputman'


if len(sys.argv) < 6:
    print("   You did not supply the proper arguments!")
    print("   Usage: MicrobeBotModularPackage.py <Wikidata user name> <Wikidata Password> <domain "
          "i.e. genes/proteins/encode_genes/encode_proteins>, <run_type>  i.e full, quick")
    sys.exit()
else:
    pass


# Login to Wikidata with bot credentials
login = PBB_login.WDLogin(sys.argv[1], sys.argv[2])

# Retrieve Current Bacterial Reference Genomes from NCBI
print('Retrieving current list of NCBI Bacterial Reference Genomes')
print('Standby...')


# if reference genome information should be updated
run_type = sys.argv[5]
genome_records = MBR.get_ref_microbe_taxids(run_type)

tid = sys.argv[3]
spec_strain = genome_records.loc[genome_records['taxid'] == tid]



# Retrieve gene and protein records from UniProt and Mygene.info by taxid
print('Retrieving gene records for {} taxid:{}'.format(spec_strain.iloc[0]['organism_name'], tid))
gene_records = MBR.mgi_qg_resources(tid)  # PANDAS DataFrame

# Iterate through gene_records for reading and writing to Wikidata
print('Commencing {} bot run  for {}'.format(sys.argv[4], spec_strain.iloc[0]['organism_name']))
gene_count = 0

for record in gene_records:
    print('{}/{}'.format(gene_count, len(gene_records)))  # , spec_strain['organism_name'])
    if sys.argv[4] == 'genes':
        gene = MBG.wd_item_construction(record, spec_strain, login)
        if gene == 'success':
            gene_count += 1
    if sys.argv[4] == 'proteins':
        protein = MBP.wd_item_construction(record, spec_strain, login)
        if protein == 'success':
            gene_count += 1
    if sys.argv[4] == 'encoder':
        encoder = MBE.encodes(record, login)
        if encoder == 'success':
            gene_count += 1
