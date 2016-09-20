import requests
import sys
import pandas as pd
import pprint
import json
import time
import numpy as np

import PBB_Core
import PBB_login


def get_entrez_qid_map(prop_nr):
    query = '''
        SELECT * WHERE {{
            ?qid wdt:{} ?id .
            ?qid wdt:P703 wd:Q15978631 .
        }}
        '''.format(prop_nr)

    results = PBB_Core.WDItemEngine.execute_sparql_query(query=query)['results']['bindings']

    id_wd_map = dict()
    for z in results:
        id_wd_map[z['id']['value']] = z['qid']['value'].split('/')[-1]

    return id_wd_map


def main():
    mirtar_base = pd.read_csv('./data/hsa_MTI.csv', low_memory=False,
                                dtype={'References (PMID)': np.str, 'Target Gene (Entrez Gene ID)': np.str})
    # mirtar_base = pd.read_csv('~/Downloads/hsa_MTI.csv', header=0, sep=',')
    pprint.pprint(mirtar_base.head(2))

    pprint.pprint(mirtar_base.iloc[1, 0])
    # pprint.pprint(mirtar_base.index)

    entrez_qid_map = get_entrez_qid_map('P351')

    # create subclass of 'mature microRNA'
    # create encoded by (if this can be determined, would be important)
    # create 'encodes' property on the genes
    # create 'found in taxon' property

    login_obj = PBB_login.WDLogin(user='ProteinBoxBot', pwd='sNxvAlNtjQ24')

    data = []
    let_7b_5p = mirtar_base.loc[mirtar_base['miRNA'].values == 'hsa-let-7b-5p', :]
    let_7b_5p = let_7b_5p[['miRTarBase ID', 'miRNA', 'Target Gene (Entrez Gene ID)']].drop_duplicates()
    print(let_7b_5p.count())
    print(len(let_7b_5p['miRTarBase ID'].unique()))
    for count, mir in let_7b_5p.iterrows():

        acc = mir['miRTarBase ID']
        ncbi_id = mir['Target Gene (Entrez Gene ID)']
        mirna_label = mir['miRNA']

        if ncbi_id not in entrez_qid_map:
            continue

        print(acc, ncbi_id)

        refs = [[
            PBB_Core.WDItemID(value='Q6826951', prop_nr='P248', is_reference=True),  # stated in
            PBB_Core.WDExternalID(value=acc, prop_nr='P2646', is_reference=True),  # source element
            PBB_Core.WDItemID(value='Q1860', prop_nr='P407', is_reference=True),  # language of work
            PBB_Core.WDMonolingualText(value=mirna_label, language='en', prop_nr='P1476', is_reference=True),
            PBB_Core.WDTime(time=time.strftime('+%Y-%m-%dT00:00:00Z'), prop_nr='P813', is_reference=True)  # retrieved
        ]]

        stmnt = PBB_Core.WDItemID(value=entrez_qid_map[ncbi_id], prop_nr='P128', references=refs)
        data.append(stmnt)

    wd_item = PBB_Core.WDItemEngine(wd_item_id='Q21414101', data=data)
    print(len(data))
    wd_item.write(login_obj)

    # pprint.pprint(mirtar_base.count())

    # unique_mirs = mirtar_base['miRNA'].unique()
    #
    # for mir in unique_mirs:
    #
    #     # references = generate_refs(chembl_id)
    #
    #     curr_mir_df = unique_mirs[unique_mirs['miRNA'] == mir]
    #
    #     statements = list()
    #
    #     # mature miRNA Q23838648
    #
    #     for x in curr_mir_df.index:
    #         curr_mesh = curr_mir_df.loc[x, 'mesh_id']
    #         if pd.notnull(curr_mesh) and curr_mesh in mesh_wd_map:
    #             print(chembl_id, curr_mesh, 'found')



if __name__ == '__main__':
    sys.exit(main())