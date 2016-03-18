import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")
import PBB_Core
import pprint
import MicrobeBotWDFunctions as wdo
import time

__author__ = 'timputman'


def wd_item_construction(gene_record, login):
    """
    generate pbb_core item object based on resources pandas dataframe for each gene
    :param gene_record: pandas dataframe of combined UniProt NCBI and MyGene.info data
    :return: PBB_Core object of WD item with claims and references for Genes
    """
    gene_record = gene_record
    # Gets WD Item Id for microbe gene is found in
    strain_qid = wdo.WDSparqlQueries(string=gene_record['taxid'], prop='P685').wd_prop2qid()
    # gets WD Item label for parent microbe based on qid
    strain_label = wdo.WDSparqlQueries(qid=strain_qid).wd_qid2label()
    item_description = 'microbial gene found in {}'.format(strain_label)

    def gene_item_statements():
        """
        construct list of referenced statements to past to PBB_Core Item engine
        :return:
        """
        # creates reference object for WD gene item claim
        ncbi_gene_reference = wdo.reference_store(source='ncbi_gene', identifier=gene_record['_id'])

        # claims for datatype string.
        WD_String_CLAIMS = {'P351': str(gene_record['_id']),
                            'P2393': gene_record['locus_tag'],
                            'P644': str(gene_record['genomic_pos.start']),
                            'P645': str(gene_record['genomic_pos.end']),
                            }
        # claims for datytpe item
        WD_Item_CLAIMS = {'P703': strain_qid,
                          'P279': 'Q7187',
                          }
        # convert integer representation of strand to corresponding WD item (Forward Strand/Reverse Strand)
        if gene_record['genomic_pos.strand'] == '1':
            WD_Item_CLAIMS['P2548'] = 'Q22809680'
        else:
            WD_Item_CLAIMS['P2548'] = 'Q22809711'

        statements = []
        # process to pbb_Core data value object and append to statments for each valid item in each datatype dict
        # WDItemID datatype
        for k, v in WD_Item_CLAIMS.items():
            if v != 'nan':
                statements.append(PBB_Core.WDItemID(value=v, prop_nr=k,
                                                    references=[ncbi_gene_reference]))
        # Entrez gene id
        for k, v in WD_String_CLAIMS.items():
            if v != 'nan':
                statements.append(PBB_Core.WDString(value=v, prop_nr=k,
                                                    references=[ncbi_gene_reference]))
        return statements

    item_name = gene_record['name'] + "\t" + gene_record['locus_tag']

    # attempt to instantiate PBB_Core item object by finding the proper item in wikidata or creating a new one (Json)
    start = time.time()
    try:
        wd_item_gene = PBB_Core.WDItemEngine(item_name=item_name, domain='genes', data=gene_item_statements(),
                                             use_sparql=True)
        #pprint.pprint(wd_item_gene.get_wd_json_representation())
        wd_item_gene.set_label(item_name)
        wd_item_gene.set_description(item_description, lang='en')
        wd_item_gene.set_aliases([gene_record['symbol'], gene_record['locus_tag']])
        wd_item_gene.write(login=login)
        new_mgs = ''
        # log actions to log file
        if wd_item_gene.create_new_item:
            new_mgs = ': New item'
        PBB_Core.WDItemEngine.log('INFO', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
            main_data_id=gene_record['_id'],
            exception_type='',
            message='success{}'.format(new_mgs),
            wd_id=wd_item_gene.wd_item_id,
            duration=time.time() - start
        ))
        print('success')

    except Exception as e:
        print(e)
        PBB_Core.WDItemEngine.log('ERROR', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
            main_data_id=gene_record['_id'],
            exception_type=type(e),
            message=e.__str__(),
            wd_id='',
            duration=time.time() - start
        ))

    end = time.time()
    print('Time elapsed:', end - start)
