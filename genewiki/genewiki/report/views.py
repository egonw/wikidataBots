from django.shortcuts import render

from genewiki.report.queries import gene_counts


def counts(request):

    results_mammals = gene_counts('mammals')
    results_microbes = gene_counts('microbes')

    vals = {}
    # test data
    # chart_data = '''
    #      [new Date(2016, 06, 10), 60058, 73450, 42079,36962],
    #      [new Date(2016, 06, 13), 70000, 80000, 50000,40000],
    #      [new Date(2016, 06, 15), 70000, 80000, 50000,40000]
    # '''
    # Overview Tab
    protein_cdata = '''
          [new Date(2016, 06, 10), 60058],
          [new Date(2016, 06, 13), 70000],
          [new Date(2016, 06, 15), 74000]
    '''
    compound_cdata = '''
          [new Date(2016, 06, 10), 60058],
          [new Date(2016, 06, 13), 70000],
          [new Date(2016, 06, 15), 74000]
    '''
    ontology_cdata = '''
          [new Date(2016, 06, 10), 60058],
          [new Date(2016, 06, 13), 70000],
          [new Date(2016, 06, 15), 74000]
    '''

    table_data = '''
          ['Genes', 'Human', 100, 10000, 5000, 'SELECT'],
          ['Proteins', 'Human', 600, 25000, 50, 'SELECT'],
          ['DiseaseOntology', 'Human', 1000, 8000, 6000, 'SELECT'],
          ['Compounds', 'Human', 150,  8000, 10, 'SELECT']
    '''
    # Jenkins Error Tab
    jenkins_change_bar = '''
        ['Event', '8-16-2016', '8-24-2016'],
        ['insertions', 500, 200],
        ['deletions', 100, 30],
        ['updates', 2000, 10000]
    ''' 
    jenkins_pie = '''
        ['Error', 'Types per Run'],
        ['Script',     11],
        ['Update',      2],
        ['Mistyped',  2],
        ['Deletion', 2]
    ''' 
    jenkins_table = '''
        ['Insertion Error',  'Major', 'Value Not Found', 'INSERT ()--error', '8-10-2016'],
        ['Update Error',  'Minor', 'Mismatch datatype', 'Update ()--error', '8-10-2016'],
        ['Insertion Error',  'Minor', 'duplicate value',  'INSERT ()--error', '8-16-2016']
    ''' 
    # Curation Issues Tab
    curation_table = '''
        ['Mismerge',  'minor', 'manual curation required', 'protein gene merged', '8-10-2016', 'http://wikidata.fixme'],
        ['Misformed', 'minor', 'manual curation required', 'statements with incorrect identifiers', '8-10-2016', 'http://wikidata.fixme'],
        ['Mismerge',  'minor', 'manual curation required', 'protein gene merged', '8-12-2016', 'http://wikidata.fixme'],
    ''' 
    curation_pie = '''
        ['Error', 'Types per Run'],
        ['Script',     11],
        ['Update',      2],
        ['Mistyped',  2],
        ['Deletion', 2]
    '''
    # vals['mammal_cdata'], vals['mammal_clegend'] = line_chart(results_mammals)
    vals['mammal_cdata'], vals['mammal_clegend'] = stackedbar_chart(results_mammals)
    # vals['microbe_cdata'], vals['microbe_clegend'] = line_chart(results_microbes)
    vals['microbe_cdata'], vals['microbe_clegend'] = stackedbar_chart(results_microbes)
    # vals['cdata'] = chart_data #for testing
    vals['protein_cdata'] = protein_cdata
    vals['compound_cdata'] = compound_cdata
    vals['ontology_cdata'] = ontology_cdata
    vals['tdata'] = table_data
    vals['jbar'] = jenkins_change_bar
    vals['jpie'] = jenkins_pie
    vals['jtable'] = jenkins_table
    vals['ctable'] = curation_table
    vals['cpie'] = curation_pie 

    return render(request, 'report/stats.jade', vals)


def line_chart(results):
    curdate = ''
    chart_data = ''
    legend = ''  # creates key
    unique_species = {}  # finds unique species
    for date in sorted(results.keys()):
        # add formatted data to string for google charts
        for species in sorted(results[date].keys()):  # need to sort to be sure generating same data lines for the chart
            ccount = str(results[date][species])
            species_name = species.replace(" ", "_")
            # get the month and subtract 1 charts index months from 0
            cdate = str(date).split("-")
            cdate_mm = "0" + str(int(cdate[1]) - 1)
            cdate_yyyy = cdate[0]
            day_time = cdate[2].split(" ")
            cdate_dd = day_time[0]
            if(curdate != date):
                chart_data += '[new Date(%s, %s, %s), %s, ' % (cdate_yyyy, cdate_mm, cdate_dd, ccount)
            else:
                chart_data += '%s, ' % (ccount)
            if species not in unique_species:
                legend += '''data.addColumn('number', '%s');\n''' % (species_name)
            curdate = date
            unique_species[str(species)] = 1
        # remove last comma and space
        chart_data = chart_data[:-2]
        chart_data += '],\n'
    return chart_data, legend

def stackedbar_chart(results):
    curdate = ''
    chart_data = ''
    legend = '''['Rank','''  # creates key
    unique_species = {}  # finds unique species
    for date in sorted(results.keys()):
        # add formatted data to string for google charts
        for species in sorted(results[date].keys()):  # need to sort to be sure generating same data lines for the chart
            ccount = str(results[date][species])
            species_name = species.replace(" ", "_")
            # get the month and subtract 1 charts index months from 0
            cdate = str(date).split("-")
            cdate_mm = "0" + str(int(cdate[1]) - 1)
            cdate_yyyy = cdate[0]
            day_time = cdate[2].split(" ")
            cdate_dd = day_time[0]
            if(curdate != date):
                chart_data += '[new Date(%s, %s, %s), %s, ' % (cdate_yyyy, cdate_mm, cdate_dd, ccount)
            else:
                chart_data += '%s, ' % (ccount)
            if species not in unique_species:
                legend += ''' '%s',''' % (species_name)
            curdate = date
            unique_species[str(species)] = 1
        # remove last comma and space
        chart_data = chart_data[:-2]
        chart_data += '],\n'
    legend += '],\n'
    return chart_data, legend
