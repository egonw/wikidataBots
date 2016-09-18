__author__ = 'andra'

from SPARQLWrapper import SPARQLWrapper, JSON
import pprint

sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
sparql.setQuery("""
PREFIX wikibase: <http://wikiba.se/ontology#>
PREFIX wd: <http://www.wikidata.org/entity/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX p: <http://www.wikidata.org/prop/>
PREFIX v: <http://www.wikidata.org/prop/statement/>
SELECT * WHERE {
  ?gene wdt:P279 wd:Q7187 .
  ?gene wdt:P279 wd:Q8054 .
 }
""")
sparql.setReturnFormat(JSON)
results = sparql.query().convert()
if len(results["results"]["bindings"]) > 0:
    raise ValueError('There is a Wikidata item of subclass Protein and Gene')
print(len(results["results"]["bindings"]))



