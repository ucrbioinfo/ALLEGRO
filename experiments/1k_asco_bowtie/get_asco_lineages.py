import pandas as pd

from Bio import Entrez
import time

import urllib.request
import urllib.error
import json

guiden = 0

def get_taxonomic_id(organism_name, max_retries=5):
    """Get the taxonomic ID for a given organism name with retry logic."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={urllib.parse.quote(organism_name)}&retmode=json"
    retries = 0
    while retries < max_retries:
        try:
            with urllib.request.urlopen(url) as response:
                result = json.loads(response.read().decode())
                if result["esearchresult"]["idlist"]:
                    return result["esearchresult"]["idlist"][0]
                else:
                    return None
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            print(f"Error fetching taxonomic ID for {organism_name}: {e}")
            retries += 1
            print(f"Retrying... ({retries}/{max_retries})")
            time.sleep(2)  # Wait 2 seconds before retrying
        except Exception as e:
            print(f"Unexpected error: {e}")
            return None
    return None

def get_lineage(tax_id, max_retries=5):
    """Get the lineage for a given taxonomic ID with retry logic."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={tax_id}&retmode=xml"
    retries = 0
    while retries < max_retries:
        try:
            with urllib.request.urlopen(url) as response:
                result = Entrez.read(response)
                if "LineageEx" in result[0]:
                    lineage = result[0]["LineageEx"]
                    lineage_full = "cellular organisms"
                    for taxon in lineage:
                        lineage_full += f"; {taxon['ScientificName']}"
                    return lineage_full
                else:
                    return "Lineage not found"
        except (urllib.error.HTTPError, urllib.error.URLError) as e:
            print(f"Error fetching lineage for tax ID {tax_id}: {e}")
            retries += 1
            print(f"Retrying... ({retries}/{max_retries})")
            time.sleep(2)  # Wait 2 seconds before retrying
        except Exception as e:
            print(f"Unexpected error: {e}")
            return None
    return None

df = pd.read_csv('ascomycota_input_species.csv')


taxa = {
"Ustilaginomycotina",
"Pucciniomycotina",
"Agaricomycetes",
"Tremellomycetes",
"Dacrymycetes",
"Sordariomycetes",
"Eurotiomycetes",
"Dothideomycetes",
"Leotiomycetes",
"Pezizomycetes",
"Lecanoromycetes",
"Orbiliomycetes",
"Coniocybomycetes",
"Uncertain Placement",
"Candelariomycetes",
"Lichinomycetes",
"Geoglossomycetes",
"Xylonomycetes",
"Chytridiomycetes",
"Neocallimastigomycetes",
"Monoblepharidomycetes",
"Saccharomycotina",
"Taphrinomycotina",
"Mucoromycotina",
"Mortierellomycotina",
"Glomeromycotina",
"Uncertain Placement",
"Kickxellomycotina",
"Entomophthoromycotina",
"Zoopagomycotina",
"Cryptomycota",
"Olpidiomycota",
"Unclassified",
"Blastocladiomycota"
}

taxa_dict = dict()

for t in taxa:
    taxa_dict[t] = []

# Set your email here
Entrez.email = "lolol@gmail.com"

n = 0
for t in df['species_name'].tolist():
    # Get the taxonomic ID
    tax_id = get_taxonomic_id(t)
    if tax_id:
        lineage = get_lineage(tax_id)
        # print(f"Lineage for {t}: {lineage}")
    else:
        print(f"Taxonomic ID not found for {t}")

    s = set()
    for i in lineage.split(';'):
        s.add(i.strip())

    for i in s.intersection(taxa):
        taxa_dict[i].append(i)
        break

    if len(s.intersection(taxa)) > 1 or len(s.intersection(taxa)) == 0:
        print(s.intersection(taxa))
        
        with open(f"1000_asco_anomaly.txt", "a") as file:
            for i in s:
                file.write(i + '\n')

    n += 1
    print(n)
    time.sleep(0.5)

with open(f"1000_asco.txt", "w") as file:
    for key, value in taxa_dict.items():
        string = ''
        for v in value:
            string += v + ', '
        file.write(f"{key}: {string}\n")