import os
import re
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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

guides = list()
with open('1k_ascomycota/1k_ascomycota_library.txt', 'r') as f:
    for l in f.readlines():
        guides.append(l.strip())


guides_targets = dict()
with open('1k_ascomycota/1k_ascomycota.txt', 'r') as f:
    for l in f.readlines():
        if l.startswith('Using'):
            # Define the regex pattern to extract the guides and their lists
            pattern = r"\('(\w+)', 1.0, \[([^\]]+)\]\)"

            # Find all matches in the string
            matches = re.findall(pattern, l)

            # Create the dictionary
            guides_targets = {match[0]: [item.strip().strip("'") for item in match[1].split(',')] for match in matches}


df = pd.read_csv('ascomycota_input_species.csv')

for guide in guides:
    s = [guide + f for f in ['AGG', 'CGG', 'GGG', 'TGG']]

    targets = list()
    locs = list()

    files = df[df['species_name'].isin(guides_targets[guide])]['ortho_file_name'].tolist()

    for f in files:
        already_counted = False

        with open(f'out/{f}', 'r') as inf:
            lines = inf.readlines()

            for l in lines:
                if l.startswith(','):
                    continue

                seq = l.split(',')[1]
                
                if seq in s:
                    locs.append(l.split(',')[3])
                    locs.append(f)
                    locs.append(l.split(',')[1])
                    locs.append('----')

                    if not already_counted:
                        # targets.append(l.split(',')[3])
                        targets.append(f.split('_')[0] + ' ' + f.split('_')[1])
                        already_counted = True


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
        taxa_dict[t] = 0

    # Set your email here
    Entrez.email = "lolol@gmail.com"

    n = 0
    for t in targets[n:]:
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
            taxa_dict[i] += 1
            break

        if len(s.intersection(taxa)) > 1 or len(s.intersection(taxa)) == 0:
            print(s.intersection(taxa))
            
            with open(f"guide{guiden}_anomaly.txt", "a") as file:
                for i in s:
                    file.write(i + '\n')

        n += 1
        print(n)
        time.sleep(0.5)

    with open(f"guide{guiden}.txt", "w") as file:
        for key, value in taxa_dict.items():
            file.write(f"{key}: {value}\n")
    guiden += 1