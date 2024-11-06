import pandas as pd
from io import StringIO
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


guides = list()
with open('1k_ascomycota_library.txt', 'r') as f:
    for l in f.readlines():
        guides.append(l.strip())

with open('reads/reads.txt', 'w') as f:
        for g in guides:
            for n in ['AGG', 'CGG', 'TGG', 'GGG']:
                f.write(f'>{g}{n}\n{g}{n}\n')

df = pd.read_csv('ascomycota_input_species.csv')
    
delimited_cds_list = df['ortho_file_name'].tolist()

for file_path in delimited_cds_list:
    base_path = os.getcwd()
    name = file_path[:file_path.find('_cds')]
    output_align_temp_path = 'reads/reads.txt'
    index_basename = os.path.join('indices/'+ name + '_idx')

    bowtie_build_command = ['bowtie-build', '-f', base_path + '/../../../ALLEGRO_Fungi_Downloader/data/fourdbs_concat/cds/' + file_path, base_path + '/' + index_basename]

    process = subprocess.Popen(bowtie_build_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, stderr = process.communicate()

    if stderr:
        print(f'{stderr.decode()}')
        print('Exiting.')

def process_file(file_path):
    base_path = os.getcwd()
    name = file_path[:file_path.find('_cds')]
    index_basename = os.path.join('indices', name + '_idx')

    bowtie_build_command = [
        'bowtie-build', 
        '-f', 
        os.path.join(base_path, '../../../ALLEGRO_Fungi_Downloader/data/fourdbs_concat/cds/', file_path), 
        os.path.join(base_path, index_basename)
    ]

    process = subprocess.Popen(bowtie_build_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, stderr = process.communicate()
    
    if process.returncode != 0:
        print(f"Error processing {file_path}: {stderr.decode('utf-8')}")
    else:
        print(f"Successfully processed {file_path}")

# Use ThreadPoolExecutor to process files in parallel
with ThreadPoolExecutor(max_workers=64) as executor:
    futures = {executor.submit(process_file, file_path): file_path for file_path in delimited_cds_list}
    
    for future in as_completed(futures):
        file_path = futures[future]
        try:
            future.result()
        except Exception as exc:
            print(f"Error processing {file_path}: {exc}")

for file_path in delimited_cds_list:
    base_path = os.getcwd()
    name = file_path[:file_path.find('_cds')]
    output_align_temp_path = 'reads/reads.txt'
    index_basename = os.path.join('indices/'+ name + '_idx')

    bowtie_command = ['bowtie', '-a', '-v', '0', '--suppress', '5,6,7,8', '-f', base_path + '/' + index_basename, base_path + '/' + output_align_temp_path, '--threads', '64']
    process = subprocess.Popen(bowtie_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    if stderr:
        print(f'{stderr.decode()}')
        print('Exiting.')

    df2 = pd.read_csv(StringIO(stdout.decode()), sep='\t',
                names=['query_name', 'strand', 'reference_name', 'start_position'],
                dtype={'reference_name': str, 'query_name': str})
    
    df2.to_csv('out/' + file_path)

# def process_file(file_path):
#     # if os.path.exists(os.path.join('out', file_path)):
#     #     return

#     base_path = os.getcwd()
#     name = file_path[:file_path.find('_cds')]
#     output_align_temp_path = 'reads/reads.txt'
#     index_basename = os.path.join('indices', name + '_idx')

#     bowtie_command = [
#         'bowtie', '-a', '-v', '0', '--suppress', '5,6,7,8', '-f', 
#         os.path.join(base_path, index_basename), 
#         os.path.join(base_path, output_align_temp_path)
#     ]
#     process = subprocess.Popen(bowtie_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#     stdout, stderr = process.communicate()

#     if stderr:
#         print(f'{stderr.decode()}')
#         print('Exiting.')
#         return

#     df2 = pd.read_csv(StringIO(stdout.decode()), sep='\t',
#                       names=['query_name', 'strand', 'reference_name', 'start_position'],
#                       dtype={'reference_name': str, 'query_name': str})
    
#     df2.to_csv(os.path.join('out', file_path))
#     print(f'Done with {file_path}')

# # Use ThreadPoolExecutor to process files in parallel
# with ThreadPoolExecutor(max_workers=64) as executor:
#     futures = {executor.submit(process_file, file_path): file_path for file_path in delimited_cds_list}
    
#     for future in as_completed(futures):
#         file_path = futures[future]
#         try:
#             future.result()
#         except Exception as exc:
#             print(f"Error processing {file_path}: {exc}")