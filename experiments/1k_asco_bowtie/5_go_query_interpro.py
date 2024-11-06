import requests
from Bio import SeqIO
import time
import os
import json

# Define the InterPro API endpoint
interpro_url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"

job_to_info = dict()

count = 0
for f in os.listdir("target_prot_sequences/"):
    records = SeqIO.parse(os.path.join('target_prot_sequences', f), 'fasta')

    if os.path.exists(f'interpro_results/{f}'):
        with open(f'interpro_results/{f}', 'r') as json_file:
            existing_data = json.load(json_file)
            
            for r in records:
                if r.id in existing_data:
                    continue

    for r in records:
        # Define the protein sequence you want to query
        protein_sequence = str(r.seq).replace('*', '')

        # Define the parameters for the API request
        params = {
            'email': 'nitothecat@gmail.com',
            'goterms': 'true',
            'stype': 'p',
            'sequence': protein_sequence,
        }

        for attempt in range(3):
            response = requests.post(interpro_url, data=params)
            if response.status_code == 200:
                job_id = response.text
                print(f"Job submitted successfully. Job ID: {job_id}")
                job_to_info[job_id] = (f, r.id)
                break
            else:
                print(f"Error: Unable to submit job to InterPro API. Status code: {response.status_code}. Attempt {attempt + 1} of {3}")
                time.sleep(10)  # Wait for 10 seconds before retrying

            if attempt == 2:
                with open('failed.txt', 'a') as failf:
                    failf.write(f'{f}, {r.id}\n')

                raise Exception("Maximum retry attempts reached. Unable to submit job.")

        time.sleep(2)

        count += 1

        if count > 10:
            for job_id, job_info in job_to_info.items():
                # Define the endpoint to check job status
                status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"

                # Poll the API for job status
                while True:
                    status_response = requests.get(status_url)
                    if status_response.status_code == 200:
                        status = status_response.text
                        if status == 'FINISHED':
                            print("Job completed successfully.")
                            break
                        elif status in ['RUNNING', 'PENDING']:
                            print(f"Job status: {status}. Checking again in 30 seconds...")
                            time.sleep(30)
                        else:
                            print(f"Job status: {status}. Exiting.")
                            break
                    else:
                        print(f"Error: Unable to check job status. Status code: {status_response.status_code}")
                        break

                # Define the endpoint to retrieve job results
                result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"

                # Retrieve the job results
                result_response = requests.get(result_url)

                # Check if the request was successful
                if result_response.status_code == 200:
                    results = result_response.json()

                    # Load existing data from the JSON file if it exists
                    if os.path.exists(f'interpro_results/{job_info[0]}'):
                        with open(f'interpro_results/{job_info[0]}', 'r') as json_file:
                            existing_data = json.load(json_file)
                    else:
                        existing_data = {}

                    # Add the new data under a new key
                    new_key = job_info[1]
                    existing_data[new_key] = results

                    # Write the updated data back to the JSON file
                    with open(f'interpro_results/{job_info[0]}', 'w') as json_file:
                        json.dump(existing_data, json_file, indent=4)
                    
                count -= 1
            job_to_info = dict()

time.sleep(2)

for job_id, job_info in job_to_info.items():
    # Define the endpoint to check job status
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"

    # Poll the API for job status
    while True:
        status_response = requests.get(status_url)
        if status_response.status_code == 200:
            status = status_response.text
            if status == 'FINISHED':
                print("Job completed successfully.")
                break
            elif status in ['RUNNING', 'PENDING']:
                print(f"Job status: {status}. Checking again in 30 seconds...")
                time.sleep(30)
            else:
                print(f"Job status: {status}. Exiting.")
                break
        else:
            print(f"Error: Unable to check job status. Status code: {status_response.status_code}")
            break

    # Define the endpoint to retrieve job results
    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/json"

    # Retrieve the job results
    result_response = requests.get(result_url)

    # Check if the request was successful
    if result_response.status_code == 200:
        results = result_response.json()

        # Load existing data from the JSON file if it exists
        if os.path.exists(f'interpro_results/{job_info[0]}'):
            with open(f'interpro_results/{job_info[0]}', 'r') as json_file:
                existing_data = json.load(json_file)
        else:
            existing_data = {}

        # Add the new data under a new key
        new_key = job_info[1]
        existing_data[new_key] = results

        # Write the updated data back to the JSON file
        with open(f'interpro_results/{job_info[0]}', 'w') as json_file:
            json.dump(existing_data, json_file, indent=4)
        
    count -= 1
