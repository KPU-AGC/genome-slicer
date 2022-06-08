import subprocess
import pandas as pd
from math import floor
import multiprocessing
import argparse
import pathlib
import io
import json
import pprint

class BlastHandler: 
    def __init__(self, blastdb, task):
        self.blastdb = blastdb
        self.blastdb_len = self._get_blastdb_len()
        self.task = task
        self.NUM_POOL = multiprocessing.cpu_count()
    def _get_blastdb_len(self): 
        blastdb_file = open(self.blastdb)
        data = blastdb_file.read()
        return data.count('>')
    def blast(self, query, output_path, tag): 
        #Generate the oligo temporary file
        json_path = output_path.joinpath(f'{tag}_blast-output.json')
        args = [
            "blastn",
            "-task",
            str(self.task),
            "-db",
            str(self.blastdb),
            "-num_alignments",
            str(self.blastdb_len),
            "-outfmt",
            "13",
            "-query",
            str(query),
            "-out",
            json_path
        ]
        subprocess.run(args)
    def multi_blast(self, queries): 
        def job_allocator(queries, NUM_GROUPS):
            list_size = floor(len(queries)/NUM_GROUPS)
            remainder = len(queries)%NUM_GROUPS
            job_list = []
            for group in range(NUM_GROUPS-1): 
                job_list.append(queries[0+(list_size*group):list_size+(list_size*group)])
            job_list.append(queries[(list_size*(NUM_GROUPS-1)):(list_size*NUM_GROUPS+remainder)])
            return job_list
        #Allocate the jobs
        job_list = job_allocator(queries, self.NUM_POOL)
        #Run the BLAST
        pool = multiprocessing.Pool(self.NUM_POOL)
        results = pool.map(self.blast, job_list)
        return results

def parse_args(): 
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'query_path',
        action='store',
        type=pathlib.Path,
        help='Path to query'
    )
    parser.add_argument(
        '--db',
        dest='db',
        action='store',
        type=pathlib.Path,
        help='BLAST database path'
    )
    parser.add_argument(
        '--task',
        dest='task',
        action='store',
        type=str,
        default='megablast',
        help='BLAST task'
    )
    parser.add_argument(
        '--output',
        dest='output_path',
        action='store',
        type=pathlib.Path,
        default=None,
        help='',
    )
    parser.add_argument(
        '--tag',
        dest='tag',
        action='store',
        type=str,
        default='',
    )
    args = parser.parse_args()
    if not args.output_path:
        args.output_path = args.query_path.parent
    return args.query_path, args.db, args.task, args.output_path, args.tag

def json_output(data, output_path): 
    # write to screened JSON file
    prettied_json =  pprint.pformat(data, compact=True, sort_dicts=False, width=160).replace("'", '"')
    with open(output_path.joinpath(f'output_json.json'), 'w') as new_JSON:
        new_JSON.write(prettied_json)

def get_json_output_list(query_path, result_path, tag): 
    with open(query_path, 'r') as query_file:
        num_seq = query_file.read().count('>')
    list_json_path = []
    for index in range(num_seq):
        json_path = result_path.joinpath(f'{tag}_blast-output_{index+1}.json')
        list_json_path.append(json_path)
    return list_json_path

def get_processed_blast_data(blast_data):
    processed_data = []
    for blast_hit in blast_data['BlastOutput2']['report']['results']['search']['hits']:
        #Main places to get the data from
        blast_description = blast_hit['description'][0]
        blast_hsp = blast_hit['hsps'][0]
        #All of the important fields  
        accession = blast_description['accession']
        sequence = blast_hsp['hseq']
        taxid = blast_description['taxid']
        sciname = blast_description['sciname']
        processed_data.append(
            {
                'accession':accession,
                'sequence':sequence,
                'taxid':taxid,
                'sciname':sciname,
            }
        )
    return processed_data

def output_fasta(blast_data, file_name, output_path): 
    fasta_path = output_path.joinpath(f'{file_name}.fasta')
    with open(fasta_path, 'r') as fasta_file:
        for sequence in blast_data:
            accession = sequence['accession']
            taxid = sequence['taxid']
            sci_name = sequence['sciname']
            fasta_file.write(f'>{accession}_{taxid}_{sci_name}\n')
            fasta_file.write(sequence['sequence']+'\n')

def main(): 
    query_path, db_path, task, output_path, tag = parse_args()
    blastHandler = BlastHandler(db_path, task)
    #Run the BLAST jobs
    blastHandler.blast(query_path, output_path, tag)
    #open the main json --> this shows you where the 
    #json_file = open(json_path, 'r')
    #json_data = json.load(json_file)
    #Note that the output files are going to be
    #tag_blast-output_#.json
    output_dict = {}
    list_json_output = get_json_output_list(query_path, output_path, tag)
    for json_path in list_json_output:
        with open(json_path, 'r') as json_file:
            blast_data = json.load(json_file)
            processed_data = get_processed_blast_data(blast_data)
        output_dict[json_path.stem] = processed_data
    #output_json_path = output_path.joinpath('output_json.json')
    json_output(output_dict, output_path)
    for key in output_dict: 
        output_fasta(output_dict[key], key, output_path)

if __name__ == '__main__': 
    main() 

