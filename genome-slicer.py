import subprocess
import pandas as pd
from math import floor
import multiprocessing
import argparse
import pathlib
import io

class Blast: 
    def __init__(self, blastdb, task):
        self.blastdb = blastdb
        self.blastdb_len = self._get_blastdb_len()
        self.task = task
        self.NUM_POOL = multiprocessing.cpu_count()
    def _get_blastdb_len(self): 
        blastdb_file = open(self.blastdb)
        data = blastdb_file.read()
        return data.count('>')
    def blast(self, query, output_path): 
        #Generate the oligo temporary file
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
            output_path.joinpath('output.json')
        ]
        subprocess.run(args)

    def multi_blast(self, oligos): 
        def job_allocator(oligos, NUM_GROUPS):
            list_size = floor(len(oligos)/NUM_GROUPS)
            remainder = len(oligos)%NUM_GROUPS
            job_list = []
            for group in range(NUM_GROUPS-1): 
                job_list.append(oligos[0+(list_size*group):list_size+(list_size*group)])
            job_list.append(oligos[(list_size*(NUM_GROUPS-1)):(list_size*NUM_GROUPS+remainder)])
            return job_list
        #Allocate the jobs
        job_list = job_allocator(oligos, self.NUM_POOL)
        #Run the BLAST
        pool = multiprocessing.Pool(self.NUM_POOL)
        results = pool.map(self.blast, job_list)
        return results

def parse_args(): 
    parser = argparse.ArgumentParser('Extract the sequence from the BLAST results')
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
    args = parse_args()
    if not args.output_path:
        args.output_path = args.query_path.parent
    return args.query_path, args.db, args.task, args.output_path

def main(): 
    query_path, db_path, task, output_path = parse_args()
    blastHandler = Blast(db_path, task)
    data = blastHandler.blast(query_path)
    print(data)

if __name__ == '__main__': 
    main() 

