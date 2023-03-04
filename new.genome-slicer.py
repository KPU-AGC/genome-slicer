#!/usr/bin/env python3
__description__ =\
"""
Purpose: Extract a region from a set of genomes using BLAST.
"""
__author__ = "Michael Ke; Erick Samera"
__version__ = "1.5.0"
__comments__ = "stable enough"
# --------------------------------------------------
from argparse import (
    Namespace,
    ArgumentParser,
    RawTextHelpFormatter)
from pathlib import Path
# --------------------------------------------------
import multiprocessing
import subprocess
import logging
import pandas as pd
import math
# --------------------------------------------------
class BlastHandler: 
    def __init__(self, blastdb: str, task: str) -> None:
        """
        """
        self.blastdb: str = blastdb
        self.blastdb_len: int = self._get_blastdb_len()
        self.task: str = task
        self.num_pool: int = multiprocessing.cpu_count()
        return None
    def _get_blastdb_len(self):
        """
        Function returns the length of the blastdb.
        
        ### Parameters
            None

        ### Returns
            blastdb_len: int
                the length of the blast database
        """
        with open(self.blastdb) as blastdb_file: return blastdb_file.read().count('>')
    def blast(self, query: Path) -> list:
        """
        Function performs BLAST using subprocess.

        ### Parameters:
            query: Path
                path of the query fasta file.
            output_path: Path
                path of directory for output
            tag: str
                string for tagging the blast-output json? 
        
        ### Returns:
            None
        """
        args = [
            "blastn",
            "-task", str(self.task),
            "-db", str(self.blastdb),
            "-num_alignments", str(self.blastdb_len),
            "-outfmt", "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp",
            "-query", str(query),
        ]
        raw_blast_results = subprocess.run(args, capture_output=True)
        blast_results: list = []
        try: blast_results: list = [result for result in raw_blast_results.stdout.decode().split('\n')]
        except: pass
        return blast_results
    def multi_blast(self, queries):
        """
        """
        def job_allocator(queries, NUM_GROUPS):
            list_size = math.floor(len(queries)/NUM_GROUPS)
            remainder = len(queries)%NUM_GROUPS
            job_list = []
            for group in range(NUM_GROUPS-1): 
                job_list.append(queries[0+(list_size*group):list_size+(list_size*group)])
            job_list.append(queries[(list_size*(NUM_GROUPS-1)):(list_size*NUM_GROUPS+remainder)])
            return job_list
        
        #Allocate the jobs
        job_list = job_allocator(queries, self.num_pool)
        
        #Run the BLAST
        pool = multiprocessing.Pool(self.num_pool)
        results = pool.map(self.blast, job_list)
        print(results)
        return results
# --------------------------------------------------
def _process_blast_data(_blast_results: list, _use_names: bool, _tag: str, _output_path: Path) -> dict:
    """
    """
    processed_results = {}

    for result in _blast_results:
        if not result: continue
        # just every single key possible, for future-proofing

        keys = ("qseqid","qgi","qacc","qaccver","qlen","sseqid","sallseqid","sgi","sallgi","sacc","saccver","sallacc","slen","qstart","qend","sstart","send","qseq","sseq","evalue","bitscore","score","length","pident","nident","mismatch","positive","gapopen","gaps","ppos","frames","qframe","sframe","btop","staxids","sscinames","scomnames","sblastnames","sskingdoms","stitle","salltitles","sstrand","qcovs","qcovhsp")
        result_dict = {key: value for key, value  in zip(keys, result.split('\t'))}

        dict_key = result_dict['qseqid']
        if dict_key not in processed_results:
            processed_results[dict_key] = []


        # if use names, raise a warning if the key is empty
        if _use_names:
            use_name_keys = ('sacc', 'staxids', 'sscinames')
            for key in use_name_keys:
                if result_dict[key] in ('N/A', '', '0'): logging.warning(f"{result_dict['saccver']} Missing data: {key}. Consider remaking blastdb w/ updated taxdb")

        processed_results[dict_key].append(result_dict)
    
    for i, (key, blast_results) in enumerate(processed_results.items()):
        prefix = '_'.join([_tag, str(i)])
        pd.DataFrame(blast_results).to_csv(_output_path.joinpath(f'{prefix}_blast-output.csv'))
    return processed_results
def _output_fasta(_processed_blast_data: dict, _output_path: Path, use_names: bool):
    """
    Function generates FASTA output from BLAST data.

    ### Parameters:
        blast_data: list
            list of data extracted from blast result
        filename: str
            filename FASTA file output
        output_path: Path
            output path of FASTA file
        use_names: bool
            boolean value for whether names are used in the header.

    ### Returns:
        None
    """
    for key, blast_data in _processed_blast_data.items():
        fasta_path = _output_path.joinpath(f'{key}_blast-output.fasta')
        with open(fasta_path, mode='w') as fasta_file:
            for blast_result in blast_data:
                id: str = ''
                description: str = ''

                if use_names:
                    sci_name = blast_result['sscinames'].replace(' ', '-')
                    accession = blast_result['sacc']
                    taxid = blast_result['staxids']
                    
                    id: str = f'{sci_name}'
                    description: str = f'_{taxid}_{accession}'
                else:
                    seq_id: str = blast_result['saccver']
                    id: str = f'{seq_id}'

                ungap_sequence = blast_result['sseq'].replace('-', '')
                fasta_file.write(f'>{id}{description}\n')
                fasta_file.write(f'{ungap_sequence}\n')
# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description=__description__,
        epilog=f"v{__version__} : {__author__} | {__comments__}",
        formatter_class=RawTextHelpFormatter)

    parser.add_argument('query_path',
        action='store',
        type=Path,
        help='path to query')
    parser.add_argument('--db', dest='db',
        action='store',
        type=Path,
        required=True,
        help='BLAST database path')
    parser.add_argument('--task', dest='task',
        action='store',
        type=str,
        default='megablast',
        help="BLAST task (default='megablast')")
    parser.add_argument('-o', '--output', dest='output_path',
        action='store',
        type=Path,
        default=None,
        help="path to output (default=in-place)",)
    parser.add_argument('--tag', dest='tag',
        action='store',
        type=str,
        default='',)   
    parser.add_argument('--use-name', dest='name',
        action='store_true',
        help="",)

    args = parser.parse_args()

    # parser errors and processing
    # --------------------------------------------------
    if not args.output_path: args.output_path = args.query_path.parent
    args.output_path.mkdir(parents=True, exist_ok=True)

    return args
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    logging.basicConfig(
        level=logging.INFO,
        handlers=[
            logging.FileHandler(args.output_path.joinpath('probesearch.log')),
            logging.StreamHandler()],
        format='%(asctime)s:%(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S',)

    # initiate blasthandler then perform BLAST
    blastHandler = BlastHandler(args.db, args.task)
    queries = open(args.query_path).read().count('>')

    blast_results = blastHandler.blast(args.query_path) if queries < 2 else blastHandler.multi_blast(args.query_path) 
    
    # process the blast results into dictionaries
    processed_results = _process_blast_data(blast_results, args.name, args.tag, args.output_path)

    # output the multi-FASTA of BLAST output
    _output_fasta(processed_results, args.output_path, args.name)

    return None
# --------------------------------------------------
if __name__ == '__main__':
    main()