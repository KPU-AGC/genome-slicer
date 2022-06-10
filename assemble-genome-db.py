import argparse
import pathlib
import csv
import json

class NCBIDataHandler():
    def __init__(self, base_path) -> None:
        #Paths
        self.base_path = self._check_base_path(base_path)
        self.catalog_path = self.base_path.joinpath('dataset_catalog.json')
        self.assembly_data_path = self.base_path.joinpath('assembly_data_report.jsonl')
        #Input data
        self.dataset_catalog = self._load_json(self.catalog_path)
        self.assembly_data = self._load_jsonl(self.assembly_data_path)
        #Processed data
    def _check_base_path(self, base_path) -> pathlib.Path:
        """
        Deal with the fact that .zip files often unzip to an a folder
        """
        if base_path.joinpath('ncbi_dataset').exists(): 
            return base_path.joinpath('ncbi_dataset/data')
        else:
            return base_path.joinpath('data')
    def _load_json(self, path) -> dict:
        '''
        json files loading
        '''
        print(path)
        with open(path, 'r') as file: 
            data = json.load(file)
        return data
    def _load_jsonl(self, path) -> list: 
        '''
        jsonl files have separate json structures on each line.
        '''
        list_data = []
        print(path)
        with open(path, 'r') as file: 
            for line in file: 
                list_data.append(json.loads(line))
        return list_data
    def generate_taxid_map(self, output_path) -> None:
        def _get_taxids(assembly_data) -> list:
            '''
            Return tuples of (accession, taxid)
            '''
            list_taxid_accession = []
            for assembly in assembly_data: 
                list_taxid_accession.append((assembly['assemblyInfo']['assemblyAccession'], assembly['taxId']))
            return list_taxid_accession
        def _get_accessions(dataset_catalog):
            """
            Return dictionary containing accessions and subaccessions
            """
            def _get_gb_accessions(jsonl_path): 
                list_gb_accessions = []
                with open(jsonl_path, 'r') as seq_report_file: 
                    for line in seq_report_file: 
                        sequence = json.loads(line)
                        list_gb_accessions.append(sequence['genbankAccession'])
                return list_gb_accessions
            accessions = {}
            for entry in dataset_catalog['assemblies']: 
                if 'accession' not in entry.keys(): 
                    pass
                else:
                    accession = entry['accession']
                    #Instead of looping through the paths, just figure it out
                    seq_report_path = self.base_path.joinpath(f'{accession}/sequence_report.jsonl')
                    accessions[accession] = _get_gb_accessions(seq_report_path)
            return accessions
        #Get necessary data from the summary files
        list_taxid_accessions = _get_taxids(self.assembly_data)
        accessions = _get_accessions(self.dataset_catalog)
        #Setup output file
        output_file_path = output_path.joinpath('taxid_map.txt')
        with open(output_file_path, 'w') as output_file: 
            #Looping
            for data in list_taxid_accessions: 
                accession = data[0]
                taxid = data[1]
                for gb_accession in accessions[accession]:
                    output_file.write(f'{gb_accession}\t{taxid}\n')
        
def parse_args(): 
    parser = argparse.ArgumentParser('Program to process NCBI datasets')
    parser.add_argument(
        'ncbi_path',
        action='store',
        type=pathlib.Path,
        help='Path to ncbi_dataset directory'
    )
    parser.add_argument(
        '--output',
        dest='output_path',
        action='store',
        type=pathlib.Path,
        default=pathlib.Path.cwd(),
        help='Path to output taxid_map.txt'
    )
    args = parser.parse_args()
    return args.ncbi_path, args.output_path

def main():
    # 
    ncbi_path, output_path = parse_args()
    ncbi_data = NCBIDataHandler(ncbi_path)
    #for child in ncbi_data.base_path.iterdir(): 
    #    print(child)
    #print(type(ncbi_data.dataset_catalog['assemblies']))
    #for entry in ncbi_data.dataset_catalog['assemblies']: 
    #   print(entry)
    #print('DATASET CATALOG')
    #print(ncbi_data.dataset_catalog)
    #print('ASSEMBLY_DATA')
    #print(ncbi_data.assembly_data)
    #for item in ncbi_data.assembly_data: 
    #   print(item.keys())
    ncbi_data.generate_taxid_map(output_path)

if __name__ == '__main__': 
    main()