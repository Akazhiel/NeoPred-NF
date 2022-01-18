#! /usr/bin/env python
"""
This tools uses the output of merge_results.py and HLAs predicted with either
the dna_pipeline.py and/or the rna_pipeline.py to compute neoantigens affinity
binding score.
The tool extracts the WT and MUT peptides and makes affinity binding
predictions for the HLAs (class I). 
The tools uses MHC-flurry for the predictions.
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from collections import Counter
from collections import defaultdict
import json
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys


def main(hla, overlap_final, alleles_file, mode, cutoff):

    # TODO perform sanity check on input parameters
    
    HLA_dict = defaultdict(list)
    for f in hla:
        with open(f) as f:
            next(f)
            for line in f.readlines():
                columns = line.split('\t')[1:7]
                for i in range(len(columns)):
                    key = columns[i].split("*")[0]
                    HLA_dict[key].append(columns[i])

    # Filter HLAs by occurrences
    filtered_hla = []
    for alleles in HLA_dict.values():
        f = Counter(alleles)
        items = f.most_common()
        filtered_hla += [y[0] for y in items if int(y[1]) >= int(cutoff)]
    # MHCflurry format (HLA-A*02:01)
    filtered_hla = ['HLA-{}'.format(':'.join(x.split(':')[0:2])) for x in filtered_hla]

    # Filter HLAs by allowed alleles in MHCflurry
    allowed_alleles = set()
    with open(alleles_file) as f:
        for x in f.readlines():
            allowed_alleles.add(x.strip())
    filtered_hla = [x for x in filtered_hla if x in allowed_alleles]

    with open('allowed_alleles.txt', 'w') as f:
        f.write(','.join(filtered_hla))

    if len(filtered_hla) == 0:
        sys.stderr.write("Error, list of HLAs is empty.\n")
        sys.exit(1)

    # Create protein FASTA file
    with open('protein_sequences_mu.fasta', 'w') as fwrite_mu:
        with open('protein_sequences_wt.fasta', 'w') as fwrite_wt:
            with open(overlap_final) as fread:
                lines = fread.readlines()
                header = lines.pop(0).strip().split('\t')
                added_proteins_mu = set()
                added_proteins_wt = set()
                for line in lines:
                    columns = line.strip().split('\t')
                    pass_somatic = int(columns[header.index('Number of DNA samples (passing)')]) > 0
                    pass_germline = int(columns[header.index('Number of RNA samples (passing)')]) > 0
                    if (mode == "both" and pass_somatic and pass_germline) or (mode == "dna" and pass_somatic)\
                        or (mode == "rna" and pass_germline) or (mode == "either" and (pass_somatic or pass_germline)):
                        protein_name = '{}_{}_{}'.format(''.join(columns[header.index('Variant key')].split()),
                                                         columns[header.index('cDNA change')],
                                                         columns[header.index('AA change')])
                        protein_seq_mu = columns[header.index('Mut Epitope')].strip()
                        protein_seq_wt = columns[header.index('Wt Epitope')].strip()
                        # Should probably make sure that all the letters in the protein seq are alpha (isalpha())
                        if protein_seq_mu != '-' and '.' not in protein_seq_mu and protein_seq_mu not in added_proteins_mu:
                            fwrite_mu.write('>{}\n{}\n'.format(protein_name, protein_seq_mu))
                            added_proteins_mu.add(protein_seq_mu)
                        # Should probably make sure that all the letters in the protein seq are alpha (isalpha())
                        if protein_seq_wt != '-' and '.' not in protein_seq_wt and protein_seq_wt not in added_proteins_wt:
                            fwrite_wt.write('>{}\n{}\n'.format(protein_name, protein_seq_wt))
                            added_proteins_wt.add(protein_seq_wt)
                del added_proteins_mu
                del added_proteins_wt

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--hla', nargs='+', default=None, required=True,
                        help='A file or files containing predicted HLAs from DNA/RNA (OptiType TSV format)')
    parser.add_argument('--variants', default=None, required=True,
                        help='A file with the final variants generated with merge_results.py (table format)')
    parser.add_argument('--alleles', default=None, required=True,
                        help='A file containing the allowed alleles in MHCflurry')
    parser.add_argument('--mode', default='either',
                        help='Mode to use to extract peptides from the variants (both, dna (default), rna, either)',
                        choices=['both', 'dna', 'rna', 'either'])
    parser.add_argument('--cutoff', default=1,
                        help='Cutoff to filter in how many samples each HLA has to be present to be kept (default = 1)')
    args = parser.parse_args()
    main(args.hla, args.variants, args.alleles, args.mode, args.cutoff)