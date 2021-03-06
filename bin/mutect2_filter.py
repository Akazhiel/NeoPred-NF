#!/usr/bin/env python

import argparse
import sys

def parse_args(args=None):
    Description = "Filter and reformat the vcf generated by Mutect2."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT> <TUMOR_ID> <NORMAL_ID>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input vcf file.")
    parser.add_argument("FILE_OUT", help="Output filtered vcf file.")
    parser.add_argument("TUMOR_ID", help="ID of Tumor sample to be replaced.")
    parser.add_argument("NORMAL_ID", help="ID of Normal sample to be replaced.")
    return parser.parse_args(args)

def index_column_substring(your_list, substring):
    for i, s in enumerate(your_list):
        if substring in s:
            return i
    return -1


def mutect2_filter(input, output, sample1_ID, sample2_ID):
    """
    Filters a Mutect2 VCF to only keep variants that PASS
    Tumor and Normal IDs are also replaced by the words TUMOR and NORMAL
    :param input: the input VCF generated with Mutect2
    :param output: the name of filtered VCF file
    :param sample1_ID: ID that identifies the tumor sample
    :param sample2_ID: ID that identifies the normal sample
    :return: None
    """
    filtered_vcf = open(output, 'w')
    vcf = open(input)
    for line in vcf:
        if line.startswith('#') and not line.startswith('#CHROM'):
            if sample1_ID in line:
                filtered_vcf.write(line.replace(sample1_ID, "TUMOR"))
            elif sample2_ID in line:
                filtered_vcf.write(line.replace(sample2_ID, "NORMAL"))
            else:
                filtered_vcf.write(line)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line.replace(sample1_ID, "TUMOR").replace(sample2_ID, "NORMAL"))
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if 'PASS' in columns[headers.index('FILTER')]:
                filtered_vcf.write(line)
        else:
            filtered_vcf.write(line)
    filtered_vcf.close()
    vcf.close()

def main(args=None):
    args = parse_args(args)
    mutect2_filter(args.FILE_IN, args.FILE_OUT, args.TUMOR_ID, args.NORMAL_ID)

if __name__ == "__main__":
    sys.exit(main())
