#!/usr/bin/env python3

import pandas as pd


def check_deletions_in_bed(deletions, bed):
    """
    Check if deletions are in the reference bed file

    Parameters
    ----------
    deletions : pd.DataFrame
        Deletions table, should contain the columns 'CHROM', 'POS', 'END' and 'interval'
    bed : pd.DataFrame
        Bed file, should contain the columns 'chrom', 'start', 'end', 'locus_tag', 'gene', 'drugs' and 'interval'
    
    Returns
    -------
    pd.DataFrame
        Annotated deletions table, containing the columns 'chrom', 'start', 'end', 'locus_tag', 'gene' and 'drugs'
    
    Notes
    -----
    Compares genomic ranges stored as pd.Interval objects, with closed bounds

    """
    deletions_annotated = []
    for index, row in deletions.iterrows():
        for index2, row2 in bed.iterrows():
            if row['CHROM'] == row2['chrom']:
                if row['interval'].overlaps(row2['interval']):
                    dictionary_overlap = {
                        'chrom': row['CHROM'],
                        'start': row['POS'],
                        'end': row['END'],
                        'locus_tag': row2['locus_tag'],
                        'gene': row2['gene'],
                        'drugs': row2['drugs']
                    }
                    deletions_annotated.append(dictionary_overlap)
    
    return pd.DataFrame(deletions_annotated)                


def main(args):
    # Read the deletion table
    deletions = pd.read_csv(args.input, sep='\t', header=0)

    if deletions.shape[0] == 0:
        with open(args.output, 'w') as f:
            f.write('chrom\tstart\tend\tlocus_tag\tgene\tdrugs\n')

    else:
        deletions['interval'] = deletions.apply(lambda x: pd.Interval(x['POS'], x['END'], closed='both'), axis=1)

        # Read the reference bed file
        bed = pd.read_csv(args.bed, sep='\t', header=None, names=['chrom', 'start', 'end', 'locus_tag', 'gene', 'drugs'])
        bed['interval'] = bed.apply(lambda x: pd.Interval(x['start'], x['end'], closed='both'), axis=1)

        annotated_deletions = check_deletions_in_bed(deletions, bed)

        # Write the annotatetd deletion table
        annotated_deletions.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Filter deletions that are not in the reference genome')
    parser.add_argument('--input', help='Deletion table', required=True)
    parser.add_argument('--bed', help='Filtered deletion table', required=True)
    parser.add_argument('--output', help='Output file', default='deletions_annotated.tsv')
    args = parser.parse_args()

    main(args)
