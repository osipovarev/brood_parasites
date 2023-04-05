#!/usr/bin/env python3
#
# This script parses degeneracy-all.sites.bed file 
# and outputs bed file with positions that are clear $d sites in all corresponding transripts

import argparse
from collections import defaultdict


__author__ = "Ekaterina Osipova, 2023."


def main():
    ## Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed_degen', type=str, help='degeneracy bed file to parse: 4-bed format with site degeneracy in the 5th column')
    args = parser.parse_args()

    ## Read input bed file line by line; store sites with their degeneracy in a dict
    sites_degen_dict = defaultdict(list)
    with open(args.bed_degen, 'r') as inf:
        for line in inf.readlines(): 
            site = tuple(line.split()[:3])
            degeneracy = line.split()[4]
            sites_degen_dict[site].append(degeneracy)
    
    ## Go through the dictionary and output sites that are clear 4D in every overlapping transcript
    for site in sites_degen_dict:
        degeneracy = sites_degen_dict[site]
        if (len(set(degeneracy)) == 1) and (set(degeneracy) == {'4'}):
            print('{}'.format('\t'.join(list(site))))


if __name__ == "__main__":
    main()

