#!/usr/bin/env python3
#
# This script parses an annotated with snpEff VCF file and extracts allele frequencies 
# separating them by class: synonymous / nonsynonymous


import argparse
import pysam

__author__ == "Ekaterina Osipova, 2023."


vcf_data = pysam.VariantFile(dir_path + vcf_name)

AF_nonsyn = []
AF_syn = []
for i in vcf_data.fetch():
    POS = i.pos
    REF = i.ref
    ALT = i.alts[0]
    AC = i.info['AC'][0]
    AN = i.info['AN']
    AA = i.info['AA']
    
    if AA == ALT:
        AC = AN - AC
    elif (AA != ALT) and (AA != REF):
        AC = 0
    AN = AN - 2
    AF = AC / AN
    
    if (AF != 1) and (AF != 0):
        if 'missense_variant' in i.info['ANN'][0]:
            AF_nonsyn.append(AF)
        elif 'synonymous_variant' in i.info['ANN'][0]:
            AF_syn.append(AF)

print(AF_syn)
