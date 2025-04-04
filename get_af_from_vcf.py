#!/usr/bin/env python3
#
# This script parses a polarized annotated with snpEff VCF file and
# extracts allele frequencies separating them by class:
# synonymous / nonsynonymous


import argparse
import pysam


__author__ = "Ekaterina Osipova, 2023."


def extract_vcf_entires(vcf_data, region):
    ## Parses vcf pysam object and extracts snps from the requested region;
    ## if no region provided, goes through the entire vcf

    if region == '':
        vcf_entires = vcf_data.fetch()
    else:
        contig = region.split('|')[0]
        start = int(region.split('|')[1])
        end = int(region.split('|')[2])
        vcf_entires = vcf_data.fetch(contig, start, end)
    return vcf_entires


def get_af_by_class(vcf_entires):
    ## Extracts syn and nonsyn allele frequencies

    AF_nonsyn = []
    AF_syn = []    
    for i in vcf_entires:
        POS = i.pos
        REF = i.ref
        ALT = i.alts[0]
        AC = i.info['AC'][0]
        AN = i.info['AN']
        AA = i.info['AA']
        
        if AA != REF:
            # check if snp wasn't introduced with the outgroup
            if (AA == ALT) and (AC != 2):
                # flip AC
                AC = AN - AC
            # exlude if it is new or 3-way to resolve
            else:
                AC = 0
        # remove outgroup from AN
        AN = AN - 2
        AF = AC / AN

        if (AF != 0) and (AF != 1):
            if 'missense_variant' in i.info['ANN'][0]:
                AF_nonsyn.append(AF)
            elif 'synonymous_variant' in i.info['ANN'][0]:
                AF_syn.append(AF)
    return AF_nonsyn, AF_syn


def main():
    ## Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument(
                        '-v',
                        '--vcf', 
                        type=str, 
                        help='input VCF file'
                        )
    parser.add_argument(
                        '-r', 
                        '--region', 
                        type=str, 
                        default='', 
                        help='bed line with a region to extract from VCF: chr_start_end'
                        )

    args = parser.parse_args()
    
    vcf_data = pysam.VariantFile(args.vcf)
    region = args.region

    ## Read input VCF and extract AF by class
    vcf_entires = extract_vcf_entires(vcf_data, region)
    AF_nonsyn, AF_syn = get_af_by_class(vcf_entires)

    ## Output AF
    AF_nonsyn_line = ','.join([str(round(i, 3)) for i in AF_nonsyn])
    AF_syn_line = ','.join([str(round(i, 3)) for i in AF_syn])
    print('__pN__:{}\t__pS__:{}'.format(AF_nonsyn_line, AF_syn_line))


if __name__ == '__main__':
    main()
