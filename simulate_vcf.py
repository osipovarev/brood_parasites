#!/usr/bin/env python3
#
# This script generates a polarized VCF of N diploid individuals 
import argparse
import msprime
import sys
import numpy as np
import random


__author__ = "Ekaterina Osipova, 2024."


def simulate_vcf(args):
    # Simulate the population
    tree_sequence = msprime.simulate(
        sample_size=2 * args.num_individuals,
        Ne=args.population_size,
        length=args.chromosome_length,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate
    )
    return(tree_sequence)


def simulate_bottleneck_vcf(args):
    # Define the demographic events to simulate a bottleneck
    demographic_events = [

        # At the end of the bottleneck period, return to original size
        msprime.PopulationParametersChange(
            time=args.bottleneck_end, initial_size=args.population_size),

        # At time of this # of generations ago, population size reduced requested size
        msprime.PopulationParametersChange(
            time=args.bottleneck_start, initial_size=args.bottleneck_size)
    ]

    # Simulate the tree sequence with the bottleneck
    tree_sequence = msprime.simulate(
        sample_size=2 * args.num_individuals, 
        Ne=args.population_size,
        length=args.chromosome_length,
        mutation_rate=args.mutation_rate,
        recombination_rate=args.recombination_rate,
        demographic_events=demographic_events
    )
    return(tree_sequence)



# def simulate_bottleneck_and_outgroup(args):
#     # Define the demographic events to simulate a bottleneck
#     demographic_events = [
#         msprime.PopulationParametersChange(
#             time=args.bottleneck_end, initial_size=args.population_size),
#         msprime.PopulationParametersChange(
#             time=args.bottleneck_start, initial_size=args.bottleneck_size)
#     ]

#     # An additional sample for the outgroup
#     num_individuals_with_outgroup = 2 * args.num_individuals + 2  # 20 individuals + 1 outgroup individual

#     # Simulate the tree sequence with the bottleneck and one outgroup individual
#     tree_sequence = msprime.simulate(
#         sample_size=num_individuals_with_outgroup,  # 42 haplotypes: 40 from main population, 2 from outgroup
#         Ne=args.population_size,
#         length=args.chromosome_length,
#         mutation_rate=args.mutation_rate,
#         recombination_rate=args.recombination_rate,
#         demographic_events=demographic_events
#     )

#     # Assume the outgroup is the last two haplotypes
#     outgroup_nodes = [num_individuals_with_outgroup - 2, num_individuals_with_outgroup - 1]

#     # Polarize the VCF based on the outgroup
#     tree_sequence.write_vcf(sys.stdout,
#                             ploidy=2, 
#                             ancestral_alleles={
#                                 variant.site.id: variant.alleles[
#                                     variant.genotypes[outgroup_nodes[0]]
#                                 ]
#                                 for variant in tree_sequence.variants()
#                             },
#                             individual_mapping=range(2 * args.num_individuals)  # Map only the intended individuals, exclude outgroup
#     )


def write_vcf(tree_sequence, args):
    # Write the VCF with or without ATGC genotypes in REF and ALT fileds

    print("##fileformat=VCFv4.2")
    print(f"##contig=<ID={args.chrom},length={args.chromosome_length}>")
    print(f"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    print(f"##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">")
    chrom_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
    print(chrom_line + '\t'.join(["ingroup_" + str(i) for i in range(args.num_individuals)]))

    # Randomly assign nucleotide bases A, T, G, C to alleles
    allele_bases = ['A', 'T', 'G', 'C']

    for variant in tree_sequence.variants():
        pos = int(variant.site.position)
        
        if args.atgc:
            # Randomize assignment of nucleotide bases to alleles, ensuring they are different
            alleles = random.sample(allele_bases, len(variant.alleles))
            ref_allele = alleles[0]
            alt_alleles = alleles[1:]
        else:
            ref_allele = '0'
            alt_alleles = ['1']

        # Write VCF entry
        genotypes = [
            f"{variant.genotypes[2 * i]}|{variant.genotypes[2 * i + 1]}"
            for i in range(args.num_individuals)
        ]
        print(f"{args.chrom}\t{pos + 1}\t.\t{ref_allele}\t{','.join(alt_alleles)}\t.\t.\t.\tGT\t" + "\t".join(genotypes))


def main():
    parser = argparse.ArgumentParser(description="Simulate genetic data and write to VCF")
    
    ## Add srguments to define population parameters
    parser.add_argument("--population-size", type=int, default=100_000, help="Population size; default: 100_000")
    parser.add_argument("--num-individuals", type=int, default=20, help="Number of diploid individuals; default: 20")
    parser.add_argument("--mutation-rate", type=float, default=1e-8, help="Mutation rate; default: 1e-8")
    parser.add_argument("--recombination-rate", type=float, default=1e-8, help="Recombination rate; default: 1e-8")
    
    ## Add arguments to define bottleneck; if you want to simulate expansion, change bottleneck_size to a large number
    parser.add_argument("--bottleneck", action='store_true', help="specify if you want to simulate recent bottleneck")
    parser.add_argument("--bottleneck_start", type=int, default=100_000, help="when bottleneck started, in generations; default: 100_000")
    parser.add_argument("--bottleneck_size", type=int, default=1000, help="population size at bottleneck; default: 1000; \
                         if you want to simulate recent population expansion, just change bottleneck_size to a large number!")
    parser.add_argument("--bottleneck_end", type=int, default=50_000, help="when bottleneck ended, in generations; default: 100_000")
    
    ## Add arguments to define length and name of generated chromosome
    parser.add_argument("--chrom", type=str, default="chr1", help="chromosome name; default: chr1")
    parser.add_argument("--chromosome-length", type=int, default=2_000_000, help="Chromosome length; default: 2_000_000")

    ## Add argument to generate real genotypes: ATGC 
    parser.add_argument("--atgc", action='store_true', help="specify if you want to generate real genotypes (A,T,G,C not just 0 or 1) in REF and ATL fileds")

    ## Add randomness
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed for simulation; default: 42")
    # parser.add_argument("--output-vcf", type=str, default="simulated.vcf", help="Output VCF file name; default: simulated.vcf")


    args = parser.parse_args()
    if args.bottleneck:
        tree_sequence = simulate_bottleneck_vcf(args)
    else:
        tree_sequence = simulate_vcf(args)
    write_vcf(tree_sequence, args)


if __name__ == "__main__":
    main()
