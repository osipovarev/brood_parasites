#!/usr/bin/env python3
#
# This script generates a polarized VCF of N diploid individuals and an outgroup \n 
# Under scenario of recent population bottleneck
import argparse
import msprime


__author__ = "Ekaterina Osipova, 2024."


def simulate_and_write_vcf(args):
    # Setup demography with a recent bottleneck
    demography = msprime.Demography()
    demography.add_population(name="ingroup", initial_size=args.population_size)
    demography.add_population(name="outgroup", initial_size=args.population_size)
    
    # Define a bottleneck event: small size at recent time
    bottleneck_time = args.bottleneck_time  # Time when bottleneck starts
    bottleneck_size = args.bottleneck_size  # Population size during bottleneck
    
    # Set the population size change events
    demography.add_population_parameters_change(
        time=bottleneck_time, 
        initial_size=bottleneck_size, 
        population="ingroup"
    )
    
    # Recover to original size
    demography.add_population_parameters_change(
        time=bottleneck_time + args.bottleneck_duration,  # assuming duration is set
        initial_size=args.population_size, 
        population="ingroup"
    )
    
    # Add a population split
    demography.add_population_split(
        time=args.divergence_time,
        derived=["ingroup"],
        ancestral="outgroup"
    )

    # Simulate the tree sequence
    tree_sequence = msprime.sim_ancestry(
        samples={"ingroup": args.num_individuals, "outgroup": args.outgroup_num},
        ploidy=2,  
        demography=demography,
        sequence_length=args.chromosome_length,
        recombination_rate=args.recombination_rate,
        random_seed=args.random_seed
    )

    # Introduce mutations on the tree sequence
    mutated_ts = msprime.sim_mutations(
        tree_sequence,
        rate=args.mutation_rate,
        random_seed=args.random_seed
    )

    # Write the VCF
    with open(args.output_vcf, "w") as vcf_file:
        # Write the VCF header
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write(f"##contig=<ID=chr1,length={args.chromosome_length}>\n")
        vcf_file.write(f"##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        
        for i in range(args.num_individuals):
            vcf_file.write(f"\tingroup_{i}")
        vcf_file.write("\toutgroup_0\n")

        # Process and flip each variant; add ancestral allele
        for variant in mutated_ts.variants():
            original_ref = variant.alleles[0]
            original_alts = variant.alleles[1:]
            
            outgroup_allele_index = variant.genotypes[-1]  # Assuming last haplotype index
            ancestral_allele = variant.alleles[outgroup_allele_index]
            
            # Flip REF and first ALT
            if original_alts:
                new_ref = original_alts[0]
                new_alts = [original_ref] + list(original_alts[1:])
                # Write out variant line
                vcf_file.write(f"{args.chrom}\t{int(variant.site.position) + 1}\t.\t{new_ref}\t{','.join(new_alts)}\t.\t.\tAA={ancestral_allele}\tGT")
                for i in range(len(variant.genotypes) // 2):
                    # Adjust genotype based on flipped alleles
                    g1 = variant.genotypes[2 * i]
                    g2 = variant.genotypes[2 * i + 1]
                    gt1 = 0 if g1 == 1 else 1
                    gt2 = 0 if g2 == 1 else 1
                    vcf_file.write(f"\t{gt1}|{gt2}")
                vcf_file.write("\n")

    print(f"VCF written to {args.output_vcf}")


def main():
    parser = argparse.ArgumentParser(description="Simulate genetic data and write to VCF")

    ## Add srguments to define population parameters
    parser.add_argument("--population-size", type=int, default=100_000, help="Population size; default: 100_000")
    parser.add_argument("--outgroup-num", type=int, default=1, help="Number of outgroup individuals; default: 1")
    parser.add_argument("--num-individuals", type=int, default=20, help="Number of diploid individuals; default: 20")
    parser.add_argument("--mutation-rate", type=float, default=1e-8, help="Mutation rate; default: 1e-8")
    parser.add_argument("--recombination-rate", type=float, default=1e-8, help="Recombination rate; default: 1e-8")
    parser.add_argument("--divergence-time", type=int, default=10_000_000, help="Time of divergence in generations for the outgroup; default: 10_000_000")
    
    ## Add arguments to define bottleneck 
    parser.add_argument("--bottleneck_time", type=int, default=100_000, help="when bottleneck started, in generations; default: 100_000")
    parser.add_argument("--bottleneck_size", type=int, default=1000, help="population size at bottleneck; default: 1000")
    parser.add_argument("--bottleneck_duration", type=int, default=10_000, help="duration of bottleneck, in generations; default: 10_000")
    
    ## Add arguments to define length and name of generated chromosome
    parser.add_argument("--chrom", type=str, default="chr1", help="chromosome name; default: chr1")
    parser.add_argument("--chromosome-length", type=int, default=2_000_000, help="Chromosome length; default: 2_000_000")
    
    ## Add randomness and output file name
    parser.add_argument("--random-seed", type=int, default=42, help="Random seed for simulation; default: 42")
    parser.add_argument("--output-vcf", type=str, default="simulated.vcf", help="Output VCF file name; default: simulated.vcf")
    
    args = parser.parse_args()
    simulate_and_write_vcf(args)

if __name__ == "__main__":
    main()
