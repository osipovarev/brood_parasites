#!/usr/bin/env python3

'''
This script prepares samples.csv file specific for snpArcher input:
BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refPath
'''

import argparse
import sys
import os



def parse_dir_files(dir_path, extensions):
    ## Makes a list of file names in the requested dir,
    ## includes only files with extensions provided

    file_list = []
    for filename in os.listdir(dir_path):
        f = os.path.join(dir_path, filename)
        if (os.path.isfile(f)) and ('.'+f.split('.')[-1] in extensions):
            file_list.append(f)
    return file_list


def get_ids_from_sample_list(file_list, extensions):
    ## Parses the list of file names and returns lists of IDs
    ## (removes extensions)

    id_list = file_list
    for ext in extensions:
        id_list = [i.replace(ext, '') for i in id_list]
    return id_list


def write_samples_file(id_list, file_list, ref_path, name_sci, project, assembly_ref, srx):
    ## Writes to stdout; uses format:
    ## BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refGenome

    # write header
    print('BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refPath')

    # write samples
    for i in id_list:
        fq_files = ','.join([f for f in file_list if f.startswith(i)])
        samples_info = '{},{}'.format(i, i)
        ref_info = '{},{},{},{}'.format(assembly_ref, srx, name_sci, project)
        full_sample_line = '{},{},{},{}'.format(samples_info, ref_info, fq_files, ref_path)
        print(full_sample_line)



def main():
    parser = argparse.ArgumentParser(description='Write sample files.')
    parser.add_argument('-f', '--fastq_dir', required=True, help="Path to dir with fastq files")
    parser.add_argument('-r', '--ref_path', required=True, help="Path to reference genome")
    parser.add_argument('-n', '--name_sci', required=True, help="Scientific name of the organism")
    parser.add_argument('-a', '--assembly_ref', default='GCA_xxxxxxxxx', help="NCBI ref assembly ID \
                        (like GCA_013399945.1) if known; default: GCA_xxxxxxxxx")
    parser.add_argument('-p', '--project', default='PRJNAxxxxxx', help="Project ID (like PRJNA839346) if known;\
                        default: PRJNAxxxxxx")
    parser.add_argument('-s', '--srx', default='SRXxxxxxxxx', help="SRX ID of the first sample (like SRX15327220) if known;\
                        default: SRXxxxxxxxx")

    args = parser.parse_args()

    ## Parse arguments
    fastq_dir = args.fastq_dir
    ref_path = args.ref_path
    name_sci = args.name_sci
    assembly_ref = args.assembly_ref
    project = args.project
    srx = args.srx

    ## Parse fastq files -> get sample IDs
    extensions = ['.gz', '.fastq']
    file_list = parse_dir_files(fastq_dir, extensions)
    id_list = get_ids_from_sample_list(file_list, extensions)

    ## Make symbolic links to reseq fastq files in data/

    ## Make symbolic links to ref genomic fasta in data/

    ## Write samples lines in snpArcher-prefered format
    write_samples_file(id_list, file_list, ref_path, name_sci, project, assembly_ref, srx)


if __name__ == "__main__":
    main()