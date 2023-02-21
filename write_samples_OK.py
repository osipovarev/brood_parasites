#!/usr/bin/env python3

'''
This script prepares samples.csv file specific for snpArcher input:
BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refPath
'''

import argparse
import re
import sys
import os
import subprocess



def parse_dir_files(dir_path, extensions):
    ## Makes a list of file names in the requested dir,
    ## includes only files with extensions provided

    file_list = []
    for filename in os.listdir(dir_path):
        f = os.path.join(dir_path, filename)
        if '.' + f.split('.')[-1] in extensions:
            file_list.append(f)
    return file_list


def make_file_links(file_list, datadir):
    ## Makes links to files in requested directory 

    wdir = os.getcwd()

    # create requested datadir if does not exist
    if not os.path.exists(datadir):
        os.makedirs(datadir)
        print("Created directory: {}".format(datadir))

    # go to the requested dir to put links into
    subprocess.run(['cd', '{}'.format(datadir)])
    print('Now in: '.format(os.getcwd()))

    # create a symlink for each file
    for f in file_list:
        f_name = f.split('/')[-1]
        print(f)
        print(f_name)
        # subprocess.run(['ln', '-s', f, f_name])

    # got back to wdir
    subprocess.run(['cd', '{}'.format(wdir)])
    print('Now in: '.format(os.getcwd()))



def get_ids_from_sample_list(file_list, extensions):
    ## Parses the list of file names and returns lists of IDs
    ## (removes extensions)

    id_list = [f.split('/')[-1] for f in file_list]
    for ext in extensions:
        id_list = [i.replace(ext, '') for i in id_list]
    id_list = list(set([re.sub('_2$', '', re.sub('_1$', '', i)) for i in id_list]))
    return id_list


def write_samples_file(id_list, file_list, ref_path, name_sci, project, assembly_ref, srx):
    ## Writes to stdout; uses format:
    ## BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refGenome

    # write header
    print('BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refPath')

    # write samples
    for i in range(len(id_list)):        
        sample_id = id_list[i]
        fq_files = ','.join(sorted([f for f in file_list if sample_id in f]))
        samples_info = '{},{}'.format(sample_id, sample_id)
        ref_info = '{},{},{},{}'.format(assembly_ref, srx+str(i), name_sci, project)
        full_sample_line = '{},{},{},{}'.format(samples_info, ref_info, fq_files, ref_path)
        print(full_sample_line)



def main():
    parser = argparse.ArgumentParser(description='Write sample files.')
    parser.add_argument(
                        '-f',
                        '--fastq_dir',
                        required=True,
                        help="Full path to dir with fastq files"
                         )
    parser.add_argument(
                        '-r',
                        '--ref_path',
                        required=True,
                        help="Full path to reference genome"
                           )
    parser.add_argument(
                        '-n',
                        '--name_sci',
                        required=True,
                        help="Scientific name of the organism"
                           )
    parser.add_argument(
                        '-a',
                        '--assembly_ref',
                        default='GCA_000000000',
                        help="NCBI ref assembly ID (like GCA_013399945.1) if known; default: GCA_000000000"
                        )
    parser.add_argument(
                        '-p',
                        '--project',
                        default='PRJNA000000',
                        help="Project ID (like PRJNA839346) if known; default: PRJNA000000"
                        )
    parser.add_argument(
                        '-s',
                        '--srx',
                        default='SRX0000000',
                        help="SRX ID of the first sample (like SRX15327220) if known; default: SRX0000000$i"
                        )
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


    ## Make symbolic links to reseq fastq files in data/local_fastq
    datadir='data/local_fastq'
    make_file_links(file_list, datadir)


    ## Make symbolic links to ref genomic fasta in data/local_genome
    datadir='data/local_genome'
    make_file_links(ref_path, datadir)


    ## Write samples lines in snpArcher-prefered format
    write_samples_file(id_list, file_list, ref_path, name_sci, project, assembly_ref, srx)


if __name__ == "__main__":
    main()
