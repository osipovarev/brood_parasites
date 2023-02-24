#!/usr/bin/env python3

'''
This script prepares samples.csv file specific for snpArcher input:
BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refPath
'''

import argparse
import re
import sys
import os



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

    # go to the requested dir to put links into
    os.chdir(datadir)

    # create a symlink for each file
    for f in file_list:
        f_name = f.split('/')[-1]
        if not os.path.islink(f_name):
            os.symlink(f, f_name)

    # got back to wdir
    os.chdir(wdir)


def get_ids_from_sample_list(file_list, extensions):
    ## Parses the list of file names and returns lists of IDs
    ## (removes extensions)

    id_list = [f.split('/')[-1] for f in file_list]
    for ext in extensions:
        id_list = [i.replace(ext, '') for i in id_list]
    id_list = list(set([re.sub('_[a-zA-Z]?2$', '', re.sub('_[a-zA-Z]?1$', '', i)) for i in id_list]))
    return id_list


def write_samples_file(id_list, file_list, ref_path, name_sci, project, assembly_ref, srx):
    ## Writes to stdout; uses format:
    ## BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refGenome

    # write header
    print('BioSample,LibraryName,refGenome,Run,Organism,BioProject,fq1,fq2,refPath')

    # define SRX format
    n_digit = len(srx.replace('SRX', ''))
    srx_1 = int(re.sub('SRX0+', '', srx))

    # write samples
    for i in range(len(id_list)):        
        sample_id = id_list[i]
        fq_files = ','.join(sorted([f for f in file_list if sample_id in f]))
        samples_info = '{},{}'.format(sample_id, sample_id)
        srx_i = 'SRX' + '0' * (n_digit - len(str(srx_1 + i))) + str(srx_1 + i)
        ref_info = '{},{},{},{}'.format(assembly_ref, srx_i, name_sci, project)
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
                        '-lf',
                        '--local_fastq',
                        type=str,
                        default='data/local_fastq/',
                        help="Relative path to dir to put links to fastq files; default: data/local_fastq/"
                         )
    parser.add_argument(
                        '-lr',
                        '--local_ref',
                        type=str,
                        default='data/local_genome',
                        help="Relative path to dir to put links to reference genome; default: data/local_genome/"
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
                        default='SRX00000001',
                        help="SRX ID of the first sample (like SRX15327220) if known; default: SRX0000000$i"
                        )
    args = parser.parse_args()


    ## Parse arguments
    fastq_dir = args.fastq_dir
    ref_path = args.ref_path
    local_fastq = args.local_fastq
    local_ref = args.local_ref
    name_sci = args.name_sci
    assembly_ref = args.assembly_ref
    project = args.project
    srx = args.srx


    ## Parse fastq files -> get sample IDs
    extensions = ['.gz', '.fastq']
    file_list = parse_dir_files(fastq_dir, extensions)
    id_list = get_ids_from_sample_list(file_list, extensions)


    ## Make symbolic links to reseq fastq files in data/local_fastq
    make_file_links(file_list, local_fastq)


    ## Make symbolic links to ref genomic fasta in data/local_genome
    make_file_links([ref_path], local_ref)


    ## Write samples lines in snpArcher-prefered format
    local_file_list = [os.getcwd() + '/' + local_fastq + f.split('/')[-1] for f in file_list]
    local_ref_path = os.getcwd() + '/' + local_ref + '/' + ref_path.split('/')[-1]
    write_samples_file(id_list, local_file_list, local_ref_path, name_sci, project, assembly_ref, srx)


if __name__ == "__main__":
    main()
