#!/usr/bin/env python3
#
# This script parses bed file with chr \t position
# and outputs bed file with continous intervals of requested size

import argparse


__author__ = "Ekaterina Osipova, 2023."


def main():
    ## Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', type=str, help='bed with two columns to parse: chr \t position')
    parser.add_argument('-s', '--size', type=int, help='max size of the window WITHOUT bed entires: split points')
    args = parser.parse_args()

    ## Read input bed file line by line;
    count = 0
    with open(args.bed) as inf:
        for line in inf:
            
            scaffold_next = line.split()[0]
            pos_next = int(line.split()[1])
            
            if count == 0:
                pos_0 = pos_next
                pos_current = pos_0
                scaffold_current = scaffold_next
            else:
                dist = pos_next - pos_current
                if (dist >= args.size) or (scaffold_next != scaffold_current):
                    ## out of the window! => print if interval big enough
                    if pos_current - pos_0 > args.size:
                        print('{}\t{}\t{}'.format(scaffold_current, pos_0, pos_current))
                    scaffold_current = scaffold_next
                    pos_0 = pos_next
                    pos_current = pos_next
                else:
                    ## within window; continue
                    pos_current = pos_next
            count += 1



if __name__ == "__main__":
    main()
