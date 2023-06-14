#!/usr/bin/env python3
#
# This script parses a file with position \t
# and outputs a file with continious positions, 
# assigning the values from the previously encoutered position  


import argparse


__author__ = "Ekaterina Osipova, 2023."


def main():
    ## Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filein', type=str, help='input file with two columns to parse: position\t value')
    args = parser.parse_args()

    ## Read input bed file line by line;
    count = 0
    with open(args.filein) as inf:
        for line in inf:
            position_curr = int(line.split()[0])  
            value_curr = line.split()[1]

            if count == 0:
            	position_previous = position_curr
            	value_previous = value_curr
            	print('{}\t{}'.format(position_previous, value_previous))

            else:
	            if position_curr > position_previous:
	            	for i in range(position_curr - position_previous):
	            		print('{}\t{}'.format(position_previous + i, value_previous))
	            	position_previous = position_curr
	            	value_previous = value_curr
            count += 1



if __name__ == "__main__":
    main()
