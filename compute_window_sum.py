#!/usr/bin/env python3
#
# This script parses a large file: chr \t position \t value1 \t value2 \t value3..
# and computes window-based sum or mean over all values


import argparse
from operator import itemgetter


__author__ = "Ekaterina Osipova, 2023."


def main():
    ## Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument(
                        '-i',
                        '--input', 
                        type=str, 
                        help='input file with at least a numeric column'
                        )
    parser.add_argument(
                        '-w', 
                        '--window', 
                        type=int, 
                        default=1000, 
                        help='window size; default=1000'
                        )
    parser.add_argument(
                        '-c', 
                        '--columns', 
                        type=str, 
                        default='3', 
                        help='comma-separated string of column numbers to compute sum/mean across'
                        )
    parser.add_argument(
                        '-m', 
                        '--mean', 
                        action='store_true', 
                        help='specify if you want to compute mean instead of sum'
                        )
    parser.add_argument(
                        '-p', 
                        '--position', 
                        type=str, 
                        default='1,2', 
                        help='comma-separated string of column numbers specifying coordinates'
                        )
    args = parser.parse_args()

    ## Read input line by line
    window = args.window
    line_count = 1
    window_value = 0
    check_columns = [int(c) - 1 for c in args.columns.split(',')]
    position_columns = [int(c) - 1 for c in args.position.split(',')]

    with open(args.input) as inf:
        for line in inf:
            
            if len(check_columns) == 1:
                line_values = float(line.split()[check_columns[0]])
                line_total = line_values
            elif len(check_columns) > 1:
                line_values = [float(i) for i in itemgetter(*check_columns)(line.split())]
                line_total = sum(line_values)
            
            if len(position_columns) == 1:
                coordinate = line.split()[position_columns[0]]
            elif len(position_columns) > 1:
                coordinate = '\t'.join([i for i in itemgetter(*position_columns)(line.split())])

            if args.mean:
                window_value += line_total / window
            else:
                window_value += line_total

            ## check if we went beyond the window
            if line_count >= window:
                print('{}\t{}'.format(coordinate, window_value))
                line_count = 0
                window_value = 0
            line_count += 1



if __name__ == "__main__":
    main()

