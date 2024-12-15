#!/usr/bin/env python3

import sys
import os
import multiprocessing
import argparse
import time

# start print time taken
start_time = time.time()

# read GFA file, reads file, number of threads, prefix file name and chunk size as arguments
parser = argparse.ArgumentParser(description='Run PHI2')
# use -g to specify the GFA file
parser.add_argument('-g', type=str, help='GFA file')
# use -r to specify the reads file
parser.add_argument('-r', type=str, help='reads file')
# use -t to specify the number of threads
parser.add_argument('-t', type=int, help='number of threads')
# use -p to specify the prefix file name
parser.add_argument('-p', type=str, help='prefix file name')
# use -c to specify the chunk size
parser.add_argument('-c', type=int, help='chunk size')
# use -b to specify the number of batches
parser.add_argument('-b', type=int, help='number of batches')
args = parser.parse_args()

# throw help message if no arguments are provided
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

# print the arguments with their values
print("GFA file: " + args.g)
print("Reads file: " + args.r)
print("Number of threads: " + str(args.t))
print("Prefix file name: " + args.p)
print("Chunk size: " + str(args.c))
print("Number of batches: " + str(args.b))


# run divide_gfa_n.py
# ./divide_gfa_n.py GFA_file n prefix
os.system("./divide_gfa_n.py " + args.g + " " + str(args.c) + " " + "GFA_" + args.p)

# run PHI2 as ./PHI2 -t {threads} -g GFA -r reads -o prefix_{i}
# where i is the chunk number 
# run in parallel with the number of threads specified with multiprocessing

# if logs folder does not exist, create it
if not os.path.exists("logs"):
    os.makedirs("logs")

def run_PHI2(i):
    os.system("./PHI2 -t " + str(args.t) + " -g GFA_" + args.p + "_" + str(i + 1) + ".gfa" + " -r " + args.r + " -o " + args.p + "_" + str(i + 1) + " 2> logs/log_" + args.p + "_" + str(i + 1) + ".txt")
    
pool = multiprocessing.Pool(args.t)
pool.map(run_PHI2, range(args.b))
pool.close()
pool.join()

# run join_haps.py prefix n
os.system("./join_haps.py " + args.p + " " + str(args.c))

# run PHI2 agian with prefix.gfa and reads file
os.system("./PHI2 -R0 -t " + str(args.t) + " -g " + args.p + ".gfa" + " -r " + args.r + " -o " + args.p + " 2> logs/log_" + args.p + "_join.txt")

# end print time taken
end_time = time.time()
print("Time taken: " + str(end_time - start_time) + " seconds")