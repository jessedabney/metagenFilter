#!/usr/bin/env python

import argparse
import sys
import os
import getopt
import glob
import subprocess
import shlex
from datetime import datetime

####################################################################
##This script is generic. Modify!
####################################################################


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest,
            os.path.abspath(os.path.expanduser(values)))

def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def is_dir(dirname):
    """Checks if a path is a directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname
        
def is_file(filename):
    """Checks if a file exists"""
    if not os.path.isfile(filename):
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename

def get_arguments(): 
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Merge individual sample VCFs into a project VCF")
    parser.add_argument('-p', '--prefix', required=True,
        help = 'sample prefix, eg ERR#, ERR1.noIndel.mpileup')
    parser.add_argument('-c', '--coverage', default=50,
        help = 'desired coverage for subsampling')
    parser.add_argument('-r', '--removereg', 
        help = 'gtf file containing regions to be removed',
        type = is_file,
        required=True)
    parser.add_argument('-g', '--gtf', 
        help = 'gtf file containing genes',
        type = is_file,
        required=True)
    return parser.parse_args()

args = get_arguments()

current_datetime = datetime.today().strftime("%d-%m-%Y-%H%M")

#calls a system command with the subprocess module
#redirects both stdout and stderr to a log file
def call_with_log(cmd):
    cmd = cmd.format(**(kvmap))

    logfile = open(current_datetime+".log", "a+")
    logfile.write("Executing command: " + cmd + "\n")
    logfile.flush()
    ret = subprocess.call(shlex.split(cmd), stdout=logfile, stderr=logfile)
    if(ret != 0):
        print("Script did not complete successfully. \n Command : \n\n" + cmd + "\n\n returned with non-zero code: " + str(ret))
        logfile.write("Script did not complete successfully. \n Command : \n\n" + cmd + "\n\n returned with non-zero code: " + str(ret))
        logfile.close()
        sys.exit(-1)
    logfile.close()

def remReg():
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {0}.noIndel.mpileup --gtf {1} --output {0}.filtered.mpileup".format(args.prefix, args.removereg))
    print("executing removal of rep regions")

def SubSample(i):
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --input {0}.filtered.mpileup --output {0}_rand{1}.mpileup --target-coverage {2} --max-cov 1000000 --min-qual 20 --fastq-type sanger --method withoutreplace".format(args.prefix, i, args.coverage))
    print("executing random subsampling {0}".format(i))

def pi(i):
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/Variance-sliding.pl --measure pi --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --window-size 100000 --step-size 10000 --input {0}_rand{1}.mpileup --output {0}_rand{1}_w100K_n10K.pi --snp-output {0}_rand{1}.snps".format(args.prefix, i))
    print("executing pi calculation on subsample {0}".format(i))

def theta(i):
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/Variance-sliding.pl --measure theta --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --window-size 100000 --step-size 10000 --input {0}_rand{1}.mpileup --output {0}_rand{1}_w100K_n10K.theta".format(args.prefix, i))
    print("executing theta calculation on subsample {0}".format(i))

def td(i):
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/Variance-sliding.pl --measure D --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --window-size 100000 --step-size 10000 --input {0}_rand{1}.mpileup --output {0}_rand{1}_w100K_n10K.td ".format(args.prefix, i))
    print("executing td calculation on subsample {0}".format(i))

def gene_pi(i):
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --pileup {0}_rand{1}.mpileup --output {0}_rand{1}_gene.pi --snp-output {0}_rand{1}.gene.snps --gtf {2}".format(args.prefix, i, args.gtf))
    print("executing pi gene-wise calculation on subsample {0}".format(i))

def gene_theta(i):
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/Variance-at-position.pl --measure theta --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --pileup {0}_rand{1}.mpileup --output {0}_rand{1}_gene.theta --gtf {2}".format(args.prefix, i, args.gtf))
    print("executing theta gene-wise calculation on subsample {0}".format(i))

def syn_nonsyn(i):
    call_with_log("perl /opt/PepPrograms/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --codon-table /opt/PepPrograms/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table /opt/PepPrograms/popoolation_1.2.2/syn-nonsyn/nsl_p1.txt --pileup {0}_rand{1}.mpileup --output {0}_rand{1}_syn-nonsyn.pi --snp-output {0}_rand{1}_syn-nonsyn.snps --gtf {2}".format(args.prefix, i, args.gtf))
    print("executing pi syn-nonsyn calculation on subsample {0}".format(i))

kvmap= {'prefix':args.prefix}

remReg()
for i in range(10) :
    #SubSample(i)
    #pi(i)
    #theta(i)
    td(i)
    gene_pi(i)
    gene_theta(i)
    syn_nonsyn(i)
