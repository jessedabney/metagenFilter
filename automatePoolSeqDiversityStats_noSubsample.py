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
        required=False)
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
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/basic-pipeline/filter-pileup-by-gtf.pl --input {0}.noIndel.mpileup --gtf {1} --output {0}.filtered.mpileup".format(args.prefix, args.removereg))
    print("executing removal of rep regions")

def SubSample(i):
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/basic-pipeline/subsample-pileup.pl --input {0}_classSeqs_noIndel_remRegs_noSB.mpileup --output {0}_rand{1}.mpileup --target-coverage {2} --max-cov 1000000 --min-qual 20 --fastq-type sanger --method withoutreplace".format(args.prefix, i, args.coverage))
    print("executing random subsampling {0}".format(i))

def pi():
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/Variance-sliding.pl --measure pi --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --window-size 100000 --step-size 10000 --input {0}_classSeqs_noIndel_remRegs_noSB.mpileup --output {0}_noSubsample_w100K_n10K.pi --snp-output {0}_noSubsample.snps".format(args.prefix))
    print("executing pi calculation")

def theta():
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/Variance-sliding.pl --measure theta --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --window-size 100000 --step-size 10000 --input {0}_classSeqs_noIndel_remRegs_noSB.mpileup --output {0}_noSubsample_w100K_n10K.theta".format(args.prefix))
    print("executing theta calculation")

def td():
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/Variance-sliding.pl --measure D --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --window-size 100000 --step-size 10000 --input {0}_classSeqs_noIndel_remRegs_noSB.mpileup --output {0}_noSubsample_w100K_n10K.td ".format(args.prefix))
    print("executing td calculation")

def gene_pi():
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --pileup {0}_classSeqs_noIndel_remRegs_noSB.mpileup --output {0}_noSubsample_gene.pi --snp-output {0}_noSubsample.gene.snps --gtf {1}".format(args.prefix, args.gtf))
    print("executing pi gene-wise calculation")

def gene_theta():
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/Variance-at-position.pl --measure theta --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --pileup {0}_classSeqs_noIndel_remRegs_noSB.mpileup --output {0}_noSubsample_gene.theta --gtf {1}".format(args.prefix, args.gtf))
    print("executing theta gene-wise calculation")

def syn_nonsyn():
    call_with_log("perl /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/Variance-at-position.pl --measure pi --pool-size 10000 --fastq-type sanger --min-count 2 --min-covered-fraction 0.5 --codon-table /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/syn-nonsyn/codon-table.txt --nonsyn-length-table /mnt/PepPop_export/PepPrograms/popoolation_1.2.2/syn-nonsyn/nsl_p1.txt --pileup {0}_classSeqs_noIndel_remRegs_noSB.mpileup --output {0}_noSubsample_syn-nonsyn.pi --snp-output {0}_noSubsample_syn-nonsyn.snps --gtf {1}".format(args.prefix, args.gtf))
    print("executing pi syn-nonsyn calculation")

kvmap= {'prefix':args.prefix}

#remReg()
#for i in range(10) :
    #SubSample(i)
pi()
theta()
td()
gene_pi()
gene_theta()
    #syn_nonsyn(i)
