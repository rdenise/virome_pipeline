from Bio import SeqIO
import argparse
import os
import subprocess
import io
import multiprocessing as mp
import shutil
import pandas as pd
from pathlib import Path

# Put error and out into the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

###########################################################
###########################################################

def create_folder(mypath):

    """
    Created the folder that I need to store my result if it doesn't exist
    :param mypath: path where I want the folder (write at the end of the path)
    :type: string
    :return: Nothing
    """

    try:
        os.makedirs(mypath)
    except OSError:
        pass

    return

###########################################################

def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("--in_fasta_file",
                        help="input_fasta", 
                        type=str, 
                        default=snakemake.input.contig)
    parser.add_argument('--outfile', 
                        help='output_folder', 
                        type=str, 
                        default=snakemake.output.blast_out)
    parser.add_argument('--tmp_dir', 
                        help='tmp folder', 
                        type=str, 
                        default=snakemake.params.tmp_output)
    parser.add_argument('--database', 
                        help='The name of the database or empty if remote blast', 
                        type=str, 
                        default=snakemake.params.database)
    parser.add_argument('--contigs_per_file', 
                        help='The number of contigs in each slice', 
                        type=int, 
                        default=1)
    parser.add_argument('--job_number', 
                        type=int, 
                        default=snakemake.threads)
    # Parse arguments
    args = parser.parse_args()

    return args

###########################################################

def main(args):
    
    #split the input file
    record_iter = SeqIO.parse(open(args.in_fasta_file), "fasta")
    files_to_run = []


    for i, batch in enumerate(batch_iterator(record_iter, args.contigs_per_file)):
        try:
            group_name=f'group_{i+1}'
            dir_name=f'{args.tmp_dir}/{group_name}'
            filename = f"{dir_name}/{group_name}.fasta"

            create_folder(f"{args.tmp_dir}/{group_name}")

            with open(filename, "w") as handle:
                count = SeqIO.write(batch, handle, "fasta")
                print(f"Wrote {count} records to {filename}")
                files_to_run.append((group_name, dir_name, filename))
        except IOError as ioerror:
            print(ioerror)

    if args.database:
        if not (Path(args.database + ".00.nsq").exists() or Path(args.database + ".nsq").exists()):
            cmd_str = f"makeblastdb -dbtype nucl -in '{args.database}'"
            stdout, stderr = execute(cmd_str)
            print(f"----makeblastdb - stdout----\n{stdout.decode('utf8')}\n----makeblastdb - stderr----\n{stderr.decode('utf8')}\n")            

    pool = mp.Pool(args.job_number)
    results = pool.map(run_job, files_to_run)
    pool.close()

    df = pd.concat(results)
    df.to_csv(args.outfile, sep='\t', index=False, header=False)

    print(f'{bcolors.OKBLUE} ------ DONE! ----------- {bcolors.ENDC}')
    shutil.rmtree(args.tmp_dir)

###########################################################

def run_job(group_tuple):
    if args.database:
        blast_remote=''
        blast_database = args.database
    else:
        blast_remote = "-task blastn -remote"
        blast_database = "nt"

    job_str=f"blastn -query {group_tuple[2]} -out {group_tuple[1]}/blast-output.txt " \
            f"-db '{blast_database}' -evalue 0.0001 -taxids 9606 {blast_remote} " \
            "-outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' " \
            "-word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes " \
            "-min_raw_gapped_score 100 -perc_identity 90 -soft_masking true -max_target_seqs 10 "

    if os.path.isfile(f'{group_tuple[1]}/blast-output.txt'):
        stdout = 'File already exists' 
        stderr = ''
    else: 
        stdout, stderr = execute(job_str)
    print(f"----BLASTn - stdout----\n{stdout.decode('utf8')}\n----BLASTn - stderr----\n{stderr.decode('utf8')}\n")
    if os.path.isfile(f'{group_tuple[1]}/blast-output.txt') and os.path.getsize(f'{group_tuple[1]}/blast-output.txt'):
        df = pd.read_csv(f'{group_tuple[1]}/blast-output.txt', sep='\t', header=None)
    else:
        df = pd.DataFrame()
    return df

###########################################################

def execute(command):
    print(f'Executing {command}')
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    stdout, stderr = process.communicate()
    return stdout, stderr

###########################################################

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

###########################################################

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

###########################################################

if __name__ == "__main__":
    args = parse_args()
    main(args)
