import os
import argparse
import sys
from multiprocessing import Pool

def process_file(args):
     # extract arguments:
    input_fasta_dir, input_fasta_suffix, output_dir, output_file_suffix, file = args
    # base file name doesn't have path or extension
    base_file_name = os.path.basename(file).replace('.'+str(input_fasta_suffix), '')
    # this is just trf command but in parallel
    os.system(f'trf {os.path.join(input_fasta_dir, base_file_name)}.{input_fasta_suffix} \
              2 7 7 80 10 50 500 -f -d -m -h -ngs \
              > {os.path.join(output_dir, base_file_name)}.{output_file_suffix}')

if __name__ == "__main__":
    if len(sys.argv) <7:
        print("Usage: python3 parallel_trf.py <input_file_dir> <input_file_suffix> \n"
              "\t <trf_output_dir> <trf_output_file_suffix> <num_processes> \n"
              "\t<file1> <file2> ...")
        sys.exit(1)
    parser = argparse.ArgumentParser(description='Process some FASTA files.')
    # directory with temporary FASTA files:
    parser.add_argument('--input-file-dir', required=True, help='the input directory')
    # input file suffix (trimmed.fasta):
    parser.add_argument('--input-file-suffix', required=True, help='the input file suffix (for example, trimmed.fasta)')
    # directory of TRF output files that contain names of reads with repeats:
    parser.add_argument('--trf-output-dir', required=True, help='the output directory')
    # output file suffix (trf_out.dat):
    parser.add_argument('--trf-output-file-suffix', required=True, help='the output file suffix (for example, trf_out.dat)')
    # number of subprocesses (CPU cores):
    parser.add_argument('--nthreads', type=int, default=1, help='the number of threads to use')
    # file names: object is python list, i.e a list of file names:
    parser.add_argument('--files', nargs='+', required=True, help='the FASTA files to process')
    args = parser.parse_args()
    
    with Pool(args.nthreads) as p:
        p.map(process_file, [(args.input_file_dir, args.input_file_suffix, \
                              args.trf_output_dir, args.trf_output_file_suffix, file) for file in args.files])
