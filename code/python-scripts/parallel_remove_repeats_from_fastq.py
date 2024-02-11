import os
import argparse
import sys
from multiprocessing import Pool

def process_file(args):
    # extract arguments
    python_script,input_fastq_dir,input_fastq_suffix, \
        trf_dat_dir, trf_dat_dir_suffix, \
            output_fastq_dir,output_fastq_suffix, file = args
    # base file name doesn't have path or extension
    base_file_name = os.path.basename(file).replace('.'+str(input_fastq_suffix), '')
    # this is just remove_repeats_from_fastq.py but in parallel
    os.system(f'python3 {python_script} \
              {os.path.join(input_fastq_dir, base_file_name)}.{input_fastq_suffix} \
              {os.path.join(trf_dat_dir, base_file_name)}.{trf_dat_dir_suffix} \
              {os.path.join(output_fastq_dir, base_file_name)}.{output_fastq_suffix}')

if __name__ == "__main__":
    if len(sys.argv) <10:
        print("Usage: python3 parallel_remove_repeats_from_fastq.py \n"
              "\t<remove_repeats.py> <input_file_dir> <input_file_suffix> \n"
              "\t<trf_file_dir> <trf_file_suffix> <output_dir> <output_file_suffix> \n"
              "\t<nthreads> <file1> <file2> ...")
        sys.exit(1)
    parser = argparse.ArgumentParser(description='Process some FASTQ/gzipped FASTQ files.')
    # script that removes repeats (not parallel):
    parser.add_argument('--remove-repeats-script', required=True, help='script that removes repeats, not parallel')
    # directory with trimmed FASTQ/gzipped FASTQ files from which we remove repeats:
    parser.add_argument('--input-file-dir', required=True, help='the input directory')
    # input file suffix (fastq.gz):
    parser.add_argument('--input-file-suffix', required=True, help='the input file suffix (for example, fastq.gz)')
    # directory of TRF dat files that contain names of reads with repeats:
    parser.add_argument('--trf-file-dir', required=True, help='the directory with TRF dat files')
    # TRF dat file suffix (trf_out.dat):
    parser.add_argument('--trf-file-suffix', required=True, help='the TRF dat file suffix (for example, trf_out.dat)')
    # directory of clean FASTQ/gzipped FASTQ files files without repeats:
    parser.add_argument('--output-dir', required=True, help='the output directory')
    # output file suffix (clean.fastq.gz)
    parser.add_argument('--output-file-suffix', required=True, help='the input file suffix (for example, clean.fastq.gz)')
    # number of subprocesses (CPU cores):
    parser.add_argument('--nthreads', type=int, default=1, help='the number of threads to use')
    # file names: object is python list, i.e a list of file names:
    parser.add_argument('--files', nargs='+', required=True, help='the FASTQ/gzipped FASTQ files to process')
    args = parser.parse_args()

    with Pool(args.nthreads) as p:
        p.map(process_file, 
              [(args.remove_repeats_script,args.input_file_dir, args.input_file_suffix, \
                args.trf_file_dir,args.trf_file_suffix, \
                    args.output_dir,args.output_file_suffix,file) for file in args.files])
        