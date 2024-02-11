import sys
import logging
import gzip
import shutil
import os

logger = logging.getLogger(__name__)

def open_file(file_path, mode='r'):
    """Open a file with gzip support if the file has a .gz extension."""
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)

def determine_output_format(input_file):
    """Determine the output file format based on the input file format."""
    if input_file.endswith('.gz'):
        return 'fastq.gz'
    else:
        return 'fastq'
        
def check_file_existence(file_path):
    """Check if a file exists."""
    return os.path.exists(file_path)

def remove_repeats_from_fastq(input_fastq, trf_output, output_fastq):
    # Check if input files exist
    if not check_file_existence(input_fastq) or not check_file_existence(trf_output):
        logger.error("Input file not found. Please check the file paths.")
        sys.exit(1)
    """ Remove the sequences from TRF that contain repeats from the output files """
    
    """This code reads the TRF output file and extracts lines that start with "@" 
    (indicating the beginning of a sequence in a FASTQ file). 
    It adds these lines to a set called sequences_with_repeats. 
    This set will contain the identifiers of sequences that were identified 
    by TRF as having repeats.
    """
    sequences_with_repeats = set()
    try:
        with open_file(trf_output) as file_handle:
            for line in file_handle:
                # sequences start with "@"
                if line[0] == "@":
                    sequences_with_repeats.add(line)
    except EnvironmentError:
        pass
                
    try:
        """
        This code tries to open the output file (output_fastq) for writing. 
        If it fails, it exits the program with an error message.
        The script then reads the input FASTQ file in chunks of four lines 
        (each set representing a single entry in the FASTQ file). 
        It checks if the sequence identifier (the first line in the set, 
        indicated by lines[0]) is in the sequences_with_repeats set. 
        If it is, the sequence is considered to have repeats, and removed_sequences 
        is incremented. 
        Otherwise, the non-repeated sequence is written to the output file.
        """
        output_format = determine_output_format(output_fastq)
        output_temp = output_fastq + '.tmp'
        with open(output_temp, 'w') as file_handle_write:
            removed_sequences = 0
            for lines in read_file_n_lines(input_fastq, 4):
                if lines[0] not in sequences_with_repeats:
                    file_handle_write.write("".join(lines))
                else:
                    removed_sequences += 1
            
            logger.info("Total number of sequences with repeats removed from file (%s): %d",
                        input_fastq, removed_sequences)
        # Compress the output file if needed
        if output_format == 'fastq.gz':
            with open_file(output_fastq, 'wb') as gz_file, open(output_temp, 'rb') as plain_file:
                shutil.copyfileobj(plain_file, gz_file)
            
            # Remove the temporary plain text file
            os.remove(output_temp)
        else:
            # Rename the temporary plain text file to the final output file
            os.rename(output_temp, output_fastq)

    except EnvironmentError:
        logger.error("Unable to open or write to file: %s", output_fastq)


def read_file_n_lines(file, n):
    """Read a file n lines at a time, using a generator to yield each set of lines."""
    line_set = []
    with open_file(file, 'rb') as file_handle:
        for line in file_handle:
            if isinstance(line, bytes):
                line = line.decode('utf-8')  # Decode bytes to string
            line_set.append(line)
            if len(line_set) == n:
                yield line_set
                line_set = []


if __name__ == '__main__':
    # Set up logging configuration
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    
    # Check if the required number of command-line arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python remove_repeats_from_fastq.py <input_fastq> <trf_output> <output_fastq>")
        sys.exit(1)  # Exit with a non-zero status to indicate an error
    
    # Specify command-line arguments
    fastq_in = sys.argv[1]
    trf_out = sys.argv[2]
    fastq_out = sys.argv[3]
    
    # Call the function
    remove_repeats_from_fastq(fastq_in, trf_out, fastq_out)
