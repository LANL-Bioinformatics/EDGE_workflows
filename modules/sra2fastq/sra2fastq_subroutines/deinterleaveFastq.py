
def deinterleaveFastq(input_file: str, output_file1: str, output_file2: str, compress: bool):
    if compress:
        file1 = open(output_file1, 'wt', compresslevel=1)
        file2 = open(output_file2, 'wt', compresslevel=1)
    else:
        file1 = open(output_file1, 'w')
        file2 = open(output_file2, 'w')

    with open(input_file, 'r') as input_fh:
        while True:
            try:
                line1 = next(input_fh)
                line2 = next(input_fh)
                line3 = next(input_fh)
                line4 = next(input_fh)

                file1.write(line1)
                file1.write(line2)
                file1.write(line3)
                file1.write(line4)

                line5 = next(input_fh)
                line6 = next(input_fh)
                line7 = next(input_fh)
                line8 = next(input_fh)

                file2.write(line5)
                file2.write(line6)
                file2.write(line7)
                file2.write(line8)
            except StopIteration:
                break

    file1.close()
    file2.close()

# # Usage example
# input_file = 'interleaved.fastq'
# output_file1 = 'f.fastq'
# output_file2 = 'r.fastq'
# compress_output = False  # Set to True if you want to compress the output files using pigz

# deinterleave_fastq(input_file, output_file1, output_file2, compress=compress_output)
