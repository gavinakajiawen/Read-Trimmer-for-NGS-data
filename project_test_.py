#!/usr/bin/env python3 
#--------------------------------------------------------------------
# Import Libraries
#--------------------------------------------------------------------
import gzip 
import re 
import sys
import argparse
import os 

#--------------------------------------------------------------------
# Declaration of  colours. 
#--------------------------------------------------------------------
bright_green = "\033[0;32m"
bright_cyan = "\033[0;96m"
bright_yellow = "\033[0;33m"
brightish = "\033[0;36m"
bright_purple ="\033[0;35m" 
reset = "\033[0m\n"
#--------------------------------------------------------------------
# INTRODUCTION
#--------------------------------------------------------------------
print(bright_purple + "#--------------------------------------------------------------------")
print(bright_green + "#--------------------------------------------------------------------" +reset)
print("Genvej: Read Trimmer for Next Generation Sequencing Data")
print("version 2.0.01 for Mac OS X 64-bit built November 22 2022")
print("Developed by M.O.L.S, Jiawen Wu")
print(bright_green + "#--------------------------------------------------------------------")
print(bright_purple + "#--------------------------------------------------------------------" + reset)
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# ARGPARSE 
#--------------------------------------------------------------------
# Create an argument parser. 
parser = argparse.ArgumentParser(prog = 'Genvej',  
    description = '''This program trims reads from compressed and uncompressed FASTQ files.''',
    
    epilog='''The reads can be trimmed from both the 3' and 5' ends.''')

# Change the name of the positional arguments.
parser._positionals.title = 'Positional arguments'

# Change the name of the optional arguments. 
parser._optionals.title = '**Optional arguments**'

#--------------------------------------------------------------------
# Add a non-optional argument :required = True
parser.add_argument("--filename", type= argparse.FileType('r'), required = False,
    help= "Enter the name of a FASTQ file for read trimming.")


parser.add_argument("--shear", type = int, required = False,
    help="Trims 3' prime end of a read by (x) amount of nucleotides")


parser.add_argument("--snip", type=int, required= False,
    help="Trims 5' prime end of a read by (x) amount of nucleotides")

#--------------------------------------------------------------------
# Add an optional argument
parser.add_argument('--name', help = "Receive a welcome message", 
    type=str, required=False)


# Add an optional argument
parser.add_argument('--version', 
    # nargs means have 0 or 1 argments 
    nargs='?', 
    # const means set the defalt when there are 0 arguments. 
    const=1,
    type=str,
    help="Display the version of the program",
    required=False)


parser.add_argument('--cite', 
    # nargs means have 0 or 1 argments 
    nargs='?', 
    # const means set the defalt when ther are 0 arguments. 
    const=1,
    type=str,
    help="Display the citation for the paper associated with the program",
    required=False)

#--------------------------------------------------------------------
# Take the argument from the command line.
args = parser.parse_args()

#--------------------------------------------------------------------
# Decide what to with the argument based on the input

if args.filename:
    print("The FASTQ file entered is:" + args.filename.name)
    
    
if args.shear:
    print("The number of nucleotides to trim from the 3' end: "+ str(args.shear))
    
if args.snip:
    print("The number of nucleotides to trim from the 5' end: "+ str(args.snip))
    
if args.version:
    print('The version of the program is version 2.0.01\n')
    
if args.cite:
    print('The citation is: \'Genvej: Read Trimmer for Next Generation Sequencing Data\'\n')
    
if args.name:
    print('Hello '+ args.name.capitalize() + 
        ", Welcome to Genvej: Read Trimmer for Next Generation Sequencing Data \n")	
    
    #--------------------------------------------------------------------
log_file = open("log.dat", "w")
#--------------------------------------------------------------------


#--------------------------------------------------------------------
# Does the filename end with .gz? 
#-------------------------------------------------------------------
try:
    search_gz = re.search(r'.gz$', args.filename.name)
except:
    # exit the program if there is no input file provided.
    pass
    sys.exit()
    
    
if search_gz is not None:
    print("This is a gzipped, compressed FASTQ file")
else:
    print("The FASTQ file is uncompressed.")
    
    
    
    #-------------------------------------------------------------------
    # uncompressed_file
    #-------------------------------------------------------------------
    
    # If the filename does not end in gz, open the file. 
if search_gz is None:
    try: 
        infile= open(args.filename.name, "r") 
    except FileNotFoundError as error_F:
        print(bright_green + "There is a minor problem with opening the FASTQ file.", reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:"+ args.filename.name)
        
        #-------------------------------------------------------------------
        
        #-------------------------------------------------------------------
        # shear and snip 
        #-------------------------------------------------------------------
        # Set the 3' and 5' ends to trim. 
three_prime_end = args.shear

five_prime_end = args.snip
#------------------------------------------------------------------

#--------------------------------------------------------------------
# compressed_file
#-------------------------------------------------------------------
# If the filename ends in gz, open the file using gzip. 
if search_gz is not None:
    try: 
        infile = gzip.open(args.filename.name, "r")
    except FileNotFoundError as error_F:
        print(bright_cyan + "There is a minor problem with opening the FASTQ file."+ reset)
        sys.exit(1)
    except gzip.BadGzipFile as error_geez:
        print(bright_cyan + "This file has the suffix of .gz, but it is actually not compressed."+ reset)
        print(bright_cyan+ "Rename the file by removing the .gz part. "+ reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:", str(args.filename.name))
        
        #-------------------------------------------------------------------
        # core functions
        #-------------------------------------------------------------------
def convert_thirty_three(s):
    ''' Take a string a convert to a phredd 33 quality score list'''
    M = []
    
    for letter in s:
        M.append(ord(letter) - 33)
        
    return M

def convert_thirty_three(s):
    ''' Take a string a convert to a phredd 33 quality score list'''
    M = []
    
    for letter in s:
        M.append(ord(letter) - 33)
        
    return M


def convert_back_thirty_three(number):
    '''Takes an integer and changes it to symbol letters'''
    M  = str(chr(number + 33))
    
    return M
#-------------------------------------------------------------------

print(bright_purple + "#--------------------------------------------------------------------#")
print(" STAGE 1")
print(bright_purple + "#--------------------------------------------------------------------#" + reset)

#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Handle an uncompressed file and determine the encoding. 

line_number = 0
phred_line = ""
seq_line = ""

if search_gz is None: 
    
    for line in infile:
        line_number += 1
        
        if line_number == 2:
            seq_line += line.strip("\n")
            
        if line_number == 4:
            phred_line += line.strip("\n")
            
infile.close()
#-------------------------------------------------------------------
#------------------------------------------------------------------
if search_gz is None:
    print("Currently", str(int(line_number/4)), "read(s) in this file.")
    print("The length of the untrimmed sequence is:", len(seq_line))
    print("The phred line is:", phred_line)
    
    # search for characters only in phred33 
    search_phred33 = re.search(r'([#=!?>;$&%+()/.*":<]|[0-9])', phred_line)
    
    if search_phred33 is not None:
        #print(search_phred33.group(1))
        print("This FASTQ file called has PHRED 33 encoding (SANGER)")
        
        # search for characters only in phred64
    search_phred64 = re.search(r'([L-Z]|[a-j])', phred_line)
    if search_phred64 is not None:
        #print(search_phred64.group(1))
        print("This FASTQ file has PHRED 64 encoding (OLD ILLUMINA)")
        
        #------------------------------------------------------------------
        #--------------------------------------------------------------------
        # phred 33 and uncompressed.
        # read each line of the file as a string into a list. 
        # Every fourth line is the quality line. 
        
header_list=[]
dna_sequences=[]
plus_signs=[]
quality_sequence=[]

line_number= 0

if search_gz is None and search_phred33 is not None:
    
    infile = open(args.filename.name, "r") 
    
    for line in infile:
        line_number += 1
        
        if line.startswith("@HWI"):
            header_list.append(line.strip())
            
            dna_sequences.append("")
            
        search_dna = re.search(r'([ACTGN]{101})', line)
        
        if search_dna is not None:
            dna_sequences[-1] += search_dna[1]
            
            plus_signs.append('')
            
        search_plus = re.search(r'^\+\s', line)
        
        if search_plus is not None:
            plus_signs[-1] += line
            
            quality_sequence.append('')
            
        if line_number % 4 == 0:
            search_quality = re.search(r'(.{101})', line, flags = re.ASCII)
            if search_quality is not None:
                quality_sequence[-1] += search_quality[1]
                
infile.close()
#--------------------------------------------------------------------
#  TRIM AND WRITE 
#--------------------------------------------------------------------
# Make use of the lists to perform the trimming 
# Append the trimmed sequences to new lists.
if search_gz is None and search_phred33 is not None:
    
    trimmed_dna_seq =[]
    trimmed_qual_seq = []
    
    for d in dna_sequences:
        trimmed_dna_seq.append(d[three_prime_end:-five_prime_end])
        
    for q in quality_sequence:
        trimmed_qual_seq.append(q[three_prime_end:-five_prime_end])
        
        # The trims are specified in the console by the user to trim each end. 
    console_reads = list(zip(header_list, trimmed_dna_seq, plus_signs, trimmed_qual_seq))
    #--------------------------------------------------------------------
    # Write to an uncompressed fastq file.
    # specification to write:
    
    if search_gz is None and search_phred33 is not None:
        outfile = open("trimmed_uncompressed_33.fq", "w")
        
        print(bright_purple + "\nUNCOMPRESSED PHRED 33" + reset)
        for h,s,p,q in console_reads:
            #print to screen. 
            print(h+"\n"+s+"\n"+p+q)
            
            # print to file.
            print(h+"\n"+s+"\n"+p+q, file = outfile)
            
        outfile.close()
        print(bright_purple + "The trimmed uncompressed PHRED33 FASTQ file is on your computer." + reset)
        print(bright_purple + str(int(line_number/4)), "reads were trimmed." + reset)
        
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        
        
        #--------------------------------------------------------------------
        # Trim each read from the 3' end based on quality (not based on input).   
        # *** Add in code here for ((moving window part)).***
        #--------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        
        
        #--------------------------------------------------------------------
print(bright_cyan + "#--------------------------------------------------------------------#")
print(" STAGE 2")
print(bright_cyan + "#--------------------------------------------------------------------#" + reset)
#-------------------------------------------------------------------
# Decode the lines from bytes to ASCII characters in a compressed file

stage_two_line = 0
second_phred_line = ""
second_seq_line = ""

if search_gz is not None:
    try: 
        infile_two = gzip.open(args.filename.name, "r")
    except FileNotFoundError as error_F:
        print(bright_cyan + "There is a minor problem with opening the FASTQ file.", reset)
        sys.exit(1)
    except gzip.BadGzipFile as error_geez:
        print(bright_cyan + "This file has the suffix of .gz, but it is actually not compressed."+ reset)
        print(bright_cyan+ "Rename the file by removing the .gz part. "+ reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:", str(args.filename.name))
        
        
if search_gz is not None: 
    
    for line in infile_two:
        stage_two_line += 1
        
        if stage_two_line == 2:
            second_seq_line = line.decode('ASCII')
            second_seq_line = second_seq_line.strip()
            print("The sequence line:", second_seq_line)
            
        if stage_two_line == 4:
            second_phred_line = line.decode('ASCII')
            second_phred_line = second_phred_line.strip()
            print(second_phred_line)
    
    infile_two.close()

#------------------------------------------------------------------
#--------------------------------------------------------------------
if search_gz is not None:
    print("Currently", str(int(stage_two_line/4)), "read(s) in this file.")
    print("The length of the untrimmed sequence is:", len(second_seq_line))
    print("The phred line is:", second_phred_line)
    
    search_phred33 = re.search(r'([#=!?>;$&%+()/.*":<]|[0-9])', second_phred_line)
    
    if search_phred33 is not None:
        #print(search_phred33.group(1))
        print("This FASTQ file called has PHRED 33 encoding (SANGER)")
        
    search_phred64 = re.search(r'([L-Z]|[a-j])', second_phred_line)
    if search_phred64 is not None:
        #print(search_phred64.group(1))
        print("This FASTQ file has PHRED 64 encoding (OLD ILLUMINA)")
        
        #--------------------------------------------------------------------
        #--------------------------------------------------------------------
        # phred 33 and compressed.
        # Read each line of the file as a string into a list. 
        # Decode the lines from bytes
        # Ever fourth line is the quality line. 
        
header_list = []
dna_sequences = []
plus_signs = []
quality_sequence=[]

new_line_number = 0

if search_gz is not None and search_phred33 is not None:
    try: 
        infile_two = gzip.open(args.filename.name, "r")
    except FileNotFoundError as error_F:
        print(bright_cyan + "There is a minor problem with opening the FASTQ file.", reset)
        
    for line in infile_two:
        new_line_number += 1
        line = line.decode('ASCII')
        
        if line.startswith("@HWI"):
            header_list.append("")
            header_list[-1] += line.strip("\n")
            dna_sequences.append("")
            
        search_dna = re.search(r'([ACTGN]{101})', line)
        
        if search_dna is not None:
            dna_sequences[-1] += search_dna[1]
            plus_signs.append("")
            
        search_plus = re.search(r'^\+\s', line)
        
        if search_plus is not None:
            plus_signs[-1] += line
            quality_sequence.append('')
            
        search_T = re.search(r'(T)', line)
        search_quality = re.search(r'(.{101})', line, flags = re.ASCII)
        
        if new_line_number % 4 == 0:
            if search_T is None:
                if search_quality is not None:
                    quality_sequence[-1] += line
                    
    infile_two.close()

#--------------------------------------------------------------------
#  TRIM AND WRITE 
#--------------------------------------------------------------------
# Make use of the lists to perform the trimming(zipping is slipping)

second_trimmed_dna_seq = []
second_trimmed_qual_seq = []

if search_gz is not None and search_phred33 is not None:
    for bases in dna_sequences:
        second_trimmed_dna_seq.append(bases[three_prime_end:-five_prime_end])
        
    for qua in quality_sequence:
        second_trimmed_qual_seq.append(qua[three_prime_end:-five_prime_end])
        
    console_reads = list(zip(header_list, second_trimmed_dna_seq, plus_signs, second_trimmed_qual_seq))
    
    #--------------------------------------------------------------------
    # Write to an uncompressed fastq file.
    # specification to write:
    
    if search_gz is not None and search_phred33 is not None:
        outfile = open("trimmed_compressed_33.fq.gz", "w")
        
        print(bright_cyan + "\nCOMPRESSED PHRED 33" + reset)
        for h,s,p,q in console_reads:
            #print to screen. 
            print(h+"\n"+s+"\n"+p+q)
            
            # print to file.
            print(h+"\n"+s+"\n"+p+q, file = outfile)
            
        outfile.close()
        print(bright_cyan + "The trimmed compressed PHRED33 FASTQ file is on your computer." + reset)
        print(bright_cyan + str(int(new_line_number/4)), "reads were trimmed." + reset)
        
        #--------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        
        '''
'''
        
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        #-------------------------------------------------------------------
        
# DEBUG FROM HERE...
#--------------------------------------------------------------------
print(bright_yellow + "#--------------------------------------------------------------------#")
print(" STAGE 3")
print(bright_yellow + "#--------------------------------------------------------------------#" + reset)
# Decode the lines from bytes to ASCII characters in a compressed file
# The encoding has already been identified. 

th_header_list = []
th_dna_sequences = []
th_plus_signs = []
th_quality_sequence =[]

th_line_number = 0

if search_gz is not None and search_phred64 is not None:
    try:
        infile_three = gzip.open(args.filename.name, "r") 
    except FileNotFoundError as error_F:
        print(bright_cyan + "There is a minor problem with opening the file.", reset)
        sys.exit(1)
    except gzip.BadGzipFile as error_geez:
        print(bright_cyan + "This file has the suffix of .gz, but it is actually not compressed."+ reset)
        print(bright_cyan+ "Rename the file by removing the .gz part. "+ reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:", str(args.filename.name))
        
    for line in infile_three:
        line = line.decode('ASCII')
        
        # Find how many lines are in the total file. 
        th_line_number += 1
        if line.startswith("@ILLUMINA-"):
            th_header_list.append("")
            th_header_list[-1] += line.strip()
            th_dna_sequences.append("")
            
        search_dna = re.search(r'([ACTGN]{,101})', line)
        
        if search_dna is not None:
            th_dna_sequences[-1] += search_dna[1]
            th_plus_signs.append("")
            
        search_plus = re.search(r'(^\+ILLUMINA-)', line)
        
        if search_plus is not None:
            th_plus_signs[-1] += line
            th_quality_sequence.append("")
            
        if th_line_number % 4 == 0:
            th_quality_sequence[-1] += line
            
    infile_three.close()
#--------------------------------------------------------------------
#  TRIM AND WRITE 
#--------------------------------------------------------------------
# Make use of the lists to perform the trimming here based on input from the command line. 
th_trimmed_dna_seq = []
th_trimmed_qual_seq = []

if search_gz is not None and search_phred64 is not None:
    for d in th_dna_sequences:
        th_trimmed_dna_seq.append(d[three_prime_end:-five_prime_end])
        
    for q in th_quality_sequence:
        th_trimmed_qual_seq.append(q[three_prime_end:-five_prime_end])
        
    console_reads = list(zip(th_header_list, th_trimmed_dna_seq, th_plus_signs, th_trimmed_qual_seq))
    
    #--------------------------------------------------------------------
    # Write to an compressed fastq file.
    # specification to write:
    
    if search_gz is not None and search_phred64 is not None:
        outfile_three = open("trimmed_compressed_64.fq.gz", "w")
        
        print(bright_yellow + "\nCOMPRESSED PHRED 64" + reset)
        
        for h,d,p,q in console_reads:
            print(h.strip() + "\n" + d + "\n" + p +"\n"+q)
            print(h.strip() + "\n" + d + "\n" + p +"\n"+q, file = outfile_three)
        
        outfile_three.close()
        print(bright_yellow + "The trimmed compressed PHRED64 FASTQ file is on your computer." + reset)
    
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
    
    
    
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
    
#--------------------------------------------------------------------
print(bright_green + "#--------------------------------------------------------------------#")
print("STAGE 4")
print(bright_green + "#--------------------------------------------------------------------#" + reset)
#--------------------------------------------------------------------
# Decode the lines from bytes to ASCII characters in an uncompressed file: phred 64 uncompressed
# The encoding has already been identified above. 


file_there = os.path.isfile('trimmed_compressed_64.fq.gz')


f_header_list= []
f_dna_sequences= []
f_plus_signs= []
f_quality_sequence=[]

f_line_number= 0

if search_gz is None and search_phred64 is not None and file_there is False:
    try:
        infile_four = open(args.filename.name, "r") 
    except FileNotFoundError as error_F:
        print(bright_cyan + "There is a minor problem with opening the file.", reset)
        sys.exit(1)
        
    for line in infile_four:
        f_line_number += 1
        
        if line.startswith("@ILLUMINA-"):
            f_header_list.append(line)
            
            f_dna_sequences.append("")
            
        search_dna = re.search(r'([ACTGN]{90,101})', line)
        
        if search_dna is not None:
            f_dna_sequences[-1] += line
            
            f_plus_signs.append("")
            
        search_plus = re.search(r'^\+ILLUMINA-', line)
        
        if search_plus is not None:
            f_plus_signs[-1] += line
            
            f_quality_sequence.append("")
            
        if f_line_number % 4 == 0:
            f_quality_sequence[-1] += line.strip("\n")
                
    infile_four.close()

#--------------------------------------------------------------------
f_trimmed_dna_seq = []
f_trimmed_qual_seq = []

if search_gz is None and search_phred64 is not None:
    
    for d in f_dna_sequences:
        f_trimmed_dna_seq.append(d[three_prime_end:-five_prime_end])
        
    for q in f_quality_sequence:
        f_trimmed_qual_seq.append(q[three_prime_end:-five_prime_end])
    
    console_reads = list(zip(f_header_list,f_trimmed_dna_seq, f_plus_signs, f_trimmed_qual_seq))
    #--------------------------------------------------------------------
    # Write to an uncompressed fastq file.
    # Make a specification for writing the file. 
    
    if search_gz is None and search_phred64 is not None:
        print(bright_green + "\nUNCOMPRESSED PHRED 64" + reset)
        
        outfile = open("trimmed_uncompressed_64.fq", "w")
    
        for t,d,p,q in console_reads:
            print(t.strip() +"\n"+ d+"\n"+ p + q)
            print(t.strip() +"\n"+ d+"\n"+ p + q, file = outfile)
            
        outfile.close()
        print(bright_green + "The trimmed uncompressed PHRED64 FASTQ file is on your computer" + reset)

#--------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
    


#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#-------------------------------------------------------------------
