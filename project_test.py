#!/usr/bin/env python3 
import gzip 
import re 
import sys
import argparse


bright_green = "\033[0;32m"
bright_cyan = "\033[0;96m"
bright_yellow = "\033[0;33m"
brightish = "\033[0;36m"
bright_purple ="\033[0;35m" 
reset = "\033[0m\n"


print(bright_green + "#--------------------------------------------------------------------#")
print(" Read Trimmer for Next Generation Sequencing Data")
print(bright_green + "#--------------------------------------------------------------------#" + reset)




#--------------------------------------------------------------------
log_file = open("log.dat", "w")
#--------------------------------------------------------------------




#--------------------------------------------------------------------
# Does the filename end with .gz? 
search_gz = re.search(r'.gz$', sys.argv[1])

if search_gz is not None:
    print("This is a gzipped, compressed FASTQ file")
else:
    print("The FASTQ file is uncompressed.")
#--------------------------------------------------------------------
    




#-------------------------------------------------------------------
# uncompressed_file
# If the filename does not end in gz, open the file. 
if search_gz is None:
    try: 
        infile= open(sys.argv[1], "r") 
    except FileNotFoundError as error_F:
        print(bright_green + "There is a minor problem with opening the FASTQ file.", reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:", str(sys.argv[1]))
#-------------------------------------------------------------------
        
        
        
        
        
        
#-------------------------------------------------------------------
# THIS VALUE MUST BE MORE THAN 0 or it will not work!
        
# Set the 3' and 5' ends to trim. 
three_prime_end = int(sys.argv[2])


five_prime_end = int(sys.argv[3])


#-------------------------------------------------------------------
                
        
        
        
        
#--------------------------------------------------------------------
# compressed_file
# If the filename ends in gz, open the file using gzip. 
        
if search_gz is not None:
    try: 
        infile = gzip.open(sys.argv[1], "r")
    except FileNotFoundError as error_F:
            print(bright_cyan + "There is a minor problem with opening the FASTQ file."+ reset)
            sys.exit(1)
    except gzip.BadGzipFile as error_geez:
        print(bright_cyan + "This file has the suffix of .gz, but it is actually not compressed."+ reset)
        print(bright_cyan+ "Rename the file by removing the .gz part. "+ reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:", str(sys.argv[1]))
#-------------------------------------------------------------------



print(bright_purple + "#--------------------------------------------------------------------#")
print(" STAGE 1")
print(bright_purple + "#--------------------------------------------------------------------#" + reset)

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
#------------------------------------------------------------------
if search_gz is None:
    print("There are", str(int(line_number/4)), "reads in this file.")
    print("The length of the untrimmed sequence is:", len(seq_line))
    print("The phred line is:", phred_line)
    
    # search for characters only in phred33 
    search_phred33 = re.search(r'([#=!?>;$&%+()/.*":<]|[0-9])', phred_line)
        
    if search_phred33 is not None:
        print(search_phred33.group(1))
        print("This FASTQ file called has PHRED 33 encoding (SANGER)")
            
    # search for characters only in phred64
    search_phred64 = re.search(r'([L-Z]|[a-j])', phred_line)
    if search_phred64 is not None:
        print(search_phred64.group(1))
        print("This FASTQ file has PHRED 64 encoding (OLD ILLUMINA)")
        
#--------------------------------------------------------------------

# I now now that it is phred 33 and uncompressed.
# read each line of the file as a string into a list. 
# if line number is divisible by 4, it is the quality line. 
    
header_list= []
dna_sequences= []
plus_signs= []
quality_sequence= []

line_number= 0

if search_gz is None and search_phred33 is not None:
    infile = open(sys.argv[1], "r") 
    for line in infile:
        line_number += 1
        
        if line.startswith("@HWI"):
            header_list.append(line.strip("\n"))
            dna_sequences.append("")
            
        search_dna = re.search(r'([ACTGN]{101})', line)
        
        if search_dna is not None:
            dna_sequences[-1] += search_dna[1].strip("\n")
            plus_signs.append("")
        
        search_plus = re.search(r'^\+\s', line)
        
        if search_plus is not None:
            plus_signs[-1] += line.strip("\n")
            quality_sequence.append("")
            
        if line_number % 4 == 0:
            quality_sequence[-1] += line.strip("\n")
        
infile.close()

#--------------------------------------------------------------------
# Make use of the lists to perform the trimming 
# based on input from the command line. 
# Append the trimmed sequences to new lists.


trimmed_dna_seq =[]
trimmed_qual_seq = []

#for h in header_list:
#  print(h, end = "\n")
    
for d in dna_sequences:
    trimmed_dna_seq.append(d[three_prime_end:-five_prime_end])
    
#for p in plus_signs:
#    print(p, end = "\n")
    
for q in quality_sequence:
    trimmed_qual_seq.append(q[three_prime_end:-five_prime_end])
        
#--------------------------------------------------------------------
# Write to an uncompressed fastq file.
reads = zip(header_list, trimmed_dna_seq, plus_signs, trimmed_qual_seq)

# specification to write:
if search_gz is None and search_phred33 is not None:
    outfile = open("trimmed_uncompressed_33.fq", "w")
        
    for heading, seq, plus, quality in reads:
        #print(heading.strip()+"\n"+seq+"\n"+ plus +"\n"+quality)
        print(heading.strip()+"\n"+seq+"\n"+ plus +"\n"+quality, file = outfile)
    
    outfile.close()
    print(bright_purple + "The trimmed uncompressed PHRED33 FASTQ file is on your computer." + reset)

#--------------------------------------------------------------------
        
        
        
        


#--------------------------------------------------------------------
print(bright_purple + "#--------------------------------------------------------------------#")
print(" STAGE 2")
print(bright_purple + "#--------------------------------------------------------------------#" + reset)

# Decode the lines from bytes to ASCII characters in a compressed file
# Identify the phred encoding.
line_number = 0
phred_line = ""
seq_line = ""

if search_gz is not None:
    try: 
        infile = gzip.open(sys.argv[1], "r")
    except FileNotFoundError as error_F:
            print(bright_cyan + "There is a minor problem with opening the FASTQ file.", reset)
            sys.exit(1)
    except gzip.BadGzipFile as error_geez:
        print(bright_cyan + "This file has the suffix of .gz, but it is actually not compressed."+ reset)
        print(bright_cyan+ "Rename the file by removing the .gz part. "+ reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:", str(sys.argv[1]))
    
    
    for line in infile:
        line = line.decode('ASCII')
        
        line_number += 1
        
        if line_number == 2:
            seq_line += str(line).strip("\n")
            print(seq_line)
            
        if line_number == 4:
            phred_line += str(line).strip("\n")
            print(phred_line)
        
infile.close()

if search_gz is not None:
    print("There are", str(int(line_number/4)), "reads in this file.")
    print("The length of the untrimmed sequence is:", len(seq_line))
    print("The phred line is:", phred_line)
    
    search_phred33 = re.search(r'([#=!?>;$&%+()/.*":<]|[0-9])', phred_line)
    
    if search_phred33 is not None:
        print(search_phred33.group(1))
        print("This FASTQ file called has PHRED 33 encoding (SANGER)")
            
    search_phred64 = re.search(r'([L-Z]|[a-j])', phred_line)
    if search_phred64 is not None:
        print(search_phred64.group(1))
        print("This FASTQ file has PHRED 64 encoding (OLD ILLUMINA)")
#--------------------------------------------------------------------
    
# I now know that it is phred 33 and compressed.
# Read each line of the file as a string into a list. 
# Decode the lines from bytes to ASCII characters. 
# if line number is divisible by 4, it is the quality line. 

header_list = []
dna_sequences= []
plus_signs = []
quality_sequence=[]

line_number= 0

if search_gz is not None and search_phred33 is not None:
    infile = gzip.open(sys.argv[1], "r") 
    
    for line in infile:
        line = line.decode('ASCII')
        
        line_number += 1
        
        if line.startswith("@HWI"):
            header_list.append(line.strip("\n"))
            dna_sequences.append("")
    
        search_dna = re.search(r'([ACTGN]{101})', line)
        
        if search_dna is not None:
            dna_sequences[-1] += search_dna[1].strip("\n")
            plus_signs.append("")
            
        search_plus = re.search(r'^\+\s', line)
        
        if search_plus is not None:
            plus_signs[-1] += line.strip("\n")
            quality_sequence.append("")
            
        if line_number % 4 == 0:
            quality_sequence[-1] += line.strip("\n")
            
infile.close()
#--------------------------------------------------------------------
# Make use of the lists to perform the trimming
# based on input from the command line



trimmed_dna_seq =[]
trimmed_qual_seq = []

for d in dna_sequences:
    trimmed_dna_seq.append(d[three_prime_end:-five_prime_end])

for q in quality_sequence:
    trimmed_qual_seq.append(q[three_prime_end:-five_prime_end])
    
#--------------------------------------------------------------------
# Write to a compressed.gz fastq file.
# Make a specification for writing the file. 
    
reads = zip(header_list, trimmed_dna_seq, plus_signs, trimmed_qual_seq)

# specification to write:
if search_gz is not None and search_phred33 is not None:
    
    outfile = open("trimmed_compressed_33.fq.gz", "w")
            
    for heading, seq, plus, quality in reads:
        print(heading.strip()+"\n"+seq+"\n"+ plus +"\n"+quality, file = outfile)
        
    outfile.close()
    print(bright_purple + "The trimmed compressed PHRED33 FASTQ file is on your computer" + reset)
#--------------------------------------------------------------------

    





#--------------------------------------------------------------------
print(bright_purple + "#--------------------------------------------------------------------#")
print(" STAGE 3")
print(bright_purple + "#--------------------------------------------------------------------#" + reset)
# Decode the lines from bytes to ASCII characters in a compressed file
# The encoding has already been identified. 

header_list= []
dna_sequences= []
plus_signs= []
quality_sequence=[]

line_number = 0

if search_gz is not None and search_phred64 is not None:
    try:
        infile = gzip.open(sys.argv[1], "r") 
    except FileNotFoundError as error_F:
        print(bright_cyan + "There is a minor problem with opening the file.", reset)
        sys.exit(1)
    except gzip.BadGzipFile as error_geez:
        print(bright_cyan + "This file has the suffix of .gz, but it is actually not compressed."+ reset)
        print(bright_cyan+ "Rename the file by removing the .gz part. "+ reset)
        sys.exit(1)
    finally:
        print("FASTQ File Information: \nThe name of the file is:", str(sys.argv[1]))
        
        for line in infile:
            line = line.decode('ASCII')
            
            line_number += 1
            
            if line.startswith("@ILLUMINA-"):
                header_list.append(line.strip("\n"))
                dna_sequences.append("")
                
            search_dna = re.search(r'([ACTGN]{101})', line)
            
            if search_dna is not None:
                dna_sequences[-1] += search_dna[1].strip("\n")
                plus_signs.append("")
                
            search_plus = re.search(r'^\+ILLUMINA-', line)
            
            if search_plus is not None:
                plus_signs[-1] += line.strip("\n")
                quality_sequence.append("")
                
            if line_number % 4 == 0:
                quality_sequence[-1] += line.strip("\n")
            
            #****THIS NEEDS TO BE REMOVED*******
            if line_number == 40:
                break
        
infile.close()

#--------------------------------------------------------------------
# Make use of the lists to perform the trimming here based on input from the command line. 

trimmed_dna_seq = []
trimmed_qual_seq = []

for d in dna_sequences:
    trimmed_dna_seq.append(d[three_prime_end:-five_prime_end])

for q in quality_sequence:
    trimmed_qual_seq.append(q[three_prime_end:-five_prime_end])
    
#--------------------------------------------------------------------
# Write to a compressed.gz fastq file.
# Make a specification for writing the file. 
reads_64 = zip(header_list, trimmed_dna_seq, plus_signs, trimmed_qual_seq)

# Specification to write:
if search_gz is not None and search_phred64 is not None:
    outfile = open("trimmed_compressed_phred64_fastq.gz", "w")
    
    for hea, siq, plud, qwal in reads_64:
        print(hea.strip() + "\n" +
              siq + "\n" +
              plud +"\n"+ 
              qwal, file = outfile)
    
    outfile.close()
    
    print(bright_purple + "The trimmed compressed PHRED64 FASTQ file is on your computer" + reset)
#--------------------------------------------------------------------





#--------------------------------------------------------------------
print(bright_purple + "#--------------------------------------------------------------------#")
print("STAGE 4")
print(bright_purple + "#--------------------------------------------------------------------#" + reset)
#--------------------------------------------------------------------
# Decode the lines from bytes to ASCII characters in an uncompressed file: phred 64 uncompressed
# The encoding has already been identified above. 

header_list= []
dna_sequences= []
plus_signs= []
quality_sequence=[]

line_number= 0


if search_gz is None:
    if search_phred64 is not None:
        try:
            infile = open(sys.argv[1], "r") 
        except FileNotFoundError as error_F:
            print(bright_cyan + "There is a minor problem with opening the file.", reset)
            sys.exit(1)
                
        for line in infile:
            line_number += 1
            
            if line.startswith("@ILLUMINA-"):
                header_list.append(line.strip("\n"))
                dna_sequences.append("")
            
            search_dna = re.search(r'([ACTGN]{90,101})', line)
            
            if search_dna is not None:
                dna_sequences[-1] += line.strip("\n")
                plus_signs.append("")
            
            search_plus = re.search(r'^\+ILLUMINA-', line)
            
            if search_plus is not None:
                plus_signs[-1] += line.strip("\n")
                quality_sequence.append("")
            
            if line_number % 4 == 0:
                quality_sequence[-1] += line.strip("\n")
            
infile.close()

#--------------------------------------------------------------------
# Make use of the lists to perform the trimming here based on input from the command line. 
trimmed_dna_seq = []
trimmed_qual_seq = []

for d in dna_sequences:
    trimmed_dna_seq.append(d[three_prime_end:-five_prime_end])
    
for q in quality_sequence:
    trimmed_qual_seq.append(q[three_prime_end:-five_prime_end])
#--------------------------------------------------------------------
# Write to an uncompressed fastq file.
# Make a specification for writing the file. 
    
reads_comp64 = zip(header_list, trimmed_dna_seq, plus_signs, trimmed_qual_seq)

# specification to write:
if search_gz is None and search_phred64 is not None:
    outfile = open("trimmed_fastq64_uncomp.fq", "w")
    
    for t,d,p,q in reads_comp64:
        print(t.strip()+"\n"+
              d + "\n"+ 
              p + "\n"+
              q, file = outfile)
    
    outfile.close()
    
    print(bright_purple + "The trimmed uncompressed PHRED64 FASTQ file is on your computer" + reset)
#--------------------------------------------------------------------
#turn the line into unicode. (* this trims a sequence line based on the value of the character)
#for letter in line:
#if ord(letter) - 33 > 20:
#print(ord(letter) - 33, end  = ",")
