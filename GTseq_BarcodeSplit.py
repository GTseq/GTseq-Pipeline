#!/usr/bin/python
# GTseq_BarcodeSplit.py
# by Nate Campbell
# Use this script to split out individual files by barcode combination rather than using "grep method".
# Input file 1 is a .csv file containing sample and barcode information [individual sample names, PlateID,i7_name,i7_sequence,i5_name,i5_sequence] for each sample.
# example: Sample,PlateID,i7_name,i7_sequence,i5_name,i5_sequence
#          Sample123,P1234,i001,ACCGTA,25,CCCTAA
#          Sample321,P1234,i001,ACCGTA,26,GGCACA
#          ....
# Note: The header line is ignored while executing the script so always include it to avoid missing data.
# Input file 2 is the .fastq file containing the sequences.
# Output is a set of fastq files for all individuals included in the library named in the format [i7name]_[i5name]_[PlateID]_[SampleID].fastq 
# This script uses only a single processing thread and was clocked at about 1.5 hours to complete barcode splitting for one lane of data.
# Speed is comparable to "grep" method while using fewer compute resources.
# Using this script for barcode splitting and the GTseq_Genotyper_v2.pl script for genotyping allows the GTseq pipeline to be run on a linux
# desktop computer.  (Raw data to genotypes in approximately 3 hours).

import threading

print('type the path to input file\n Use single quotes if using python2 Format= \'/home/user/...\'')
path1 = input()

print('type the path to the fastq file to split\n Use single quotes if using python2 Format= \'/home/user/...\'')
path2 = input()


def split_file(individual_list):
     fq = open(path2, 'r')
     BC_Dict = {}
     File_Dict = {}
     handle_dict = {}

     for lines in individual_list:
          stuff = lines.split(',')
          name = stuff[2] + '_' + stuff[4] + '_' + stuff[1] + '_' + stuff[0] + '.fastq'
          BC_Dict[name] = stuff[3] + stuff[5]
          File_Dict[stuff[3] + stuff[5]] = name
          handle_dict[stuff[3] + stuff[5]] = open(File_Dict[stuff[3] + stuff[5]], 'a')

     lineNo2 = 0
     writelines = 0

     for line in fq:
          lineNo2 = lineNo2 + 1
          if writelines < 5 and writelines > 0:
               f_out.write(line)
               writelines = writelines + 1
               if writelines == 4:
                    writelines = 0

          if '@HISEQ' in line:
               info = line.split(':')
               BC = info[9]
               if BC in File_Dict:
                    f_out = handle_dict[BC]
                    f_out.write(line)
                    writelines = 1
     fq.close()
     return

individuals = 0
start = 1  # skip line 1 of input file (header line)...
end = 0
f = open(path1, 'r')
lineNo = 0

# Determine the number of individuals in the library...
for line in f:
     lineNo = lineNo + 1
file_end = lineNo
individuals = lineNo - 1
if individuals > 500 + start:
     end = 500 + start
else:
     end = file_end

sets = 0
Samples = True
starting_ind = 1

while Samples == True:
     lineNo = 0

     f = open(path1, 'r')
     list1 = []     # set empty list...
     for line in f:
          lineNo = lineNo + 1
          if lineNo > start and lineNo <= end: # Populate list in sets of 500 individuals ...
               list1.append(line)
     sets = sets + 1
     size = len(list1)
     ending_ind = end - 1
     print('Working on individuals %s to %s' % (starting_ind, ending_ind))
     print('Process %s is splitting %s samples' % (sets, size))

     t = threading.Thread(target=split_file(list1))  # process each set of 500 individuals by pushing onto thread...
     t.start()

     if end == file_end:
          Samples = False
          break
     elif start > file_end:
          Samples = False
          break

     start = end
     end = start + 500
     starting_ind = starting_ind + 500

     if end > file_end:
          end = file_end
          ending_ind = file_end - 1
     #print('Next set of individuals begins at line %s and ends at %s of input file' % (start, end))

     f.close()

