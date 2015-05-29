#!/usr/local/bin/python3.4
# GTseq_BarcodeSplit_MP.py
# by Nate Campbell
# Use this script to split out individual files by barcode combination rather than using "grep method".
# This version will process up to 5,000 individual samples and simultaneously run up to 10 processing cores for faster parsing of individual sequences.
# Input file 1 is a .csv file containing sample and barcode information [individual sample names, PlateID,i7_name,i7_sequence,i5_name,i5_sequence] for each sample.
# example: Sample,PlateID,i7_name,i7_sequence,i5_name,i5_sequence
#          Sample123,P1234,i001,ACCGTA,25,CCCTAA
#          Sample321,P1234,i001,ACCGTA,26,GGCACA
#          ....
# Note: The header line is ignored while executing the script so always include it to avoid missing data.  Also, note that output files are appended
# and therefore not overwritten if the script is called multiple times.  If for some reason the script needs to be run more than once for the same set of samples, 
# the original output files will need to be deleted to avoid having files with duplicated sequences.
# Input file 2 is the .fastq file containing the sequences.
# Output is a set of fastq files for all individuals included in the library named in the format [i7name]_[i5name]_[PlateID]_[SampleID].fastq 
# This script uses multiple processors (1 for every 500 samples) and was clocked at about 20 minutes to complete barcode splitting for one lane of data.
# Speed is faster than grep method while using fewer compute resources.
# Using this script for barcode splitting and the GTseq_Genotyper_v2.pl script for genotyping allows the GTseq pipeline to be run on a linux
# desktop computer.  (Raw data to genotypes in less than 1 hour).

from multiprocessing import Process

print('type the path to input file\nFormat= /home/user/...')
path1 = input()

print('type the path to the fastq file to split\nFormat= /home/user/...')
path2 = input()


def split_file(individual_list):
     if individual_list == 'list1':
          individual_list = list1
     elif individual_list == 'list2':
          individual_list = list2
     elif individual_list == 'list3':
          individual_list = list3
     elif individual_list == 'list4':
          individual_list = list4
     elif individual_list == 'list5':
          individual_list = list5
     elif individual_list == 'list6':
          individual_list = list6
     elif individual_list == 'list7':
          individual_list = list7
     elif individual_list == 'list8':
          individual_list = list8
     elif individual_list == 'list9':
          individual_list = list9
     else:
          individual_list = list10

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

def Main():
     if len(list1) > 0:
          P1 = Process(target=split_file, args=('list1',))
          print('Process initiated for 1st set of samples')
          P1.start()
     if len(list2) > 0:
          P2 = Process(target=split_file, args=('list2',))
          print('Process initiated for 2nd set of samples')
          P2.start()
     if len(list3) > 0:
          P3 = Process(target=split_file, args=('list3',))
          print('Process initiated for 3rd set of samples')
          P3.start()
     if len(list4) > 0:
          P4 = Process(target=split_file, args=('list4',))
          print('Process initiated for 4th set of samples')
          P4.start()
     if len(list5) > 0:
          P5 = Process(target=split_file, args=('list5',))
          print('Process initiated for 5th set of samples')
          P5.start()
     if len(list6) > 0:
          P6 = Process(target=split_file, args=('list6',))
          print('Process initiated for 6th set of samples')
          P6.start()
     if len(list7) > 0:
          P7 = Process(target=split_file, args=('list7',))
          print('Process initiated for 7th set of samples')
          P7.start()
     if len(list8) > 0:
          P8 = Process(target=split_file, args=('list8',))
          print('Process initiated for 8th set of samples')
          P8.start()
     if len(list9) > 0:
          P9 = Process(target=split_file, args=('list9',))
          print('Process initiated for 9th set of samples')
          P9.start()
     if len(list10) > 0:
          P10 = Process(target=split_file, args=('list10',))
          print('Process initiated for 10th set of samples')
          P10.start()
     print('Keepin it 100!  Your files will be ready shortly...')


list1 = []
list2 = []
list3 = []
list4 = []
list5 = []
list6 = []
list7 = []
list8 = []
list9 = []
list10 = []

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
     list0 = []     # set empty list...
     for line in f:
          lineNo = lineNo + 1
          if lineNo > start and lineNo <= end: # Populate list in sets of 500 individuals ...
               list0.append(line)
     sets = sets + 1
     size = len(list0)
     if sets == 1:
          list1 = tuple(list0)
     elif sets == 2:
          list2 = tuple(list0)
     elif sets == 3:
          list3 = tuple(list0)
     elif sets == 4:
          list4 = tuple(list0)
     elif sets == 5:
          list5 = tuple(list0)
     elif sets == 6:
          list6 = tuple(list0)
     elif sets == 7:
          list7 = tuple(list0)
     elif sets == 8:
          list8 = tuple(list0)
     elif sets == 9:
          list9 = tuple(list0)
     else:
          list10 = tuple(list0)

     ending_ind = end - 1

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

if __name__ == '__main__':
     Main()

