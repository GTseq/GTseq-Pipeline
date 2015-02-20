#!/usr/bin/python
# GTseq_GenoCorrect_v2.py
# by Nate Campbell
# Requires the matplotlib python module for graphics output.
# interactive script to correct for psv signal in specified loci...
# This version uses ".genos" output files from GTseq_Genotyper_v2.pl script and modifies allele counts using user input.
# Displays xy cartesian graph of specified loci using allele1 and allele2 correction values.
# Script assumes the presence of a single alternate paralogous sequence variant.  The total read counts for each sample at the 
# given locus are divided by four.  This number is then multiplied by the users correction factor and subtracted from the reads
# for the given allele.  Correction values should be between 0 and 3 and only 1 allele should be corrected at a given locus.
# correction values of 0 make no changes to the original data.


from os import listdir
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys

print('type the path to directory containing .genos files for library *use single quotes*\nFormat= \'/home/user/...\'')
path = input()

print('type the library name *Use single quotes*')
Lname = input()

fname = Lname + '_Corrections.pdf'
print(fname)

pp=PdfPages(fname)

list = listdir( path )
flist = []
assaylist = []

#filter for .genos files in directory and add to "flist"...
for i in list:
  if '.genos' in i:
    j = path + '/' + i
    flist.append(j)

#open top file and create assay list...
f = open(flist[0])
lineNo = 0
for line in f:
  lineNo = lineNo + 1
  if lineNo > 1:
    stuff = line.split(',')
    assaylist.append(stuff[0])

f.close()

Hit = 0
while True:
 print('Enter Assay name for correction *Use single quotes*')
 Locus = input()
 for loci in assaylist:
  if Locus in loci:
   print('is this the assay you want to examine? %s (Y/N) *Use single quotes*' % loci)
   Answer = input()
   if Answer == 'Y':
    Hit = Hit + 1
    Locus = loci
    break
  else:
   pass

#Plot A1 vs A2 counts given correction factor...
 while Hit == 1:
     xmax = 0
     ymax = 0
     fmax = 0
     gt_per = float(0)
     gt_inds = float(0)
     inds = float(len(flist))
     print('Enter correction factor for Allele 1 (enter 0 for no change)')
     A1corr = float(input())
     print('Enter correction factor for Allele 2 (enter 0 for no change)')
     A2corr = float(input())
     print('Working on it...\nKeepin it 100...')
     for genos in flist:
      g = open(genos)
      for line in g:
       if Locus in line:
        info = line.split(',')
        xarr = info[1].split('=')
        yarr = info[2].split('=')
        x = int(round(float(xarr[1])))
        y = int(round(float(yarr[1])))
        sum_xy = x + y
        x = x - (sum_xy / 4 * A1corr)
        if x < 0:
         x = 0
        y = y - (sum_xy / 4 * A2corr)
        if y < 0:
         y = 0
        flx = 0
        if x == 0:
         flx = float(0.01)
        else:
         flx = float(x)
        fly = 0
        if y == 0:
         fly = float(0.01)
        else:
         fly = float(y)
        ratio = flx / fly
        sum_xy = x + y
        scale = 100.0
        if x > xmax:
         xmax = x
        if y > ymax:
         ymax = y
        if sum_xy < 10:
         color = 'yellow'
        elif sum_xy == 0:
         color = 'yellow'
        elif ratio > 5:
         color = 'red'
         gt_inds = gt_inds + 1
        elif ratio < 0.2:
         color = 'blue'
         gt_inds = gt_inds + 1
        elif ratio < 2 and ratio > 0.5:
         color = 'purple'
         gt_inds = gt_inds + 1
        else:
         color = 'yellow'
        plt.scatter(x, y, c=color, s=scale, label=color,alpha=0.4, edgecolors='none')
     if xmax > ymax:
      fmax = xmax
     else:
      fmax = ymax
     gt_per = gt_inds / inds * 100
     gt_per = str(gt_per)
     A1corr = str(A1corr)
     A2corr = str(A2corr)
     text = ': Corrections [' + A1corr[:3] + ' : ' + A2corr[:3] + ']: ' + 'GT% = ' + gt_per[:4]
     plt.grid(True)
     plt.axis([-5, fmax, -5, fmax])
     plt.plot([0, 10], [10, 0], 'y-', linewidth=2.0)
     plt.plot([8, 10000], [2, 2000], 'r-', linewidth=2.0)
     plt.plot([2000, 2], [10000, 8], 'b-', linewidth=2.0)
     plt.plot([5, 10000], [5, 10000], 'm-', linewidth=2.0)
     plt.plot([6.6, 10000], [3.3, 5000], 'k-', linewidth=2.0)
     plt.plot([3.3, 5000], [6.6, 10000], 'k-', linewidth=2.0)
     plt.title(Locus + text)
     plt.xlabel('A1 counts')
     plt.ylabel('A2 counts')

     plt.show(block=False)
	
     print('Does this look better? (Y/N) *Use single quotes* *Y prints current figure to pdf*')
     answer = input()

     if answer == 'Y':
      plt.savefig(pp, format='pdf')
      plt.close()
      break

     elif answer == 'N':
      plt.close()

 if Hit == 0:
  print('Sorry no loci matching that name were found...')
 elif Hit == 1:
  print('All done? (Y/N)')
  ans = str(input())
  if ans == 'Y':
   pp.close()
   break
  else:
   Hit = 0
   pass

