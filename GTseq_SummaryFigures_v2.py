#!/usr/bin/python
# GTseq_SummaryFigures_v2.py
# by Nate Campbell
# produce summary figures for GTseq libraries using GTseq_Genotyper_v2 output formatted files.
# Requires matplotlib module for plotting summary figures.
# Also outputs summary data in text format for further analysis.

import math
from os import listdir
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

print('type the path to directory containing .genos files for library *use single quotes*\nFormat= \'/home/user/...\'')
path = input()

print('type the library name *use single quotes*')
Lname = input()

fname = Lname + '.pdf'
print(fname)

fout_name = Lname + '_SummaryData.txt'
f_out = open(fout_name, 'w')
f_out.write(Lname + '\t')
f_out.write('GTseq Summary Data\n\n')

pp=PdfPages(fname)

list = listdir( path )
flist = []
assaylist = []

#filter for .genos files in directory and add to "flist"...
for i in list:
 if '.genos' in i:
  j = path + '/' + i
  flist.append(j)

#open top file and create assay list...initialize dictionary of loci with their corresponding percentage of on-target reads
OT_Dict = {}
StDEV_Dict = {}
OTP_Dict = {}
StDEV2_Dict = {}
f = open(flist[0])
lineNo = 0
for line in f:
 lineNo = lineNo + 1
 if lineNo > 1:
  stuff = line.split(',')
  assaylist.append(stuff[0])
  OT_Dict[stuff[0]] = float(0)
  StDEV_Dict[stuff[0]] = float(0)
  OTP_Dict[stuff[0]] = float(0)
  StDEV2_Dict[stuff[0]] = float(0)

f.close()

#Initialize variables...
AssayNum = float(len(assaylist))
xmax = 0
num90 = float(0)
per90 = float(0)
inds = float(len(flist))
aveOTP = float(0)

AssayNum = float(len(assaylist))
inds = float(len(flist))

for loci in assaylist:
 for genos in flist:
  g = open(genos)
  for line in g:
   if loci in line:
    info = line.split(',')
    OT_Dict[info[0]] = OT_Dict[info[0]] + float(info[10])
    OTP_Dict[info[0]] = OTP_Dict[info[0]] + float(info[9])

#Get the mean of percenage of OT reads for each locus and mean of percentage
#forward primer reads containing probe sequences for each locus...
for loci in assaylist:
 OT_Dict[loci] = OT_Dict[loci] / inds
 OTP_Dict[loci] = OTP_Dict[loci] / inds

#Calculate standard deviations at each locus...
for loci in assaylist:
 for genos in flist:
  g = open(genos)
  for line in g:
   if loci in line:
    info = line.split(',')
    variance = (float(info[10]) - OT_Dict[loci])**2
    variance2 = (float(info[9]) - OTP_Dict[loci])**2
    StDEV_Dict[loci] = StDEV_Dict[loci] + variance
    StDEV2_Dict[loci] = StDEV2_Dict[loci] + variance2
    

for loci in assaylist:
 StDEV_Dict[loci] = (math.sqrt(StDEV_Dict[loci] / inds))
 StDEV2_Dict[loci] = (math.sqrt(StDEV2_Dict[loci] / inds))

#plot read distribution bar graph using means and standard deviation for each locus...
bar_width = 1
opacity = 0.4
left = 0
L_AvOTP = 100 / AssayNum

print('Read Distribution data (sorted by locus name)\n')
f_out.write('Read Distribution data (sorted by locus name)\n')
for loci in assaylist:
 plt.bar(left, OT_Dict[loci], width=bar_width, bottom=None, hold=None, 
 alpha=opacity, yerr=StDEV_Dict[loci], ecolor='black',
 color='b', capsize=1)
 print(loci, OT_Dict[loci], StDEV_Dict[loci])
 f_out.write(loci + '\t' + str(OT_Dict[loci]) + '\t' + str(StDEV_Dict[loci]) + '\n')
 left = left + bar_width

plt.plot([0, left], [L_AvOTP, L_AvOTP], 'r-', linewidth=1.0)
plt.xlabel('Loci')
plt.ylabel('Percentage of On-Target Reads')
plt.title(Lname + ' : Read distribution among loci (unsorted)')
plt.tight_layout()
plt.savefig(pp, format='pdf')
plt.clf()

Sorted_OT = sorted(OT_Dict.values())
Sorted_OTkeys = sorted(OT_Dict, key=OT_Dict.get)
Sorted_stDEV = []
Sorted_OTP = sorted(OTP_Dict.values())
Sorted_OTPkeys = sorted(OTP_Dict, key=OTP_Dict.get)
Sorted_stDEV2 = []

#Get sorted list of standard deviations for error bars...
for loci in Sorted_OTkeys:
 Sorted_stDEV.append(StDEV_Dict[loci])
 Sorted_stDEV2.append(StDEV2_Dict[loci])

print('Read Distribution data (sorted by value)\n')
f_out.write('\nRead Distribution data (sorted by value)\n')
#populate read distribution bar graph using sorted data...
left = 0
Assays = int(AssayNum)
for x in range(0, Assays):
 plt.bar(left, Sorted_OT[x], width=bar_width, bottom=None, hold=None, 
 alpha=opacity, yerr=Sorted_stDEV[x], ecolor='black',
 color='green', capsize=1)
 print(Sorted_OTkeys[x], Sorted_OT[x], Sorted_stDEV[x])
 f_out.write(Sorted_OTkeys[x] + '\t' + str(Sorted_OT[x]) + '\t' + str(Sorted_stDEV[x]) + '\n')
 left = left + bar_width

plt.plot([0, left], [L_AvOTP, L_AvOTP], 'r-', linewidth=1.0)
plt.xlabel('Loci')
plt.ylabel('Percentage of On-Target Reads')
plt.title(Lname + ' : Read distribution among loci (Sorted)')
plt.tight_layout()
plt.savefig(pp, format='pdf')
plt.clf()

print('Primer Reads On-Target (reads with forward primer AND probe / reads with fwd primer)*100\n')
f_out.write('\nPrimer Reads On-Target (reads with forward primer AND probe / reads with fwd primer)*100\n')
#Create bar graph of percentage (reads with forward primer AND probe / reads with fwd primer)...
left = 0
for loci in assaylist:
 plt.bar(left, OTP_Dict[loci], width=bar_width, bottom=None, hold=None, 
 alpha=opacity, yerr=StDEV2_Dict[loci], ecolor='black',
 color='b', capsize=1)
 print(loci, OTP_Dict[loci], StDEV2_Dict[loci])
 f_out.write(loci + '\t' + str(OTP_Dict[loci]) + '\t' + str(StDEV2_Dict[loci]) + '\n')
 left = left + bar_width

plt.xlabel('Loci')
plt.ylabel('Percentage On-Target Primers')
plt.title(Lname + ' : Primers On-Target (unsorted)')
plt.tight_layout()
plt.savefig(pp, format='pdf')
plt.clf()

print('Primer Reads On-Target (sorted) (reads with forward primer AND probe / reads with fwd primer)*100\n')
f_out.write('\nPrimer Reads On-Target (sorted) (reads with forward primer AND probe / reads with fwd primer)*100\n')
#Create sorted bar graph of percentage (reads with forward primer AND probe / reads with fwd primer)...
left = 0
for x in range(0, Assays):
 plt.bar(left, Sorted_OTP[x], width=bar_width, bottom=None, hold=None, 
 alpha=opacity, yerr=Sorted_stDEV2[x], ecolor='black',
 color='red', capsize=1)
 print(Sorted_OTPkeys[x], Sorted_OTP[x], Sorted_stDEV2[x])
 f_out.write(Sorted_OTPkeys[x] + '\t' + str(Sorted_OTP[x]) + '\t' + str(Sorted_stDEV2[x]) + '\n')
 left = left + bar_width

plt.xlabel('Loci')
plt.ylabel('Percentage On-Target Primers')
plt.title(Lname + ' : Primers On-Target (sorted)')
plt.tight_layout()
plt.savefig(pp, format='pdf')
plt.clf()

print('Library Summary\n')
f_out.write('\nLibrary Summary\n')
#open each .genos file and define x and y coordinates for plotting Raw-reads vs. GT% graph...

for genos in flist:
 g = open(genos)
 genper = float(0)
 otreads = 0
 rawreads = 0
 NA = 0
 lineNo2 = 0
 for line in g:
  lineNo2 = lineNo2 + 1
  if lineNo2 == 1:
   info = line.split(',')
   RRarr = info[1].split(':')
   OTarr = info[2].split(':')
   OTParr = info[3].split(':')
   otreads = int(OTarr[1])
   rawreads = int(RRarr[1])
   aveOTP = aveOTP + float(OTParr[1])
  elif lineNo2 > 1:
   info2 = line.split(',')
   if 'NA' in info2[5]:
    NA = NA + 1
 if rawreads > xmax:
  xmax = rawreads
 NA = float(NA)
 genper = (1 - (NA / AssayNum)) * 100
 scale = 100.0
 if genper >= 90:
  num90 = num90 + 1
 #print(genos,otreads,genper)
 plt.scatter(rawreads, genper, c='orange', s=scale, label='black',alpha=0.4, edgecolors='none')
per90 = num90 / inds * 100
aveOTP = aveOTP / inds

perSTR = str(per90)
aveOTPSTR = str(aveOTP)

text1 = 'Samples in library: ' + str(int(inds))
text2 = 'Samples over 90% GT: ' + str(int(num90))
text3 = 'Percentage over 90% GT: ' + perSTR[:4] + '%'
text4 = 'Average OT-Percentage: ' + aveOTPSTR[:4] + '%'
print(text1, text2, text3, text4)
f_out.write(text1 + '\n')
f_out.write(text2 + '\n')
f_out.write(text3 + '\n')
f_out.write(text4 + '\n')

font = {'family' : 'serif',
 'color'  : 'darkred',
 'weight' : 'normal',
 'size'   : 16,
 }

plt.plot([0, xmax], [90, 90], 'r-', linewidth=2.0)
plt.grid(True)
plt.axis([-5, xmax, -5, 105])
plt.title(Lname)
plt.xlabel('Raw Reads')
plt.ylabel('Genotyping Percentage')
plt.text(25000, 50, text1, fontdict=font)
plt.text(25000, 40, text2, fontdict=font)
plt.text(25000, 30, text3, fontdict=font)
plt.text(25000, 20, text4, fontdict=font)
plt.savefig(pp, format='pdf')
plt.clf()

#Re-initialize variable at zero...

xmax = 0
num90 = float(0)
per90 = float(0)
aveOTP = float(0)

#open each .genos file and define x and y coordinates for plotting OT-reads vs. GT% graph...

for genos in flist:
 g = open(genos)
 genper = float(0)
 otreads = 0
 rawreads = 0
 NA = 0
 lineNo2 = 0
 for line in g:
  lineNo2 = lineNo2 + 1
  if lineNo2 == 1:
   info = line.split(',')
   RRarr = info[1].split(':')
   OTarr = info[2].split(':')
   OTParr = info[3].split(':')
   otreads = int(OTarr[1])
   rawreads = int(RRarr[1])
   aveOTP = aveOTP + float(OTParr[1])
  elif lineNo2 > 1:
   info2 = line.split(',')
   if 'NA' in info2[5]:
    NA = NA + 1
 if otreads > xmax:
  xmax = otreads
 NA = float(NA)
 genper = (1 - (NA / AssayNum)) * 100
 scale = 100.0
 if genper >= 90:
  num90 = num90 + 1
 #print(genos,otreads,genper)
 plt.scatter(otreads, genper, c='black', s=scale, label='black',alpha=0.4, edgecolors='none')
per90 = num90 / inds * 100
aveOTP = aveOTP / inds

perSTR = str(per90)
aveOTPSTR = str(aveOTP)

text1 = 'Samples in library: ' + str(int(inds))
text2 = 'Samples over 90% GT: ' + str(int(num90))
text3 = 'Percentage over 90% GT: ' + perSTR[:4] + '%'
text4 = 'Average OT-Percentage: ' + aveOTPSTR[:4] + '%'

font = {'family' : 'serif',
 'color'  : 'darkred',
 'weight' : 'normal',
 'size'   : 16,
 }

plt.plot([0, xmax], [90, 90], 'r-', linewidth=2.0)
plt.grid(True)
plt.axis([-5, xmax, -5, 105])
plt.title(Lname)
plt.xlabel('On-Target Reads')
plt.ylabel('Genotyping Percentage')
plt.text(25000, 50, text1, fontdict=font)
plt.text(25000, 40, text2, fontdict=font)
plt.text(25000, 30, text3, fontdict=font)
plt.text(25000, 20, text4, fontdict=font)
plt.savefig(pp, format='pdf')
plt.clf()

print('Per-locus Summary\n')
f_out.write('\nPer-locus Summary\n')
#Create XY scatter plot of each locus attempted...
for loci in assaylist:
 xmax = 0
 ymax = 0
 fmax = 0
 A1_corr = float(0)
 A2_corr = float(0)
 gt_per = float(0)
 gt_inds = float(0)
 inds = float(len(flist))
 for genos in flist:
  g = open(genos)
  for line in g:
   if loci in line:
    info = line.split(',')
    size = len(info)
    if size == 11:
     A1_corr = float(info[6])
     A2_corr = float(info[7])
    xarr = info[1].split('=')
    yarr = info[2].split('=')
    x = int(round(float(xarr[1])))
    y = int(round(float(yarr[1])))
    ratio = float(info[3])
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
 A1_corr = str(A1_corr)
 A2_corr = str(A2_corr)
 text = ' : Corrections [' + A1_corr[:3] + ', ' + A2_corr[:3] + '] : GT% = ' + gt_per[:5]
 print(loci + text)
 f_out.write(loci + text + '\n')
 plt.grid(True)
 plt.axis([-5, fmax, -5, fmax])
 plt.plot([0, 10], [10, 0], 'y-', linewidth=2.0)
 plt.plot([8, 10000], [2, 2000], 'r-', linewidth=2.0)
 plt.plot([2000, 2], [10000, 8], 'b-', linewidth=2.0)
 plt.plot([5, 10000], [5, 10000], 'm-', linewidth=2.0)
 plt.plot([6.6, 10000], [3.3, 5000], 'k-', linewidth=2.0)
 plt.plot([3.3, 5000], [6.6, 10000], 'k-', linewidth=2.0)
 plt.title(loci + text)
 plt.xlabel('A1 counts')
 plt.ylabel('A2 counts')
 plt.savefig(pp, format='pdf')
 plt.clf()

pp.close()
f_out.close()
