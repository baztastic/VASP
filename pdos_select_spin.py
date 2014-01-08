#!/usr/bin/env python

# Based on pdos_select.py by Isaac Tamblyn
# https://github.com/itamblyn/VASP
# This version works for spin-polarised VASP calculations, outputting the up and down states separately

import sys
import commands
import numpy

natom_skip = 0 # in this case, the atoms of interest are at the top of the file
natom = 60      # how many atoms (sequentially) after natom_skip you would like to include

def main():

  try:
    input_filename = sys.argv[1]
  except IndexError:
    print '\nusage: ' + sys.argv[0] + ' input_filename'
    print 'Since no file was provided, I will try DOSCAR in the run dir'
    input_filename = 'DOSCAR'

  core_valence_divide = -1000 # <- all valence 

  command_line_counter = commands.getoutput('wc -l ' + input_filename).split()

  if len(command_line_counter) != 2:
    print 'Error determining file size'
  else:
    number_of_lines = int(command_line_counter[0])

  outputFile_dos = open('dos.dat', 'w')
  outputFile_pdos = open('pdos.dat', 'w')

  inputFile = open(input_filename, 'r')

#  If you wanted to read natom from the file, here it is. For this code though we want to specify
#  so we will simply discard this line
  null = inputFile.readline()
#  natom = int(inputFile.readline().split()[0])

  fermi_subtract = True

  print '!!! natom = ' + str(natom)

  for i in range(4):
      line = inputFile.readline()

  inline_total = inputFile.readline()

  emin_total = float(inline_total.split()[1])
  emax_total = float(inline_total.split()[0])
  enum_total = int(inline_total.split()[2])
  efermi_total = float(inline_total.split()[3])

  dos = numpy.zeros((enum_total, 5), dtype = numpy.float) # E, dos up, dos down, idos up, idos down


  for i in range(enum_total):
      inline_total = inputFile.readline().split()
      if fermi_subtract == False: 
          dos[i][0] = float(inline_total[0])  
      else:
          dos[i][0] = float(inline_total[0]) - efermi_total 
      dos[i][1], dos[i][2], dos[i][3], dos[i][4] = float(inline_total[1]), float(inline_total[2]), float(inline_total[3]), float(inline_total[4])


  for row in dos:
      for element in row:
          outputFile_dos.write(str(element) + ' ')
      outputFile_dos.write('\n')

############


  if fermi_subtract == True: print 'Fermi level was subtracted'
  else: print 'Fermi level was NOT subtracted'

  enum_project = enum_total

  spacing = (emax_total - emin_total)/enum_total
  efermi_bin = int((efermi_total - emin_total)/spacing)

  pdos = numpy.zeros((enum_project, 1 + 3*3), dtype = numpy.float) # E, s up, s down, p up, p down, d up, d down, integrated

  for i in range(len(pdos)):

     pdos[i][0] = dos[i][0]  # fills in the E column of pdos


  # skip read the atoms you would like to look at
  for atoms in range(natom_skip):
      inline_project = inputFile.readline()

      for i in range(enum_project):
        inline_project = inputFile.readline().split()

  # read the atoms you care about
  for atoms in range(natom):
      inline_project = inputFile.readline()

      for i in range(enum_project):
        inline_project = inputFile.readline().split()
        pdos[i][1] += float(inline_project[1]) # s up
        pdos[i][2] += float(inline_project[2]) # s down
        pdos[i][3] += float(inline_project[3]) + float(inline_project[5]) + float(inline_project[7]) # py pz px up
        pdos[i][4] += float(inline_project[4]) + float(inline_project[6]) + float(inline_project[8]) # py pz px down
        pdos[i][5] += float(inline_project[9]) + float(inline_project[11]) + float(inline_project[13]) + float(inline_project[15]) + float(inline_project[17]) # dxy ... up
        pdos[i][6] += float(inline_project[10]) + float(inline_project[12]) + float(inline_project[14]) + float(inline_project[16]) + float(inline_project[18]) # dxy ... down


  for i in numpy.arange(1,enum_project):
      pdos[i][7] = pdos[i - 1][7] + pdos[i][1] + pdos[i][2]
      pdos[i][8] = pdos[i - 1][8] + pdos[i][3] + pdos[i][4]
      pdos[i][9] = pdos[i - 1][9] + pdos[i][5] + pdos[i][6]


  energies  = numpy.transpose(pdos)[0]
  s_up_project   = numpy.transpose(pdos)[1]
  p_up_project   = numpy.transpose(pdos)[3]
  d_up_project   = numpy.transpose(pdos)[5]
  s_down_project = numpy.transpose(pdos)[2]
  p_down_project = numpy.transpose(pdos)[4]
  d_down_project = numpy.transpose(pdos)[6]
#energies_core  =  energies[0:core_valence_divide]
#s_project_core = s_project[0:core_valence_divide]
#p_project_core = p_project[0:core_valence_divide]
#d_project_core = d_project[0:core_valence_divide]


  energies_valence  =  energies[core_valence_divide:efermi_bin]
  s_up_project_valence = s_up_project[core_valence_divide:efermi_bin]
  p_up_project_valence = p_up_project[core_valence_divide:efermi_bin]
  d_up_project_valence = d_up_project[core_valence_divide:efermi_bin]
  s_down_project_valence = s_down_project[core_valence_divide:efermi_bin]
  p_down_project_valence = p_down_project[core_valence_divide:efermi_bin]
  d_down_project_valence = d_down_project[core_valence_divide:efermi_bin]

#spd_sum_core    = sum(s_project_core)    + sum(p_project_core)    + sum(d_project_core)
  spd_sum_valence = sum(s_up_project_valence) + sum(p_up_project_valence) + sum(d_up_project_valence) + sum(s_down_project_valence) + sum(p_down_project_valence) + sum(d_down_project_valence)

#s_fraction_core = sum(s_project_core)/spd_sum_core
#p_fraction_core = sum(p_project_core)/spd_sum_core

  s_fraction_valence = (sum(s_up_project_valence)+sum(s_down_project_valence))/spd_sum_valence
  p_fraction_valence = (sum(p_up_project_valence)+sum(p_down_project_valence))/spd_sum_valence
  d_fraction_valence = (sum(d_up_project_valence)+sum(d_down_project_valence))/spd_sum_valence

#print
#print "%s core", s_fraction_core*100
#print "%p core", p_fraction_core*100
  print
  print "%s valence", s_fraction_valence*100
  print "%p valence", p_fraction_valence*100
  print "%d valence", d_fraction_valence*100
  print
  print "norm = ", str(s_fraction_valence + p_fraction_valence + d_fraction_valence)
  print
  print "efermi_total = ", efermi_total
  print "spacing = ", spacing
  print

  for row in pdos:
      for element in row:
          outputFile_pdos.write(str(element) + ' ')
      outputFile_pdos.write('\n')

 
if __name__ == '__main__':
     main()
