#!/usr/bin/python
import sys
# Run with AP.py ALA-ALA etc.
PEPTIDE=sys.argv[1]

# Open file to read.
initial_sasa_file=open("{}_sasa_init.xvg".format(PEPTIDE), "r")
# Go through all the lines in the file. We're only interested in the last.
for initial_line in initial_sasa_file:
	pass
# The last line is left after the end for the loop. Split the line at spaces.
initial_line_split = initial_line.split(' ')
# The SASA value is the last item in the line, and is a string.
sasa_value_string = initial_line_split[-1]
# Converting it to a floating point (computer represendation of decimal) number
initial_sasa = float(sasa_value_string)

# The same for the final sasa value.
final_sasa_file=open("{}_sasa_final.xvg".format(PEPTIDE), "r")
for final_line in final_sasa_file:
	pass
final_sasa = float(final_line.split(' ')[-1])

# Calculate AP
AP = initial_sasa / final_sasa
# Output the number as eg 1.002 or 2.050
print('{:.3f}'.format(AP))
