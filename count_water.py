#!/usr/bin/python
# Count the MARTINI waters in a file.
# Should be Py2/3 compatible.
import sys
water_residues, water_count = set(('W', 'WF')), 0
with open(sys.argv[1]) as gro_file:
	for line in gro_file:
		try:
			resname = line[5:10].rstrip(' ')
		except ValueError:
			pass
		else:
			if resname in water_residues:
				water_count += 1
print(water_count)
