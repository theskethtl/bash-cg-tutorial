package require psfgen
topology ../top_all36_prot.rtf
segment PEPTIDE {pdb PEPTIDE_aa.pdb}
coordpdb PEPTIDE_aa.pdb PEPTIDE
guesscoord
writepdb PEPTIDE_aa.pdb
exit
