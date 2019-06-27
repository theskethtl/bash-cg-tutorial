#!/bin/bash
# Written by Travis Hesketh, 2019. (C) TuttleLab, University of Strathclyde.
# This script is the basic version of the AP calculation script.
# The comments in this file should help you to understand what every line does.

# Loading the GROMACS module. We need this for the SASA command.
module load gromacs/intel-2018.2/2016.5-single-wee-archie
AP_FILE="dipep_ap.txt"

# The amino acids we're interested in. Usually the full set. We use HSE instead of histidine
# because this is how it is represented in the CHARMM forcefield.
ACID_CODES=(ALA CYS GLY HSE ILE LEU MET ASN PRO GLN SER THR VAL ASP GLU LYS ARG PHE TRP TYR)
ONE_LETTER_CODES=(A C G H I L M N P Q S T V D E K R F W Y)

# These 'for' loops go through all the items in the array above. Because we're "nesting" them,
# the second loop goes through the whole loop each time the first goes round. In other words,
# we get all possible combinations:
# ALA-CYS, ALA-GLY, ALA-HSE, ... TYR-TYR
NUM_ACIDS=$((${#ACID_CODES[@]}-1))
for INDEX_ONE in $(seq 0 ${NUM_ACIDS}); do
	for INDEX_TWO in $(seq 0 ${NUM_ACIDS}); do
		# The name of our peptide (e.g. ALA-ALA)
		ACID_ONE=${ACID_CODES[${INDEX_ONE}]}
		ACID_TWO=${ACID_CODES[${INDEX_TWO}]}

		# These aren't used because WF, CL and NA are reserved names, but helpful for tri/tetrapeptides.
		ACID_ONE_OLC=${ONE_LETTER_CODES[${INDEX_ONE}]}
		ACID_TWO_OLC=${ONE_LETTER_CODES[${INDEX_TWO}]}

		PEPTIDE=${ACID_ONE}-${ACID_TWO} # For tri/tetrapeptides, use PEPTIDE=${ACID_ONE_OLC}${ACID_TWO_OLC}${ACID_THREE_OLC} etc.

		# If the peptide directory has not been made, there's no point in running this script.
		if [ ! -d "${PEPTIDE}" ]; then
			echo "Peptide ${PEPTIDE} hasn't been created yet."
			echo "You probably need to run the submission script first."
			exit
		fi

		cd ${PEPTIDE}

		# If the peptide is unfinished, we don't need to calulate the AP score.
		# We can check to see if it's unfinished by testing to see if the file does not exist ([ ! -f FILE ])
		# OR (||) the log does not contain "Finished mdrun".
		if [ ! -f "${PEPTIDE}_eq.log" ] || ! grep "Finished mdrun" ${PEPTIDE}_eq.log > /dev/null; then
			echo "Peptide ${PEPTIDE} not finished. Not calculating AP score."
			cd ../
			continue
		fi
		# If there's an AP file [ -f FILE ] AND (&&) the peptide has already had the AP calculated we don't need to do it again
		if [ -f ../"${AP_FILE}" ] && grep ${PEPTIDE} ../${AP_FILE} > /dev/null; then
			echo "Peptide ${PEPTIDE} AP already calculated."
			cd ../
			continue
		fi

		echo "Calculating AP for ${PEPTIDE}"
		# We only need to do each of these if the file hasn't been made yet, so we have tests for each one.

		# Calculate the initial SASA score from the minimised GRO file.
		if [ ! -f "${PEPTIDE}_sasa_init.xvg" ]; then
			echo Protein | gmx_mpi sasa -f ${PEPTIDE}_min.gro -s ${PEPTIDE}_min.tpr -o ${PEPTIDE}_sasa_init.xvg
		fi

		if [ ! -f "${PEPTIDE}_eq_centred.gro" ] || [ ! -f "${PEPTIDE}_sasa_final.xvg"]; then
			# Cluster and center the final gro file of the equilibrated protein and write out the coordinates of the whole system.
			echo Protein Protein System | gmx_mpi trjconv -f ${PEPTIDE}_eq.gro -s ${PEPTIDE}_eq.tpr -pbc cluster -center -o ${PEPTIDE}_eq_centred.gro
			# And calculate the final SASA score from the clustered coordinates.
			echo Protein | gmx_mpi sasa -f ${PEPTIDE}_eq_centred.gro -s ${PEPTIDE}_eq.tpr -o ${PEPTIDE}_sasa_final.xvg
		fi
		# Use AP.py to calculate AP score.
		AP=$(python ../AP.py ${PEPTIDE})
		# And add it to the AP file.
		echo "${PEPTIDE} ${AP}" >> ../${AP_FILE}
		rm ./\#* *.edr *.mdp > /dev/null 2>&1
		cd ../
	done
done

echo ""
echo "All molecules done. Sorting AP file."

# Sorting by AP.
sort -rk2 -o ${AP_FILE} ${AP_FILE}
