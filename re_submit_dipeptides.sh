#!/bin/bash
# Written by Travis Hesketh, 2019. (C) TuttleLab, University of Strathclyde.
# This script is the basic version of the dipeptide re-submission script.
# The comments in this file should help you to understand what every line does.

# Loading the VMD and GROMACS modules. We need these for the simulation setup.
module load vmd/1.9.1
module load gromacs/intel-2018.2/2016.5-single

# The amino acids we're interested in. Usually the full set. We use HSE instead of histidine
# because this is how it is represented in the CHARMM forcefield.
ACID_CODES=(ALA CYS GLY HSE ILE LEU MET ASN PRO GLN SER THR VAL ASP GLU LYS ARG PHE TRP TYR)

# These 'for' loops go through all the items in the array above. Because we're "nesting" them,
# the second loop goes through the whole loop each time the first goes round. In other words,
# we get all possible combinations:
# ALA-CYS, ALA-GLY, ALA-HSE, ... TYR-TYR

# We're assuming the coarse graining has gone okay (otherwise something has gone seriously wrong).
# Assign the contents of our queue to a variable.
QUEUE_CONTENTS=$(squeue -u ${USER} -o '%.100j')

for ACID_ONE in ${ACID_CODES[@]}; do
	for ACID_TWO in ${ACID_CODES[@]}; do
		PEPTIDE=${ACID_ONE}-${ACID_TWO}

		# We need to use some test logic to determine the best course of action.
		# If the peptide directory has not been made, there's no point in running this script.
		if [ ! -d "${PEPTIDE}" ]; then
			echo "Peptide ${PEPTIDE} hasn't been created yet."
			echo "You probably need to run the submission script first."
			exit
		fi

		# If the peptide is still in the queue, it hasn't failed. We don't need to re-submit it.
		# We can use squeue to get our jobs (-o '%.100j' means give us 100 characters of job name)
		# grep checks for our peptide.
		if echo ${QUEUE_CONTENTS} | grep "${PEPTIDE}" > /dev/null; then
			echo "${PEPTIDE} is still running."
			# This continue statement means 'go to the next peptide in the loop'.
			continue
		fi

		# Okay, so the folder exists and the peptide isn't running. Change into the peptide folder.
		cd ${PEPTIDE}

		# If there's an equilibration log, we can check if it's finished.
		if [ -f "${PEPTIDE}_eq.log" ]; then
			# If the job is finished, we don't need to do anything. Go back to root folder and on to next peptide
			if grep "Finished mdrun" ${PEPTIDE}_eq.log > /dev/null; then
				echo "Peptide ${PEPTIDE} completed equilibration successfully."
				cd ../
				continue
			# Otherwise, something must be wrong with the equilibration.
			# Let's see if the time limit was the issue first.
			elif grep "TIME LIMIT" ${PEPTIDE}_slurm.log > /dev/null; then
				echo "Peptide ${PEPTIDE} failed due to time limit."
				# If we haven't changed the time limit already, we need to exit here. Change time limit.
				if [ ! "$1" == "-i" ]; then
					# Run script with "./re_submit_failed.sh -i" to skip this exit command.
					echo "Change the time limit in configuration_dipeptides.sh and restart this script with the -i option."
					cd ../
					exit
				fi
			fi
		fi
		echo "Peptide ${PEPTIDE} failed. Trying to re-submit..."

		# We'll try to go through the whole procedure again. New starting positions can make it more likely a job will work.
		# These commands are the same as submit_dipeptides.sh
		gmx_mpi insert-molecules -f ../water_125A.gro -ci ${PEPTIDE}.pdb -radius 0.25 -replace W -nmol 300 -o ${PEPTIDE}_box.gro
		N_WATER=$(python2 ../count_water.py ./${PEPTIDE}_box.gro)
		sed "s/PEPTIDE/${PEPTIDE}/g; s/N_MOLECULES/300/g; s/N_WATER/${N_WATER}/g" ../template.top > ${PEPTIDE}.top
 		gmx_mpi grompp -f ../min_parameter_file.mdp -c ${PEPTIDE}_box.gro -p ${PEPTIDE}.top -o ${PEPTIDE}_min.tpr
		echo Water | gmx_mpi genion -s ${PEPTIDE}_min.tpr -neutral -pname NA -nname CL -p ${PEPTIDE}.top -o ${PEPTIDE}_box.gro
		mpirun -np 1 gmx_mpi mdrun -ntomp 1 -deffnm ${PEPTIDE}_min
		gmx_mpi grompp -f ../eq_parameter_file.mdp -c ${PEPTIDE}_min.gro -p ${PEPTIDE}.top -o ${PEPTIDE}_eq.tpr -maxwarn 2
		sed "s/PEPTIDE/${PEPTIDE}/g" ../slurm_template.sh > ${PEPTIDE}.sh
		sbatch ${PEPTIDE}.sh
		rm ./\#* *.edr *.mdp *_prev.cpt > /dev/null 2>&1
		cd ../
	done
done

echo ""
echo "All molecules done."
