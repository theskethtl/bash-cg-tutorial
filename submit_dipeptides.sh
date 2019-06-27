#!/bin/bash
# Written by Travis Hesketh, 2019. (C) TuttleLab, University of Strathclyde.
# This script is the basic version of the dipeptide submission script.
# The comments in this file should help you to understand what every line does.

# To adapt this file for tripeptides/tetrapeptides etc.:
#	Add additional for loops for each extra acid [after line 43] and acid definitions [line 47]
#		- e.g. for INDEX_THREE in $(seq 0 ${NUM_ACIDS}); do [line 43] and ACID_THREE=${ACID_CODES[${INDEX_THREE}]} [line 47]
#		- Remember to add additional 'done' statements [before line 133] for each new loop.
#		- It is helpful for troubleshooting to indent lines to match the new logic levels.
#	Change the definition of PEPTIDE to include your new acids. [line 51]
#		- For longer peptides we can use one letter codes instead of three letter codes
#		  (WF, NA and CL all clash with reserved MARTINI names)
#	Alter WX.pdb to have backbone atoms for your new acids.
#		- There are templates for tri and tetrapeptides in the dropbox (WXY.pdb, WXYZ.pdb)
#	Adjust the sed statement [line 67] to account for your new acids and template.
#		- e.g. sed "s/WWW/${ACID_ONE}/g; s/XXX/${ACID_TWO}/g; s/YYY/${ACID_THREE}/g" ../WXY.pdb > ${PEPTIDE}_aa.pdb
#	Update the secondard strucure in the martinize.py call [line 68] to contain extra E characters for each acid (e.g. -ss EEE)
#	Remember to update the resubmission script and AP calculation scripts.

# Loading the VMD and GROMACS modules. We need these for the simulation setup.
module load vmd/1.9.1
module load gromacs/intel-2018.2/2016.5-single

# The amino acids we're interested in. Usually the full set. We use HSE instead of HIS
# because this is how histidine is represented in the CHARMM forcefield.
# If you're interested in generating a subset of dipeptides, you could define multiple arrays here.
#ACID_CODES=(ALA CYS GLY HSE ILE LEU MET ASN PRO GLN SER THR VAL ASP GLU LYS ARG PHE TRP TYR)
#ONE_LETTER_CODES=(A C G H I L M N P Q S T V D E K R F W Y)
ACID_CODES=(ALA)
ONE_LETTER_CODES=(A)

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

		PEPTIDE=${ACID_ONE}-${ACID_TWO}  # For tri/tetrapeptides, use PEPTIDE=${ACID_ONE_OLC}${ACID_TWO_OLC}${ACID_THREE_OLC} etc.

		# Making a directory called the name of our peptide, and moving into it. This is to keep all
		# the files relevant to each system separated.
		mkdir ${PEPTIDE}
		cd ${PEPTIDE}
		# Note: this now means all the files in the root folder are up one level (need to be prefaced with ../)

		# First, we build an atomistic structure for our peptide
		# We can use VMD to fill in the missing coordinates from our backbone pdb file (WX.pdb)
		# VMD's psfgen package knows how to fill them because the CHARMM topology file (top_all36_prot.rtf),
		# contains all the necessary angles and bond lengths to calculate the structure.

		# First, we need to use the sed command to swap "WWW" and "XXX" for the three letter codes of
		# our acids, and write this to the name of our peptide (${PEPTIDE}_aa.pdb, where aa=all atoms). This is our backbone template.
		sed "s/WWW/${ACID_ONE}/g; s/XXX/${ACID_TWO}/g" ../WX.pdb > ${PEPTIDE}_aa.pdb
		# Then, we need to do the same for the set of VMD commands, so that they refer to this particular peptide.
		sed "s/PEPTIDE/${PEPTIDE}/g" ../structure_gen.tcl > ${PEPTIDE}.pgn
		# And finally we can use VMD to fill in the missing coordinates.
		vmd -dispdev text -e ${PEPTIDE}.pgn
		# We now have a full atomistic structure for our peptide (${PEPTIDE}_aa.pdb)

		# We need to convert this to a coarse grained structure using martinize.py.
		# MARTINI expects histidine to be called 'HIS' rather than 'HSE', so first we need to swap HSE for HIS.
		# The -i flag means "in place" (i.e. don't write the output to a new file or to the terminal).
		sed -i "s/HSE/HIS/g" ${PEPTIDE}_aa.pdb
		python2 ../martinize.py -ff martini22 -f ${PEPTIDE}_aa.pdb -x ${PEPTIDE}.pdb -o ${PEPTIDE}.top -name ${PEPTIDE} -ss EE
		# And the coarse graining is done.
		# This produces an itp file (${PEPTIDE}.itp), a top file (${PEPTIDE}.top), and a coarse grained pdb file (${PEPTIDE}.pdb).

		# In the past, we built up our system by adding N_MOLECULES protein molecules to an empty box and adding water
		# This was error prone (it's harder to solvate to the right density), so now we take a box of water of known density
		# and add protein, replacing water in the process. This is a little bit more complicated where the top file is concerned,
		# but is more consistent.
		# If you are unsure about what each option for a GROMACS commadn does (-f -ci etc.), load the GROMACS module
		# and use the '-h' flag after the name of the command.
		# For insert-molecules, that would be like this:
		# 			gmx_mpi insert-molecules -h
		# This will produce a help message with information about what the command and all its options do.
		gmx_mpi insert-molecules -f ../water_125A.gro -ci ${PEPTIDE}.pdb -radius 0.25 -replace W -nmol 300 -o ${PEPTIDE}_box.gro

		# Great. We now we have a water box with $N_MOLECULES molecules in it (${PEPTIDE}_box.gro)! Unfortunately, our '.top' file
		# no longer matches the system. We could add the number of water molecules to the existing file, but because the files need to be
		# in order (i.e. the number of W molecules needs to be above the peptide) it's actually easier to just remake them from a template.

		# We have a Python script (count_water.py) which counts the number of molecules in a MARTINI gro file.
		# Let's first use that to get the number of water molecules
		N_WATER=$(python2 ../count_water.py ./${PEPTIDE}_box.gro)
		# And then use sed to fill in our topology template, swapping in the peptide we're using and the correct number
		# of each type of molecule.
		sed "s/PEPTIDE/${PEPTIDE}/g; s/N_MOLECULES/300/g; s/N_WATER/${N_WATER}/g" ../template.top > ${PEPTIDE}.top

		# Now we can set up the minimisation using our minimisation parameter file (min_parameter_file.mdp)
 		gmx_mpi grompp -f ../min_parameter_file.mdp -c ${PEPTIDE}_box.gro -p ${PEPTIDE}.top -o ${PEPTIDE}_min.tpr

		# But before we run it, we might need to replace water with some ions so that our system is neutral.
		# Unlike insert-molecules, this automatically updates the top file.
		# The "echo Water |" at the start tells genion which group to replace with ions (it would stop and ask us otherwise).
		echo Water | gmx_mpi genion -s ${PEPTIDE}_min.tpr -neutral -pname NA -nname CL -p ${PEPTIDE}.top -o ${PEPTIDE}_box.gro

		# The following command runs the minimisation. This can be done locally as it's quick and not computationally demanding.
		# We need to preface this with mpirun as the version of GROMACS we're using is compiled to work with
		# the 'message passing interface' so that multiple processes can communicate. In this case, we just want one core.
		mpirun -np 1 gmx_mpi mdrun -ntomp 1 -deffnm ${PEPTIDE}_min

		# After this, we can set up the equilibration. This is similar to the minimisation setup,
		# but uses the minimised coordinates and the equilibration parameters instead.
		gmx_mpi grompp -f ../eq_parameter_file.mdp -c ${PEPTIDE}_min.gro -p ${PEPTIDE}.top -o ${PEPTIDE}_eq.tpr -maxwarn 2

		# Unfortunately, we can't run this locally. We need to swap out the slurm submission template and submit it to the cluster.
		sed "s/PEPTIDE/${PEPTIDE}/g" ../slurm_template.sh > ${PEPTIDE}.sh
		# The following command submits the job to the cluster.
		# We're using command line flags so that we don't need to sed 4 different things into the submission template.
		# but these comand line options can also be given with e.g. "#SBATCH --partition=standard" in the sh file.
		sbatch ${PEPTIDE}.sh
		# This rm command removes awkward temporary files (GROMACS generates a lot of them), which help keep the folder neat.
		# The redirection of the output (>) and error (2>) streams to /dev/null means we won't see an output complaining
		# if the files are missing.
		rm ./\#* *.edr *.mdp > /dev/null 2>&1
		# And changing back up to the root folder (because we moved into the peptide directory after we created it).
		cd ../
	done
done

echo ""
echo "All molecules done."
