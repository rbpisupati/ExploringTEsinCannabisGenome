# Make a soft link of the main and the control file 
# Copy the script to the location contig_positions
# Run the script directly without any inputs

# Generating control file
# running main file

# Get working directory
workDir=`(pwd)`

for workFile in CL*Contig*fasta
do

# Sed the variable
  workFile=$(sed 's|\.fasta||' <<< $workFile)

# Make a new directory for each contig and create symbolic link
  mkdir $workFile
  mv "$workFile".fasta $workFile/
  cd $workFile
  ln -sf ~/AnteR/main "$workDir"/"$workFile"/main

# Edit control file each time you change the input file
  echo "SEQUENCE_FILE "$workFile".fasta
NUM_ROUNDS 3" > control.txt

# Run the main script of AnteR for all the contigs
  ./main > log.main.txt 2> log.main.err

  cd ..
  
done
