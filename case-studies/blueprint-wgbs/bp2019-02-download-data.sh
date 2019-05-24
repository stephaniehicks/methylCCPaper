cat /users/shicks1/data/DNAm/blueprint_ihec/blueprint_blood_ftp_paths.csv | while read line
do
wget $line
done
