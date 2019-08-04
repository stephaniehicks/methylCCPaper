#$ -cwd
#$ -l mem_free=5G,h_vmem=5G

mkdir files_bigwig
mkdir files_bsseq

cd files_bigwig
cat /users/shicks1/data/DNAm/blueprint_ihec/blueprint_blood_ftp_paths.csv | while read line
do
wget $line
done
