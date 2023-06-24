#!/bin/bash


#inputs:

#-p $barcode_path parent dir of barcode folders
#-h $host_ref bacterial host reference file path
#-a $assembly_name desired prefix for assemblies; ex "CanuAssembly_MRL300_MOL200"


barcode_path=''
host_ref=''
assembly_name=''

print_usage() {
  printf "Usage: ..."
}

while getopts p:h:a: flag
do
    case "${flag}" in
        p) barcode_path=${OPTARG};;
        h) host_ref=${OPTARG};;
 a) assembly_name=${OPTARG};;


    esac
done




#ok lets align & remove host reads
#khmer trim? skip for now


#1 map to host with minimmap2
#minimap is in conda env seqqc

#host: p.savastanoi.ref.nbrks.fasta
for barcode in ${barcode_path}/barcode*;
do bn=$(basename $barcode);
 echo "barcode path: ${barcode_path}, (loop index) barcode:  ${barcode}; bn: ${bn}"
 minimap2 -x map-ont -k15 -a -t 8 ${host_ref} ${barcode}/*all*fasta > ${barcode}/${bn}_host_minimap;
 samtools view --threads 24 -b -o ${barcode}/${bn}_host.bam ${barcode}/${bn}_host_minimap;
 samtools sort -o ${barcode}/${bn}_host_sorted.bam -T ${barcode}/${bn}_samtools_temp_dir --threads 24 ${barcode}/${bn}_host.bam;
 #view mapping info for reads
 samtools flagstat ${barcode}/${bn}_host_sorted.bam >> ${barcode}/${bn}_host_BAMflagstat.txt;
 #filter mapped (host) reads with samtools
 samtools fasta -n -f 4 ${barcode}/${bn}_host_sorted.bam > ${barcode}/${bn}.unmapped.fasta;
done



#2 assemble!

#reexport canu to path (no conda env needed)
export PATH=$PATH:~/canu-2.2/bin


#non host read for each barcode now in *unmapped.fasta
#with shorter allowed reads and overlaps
#this one!!
for barcode in ${barcode_path}/barcode*; do bn=$(basename $barcode); echo; echo ${bn}; echo; canu -p ${bn}_canu_assembly -d ${barcode}/${bn}_${assembly_name} genomeSize=5k -minReadLength=300 -minOverlapLength=200 -nanopore ${barcode}/${bn}.unmapped.fasta; done



#3. polishing with medaka:
#conda env: medaka

for barcode in ${barcode_path}/barcode*; do bn=$(basename $barcode);
medaka_consensus -i ${barcode}/${bn}.all.fastq -d ${barcode}/${bn}_${assembly_name}/${bn}_canu_assembly.contigs.fasta -o ${barcode}/${bn}_${assembly_name}/medaka -m r941_min_hac_g507 -t 8; done


#name medaka consensus files with barcode name

for barcode in ${barcode_path}/barcode*; do bn=$(basename $barcode); mv ${barcode}/${bn}_${assembly_name}/medaka/consensus.fasta ${barcode}/${bn}_${assembly_name}/medaka/${bn}_medaka_consensus.fasta; done
