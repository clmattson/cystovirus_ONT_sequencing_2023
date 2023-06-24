#!/bin/bash


#inputs:

#-p $parent_dir directory to set as working
#-c $contig_suffix common, specific suffix of complete phylogentics sequence files (ex: _ALL_named_contigs_blast_len.fasta.fasta)


barcode_path=''
host_ref=''
assembly_name=''

print_usage() {
  printf "Usage: ..."
}

while getopts p:c: flag
do
    case "${flag}" in
        p) parent_dir=${OPTARG};;
        c) contig_suffix=${OPTARG};;

    esac
done


#set working dir
cd ${parent_dir} || exit
echo "working shell directory:"
pwd
echo


#run prokka!
for seqs_file in *"${contig_suffix}";
do prokka --kingdom Viruses --prefix ${seqs_file}_prokka ${seqs_file};
done
