#!/bin/bash


#inputs:

#-p $parent_dir directory to set as workign
#-c $contig_suffix common, specific suffix of contig files (M"_XYZ.fasta")
#-x $Phi6_S fasta with phi 6 S segment ref only
#-y $Phi6_M fasta with phi 6 M segment ref only
#-z $Phi6_L fasta with phi 6 L segment ref only
#-s $S_db database of S segment refs
#-m $M_db database of M segment refs
#-l #L_db database of L segment refs




barcode_path=''
host_ref=''
assembly_name=''

print_usage() {
  printf "Usage: ..."
}

while getopts p:c:x:y:z:s:m:l: flag
do
    case "${flag}" in
        p) parent_dir=${OPTARG};;
        c) contig_suffix=${OPTARG};;
 x) Phi6_S=${OPTARG};;
 y) Phi6_M=${OPTARG};;
 z) Phi6_L=${OPTARG};;
 s) S_db=${OPTARG};;
 m) M_db=${OPTARG};;
 l) L_db=${OPTARG};;

    esac
done

#set working dir
cd ${parent_dir} || exit
echo "working shell directory:"
pwd
echo

#ex: ${contig_suffix}=="_named_contigs_BL.fasta"

#add *phi6 references* to beginning of each file of sequences

cat ${Phi6_S} [sS]"${contig_suffix}" >> S_P6"${contig_suffix}"
cat ${Phi6_M} [mM]"${contig_suffix}" >> M_P6"${contig_suffix}"
cat ${Phi6_L} [lL]"${contig_suffix}" >> L_P6"${contig_suffix}"


#mafft-linsi - dont do for now?
#for contigs_file in *_P6"${contig_suffix}";
#do mafft-linsi --adjustdirectionaccurately --maxiterate 1000 --localpair ${contigs_file} > "${contigs_file}".mafft.fasta;
#done


#mafft einsi, DO THIS FOR NOW

mafft-einsi --genafpair --adjustdirectionaccurately --maxiterate 1000 S_P6"${contig_suffix}" > S_P6"${contig_suffix}".mafft.fasta
mafft-einsi --genafpair --adjustdirectionaccurately --maxiterate 1000 M_P6"${contig_suffix}" > M_P6"${contig_suffix}".mafft.fasta
mafft-einsi --genafpair --adjustdirectionaccurately --maxiterate 1000 L_P6"${contig_suffix}" > L_P6"${contig_suffix}".mafft.fasta



#add db refs:

mafft --add ${S_db} --maxiterate 1000 --oldgenafpair S_P6"${contig_suffix}".mafft.fasta > S_ALL"${contig_suffix}".mafft.fasta
mafft --add ${M_db} --maxiterate 1000 --oldgenafpair M_P6"${contig_suffix}".mafft.fasta > M_ALL"${contig_suffix}".mafft.fasta
mafft --add ${L_db} --maxiterate 1000 --oldgenafpair L_P6"${contig_suffix}".mafft.fasta > L_ALL"${contig_suffix}".mafft.fasta

#or  ?
#mafft --seed --oldgenafpair mafft/l_seqs_blast_len.fasta.mafft.fasta /mnt/data0/cysto_refs/cystovirus_L_db.fasta > mafft/l_seqs_seed_test.mafft.fasta

#make regular FASTA of all aligned seqs:
cat S_P6"${contig_suffix}" ${S_db} >> S_ALL"${contig_suffix}"
cat M_P6"${contig_suffix}" ${M_db} >> M_ALL"${contig_suffix}"
cat L_P6"${contig_suffix}" ${L_db} >> L_ALL"${contig_suffix}"

#move new files to new dir
#mkdir ./mafft_"${contig_suffix}"
#mv *ALL* ./mafft_"${contig_suffix}"
#mv *P6* ./mafft_"${contig_suffix}"
