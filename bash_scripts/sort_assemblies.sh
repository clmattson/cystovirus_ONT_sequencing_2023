#!/bin/bash

#get input:

#-p ${parent_dir} directory containing barcode folders
#-n ${assembly_suffix} assembly file name suffix (follows 'barcodeXX', ex: _medaka_consensus.fasta)
#-s ${host_s} S segment reference file & path
#-m ${host_m} M segment reference file & path
#-l ${host_l} L segment reference file & path


parent_dir=''
assembly_suffix=''
host_s=''
host_m=''
host_l=''


print_usage() {
  printf "Usage: ..."
}

while getopts p:n:s:m:l: flag
do
    case "${flag}" in
        p) parent_dir=${OPTARG};;
 n) assembly_suffix=${OPTARG};;
 s) host_s=${OPTARG};;
 m) host_m=${OPTARG};;
 l) host_l=${OPTARG};;

    esac
done


cd ${parent_dir} || exit

echo "working shell directory:"
pwd
echo

#get assemblies with 3 contigs only, put in an array:
#for file in *"${assembly_suffix}";

#       do num_contigs="$(grep -c ">" "$file")";
#        if [ "${num_contigs}" -eq 3 ]; then

#               contigs_array+=("$file");
#       fi;
#done


#blast all assemblies against phi6 M, S, and L separately:

#S segment
#for assembly in "${contigs_array[@]}";
for assembly in *"${assembly_suffix}";

 #get file name only:
 do filename=$(basename ${assembly});
        #gets only barcode name from filename (hard coded, gets first 9 characters only)
        barcode="${filename:0:9}";
        #add barcode name to contig headers if not already there
 header_check="$(grep -c "barcode" "$assembly")"
 if [ "${header_check}" -eq 0 ]; then
        sed -i "s/>/>${barcode}\//I" ${assembly};
 fi;
        #blast the current assembly against phi6 S segment reference
        blastn -query ${assembly} -subject ${host_s} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${barcode}_S_blast_out.csv;
        blast_check="$(cat ${barcode}_S_blast_out.csv)"
        if [ -n "${blast_check}" ]; then
                echo "${barcode} ${blast_check} not empty"
                #makes file with contig names only
                awk -F"," '{print $1}' ${barcode}_S_blast_out.csv > ${barcode}_S_blast_readlist.txt;
                #gets most frequent read
                s_top="$(cut -d' ' -f1 ${barcode}_S_blast_readlist.txt | sort | uniq -c | sort -rn | head -n1 | awk '{print $2}')";
                echo "${barcode}_${s_top}_"
                echo ${s_top} >> s_contigs_blast.txt;
                #if [ -n "{$s_top}" ]; then
                #echo "${barcode}_${s_top}_not empty!";
                cat s_contigs_blast.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${assembly} >> s_contigs_blast.fasta;
        fi;
done


#M segment
#for assembly in "${contigs_array[@]}";
for assembly in *"${assembly_suffix}";

        #get file name only:
 do filename=$(basename ${assembly});
        #gets only barcode name from filename (hard coded, gets first 9 characters only)
        barcode="${filename:0:9}";
        #can add line here to remove contig header whitespace if necessary
        #blast the current assembly against phi6 M segment reference
        blastn -query ${assembly} -subject ${host_m} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${barcode}_M_blast_out.csv;
        blast_check="$(cat ${barcode}_M_blast_out.csv)"
        if [ -n "${blast_check}" ]; then
                echo "${barcode} ${blast_check} not empty"
                #makes file with contig names only
                awk -F"," '{print $1}' ${barcode}_M_blast_out.csv > ${barcode}_M_blast_readlist.txt;
                #gets most frequent read
                m_top="$(cut -d' ' -f1 ${barcode}_M_blast_readlist.txt | sort | uniq -c | sort -rn | head -n1 | awk '{print $2}')";
                echo "${barcode} _${m_top}_"
                echo ${m_top} >> m_contigs_blast.txt;
                #if [ -n "{$m_top}" ]; then
                #echo "${barcode}_${m_top}_not empty!";
                cat m_contigs_blast.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${assembly} >> m_contigs_blast.fasta;
        fi;
done


#L segment
#or assembly in "${contigs_array[@]}";
for assembly in *"${assembly_suffix}";

        #get file name only:
 do filename=$(basename ${assembly});
        #gets only barcode name from filename (hard coded, gets first 9 characters only)
        barcode="${filename:0:9}";
        #can add line here to remove contig header whitespace if necessary
        #blast the current assembly against phi6 L segment reference
        blastn -query ${assembly} -subject ${host_l} -outfmt '6 delim=,' -qcov_hsp_perc 20 >> ${barcode}_L_blast_out.csv;
        blast_check="$(cat ${barcode}_L_blast_out.csv)"
        if [ -n "${blast_check}" ]; then
                echo "${barcode} ${blast_check} not empty"
                #makes file with contig names only
                awk -F"," '{print $1}' ${barcode}_L_blast_out.csv > ${barcode}_L_blast_readlist.txt;
                #gets most frequent read
                l_top="$(cut -d' ' -f1 ${barcode}_L_blast_readlist.txt | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}')";
                echo "${barcode}_${l_top}_"
                echo ${l_top} >> l_contigs_blast.txt;
                #if [ -n "{$l_top}" ]; then
                #echo "${barcode}_${l_top}_not empty!";
                cat l_contigs_blast.txt | awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - ${assembly} >> l_contigs_blast.fasta;
        fi;
done


#remove sequences shorter than:
# 2000 S
# 3500 M
# 5000 L

seqkit seq -m 2000 --remove-gaps s_contigs_blast.fasta > s_contigs_blast_len.fasta

seqkit seq -m 3500 --remove-gaps m_contigs_blast.fasta > m_contigs_blast_len.fasta

seqkit seq -m 5000 --remove-gaps l_contigs_blast.fasta > l_contigs_blast_len.fasta
