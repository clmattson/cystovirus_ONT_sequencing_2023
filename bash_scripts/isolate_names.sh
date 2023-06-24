#!/bin/bash


#inputs:

#-p $parent_dir directory to set as workign
#-c $contig_suffix common, specific suffix of contig files (.fasta)
#-s $s_key tab sep barcode - isolate key files with S suffix, no path needed
#-m $m_key " "M
#-l $l_key " "L



barcode_path=''
host_ref=''
assembly_name=''

print_usage() {
  printf "Usage: ..."
}

while getopts p:c:s:m:l: flag
do
    case "${flag}" in
        p) parent_dir=${OPTARG};;
        c) contig_suffix=${OPTARG};;
 s) s_key=${OPTARG};;
 m) m_key=${OPTARG};;
 l) l_key=${OPTARG};;

    esac
done

#set working dir
cd ${parent_dir} || exit
echo "working shell directory:"
pwd
echo


#small seg assembly names:

awk '
FNR==NR {f2[$1]=$2;next}
/^>/ {
  for (i in f2) {
    if (index(substr($1,2), i)) {
      print ">" f2[i]; next
    }
  }
}1' ${s_key} [sS]"${contig_suffix}" >> S_named"${contig_suffix}"



#med seg assembly names:

awk '
FNR==NR {f2[$1]=$2;next}
/^>/ {
  for (i in f2) {
    if (index(substr($1,2), i)) {
      print ">" f2[i]; next
    }
  }
}1' ${m_key} [mM]"${contig_suffix}" >> M_named"${contig_suffix}"



#large seg assembly names:

awk '
FNR==NR {f2[$1]=$2;next}
/^>/ {
  for (i in f2) {
    if (index(substr($1,2), i)) {
      print ">" f2[i]; next
    }
  }
}1' ${l_key} [lL]"${contig_suffix}" >> L_named"${contig_suffix}"


