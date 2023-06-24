#!/bin/bash

#bash script for initial processing seq data

#input to gather:
#p $fast5_pass_path path where the fast5_pass is
#f $fast5_fail_path path where the fast5_fail is

#b $pass_bc_path path for passed fast5s once basecalled
#l $fail_bc_path path for failed fast5s once basecalled

#d $demux_pass_all path for all passed, demultiplexed, trimmed fastq's

#c $cysto_blast_db for the desired cystovirus blast database, for now: mnt/data0/wild.cysto.ref.genomes.fasta
#h $pp_blast_db for the desired pp blast database, for now: /mnt/data0/pp_ref.fasta


fast5_pass_path=''
fast5_fail_path=''
pass_bc_path=''
fail_bc_path=''
demux_pass_all=''
cysto_blast_db=''
host_blast_db=''

print_usage() {
        printf "Usage: #input to gather"
}

while getopts p:f:b:l:d:c:h: flag
do
    case "${flag}" in
 p) fast5_pass_path=${OPTARG};;
 f) fast5_fail_path=${OPTARG};;
 b) pass_bc_path=${OPTARG};;
 l) fail_bc_path=${OPTARG};;
 d) demux_pass_all=${OPTARG};;
 c) cysto_blast_db=${OPTARG};;
 h) host_blast_db=${OPTARG};;
    esac
done


#Basecalling and demultiplexing have to be very careful/specific as we want to re-basecall all (passed and failed) fast5's but filenames are the same between the fast5_pass and fast5_fail dirs. we must not accidentally write over our data when processing the second group

#make log files
echo > ${demux_pass_all}/output/basecalling.log
echo > ${demux_pass_all}/output/barcoding.log

#basecalling

#basecall fast5_pass:
~/ont-guppy/bin/guppy_basecaller --input_path $fast5_pass_path --save_path $pass_bc_path --config dna_r9.4.1_450bps_hac.cfg -x cuda:0 --num_callers 8 --verbose_logs

#basecall fast5_fail:
~/ont-guppy/bin/guppy_basecaller --input_path $fast5_fail_path --save_path $fail_bc_path --config dna_r9.4.1_450bps_hac.cfg -x cuda:0 --num_callers 8 --verbose_logs


#demultiplexing:

#demux passed fastq's from fast5_pass

#require both barcodes! - OFF
#~/ont-guppy/bin/guppy_barcoder --require_barcodes_both_ends -i ${pass_bc_path}/pass -s ${demux_pass_all} --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg" --trim_barcodes -x cuda:0 --worker_threads 8 --verbose_logs

#Require one barcode!
~/ont-guppy/bin/guppy_barcoder -i ${pass_bc_path}/pass -s ${demux_pass_all} --config configuration.cfg --barcode_kits "EXP-NBD104 EXP-NBD114 EXP-NBD196" --trim_adapters --detect_mid_strand_barcodes --enable_trim_barcodes --verbose_logs -x cuda:0 --worker_threads 8

#demux passed fastq's from fast5_fail

#first change name of each fastq:
for file in ${fail_bc_path}/pass/*; do mv -v ${file} ${file}_from_fail.fastq; done

#now demux, cant output to final dest yet bc its going to make barcode folders and dont wanna overwrite the barcode folders from before:
#require both barcodes! - OFF
#~/ont-guppy/bin/guppy_barcoder --require_barcodes_both_ends -i ${fail_bc_path}/pass -s ${fail_bc_path}/pass/demux_temp --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg" --trim_barcodes --verbose_logs -x cuda:0 --worker_threads 8

#Require one barcode!
~/ont-guppy/bin/guppy_barcoder -i ${fail_bc_path}/pass -s ${fail_bc_path}/pass/demux_temp --config configuration.cfg --barcode_kits "EXP-NBD104 EXP-NBD114 EXP-NBD196" --trim_adapters --detect_mid_strand_barcodes --enable_trim_barcodes --verbose_logs -x cuda:0 --worker_threads 8

#apparently file name gets changed back to the base ont name (no "from fail" tag) -.-
for file in ${fail_bc_path}/pass/demux_temp/barcode*/*.fastq; do mv -v ${file} ${file}_from_fail.fastq; done

#ok nowwww move output to final demultiplexed read location:
for barcode in ${fail_bc_path}/pass/demux_temp/barcode*; do bn=$(basename $barcode); for file in $barcode/*; do file_bn=$(basename $file); mv $barcode/${file_bn} ${demux_pass_all}/${bn}/; done; done

#ok, now get # reads/ barcode:
#for each barcode folder get all reads into one file, they are currently spread across many:
for barcode in ${demux_pass_all}/barcode*; do bn=$(basename $barcode); cat ${barcode}/*.fastq >> ${barcode}/${bn}.all.fastq; done


#make output dir
mkdir ${demux_pass_all}/output

#get readcount, print and pass to readcounts.txt file
wc -l ${demux_pass_all}/barcode*/*.all.fastq | awk '{print $2, $1/4}'
wc -l ${demux_pass_all}/barcode*/*.all.fastq | awk '{print $2, $1/4}' >> ${demux_pass_all}/output/readcounts.txt


#get read lengths and info

for i in ${demux_pass_all}/barcode*/*.all.fastq; do awk 'BEGIN {t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n" ,n,m,sq/n-m*m);}' ${i}; done
for i in ${demux_pass_all}/barcode*/*.all.fastq; do awk 'BEGIN {t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("total %d avg=%f stddev=%f\n" ,n,m,sq/n-m*m);}' ${i} >> ${demux_pass_all}/output/readstats.txt; done

#convert fastqs to fastas

for i in ${demux_pass_all}/barcode*/*.all.fastq; do sed -n '1~4s/^@/>/p;2~4p' $i > $i.fasta; done

#blast to get # of reads that appear to be cystovirus
for i in ${demux_pass_all}/barcode*/*.fasta; do blastn -query $i -subject ${cysto_blast_db} -outfmt '6 delim=,' -max_target_seqs 1 -max_hsps 1 >> ${demux_pass_all}/output/cysto_blast_out.csv; done

#number hsps?
echo "number cysto blast hits:"
wc -l ${demux_pass_all}/output/cysto_blast_out.csv

#for i in ${demux_pass_all}/barcode*/*.fasta; do blastn -query $i -subject ${host_blast_db} -outfmt '6 delim=,' -max_target_seqs 1 -max_hsps 1 >> ${demux_pass_all}/output/host_blast_out.csv; done

#blast to get # of reads that appear to be PP
for i in ${demux_pass_all}/barcode*/*.fasta; do blastn -query $i -subject ${host_blast_db} -outfmt '6 delim=,' -max_target_seqs 1 -max_hsps 1 >> ${demux_pass_all}/output/pp_blast_out.csv; done


#number hsps?
echo  "number pp blast hits:"
wc -l ${demux_pass_all}/output/pp_blast_out.csv

echo >  blast_counts.txt

echo  "number cysto blast hits:" > blast_counts.txt
echo wc -l ${demux_pass_all}/output/cysto_blast_out.csv > blast_counts.txt
echo  "number pp blast hits:" > blast_counts.txt
echo wc -l ${demux_pass_all}/output/pp_blast_out.csv > blast_counts.txt
