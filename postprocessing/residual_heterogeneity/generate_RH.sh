#!/bin/bash
set -e

#This wrapper script will take a maps2 output file and generate maps2 files containing only RH and non-RH mutations, along with all of the intermediate files.
#This should be run in the directory containing the files. Lengths file must be absolute path or in the directory, do not give a relative path
#Usage:bash -x generate_RH.sh [glob for maps2 files to process, e.g. "*.deduped.tsv"] reference_length_file.length
#Don't forget quotes on the file glob!
#This can also be run with single files, no quotes required

if [ "$#" -ne 2 ]; then
	echo "Usage: bash -x $0 MAPSfile.tsv lengthfile.tsv" >&2
	exit 1
fi

#path to scripts
path="./"  ##Modify for path to folder
files=$1
f="${files%.*}"
extension="${files##*.}"

mkdir -p intermediate_RH_files
#Correct Het to Hom
if [ ${extension} == "gz" ] 
	then
		for i in ${files}; do gzcat ${i} > "${i%.*}"; done
		files=$f
		f="${f%.*}"
		for i in ${files}; do awk -v OFS='\t' '{if ( $9 ~ /[0-9]/ && $9 != 0 && $9/($9+$10) < 0.15) {$8 = "hom"; print $0, "Warning: corrected to hom";} else print $0}' ${i} > intermediate_RH_files/"${i%.*}".corrected.tsv; done

	else
		for i in ${files}; do awk -v OFS='\t' '{if ( $9 ~ /[0-9]/ && $9 != 0 && $9/($9+$10) < 0.15) {$8 = "hom"; print $0, "Warning: corrected to hom";} else print $0}' ${i} > intermediate_RH_files/"${i%.*}".corrected.tsv; done
fi

cd intermediate_RH_files
ln -s ../${2} .

#Generate binned MAPS calls and lengths file
for i in ${f}.corrected.tsv; do python ${path}/binMAPSdata.py -m ${i} -b 10000 -l ${2} -o "${i%.*}".10kb_bins.tsv -L "${i%.*}".10kb_bins.length; done

#Generate byContig files
for i in ${f}.corrected.10kb_bins.tsv; do python ${path}/detectHeteroZonesV1.5.py -t 1 -p 1 -m ${i} -l "${i%.*}".length "${i%.*}".RH ; done

#Add MI to byContig files
for i in ${f}.corrected.10kb_bins.tsv; do python ~/bin/tilling_project/processing_scripts/addMItoRHv2.py -m ${i} -r "${i%.*}".RH.byContig.tsv -z -o "${i%.*}".RH.byContig.MI.tsv; done

#Call RH regions
for i in ${f}.corrected.10kb_bins.RH.byContig.MI.tsv; do python ${path}/callRHcontigs.py -i ${i} -r "${i%.*}".RH_only.tsv -n "${i%.*}".No_RH.tsv; done

#Get MAPS outputs
for i in ${f}.corrected.10kb_bins.tsv; do awk 'NR==FNR {a[$6","$1];next} ($1","$7 in a) {print $0}' "${i%.*}".RH.byContig.MI.RH_only.tsv ${i} > "${i%.*}".RH.byContig.MI.RH_only.maps.tsv; done
for i in ${f}.corrected.10kb_bins.tsv; do awk 'NR==FNR {a[$6","$1];next} ($1","$7 in a) {print $0}' "${i%.*}".RH.byContig.MI.No_RH.tsv ${i} > "${i%.*}".RH.byContig.MI.No_RH.maps.tsv; done

#Add header back to maps format files and remove bin information
for i in ${f}.corrected.10kb_bins.tsv; do cat <(head -n 1 ${i}) <(sed 's/_bin_[0-9]*-[0-9]*\t/\t/' "${i%.*}".RH.byContig.MI.RH_only.maps.tsv) > "${i%.*}".RH.byContig.MI.RH_only.maps.tsv.tmp > "${i%.*}".RH.byContig.MI.RH_only.maps.tsv.tmp; mv "${i%.*}".RH.byContig.MI.RH_only.maps.tsv.tmp ../"${i%.*}".RH.byContig.MI.RH_only.maps.tsv; done
for i in ${f}.corrected.10kb_bins.tsv; do cat <(head -n 1 ${i}) <(sed 's/_bin_[0-9]*-[0-9]*\t/\t/' "${i%.*}".RH.byContig.MI.No_RH.maps.tsv) > "${i%.*}".RH.byContig.MI.No_RH.maps.tsv.tmp; mv "${i%.*}".RH.byContig.MI.No_RH.maps.tsv.tmp ../"${i%.*}".RH.byContig.MI.No_RH.maps.tsv; done

#Remove all intermediate files
rm lowcontigs.forEnsembl.bed

#cd ..
#rm -r intermediate_RH_files
