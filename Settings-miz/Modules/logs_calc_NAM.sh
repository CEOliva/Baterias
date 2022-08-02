#!/bin/bash

mkdir Info_files
rm Info_files/chemaxon_logs.csv
i=1
while IFS="," read -r smi1 smi2 ph1 ph2; do
      #echo "$i"
      #echo "$ph"
      #cxcalc $smi logs -H $ph1 logs -H $ph2 >> Info_files/chemaxon_logs.csv
      cxcalc $smi1 logs -H $ph1 > Info_files/chemaxon_logs.csv
      cxcalc $smi1 logs -H $ph2 > Info_files/chem-dft_logs.csv
      cxcalc $smi1 logd -H $ph1 > Info_files/chemaxon_logd.csv
      cxcalc $smi1 logd -H $ph2 > Info_files/chem-dft_logd.csv
      #sed -i '/logS/d' Info_files/chemaxon_logs.csv 
      i=$((i+1))
done < "Info_files/chemaxon_logs_params.csv"



