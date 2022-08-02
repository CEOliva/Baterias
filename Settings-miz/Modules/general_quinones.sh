#!/bin/bash

read -r idmol < uuidmol

#---------------------------------------------------------------------#
# Conformation analisys for Nicotinamides (NAM) and HNAM+ with obabel #
#---------------------------------------------------------------------#

# Cores:
# [He]C1=C([Ne])C2=C(C(Br)=C1I)C(=O)C(Cl)=C(F)C2=O --> 1,4-naphtoquinone core
# C9($R7)=C($R8)C8=C(C($R5)=C9($R6))C(=O)C($R3)=C($R2)C8=O --> Substitutions
# C9(CC)=C(O)C8=C(C(CN)=C9(N))C(=O)C(C)=C(CO)C8=O --> Example
# C9($R7)=C($R8)C8=C(C(O)=C($R3)C($R2)=C8(O))C($R5)=C9($R6) --> 1,4-hidroquinone substitutions
# C9(CC)=C(O)C8=C(C(O)=C(C)C(CO)=C8(O))C(CN)=C9(N) --> Examples

start_time=$(date +%s)
filename="subs.csv"
dir="Conformers"

mkdir Info_files
while IFS=',' read -r R2 R3 R5 R6 R7 R8 ; do
      #Conformer search and minimizations for BQ
      mol1="NQ"
      mkdir -p ${dir}/${mol1}
      cat > ${dir}/${mol1}/smiles.smi << !
C9($R7)=C($R8)C8=C(C($R5)=C9($R6))C(=O)C($R3)=C($R2)C8=O
!
      obabel ${dir}/${mol1}/smiles.smi -O ${dir}/${mol1}/coord_init.xyz -h -nowarn --gen3D
      echo "${idmol} ${mol1} conformer search..."
      obconformer 50 500 ${dir}/${mol1}/coord_init.xyz > ${dir}/${mol1}/coord_conf.xyz 2> ${dir}/${mol1}/info_conf.txt 
      sed -i '/WARNING/d' ${dir}/${mol1}/coord_conf.xyz
      echo "${idmol} ${mol1} energy minimization..."
      obminimize -ff mmff94 -sd -n 10000 -ixyz ${dir}/${mol1}/coord_conf.xyz -oxyz > ${dir}/${mol1}/coord_opt.xyz 2> ${dir}/${mol1}/info_min.txt
      pos=$(awk -F+ '{print NF-1}' < ${dir}/${mol1}/smiles.smi)
      neg=$(awk -F- '{print NF-1}' < ${dir}/${mol1}/smiles.smi)
      mol_charge=$(($pos - $neg))
      echo "$mol_charge" > ${dir}/${mol1}/mol_charge.txt

      #Conformer search and minimizations for H2BQ
      mol1="H2NQ"
      mkdir -p ${dir}/${mol1}
      cat > ${dir}/${mol1}/smiles.smi << !
C9($R7)=C($R8)C8=C(C(O)=C($R3)C($R2)=C8(O))C($R5)=C9($R6)
!
      obabel ${dir}/${mol1}/smiles.smi -O ${dir}/${mol1}/coord_init.xyz -h -nowarn --gen3D
      echo "${idmol} ${mol1} conformer search..."
      obconformer 50 500 ${dir}/${mol1}/coord_init.xyz > ${dir}/${mol1}/coord_conf.xyz 2> ${dir}/${mol1}/info_conf.txt
      sed -i '/WARNING/d' ${dir}/${mol1}/coord_conf.xyz
      echo "${idmol} ${mol1} energy minimization..."
      obminimize -ff mmff94 -sd -n 10000 -ixyz ${dir}/${mol1}/coord_conf.xyz -oxyz > ${dir}/${mol1}/coord_opt.xyz 2> ${dir}/${mol1}/info_min.txt
      pos=$(awk -F+ '{print NF-1}' < ${dir}/${mol1}/smiles.smi)
      neg=$(awk -F- '{print NF-1}' < ${dir}/${mol1}/smiles.smi)
      mol_charge=$(($pos - $neg))
      echo "$mol_charge" > ${dir}/${mol1}/mol_charge.txt
done < "$filename"

finish_time=$(date +%s)
echo "Time duration: $((finish_time - start_time)) secs." > Info_files/time_conf.txt

#-----------------------------------------------------------------------------#
# General settings for the input files in Gaussian for NQs                    #
#-----------------------------------------------------------------------------#

Functional="B3LYP"       #Options=PM7,B3LYP,wB97XD,... --> DFT functional  
Base="6-311+G(d,p)"          #Basis set accorging to gaussian options
SCRF="SMD"               #Options=PCM,SMD,... --> Solvation model. Use SCRF="" for gas phase calculations. 
Solvent="water"          #Options=water,DiMethylSulfoxide,Acetonitrile,... --> Solvent
SP_Opt="Opt"             #Options=SP,Opt --> Single Point (SP) or Optimization (Opt)
memory="64"              #Memory for the input file in GB (32Gb for q_residual, 64Gb for q_htc, 240 Gb for q_16_256g)
queue="htc"

#Other general settings...
if [[ "${SCRF}" = "" ]]; then
   method="${Functional}_${Base//[,()]/}_${SP_Opt}"
else
   method="${Functional}_${Base//[,()]/}_${SP_Opt}_${SCRF}"
fi

if [[ "$SP_Opt" = "Opt" ]] ; then
   SP_Opt2="Opt=(maxcyc=800)"  #calcFC
else
   SP_Opt2=""
fi

IodineFoot() {
#Function to read the smiles and generate a string with the atoms for the input foot different from I.
if [ `echo $smiles | grep "O" | wc -l` -ge 1 ] ; then
    O_atom=" O"
else
    O_atom=""
fi
if [ `echo $smiles | grep "S" | wc -l` -ge 1 ] ; then
   S_atom=" S"
else
   S_atom=""
fi
if [ `echo $smiles | grep "F" | wc -l` -ge 1 ] ; then
   F_atom=" F"
else
   F_atom=""
fi
if [ `echo $smiles | grep "Cl" | wc -l` -ge 1 ] ; then
   Cl_atom=" Cl"
else
   Cl_atom=""
fi
if [ `echo $smiles | grep "Br" | wc -l` -ge 1 ] ; then
   Br_atom=" Br"
else
   Br_atom=""
fi
if [ `echo $smiles | grep "Si" | wc -l` -ge 1 ] ; then
   Si_atom=" Si"
else
   Si_atom=""
fi
string="C H N$O_atom$S_atom$F_atom$Cl_atom$Br_atom$Si_atom 0"
}

CheckJobs() {
#function to check for jobs done
job_run=`bjobs $job_id | grep -c 'RUN'`
job_pend=`bjobs $job_id | grep -c 'PEND'`
while [[ ${job_run} = "1" || ${job_pend} = "1" ]] ; do
      #if [[ ${job_run} = "1" ]] ; then
      #   echo "Job $job_id RUN"
      #elif [[ ${job_pend} = "1" ]] ; then
      #   echo "Job $job_id PEND"
      #fi
      sleep 5
      job_run=`bjobs $job_id | grep -c 'RUN'`
      job_pend=`bjobs $job_id | grep -c 'PEND'`
done
echo "${idmol} ${Molecule}: Job ${job_id} done!"
}

BsubBatch() {
#function to submit jobs in batches of 10 
jobs_status=0
while [ "$jobs_status" -eq "0" ] ; do
      jobs_num=`bjobs | grep -c 'mflol_a'`
      if [ "$jobs_num" -lt "10" ] ; then
         cd ${dir_writing}/${Molecule}
         job_id=`bsub < g16script | cut -d \< -f 2 | cut -d \> -f 1`  #Mizlti
         #g16 < $j/str.gjf > $j/str.log  #Erikuchis
         cd ../../
         sleep 5
         jobs_status=1
      else
         sleep 30
         jobs_status=0
      fi
done
}

HeaderFootScript() {
#A header and a foot are defined for the input file
if [[ "$SCRF" = "PCM" || "$SCRF" = "SMD" ]]; then
   if [ $iodine_num -ge 1 ] ; then
      #Header for iodine in solution
      cat > ${dir_writing}/${Molecule}/header << !
%chk=str.chk
%mem=${memory}Gb
%nproc=16
# ${Functional}/Gen ${SP_Opt2} scf=xqc scrf=(${SCRF},solvent=${Solvent}) Freq output=wfn 

Title

${Charge} ${Mult}
!
   else
      #Header for non-iodine in solution
      cat > ${dir_writing}/${Molecule}/header << !
%chk=str.chk
%mem=${memory}Gb
%nproc=16
# ${Functional}/${Base} ${SP_Opt2} scf=xqc scrf=(${SCRF},solvent=${Solvent}) Freq output=wfn 

Title

${Charge} ${Mult}
!
   fi
elif [[ "$SCRF" = "" ]] ; then
   if [ $iodine_num -ge 1 ] ; then
      #Header for iodine in gas phase
      cat > ${dir_writing}/${Molecule}/header << !
%chk=str.chk
%mem=${memory}Gb
%nproc=16
# ${Functional}/Gen ${SP_Opt2} scf=xqc Freq output=wfn 

Title

${Charge} ${Mult}
!
   else
      #Header for non-iodine in gas phase
      cat > ${dir_writing}/${Molecule}/header << !
%chk=str.chk
%mem=${memory}Gb
%nproc=16
# ${Functional}/${Base} ${SP_Opt2} scf=xqc Freq output=wfn

Title

${Charge} ${Mult}
!
   fi
fi
#Foot for the input file
if [ $iodine_num -ge 1 ] ; then
   IodineFoot
   #foot for iodine
   cat > ${dir_writing}/${Molecule}/foot << !

$string
$Base
****
I 0
3-21G
****

str.wfn

!
else
   #regular foot
   cat > ${dir_writing}/${Molecule}/foot << !

str.wfn

!
fi
#Script for gaussian16
cat > ${dir_writing}/${Molecule}/g16script << !
#!/bin/bash
#BSUB -oo str.o
#BSUB -eo str.e
#BSUB -q q_htc
#BSUB -n 16
#BSUB -R "span[hosts=1]"

module load gaussian/2016.03
export GAUSS_SCRDIR=$TMPU/g16scrdir

g16 < str.gjf > str.log &

!
sed '1,2d' ${dir_writing}/${Molecule}/coord_init.xyz > ${dir_writing}/${Molecule}/coord #Removes the two first rows in the .xyz file
cat ${dir_writing}/${Molecule}/header ${dir_writing}/${Molecule}/coord ${dir_writing}/${Molecule}/foot > ${dir_writing}/${Molecule}/str.gjf
}

CheckErrors() {
#looking for common gaussian errors!!!
galloc_error=`grep -c "galloc" ${dir_writing}/${Molecule}/str.e`
opt_error=`grep -c " Optimized" ${dir_writing}/${Molecule}/str.log`
if [[ "$galloc_error" -ge "0" ]] && [[ "$opt_error" -eq "0" ]] ; then
   #memory error --> galloc:  could not allocate memory
   #solution --> decrease the amount of memory and submit the calculation again
   if [ "$galloc_error" -gt "0" ] ; then
      echo "${Molecule}: galloc:  could not allocate memory" >> Info_files/error_${method}.txt
      rm ${dir_writing}/${Molecule}/core.*
      mem_new=$((memory-32))
      sed -i "s/%mem=${memory}Gb/%mem=${mem_new}Gb/g" ${dir_writing}/${Molecule}/str.gjf
   elif [[ "$galloc_error" -eq "0" ]] && [[ "$opt_error" -eq "0" ]] ; then
      nonopt_error=`grep -c "Non-Optimized" ${dir_writing}/${Molecule}/str.log`
      #error --> Non-Optimized Parameters
      #solution --> restart the calculation from the str.chk file 
      if [ "$nonopt_error" -gt "0" ] ; then
         echo "${Molecule}: Non-Optimized Parameters" >> Info_files/error_${method}.txt
         rm ${dir_writing}/${Molecule}/core.*
         cp ${dir_writing}/${Molecule}/str.gjf ${dir_writing}/${Molecule}/errstr.gjf
         sed -i '9,/^$/d' ${dir_writing}/${Molecule}/str.gjf
         sed -i '9i\\' ${dir_writing}/${Molecule}/str.gjf
         sed -i "s/${SP_Opt2}/${SP_Opt2} Geom=Checkpoint/g" ${dir_writing}/${Molecule}/str.gjf
      else
         #error --> some other error
         #solution --> send the calculation again
         echo "${Molecule}: Unknown error" >> Info_files/error_${method}.txt
         rm ${dir_writing}/${Molecule}/core.*
      fi
   fi
   #re-submit calculation
   cd ${dir_writing}/${Molecule}
   job_id=`bsub < g16script | cut -d \< -f 2 | cut -d \> -f 1`  #Mizlti
   cd ../../
   echo "${idmol} ${Molecule}: Job ${job_id} re-submitted..."
   echo "${Molecule}: calculation re-submitted" >> Info_files/error_${method}.txt
   sleep 5
elif [[ "$galloc_error" -eq "0" ]] && [[ "$opt_error" -gt "0" ]] ; then
   #error --> Optimized Parameters but no Free Energies
   gibbs=`grep -c "Free Energies" ${dir_writing}/${Molecule}/str.log`
   if [ ${gibbs} -eq "0" ] ; then
      echo "${Molecule}: Optimized Parameters but no Free Energies" >> Info_files/error_${method}.txt
      #solution1 --> restart freq calculations from str.rwf file if exists
      rwf_file=`ls ${dir_writing}/${Molecule} | grep -c "str.rwf"`
      if [[ ${rwf_file} -gt "0" ]] ; then
         echo "Frequencies calculation restarted from str.rwf" >> Info_files/error_${method}.txt
         cat > ${dir_writing}/${Molecule}/strfreq.gjf << !
%rwf=str.rwf
%NoSave
%chk=str.chk
%mem=${memory}Gb
%nproc=16
# Restart 

!
         sed -i "s/str/strfreq/g" ${dir_writing}/${Molecule}/g16script
      #solution2 --> extract the optimized coord and resubmit the freq calculations
      else
         echo "Frequencies calculation resubmitted from optimized geometry" >> Info_files/error_${method}.txt
         obabel -ilog ${dir_writing}/${Molecule}/str.log -oxyz > ${dir_writing}/${Molecule}/coord_opt.xyz
         sed '1,2d' ${dir_writing}/${Molecule}/coord_opt.xyz > ${dir_writing}/${Molecule}/coord
         cat ${dir_writing}/${Molecule}/header ${dir_writing}/${Molecule}/coord ${dir_writing}/${Molecule}/foot > ${dir_writing}/${Molecule}/strfreq.gjf
         sed -i "s/str/strfreq/g" ${dir_writing}/${Molecule}/g16script
      fi
      #re-submit calculation
      cd ${dir_writing}/${Molecule}
      job_id=`bsub < g16script | cut -d \< -f 2 | cut -d \> -f 1`  #Mizlti
      cd ../../
      echo "${idmol} ${Molecule}: Job ${job_id} re-submitted..."
      echo "${Molecule}: calculation re-submitted" >> Info_files/error_${method}.txt
      sleep 5
   else
      #Optimized Parameters and Free Energies
      echo "${Molecule}: Optimized Parameters and Free Energies" >> Info_files/error_${method}.txt
   fi
fi
}

#----------------------------------------------------------------------------------#
# This section creates the input files in Gaussian and sumbits the jobs for NQ     #
#----------------------------------------------------------------------------------#

### Additional settings for the gaussian input file for NQ  ###
dir_reading="Conformers" #change the reading file if necessary!
Molecule="NQ"            #Label of the molecule (e.g. Bpy, HBpy+, H2Bpy2+,...)
Mult="1"                 #Multiplicity

dir_writing="DFT_${method}"
start_time=$(date +%s)

### Parsing the data from Conformers ###
mkdir -p ${dir_writing}/${Molecule}
cp ${dir_reading}/${Molecule}/coord_opt.xyz ${dir_writing}/${Molecule}/coord_init.xyz
cp ${dir_reading}/${Molecule}/smiles.smi ${dir_writing}/${Molecule}/smiles.smi
cp ${dir_reading}/${Molecule}/mol_charge.txt ${dir_writing}/${Molecule}/mol_charge.txt

### This section creates the input files and submits the jobs ###
read -r Charge < ${dir_writing}/${Molecule}/mol_charge.txt
read -r smiles < ${dir_writing}/${Molecule}/smiles.smi
iodine_num=`echo $smiles | grep "I" | wc -l`

HeaderFootScript

#Submit the jobs... 
cd ${dir_writing}/${Molecule}
job_id=`bsub < g16script | cut -d \< -f 2 | cut -d \> -f 1`  #Mizlti
cd ../../
echo "${idmol} ${Molecule}: Job ${job_id} submitted..."
sleep 5

### This section checks for failed jobs and jobs done before send the next molecule ###
CheckJobs
CheckErrors
CheckJobs

#remove temporal files head,foot,coord and extract coord_opt.yxz
rm ${dir_writing}/${Molecule}/header ${dir_writing}/${Molecule}/foot ${dir_writing}/${Molecule}/coord
obabel -ilog ${dir_writing}/${Molecule}/str.log -oxyz > ${dir_writing}/${Molecule}/coord_opt.xyz
obabel -ilog ${dir_writing}/${Molecule}/str.log -omol > ${dir_writing}/${Molecule}/coord_opt.mol

#Time estimation of the calculation for NAM+
finish_time=$(date +%s)
echo "${Molecule} calculation time: $((finish_time - start_time)) secs." >> Info_files/time_${method}.txt

#-----------------------------------------------------------------------------------#
# This section creates the input files in Gaussian and sumbit the jobs for H2NQ     #
#-----------------------------------------------------------------------------------#

### Additional settings for the gaussian input file for H2NQ  ###
dir_reading="Conformers" #change the reading file if necessary!
Molecule="H2NQ"          #Label of the molecule (e.g. Bpy, HBpy+, H2Bpy2+,...) 
Mult="1"                 #Multiplicity

dir_writing="DFT_${method}"
start_time=$(date +%s)

### Parsing the data from DFT-NQ ### 
mkdir -p ${dir_writing}/${Molecule}
cp ${dir_reading}/${Molecule}/coord_opt.xyz ${dir_writing}/${Molecule}/coord_init.xyz
cp ${dir_reading}/${Molecule}/smiles.smi ${dir_writing}/${Molecule}/smiles.smi
cp ${dir_reading}/${Molecule}/mol_charge.txt ${dir_writing}/${Molecule}/mol_charge.txt

### This section creates the input files and submits the jobs ###
#read -r q < ${dir_writing}/${Molecule}/mol_charge.txt
#Charge=$(($q - 1))
#echo "$Charge" > ${dir_writing}/${Molecule}/mol_charge.txt
read -r Charge < ${dir_writing}/${Molecule}/mol_charge.txt
read -r smiles < ${dir_writing}/${Molecule}/smiles.smi
iodine_num=`echo $smiles | grep "I" | wc -l`

HeaderFootScript

#Submit the jobs...
cd ${dir_writing}/${Molecule}
job_id=`bsub < g16script | cut -d \< -f 2 | cut -d \> -f 1`  #Mizlti
cd ../../
echo "${idmol} ${Molecule}: Job ${job_id} submitted..."
sleep 5

### This section checks for failed jobs and jobs done before send the next molecule ###
CheckJobs
CheckErrors
CheckJobs

#remove temporal files head,foot,coord and extract coord_opt.yxz
rm ${dir_writing}/${Molecule}/header ${dir_writing}/${Molecule}/foot ${dir_writing}/${Molecule}/coord
obabel -ilog ${dir_writing}/${Molecule}/str.log -oxyz > ${dir_writing}/${Molecule}/coord_opt.xyz
obabel -ilog ${dir_writing}/${Molecule}/str.log -omol > ${dir_writing}/${Molecule}/coord_opt.mol

#Time estimation of the calculation for NAM_r
finish_time=$(date +%s)
echo "${Molecule} calculation time: $((finish_time - start_time)) secs." >> Info_files/time_${method}.txt

#---------------------------------------------------------------------------------------#
# This section reads the energies from the output files for reaction 1a-4d: NQ --> H2NQ #
#---------------------------------------------------------------------------------------#
### Settings to read energies for reaction 1a-4d ###
Reaction="1a-4d"
Reactant="NQ"  #Protonated/Oxidized
Product="H2NQ"    #Deprotonated/Reduced
dir_reading="DFT_${method}"  #Change the reading directory!

rm Energy_files/${dir_reading}/rxn${Reaction}_*.txt
rm Info_files/info_${dir_reading}_rxn${Reaction}.txt
rm -r Coord_opt/${dir_reading}_${Reactant}
rm -r Coord_opt/${dir_reading}_${Product}
mkdir -p Energy_files/${dir_reading}
mkdir -p Coord_opt/${dir_reading}_${Reactant}
mkdir -p Coord_opt/${dir_reading}_${Product}

EnergiesReactants() {
#Energies for the reactants
error=`grep -c " Optimized" ${dir_reading}/${Reactant}/str.log`    # Optimized/Normal termination
obabel -ilog ${dir_reading}/${Reactant}/str.log -oxyz > ${dir_reading}/${Reactant}/coord_opt.xyz
obabel -ilog ${dir_reading}/${Reactant}/str.log -omol > ${dir_reading}/${Reactant}/coord_opt.mol
cp ${dir_reading}/${Reactant}/coord_opt.* Coord_opt/${dir_reading}_${Reactant}
cp ${dir_reading}/${Reactant}/smiles.smi Coord_opt/${dir_reading}_${Reactant}
if [[ ${error} != 0 ]] ; then
   echo "${dir_reading}/${Reactant}/str.log Optimized Parameters!" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
   grep " SCF Done" ${dir_reading}/${Reactant}/str.log | tail -n 1 >> Energy_files/${dir_reading}/rxn${Reaction}_E.txt
   gibbs=`grep -c "Free Energies" ${dir_reading}/${Reactant}/str.log`
   if [[ ${gibbs} != 0 ]] ; then
      echo "${dir_reading}/${Reactant}/str.log Free Energies done!" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
      grep "Free Energies" ${dir_reading}/${Reactant}/str.log >> Energy_files/${dir_reading}/rxn${Reaction}_G.txt
   else
      echo "Error. Free Energies not found! ${dir_reading}/${Reactant}/str.log" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
      echo " Sum of electronic and thermal Free Energies= NaN" >> Energy_files/${dir_reading}/rxn${Reaction}_G.txt
   fi
else
   echo "Error. Optimized Parameters not found! ${dir_reading}/${Reactant}/str.log" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
   echo " SCF Done:  E(RB3LYP) = NaN" >> Energy_files/${dir_reading}/rxn${Reaction}_E.txt
   echo " Sum of electronic and thermal Free Energies= NaN" >> Energy_files/${dir_reading}/rxn${Reaction}_G.txt
fi
}
EnergiesReactants

EnergiesProducts() {
#Energies for the products 
error=`grep -c " Optimized" ${dir_reading}/${Product}/str.log`    # Optimized/Normal termination
obabel -ilog ${dir_reading}/${Product}/str.log -oxyz > ${dir_reading}/${Product}/coord_opt.xyz
obabel -ilog ${dir_reading}/${Product}/str.log -omol > ${dir_reading}/${Product}/coord_opt.mol
cp ${dir_reading}/${Product}/coord_opt.* Coord_opt/${dir_reading}_${Product}
cp ${dir_reading}/${Product}/smiles.smi Coord_opt/${dir_reading}_${Product}
if [[ ${error} != 0 ]] ; then
   echo "${dir_reading}/${Product}/str.log Optimized Parameters!" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
   grep " SCF Done" ${dir_reading}/${Product}/str.log | tail -n 1 >> Energy_files/${dir_reading}/rxn${Reaction}_E.txt
   gibbs=`grep -c "Free Energies" ${dir_reading}/${Product}/str.log`
   if [[ ${gibbs} != 0 ]] ; then
      echo "${dir_reading}/${Product}/str.log Free Energies done!" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
      grep "Free Energies" ${dir_reading}/${Product}/str.log >> Energy_files/${dir_reading}/rxn${Reaction}_G.txt
   else
      echo "Error. Free Energies not found! ${dir_reading}/${Product}/str.log" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
      echo " Sum of electronic and thermal Free Energies= NaN" >> Energy_files/${dir_reading}/rxn${Reaction}_G.txt
   fi
else
   echo "Error. Optimized Parameters not found! ${dir_reading}/${Product}/str.log" >> Info_files/info_${dir_reading}_rxn${Reaction}.txt
   echo " SCF Done:  E(RB3LYP) = NaN" >> Energy_files/${dir_reading}/rxn${Reaction}_E.txt
   echo " Sum of electronic and thermal Free Energies= NaN" >> Energy_files/${dir_reading}/rxn${Reaction}_G.txt
fi
}
EnergiesProducts


echo  "$idmol: dft calculations completed!"

