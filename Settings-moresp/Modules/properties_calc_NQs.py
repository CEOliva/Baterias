#!/usr/bin/python
# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import os
import PCET_SMD_NQs as pcet
import scscore as scs

###---Redox potentials calculations
file_dir = 'Energy_files/DFT_B3LYP_6-311+Gdp_Opt_SMD'
solv_model='SMD'         #Solvation Model (PCM o SMD)
reaction_redox=['1a-4d']       #Reaction to calculate according to PCET scheme
electronic_energy='Yes'  #Electronic energy? --> Yes/No
gibbs_energy='Yes'       #Gibbs free energy avail --> Yes/No

def RedoxPotentials(reaction_redox,electronic_energy,gibbs_energy):
    if electronic_energy == 'Yes':
       file=open('%s/rxn%s_E.txt' %(file_dir,reaction_redox),'r')  #%d for integers
       lines=file.readlines()
       elec_energies = pcet.get_energies(lines)
       file.close()
       molecules = elec_energies[0]
       energies_ox_s = elec_energies[1]
       energies_red_s = elec_energies[2]
    else:
       molecules = []
       energies_ox_s=[]
       energies_red_s=[]
       for i in molecules:
           molecules.append(np.nan)
           energies_ox_s.append(np.nan)
           energies_red_s.append(np.nan)

    if gibbs_energy == 'Yes':
       file=open('%s/rxn%s_G.txt' %(file_dir,reaction_redox),'r') 
       lines=file.readlines()
       gibbs_energies = pcet.get_energies(lines)
       file.close()
       molecules_G = gibbs_energies[0]
       gibbs_ox_s=gibbs_energies[1]
       gibbs_red_s=gibbs_energies[2]
    else:
       molecules_G = []
       gibbs_ox_s=[]
       gibbs_red_s=[]
       for i in molecules:
           molecules_G.append(np.nan)
           gibbs_ox_s.append(np.nan)
           gibbs_red_s.append(np.nan)
    #print(molecules)
    #print(molecules_G)
    print(energies_ox_s)
    print(energies_red_s)
    print(gibbs_ox_s)
    print(gibbs_red_s)

    Energies_redox_DataSet=list(zip(molecules,energies_ox_s,energies_red_s,gibbs_ox_s,gibbs_red_s))
    df_Energies_redox = pd.DataFrame(data = Energies_redox_DataSet, columns=['Molecule','energy_ox_S','energy_red_S','gibbs_ox_S','gibbs_red_S'])

    if electronic_energy == 'Yes':
       potentials_from_E = pcet.Potentials_from_E(reaction_redox,energies_ox_s,energies_red_s)
       DE_redox_S=potentials_from_E[0]
       E_calc_elec=potentials_from_E[1]
    else:
       DE_redox_S=[]
       E_calc_elec=[]
       for i in moleculas:
           DE_redox_S.append(np.nan)
           E_calc_elec.append(np.nan)
    if gibbs_energy == 'Yes':
       potentials_from_G = pcet.Potentials_from_G(reaction_redox,gibbs_ox_s,gibbs_red_s)
       DG_redox_S=potentials_from_G[0]
       E_calc_gibbs=potentials_from_G[1]
    else:
       DG_redox_S=[]
       E_calc_gibbs=[]
       for i in moleculas:
           DG_redox_S.append(np.nan)
           E_calc_gibbs.append(np.nan)

    print(E_calc_elec)
    print(E_calc_gibbs)
    return [molecules,reaction_redox,energies_ox_s,energies_red_s,gibbs_ox_s,gibbs_red_s,DE_redox_S,DG_redox_S,E_calc_elec,E_calc_gibbs]

Potential_rxn1a4d = RedoxPotentials(reaction_redox[0],electronic_energy,gibbs_energy)
print(Potential_rxn1a4d)

###---pKa calculations
reaction_pka=['A','B']
electronic_energy='Yes'  #Electronic energy? --> Yes/No
gibbs_energy='Yes'       #Gibbs free energy avail --> Yes/No

def pkaCalc(reaction_pka,electronic_energy,gibbs_energy):
    if electronic_energy == 'Yes':
       file=open('%s/rxn%s_E.txt' %(file_dir,reaction_pka),'r')  #%d for integers
       lines=file.readlines()
       elec_energies = pcet.get_energies(lines)
       file.close()
       molecules_pka = elec_energies[0]
       energies_prot_s = elec_energies[1]
       energies_deprot_s = elec_energies[2]
    else:
       molecules_pka = []
       energies_prot_s=[]
       energies_deprot_s=[]
       for i in molecules:
           molecules_pka.append(np.nan)
           energies_prot_s.append(np.nan)
           energies_deprot_s.append(np.nan)

    if gibbs_energy == 'Yes':
       file=open('%s/rxn%s_G.txt' %(file_dir,reaction_pka),'r') 
       lines=file.readlines()
       gibbs_energies = pcet.get_energies(lines)
       file.close()
       molecules_pka_G = gibbs_energies[0]
       gibbs_prot_s=gibbs_energies[1]
       gibbs_deprot_s=gibbs_energies[2]
    else:
       molecules_pka_G = []
       gibbs_prot_s=[]
       gibbs_deprot_s=[]
       for i in molecules:
           molecules_pka_G.append(np.nan)
           gibbs_prot_s.append(np.nan)
           gibbs_deprot_s.append(np.nan)

    #print(molecules_pka)
    #print(molecules_pka_G)
    #print(energies_prot_s)
    #print(gibbs_deprot_s)

    Energies_pka_DataSet=list(zip(molecules_pka,energies_prot_s,energies_deprot_s,gibbs_prot_s,gibbs_deprot_s))
    df_Energies_pka = pd.DataFrame(data = Energies_pka_DataSet, columns=['Molecule_pka','energy_prot_S','energy_deprot_S','gibbs_prot_S','gibbs_deprot_S'])

    if electronic_energy == 'Yes':
       pka_E = pcet.pka_from_E_simple(reaction_pka,energies_prot_s,energies_deprot_s)
       DE_Solution=pka_E[0]
       pka_calc_elec=pka_E[1]
    else:
       DE_Solution=[]
       pka_calc_elec=[]
       for i in molecules_pka:
           DE_Solution.append(np.nan)
           pka_calc_elec.append(np.nan)

    if gibbs_energy == 'Yes':
       pka_G = pcet.pka_from_G_simple(reaction_pka,gibbs_prot_s,gibbs_deprot_s)
       DG_Solution=pka_G[0]
       pka_calc_gibbs=pka_G[1]
    else:
       DG_Solution=[]
       pka_calc_gibbs=[]
       for i in molecules_pka:
           DG_Solution.append(np.nan)
           pka_calc_gibbs.append(np.nan)

    #print(pka_calc_elec)
    #print(pka_calc_gibbs)
    return [molecules_pka,reaction_pka,energies_prot_s,energies_deprot_s,gibbs_prot_s,gibbs_deprot_s,DE_Solution,DG_Solution,pka_calc_elec,pka_calc_gibbs]

#pka_rxnA = pkaCalc(reaction_pka[0],electronic_energy,gibbs_energy)
#print(pka_rxnB)

###---Calculating pka with chemaxon
#os.system('mkdir Info_files')
#os.system('cxcalc Coord_opt/DFT_B3LYP_6-311G+dp_Opt_SMD_NAM/smiles.smi pka -a 3 -b 3 > Info_files/chemaxon_pkas.csv')

###---Reading pkas (chemaxon and DFT) to get the ph values to calculate solubilities
#pka_chemaxon = pd.read_csv('Info_files/chemaxon_pkas.csv',delimiter='\t')
#pka_chemaxon['bpKa1'] = np.nan
#ph_chemaxon  = pd.DataFrame(data=[7.0],columns=['ph_chemaxon'])
#pka_dft = pd.DataFrame(pka_calc_gibbs)
#ph_dft = pd.DataFrame(data=[7.0],columns=['ph_dft'])
smiles_ox = pd.read_csv('Coord_opt/DFT_B3LYP_6-311+Gdp_Opt_SMD_NQ/smiles.smi',header=None,sep='\t',usecols=[0],names=['smiles_ox'])
smiles_red = pd.read_csv('Coord_opt/DFT_B3LYP_6-311+Gdp_Opt_SMD_H2NQ/smiles.smi',header=None,sep='\t',usecols=[0],names=['smiles_red'])
#smiles_HNAM = pd.read_csv('Coord_opt/DFT_B3LYP_6-311G+dp_Opt_SMD_HNAM/smiles.smi',header=None,sep='\t',usecols=[0],names=['smiles_HNAM'])

###---Solubily calculation (logs and logd) with chemaxon using chemaxon and dft pH values
#logs_params = pd.concat([smiles_NAM+,ph_chemaxon,ph_dft],axis=1)
#logs_params.to_csv('Info_files/chemaxon_logs_params.csv', index=False, header=False)
#os.system('bash ../../Modules/logs_calc.sh')
#os.system('cxcalc Coord_opt/DFT_B3LYP_6-31Gd_Opt_SMD_Diquat/smiles.smi logs -H 7.0 > Info_files/chemaxon_logs.csv')

#Calculating SCScore
#model = scs.SCScorer()
#model.restore(scs.WEIGHTS_FILE)
#scscore = [model.get_score_from_smi(m)[1] for m in smiles_NAM+['smiles']]
#smiles_NAM+['sc_score'] = scscore

###---Read de functional groups and logS files
df_func_groups = pd.read_csv('subs.csv',header=None,names=['R2','R3','R5','R6','R7','R8'])
df_func_groups['core'] = 'C9($R7)=C($R8)C8=C(C($R5)=C9($R6))C(=O)C($R3)=C($R2)C8=O'
df_func_groups['type'] = '1,4-Naphtoquinones'
#logs_chemaxon  = pd.read_csv('Info_files/chemaxon_logs.csv',header=0,sep='\t',usecols=[1],names=['logs_chemaxon'])
#logs_chemdft   = pd.read_csv('Info_files/chem-dft_logs.csv',header=0,sep='\t',usecols=[1],names=['logs_chemaxon-dft'])
#logd_chemaxon  = pd.read_csv('Info_files/chemaxon_logd.csv',header=0,sep='\t',usecols=[1],names=['logd_chemaxon'])
#logd_chemdft   = pd.read_csv('Info_files/chem-dft_logd.csv',header=0,sep='\t',usecols=[1],names=['logd_chemaxon-dft'])
#print(logs_chemaxon)

###---Generate and save the dataset
DataSet_pot_pka=list(zip(Potential_rxn1a4d[6],Potential_rxn1a4d[7],Potential_rxn1a4d[8],Potential_rxn1a4d[9]))
labels_dataset_pot_pka = ['DE_redox_1a-4d','DG_redox_1a-4d','pot_elec_1a-4d','pot_gibbs_1a-4d']

df_pot_pka=pd.DataFrame(data=DataSet_pot_pka, columns=labels_dataset_pot_pka)

df_general=pd.concat([df_func_groups,smiles_ox,smiles_red,df_pot_pka],axis=1)
labels_general = ['R2','R3','R5','R6','R7','R8','core','type','smiles_ox','smiles_red',
                  'DE_redox_1a-4d','DG_redox_1a-4d','pot_elec_1a-4d','pot_gibbs_1a-4d']

df_general.to_csv('results.csv', index=False, header=True, columns=labels_general) #Save file







