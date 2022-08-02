#!/usr/bin/python
# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import os
import PCET_SMD_NAM as pcet
import scscore as scs

###---Redox potentials calculations
file_dir = 'Energy_files/DFT_B3LYP_6-311+Gdp_Opt_SMD'
solv_model='SMD'         #Solvation Model (PCM o SMD)
reaction_redox=['1','2','3','1a','2b']       #Reaction to calculate according to PCET scheme
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
    #print(energies_ox_s)
    #print(gibbs_ox_s)

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

    #print(E_calc_elec)
    #print(E_calc_gibbs)
    return [molecules,reaction_redox,energies_ox_s,energies_red_s,gibbs_ox_s,gibbs_red_s,DE_redox_S,DG_redox_S,E_calc_elec,E_calc_gibbs]

Potential_rxn1 = RedoxPotentials(reaction_redox[0],electronic_energy,gibbs_energy)
Potential_rxn2 = RedoxPotentials(reaction_redox[1],electronic_energy,gibbs_energy)
Potential_rxn3 = RedoxPotentials(reaction_redox[2],electronic_energy,gibbs_energy)
Potential_rxn1a = RedoxPotentials(reaction_redox[3],electronic_energy,gibbs_energy)
Potential_rxn2b = RedoxPotentials(reaction_redox[4],electronic_energy,gibbs_energy)
print(Potential_rxn1)
print(Potential_rxn2)
print(Potential_rxn3)
print(Potential_rxn1a)
print(Potential_rxn2b)

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
       pka_E = pcet.pka_from_E_ref(reaction_pka,energies_prot_s,energies_deprot_s)
       DE_Solution_simple=pka_E[0]
       pka_calc_elec_simple=pka_E[1]
       DE_Solution_ref=pka_E[2]
       pka_calc_elec_ref=pka_E[3]
    else:
       DE_Solution_simple=[np.nan for i in molecules_pka]
       pka_calc_elec_simple=[np.nan for i in molecules_pka]
       DE_Solution_ref=[np.nan for i in molecules_pka]
       pka_calc_elec_ref=[np.nan for i in molecules_pka]

    if gibbs_energy == 'Yes':
       pka_G = pcet.pka_from_G_ref(reaction_pka,gibbs_prot_s,gibbs_deprot_s)
       DG_Solution_simple=pka_G[0]
       pka_calc_gibbs_simple=pka_G[1]
       DG_Solution_ref=pka_G[2]
       pka_calc_gibbs_ref=pka_G[3]
    else:
       DG_Solution_simple=[np.nan for i in molecules_pka]
       pka_calc_gibbs_simple=[np.nan for i in molecules_pka]
       DG_Solution_ref=[np.nan for i in molecules_pka]
       pka_calc_gibbs_ref=[np.nan for i in molecules_pka]

    #print(pka_calc_elec)
    #print(pka_calc_gibbs)
    return [molecules_pka,reaction_pka,energies_prot_s,energies_deprot_s,gibbs_prot_s,gibbs_deprot_s,
            DE_Solution_simple,DG_Solution_simple,pka_calc_elec_simple,pka_calc_gibbs_simple,
            DE_Solution_ref,DG_Solution_ref,pka_calc_elec_ref,pka_calc_gibbs_ref]

pka_rxnA = pkaCalc(reaction_pka[0],electronic_energy,gibbs_energy)
pka_rxnB = pkaCalc(reaction_pka[1],electronic_energy,gibbs_energy)
print(pka_rxnA)
print(pka_rxnB)

###---Calculating pka with chemaxon
#os.system('mkdir Info_files')
#os.system('cxcalc Coord_opt/DFT_B3LYP_6-311G+dp_Opt_SMD_NAM/smiles.smi pka -a 3 -b 3 > Info_files/chemaxon_pkas.csv')

###---Reading pkas (chemaxon and DFT) to get the ph values to calculate solubilities
#pka_chemaxon = pd.read_csv('Info_files/chemaxon_pkas.csv',delimiter='\t')
#pka_chemaxon['bpKa1'] = np.nan
#ph_chemaxon  = pd.DataFrame(data=[7.0],columns=['ph_chemaxon'])
#pka_dft = pd.DataFrame(pka_calc_gibbs)
#ph_dft = pd.DataFrame(data=[7.0],columns=['ph_dft'])
smiles_NAM = pd.read_csv('Coord_opt/DFT_B3LYP_6-311+Gdp_Opt_SMD_NAM+/smiles.smi',header=None,sep='\t',usecols=[0],names=['smiles_NAM+'])
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
df_func_groups = pd.read_csv('subs.csv',header=None,names=['R1','R2','R3','R4','R5','R6'])
df_func_groups['core'] = 'C9($R2)=C($R3)C($R4)=C($R5)C($R6)=[N+]9($R1)'
df_func_groups['type'] = 'Nicotinamides'
#logs_chemaxon  = pd.read_csv('Info_files/chemaxon_logs.csv',header=0,sep='\t',usecols=[1],names=['logs_chemaxon'])
#logs_chemdft   = pd.read_csv('Info_files/chem-dft_logs.csv',header=0,sep='\t',usecols=[1],names=['logs_chemaxon-dft'])
#logd_chemaxon  = pd.read_csv('Info_files/chemaxon_logd.csv',header=0,sep='\t',usecols=[1],names=['logd_chemaxon'])
#logd_chemdft   = pd.read_csv('Info_files/chem-dft_logd.csv',header=0,sep='\t',usecols=[1],names=['logd_chemaxon-dft'])
#print(logs_chemaxon)

###---Generate and save the dataset
DataSet_pot_pka=list(zip(Potential_rxn1[6],Potential_rxn1[7],Potential_rxn1[8],Potential_rxn1[9],
                         Potential_rxn2[6],Potential_rxn2[7],Potential_rxn2[8],Potential_rxn2[9],
                         Potential_rxn3[6],Potential_rxn3[7],Potential_rxn3[8],Potential_rxn3[9],
                         Potential_rxn1a[6],Potential_rxn1a[7],Potential_rxn1a[8],Potential_rxn1a[9],
                         Potential_rxn2b[6],Potential_rxn2b[7],Potential_rxn2b[8],Potential_rxn2b[9],
                         pka_rxnA[6],pka_rxnA[7],pka_rxnA[8],pka_rxnA[9],
                         pka_rxnB[6],pka_rxnB[7],pka_rxnB[8],pka_rxnB[9],
                         pka_rxnB[10],pka_rxnB[11],pka_rxnB[12],pka_rxnB[13]))
labels_dataset_pot_pka = ['DE_redox_1','DG_redox_1','pot_elec_1','pot_gibbs_1',
                          'DE_redox_2','DG_redox_2','pot_elec_2','pot_gibbs_2',
                          'DE_redox_3','DG_redox_3','pot_elec_3','pot_gibbs_3',
                          'DE_redox_1a','DG_redox_1a','pot_elec_1a','pot_gibbs_1a',
                          'DE_redox_2b','DG_redox_2b','pot_elec_2b','pot_gibbs_2b',
                          'DE_rxn_A','DG_rxn_A','pka_elec_A','pka_gibbs_A',
                          'DE_rxn_B','DG_rxn_B','pka_elec_B','pka_gibbs_B',
                          'DE_rxn_B_ref','DG_rxn_B_ref','pka_elec_B_ref','pka_gibbs_B_ref']
df_pot_pka=pd.DataFrame(data=DataSet_pot_pka, columns=labels_dataset_pot_pka)

df_general=pd.concat([df_func_groups,smiles_NAM,df_pot_pka],axis=1)
labels_general = ['R1','R2','R3','R4','R5','R6','core','type','smiles_NAM+',
                  'DE_redox_1','DG_redox_1','pot_elec_1','pot_gibbs_1',
                  'DE_redox_2','DG_redox_2','pot_elec_2','pot_gibbs_2',
                  'DE_redox_3','DG_redox_3','pot_elec_3','pot_gibbs_3',
                  'DE_redox_1a','DG_redox_1a','pot_elec_1a','pot_gibbs_1a',
                  'DE_redox_2b','DG_redox_2b','pot_elec_2b','pot_gibbs_2b',
                  'DE_rxn_A','DG_rxn_A','pka_elec_A','pka_gibbs_A',
                  'DE_rxn_B','DG_rxn_B','pka_elec_B','pka_gibbs_B',
                  'DE_rxn_B_ref','DG_rxn_B_ref','pka_elec_B_ref','pka_gibbs_B_ref']
df_general.to_csv('results.csv', index=False, header=True, columns=labels_general) #Save file







