#!/usr/bin/python
# -*- coding: latin-1 -*-

import pandas as pd
import numpy as np
import os
import PCET_SMD
import scscore as scs

###---Redox potentials calculations
file_dir = 'Energy_files/DFT_B3LYP_6-31Gd_Opt_SMD'
solv_model='SMD'         #Solvation Model (PCM o SMD)
reaction_redox='1'       #Reaction to calculate according to PCET scheme
electronic_energy='Yes'  #Electronic energy? --> Yes/No
gibbs_energy='Yes'       #Gibbs free energy avail --> Yes/No

if electronic_energy == 'Yes':
    file=open('%s/rxn%s_E.txt' %(file_dir,reaction_redox),'r')  #%d for integers
    lines=file.readlines()
    elec_energies = PCET_SMD.get_energies(lines)
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
    gibbs_energies = PCET_SMD.get_energies(lines)
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
    potentials_from_E = PCET_SMD.Potentials_from_E(reaction_redox,energies_ox_s,energies_red_s)
    DE_redox_S=potentials_from_E[0]
    E_calc_elec=potentials_from_E[1]
else:
    DE_redox_S=[]
    E_calc_elec=[]
    for i in moleculas:
        DE_redox_S.append(np.nan)
        E_calc_elec.append(np.nan)
if gibbs_energy == 'Yes':
    potentials_from_G = PCET_SMD.Potentials_from_G(reaction_redox,gibbs_ox_s,gibbs_red_s)
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

###---Fitting the calculated redox potentials. Linear regression with the experimental values
if solv_model == 'PCM':
    if electronic_energy == 'Yes':
        E_fit_elec  = [np.nan for x in E_calc_elec]
    else:
        E_fit_elec  = [np.nan for x in E_calc_elec]
    
    if gibbs_energy == 'Yes':
        E_fit_gibbs = [np.nan for x in E_calc_gibbs]
    else:
        E_fit_gibbs = [np.nan for x in E_calc_gibbs]

if solv_model == 'SMD':
    if electronic_energy == 'Yes':
        E_fit_elec  = [((0.809909*x) + 0.262269) for x in E_calc_elec] #r2=0.958951
    else:
        E_fit_elec  = [np.nan for x in E_calc_elec]

    if gibbs_energy == 'Yes':
        E_fit_gibbs = [((0.804430*x) + 0.181826) for x in E_calc_gibbs] #r2=0.986857
    else:
        E_fit_gibbs = [np.nan for i in E_calc_gibbs]

#print('redox potential fitted values')
#print(E_fit_elec)
#print(E_fit_gibbs)

###---pKa calculations
reaction_pka='A'
electronic_energy='No'  #Electronic energy? --> Yes/No
gibbs_energy='No'       #Gibbs free energy avail --> Yes/No

if electronic_energy == 'Yes':
    file=open('%s/rxn%s_E.txt' %(file_dir,reaction_pka),'r')  #%d for integers
    lines=file.readlines()
    elec_energies = PCET_SMD.get_energies(lines)
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
    gibbs_energies = PCET_SMD.get_energies(lines)
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
    pka_E = PCET_SMD.pka_from_E(reaction_pka,energies_prot_s,energies_deprot_s)
    DE_Solution_ref=pka_E[0]
    DE_Solution_simple=pka_E[2]
    pka_calc_elec=pka_E[1]
else:
    DE_Solution_ref=[]
    DE_Solution_simple=[]
    pka_calc_elec=[]
    for i in molecules_pka:
        DE_Solution_ref.append(np.nan)
        DE_Solution_simple.append(np.nan)
        pka_calc_elec.append(np.nan)

if gibbs_energy == 'Yes':
    pka_G = PCET_SMD.pka_from_G(reaction_pka,gibbs_prot_s,gibbs_deprot_s)
    DG_Solution_ref=pka_G[0]
    DG_Solution_simple=pka_G[2]
    pka_calc_gibbs=pka_G[1]
else:
    DG_Solution_ref=[]
    DG_Solution_simple=[]
    pka_calc_gibbs=[]
    for i in molecules_pka:
        DG_Solution_ref.append(np.nan)
        DG_Solution_simple.append(np.nan)
        pka_calc_gibbs.append(np.nan)

#print(pka_calc_elec)
#print(pka_calc_gibbs)

###---Fitting the calculated pkas. Linear regression with the experimental values
if solv_model == 'PCM':
    if electronic_energy == 'Yes':
        pka_fit_elec = [np.nan for x in pka_calc_elec]
    else:
        pka_fit_elec = [np.nan for x in pka_calc_elec]
    
    if gibbs_energy == 'Yes':
        pka_fit_gibbs = [np.nan for x in pka_calc_egibbs]
    else:
        pka_fit_gibbs = [np.nan for x in pka_calc_egibbs]

if solv_model == 'SMD':
    if electronic_energy == 'Yes':
        pka_fit_elec = [np.nan for x in pka_calc_elec]
        pka_fit_DE_simple = [np.nan for x in DE_Solution_simple]
    else:
        pka_fit_elec = [np.nan for x in pka_calc_elec]
        pka_fit_DE_simple = [np.nan for x in DE_Solution_simple]
    
    if gibbs_energy == 'Yes':
        pka_fit_gibbs = [np.nan for x in pka_calc_gibbs]
        pka_fit_DG_simple = [np.nan for x in DG_Solution_simple]
    else:
        pka_fit_gibbs = [np.nan for x in pka_calc_gibbs]
        pka_fit_DG_simple = [np.nan for x in DG_Solution_simple]

#print('pka fitted values')
#print(pka_fit_elec)
#print(pka_fit_DE_simple)
#print(pka_fit_gibbs)
#print(pka_fit_DG_simple)

###---Calculating pka with chemaxon
os.system('mkdir Info_files')
os.system('cxcalc Coord_opt/DFT_B3LYP_6-31Gd_Opt_SMD_Diquat/smiles.smi pka -a 3 -b 3 > Info_files/chemaxon_pkas.csv')

###---Reading pkas (chemaxon and DFT) to get the ph values to calculate solubilities
pka_chemaxon = pd.read_csv('Info_files/chemaxon_pkas.csv',delimiter='\t')
pka_chemaxon['bpKa1'] = np.nan
ph_chemaxon  = pd.DataFrame(data=[7.0],columns=['ph_chemaxon'])
pka_dft = pd.DataFrame(pka_fit_gibbs)
ph_dft = pd.DataFrame(data=[7.0],columns=['ph_dft'])
smiles = pd.read_csv('Coord_opt/DFT_B3LYP_6-31Gd_Opt_SMD_Diquat/smiles.smi',header=None,sep='\t',usecols=[0],names=['smiles'])

###---Solubily calculation (logs and logd) with chemaxon using chemaxon and dft pH values
logs_params = pd.concat([smiles,ph_chemaxon,ph_dft],axis=1)
logs_params.to_csv('Info_files/chemaxon_logs_params.csv', index=False, header=False)
os.system('bash ../../Modules/logs_calc.sh')
#os.system('cxcalc Coord_opt/DFT_B3LYP_6-31Gd_Opt_SMD_Diquat/smiles.smi logs -H 7.0 > Info_files/chemaxon_logs.csv')

#Calsulating SCScore
model = scs.SCScorer()
model.restore(scs.WEIGHTS_FILE)
scscore = [model.get_score_from_smi(m)[1] for m in smiles['smiles']]
smiles['sc_score'] = scscore

###---Read de functional groups and logS files
df_func_groups = pd.read_csv('subs.csv',header=None,names=['R1','R2','R3','R4','R5','R6','R7','R8'])
df_func_groups['core'] = 'C1([H])=[N+]2CC[N+]3=C(C([H])=C([H])C([H])=C3([H]))C2=C([H])C([H])=C1([H])'
df_func_groups['type'] = 'Diquat'
logs_chemaxon  = pd.read_csv('Info_files/chemaxon_logs.csv',header=0,sep='\t',usecols=[1],names=['logs_chemaxon'])
logs_chemdft   = pd.read_csv('Info_files/chem-dft_logs.csv',header=0,sep='\t',usecols=[1],names=['logs_chemaxon-dft'])
logd_chemaxon  = pd.read_csv('Info_files/chemaxon_logd.csv',header=0,sep='\t',usecols=[1],names=['logd_chemaxon'])
logd_chemdft   = pd.read_csv('Info_files/chem-dft_logd.csv',header=0,sep='\t',usecols=[1],names=['logd_chemaxon-dft'])
#print(logs_chemaxon)

###---Generate and save the dataset
DataSet_General=list(zip(molecules,E_calc_elec,E_calc_gibbs,E_fit_elec,E_fit_gibbs,DE_redox_S,DG_redox_S,
                         molecules_pka,pka_calc_elec,pka_calc_gibbs,pka_fit_elec,pka_fit_gibbs,DE_Solution_ref,DG_Solution_ref))
values_pot_pka = ['molecule','pot_calc_elec','pot_calc_gibbs','pot_fit_elec','pot_fit_gibbs','DE_redox_sol','DG_redox_sol',
                  'molecule_pka','pka_calc_elec','pka_calc_gibbs','pka_fit_elec','pka_fit_gibbs','DE_acid_sol','DG_acid_sol']
df_pot_pka=pd.DataFrame(data=DataSet_General, columns=values_pot_pka)

df_general=pd.concat([df_func_groups,smiles,df_pot_pka,ph_chemaxon,ph_dft,logs_chemaxon,logs_chemdft,logd_chemaxon,logd_chemdft],axis=1)

values_general = ['R1','R2','R3','R4','R5','R6','R7','R8','core','type','smiles','sc_score',
                  'pot_calc_elec','pot_calc_gibbs','pot_fit_elec','pot_fit_gibbs','DE_redox_sol','DG_redox_sol',
                  'pka_calc_elec','pka_calc_gibbs','pka_fit_elec','pka_fit_gibbs','DE_acid_sol','DG_acid_sol',
                  'ph_chemaxon','ph_dft','logs_chemaxon','logs_chemaxon-dft','logd_chemaxon','logd_chemaxon-dft']
df_general.to_csv('results.csv', index=False, header=True, columns=values_general) #Save file







