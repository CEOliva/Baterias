#!/usr/bin/python
# -*- coding: latin-1 -*-

## This module contains some functions to get the energies and
## calculate redox potentials and pKa according to the PCET scheme
## It is a modifications of the original version that works only for SMD

#Function to get the energies from the rxn*.txt files
def get_energies(lines):
    molecules=[]
    energies_ox_s=[]
    energies_red_s=[]
    i=0
    while i < (len(lines)/2):              #Se divide entre 4 porque hay 2 energías --> se forman 2 bloques
        molecules.append((i+1))
        start=lines[0].find('=')
        end=lines[0].find('.')            #Encuentra la posición del símbolo'=' y lee los dígitos posteriores
        if end <= 0:
        	end=start+5
        energies_ox_s.append(float(lines[i][start+1:end+11]))
        energies_red_s.append(float(lines[i+int(len(lines)/2)][start+1:end+11]))
        i=i+1
    return [molecules, energies_ox_s, energies_red_s]

#Options to calculate redox potentials and pka's according to the PCET scheme
electron_transfer=['1','2','3','4','5','6']
two_electron_transfer=['12','34','56']
proton_transfer1=['A','C','E']
proton_transfer2=['B','D','F']
proton_electron_transfer=['1A','2C','3B','4D']
two_proton_electron_transfer=['1A4D']

#Función para calcular potenciales redox a partir de energías electrónicas
def Potentials_from_E(reaction,energies_ox_s,energies_red_s):
    F=96.4853                      #Faraday constant (kJ V-1 mol-1)
    SHE=4.44                       #redox potential (V) for SHE in water
    H2_GP=-1.17548236813           #Electronic energy (Ha) of H2 optimized in gas phase at B3LYP/6-31G(d) level of theory
    H2_S_PCM=-1.17556856786        #Electronic energy (Ha) of H2 optimized in water at B3LYP/6-31G(d)/PCM level of theory
    H2_S_SMD=-1.17268376198        #Electronic energy (Ha) of H2 optimized in water at B3LYP/6-31G(d)/SMD level of theory
    DE_redox_S=[]
    E_calc=[]
    if reaction in electron_transfer:
        n=1.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using electronic energies for one electron transfer')
        while i < len(energies_ox_s):
            DE_redox_S.append((energies_red_s[i]-energies_ox_s[i])*2625.49962) #Usando solo la diferencia de las especies en disolución
            E_calc.append((DE_redox_S[i]/(-n*F))-SHE)
            i=i+1
    elif reaction in proton_electron_transfer:
        n=1.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using electronic energies for one proton-electron transfer')
        while i < len(energies_ox_s):
            DE_redox_S.append((energies_red_s[i]-energies_ox_s[i]-(H2_S_SMD/2))*2625.49962)
            E_calc.append((DE_redox_S[i]/(-n*F)))
            i=i+1
    elif reaction in two_electron_transfer:
        n=2.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using electronic energies for two electron transfer')
        while i < len(energies_ox_s):
            DE_redox_S.append((energies_red_s[i]-energies_ox_s[i])*2625.49962)
            E_calc.append((DE_redox_S[i]/(-n*F))-SHE)
            i=i+1
    elif reaction in two_proton_electron_transfer:
        n=2.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using electronic energies for two proton-electron transfer')
        while i < len(energies_ox_s):
            DE_redox_S.append((energies_red_s[i]-energies_ox_s[i]-H2_sol_SMD)*2625.49962)
            E_calc.append((DE_redox_S[i]/(-n*F)))
            i=i+1
    return [DE_redox_S, E_calc]

def Potentials_from_G(reaction,gibbs_ox_s,gibbs_red_s):
    F=96.4853                      #Faraday constant (kJ V-1 mol-1)
    SHE=4.44                       #redox potential (V) for SHE in water
    H2_GP_G=-1.176819              #Electronic energy (Ha) of H2 optimized in gas phase at B3LYP/6-31G(d) level of theory
    H2_S_PCM_G=-1.176923           #Electronic energy (Ha) of H2 optimized in water at B3LYP/6-31G(d)/PCM level of theory
    H2_S_SMD_G=-1.174035           #Electronic energy (Ha) of H2 optimized in water at B3LYP/6-31G(d)/SMD level of theory
    DG_redox_S=[]
    E_calc_G=[]
    if reaction in electron_transfer:
        n=1.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using Gibbs energies for one electron transfer')
        while i < len(gibbs_ox_s):
            DG_redox_S.append((gibbs_red_s[i]-gibbs_ox_s[i])*2625.49962) #Usando solo la diferencia de las especies en disolución
            E_calc_G.append((DG_redox_S[i]/(-n*F))-SHE)
            i=i+1
    elif reaction in proton_electron_transfer:
        n=1.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using Gibbs energies for one proton-electron transfer')
        while i < len(gibbs_ox_s):
            DG_redox_S.append((gibbs_red_s[i]-gibbs_ox_s[i]-(H2_S_SMD_/2))*2625.49962)
            E_calc_G.append((DG_redox_S[i]/(-n*F)))
            i=i+1
    elif reaction in two_electron_transfer:
        n=2.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using Gibbs energies for two electron transfer')
        while i < len(gibbs_ox_s):
            DG_redox_S.append((gibbs_red_s[i]-gibbs_ox_s[i])*2625.49962)
            E_calc_G.append((DG_redox_S[i]/(-n*F))-SHE)
            i=i+1
    elif reaction in two_proton_electron_transfer:
        n=2.0     #number of e- exchanged
        i=0
        print('Redox potentials calculated using Gibbs energies for two proton-electron transfer')
        while i < len(gibbs_ox_s):
            DG_redox_S.append((gibbs_red_s[i]-gibbs_ox_s[i]-H2_sol_SMD_G)*2625.49962)
            E_calc_G.append((DG_redox_S[i]/(-n*F)))
            i=i+1
    return [DG_redox_S, E_calc_G]

def pka_from_E(reaction,energies_prot_s,energies_deprot_s):
    E_H2Ref_gas = -495.998195          #Electronic energy (Ha) of 2,2'-H2Bpy2+ (H2Ref2+) optimized in gas phase
    E_HRef_gas = -495.777608           #Electronic energy (Ha) of 2,2'-HBpy+ (HRef+) optimized in gas phase
    E_Ref_gas = -495.376530            #Electronic energy (Ha) of 2,2'-Bpy- (Ref) optimized in gas phase
    E_H2Ref_solution_PCM = -496.266693 #Electronic energy (Ha) of 2,2'-H2Bpy2+ (H2Ref2+) optimized in water at B3LYP/6-31G(d)/PCM level of theory
    E_HRef_solution_PCM = -495.847175  #Electronic energy (Ha) of 2,2'-HBpy+ (HRef+) optimized in water at B3LYP/6-31G(d)/PCM level of theory
    E_Ref_solution_PCM = -495.387254   #Electronic energy (Ha) of 2,2'-Bpy- (Ref) optimized in water at B3LYP/6-31G(d)/PCM level of theory
    E_H2Ref_solution_SMD = -496.294652 #Electronic energy (Ha) of 2,2'-H2Bpy2+ (H2Ref2+) optimized in water at B3LYP/6-31G(d)/SMD level of theory
    E_HRef_solution_SMD = -495.856530  #Electronic energy (Ha) of 2,2'-HBpy+ (HRef+) optimized in water at B3LYP/6-31G(d)/SMD level of theory
    E_Ref_solution_SMD = -495.391348   #Electronic energy (Ha) of 2,2'-Bpy- (Ref) optimized in water at B3LYP/6-31G(d)/SMD level of theory
    R =0.008314472                     #Constant of gases (kJ mol-1 k-1)
    T = 298.15                         #Temperature (K)
    pka1_HRef = 4.3                    #pKa1 experimental of 2,2'-HBpy+ --> 2,2'-Bpy- + H+
    pka2_H2Ref = -0.2                  #pKa2 experimental of 2,2'-H2Bpy2+ --> 2,2'-HBpy+ + H+
    DE_solution=[]                     #Energy diference of the reaction including the reference system
    pka_calc_elec=[]
    DE_simple=[]                       #Energy diference of the reaction without reference system
    if reaction in proton_transfer1:
        i=0
        print('pka1 calculated using electronic energies for one proton transfer')
        while i < len(energies_prot_s):
            DE_simple.append((energies_deprot_s[i] - energies_prot_s[i])*627.50961)
            DE_solution.append((energies_deprot_s[i] + E_HRef_solution_SMD - energies_prot_s[i] - E_Ref_solution_SMD)*2625.49962)
            pka_calc_elec.append((DE_solution[i]/(R*T*2.303))+pka1_HRef)
            i=i+1
    elif reaction in proton_transfer2:
        i=0
        print('pka2 calculated using electronic energies for one proton transfer')
        while i < len(energies_prot_s):
            DE_simple.append((energies_deprot_s[i] - energies_prot_s[i])*627.50961)
            DE_solution.append((energies_deprot_s[i] + E_H2Ref_solution_SMD - energies_prot_s[i] - E_HRef_solution_SMD)*2625.49962)
            pka_calc_elec.append((DE_solution[i]/(R*T*2.303))+pka2_H2Ref)
            i=i+1
    return [DE_solution, pka_calc_elec, DE_simple]

def pka_from_G(reaction,gibbs_prot_s,gibbs_deprot_s):
    G_H2Ref_gas = -495.847911          #Gibbs free energy (Ha) of 2,2'-H2Bpy2+ (H2Ref2+) optimized in gas phase
    G_HRef_gas = -495.639542           #Gibbs free energy (Ha) of 2,2'-HBpy+ (HRef+) optimized in gas phase
    G_Ref_gas = -495.253176            #Gibbs free energy (Ha) of 2,2'-Bpy- (Ref) optimized in gas phase
    G_H2Ref_solution_PCM = -496.115371 #Gibbs free energy (Ha) of 2,2'-H2Bpy2+ (H2Ref2+) optimized in water at B3LYP/6-31G(d)/PCM level of theory
    G_HRef_solution_PCM = -495.709946  #Gibbs free energy (Ha) of 2,2'-HBpy+ (HRef+) optimized in water at B3LYP/6-31G(d)/PCM level of theory
    G_Ref_solution_PCM = -495.263602   #Gibbs free energy (Ha) of 2,2'-Bpy- (Ref) optimized in water at B3LYP/6-31G(d)/PCM level of theory
    G_H2Ref_solution_SMD = -496.141776 #Gibbs free energy (Ha) of 2,2'-H2Bpy2+ (H2Ref2+) optimized in water at B3LYP/6-31G(d)/SMD level of theory
    G_HRef_solution_SMD = -495.718138  #Gibbs free energy (Ha) of 2,2'-HBpy+ (HRef+) optimized in water at B3LYP/6-31G(d)/SMD level of theory
    G_Ref_solution_SMD = -495.267150   #Gibbs free energy (Ha) of 2,2'-Bpy- (Ref) optimized in water at B3LYP/6-31G(d)/SMD level of theory
    R =0.008314472                     #Constant of gases (kJ mol-1 k-1)
    T = 298.15                         #Temperature (K)
    pka1_HRef = 4.3                    #pKa1 experimental of 2,2'-HBpy+ --> 2,2'-Bpy- + H+
    pka2_H2Ref = -0.2                  #pKa2 experimental of 2,2'-H2Bpy2+ --> 2,2'-HBpy+ + H+
    DG_solution=[]
    pka_calc_gibbs=[]
    DG_simple=[]
    if reaction in proton_transfer1:
        i=0
        print('pka1 calculated using Gibbs energies for one proton transfer')
        while i < len(gibbs_prot_s):
            DG_simple.append((gibbs_deprot_s[i] - gibbs_prot_s[i])*627.50961)
            DG_solution.append((gibbs_deprot_s[i] + G_HRef_solution_SMD - gibbs_prot_s[i] - G_Ref_solution_SMD)*2625.49962)
            pka_calc_gibbs.append((DG_solution[i]/(R*T*2.303))+pka1_HRef)
            i=i+1
    elif reaction in proton_transfer2:
        i=0
        print('pka2 calculated using Gibbs energies for one proton transfer')
        while i < len(gibbs_prot_s):
            DG_simple.append((gibbs_deprot_s[i] - gibbs_prot_s[i])*627.50961)
            DG_solution.append((gibbs_deprot_s[i] + G_H2Ref_solution_SMD - gibbs_prot_s[i] - G_HRef_solution_SMD)*2625.49962)
            pka_calc_gibbs.append((DG_solution[i]/(R*T*2.303))+pka2_H2Ref)
            i=i+1
    return [DG_solution, pka_calc_gibbs, DG_simple]
