
# coding: utf-8

# In[4]:


import math
import numpy as np
import glob
import os
import fractions
from random import choice
from pymatgen import MPRester
from pymatgen import Structure
from pymatgen.io.vasp.sets import MPRelaxSet


# In[5]:


def run_calculations():
    #assume that calc files are written in the same directory as this script
    folders = glob.glob('./*')
    folders = [folder for folder in folders if "tag" in folder] # filters the files to retrieve relevant ones

    #program to write script
    for folder in folders:
        directory = "%s/run.slurm" % folder
        script = """#!/bin/bash -l \n
                    #SBATCH -J %s \n
                    #SBATCH -q regular \n
                    #SBATCH -N 2 \n#SBATCH -t 02:30:00 \n
                    #SBATCH -C haswell \n\nmodule load vasp \n
                    srun -n 64 ../vasp_std \nmodule load vasp \n\n
                    srun -n 64 ../vasp_std" % (folder[2:])"""
        #creates and writes run.slurm in the current file
        f= open(directory,"w+")
        f.write(script)
        f.close()
        #executes run.slurm in the current file
        os.chroot(folder)
        os.system("sbatch run.slurm")


# In[16]:


VASP_LIMIT = 160
#Connect to MP database
VAC_FRAC = 0
mpr = MPRester('sha1JueA3CU5QzB8')


# In[17]:


#Get primitive cell from MP
primitive_str = mpr.get_structure_by_material_id('mp-765279')
test_str = primitive_str.copy()
print(test_str.composition.as_dict()['Li'])


# In[18]:


#Create supercell
#Q: how to decide the supercell size?
primitive_str.make_supercell([4, 1, 1])
structure = primitive_str


# In[19]:


#Create vacancies in lithium layer. Concentrations are 0, 0.25, 0.5, 0.6, 0.75, 1
#This arrangement may be suitable for all layered-structure lithium cathodes.
#Q: we should hardcode vacancies arrangement for all kinds of materials?
"""This function hardcode the vacancy position of different concentration."""
def vacancy_position(concentration, structure):
    if concentration == 0.75:
        indice_list = np.arange(0,19,2)
        indice_list = sorted(indice_list,reverse=True)
        return indice_list
    elif concentration == 0.5:
        indice_list = np.arange(0,20)
        indice_list = sorted(indice_list,reverse=True)
        return indice_list
    elif concentration == 0.25:
        indice_list1 = [21,23,24,26,29,31,33,35,37,39]
        indice_list2 = np.arange(0,20)
        indice_list = indice_list1 + indice_list2
        indice_list = sorted(indice_list,reverse=True)
        return indice_list
    elif concentration == 0.6:
        indice_list = np.arange(24,40)
        indice_list = sorted(indice_list,reverse=True)
        return indice_list
    elif concentration == 0:
        indice = np.arange(0,40)
        indice_list = sorted(indice_list,reverse=True)
        return indice_list

"""This function remove the vacancies according to indice list, meanwhile record the coordinates.
   Return a dictionary containing all coordinates."""
def create_vacancy(vac_indice_list, structure):
    vacancy_frac_coords_list =[]
    for i in vac_indice_list:
        site = structure.pop(i)
        vacancy_frac_coords_list.append(site.frac_coords)
    return vacancy_frac_coords_list        


# In[20]:


"""This function remove the vacancies according to indice list, meanwhile record the coordinates.
   Return a dictionary containing all coordinates."""
def create_vacancy(vac_indice_list, structure):
    vacancy_frac_coords_list =[]
    for i in vac_indice_list:
        site = structure.pop(i)
        vacancy_frac_coords_list.append(site.frac_coords)
    return vacancy_frac_coords_list
"""This function swap the species at two given indice."""
def swap_element(indice1, indice2, structure):
    specie1 = structure[indice1].specie
    specie2 = structure[indice2].specie
    structure.replace(indice2,specie1)
    structure.replace(indice1,specie2)
    return structure

"""Get the distance between TM atom and all Li atoms. Returnt the indice of the three nearest ones."""
def find_nearest(indice_TM, structure):
    distance = [];
    indice_Li = structure.indices_from_symbol('Li')
    for i in indice_Li:
        d = structure.get_distance(indice_TM,i)
        distance.append(d)
    dictionary = dict(zip(indice_Li,distance))
    dictionary_indice = np.array(sorted(dictionary.items(), key = lambda item:item[1]))
    three_nearest_Li = dictionary_indice[:3, 0].astype(int)
    return three_nearest_Li


"""This function swap TM atom and Li atom.
   Then generate the input files.
   Attention: specie should be a string or character with '' !"""
def swap_Li(specie, structure):
    indice_list_specie = structure.indices_from_symbol(specie)
    print(indice_list_specie)
    indice_specie = int(choice(indice_list_specie))
    indice_Li_list = find_nearest(indice_specie, structure)
    for i in indice_Li_list:
        new_structure = structure.copy()
        new_structure = swap_element(indice_specie, i, new_structure)
        file_input = MPRelaxSet(new_structure)
        TM_string = str(indice_specie)
        Li_string = str(i)
        #insert an id tag so that glob can easily identify the proper folders, can replace 'tag' with another word
        file_input.write_input("tag_" + specie+TM_string+'Li'+Li_string)

"""This function find the three nearest vacancies around the TM atom."""

def find_nearest_vac(indice_TM, structure, vacancy_frac_coords_list):
    vacancy_distance = []
    index = np.arange(0,len(vacancy_frac_coords_list))
    site_TM = structure[indice_TM]
    xcoord = site_TM.frac_coords[0]
    ycoord = site_TM.frac_coords[1]
    zcoord = site_TM.frac_coords[2]
    for i in range(0,len(vacancy_frac_coords_list)):
        vx = vacancy_frac_coords_list[i][0]
        vy = vacancy_frac_coords_list[i][1]
        vz = vacancy_frac_coords_list[i][2]
        distance_square = math.pow((xcoord-vx),2)+math.pow((ycoord-vy),2)+math.pow((zcoord-vz),2)
        vacancy_distance.append(distance_square)
    vac_dict = dict(zip(index,vacancy_distance))
    vac_dict = np.array(sorted(vac_dict.items(), key = lambda item:item[1]))
    three_index = vac_dict[:3, 0].astype(int)
    three_nearest_vac = [vacancy_frac_coords_list[three_index[0]],vacancy_frac_coords_list[three_index[1]],vacancy_frac_coords_list[three_index[2]]]
    return three_nearest_vac

"""This function swap TM atom and vacancy position.
   Then generate input files."""
def swap_vac(specie, Li_concentration, structure):
    vac_all_indice = vacancy_position(Li_concentration, structure)
    vac_position = create_vacancy(vac_all_indice,structure)
    indice_list_specie = structure.indices_from_symbol(specie)
    indice_specie = choice(indice_list_specie)
    vac_coord_list = find_nearest_vac(indice_specie, structure, vac_position)

    for i in range(0,len(vac_coord_list)):
        new_structure = structure.copy()
        new_structure.pop(indice_specie)
        new_structure.append(specie, vac_coord_list[i])
        file_input = MPRelaxSet(new_structure)
        TM_string = str(indice_specie)
        vac_string = str(i)
        file_input.write_input(specie+TM_string+'Vac'+vac_string)

def lcm(a,b): return abs(a * b) / fractions.gcd(a,b) if a and b else 0

def get_cell_limit(species, vacancy_conc, structure):
    #get number of species
    site = len(structure.sites)
    num = test_str.composition.as_dict()[species]
    
    #no need for supercell
    if (vacancy_conc % num)== 0:
        return 1
    
    #find lcm between inverse of vacancy_conc and num
    lcm_val = lcm((1/vacancy_conc), num)
    
    #if lcm is within limit, build the supercell
    if (lcm_val / num * site <= VASP_LIMIT):
        return lcm_val / num
    
    #otherwise, change global vacancy variable
    
    num2 = vacancy_conc/num
    num2 = math.floor(num2*10)/10
    global VAC_FRAC
    VAC_FRAC = num2
    return 1
    
def make_custom_supercell(structure):
    #see how many cells can be put together
    res_struct = structure.copy()
    supercell_params = [1, 1, 1]
    cell_limit = get_cell_limit("Li", VAC_FRAC, res_struct) - 1 #number of cells that can be added to primitive cell
    struct_info = structure.as_dict().get('lattice', None)
    parameters = np.array([struct_info.get('a', None), struct_info.get('b', None), struct_info.get('c', None)])
    final_param = parameters.copy()

    while(cell_limit >= 1):
        increment = 1
        if (cell_limit >= 1):
            min_index = np.argmin(parameters)
            max_index = np.argmax(parameters)
            #adds 1 row of cells along the shortest lattice direction, until it is over the longest lattice + short lattice unit
            if((final_param[min_index] + parameters[min_index]) > (final_param[max_index] + parameters[min_index])):
                while(cell_limit >= np.product(supercell_params)):
                    #if structure is relatively cubic, just add 1 in each direction until total cubic limit is met
                    cell_limit -= np.product(supercell_params)
                    supercell_params = [x+1 for x in supercell_params]
                break

            while ((final_param[min_index] + parameters[min_index]) < (final_param[max_index] + parameters[min_index]) ):
                tmp = supercell_params.copy()
                del tmp[min_index]
                increment *= np.product(tmp)
                if (cell_limit < increment):
                    res_struct.make_supercell(supercell_params)
                    return res_struct

                final_param[min_index] = final_param[min_index] + parameters[min_index]
                supercell_params[min_index] = supercell_params[min_index] + 1
                cell_limit -= increment
        else:
            break
    res_struct.make_supercell(supercell_params)
    return res_struct

#Get primitive cell from MP
primitive_str = mpr.get_structure_by_material_id('mp-765279')
res_struct = make_custom_supercell(primitive_str)


# In[21]:


print(res_struct)


# In[87]:




