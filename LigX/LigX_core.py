import numpy as np
import math
import warnings
import csv
import os
import shutil
import sys
import time
import statistics
import subprocess
from copy import deepcopy
import logging
from logging import getLogger, StreamHandler, FileHandler, Formatter
import copy
import re
import itertools
from importlib import resources
from ._version import __version__
import platform

class LigXRun:
    """Core class in LigX<br>

    ------------------
    ## Key Parameters
    **linknum** (*int*):<br>
        Maximum number of fragment linkings. (default: 5)<br>
    **outnum** (*int*):<br>
        Number of structures to generate in this run. (default: 2000)<br>
    **ex_frag** (*str*):<br>
        Relative path to the external fragment set specified by the user.
        <br>(example: 'external_fragment/fragment')<br>
    **ex_rollfrag** (*str*):
        Relative path to the external rollout fragment set specified by the user.
        <br>(example: 'external_fragment/rollout_fragment')<br>
    ------------------

    ## Key Attributes

    **smiles_filter** (*bool*): <br>
    Set to `False` when explicitly turning off the SMILES filter.
    LigX program automatically sets the value to `False`
    in environments where RDkit is unavailable.(default: `True`)<br>
    **result_dir** (*str*):<br>
    Directory name where output is saved. (default: 'result')<br>
    **ligand_xyz** (*str*):<br>
    Directory name where the ligand (initiator) structure used for
    input is stored. (default: 'ligand_xyz')<br>
    **dir_protein_xyz** (*str*):<br>
    Directory name where the protein structure used for input
    is stored. (default: 'protein_xyz')<br>
    **cutoff_score_init** (*float*):<br>
    Initial cutoff value for ligand selection. (default: 15)<br>
    **cutoff_score_fin** (*float*):<br>
    The maximum cutoff value that increases with the number of
    generated structures. Increases linearly from cutoff_score_init.
    (default: 15)<br>
    **sigma_sum_cut** (*float*):<br>
    Collision detection value for stetic clash
    between fragments and initiators. Reducing this value increases
    structural flexibility, but also makes it easier to generate
    unreasonable structures. (default: 10)<br>
    **rollout_num** (*int*):<br>
    The number of fragments for which rollout is performed after evaluating the interaction of fragments. (default: 5) <br>
    **clean** (*bool*):<br>
    For MacOS only. LigX detects the operating system. If running on macOS, clean will set to `True` and the program will execute clean.sh to remove automatically generated files such as .DS_Store and AppleDouble. <br>
    Set this option to `False` if you do not allow the execution of clean.sh. <br>There is no need to configure this setting on Windows or Linux systems. (default: `False`) <br>

    ------------------
    """
    def __init__(self, linknum = None, outnum = None, ex_frag = None, ex_rollfrag = None, logfile_root = None, result_dir_root = None, paramlog_root = None, ligand_xyz_root = None, dir_protein_xyz_root = None, desired_frag_root = None, desired_frag_rollout_root = None):

        if logfile_root is None:
            self.logfile = 'logfile.log'
        else:
            self.logfile = logfile_root
        if paramlog_root is None:
            self.paramlog = 'Param.log'
        else:
            self.paramlog = paramlog_root
        if result_dir_root is None:
            self.result_dir = 'result'
        else:
            self.result_dir = result_dir_root
        if ligand_xyz_root is None:
            self.ligand_xyz = 'ligand_xyz'
        else:
            self.ligand_xyz = ligand_xyz_root
        if dir_protein_xyz_root is None:
            self.dir_protein_xyz = 'protein_xyz'
        else:
            self.dir_protein_xyz = dir_protein_xyz_root
        if ex_frag is None:
            self.fragment_set = None
        else:
            self.fragment_set = ex_frag
        if ex_rollfrag is None:
            self.rollout_set = None
        else:
            self.rollout_set = ex_rollfrag
        if desired_frag_root is None:
            self.desired_frag = None
        else:
            self.desired_frag = desired_frag_root
        if desired_frag_rollout_root is None:
            self.desired_frag_rollout = None
        else:
            self.desired_frag_rollout = desired_frag_rollout_root

        if linknum is None:
            self.cycle = 5
        else:
            self.cycle = linknum

        if outnum is None:
            self.max_conf = 2000
        else:
            self.max_conf = outnum

        self.des_frag_cycle = -1
        self.max_branch_core = 0
        self.max_branch_frag = 0
        self.result_conf = self.max_conf
        self.evaluation_param = 'LE_interact_E'
        self.evaluation_method = 'highest'
        self.cutoff_score_init = 15
        self.cutoff_score_fin = 15
        self.rollout_num = 5
        self.rollout_frag_max = 3
        self.image_dir_name = 'image_ligx'
        self.count_for_data_processing = -1
        self.rollout_conf_inner = 1
        self.rollout_conf_inner_des = 20
        self.sigma_sum_cut = 10
        self.clean = False
        if platform.system() == 'Darwin':
            self.clean = True
        else:
            self.clean = False



    def exec_ligx(self):

        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
        except Exception as e:
            self.smiles_filter = False
            red_print("SMILES filter disabled (RDKit import failed): " +str(e))
        else:
            self.smiles_filter = True
            from rdkit.Chem import AllChem
            from rdkit.Chem import Draw
            from rdkit import RDLogger
            RDLogger.DisableLog('rdApp.*')
            red_print("SMILES filter enabled")

        logfile = os.path.join(os.getcwd(),self.logfile)
        paramlog = os.path.join(os.getcwd(),self.paramlog)
        result_dir = os.path.join(os.getcwd(),self.result_dir)

        dir_parent = os.path.join(os.getcwd(),self.ligand_xyz)
        dir_protein_xyz = os.path.join(os.getcwd(),self.dir_protein_xyz)

        if self.fragment_set == None:
            fragment_set = str(resources.files(__package__ + ".fragment"))
        else:
            fragment_set = os.path.join(os.getcwd(),self.fragment_set)

        if self.rollout_set == None:
            rollout_set = str(resources.files(__package__ + ".rollout_fragment"))
        else:
            rollout_set = os.path.join(os.getcwd(),self.rollout_set)

        if self.desired_frag == None:
            desired_frag_root = str(resources.files(__package__ + ""))
        else:
            desired_frag_root = os.path.join(os.getcwd(),self.desired_frag)

        if self.desired_frag_rollout == None:
            desired_frag_rollout_root = str(resources.files(__package__ + ""))
        else:
            desired_frag_rollout_root = os.path.join(os.getcwd(),self.desired_frag_rollout)

        cycle = self.cycle
        #desired_fragを入れるcycle（0にしてはだめ,-1だとdesired_fragmentは入らない）
        des_frag_cycle = self.des_frag_cycle

        if des_frag_cycle >= 1:
            desired_frag = [desired_frag_root]
            desired_frag_rollout = [desired_frag_rollout_root]
        else:
            desired_frag = ['']
            desired_frag_rollout = ['']

        #clean_path = os.path.join(os.getcwd(), "clean.sh")
        clean_path = os.path.join(resources.files(__package__),'clean.sh')

        if self.clean == True:
            subprocess.call([clean_path])
        max_branch_core = self.max_branch_core

        max_branch_frag = self.max_branch_frag
        result_conf = self.result_conf

        max_conf = self.max_conf

        cutoff_score_init = self.cutoff_score_init
        cutoff_score_fin = self.cutoff_score_fin

        rollout_num = self.rollout_num

        rollout_frag_max = self.rollout_frag_max

        evaluation_param = self.evaluation_param

        evaluation_method = self.evaluation_method

        count_for_data_processing = self.count_for_data_processing

        image_dir_name = self.image_dir_name
        rollout_conf_inner = self.rollout_conf_inner
        rollout_conf_inner_des = self.rollout_conf_inner_des

        sigma_sum_cut = self.sigma_sum_cut

        if os.path.isdir(result_dir):
            shutil.rmtree(result_dir)

        os.mkdir(result_dir)
        image_dir = result_dir +'/' +image_dir_name
        os.mkdir(image_dir)


        if(os.path.isfile(logfile)):
            os.remove(logfile)

        if(os.path.isfile(paramlog)):
            os.remove(paramlog)
        logger=my_logfile(logfile,'INFO','a')
        logger_notime=my_logfile_notime(logfile,'INFO_notime','a')
        logger_state=my_logfile(paramlog,'Param','a')
        logger_state_notime=my_logfile_notime(paramlog,'Param_notime','a')


        logger.info('Execute LigX program')
        logger_notime.info('\nLigX: Ligand Extending Program for Medicinal Chemistry. (2025)')
        logger_notime.info('\nVersion: ' + str(__version__))

        logger_state.info('Execute LigX program')
        logger_state_notime.info('\nLigX: Ligand Extending Program for Medicinal Chemistry. (2025)')
        logger_state_notime.info('\nVersion: ' + str(__version__))
        logger_state_notime.info('\nDetailed parameters will be written in this file.')

        logger_state_notime.info('\n################################LigXLigX###############################')

        desired_list = []

        if desired_frag != ['']:
            for des in desired_frag:
                if self.clean == True:
                    subprocess.call([clean_path])
                despath_list = os.listdir(des)
                for i in despath_list:
                    if i != 'charge.txt' and i != 'param.txt':
                        desired_list.append(des+'/'+i)

                    des_category = extract_des_category(des)

        if desired_frag == ['']:
            des_category = None

        desired_list_rollout = []
        if desired_frag_rollout != ['']:
            for des in desired_frag_rollout:
                if self.clean == True:
                    subprocess.call([clean_path])
                despath_list_rollout = os.listdir(des)
                for i in despath_list_rollout:
                    if i != 'charge.txt' and i != 'param.txt':
                        desired_list_rollout.append(des+'/'+i)
                    rol_des_category = extract_des_category(des)

        if desired_frag_rollout == ['']:
            rol_des_category = None

        if self.clean == True:
            subprocess.call([clean_path])
        core_dir_list_root = os.listdir(dir_parent)

        core_dir_list = [os.path.join(dir_parent,x) for x in core_dir_list_root]

        protein_xyz_dir_list = os.listdir(dir_protein_xyz)
        protein_xyz_dir_list.sort()

        #Load protein structue
        protein_index = 0
        all_protein_cord = []
        for protein_xyz_dir in protein_xyz_dir_list:
            var_name = 'protein_cord_'+str(protein_index)
            protein_path = os.path.join(dir_protein_xyz,protein_xyz_dir)
            protein_cord,empty_name = read_xyz_data(protein_path)
            exec("{} = protein_cord".format(var_name))
            all_protein_cord.append(protein_cord)
            protein_index +=1

        num_protein_structures = len(all_protein_cord)

        count_generated = 0
        try_num = 0

        all_result = []

        all_score_strategy_0 = ['ID']
        time_init = 0

        for i in range(num_protein_structures):
            x = 'LJ['+str(i)+']'
            all_score_strategy_0.append(x)
        all_score_strategy_0.append('LJ_average')

        for i in range(num_protein_structures):
            x = 'coulomb['+str(i)+']'
            all_score_strategy_0.append(x)
        all_score_strategy_0.append('coulomb_average')

        for i in range(num_protein_structures):
            x = 'interact_E['+str(i)+']'
            all_score_strategy_0.append(x)
        all_score_strategy_0.append('interact_E_average')

        for i in range(num_protein_structures):
            x = 'LE_LJ['+str(i)+']'
            all_score_strategy_0.append(x)
        all_score_strategy_0.append('LE_LJ_average')

        for i in range(num_protein_structures):
            x = 'LE_interact_E['+str(i)+']'
            all_score_strategy_0.append(x)
        all_score_strategy_0.append('LE_interact_E_average')

        for i in range(num_protein_structures):
            x = 'LE_LJ+coulomb['+str(i)+']'
            all_score_strategy_0.append(x)
        all_score_strategy_0.append('LE_LJ+coulomb_average')
        all_score_strategy_0.append('strategy')
        all_score_strategy_0.append('SMILES')
        all_score_strategy =[all_score_strategy_0]

        transform_dict = get_transform_dict_hyd(fragment_set)

        setting_info = sum_set_info(cycle,des_frag_cycle,max_branch_core,max_branch_frag,result_conf,max_conf,cutoff_score_init,cutoff_score_fin,rollout_num,rollout_frag_max,fragment_set,rollout_set,desired_frag,desired_frag_rollout,evaluation_param,evaluation_method,core_dir_list,protein_xyz_dir_list)

        logger_notime.info('\n################################ Setting Info ###############################\n')
        logger_notime.info(setting_info)
        logger_notime.info('############################# Setting Info Ending ###########################\n')


        blue_print('Loading fragments…')
        logger_notime.info('Loading fragments…')
        logger_notime.info('\n################################LigXLigX###############################')

        #Initialization of Fragments

        Frag_F1L = Fragments('F1L',transform_dict)
        Frag_F2L = Fragments('F2L',transform_dict)
        Frag_C1L = Fragments('C1L',transform_dict)
        Frag_C2L = Fragments('C2L',transform_dict)
        Frag_dict = {'F1L': Frag_F1L, 'F2L': Frag_F2L,'C1L': Frag_C1L, 'C2L': Frag_C2L}

        if 'F3L' in os.listdir(fragment_set):
            Frag_F3L = Fragments('F3L',transform_dict)
            Frag_dict['F3L'] = Frag_F3L
        if 'F2LB' in os.listdir(fragment_set):
            Frag_F2LB = Fragments('F2LB',transform_dict)
            Frag_dict['F2LB'] = Frag_F2LB
        if 'C3L' in os.listdir(fragment_set):
            Frag_C3L = Fragments('C3L',transform_dict)
            Frag_dict['C3L'] = Frag_C3L
        if 'XHD' in os.listdir(fragment_set):
            Frag_XHD = Fragments('XHD',transform_dict)
            Frag_dict['XHD'] = Frag_XHD
        if 'XFO' in os.listdir(fragment_set):
            Frag_XFO = Fragments('XFO',transform_dict)
            Frag_dict['XFO'] = Frag_XFO

        Frag_desired = desFragments(desired_list,des_category)
        #Frag_dict = {'F1L': Frag_F1L, 'F2L': Frag_F2L, 'F3L': Frag_F3L, 'F2LB': Frag_F2LB, 'C1LB': Frag_C1LB, 'C1L': Frag_C1L, 'C2L': Frag_C2L, 'C3L': Frag_C3L,'XHD': Frag_XHD, 'XFO': Frag_XFO}
        start = time.time()
        logger_notime.info('Generating molecules…')
        logger_notime.info('\n################################LigXLigX###############################')

        count_for_top = 0
        all_state = [['init',100,1]]

        #Generation
        while count_generated < max_conf:
            cutoff_score = cutoff_score_init + (count_generated*(cutoff_score_fin-cutoff_score_init)/(max_conf))

            current_branch_core = 0
            current_branch_frag = 0
            current_cycle = 1
            depth = 0

            EOF_value = False
            blue_print('Generating molecules…')
            while EOF_value == False:

                if depth == 0:
                    if self.clean == True:
                        subprocess.call([clean_path])
                    each_strategy = []
                    infomation = 'Trying Sequence No.' +str(try_num+1)
                    logger.info(infomation)

                    for dir_core in core_dir_list:

                        current_state = []
                        initial_remaining_data,S,Q,P,total_charge,transform_log_main,transform_log_side,c_for_next,param_for_next,added_frag,ref_frag_list_next = import_core(dir_core)

                        current_state.append(dir_core)

                        state_for_rollout = current_state

                        first_elements = [sublist[0] for sublist in all_state]
                        add_cycle = 0
                        core = Core_mol(initial_remaining_data,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next)
                        Q_atom = core.Q_atom
                        if not any(current_state == item for item in first_elements):

                            if current_cycle != des_frag_cycle:

                                all_atom = [x for i in core.operation for x in Frag_dict[i].atom]
                                all_coord = [x for i in core.operation for x in Frag_dict[i].coord]
                                all_charge = [x for i in core.operation for x in Frag_dict[i].charge]
                                all_each_frag = [x for i in core.operation for x in Frag_dict[i].each_frag]
                                all_len_list = [x for i in core.operation for x in Frag_dict[i].len_list]
                                all_central_atom = [x for i in core.operation for x in Frag_dict[i].central_atom]
                                all_category = [x for i in core.operation for x in Frag_dict[i].category]
                                all_param = [x for i in core.operation for x in Frag_dict[i].param]
                                all_frag = [x for i in core.operation for x in Frag_dict[i].frag_name]
                                all_add_branch_core = [x for i in core.operation for x in Frag_dict[i].add_branch_core]
                                all_add_branch_frag = [x for i in core.operation for x in Frag_dict[i].add_branch_frag]
                                all_ref_frag_list = [x for i in core.operation for x in Frag_dict[i].ref_frag_list]
                                split_ind_list_1,split_ind_list_2 = create_split_ind_list(all_len_list)
                                possibility = True

                            elif current_cycle == des_frag_cycle and any(x == Frag_desired.frag_category_val[0] for x in core.operation):
                                all_atom = [x for x in Frag_desired.atom]
                                all_coord = [x for x in Frag_desired.coord]
                                all_charge = [x for x in Frag_desired.charge]
                                all_each_frag = [x for x in Frag_desired.each_frag]
                                all_len_list = [x for x in Frag_desired.len_list]
                                all_central_atom = [x for x in Frag_desired.central_atom]
                                all_category = [x for x in Frag_desired.category]
                                all_param = [x for x in Frag_desired.param]
                                all_frag = [x for x in Frag_desired.frag_name]
                                all_add_branch_core = [x for x in Frag_desired.add_branch_core]
                                all_add_branch_frag = [x for x in Frag_desired.add_branch_frag]
                                all_ref_frag_list = [x for x in Frag_desired.ref_frag_list]
                                split_ind_list_1,split_ind_list_2 = create_split_ind_list(all_len_list)
                                possibility = True

                            elif current_cycle == des_frag_cycle and not any(x == Frag_desired.frag_category_val[0] for x in core.operation):

                                escape_state = deepcopy(to_frons_state)
                                first_list = [sublist[0] for sublist in all_state]
                                escape_ind= first_list.index(escape_state[0])
                                all_state[escape_ind][1] = 100
                                all_state[escape_ind][2] = 10000000
                                possibility = False

                            if possibility:
                                score = execute_rollout_v2(initial_remaining_data,S,Q,P,total_charge,transform_log_main,transform_log_side,c_for_next,param_for_next,added_frag,ref_frag_list_next,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,evaluation_param,evaluation_method,desired_list_rollout,rollout_set,all_protein_cord,add_cycle,state_for_rollout,rol_des_category,des_frag_cycle,dir_core,cutoff_score_init,rollout_conf_inner,rollout_conf_inner_des,sigma_sum_cut)

                                all_state.append([current_state,score,1])

                    core_state = [item for item in all_state if len(item[0]) == 1]
                    core_state = sorted(core_state, key=lambda x: (x[2],x[1]))
                    selected_core_state = []
                    for i in core_state:
                        if i[1] <= cutoff_score and i[2]<10000000:
                            selected_core_state.append(i)

                    if selected_core_state ==[]:
                        if count_generated == 0:

                            print('no possibility in this core')
                            logger.info('no possibility in this core')
                            logger_notime.info('\n################################LigXLigX###############################')
                            sys.exit()
                        else:
                            print('the search is terminated because of no possibility.')
                            logger.info('no possibility in this core')
                            logger_notime.info('\n################################LigXLigX###############################')
                            max_conf = count_generated
                            break

                    else:
                        ubc_ind = calc_UCB(selected_core_state,cutoff_score)
                        core = selected_core_state[ubc_ind][0][0]

                        initial_remaining_data,S,Q,P,total_charge,transform_log_main,transform_log_side,c_for_next,param_for_next,added_frag,ref_frag_list_next = import_core(dir_core)
                        core = Core_mol(initial_remaining_data,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next)


                        first_list = [sublist[0] for sublist in all_state]

                        ind= first_list.index(selected_core_state[ubc_ind][0])

                        all_state[ind][2] +=1

                        depth = 1
                        next_state = selected_core_state[ubc_ind]


                #depth != 0
                else:
                    temp_current_state = []




                    to_frons_state = next_state

                    temp_state = next_state

                    current_state = []

                    if current_cycle != des_frag_cycle:
                        #その時点で連結可能なすべてのフラグメントを処理
                        all_atom = [x for i in core.operation for x in Frag_dict[i].atom]
                        all_coord = [x for i in core.operation for x in Frag_dict[i].coord]
                        all_charge = [x for i in core.operation for x in Frag_dict[i].charge]
                        all_each_frag = [x for i in core.operation for x in Frag_dict[i].each_frag]
                        all_len_list = [x for i in core.operation for x in Frag_dict[i].len_list]
                        all_central_atom = [x for i in core.operation for x in Frag_dict[i].central_atom]
                        all_category = [x for i in core.operation for x in Frag_dict[i].category]
                        all_param = [x for i in core.operation for x in Frag_dict[i].param]
                        all_frag = [x for i in core.operation for x in Frag_dict[i].frag_name]
                        all_add_branch_core = [x for i in core.operation for x in Frag_dict[i].add_branch_core]
                        all_add_branch_frag = [x for i in core.operation for x in Frag_dict[i].add_branch_frag]
                        all_ref_frag_list = [x for i in core.operation for x in Frag_dict[i].ref_frag_list]
                        all_frag_category_val = [x for i in core.operation for x in Frag_dict[i].frag_category_val]
                        split_ind_list_1,split_ind_list_2 = create_split_ind_list(all_len_list)
                        frag_c_for_next = [x for i in core.operation for x in Frag_dict[i].c_for_next]
                        possibility = True
                        desired_check = False


                    elif current_cycle == des_frag_cycle and any(x == Frag_desired.frag_category_val[0] for x in core.operation):

                        all_atom = [x for x in Frag_desired.atom]
                        all_coord = [x for x in Frag_desired.coord]
                        all_charge = [x for x in Frag_desired.charge]
                        all_each_frag = [x for x in Frag_desired.each_frag]
                        all_len_list = [x for x in Frag_desired.len_list]
                        all_central_atom = [x for x in Frag_desired.central_atom]
                        all_category = [x for x in Frag_desired.category]
                        all_param = [x for x in Frag_desired.param]
                        all_frag = [x for x in Frag_desired.frag_name]
                        all_add_branch_core = [x for x in Frag_desired.add_branch_core]
                        all_add_branch_frag = [x for x in Frag_desired.add_branch_frag]
                        all_ref_frag_list = [x for x in Frag_desired.ref_frag_list]
                        all_frag_category_val = [x for x in Frag_desired.frag_category_val]
                        all_frag_category_val = [x for x in Frag_desired.frag_category_val]
                        frag_c_for_next = [x for x in Frag_desired.c_for_next]
                        split_ind_list_1,split_ind_list_2 = create_split_ind_list(all_len_list)
                        possibility = True
                        desired_check = True


                    elif current_cycle == des_frag_cycle and not any(x == Frag_desired.frag_category_val[0] for x in core.operation):
                        escape_state = deepcopy(to_frons_state)
                        first_list = [sublist[0] for sublist in all_state]
                        escape_ind= first_list.index(escape_state[0])
                        all_state[escape_ind][1] = 100
                        all_state[escape_ind][2] = 10000000
                        possibility = False
                        depth = 0

                    if possibility:

                        temp_state_0 = temp_state[0]
                        first_elements_for_check = [sublist[0] for sublist in all_state if len(sublist[0])==(len(temp_state_0)+1) and sublist[0][:len(temp_state_0)]== temp_state_0]


                        int_result = link_result(core,all_atom,all_coord,all_charge,all_central_atom,split_ind_list_1,split_ind_list_2,all_category,all_param,all_len_list,all_frag,all_protein_cord,evaluation_param,evaluation_method,all_ref_frag_list,frag_c_for_next,sigma_sum_cut)
                        if first_elements_for_check == []:

                            int_result_list = int_result.evaluation.tolist()
                            current_state_0 = [[a,b,2] for a,b in zip(all_each_frag,int_result_list)]
                            current_state = deepcopy(current_state_0)
                            sorted_current_state = sorted(current_state, key=lambda x: (x[2],x[1]))
                            sorted_current_state_list = [sublist[1] for sublist in sorted_current_state if sublist[1] <= 0]
                            sorted_current_state_frag_dir = [sublist[0][0].split('/')[-2] for sublist in sorted_current_state if sublist[1] <= 0]
                            if desired_check == False:
                                if len(sorted_current_state_list) >= rollout_num:
                                    new_sorted_current_state_list = []
                                    frag_count_list = []
                                    x = 0
                                    for s,f in zip(sorted_current_state_list,sorted_current_state_frag_dir):
                                        if not any(f == a[0] for a in frag_count_list):
                                            frag_count_list.append([f,1])
                                            new_sorted_current_state_list.append(s)
                                            x += 1
                                            if x >= rollout_num:
                                                break

                                        else:
                                            index_frag = [a[0] for a in frag_count_list].index(f)
                                            if frag_count_list[index_frag][1] <= rollout_frag_max:
                                                frag_count_list[index_frag][1] += 1
                                                new_sorted_current_state_list.append(s)
                                                x += 1
                                                if x >= rollout_num:
                                                    break
                                            else:
                                                frag_count_list[index_frag][1] += 1

                                    sorted_current_state_list = new_sorted_current_state_list

                            #Manipulation of each fragment
                            for i in range(len(sorted_current_state_list)):

                                best_index = current_state.index(sorted_current_state[i])


                                next_frag_atom = int_result.atom[split_ind_list_1[best_index]:split_ind_list_2[best_index]]
                                next_frag_coord = int_result.coord[split_ind_list_1[best_index]:split_ind_list_2[best_index]]
                                EOF_value = int_result.eof_val[best_index]
                                category = int_result.temp_category[best_index]
                                temp_frag = int_result.temp_frag[best_index]
                                ref_frag_list_next = int_result.ref_frag_list_next_list[best_index]

                                first_list = [sublist[0] for sublist in all_state]

                                new_core_structure = create_new_core_structure(core,next_frag_atom,next_frag_coord,category,temp_frag,ref_frag_list_next)

                                current_cycle = core.add_cycle+core.current_cycle
                                current_branch_core = core.current_branch_core + all_add_branch_core[best_index]
                                current_branch_frag = core.current_branch_frag + all_add_branch_frag[best_index]
                                ref_frag_list_next = int_result.ref_frag_list_next_list[best_index]
                                param_for_next = int_result.param_for_next[best_index]
                                total_charge = int_result.total_charge[best_index]
                                transform_log_main = int_result.transform_log_main_list[best_index]
                                transform_log_side = int_result.transform_log_side_list[best_index]

                                if int_result.temp_frag[best_index] != 'sym24_XHD_0' and int_result.temp_frag[best_index] != 'sym24_XHD_1' and int_result.temp_frag[best_index] != 'sym24_XHD_2' and int_result.temp_frag[best_index] != 'sym24_XHD_3' and int_result.temp_frag[best_index] != 'sym24_XHD_4':
                                    c_for_next = int_result.c_for_next_list[best_index]
                                else:
                                    c_for_next = ['XCA_1']

                                if all_frag_category_val[best_index] != 'F1L' and all_frag_category_val[best_index] != 'C1L':

                                    score = execute_rollout_v2(new_core_structure,S,Q,P,total_charge,transform_log_main,transform_log_side,c_for_next,param_for_next,added_frag,ref_frag_list_next,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,evaluation_param,evaluation_method,desired_list_rollout,rollout_set,all_protein_cord,add_cycle,state_for_rollout,rol_des_category,des_frag_cycle,dir_core,cutoff_score_init,rollout_conf_inner,rollout_conf_inner_des,sigma_sum_cut)
                                else:
                                    score = current_state[best_index][1]

                                current_state[best_index][1] = score
                                current_state[best_index][2] = 1

                            each_state_0 = [[[next_state[0][0],sublist[0]]] + sublist[1:] if depth == 0 else [next_state[0]+[sublist[0]]] + sublist[1:] for sublist in current_state]
                            each_state = deepcopy(each_state_0)
                            ubc_ind = calc_UCB(current_state,cutoff_score)
                            best_for_opt = deepcopy(each_state[ubc_ind])
                            current_branch_core = core.current_branch_core + all_add_branch_core[ubc_ind]
                            current_branch_frag = core.current_branch_frag + all_add_branch_frag[ubc_ind]
                            ref_frag_list_next = int_result.ref_frag_list_next_list[ubc_ind]
                            param_for_next = int_result.param_for_next[ubc_ind]
                            total_charge = int_result.total_charge[ubc_ind]
                            transform_log_main = int_result.transform_log_main_list[ubc_ind]
                            transform_log_side = int_result.transform_log_side_list[ubc_ind]

                            if int_result.temp_frag[ubc_ind] != 'sym24_XHD_0' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_1' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_2' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_3' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_4':
                                c_for_next = int_result.c_for_next_list[ubc_ind]
                            else:
                                c_for_next = ['XCA_1']

                            if int_result.temp_frag[ubc_ind] != 'sym24_XHD_0' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_1' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_2' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_3' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_4':

                                for_opt_state = [best_for_opt]

                                plus_15_deg_frag = deepcopy(best_for_opt)
                                plus_15_deg_frag[0][-1][2] += 10
                                for_opt_state.append(plus_15_deg_frag)


                                minus_15_deg_frag = deepcopy(best_for_opt)
                                minus_15_deg_frag[0][-1][2] -= 10
                                for_opt_state.append(minus_15_deg_frag)

                                for_opt_score = []

                                for i in range(len(for_opt_state)):

                                    transforming_core_point_int = create_structure_from_strategy(for_opt_state,i,cycle,max_branch_core,max_branch_frag)
                                    if all_frag_category_val[ubc_ind] != 'F1L' and all_frag_category_val[ubc_ind] != 'C1L':
                                        score = execute_rollout_v2(transforming_core_point_int,S,Q,P,total_charge,transform_log_main,transform_log_side,c_for_next,param_for_next,added_frag,ref_frag_list_next,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,evaluation_param,evaluation_method,desired_list_rollout,rollout_set,all_protein_cord,add_cycle,state_for_rollout,rol_des_category,des_frag_cycle,dir_core,cutoff_score_init,rollout_conf_inner,rollout_conf_inner_des,sigma_sum_cut)
                                    else:
                                        atom_list = [x[0] for x in transforming_core_point_int]
                                        transforming_core_point_int_2= transforming_core_point_int[atom_list.index(core.Q_atom)+1:]
                                        LJ_array_sum,coulomb_array_sum,int_E_array_sum = calc_interaction_from_strategy(transforming_core_point_int_2,all_protein_cord)
                                        heavy_atom_list = [calc_heavy_atoms(transforming_core_point_int_2)]
                                        score = calc_evaluation_v2(LJ_array_sum,coulomb_array_sum,int_E_array_sum,evaluation_param,evaluation_method,heavy_atom_list)[0]

                                    for_opt_score.append(score)

                                opt_index = for_opt_score.index(min(for_opt_score))
                                each_state[ubc_ind][0][-1][2] = for_opt_state[opt_index][0][-1][2]
                                each_state[ubc_ind][1] = for_opt_score[opt_index]
                                current_state[ubc_ind][0][2] = for_opt_state[opt_index][0][-1][2]
                                current_state[ubc_ind][1] = for_opt_score[opt_index]

                        else:
                            int_state = [sublist for sublist in all_state if len(sublist[0])==(len(temp_state_0)+1) and sublist[0][:len(temp_state_0)]== temp_state_0]
                            current_state = [[a[0][-1],a[1],a[2]] for a in int_state]
                            each_state = [[[next_state[0][0],sublist[0]]] + sublist[1:] if depth == 0 else [next_state[0]+[sublist[0]]] + sublist[1:] for sublist in current_state]

                        current_score = [sublist[1] for sublist in current_state]
                        current_selection = [sublist[2] for sublist in current_state]
                        ubc_ind = calc_UCB(current_state,cutoff_score)


                        if current_score[ubc_ind] <= cutoff_score and current_selection[ubc_ind] < 10000:

                            #Update initiator
                            if first_elements_for_check == []:
                                all_state.extend(each_state)
                            first_list = [sublist[0] for sublist in all_state]
                            all_state_ind = first_list.index(next_state[0]+[current_state[ubc_ind][0]])

                            EOF_value = int_result.eof_val[ubc_ind]
                            next_frag_atom = int_result.atom[split_ind_list_1[ubc_ind]:split_ind_list_2[ubc_ind]]
                            next_frag_coord = int_result.coord[split_ind_list_1[ubc_ind]:split_ind_list_2[ubc_ind]]
                            category = int_result.temp_category[ubc_ind]
                            temp_frag = int_result.temp_frag[ubc_ind]
                            ref_frag_list_next = int_result.ref_frag_list_next_list[ubc_ind]
                            current_branch_core = core.current_branch_core + all_add_branch_core[ubc_ind]
                            current_branch_frag = core.current_branch_frag + all_add_branch_frag[ubc_ind]
                            ref_frag_list_next = int_result.ref_frag_list_next_list[ubc_ind]
                            param_for_next = int_result.param_for_next[ubc_ind]
                            total_charge = int_result.total_charge[ubc_ind]
                            transform_log_main = int_result.transform_log_main_list[ubc_ind]
                            transform_log_side = int_result.transform_log_side_list[ubc_ind]

                            if int_result.temp_frag[ubc_ind] != 'sym24_XHD_0' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_1' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_2' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_3' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_4':
                                c_for_next = int_result.c_for_next_list[ubc_ind]
                            else:
                                c_for_next = ['XCA_1']

                            core_name = each_state[ubc_ind][0][0]
                            frag_angle = each_state[ubc_ind][0][-1][-2]
                            score_info = "{:.2f}".format(current_score[ubc_ind])

                            infomation = '\tcore: '+str(core_name)+', depth: '+str(depth)+', fragment: '+str(temp_frag)+', dihedral: '+str(frag_angle)+', score: '+str(score_info)
                            logger_notime.info(infomation)

                            transforming_core_point_int = create_structure_from_strategy(each_state,ubc_ind,cycle,max_branch_core,max_branch_frag)

                            if EOF_value == False:
                                score = execute_rollout_v2(transforming_core_point_int,S,Q,P,total_charge,transform_log_main,transform_log_side,c_for_next,param_for_next,added_frag,ref_frag_list_next,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,evaluation_param,evaluation_method,desired_list_rollout,rollout_set,all_protein_cord,add_cycle,state_for_rollout,rol_des_category,des_frag_cycle,dir_core,cutoff_score_init,rollout_conf_inner,rollout_conf_inner_des,sigma_sum_cut)
                            elif EOF_value == True:
                                atom_list = [x[0] for x in transforming_core_point_int]
                                transforming_core_point_int_2= transforming_core_point_int[atom_list.index(core.Q_atom)+1:]
                                LJ_array_sum,coulomb_array_sum,int_E_array_sum = calc_interaction_from_strategy(transforming_core_point_int_2,all_protein_cord)
                                heavy_atom_list = [calc_heavy_atoms(transforming_core_point_int_2)]
                                score = calc_evaluation_v2(LJ_array_sum,coulomb_array_sum,int_E_array_sum,evaluation_param,evaluation_method,heavy_atom_list)[0]

                            if EOF_value == False:
                                all_state[all_state_ind][2] +=1
                                depth += 1
                                next_state = all_state[all_state_ind]

                                category = int_result.temp_category[ubc_ind]
                                temp_frag = int_result.temp_frag[ubc_ind]
                                ref_frag_list_next = int_result.ref_frag_list_next_list[ubc_ind]
                                current_cycle = core.add_cycle+core.current_cycle
                                current_branch_core = core.current_branch_core + all_add_branch_core[ubc_ind]
                                current_branch_frag = core.current_branch_frag + all_add_branch_frag[ubc_ind]
                                ref_frag_list_next = int_result.ref_frag_list_next_list[ubc_ind]
                                param_for_next = int_result.param_for_next[ubc_ind]
                                total_charge = int_result.total_charge[ubc_ind]
                                transform_log_main = int_result.transform_log_main_list[ubc_ind]

                                transform_log_side = int_result.transform_log_side_list[ubc_ind]

                                if int_result.temp_frag[ubc_ind] != 'sym24_XHD_0' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_1' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_2' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_3' and int_result.temp_frag[ubc_ind] != 'sym24_XHD_4':
                                    c_for_next = int_result.c_for_next_list[ubc_ind]

                                else:
                                    c_for_next = ['XCA_1']

                                core = Core_mol(transforming_core_point_int,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next)

                                infomation_for_param = '#LigX# @Sequence No.' +str(try_num+1)
                                logger_state.info(infomation_for_param)
                                param_log = param_for_logfile(transforming_core_point_int,current_cycle,current_branch_core,current_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next,temp_frag,next_frag_atom)
                                logger_state_notime.info('')
                                logger_state_notime.info(param_log)

                                logger_state_notime.info('\n################################LigXLigX###############################')

                            elif EOF_value == True:
                                infomation_for_param = '#LigX# @Sequence No.' +str(try_num+1)
                                logger_state.info(infomation_for_param)
                                param_log = param_for_logfile(transforming_core_point_int,current_cycle,current_branch_core,current_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next,temp_frag,next_frag_atom)
                                logger_state.info('')
                                logger_state_notime.info(param_log)
                                logger_state_notime.info('\n#LigX# EOF_value == True')


                                logger_state_notime.info('\n################################LigXLigX###############################')
                                logger_state_notime.info('################################LigXLigX###############################')
                                logger_state_notime.info('################################LigXLigX###############################')

                                count_generated +=1
                                try_num +=1
                                all_state[all_state_ind][2] +=10000
                                created = True

                        else:
                            first_list_to_escape = [sublist[0] for sublist in all_state]

                            escape_ind = first_list_to_escape.index(temp_state_0)
                            all_state[escape_ind][1] = 100
                            all_state[escape_ind][2] = 10000
                            depth = 0
                            try_num +=1
                            created = False

                            green_print('Not generated in this sequence. Node pruning')

                            logger_notime.info('Not generated in this sequence. Removed by cutoff')
                            state_for_logfile = '\tstate:'+str(each_state[ubc_ind])
                            logger_notime.info('')
                            logger_notime.info(state_for_logfile)
                            logger_notime.info('\n################################LigXLigX###############################')

                            break

                    elif possibility == False:
                        temp_state_0 = temp_state[0]
                        first_list_to_escape = [sublist[0] for sublist in all_state]

                        escape_ind = first_list_to_escape.index(temp_state_0)
                        all_state[escape_ind][1] = 100
                        all_state[escape_ind][2] = 10000
                        depth = 0
                        green_print('Not generated in this sequence. Possibility: False. Node pruning')

                        logger_notime.info('Not generated in this sequence. Possibility: False')
                        logger_notime.info('\n################################LigXLigX###############################')
                        try_num +=1
                        created = False
                        break

            #Save Data
            if created == False:
                continue
            depth = 0




            transforming_core_point_int = create_structure_from_strategy(each_state,ubc_ind,cycle,max_branch_core,max_branch_frag)
            LJ_list = [x[0] for x in LJ_array_sum]
            coulomb_list = [x[0] for x in coulomb_array_sum]
            int_E_list = [x[0] for x in int_E_array_sum]
            heavy_atom = heavy_atom_list[0]
            if heavy_atom == 0:
                heavy_atom = 1

            LE_LJ =[float(x)/int(heavy_atom) for x in LJ_list]
            LE_LJ_coulomb = []
            for a,b in zip(LE_LJ,coulomb_list):
                LE_LJ_coulomb.append(a+float(b))

            transforming_core_point =atomname_change_final(transforming_core_point_int)
            final_total_atom = [str(int(len(transforming_core_point)))]
            each_name = ['ID_'+str(count_generated)]
            temp_structure = []

            temp_structure.append(final_total_atom)
            temp_structure.append(each_name)

            for j in transforming_core_point:
                float_point =[round(value.item(), 4) if isinstance(value, np.float64) else value for value in j]
                temp_structure.append(float_point)

            smiles = get_smiles(temp_structure,total_charge,self.smiles_filter)

            if smiles != 'no_SMILES' and current_cycle-1 >= des_frag_cycle:

                all_result.append(final_total_atom)
                all_result.append(each_name)

                for j in transforming_core_point:
                    float_point =[round(value.item(), 4) if isinstance(value, np.float64) else value for value in j]
                    all_result.append(float_point)

                all_result_name = 'result_all_total_conf_'+str(count_generated)
                save_csv(all_result,all_result_name,result_dir)
                if count_generated >1:
                    remove_file = result_dir+'/result_all_total_conf_'+str(count_generated-1)+'.xyz'
                    os.remove(remove_file)

                each_result = [count_generated]

                for i in range(num_protein_structures):
                    each_result.append(LJ_list[i])

                LJ_average = statistics.mean(LJ_list)
                each_result.append(LJ_average)

                for i in range(num_protein_structures):
                    each_result.append(coulomb_list[i])
                coulomb_average = statistics.mean(coulomb_list)
                each_result.append(coulomb_average)

                for i in range(num_protein_structures):
                    each_result.append(int_E_list[i])
                int_E_average = statistics.mean(int_E_list)
                each_result.append(int_E_average)

                for i in range(num_protein_structures):
                    each_result.append(float(LJ_list[i])/int(heavy_atom))
                LJ_LE_average = statistics.mean(LJ_list)/int(heavy_atom)
                each_result.append(LJ_LE_average)

                for i in range(num_protein_structures):
                    each_result.append(float(int_E_list[i])/int(heavy_atom))
                int_E_LE_average = statistics.mean(int_E_list)/int(heavy_atom)
                each_result.append(int_E_LE_average)

                for i in range(num_protein_structures):
                    each_result.append(LE_LJ_coulomb[i])

                LE_LJ_coulomb_LJ_average = statistics.mean(LE_LJ_coulomb)
                each_result.append(LE_LJ_coulomb_LJ_average)

                each_result.append(each_state[ubc_ind])
                each_result.append(smiles)

                all_score_strategy.append(each_result)

                score_name = 'result_score'
                save_csv_2(all_score_strategy,score_name,result_dir)

                score_for_info =  "{:.3f}".format(score)
                generated_time = time.time()

                elapsed_time = generated_time-start

                performance = calc_performance(count_generated,elapsed_time)

                generated_info = '#LigX# '+str(count_generated)+' conf generated (' +performance +'), score = '+str(score_for_info)

                red_print(generated_info)

                structure_info = '\tName: '+str(each_name[0])+', SMILES: '+str(smiles)
                state_for_logfile = '\tstate:'+str(each_state[ubc_ind])
                logger_notime.info('')
                logger_notime.info(state_for_logfile)
                logger_notime.info('')
                logger.info(generated_info)
                logger_notime.info(structure_info)
                logger_notime.info('\n################################LigXLigX###############################')

                if smiles != 'without_SMILES_conversion':
                    for_fig = Chem.MolFromSmiles(smiles)
                    Draw.MolToFile(for_fig,image_dir+'/'+str(each_name[0])+'.png',size=(900, 900))

            else:
                green_print('Removed by SMILES filter')
                state_for_logfile = '\tstate:'+str(each_state[ubc_ind])
                logger_notime.info('')
                logger_notime.info(state_for_logfile)
                logger_notime.info('Removed by SMILES filter')
                logger_notime.info('\n################################LigXLigX###############################')
                count_generated-=1

            if (count_generated)//count_for_data_processing > (count_generated-1)//count_for_data_processing and smiles != 'no_SMILES':


                if os.path.isdir(result_dir+'/result_conf_'+str(count_for_top)):
                    shutil.rmtree(result_dir+'/result_conf_'+str(count_for_top))
                    del count_for_top

                count_for_top = deepcopy(count_generated)

                print('creating top at conf '+ str(count_for_top))
                os.mkdir(result_dir+'/result_conf_'+str(count_for_top))
                index_list = get_reconstract_strategy(num_protein_structures,all_score_strategy)
                all_score_strategy_cut = all_score_strategy[1:]
                if self.clean == True:
                    subprocess.call([clean_path])


                for index_name,index in index_list[1:-1]:
                    sorted_all_score_strategy_cut = sorted(all_score_strategy_cut, key=lambda x: x[index])
                    count_sel = 0
                    ranking_list = []
                    each_ranking_dir = result_dir+'/result_conf_'+str(count_generated)+'/'+index_name
                    os.mkdir(each_ranking_dir)

                    for sorted_list in sorted_all_score_strategy_cut:

                        count_sel += 1
                        if count_sel >result_conf:
                            break
                        depth = 0
                        ID_for_reult =sorted_list[0]

                        each_strategy = sorted_list[-2]
                        each_strategy_fin = [each_strategy]

                        final_charge = get_core_charge(each_strategy_fin[0][0][0])


                        for i in each_strategy_fin[0][0][1:]:
                            final_charge += int(i[-1])
                            #print(i[-1])

                        transforming_core_point = create_structure_from_strategy(each_strategy_fin,0,cycle,max_branch_core,max_branch_frag)

                        transforming_core_point =atomname_change_final(transforming_core_point)
                        final_total_atom = [str(int(len(transforming_core_point)))]
                        rank_num = str(count_sel).zfill(4)
                        each_name = ['rank_'+rank_num+'_'+index_name+'_ID_'+str(ID_for_reult)]

                        each_list = []

                        each_list.append(final_total_atom)
                        each_list.append(each_name)
                        ranking_list.append(final_total_atom)
                        ranking_list.append(each_name)

                        for j in transforming_core_point:
                            float_point =[round(value.item(), 4) if isinstance(value, np.float64) else value for value in j]
                            each_list.append(float_point)
                            ranking_list.append(float_point)

                        each_result ='result_conf_'+str(count_generated)+'/'+index_name+'/'+str(each_name[0])+'/'+str(each_name[0])
                        os.mkdir(each_ranking_dir+'/'+str(each_name[0]))
                        dir_charge = each_ranking_dir+'/'+str(each_name[0])

                        save_csv(each_list,each_result,result_dir)
                        save_charge(final_charge,dir_charge)
                    result_name ='result_conf_'+str(count_generated)+'/'+'ranking_'+index_name+'_top'
                    save_csv(ranking_list,result_name,result_dir)


        if os.path.isdir(result_dir+'/result_conf_'+str(count_for_top)):
            shutil.rmtree(result_dir+'/result_conf_'+str(count_for_top))
            del count_for_top

        count_for_top = deepcopy(count_generated)

        print('creating top at conf '+ str(count_generated))
        os.mkdir(result_dir+'/result_conf_'+str(count_generated))
        index_list = get_reconstract_strategy(num_protein_structures,all_score_strategy)
        all_score_strategy_cut = all_score_strategy[1:]
        if self.clean == True:
            subprocess.call([clean_path])

        for index_name,index in index_list[1:-1]:
            sorted_all_score_strategy_cut = sorted(all_score_strategy_cut, key=lambda x: x[index])
            count_sel = 0
            ranking_list = []
            each_ranking_dir = result_dir+'/result_conf_'+str(count_generated)+'/'+index_name
            os.mkdir(each_ranking_dir)

            for sorted_list in sorted_all_score_strategy_cut:

                count_sel += 1
                if count_sel >result_conf:
                    break
                depth = 0
                ID_for_reult =sorted_list[0]

                each_strategy = sorted_list[-2]
                each_strategy_fin = [each_strategy]

                final_charge = get_core_charge(each_strategy_fin[0][0][0])

                for i in each_strategy_fin[0][0][1:]:
                    final_charge += int(i[-1])

                transforming_core_point = create_structure_from_strategy(each_strategy_fin,0,cycle,max_branch_core,max_branch_frag)

                transforming_core_point =atomname_change_final(transforming_core_point)
                final_total_atom = [str(int(len(transforming_core_point)))]
                rank_num = str(count_sel).zfill(4)
                each_name = ['rank_'+rank_num+'_'+index_name+'_ID_'+str(ID_for_reult)]

                each_list = []

                each_list.append(final_total_atom)
                each_list.append(each_name)
                ranking_list.append(final_total_atom)
                ranking_list.append(each_name)

                for j in transforming_core_point:
                    float_point =[round(value.item(), 4) if isinstance(value, np.float64) else value for value in j]
                    each_list.append(float_point)
                    ranking_list.append(float_point)

                each_result ='result_conf_'+str(count_generated)+'/'+index_name+'/'+str(each_name[0])+'/'+str(each_name[0])
                os.mkdir(each_ranking_dir+'/'+str(each_name[0]))
                dir_charge = each_ranking_dir+'/'+str(each_name[0])

                save_csv(each_list,each_result,result_dir)
                save_charge(final_charge,dir_charge)
            result_name ='result_conf_'+str(count_generated)+'/'+'ranking_'+index_name+'_top'
            save_csv(ranking_list,result_name,result_dir)
        print('Normal termination')
        logger.info('Normal termination')
        logger_notime.info('\n################################LigXLigX###############################')

##basic class for LigX

class Fragments:
    def __init__(self,frag_category_val,transform_dict):
        pattern_in_category = transform_dict[frag_category_val]
        combined_atom_list = []
        combined_coord_list = []
        combined_charge_list = []
        each_frag_list = []
        len_list = []
        cent_atom_list = []
        param_list = []
        category_list = []
        frag_name_list = []
        add_branch_core_list = []
        add_branch_frag_list = []
        all_ref_frag_list = []
        frag_category_val_list = []
        c_for_next_list = []
        for each_frag in pattern_in_category:
            frag_cord_list,frag_name = read_xyz_data(each_frag[0])
            atom_list,coord_list = separate_atom_and_coord(frag_cord_list)
            param,category = read_param(each_frag[1])
            ref_frag_list = ref_frag_from_param(each_frag[1])
            if ref_frag_list != None:
                ref_frag = ref_frag_list[0]
            else:
                ref_frag = None
            angle = each_frag[2]
            charge = int(each_frag[3])
            if frag_name == 'sym24_H1L':
                if category == ['link_c']:
                    coord_list = [[0,0,0]]
                    atom_list = ['D']
                    key_elements = ['D']
                elif category == ['link_f']:
                    coord_list = [[0,0,0]]
                    atom_list = ['H']
                    key_elements = ['D']
            else:
                key_elements = extract_key_elements(atom_list,category,param,ref_frag)
                coord_list = initialize_fragment(atom_list,coord_list,key_elements,angle)
            for i in atom_list:
                combined_atom_list.append(i)
            for i in coord_list:
                combined_coord_list.append(i)
            combined_charge_list.append(charge)
            each_frag_list.append(each_frag)
            len_list.append(len(atom_list))
            cent_atom_list.append(key_elements[0])
            param_list.append(param)
            category_list.append(category)
            frag_name_list.append(frag_name)
            all_ref_frag_list.append(ref_frag_list)
            frag_category_val_list.append(frag_category_val)
            if frag_category_val == 'F3L':
                add_branch_core_list.append(0)
                add_branch_frag_list.append(1)
            elif frag_category_val == 'C3L':
                add_branch_core_list.append(1)
                add_branch_frag_list.append(0)
            else:
                add_branch_core_list.append(0)
                add_branch_frag_list.append(0)
            c_for_next = []
            f_atom_list = ['YCA','YO','YNA','YN3','YN2','YS','YSA']
            for x in ['XC3_','XC3g_','XC2_','XC1_','XCA_']:
                for n in range(1,4):
                    X_name_next = x+str(n)
                    if any(X_name_next in element for element in atom_list):
                        c_for_next.append(X_name_next)
            if category == ['link_c']:
                for f_atom in f_atom_list:
                    if any(f_atom in element for element in atom_list):
                        c_for_next = []
                        break
            c_for_next_list.append(c_for_next)
        self.atom = combined_atom_list
        self.coord = combined_coord_list
        self.charge = combined_charge_list
        self.each_frag = each_frag_list
        self.len_list = len_list
        self.central_atom = cent_atom_list
        self.param = param_list
        self.category = category_list
        self.ref_frag_list = all_ref_frag_list
        self.frag_name = frag_name_list
        self.add_branch_core = add_branch_core_list
        self.add_branch_frag = add_branch_frag_list
        self.frag_category_val = frag_category_val_list
        self.c_for_next = c_for_next_list



class desFragments:
    def __init__(self,desired_list,des_category):
        pattern_in_category =get_desired_list(desired_list)
        combined_atom_list = []
        combined_coord_list = []
        combined_charge_list = []
        each_frag_list = []
        len_list = []
        cent_atom_list = []
        param_list = []
        category_list = []
        frag_name_list = []
        add_branch_core_list = []
        add_branch_frag_list = []
        all_ref_frag_list = []
        frag_category_val_list = []
        c_for_next_list = []
        for each_frag in pattern_in_category:
            frag_cord_list,frag_name = read_xyz_data(each_frag[0])
            atom_list,coord_list = separate_atom_and_coord(frag_cord_list)
            param,category = read_param(each_frag[1])
            ref_frag_list = ref_frag_from_param(each_frag[1])
            if ref_frag_list != None:
                ref_frag = ref_frag_list[0]
            else:
                ref_frag = None
            angle = each_frag[2]
            charge = int(each_frag[3])
            if frag_name == 'sym24_H1L':
                if category == ['link_c']:
                    coord_list = [0,0,0]
                    atom_list = [['D']]
                    key_elements = ['D']
                elif category == ['link_f']:
                    coord_list = [[0,0,0]]
                    atom_list = ['H']
                    key_elements = ['D']
            else:
                key_elements = extract_key_elements(atom_list,category,param,ref_frag)
                coord_list = initialize_fragment(atom_list,coord_list,key_elements,angle)
            for i in atom_list:
                combined_atom_list.append(i)
            for i in coord_list:
                combined_coord_list.append(i)
            combined_charge_list.append(charge)
            each_frag_list.append(each_frag)
            len_list.append(len(atom_list))
            cent_atom_list.append(key_elements[0])
            param_list.append(param)
            category_list.append(category)
            frag_name_list.append(frag_name)
            if des_category == 'F3L':
                add_branch_core_list.append(0)
                add_branch_frag_list.append(1)
            elif des_category == 'C3L':
                add_branch_core_list.append(1)
                add_branch_frag_list.append(0)
            else:
                add_branch_core_list.append(0)
                add_branch_frag_list.append(0)
            all_ref_frag_list.append(ref_frag_list)
            frag_category_val_list.append(des_category)
            c_for_next = []
            f_atom_list = ['YCA','YO','YNA','YN3','YN2','YS','YSA']
            for x in ['XC3_','XC3g_','XC2_','XC1_','XCA_']:
                for n in range(1,4):
                    X_name_next = x+str(n)
                    if any(X_name_next in element for element in atom_list):
                        c_for_next.append(X_name_next)
            if category == ['link_c']:
                for f_atom in f_atom_list:
                    if any(f_atom in element for element in atom_list):
                        c_for_next = []
                        break
            c_for_next_list.append(c_for_next)
        self.atom = combined_atom_list
        self.coord = combined_coord_list
        self.charge = combined_charge_list
        self.each_frag = each_frag_list
        self.len_list = len_list
        self.central_atom = cent_atom_list
        self.param = param_list
        self.category = category_list
        self.ref_frag_list = all_ref_frag_list
        self.frag_name = frag_name_list
        self.add_branch_core = add_branch_core_list
        self.add_branch_frag = add_branch_frag_list
        self.frag_category_val = frag_category_val_list
        self.c_for_next = c_for_next_list


class single_Fragment:
    def __init__(self,each_frag):
        combined_atom_list = []
        combined_coord_list = []
        combined_charge_list = []
        each_frag_list = []
        len_list = []
        cent_atom_list = []
        param_list = []
        category_list = []
        frag_name_list = []
        add_branch_core_list = []
        add_branch_frag_list = []
        all_ref_frag_list = []
        frag_category_val = each_frag[0].split('/')[1]
        c_for_next_list = []
        frag_cord_list,frag_name = read_xyz_data(each_frag[0])
        atom_list,coord_list = separate_atom_and_coord(frag_cord_list)
        param,category = read_param(each_frag[1])
        ref_frag_list = ref_frag_from_param(each_frag[1])
        if ref_frag_list != None:
            ref_frag = ref_frag_list[0]
        else:
            ref_frag = None
        angle = each_frag[2]
        charge = int(each_frag[3])
        if frag_name == 'sym24_H1L':
            if category == ['link_c']:
                coord_list = [[0,0,0]]
                atom_list = ['D']
                key_elements = ['D']
            elif category == ['link_f']:
                coord_list = [[0,0,0]]
                atom_list = ['H']
                key_elements = ['D']
        else:
            key_elements = extract_key_elements(atom_list,category,param,ref_frag)
            coord_list = initialize_fragment(atom_list,coord_list,key_elements,angle)
        for i in atom_list:
            combined_atom_list.append(i)
        for i in coord_list:
            combined_coord_list.append(i)
        combined_charge_list.append(charge)
        each_frag_list.append(each_frag)
        len_list.append(len(atom_list))
        cent_atom_list.append(key_elements[0])
        param_list.append(param)
        category_list.append(category)
        frag_name_list.append(frag_name)
        all_ref_frag_list.append(ref_frag_list)
        if frag_category_val == 'F3L':
            add_branch_core_list.append(0)
            add_branch_frag_list.append(1)
        elif frag_category_val == 'C3L':
            add_branch_core_list.append(1)
            add_branch_frag_list.append(0)
        else:
            add_branch_core_list.append(0)
            add_branch_frag_list.append(0)
        c_for_next = []
        f_atom_list = ['YCA','YO','YNA','YN3','YN2','YS','YSA']
        for x in ['XC3_','XC3g_','XC2_','XC1_','XCA_']:
            for n in range(1,4):
                X_name_next = x+str(n)
                if any(X_name_next in element for element in atom_list):
                    c_for_next.append(X_name_next)
        if category == ['link_c']:
            for f_atom in f_atom_list:
                if any(f_atom in element for element in atom_list):
                    c_for_next = []
                    break
        c_for_next_list.append(c_for_next)
        self.atom = combined_atom_list
        self.coord = combined_coord_list
        self.charge = combined_charge_list
        self.each_frag = each_frag_list
        self.len_list = len_list
        self.central_atom = cent_atom_list
        self.param = param_list
        self.category = category_list
        self.ref_frag_list = all_ref_frag_list
        self.frag_name = frag_name_list
        self.add_branch_core = add_branch_core_list
        self.add_branch_frag = add_branch_frag_list
        self.c_for_next = c_for_next_list

class Core_mol:
    def __init__(self,core_coord,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next):

        atom_list,coord_list = separate_atom_and_coord(core_coord)



        if any('NQ' in element for element in atom_list):
            Q = 'NQ'
            QE ='N'
        elif any('OQ' in element for element in atom_list):
            Q = 'OQ'
            QE ='O'
        elif any('HQ' in element for element in atom_list):
            Q = 'HQ'
            QE ='H'
        elif any('DQ' in element for element in atom_list):
            Q = 'DQ'
            QE ='D'
        elif any('SQ' in element for element in atom_list):
            Q = 'SQ'
            QE ='S'
        elif any('CQ' in element for element in atom_list):
            Q = 'CQ'
            QE ='C'
        else:
            print('no Q atom')
            sys.exit()

        self.Q_atom = Q
        self.QE_atom = QE

        if any('XC5' in element for element in atom_list):
            self.add_cycle = 1
            self.operation = ['XFO']
            central_atom = 'XC5'
            for r in ['RH_1','RC_1','RO_1','RN_1','RS_1','RCl_1','RBr_1','RI_1','RP_1','RD_1','RF_']:
                if any(r in element for element in atom_list):
                    reference_atom = r
                    break
            virtual_atom = 'VC_1'


        elif any('X4' in element for element in atom_list):
            self.add_cycle = 1
            self.operation = ['XHD']
            central_atom = 'X4'
            reference_atom = 'B4'
            virtual_atom = 'A4'

        elif any('VLC' in element for element in atom_list):

            reference_atom = ref_frag_list_next[-1]
            central_atom = param_for_next[-1]
            for v in ['VLC','VLB','VLA']:
                if any(v in element for element in atom_list):
                    virtual_atom = v
                    break
                else:
                    continue

            self.add_cycle = 0
            self.operation = ['C1LB']

        elif any('VLB' in element for element in atom_list):

            reference_atom = ref_frag_list_next[-1]
            central_atom = param_for_next[-1]
            for v in ['VLC','VLB','VLA']:
                if any(v in element for element in atom_list):
                    virtual_atom = v
                    break
                else:
                    continue
            if cycle > current_cycle:
                if current_branch_core < max_branch_core:
                    self.add_cycle = 1
                    self.operation = ['C1L','C2L','C3L']

                elif current_branch_core == max_branch_core:
                    self.add_cycle = 1
                    self.operation = ['C1L','C2L']

            elif cycle == current_cycle:
                self.add_cycle = 1
                self.operation = ['C1L']

        elif any('VC_3' in element for element in atom_list):

            for n in reversed(range(1,4)):
                for x in ['XC3_','XC3g_','XC2_','XC1_','XCA_']:
                    X_name = x+str(n)

                    if any(X_name in element for element in atom_list):
                        central_atom = X_name
                        n_at_front = copy.deepcopy(n)
                        break
                    else:
                        continue
                else:
                    continue
                break

            for r in ['RH_','RC_','RO_','RN_','RS_','RCl_','RBr_','RI_','RP_','RD_','RF_']:
                R_name = r+str(n_at_front)
                if any(R_name in element for element in atom_list):
                    reference_atom = R_name
                    break
                else:
                    continue
            for vc in ['VC_3','VC_2','VC_1']:
                if any(vc in element for element in atom_list):
                    virtual_atom = vc
                    break


            self.add_cycle = 0
            self.operation = ['F2LB']

        elif any('VC_2' in element for element in atom_list):

            for n in reversed(range(1,4)):
                for x in ['XC3_','XC3g_','XC2_','XC1_','XCA_']:
                    X_name = x+str(n)

                    if any(X_name in element for element in atom_list):
                        central_atom = X_name
                        n_at_front = copy.deepcopy(n)
                        break
                    else:
                        continue
                else:
                    continue
                break

            for r in ['RH_','RC_','RO_','RN_','RS_','RCl_','RBr_','RI_','RP_','RD_','RF_']:
                R_name = r+str(n_at_front)
                if any(R_name in element for element in atom_list):
                    reference_atom = R_name
                    break
                else:
                    continue
            for vc in ['VC_3','VC_2','VC_1']:
                if any(vc in element for element in atom_list):
                    virtual_atom = vc
                    break

            if cycle > current_cycle:
                if current_branch_frag < max_branch_frag:
                    self.add_cycle = 1
                    self.operation = ['F1L','F2L','F3L']

                elif current_branch_frag == max_branch_frag:
                    self.add_cycle = 1
                    self.operation = ['F1L','F2L']

            elif cycle == current_cycle:
                self.add_cycle = 1
                self.operation = ['F1L']

        key_elements = [central_atom,reference_atom,virtual_atom]

        self.central_atom = central_atom

        inv_rotation_matrix_4th,inv_rotation_matrix_3rd, inv_rotation_matrix_2nd, inv_rotation_matrix_1st, coord_central = initialize_core(atom_list,coord_list,key_elements)

        self.matrix_4th = inv_rotation_matrix_4th
        self.matrix_3rd = inv_rotation_matrix_3rd
        self.matrix_2nd = inv_rotation_matrix_2nd
        self.matrix_1st = inv_rotation_matrix_1st
        self.central_coord = coord_central

        self.atom = atom_list
        self.coord = coord_list

        self.total_charge= total_charge
        self.transform_log_main = transform_log_main
        self.transform_log_side = transform_log_side
        self.c_for_next = c_for_next
        self.current_cycle = current_cycle
        self.current_branch_core = current_branch_core
        self.current_branch_frag = current_branch_frag
        self.param_for_next = param_for_next
        self.ref_frag_list_next = ref_frag_list_next

class link_result:
    def __init__(self,Core_mol,all_atom,all_coord,all_charge,all_central_atom,split_ind_list_1,split_ind_list_2,all_category,all_param,all_len_list,all_frag,all_protein_cord,evaluation_param,evaluation_method,all_ref_frag_list,frag_c_for_next,sigma_sum_cut):

        all_transform_log = []

        for i in all_category:

            if i == ['link_c']:
                if len(Core_mol.param_for_next)==2:
                    all_transform_log.append(Core_mol.transform_log_side)

                elif len(Core_mol.param_for_next)==1:
                    all_transform_log.append(Core_mol.transform_log_main)

            elif i == ['link_f']:
                if len(Core_mol.c_for_next)==3:
                    all_transform_log.append(Core_mol.transform_log_side)
                elif len(Core_mol.c_for_next)==2 or len(Core_mol.c_for_next)==1:
                    all_transform_log.append(Core_mol.transform_log_main)


        if all('sym24_XHD_0' not in x and 'sym24_XHD_1' not in x and 'sym24_XHD_2' not in x and 'sym24_XHD_3' not in x and 'sym24_XHD_4' not in x and 'XFO' not in x for x in all_frag):

            check_list = [check_permission_previous(frag_central_atom,Core_mol.central_atom) != False and check_permission_combination(frag_central_atom,Core_mol.central_atom,transform_log) != False for frag_central_atom,transform_log in zip(all_central_atom,all_transform_log)]

        else:
            check_list = [True for _ in all_frag]
        dist_list = [[float(0),float(0),float(get_distance(frag_central_atom,Core_mol.central_atom))] if check_permission_previous(frag_central_atom,Core_mol.central_atom) != False else [np.nan,np.nan,np.nan] for frag_central_atom in all_central_atom]

        dist_array = np.array([row for row, count in zip(dist_list, all_len_list) for _ in range(count)])


        rot_matrix = np.dot(Core_mol.matrix_4th,np.dot(Core_mol.matrix_3rd,np.dot(Core_mol.matrix_2nd,Core_mol.matrix_1st )))

        coord_array = np.dot(rot_matrix, (np.array(all_coord)-dist_array).T).T + Core_mol.central_coord

        all_EOF = []
        transform_log_main_list = []
        transform_log_side_list = []
        c_for_next_list = []
        ref_frag_list_next_list = []
        temp_frag_list = []
        temp_category_list = []
        param_for_next_list = []


        new_total_charge= (np.array(all_charge)+Core_mol.total_charge).tolist()


        for temp_param, temp_category, temp_frag,ref_frag_list,temp_c_for_next in zip(all_param,all_category,all_frag,all_ref_frag_list,frag_c_for_next):

            transform_log_main = copy.deepcopy(Core_mol.transform_log_main)
            transform_log_side = copy.deepcopy(Core_mol.transform_log_side)

            if temp_category ==['link_c']:

                if temp_frag != 'sym24_H1L':

                    ### C→C or C→N
                    if temp_param == []:
                         ### C→C　

                        if len(Core_mol.param_for_next)==2:

                            transform_log_side.append(temp_c_for_next[-1])
                        elif len(Core_mol.param_for_next)==1:
                            if len(temp_c_for_next) ==3:
                                transform_log_side.append(temp_c_for_next[-1])
                                transform_log_main.append(temp_c_for_next[-2])
                            elif len(temp_c_for_next) ==2 or len(temp_c_for_next) ==1:
                                transform_log_main.append(temp_c_for_next[-1])

                        ref_frag_list_next = copy.deepcopy(Core_mol.ref_frag_list_next)
                        param_for_next = copy.deepcopy(Core_mol.param_for_next)
                        del param_for_next[-1]

                        del ref_frag_list_next[-1]
                    else:
                        ### C→N　
                        param_for_next = temp_param[1:]

                        if len(temp_c_for_next) ==3:
                            transform_log_side.append(temp_param[-1])
                        elif len(temp_c_for_next) ==2 or len(temp_c_for_next) ==1 or len(temp_c_for_next) ==0:
                            if len(param_for_next)==2:
                                transform_log_main.append(temp_param[-2])
                                transform_log_side.append(temp_param[-1])
                            elif len(param_for_next)==1:
                                transform_log_main.append(temp_param[-1])

                        ref_frag_list_next = [param_for_next[-1]]

                    c_for_next = copy.deepcopy(temp_c_for_next)

                if temp_frag == 'sym24_H1L':
                    ref_frag_list_next = copy.deepcopy(Core_mol.ref_frag_list_next)
                    param_for_next = copy.deepcopy(Core_mol.param_for_next)
                    c_for_next = copy.deepcopy(temp_c_for_next)

                    del param_for_next[-1]
                    del ref_frag_list_next[-1]

            elif temp_category ==['link_f']:

                if temp_frag == 'XFO':
                    if len(Core_mol.c_for_next) ==3:
                        transform_log_side.append(temp_param[-1])
                    elif len(Core_mol.c_for_next) ==2 or len(Core_mol.c_for_next) ==1:
                        if len(Core_mol.param_for_next)==2:
                            transform_log_main.append(temp_param[-2])
                            transform_log_side.append(temp_param[-1])
                        elif len(Core_mol.param_for_next)==1:
                            transform_log_main.append(temp_param[-1])
                    ref_frag_list_next=[]
                    param_for_next = temp_param[1:]
                    c_for_next = copy.deepcopy(temp_c_for_next)

                else:

                    if temp_frag == 'sym24_XHD_0' or temp_frag == 'sym24_XHD_1' or temp_frag == 'sym24_XHD_2' or temp_frag == 'sym24_XHD_3' or temp_frag == 'sym24_XHD_4':
                        param_for_next = temp_param[0:]

                    else:
                        param_for_next = temp_param[1:]


                    if temp_frag == 'sym24_XHD_0' or temp_frag == 'sym24_XHD_1' or temp_frag == 'sym24_XHD_2' or temp_frag == 'sym24_XHD_3' or temp_frag == 'sym24_XHD_4':
                        # link hudrazone
                        transform_log_main.append(temp_param[-1])
                        ref_frag_list_next = copy.deepcopy(ref_frag_list)
                        c_for_next = copy.deepcopy(temp_c_for_next)

                        del ref_frag_list_next[0]


                    else:

                        if temp_param[0] =='dummy':

                        #N→C
                            if len(param_for_next)==2:
                                transform_log_side.append(temp_param[-1])
                            elif len(param_for_next)==1:
                                if len(temp_param) ==3:
                                    transform_log_side.append(temp_param[-1])
                                    transform_log_main.append(temp_param[-2])
                                elif len(temp_param) ==2 or len(temp_param) ==1:
                                    transform_log_main.append(temp_param[-1])
                            ref_frag_list_next = copy.deepcopy(ref_frag_list)

                            del param_for_next[-1]
                            del ref_frag_list_next[-1]
                            c_for_next = copy.deepcopy(temp_c_for_next)

                        else:
                        #N→N

                            if len(Core_mol.c_for_next) ==3:
                                transform_log_side.append(temp_param[-1])
                            elif len(Core_mol.c_for_next) ==2 or len(Core_mol.c_for_next) ==1:
                                if len(param_for_next)==2:
                                    transform_log_main.append(temp_param[-2])
                                    transform_log_side.append(temp_param[-1])
                                elif len(param_for_next)==1:
                                    transform_log_main.append(temp_param[-1])
                            ref_frag_list_next = copy.deepcopy(ref_frag_list)
                            del ref_frag_list_next[0]
                            c_for_next = copy.deepcopy(temp_c_for_next)



            if temp_frag == 'sym24_H1L' or temp_frag == 'carboxylate':
                all_EOF.append(True)
            else:
                all_EOF.append(False)

            transform_log_main_list.append(transform_log_main)
            transform_log_side_list.append(transform_log_side)
            c_for_next_list.append(c_for_next)
            ref_frag_list_next_list.append(ref_frag_list_next)
            temp_frag_list.append(temp_frag)
            temp_category_list.append(temp_category)
            param_for_next_list.append(param_for_next)

        #Calculation of interact E
        remaining_frag_atom = Core_mol.atom[Core_mol.atom.index(Core_mol.Q_atom)+1:]
        remaining_frag_coord = Core_mol.coord[Core_mol.atom.index(Core_mol.Q_atom)+1:]

        if split_ind_list_1 != None or split_ind_list_2 != None:

            collision_check = check_collision_lig_fin_v2_1(Core_mol.atom, Core_mol.coord, all_atom, coord_array,split_ind_list_1,split_ind_list_2,temp_frag_list,sigma_sum_cut)

            heavy_atom_remaining = sum(1 for row in remaining_frag_atom if row != 'H' and row != 'D' and row != 'VLA' and row != 'VLB' and row != 'VLC' and row != 'VC_1' and row != 'VC_2' and row != 'VC_3')
            heavy_atom_list = [heavy_atom_remaining + sum(1 for row in each if row != 'H' and row != 'D' and row != 'VLA' and row != 'VLB' and row != 'VLC' and row != 'VC_1' and row != 'VC_2' and row != 'VC_3') for each in [all_atom[a:b] for a,b in zip(split_ind_list_1,split_ind_list_2)]]

            LJ_list_remaining,coulomb_list_remaining,int_E_list_remaining = calc_interaction_remaining(remaining_frag_atom,remaining_frag_coord, all_protein_cord)
            LJ_list,coulomb_list,int_E_list = calc_interaction_additional_frag(all_atom,coord_array, all_protein_cord,split_ind_list_1,split_ind_list_2)
        else:
            collision_check = check_collision_lig_fin_v2_1(Core_mol.atom, Core_mol.coord, all_atom, coord_array,split_ind_list_1,split_ind_list_2,temp_frag_list,sigma_sum_cut)

            heavy_atom_remaining = sum(1 for row in remaining_frag_atom if row != 'H' and row != 'D' and row != 'VLA' and row != 'VLB' and row != 'VLC' and row != 'VC_1' and row != 'VC_2' and row != 'VC_3')
            heavy_atom_list = [heavy_atom_remaining + sum(1 for row in each if row != 'H' and row != 'D' and row != 'VLA' and row != 'VLB' and row != 'VLC' and row != 'VC_1' and row != 'VC_2' and row != 'VC_3') for each in [all_atom[:]]]

            LJ_list_remaining,coulomb_list_remaining,int_E_list_remaining = calc_interaction_remaining(remaining_frag_atom,remaining_frag_coord, all_protein_cord)
            LJ_list,coulomb_list,int_E_list = calc_interaction_additional_frag(all_atom,coord_array, all_protein_cord,split_ind_list_1,split_ind_list_2)

        LJ_array_sum = np.array(LJ_list_remaining)+np.array(LJ_list)
        coulomb_array_sum = np.array(coulomb_list_remaining)+np.array(coulomb_list)
        int_E_array_sum = np.array(int_E_list_remaining)+np.array(int_E_list)


        evaluation_array = calc_evaluation_v2(LJ_array_sum,coulomb_array_sum,int_E_array_sum,evaluation_param,evaluation_method,heavy_atom_list)
        evaluation_array = np.array([x if bool_val == False and check == True else 100 for x,bool_val,check in zip(evaluation_array,collision_check,check_list)])


        self.total_charge = new_total_charge
        self.eof_val = all_EOF
        self.ref_frag_list_next_list = ref_frag_list_next_list
        self.c_for_next_list = c_for_next_list
        self.transform_log_main_list = transform_log_main_list
        self.transform_log_side_list = transform_log_side_list
        self.evaluation = evaluation_array
        self.atom = all_atom
        self.coord = coord_array
        self.temp_frag = temp_frag_list
        self.temp_category = temp_category_list
        self.param_for_next = param_for_next_list

class link_result_no_calc_interaction:
    def __init__(self,Core_mol,all_atom,all_coord,all_charge,all_central_atom,split_ind_list_1,split_ind_list_2,all_category,all_param,all_len_list,all_frag,all_ref_frag_list,frag_c_for_next):

        all_transform_log = []

        for i in all_category:

            if i == ['link_c']:
                if len(Core_mol.param_for_next)==2:
                    all_transform_log.append(Core_mol.transform_log_side)

                elif len(Core_mol.param_for_next)==1:
                    all_transform_log.append(Core_mol.transform_log_main)

            elif i == ['link_f']:
                if len(Core_mol.c_for_next)==3:
                    all_transform_log.append(Core_mol.transform_log_side)
                elif len(Core_mol.c_for_next)==2 or len(Core_mol.c_for_next)==1:
                    all_transform_log.append(Core_mol.transform_log_main)

        if all('sym24_XHD_0' not in x and 'sym24_XHD_1' not in x and 'sym24_XHD_2' not in x and 'sym24_XHD_3' not in x and 'sym24_XHD_4' not in x and 'XFO' not in x for x in all_frag):

            check_list = [check_permission_previous(frag_central_atom,Core_mol.central_atom) != False and check_permission_combination(frag_central_atom,Core_mol.central_atom,transform_log) != False for frag_central_atom,transform_log in zip(all_central_atom,all_transform_log)]


        else:
            check_list = [True for _ in all_frag]
        dist_list = [[float(0),float(0),float(get_distance(frag_central_atom,Core_mol.central_atom))] if check_permission_previous(frag_central_atom,Core_mol.central_atom) != False else [np.nan,np.nan,np.nan] for frag_central_atom in all_central_atom]


        dist_array = np.array([row for row, count in zip(dist_list, all_len_list) for _ in range(count)])


        rot_matrix = np.dot(Core_mol.matrix_4th,np.dot(Core_mol.matrix_3rd,np.dot(Core_mol.matrix_2nd,Core_mol.matrix_1st )))

        coord_array = np.dot(rot_matrix, (np.array(all_coord)-dist_array).T).T + Core_mol.central_coord

        all_EOF = []
        transform_log_main_list = []
        transform_log_side_list = []
        c_for_next_list = []
        ref_frag_list_next_list = []
        temp_frag_list = []
        temp_category_list = []
        param_for_next_list = []
        new_total_charge= (np.array(all_charge)+Core_mol.total_charge).tolist()
        for temp_param, temp_category, temp_frag,ref_frag_list,temp_c_for_next in zip(all_param,all_category,all_frag,all_ref_frag_list,frag_c_for_next):

            transform_log_main = copy.deepcopy(Core_mol.transform_log_main)
            transform_log_side = copy.deepcopy(Core_mol.transform_log_side)

            if temp_category ==['link_c']:

                if temp_frag != 'sym24_H1L':

                    ### C→C or C→N
                    if temp_param == []:
                         ### C→C

                        if len(Core_mol.param_for_next)==2:

                            transform_log_side.append(temp_c_for_next[-1])
                        elif len(Core_mol.param_for_next)==1:

                            if len(temp_c_for_next) ==3:
                                transform_log_side.append(temp_c_for_next[-1])
                                transform_log_main.append(temp_c_for_next[-2])
                            elif len(temp_c_for_next) ==2 or len(temp_c_for_next) ==1:
                                transform_log_main.append(temp_c_for_next[-1])

                        ref_frag_list_next = copy.deepcopy(Core_mol.ref_frag_list_next)
                        param_for_next = copy.deepcopy(Core_mol.param_for_next)
                        del param_for_next[-1]

                        del ref_frag_list_next[-1]
                    else:
                        ### C→N
                        param_for_next = temp_param[1:]

                        if len(temp_c_for_next) ==3:
                            transform_log_side.append(temp_param[-1])
                        elif len(temp_c_for_next) ==2 or len(temp_c_for_next) ==1 or len(temp_c_for_next) ==0:
                            if len(param_for_next)==2:
                                transform_log_main.append(temp_param[-2])
                                transform_log_side.append(temp_param[-1])
                            elif len(param_for_next)==1:
                                transform_log_main.append(temp_param[-1])

                        ref_frag_list_next = [param_for_next[-1]]
                    c_for_next = copy.deepcopy(temp_c_for_next)



                if temp_frag == 'sym24_H1L':
                    ref_frag_list_next = copy.deepcopy(Core_mol.ref_frag_list_next)
                    param_for_next = copy.deepcopy(Core_mol.param_for_next)
                    c_for_next = copy.deepcopy(temp_c_for_next)


                    del param_for_next[-1]
                    del ref_frag_list_next[-1]

            elif temp_category ==['link_f']:

                if temp_frag == 'XFO':
                    if len(Core_mol.c_for_next) ==3:
                        transform_log_side.append(temp_param[-1])
                    elif len(Core_mol.c_for_next) ==2 or len(Core_mol.c_for_next) ==1:
                        if len(Core_mol.param_for_next)==2:
                            transform_log_main.append(temp_param[-2])
                            transform_log_side.append(temp_param[-1])
                        elif len(Core_mol.param_for_next)==1:
                            transform_log_main.append(temp_param[-1])
                    ref_frag_list_next=[]
                    param_for_next = temp_param[1:]
                    c_for_next = copy.deepcopy(temp_c_for_next)

                else:

                    if temp_frag == 'sym24_XHD_0' or temp_frag == 'sym24_XHD_1' or temp_frag == 'sym24_XHD_2' or temp_frag == 'sym24_XHD_3' or temp_frag == 'sym24_XHD_4':
                        param_for_next = temp_param[0:]

                    else:
                        param_for_next = temp_param[1:]

                    if temp_frag == 'sym24_XHD_0' or temp_frag == 'sym24_XHD_1' or temp_frag == 'sym24_XHD_2' or temp_frag == 'sym24_XHD_3' or temp_frag == 'sym24_XHD_4':
                        # link hudrazone
                        transform_log_main.append(temp_param[-1])
                        ref_frag_list_next = copy.deepcopy(ref_frag_list)
                        c_for_next = copy.deepcopy(temp_c_for_next)

                        del ref_frag_list_next[0]

                    else:

                        if temp_param[0] =='dummy':

                        #N→C
                            if len(param_for_next)==2:
                                transform_log_side.append(temp_param[-1])
                            elif len(param_for_next)==1:
                                if len(temp_param) ==3:
                                    transform_log_side.append(temp_param[-1])
                                    transform_log_main.append(temp_param[-2])
                                elif len(temp_param) ==2 or len(temp_param) ==1:
                                    transform_log_main.append(temp_param[-1])
                            ref_frag_list_next = copy.deepcopy(ref_frag_list)

                            del param_for_next[-1]
                            del ref_frag_list_next[-1]
                            c_for_next = copy.deepcopy(temp_c_for_next)

                        else:
                        #N→N

                            if len(Core_mol.c_for_next) ==3:
                                transform_log_side.append(temp_param[-1])
                            elif len(Core_mol.c_for_next) ==2 or len(Core_mol.c_for_next) ==1:
                                if len(param_for_next)==2:
                                    transform_log_main.append(temp_param[-2])
                                    transform_log_side.append(temp_param[-1])
                                elif len(param_for_next)==1:
                                    transform_log_main.append(temp_param[-1])
                            ref_frag_list_next = copy.deepcopy(ref_frag_list)
                            del ref_frag_list_next[0]
                            c_for_next = copy.deepcopy(temp_c_for_next)

            if temp_frag == 'sym24_H1L' or temp_frag == 'carboxylate':
                all_EOF.append(True)
            else:
                all_EOF.append(False)

            transform_log_main_list.append(transform_log_main)
            transform_log_side_list.append(transform_log_side)
            c_for_next_list.append(c_for_next)
            ref_frag_list_next_list.append(ref_frag_list_next)
            temp_frag_list.append(temp_frag)
            temp_category_list.append(temp_category)
            param_for_next_list.append(param_for_next)

        self.total_charge = new_total_charge
        self.eof_val = all_EOF
        self.ref_frag_list_next_list = ref_frag_list_next_list
        self.c_for_next_list = c_for_next_list
        self.transform_log_main_list = transform_log_main_list
        self.transform_log_side_list = transform_log_side_list

        self.atom = all_atom
        self.coord = coord_array
        self.temp_frag = temp_frag_list
        self.temp_category = temp_category_list
        self.param_for_next = param_for_next_list

##function for LigX

def get_distance(atomname_add,atomname_previous):
    pattern = re.compile(r'^(X[^_]*_?).*$')
    atomname_add_modi = pattern.sub(r'\1', atomname_add)
    atomname_previous_modi = pattern.sub(r'\1', atomname_previous)
    atom_names = atomname_add_modi+','+atomname_previous_modi
    distance = distance_dict[atom_names]
    return distance

def check_permission_previous(atomname_add,atomname_previous):
    pattern = re.compile(r'^(X[^_]*_?).*$')
    atomname_add_modi = pattern.sub(r'\1', atomname_add)
    atomname_previous_modi = pattern.sub(r'\1', atomname_previous)
    atom_names = atomname_add_modi+','+atomname_previous_modi
    if atom_names in distance_dict:
        return True
    else:
        return False

def check_permission_combination(atomname_add,atomname_previous,transform_log):
    atom_minus2 = transform_log[-2]
    pattern = re.compile(r'^(X[^_]*_?).*$')
    atomname_add_modi = pattern.sub(r'\1', atomname_add)
    atomname_previous_modi = pattern.sub(r'\1', atomname_previous)
    atom_minus2_modi = pattern.sub(r'\1', atom_minus2)
    atom_names =atomname_add_modi+','+atomname_previous_modi+','+atom_minus2_modi
    if atom_names in not_permitted_comb_dict:
        return False
    else:
        return True


def get_desired_list(desired_list):
    desired_frag_list = []
    for i in desired_list:
        temp =''
        j = i.split('/')
        for k in range(len(j)-1):
            temp +=j[k]+'/'
        frag_name = j[-1]
        charge_file = temp+'charge.txt'
        param_file = temp+'param.txt'
        charge = extract_charge(charge_file)
        angle_value = 0
        sym_num = extract_number(frag_name)
        while angle_value < 360//sym_num:
            add_list =[i,param_file,angle_value,charge]
            angle_value += 30
            desired_frag_list.append(add_list)
    return desired_frag_list

def get_desired_list_rollout(desired_list):
    desired_frag_list = []
    for i in desired_list:
        temp =''
        j = i.split('/')
        for k in range(len(j)-1):
            temp +=j[k]+'/'
        frag_name = j[-1]
        charge_file = temp+'charge.txt'
        param_file = temp+'param.txt'
        charge = extract_charge(charge_file)
        angle_value = 0
        sym_num = extract_number(frag_name)
        while angle_value < 360//sym_num:
            add_list =[i,param_file,angle_value,charge]
            angle_value += 30
            desired_frag_list.append(add_list)
        return desired_frag_list

def import_structure_2(dir_core,core_file):
    core_name = str(core_file).replace('.xyz', '')
    path_core = dir_core+'/'+str(core_file)
    with open(path_core) as core_init:
        reader = csv.reader(core_init)
        list_row = core_init.read().split('\n')

    i_arr_1 =[]
    for i in list_row:
        i_list =  str(i).split()
        i_arr_1.append(i_list)

    first_column = [item[0] for item in i_arr_1 if len(item) > 0]


    S = None
    SE = None
    P = None
    PE = None

    if any('NQ' in element for element in first_column):
        Q = 'NQ'
        QE ='N'
    elif any('OQ' in element for element in first_column):
        Q = 'OQ'
        QE ='O'
    elif any('HQ' in element for element in first_column):
        Q = 'HQ'
        QE ='H'
    elif any('DQ' in element for element in first_column):
        Q = 'DQ'
        QE ='D'
    elif any('SQ' in element for element in first_column):
        Q = 'SQ'
        QE ='S'
    elif any('CQ' in element for element in first_column):
        Q = 'CQ'
        QE ='C'
    else:
        print('no Q')
        sys.exit()


    core_atoms =[]
    fin_index = len(i_arr_1)-1
    for sep in range(2,fin_index):
        core_atoms.append(i_arr_1[sep])
    return core_atoms, S, Q, P


def move_origin_then_rotation_twice(list1, origin, rotation_point_1,rotation_point_2):
    origin_coords = np.array([point[1:4] for point in list1 if point[0] == origin], dtype=float)[0]
    translated_points = []
    for point in list1:
        name, x, y, z = point
        translated_coords = np.array([float(x), float(y), float(z)], dtype=float) - origin_coords
        translated_points.append([name, *translated_coords])

    coords_after_translation_1 = np.array([point[1:4] for point in translated_points if point[0] == rotation_point_1], dtype=float)[0]
    rho_1 = np.linalg.norm(coords_after_translation_1)
    phi_1 = math.degrees(math.atan2(coords_after_translation_1[1], coords_after_translation_1[0]))
    theta_1 = math.degrees(math.acos(coords_after_translation_1[2] / rho_1))
    phi_xy_side_2nd = np.radians(-1*phi_1)
    rotation_matrix_xy_side_2nd = np.array([[np.cos(phi_xy_side_2nd), -np.sin(phi_xy_side_2nd), 0],
                                  [np.sin(phi_xy_side_2nd), np.cos(phi_xy_side_2nd), 0],
                                  [0, 0, 1]])
    translated_points_2 = []
    for point in translated_points:
        names, xs, ys, zs = point
        original_coords_1 = np.array([xs, ys, zs])
        rotated_coords_1 = np.dot(rotation_matrix_xy_side_2nd, original_coords_1)
        translated_points_2.append([names, *rotated_coords_1])
    theta_xz_side_2nd = np.radians(-1*theta_1)
    rotation_matrix_xz_side_2nd = np.array([[np.cos(theta_xz_side_2nd), 0, np.sin(theta_xz_side_2nd)],
                                  [0, 1, 0],
                                  [-np.sin(theta_xz_side_2nd), 0, np.cos(theta_xz_side_2nd)]])
    translated_points_3 = []
    for point in translated_points_2:
        names, xs, ys, zs = point
        original_coords_2 = np.array([xs, ys, zs])
        rotated_coords_2 = np.dot(rotation_matrix_xz_side_2nd, original_coords_2)
        translated_points_3.append([names, *rotated_coords_2])

    coords_after_translation_2 = np.array([point[1:4] for point in translated_points_3 if point[0] == rotation_point_2], dtype=float)[0]
    rho_2 = np.linalg.norm(coords_after_translation_2)
    phi_2 = math.degrees(math.atan2(coords_after_translation_2[1], coords_after_translation_2[0]))
    theta_2 = math.degrees(math.acos(coords_after_translation_2[2] / rho_2))

    phi_xy_side_2nd = np.radians(-1*phi_2)
    rotation_matrix_xy_side_2nd = np.array([[np.cos(phi_xy_side_2nd), -np.sin(phi_xy_side_2nd), 0],
                                  [np.sin(phi_xy_side_2nd), np.cos(phi_xy_side_2nd), 0],
                                  [0, 0, 1]])

    translated_points_4 = []
    for point in translated_points_3:
        names, xs, ys, zs = point
        original_coords_3 = np.array([xs, ys, zs])
        rotated_coords_3 = np.dot(rotation_matrix_xy_side_2nd, original_coords_3)
        translated_points_4.append([names, *rotated_coords_3])

    return translated_points_4,[phi_1,theta_1,phi_2,theta_2],origin_coords

def read_param(file_path):
    filename = str(file_path).replace('.txt', '')
    with open(file_path) as file_init:
        reader = csv.reader(file_init)
        list_row = file_init.read().split('\n')
    txt_list =[]
    link_category = []
    for i in list_row:
        i_list =  str(i).split(' ')
        if i_list[0] == 'link_c':
            link_category.append(i_list[0])
        if i_list[0] == 'link_f':
            link_category.append(i_list[0])
        if i_list[0] == 'dihedral':
            link_category.append(i_list[0])
        if i_list[0] == 'one_atom':
            link_category.append(i_list[0])
        if i_list[0] == 'cyano':
            link_category.append(i_list[0])
        if i_list[0] == 'link_atom_1':
            txt_list.append(i_list[1])
        elif i_list[0] == 'link_atom_2':
            txt_list.append(i_list[1])
        elif i_list[0] == 'link_atom_3':
            txt_list.append(i_list[1])
    return txt_list,link_category

def ref_frag_from_param(file_path):
    with open(file_path) as file_init:
        reader = csv.reader(file_init)
        list_row = file_init.read().split('\n')
    ref_frag_list = []
    for i in list_row:
        i_list =  str(i).split(' ')
        if i_list[0] == 'ref_atom':
            ref_frag_list.append(i_list[1])
        elif i_list[0] == 'ref_atom_2':
            ref_frag_list.append(i_list[1])
        elif i_list[0] == 'ref_atom_3':
            ref_frag_list.append(i_list[1])
    if len(ref_frag_list) !=0:
        return ref_frag_list
    else:
        return None

def read_xyz_data(file_path):
    filename_1 = str(file_path).replace('.xyz', '')
    filename_2 = filename_1.split('/')
    filename = '/'.join(filename_2[-1:])
    with open(file_path) as file_init:
        reader = csv.reader(file_init)
        list_row = file_init.read().split('\n')
    txt_list =[]
    for i in list_row:
        i_list =  str(i).split()
        txt_list.append(i_list)
    txt_list_removed =[]
    fin_index = len(txt_list)-1
    for j in range(2,fin_index):
        txt_list_removed.append(txt_list[j])
    return txt_list_removed,filename

def get_transform_dict_hyd(fragment_set):
    all_unit = ['F1L','F2L','F3L','F2LB','C1LB','C1L','C2L','C3L','XHD','XFO']
    transform_dict = {}
    for transformations in all_unit:
        transform_list = []
        frag_dir = fragment_set+'/'+transformations
        try:
            os.listdir(frag_dir)
        except:
            continue
        else:
            frag_list= os.listdir(frag_dir)
        frag_list.sort()
        for each_frag in frag_list:
            each_frag_dir = frag_dir+'/'+each_frag
            frag_list_ini = os.listdir(each_frag_dir)
            param_file = each_frag_dir+'/param.txt'
            charge_file = each_frag_dir+'/charge.txt'
            charge = extract_charge(charge_file)
            frag_list_init_2 = [element for element in frag_list_ini if element != 'param.txt']
            frag_list = [element for element in frag_list_init_2 if element != 'charge.txt']
            for frag_file in frag_list:
                path_frag = each_frag_dir+'/'+str(frag_file)
                angle_value = 0
                sym_num = extract_number(frag_file)
                if transformations != 'XHD':
                    while angle_value < 360//sym_num:
                        add_list =[path_frag,param_file,angle_value,charge]
                        angle_value += 30
                        transform_list.append(add_list)
                elif transformations == 'XHD':
                        add_list =[path_frag,param_file,angle_value,charge]
                        angle_value += 180
                        transform_list.append(add_list)
        transform_dict[transformations]=transform_list
    return transform_dict

def get_transform_dict_hyd_rollout(fragment_set):
    all_unit = ['F1L','F2L','F3L','F2LB','C1LB','C1L','C2L','C3L','XHD','XFO']
    transform_dict = {}
    for transformations in all_unit:
        transform_list = []
        frag_dir = fragment_set+'/'+transformations
        try:
            os.listdir(frag_dir)
        except:
            continue
        frag_list= os.listdir(frag_dir)
        frag_list.sort()
        for each_frag in frag_list:
            each_frag_dir = frag_dir+'/'+each_frag
            frag_list_ini = os.listdir(each_frag_dir)
            param_file = each_frag_dir+'/param.txt'
            charge_file = each_frag_dir+'/charge.txt'
            charge = extract_charge(charge_file)
            frag_list_init_2 = [element for element in frag_list_ini if element != 'param.txt']
            frag_list = [element for element in frag_list_init_2 if element != 'charge.txt']
            for frag_file in frag_list:
                path_frag = each_frag_dir+'/'+str(frag_file)
                angle_value = 0
                sym_num = extract_number(frag_file)
                while angle_value < 360//sym_num:
                    add_list =[path_frag,param_file,angle_value,charge]
                    angle_value += 30
                    transform_list.append(add_list)
        transform_dict[transformations]=transform_list
    return transform_dict

def extract_number(frag_file):
    match = re.match(r'^sym(\d{1,2})(_.*)?$', frag_file)
    if match:
        number = int(match.group(1))
        if 1 <= number <= 24:
            return number
    return 1

def extract_charge(charge_file):
    (charge_file)
    with open(charge_file) as file_init:
        reader = csv.reader(file_init)
        list_row = file_init.read().split('\n')
    return list_row[0]

def save_csv(cord_list,output_name,result_dir):
    with open('output.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        for row in cord_list:
            writer.writerow(row)
    with open('output.csv') as reader:
        content = reader.read()
    content = content.replace(',', ' ')
    output_file_name =result_dir +'/'+ output_name+'.xyz'
    with open(output_file_name, 'w') as writer:
        writer.write(content)
    os.remove('output.csv')

def save_csv_2(cord_list,output_name,result_dir):
    with open('output.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        for row in cord_list:
            writer.writerow(row)
    with open('output.csv') as reader:
        content = reader.read()
    content = content.replace(',', ' ')
    output_file_name =result_dir+'/'+ output_name+'.csv'
    with open(output_file_name, 'w') as writer:
        writer.write(content)
    os.remove('output.csv')

def save_csv_3(cord_list,output_name,result_dir):
    with open('output.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        for row in cord_list:
            writer.writerow(row)
    with open('output.csv') as reader:
        content = reader.read()
    output_file_name =result_dir+'/'+ output_name+'.csv'
    with open(output_file_name, 'a') as writer:
        writer.write(content)
    os.remove('output.csv')

def save_csv_debug(cord_list,output_name):
    with open('output.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        for row in cord_list:
            writer.writerow(row)
    with open('output.csv') as reader:
        content = reader.read()
    content = content.replace(',', ' ')
    output_file_name ='debug/'+ output_name+'.xyz'
    with open(output_file_name, 'w') as writer:
        writer.write(content)
    os.remove('output.csv')

def save_charge(total_charge,dir_charge):
    with open('output.csv', 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(str(total_charge))
    with open('output.csv') as reader:
        content = reader.read()
    content = content.replace(',', '')
    output_file_name =dir_charge+'/total_charge.txt'
    with open(output_file_name, 'w') as writer:
        writer.write(content)
    os.remove('output.csv')

def get_smiles(cord_list,total_charge,smiles_filter):
    if smiles_filter == True:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        from rdkit.Chem import Draw
        from rdkit.Chem import AllChem
        with open('temp.csv', 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            for row in cord_list:
                writer.writerow(row)
        with open('temp.csv') as reader:
            content = reader.read()
        content = content.replace(',', ' ')
        output_file_name ='temp.xyz'
        with open(output_file_name, 'w') as writer:
            writer.write(content)
        os.remove('temp.csv')
        mol = Chem.MolFromXYZFile('temp.xyz')
        try:
            rdDetermineBonds.DetermineBonds(mol,charge=total_charge)
        except:
            return ('no_SMILES')
        else:
            #rdDetermineBonds.DetermineBonds(mol,charge=total_charge)
            try:
                smiles = Chem.MolToSmiles(mol)
            except ValueError:
                return ('no_SMILES')
            else:
                #smiles = Chem.MolToSmiles(mol)
                os.remove('temp.xyz')
                return smiles
    else:
        return 'without_SMILES_conversion'

def atomname_change_after_link_XHD(atom_list):
    atom_list = ['C' if x == 'X4' else x for x in atom_list]
    atom_list = ['H' if x == 'B4' else x for x in atom_list]
    atom_list = ['N' if x == 'W1' else x for x in atom_list]
    return atom_list

def atomname_change_after_link_f(atom_list,original_name_x,original_name_r,n_at_front):
    atom_list = ['CX_0' if x == original_name_x else x for x in atom_list]
    for r in ['H','C','O','N','S','Cl','Br','I','P','F','D']:
        R_name = 'R'+r+'_'+str(n_at_front)
        if R_name == original_name_r:
            atom_list = [r if x == original_name_r else x for x in atom_list]
    return atom_list

def remove_VC_after_link_f(atom_list,cord_list,VC_name_at_front,V_name_at_front):
    new_atom_list = [x for x in atom_list if x != VC_name_at_front]
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != VC_name_at_front]
    new_atom_list_2 = [x for x in new_atom_list if x != V_name_at_front]
    new_cord_list_2 = [y for x,y in zip(new_atom_list,new_cord_list) if x != V_name_at_front]
    return new_atom_list_2,new_cord_list_2

def remove_atoms_after_link_XHD(atom_list,cord_list):
    new_atom_list = [x for x in atom_list if x != 'A4']
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != 'A4']
    new_atom_list_2 = [x for x in new_atom_list if x != 'Y1']
    new_cord_list_2 = [y for x,y in zip(new_atom_list,new_cord_list) if x != 'Y1']
    new_atom_list_3 = [x for x in new_atom_list_2 if x != 'Z1']
    new_cord_list_3 = [y for x,y in zip(new_atom_list_2,new_cord_list_2) if x != 'Z1']
    return new_atom_list_3,new_cord_list_3

def atomname_change_after_link_c(atom_list,original_name_x,original_name_r):
    if original_name_x != 'XC3g_1':
        atom_list = ['C' if x == original_name_x else x for x in atom_list]
    elif original_name_x == 'XC3g_1':
        atom_list = ['XC3g_2' if x == original_name_x else x for x in atom_list]
    atom_list = ['C' if x == 'CX_0' else x for x in atom_list]
    atom_list = ['N' if x == 'NX_0' else x for x in atom_list]
    for r in ['H','C','O','N','S','Cl','Br','I','P','F','D']:
        R_name = 'R'+r+'_1'
        if R_name == original_name_r:
            atom_list = [r if x == original_name_r else x for x in atom_list]
    return atom_list

def atomname_change_after_link_c_f(atom_list,original_name_x,original_name_r):
    if original_name_x != 'XC3g_1':
        atom_list = ['C' if x == original_name_x else x for x in atom_list]
    elif original_name_x == 'XC3g_1':
        atom_list = ['XC3g_2' if x == original_name_x else x for x in atom_list]
    for r in ['H','C','O','N','S','Cl','Br','I','P','F','D']:
        R_name = 'R'+r+'_1'
        if R_name == original_name_r:
            atom_list = [r if x == original_name_r else x for x in atom_list]
    return atom_list

def atomname_change_after_link_f_to_c(atom_list):
    atom_list = ['C' if x == 'CX_0' else x for x in atom_list]
    atom_list = ['XC3g_2' if x == 'XC3g_1' else x for x in atom_list]
    return atom_list

def remove_VC_after_link_c(atom_list,cord_list,VC_name_at_front,V_name_at_front):
    new_atom_list = [x for x in atom_list if x != VC_name_at_front]
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != VC_name_at_front]
    new_atom_list_2 = [x for x in new_atom_list if x != V_name_at_front]
    new_cord_list_2 = [y for x,y in zip(new_atom_list,new_cord_list) if x != V_name_at_front]
    return new_atom_list_2,new_cord_list_2

def remove_V_only_after_link_c(atom_list,cord_list,V_name_at_front):
    new_atom_list = [x for x in atom_list if x != V_name_at_front]
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != V_name_at_front]
    return new_atom_list,new_cord_list

def remove_VC_only_after_link_c(atom_list,cord_list,VC_name_at_front):
    new_atom_list = [x for x in atom_list if x != VC_name_at_front]
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != VC_name_at_front]
    return new_atom_list,new_cord_list

def atomname_change_after_link_h(atom_list):
    atom_list = ['C' if x == 'CX_0' else x for x in atom_list]
    return atom_list

def remove_VC_after_link_h(atom_list,cord_list,V_name_at_front):
    new_atom_list = [x for x in atom_list if x != V_name_at_front]
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != V_name_at_front]
    return new_atom_list,new_cord_list

def remove_VC_after_link_h_f(atom_list,cord_list,VC_name_at_front):
    new_atom_list = [x for x in atom_list if x != VC_name_at_front]
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != VC_name_at_front]
    return new_atom_list,new_cord_list

def remove_V_name_at_front_of_frag_after_link_f_c(atom_list,cord_list,V_name_at_front):
    new_atom_list = [x for x in atom_list if x != V_name_at_front]
    new_cord_list = [y for x,y in zip(atom_list,cord_list) if x != V_name_at_front]
    return new_atom_list,new_cord_list

def atomname_change_after_link_XFO(atom_list,original_name_r):
    atom_list = ['C' if x == 'XC5' else x for x in atom_list]
    for r in ['H','C','O','N','S','Cl','Br','I','P','F','D']:
        R_name = 'R'+r+'_1'
        if R_name == original_name_r:
            atom_list = [r if x == original_name_r else x for x in atom_list]
    return atom_list

def get_initial_x(initial_remaining_data):
    first_column= [item[0] for item in initial_remaining_data if len(item) > 0]
    for x in ['XC3_2','XC3g_2','XC2_2','XC1_2','XCA_2']:
        if any(x in element for element in first_column):
            transform_log_main = ['',x]
            transform_log_side = ['',x]
            c_for_next = [x]
            param_for_next = []
            ref_frag_list_next = []
            return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('X4' in element for element in first_column):
        transform_log_main = ['','XCA_1']
        transform_log_side = ['','XCA_1']
        c_for_next = ['XCA_1']
        param_for_next = []
        ref_frag_list_next = []
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('XC5' in element for element in first_column):
        transform_log_main = ['','XCA_1']
        transform_log_side = ['','XCA_1']
        c_for_next = ['XCA_1']
        param_for_next = []
        ref_frag_list_next = []
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('YCA' in element for element in first_column):
        transform_log_main = ['XCA_1','YCA']
        transform_log_side = ['XCA_1','YCA']
        c_for_next = []
        param_for_next = ['YCA']
        if any('CX_0' in element for element in first_column):
            ref_frag_list_next = ['CX_0']
        elif any('NX_0' in element for element in first_column):
            ref_frag_list_next = ['NX_0']
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('YO' in element for element in first_column):
        transform_log_main = ['XC3_1','YO']
        transform_log_side = ['XC3_1','YO']
        c_for_next = []
        param_for_next = ['YO']
        if any('CX_0' in element for element in first_column):
            ref_frag_list_next = ['CX_0']
        elif any('NX_0' in element for element in first_column):
            ref_frag_list_next = ['NX_0']
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('YS' in element for element in first_column):
        transform_log_main = ['XC3_1','YS']
        transform_log_side = ['XC3_1','YS']
        c_for_next = []
        param_for_next = ['YS']
        ref_frag_list_next = ['CX_0']
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('YSA' in element for element in first_column):
        transform_log_main = ['XC3_1','YSA']
        transform_log_side = ['XC3_1','YSA']
        c_for_next = []
        param_for_next = ['YSA']
        ref_frag_list_next = ['YN3']
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('YNA' in element for element in first_column):
        transform_log_main = ['XCA_1','YNA']
        transform_log_side = ['XCA_1','YNA']
        c_for_next = []
        param_for_next = ['YNA']
        ref_frag_list_next = ['CX_0']
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('YN2' in element for element in first_column):
        transform_log_main = ['XCA_1','YN2']
        transform_log_side = ['XCA_1','YN2']
        c_for_next = []
        param_for_next = ['YN2']
        if any('CX_0' in element for element in first_column):
            ref_frag_list_next = ['CX_0']
        elif any('NX_0' in element for element in first_column):
            ref_frag_list_next = ['NX_0']
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next
    if any('YN3' in element for element in first_column):
        transform_log_main = ['XC3_1','YN3']
        transform_log_side = ['XC3_1','YN3']
        c_for_next = []
        param_for_next = ['YN3']
        if any('CX_0' in element for element in first_column):
            ref_frag_list_next = ['CX_0']
        elif any('NX_0' in element for element in first_column):
            ref_frag_list_next = ['NX_0']
        return transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next

def cut_core(combined_cord_2,Q):
    frag_only = []
    first_column= [item[0] for item in combined_cord_2 if len(item) > 0]
    for i in range(len(combined_cord_2)):
        if first_column[i] == Q:
            qq = i+1
            break
    for i in range(qq,len(combined_cord_2)):
        frag_only.append(combined_cord_2[i])
    return frag_only

def cut_frag(combined_cord_2,Q):
    core_only = []
    first_column= [item[0] for item in combined_cord_2 if len(item) > 0]
    for i in range(len(combined_cord_2)):
        if first_column[i] == Q:
            qq = i+1
    for i in range(qq):
        core_only.append(combined_cord_2[i])
    return core_only

def atomname_change_before_link_n(atom_list):
    atom_list = ['C' if x == 'YCA' else x for x in atom_list]
    atom_list = ['C' if x == 'YCO' else x for x in atom_list]
    atom_list = ['O' if x == 'YO' else x for x in atom_list]
    atom_list = ['S' if x == 'YS' else x for x in atom_list]
    atom_list = ['S' if x == 'YSA' else x for x in atom_list]
    atom_list = ['N' if x == 'YNA' else x for x in atom_list]
    atom_list = ['N' if x == 'YN3' else x for x in atom_list]
    atom_list = ['N' if x == 'YN2' else x for x in atom_list]
    return atom_list

def atomname_change_before_link_n_c(atom_list):
    atom_list = ['C' if x == 'YCA' else x for x in atom_list]
    atom_list = ['C' if x == 'YCO' else x for x in atom_list]
    atom_list = ['C' if x == 'CX_0' else x for x in atom_list]
    atom_list = ['O' if x == 'YO' else x for x in atom_list]
    atom_list = ['S' if x == 'YS' else x for x in atom_list]
    atom_list = ['S' if x == 'YSA' else x for x in atom_list]
    atom_list = ['N' if x == 'YNA' else x for x in atom_list]
    atom_list = ['N' if x == 'YN3' else x for x in atom_list]
    atom_list = ['N' if x == 'YN2' else x for x in atom_list]
    return atom_list

def atomname_change_final(cord_list):
    for k in range(len(cord_list)):
        if cord_list[k][0] == 'YCA':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'YCO':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'YO':
            cord_list[k][0] = 'O'
        if cord_list[k][0] == 'YS':
            cord_list[k][0] = 'S'
        if cord_list[k][0] == 'YSA':
            cord_list[k][0] = 'S'
        if cord_list[k][0] == 'YNA':
            cord_list[k][0] = 'N'
        if cord_list[k][0] == 'YN3':
            cord_list[k][0] = 'N'
        if cord_list[k][0] == 'YN2':
            cord_list[k][0] = 'N'
        if cord_list[k][0] == 'NX_0':
            cord_list[k][0] = 'N'
        if cord_list[k][0] == 'CX_0':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC3g_1':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC3g_2':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC3g_3':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC3_1':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC3_2':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC3_3':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC2_1':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC2_2':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC2_3':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC1_1':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC1_2':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XC1_3':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XCA_1':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XCA_2':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'XCA_3':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'NS':
            cord_list[k][0] = 'N'
        if cord_list[k][0] == 'OS':
            cord_list[k][0] = 'O'
        if cord_list[k][0] == 'HS':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'DS':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'SS':
            cord_list[k][0] = 'S'
        if cord_list[k][0] == 'CS':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'NQ':
            cord_list[k][0] = 'N'
        if cord_list[k][0] == 'OQ':
            cord_list[k][0] = 'O'
        if cord_list[k][0] == 'HQ':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'DQ':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'SQ':
            cord_list[k][0] = 'S'
        if cord_list[k][0] == 'CQ':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'NP':
            cord_list[k][0] = 'N'
        if cord_list[k][0] == 'OP':
            cord_list[k][0] = 'O'
        if cord_list[k][0] == 'HP':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'DP':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'SP':
            cord_list[k][0] = 'S'
        if cord_list[k][0] == 'CP':
            cord_list[k][0] = 'C'
        if cord_list[k][0] == 'D':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'DH':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'HA':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'HR':
            cord_list[k][0] = 'H'
        if cord_list[k][0] == 'HK':
            cord_list[k][0] = 'H'
        for r in ['H','C','O','N','S','Cl','Br','I','P','D','F']:
            R_name = 'R'+r+'_1'
            if cord_list[k][0] == R_name:
                cord_list[k][0] = r
        for r in ['H','C','O','N','S','Cl','Br','I','P','D','F']:
            R_name = 'R'+r+'_2'
            if cord_list[k][0] == R_name:
                cord_list[k][0] = r
    return cord_list

def get_reconstract_strategy(num_protein_structures,all_score_strategy):
    index_list =[]
    id_index = all_score_strategy[0].index('ID')
    index_list.append(['ID',id_index])
    for i in range(num_protein_structures):
        x = 'LJ['+str(i)+']'
        index_num = all_score_strategy[0].index(x)
        index_list.append([x,index_num])
    LJ_index = all_score_strategy[0].index('LJ_average')
    index_list.append(['LJ_average',LJ_index])
    for i in range(num_protein_structures):
        x = 'coulomb['+str(i)+']'
        index_num = all_score_strategy[0].index(x)
        index_list.append([x,index_num])
    coulomb_index = all_score_strategy[0].index('coulomb_average')
    index_list.append(['coulomb_average',coulomb_index])
    for i in range(num_protein_structures):
        x = 'interact_E['+str(i)+']'
        index_num = all_score_strategy[0].index(x)
        index_list.append([x,index_num])
    interact_E_index = all_score_strategy[0].index('interact_E_average')
    index_list.append(['interact_E_average',interact_E_index])
    for i in range(num_protein_structures):
        x = 'LE_LJ['+str(i)+']'
        index_num = all_score_strategy[0].index(x)
        index_list.append([x,index_num])
    LE_LJ_index = all_score_strategy[0].index('LE_LJ_average')
    index_list.append(['LE_LJ_average',LE_LJ_index])
    for i in range(num_protein_structures):
        x = 'LE_interact_E['+str(i)+']'
        index_num = all_score_strategy[0].index(x)
        index_list.append([x,index_num])
    LE_int_E_index = all_score_strategy[0].index('LE_interact_E_average')
    index_list.append(['LE_interact_E_average',LE_int_E_index])
    for i in range(num_protein_structures):
        x = 'LE_LJ+coulomb['+str(i)+']'
        index_num = all_score_strategy[0].index(x)
        index_list.append([x,index_num])
    LE_int_E_cl_index = all_score_strategy[0].index('LE_LJ+coulomb_average')
    index_list.append(['LE_LJ+coulomb_average',LE_int_E_cl_index])
    strategy_index = all_score_strategy[0].index('strategy')
    index_list.append(['strategy',strategy_index])
    return index_list

def extract_des_category(des):
    j = des.split('/')
    return j[-2]

def create_split_ind_list(len_list):
    split_ind_list_1 = [0]
    split_ind_list_2 = []
    sum_count = 0
    for i in len_list:
        sum_count += i
        split_ind_list_1.append(sum_count)
        split_ind_list_2.append(sum_count)
    split_ind_list_1.pop()
    return split_ind_list_1,split_ind_list_2

def create_split_ind_list_for_interaction(split_ind_list_2,len_core):
    split_ind_list_1st = [0]
    split_ind_list_2nd = []
    sum_count = 0
    for i in split_ind_list_2:
        split_ind_list_1st.append(i*len_core)
        split_ind_list_2nd.append(i*len_core)
    split_ind_list_1st.pop()
    return split_ind_list_1st,split_ind_list_2nd

def calc_interaction_remaining(all_atom,coord_array, all_protein_cord):
    coord_array = np.array(coord_array)

    LJ_list = []
    coulomb_list = []
    int_E_list = []

    if all_atom != []:

        cut_off_ratio =3
        cut_off_dist = 10
        electronic_ratio = 0.8

        LJ_param = {'H':[2.38760856462,0.1882800],'D':[0.40001352445,0.1924640],'C':[3.58141284692,0.2343040],'N':[3.52795892384,0.2928800],'O':[3.02905564168,0.5020800],'F':[3.02905564168,0.5020800],'P':[3.83086448800,2.4476400],'S':[3.56359487256,1.8828000],'Cl':[3.31414323148,0.962320],'Br':[3.52795892384,1.3388800],'I':[3.99122625727,2.1756800],'Mg':[2.11142996199,0.0627600],'Zn':[1.94215920555,1.0460000],'Vn':[0,0]}

        coulomb_param = {'H':0.09,'D':0.43,'C':0.07,'N':-0.4,'O':-0.5,'F':-0.19,'P':1.0,'S':0.12,'Cl':-0.17,'Br':-0.14,'I':-0.08,'Mg':2,'Zn':2,'Vn':0}

        water_epsilon = 78.5

        replace_dict = {'CQ': 'C', 'X1': 'C', 'X2': 'C', 'X3': 'C', 'X3': 'C','X4':'N', 'W1': 'N', 'A1': 'H', 'A2': 'H', 'A3': 'H', 'A3': 'H', 'B1': 'H', 'B2': 'H', 'B3': 'H', 'B4': 'H', 'NS': 'N', 'OS': 'O', 'SS': 'S', 'HS': 'H', 'CS': 'C', 'NQ': 'N', 'OQ': 'O', 'SQ': 'S', 'HQ': 'H', 'CQ': 'C', 'NP': 'N', 'OP': 'O', 'SP': 'S', 'HP': 'H', 'CP': 'C', 'VLA': 'Vn', 'VLB': 'Vn', 'VLC': 'Vn', 'VC_1': 'Vn', 'VC_2': 'Vn', 'VC_3': 'Vn', 'VC_4': 'Vn', 'Y1': 'Vn', 'Z1': 'Vn', 'RH_1': 'H', 'RH_2': 'H', 'RH_3': 'H', 'RH_4': 'H', 'RC_1': 'C', 'RC_2': 'C', 'RC_3': 'C', 'RC_4': 'C', 'RO_1': 'O', 'RO_2': 'O', 'RO_3': 'O', 'RO_4': 'O', 'RN_1': 'N', 'RN_2': 'N', 'RN_3': 'N', 'RN_4': 'N','RS_1': 'S', 'RS_2': 'S', 'RS_3': 'S', 'RS_4': 'S','YCA': 'C', 'YCO': 'C', 'YNA': 'N','YN2': 'N','YN3': 'N', 'YS': 'S', 'YSA': 'S', 'YO': 'O','XC3_1': 'C','XC3_2': 'C','XC3_3': 'C','XC3_4': 'C', 'XC3g_1': 'C','XC3g_2': 'C','XC3g_3': 'C','XC3g_4': 'C', 'XC2_1': 'C','XC2_2': 'C','XC2_3': 'C','XC2_4': 'C', 'XC1_1': 'C','XC1_2': 'C','XC1_3': 'C','XC1_4': 'C', 'XCA_1': 'C','XCA_2': 'C','XCA_3': 'C','XCA_4': 'C','CX_0': 'C','A4': 'N', 'B4': 'H','NX_0': 'N', 'DH': 'D', 'HA': 'H', 'HR': 'H', 'HK': 'H', 'CAD': 'C', 'CRD': 'C', 'CKD': 'C','HET': 'N', 'DH_1': 'D', 'HA_1': 'H', 'HR_1': 'H', 'HK_1': 'H', 'DH_2': 'D', 'HA_2': 'H', 'HR_2': 'H', 'HK_2': 'H', 'DH_3': 'D', 'HA_3': 'H', 'HR_3': 'H', 'HK_3': 'H', 'DH_4': 'D', 'HA_4': 'H', 'HR_4': 'H', 'HK_4': 'H', 'DH_5': 'D', 'HA_5': 'H', 'HR_5': 'H', 'HK_5': 'H', 'DH_6': 'D', 'HA_6': 'H', 'HR_6': 'H', 'HK_6': 'H', 'DH_7': 'D', 'HA_7': 'H', 'HR_7': 'H', 'HK_7': 'H', 'DH_8': 'D', 'HA_8': 'H', 'HR_8': 'H', 'HK_8': 'H', 'DH_9': 'D', 'HA_9': 'H', 'HR_9': 'H', 'HK_9': 'H', 'V_core': 'Vn', 'RF_1': 'F', 'RD_1': 'D'}


        for core in all_protein_cord:
            LJ_u = 0
            coulomb_u = 0




            names2, coords2 = np.array(core)[:, 0], np.array(core)[:, 1:].astype(float)
            distances = np.linalg.norm(coord_array[:, np.newaxis, :] - coords2, axis=-1)
            combinations = np.array(np.meshgrid(all_atom, names2)).T.reshape(-1, 2)

            for key, value in replace_dict.items():
                mask = combinations == key
                combinations[mask] = value
            true_elem_distance = np.column_stack((combinations, distances.flatten()))
            i,j,dist = zip(*true_elem_distance)

            sigma_list = np.array([[((LJ_param[i][0]+LJ_param[j][0])/2)/float(dist)]for i, j, dist in zip(i, j, dist) if float(dist) <=cut_off_dist])
            epsilon_list = np.array([[math.sqrt(LJ_param[i][1]*LJ_param[j][1])]for i, j, dist in zip(i, j, dist) if float(dist) <=cut_off_dist])
            electronic_list = np.array([[(138.94*electronic_ratio*coulomb_param[i]*coulomb_param[j])/(water_epsilon*float(dist))] for i, j, dist in zip(i, j, dist) if float(dist) <=cut_off_dist])

            all_pair_LJ = 4*epsilon_list*((sigma_list**12)-(sigma_list**6))

            LJ_list.append(np.sum(all_pair_LJ))
            coulomb_list.append(np.sum(electronic_list))

            total_E = (np.array(LJ_list)+np.array(coulomb_list)).tolist()

            int_E_list.append(total_E)

        return LJ_list,coulomb_list,int_E_list

    else:
         for core in all_protein_cord:
            LJ_list.append(float(0))
            coulomb_list.append(float(0))

            int_E_list.append(float(0))
         return LJ_list,coulomb_list,int_E_list

def calc_interaction_additional_frag(all_atom,coord_array, all_protein_cord,split_ind_list_1,split_ind_list_2):
    LJ_list = []
    coulomb_list = []
    int_E_list = []
    cut_off_dist = 10
    electronic_ratio = 0.8
    LJ_param = {'H':[2.38760856462,0.1882800],'D':[0.40001352445,0.1924640],'C':[3.58141284692,0.2343040],'N':[3.52795892384,0.2928800],'O':[3.02905564168,0.5020800],'F':[3.02905564168,0.5020800],'P':[3.83086448800,2.4476400],'S':[3.56359487256,1.8828000],'Cl':[3.31414323148,0.962320],'Br':[3.52795892384,1.3388800],'I':[3.99122625727,2.1756800],'Mg':[2.11142996199,0.0627600],'Zn':[1.94215920555,1.0460000],'Vn':[0,0]}
    coulomb_param = {'H':0.09,'D':0.43,'C':0.07,'N':-0.4,'O':-0.5,'F':-0.19,'P':1.0,'S':0.12,'Cl':-0.17,'Br':-0.14,'I':-0.08,'Mg':2,'Zn':2,'Vn':0}
    water_epsilon = 78.5
    replace_dict = {'CQ': 'C', 'X1': 'C', 'X2': 'C', 'X3': 'C', 'X3': 'C','X4':'N', 'W1': 'N', 'A1': 'H', 'A2': 'H', 'A3': 'H', 'A3': 'H', 'B1': 'H', 'B2': 'H', 'B3': 'H', 'B4': 'H', 'NS': 'N', 'OS': 'O', 'SS': 'S', 'HS': 'H', 'CS': 'C', 'NQ': 'N', 'OQ': 'O', 'SQ': 'S', 'HQ': 'H', 'CQ': 'C', 'NP': 'N', 'OP': 'O', 'SP': 'S', 'HP': 'H', 'CP': 'C', 'VLA': 'Vn', 'VLB': 'Vn', 'VLC': 'Vn', 'VC_1': 'Vn', 'VC_2': 'Vn', 'VC_3': 'Vn', 'VC_4': 'Vn', 'Y1': 'Vn', 'Z1': 'Vn', 'RH_1': 'H', 'RH_2': 'H', 'RH_3': 'H', 'RH_4': 'H', 'RC_1': 'C', 'RC_2': 'C', 'RC_3': 'C', 'RC_4': 'C', 'RO_1': 'O', 'RO_2': 'O', 'RO_3': 'O', 'RO_4': 'O', 'RN_1': 'N', 'RN_2': 'N', 'RN_3': 'N', 'RN_4': 'N','RS_1': 'S', 'RS_2': 'S', 'RS_3': 'S', 'RS_4': 'S','YCA': 'C', 'YCO': 'C', 'YNA': 'N','YN2': 'N','YN3': 'N', 'YS': 'S', 'YSA': 'S', 'YO': 'O','XC3_1': 'C','XC3_2': 'C','XC3_3': 'C','XC3_4': 'C', 'XC3g_1': 'C','XC3g_2': 'C','XC3g_3': 'C','XC3g_4': 'C', 'XC2_1': 'C','XC2_2': 'C','XC2_3': 'C','XC2_4': 'C', 'XC1_1': 'C','XC1_2': 'C','XC1_3': 'C','XC1_4': 'C', 'XCA_1': 'C','XCA_2': 'C','XCA_3': 'C','XCA_4': 'C','CX_0': 'C','A4': 'N', 'B4': 'H','NX_0': 'N','DH': 'D', 'HA': 'H', 'HR': 'H', 'HK': 'H','RF_1': 'F','RD_1': 'D'}
    for core in all_protein_cord:
        names2, coords2 = np.array(core)[:, 0], np.array(core)[:, 1:].astype(float)
        distances = np.linalg.norm(coord_array[:, np.newaxis, :] - coords2, axis=-1)
        combinations = np.array(np.meshgrid(all_atom, names2)).T.reshape(-1, 2)
        for key, value in replace_dict.items():
            mask = combinations == key
            combinations[mask] = value
        i,j,dist = zip(*np.column_stack((combinations, distances.flatten())))
        sigma_list = np.array([[((LJ_param[i][0]+LJ_param[j][0])/2)/float(dist)] if float(dist) <=cut_off_dist else [float(0)] for i, j, dist in zip(i, j, dist)])
        epsilon_list = np.array([[math.sqrt(LJ_param[i][1]*LJ_param[j][1])] if float(dist) <=cut_off_dist else [float(0)] for i, j, dist in zip(i, j, dist)])
        electronic_list = np.array([[(138.94*electronic_ratio*coulomb_param[i]*coulomb_param[j])/(water_epsilon*float(dist))] if float(dist) <=cut_off_dist else [float(0)] for i, j, dist in zip(i, j, dist)])
        all_pair_LJ = 4*epsilon_list*((sigma_list**12)-(sigma_list**6))
        if split_ind_list_2 != None:
            split_ind_list_1st,split_ind_list_2nd = create_split_ind_list_for_interaction(split_ind_list_2,len(names2))
            LJ_list.append([np.sum(all_pair_LJ[a:b]) for a,b in zip(split_ind_list_1st,split_ind_list_2nd)])
            coulomb_list.append([np.sum(electronic_list[a:b]) for a,b in zip(split_ind_list_1st,split_ind_list_2nd)])
        else:
            LJ_list.append([np.sum(all_pair_LJ[:])])
            coulomb_list.append([np.sum(electronic_list[:])])
        total_E = (np.array(LJ_list)+np.array(coulomb_list)).tolist()
        int_E_list= [x for x in total_E]
    return LJ_list,coulomb_list,int_E_list

def calc_heavy_atoms(core_coord):
    heavy_atom = 0
    atom_list = [names for names, xs, ys, zs in core_coord]
    for row in atom_list:
        if row != 'H' and row != 'D' and row != 'VLA' and row != 'VLB' and row != 'VLC' and row != 'VC_1' and row != 'VC_2' and row != 'VC_3':
            heavy_atom += 1
    return heavy_atom

def calc_interaction_from_strategy(core_coord, all_protein_cord):
    coord_list = [[float(xs), float(ys), float(zs)] for names, xs, ys, zs in core_coord]
    atom_list = [names for names, xs, ys, zs in core_coord]
    all_atom = atom_list
    coord_array = np.array(coord_list)
    LJ_list = []
    coulomb_list = []
    int_E_list = []
    cut_off_dist = 10
    electronic_ratio = 0.8
    LJ_param = {'H':[2.38760856462,0.1882800],'D':[0.40001352445,0.1924640],'C':[3.58141284692,0.2343040],'N':[3.52795892384,0.2928800],'O':[3.02905564168,0.5020800],'F':[3.02905564168,0.5020800],'P':[3.83086448800,2.4476400],'S':[3.56359487256,1.8828000],'Cl':[3.31414323148,0.962320],'Br':[3.52795892384,1.3388800],'I':[3.99122625727,2.1756800],'Mg':[2.11142996199,0.0627600],'Zn':[1.94215920555,1.0460000],'Vn':[0,0]}
    coulomb_param = {'H':0.09,'D':0.43,'C':0.07,'N':-0.4,'O':-0.5,'F':-0.19,'P':1.0,'S':0.12,'Cl':-0.17,'Br':-0.14,'I':-0.08,'Mg':2,'Zn':2,'Vn':0}
    water_epsilon = 78.5
    replace_dict = {'CQ': 'C', 'X1': 'C', 'X2': 'C', 'X3': 'C', 'X3': 'C','X4':'N', 'W1': 'N', 'A1': 'H', 'A2': 'H', 'A3': 'H', 'A3': 'H', 'B1': 'H', 'B2': 'H', 'B3': 'H', 'B4': 'H', 'NS': 'N', 'OS': 'O', 'SS': 'S', 'HS': 'H', 'CS': 'C', 'NQ': 'N', 'OQ': 'O', 'SQ': 'S', 'HQ': 'H', 'CQ': 'C', 'NP': 'N', 'OP': 'O', 'SP': 'S', 'HP': 'H', 'CP': 'C', 'VLA': 'Vn', 'VLB': 'Vn', 'VLC': 'Vn', 'VC_1': 'Vn', 'VC_2': 'Vn', 'VC_3': 'Vn', 'VC_4': 'Vn', 'Y1': 'Vn', 'Z1': 'Vn', 'RH_1': 'H', 'RH_2': 'H', 'RH_3': 'H', 'RH_4': 'H', 'RC_1': 'C', 'RC_2': 'C', 'RC_3': 'C', 'RC_4': 'C', 'RO_1': 'O', 'RO_2': 'O', 'RO_3': 'O', 'RO_4': 'O', 'RN_1': 'N', 'RN_2': 'N', 'RN_3': 'N', 'RN_4': 'N','RS_1': 'S', 'RS_2': 'S', 'RS_3': 'S', 'RS_4': 'S','YCA': 'C', 'YCO': 'C', 'YNA': 'N','YN2': 'N','YN3': 'N', 'YS': 'S', 'YSA': 'S', 'YO': 'O','XC3_1': 'C','XC3_2': 'C','XC3_3': 'C','XC3_4': 'C', 'XC3g_1': 'C','XC3g_2': 'C','XC3g_3': 'C','XC3g_4': 'C', 'XC2_1': 'C','XC2_2': 'C','XC2_3': 'C','XC2_4': 'C', 'XC1_1': 'C','XC1_2': 'C','XC1_3': 'C','XC1_4': 'C', 'XCA_1': 'C','XCA_2': 'C','XCA_3': 'C','XCA_4': 'C','CX_0': 'C','A4': 'N', 'B4': 'H','NX_0': 'N','DH': 'D', 'HA': 'H', 'HR': 'H', 'HK': 'H', 'RF_1': 'F', 'RD_1': 'D'}
    for core in all_protein_cord:
        names2, coords2 = np.array(core)[:, 0], np.array(core)[:, 1:].astype(float)
        distances = np.linalg.norm(coord_array[:, np.newaxis, :] - coords2, axis=-1)
        combinations = np.array(np.meshgrid(all_atom, names2)).T.reshape(-1, 2)
        for key, value in replace_dict.items():
            mask = combinations == key
            combinations[mask] = value
        i,j,dist = zip(*np.column_stack((combinations, distances.flatten())))
        sigma_list = np.array([[((LJ_param[i][0]+LJ_param[j][0])/2)/float(dist)] if float(dist) <=cut_off_dist else [float(0)] for i, j, dist in zip(i, j, dist)])
        epsilon_list = np.array([[math.sqrt(LJ_param[i][1]*LJ_param[j][1])] if float(dist) <=cut_off_dist else [float(0)] for i, j, dist in zip(i, j, dist)])
        electronic_list = np.array([[(138.94*electronic_ratio*coulomb_param[i]*coulomb_param[j])/(water_epsilon*float(dist))] if float(dist) <=cut_off_dist else [float(0)] for i, j, dist in zip(i, j, dist)])
        all_pair_LJ = 4*epsilon_list*((sigma_list**12)-(sigma_list**6))
        LJ_list.append([np.sum(all_pair_LJ[:])])
        coulomb_list.append([np.sum(electronic_list[:])])
        total_E = (np.array(LJ_list)+np.array(coulomb_list)).tolist()
        int_E_list= [x for x in total_E]
    return LJ_list,coulomb_list,int_E_list

def check_collision_lig_fin_v2_1(names1, coords1, names2_all, coords2_all,split_ind_list_1,split_ind_list_2,temp_frag_list,sigma_sum_cut_init):
    replace_dict = {'CQ': 'C', 'X1': 'C', 'X2': 'C', 'X3': 'C', 'X3': 'C','X4':'C', 'W1': 'N', 'A1': 'H', 'A2': 'H', 'A3': 'H', 'A3': 'H', 'B1': 'H', 'B2': 'H', 'B3': 'H', 'B4': 'H', 'NS': 'N', 'OS': 'O', 'SS': 'S', 'HS': 'H', 'CS': 'C', 'NQ': 'N', 'OQ': 'O', 'SQ': 'S', 'HQ': 'H', 'CQ': 'C', 'NP': 'N', 'OP': 'O', 'SP': 'S', 'HP': 'H', 'CP': 'C', 'VLA': 'Vn', 'VLB': 'Vn', 'VLC': 'Vn', 'VC_1': 'Vn', 'VC_2': 'Vn', 'VC_3': 'Vn', 'VC_4': 'Vn', 'Y1': 'Vn', 'Z1': 'Vn', 'RH_1': 'H', 'RH_2': 'H', 'RH_3': 'H', 'RH_4': 'H', 'RC_1': 'C', 'RC_2': 'C', 'RC_3': 'C', 'RC_4': 'C', 'RO_1': 'O', 'RO_2': 'O', 'RO_3': 'O', 'RO_4': 'O', 'RN_1': 'N', 'RN_2': 'N', 'RN_3': 'N', 'RN_4': 'N','RS_1': 'S', 'RS_2': 'S', 'RS_3': 'S', 'RS_4': 'S','YCA': 'C', 'YCO': 'C', 'YNA': 'N','YN2': 'N','YN3': 'N', 'YS': 'S', 'YSA': 'S', 'YO': 'O','XC3_1': 'C','XC3_2': 'C','XC3_3': 'C','XC3_4': 'C', 'XC3g_1': 'C','XC3g_2': 'C','XC3g_3': 'C','XC3g_4': 'C', 'XC2_1': 'C','XC2_2': 'C','XC2_3': 'C','XC2_4': 'C', 'XC1_1': 'C','XC1_2': 'C','XC1_3': 'C','XC1_4': 'C', 'XCA_1': 'C','XCA_2': 'C','XCA_3': 'C','XCA_4': 'C','CX_0': 'C','A4': 'Vn', 'B4': 'H', 'XC5': 'C','NX_0': 'N','DH': 'D', 'HA': 'H', 'HR': 'H', 'HK': 'H', 'RF_1': 'F', 'RD_1': 'D'}

    LJ_param = {'H':[2.38760856462,0.1882800],'D':[0.40001352445,0.1924640],'C':[3.58141284692,0.2343040],'N':[3.52795892384,0.2928800],'O':[3.02905564168,0.5020800],'F':[3.02905564168,0.5020800],'P':[3.83086448800,2.4476400],'S':[3.56359487256,1.8828000],'Cl':[3.31414323148,0.962320],'Br':[3.52795892384,1.3388800],'I':[3.99122625727,2.1756800],'Mg':[2.11142996199,0.0627600],'Zn':[1.94215920555,1.0460000],'Vn':[0,0]}
    sigma_sum_cut = sigma_sum_cut_init
    collision_list = []
    if split_ind_list_1 != None or split_ind_list_1 != None:
        for a,b,c in zip(split_ind_list_1,split_ind_list_2,temp_frag_list):
            if c == 'sym24_H1L':
                collision_list.append(False)
            else:
                if c == 'sym24_XHD_0' or c == 'sym24_XHD_1' or c == 'sym24_XHD_2' or c == 'sym24_XHD_3' or c == 'sym24_XHD_4':
                    sigma_sum_cut = 10
                names2 = names2_all[a:b]
                coords2 = coords2_all[a:b]
                distances = np.linalg.norm(np.array(coords1)[:, np.newaxis, :] - np.array(coords2), axis=-1)
                combinations = np.array(np.meshgrid(names1, names2)).T.reshape(-1, 2)
                for key, value in replace_dict.items():
                    mask = combinations == key
                    combinations[mask] = value
                i,j,dist = zip(*np.column_stack((combinations, distances.flatten())))
                sigma_list = np.array([[((LJ_param[i][0]+LJ_param[j][0])/2)/float(dist)] for i, j, dist in zip(i, j, dist)])
                epsilon_list = np.array([[math.sqrt(LJ_param[i][1]*LJ_param[j][1])] for i, j, dist in zip(i, j, dist)])
                all_pair_LJ = (4*epsilon_list*((sigma_list**12)-(sigma_list**6))).flatten()
                all_pair_LJ = np.sort(all_pair_LJ)[::-1]
                all_pair_LJ = all_pair_LJ[7:]
                if sum(all_pair_LJ) != np.nan:
                    if sum(all_pair_LJ) <= sigma_sum_cut:
                        collision_list.append(False)
                    else:
                        collision_list.append(True)
                else:
                    collision_list.append(True)
    else:
        if temp_frag_list == 'sym24_H1L':
            collision_list.append(False)
        else:
            if temp_frag_list == 'sym24_XHD_0' or temp_frag_list == 'sym24_XHD_1' or temp_frag_list == 'sym24_XHD_2' or temp_frag_list == 'sym24_XHD_3' or temp_frag_list == 'sym24_XHD_4':
                sigma_sum_cut = 10
            names2 = names2_all[:]
            coords2 = coords2_all[:]
            distances = np.linalg.norm(np.array(coords1)[:, np.newaxis, :] - np.array(coords2), axis=-1)
            combinations = np.array(np.meshgrid(names1, names2)).T.reshape(-1, 2)
            for key, value in replace_dict.items():
                mask = combinations == key
                combinations[mask] = value
            i,j,dist = zip(*np.column_stack((combinations, distances.flatten())))
            sigma_list = np.array([[((LJ_param[i][0]+LJ_param[j][0])/2)/float(dist)] for i, j, dist in zip(i, j, dist)])
            epsilon_list = np.array([[math.sqrt(LJ_param[i][1]*LJ_param[j][1])] for i, j, dist in zip(i, j, dist)])
            all_pair_LJ = (4*epsilon_list*((sigma_list**12)-(sigma_list**6))).flatten()
            all_pair_LJ = np.sort(all_pair_LJ)[::-1]
            all_pair_LJ = all_pair_LJ[7:]
            if sum(all_pair_LJ) != np.nan:
                if sum(all_pair_LJ) <= sigma_sum_cut:
                    collision_list.append(False)
                else:
                    collision_list.append(True)
            else:
                collision_list.append(True)
    return collision_list

def my_logfile(file_path, name,mode):
    logger_INFO = getLogger(name)
    logger_INFO.setLevel(logging.DEBUG)
    handler_format = Formatter('%(message)s on %(asctime)s',datefmt="%Y-%m-%d %H:%M:%S")
    file_handler = FileHandler(file_path, mode)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(handler_format)
    logger_INFO.addHandler(file_handler)
    return logger_INFO

def my_logfile_notime(file_path, name, mode):
    logger_notime = getLogger(name)
    logger_notime.setLevel(logging.DEBUG)
    handler_format_notime = Formatter('%(message)s')
    file_handler = FileHandler(file_path, mode)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(handler_format_notime)
    logger_notime.addHandler(file_handler)
    return logger_notime

def sum_set_info(cycle,des_frag_cycle,max_branch_core,max_branch_frag,result_conf,max_conf,cutoff_score_init,cutoff_score_fin,rollout_num,rollout_frag_max,fragment_set,rollout_set,desired_frag,desired_frag_rollout,evaluation_param,evaluation_method,core_dir_list,protein_xyz_dir_list):
    fragment_set_info = []
    rollout_set_info = []
    for x in os.listdir(fragment_set):
        if x != '__init__.py' and x != '__pycache__':
            fragment_set_list = os.listdir(fragment_set+'/'+x)
            for y in fragment_set_list:
                each_frag_name = x+'/'+y
                fragment_set_info.append(each_frag_name)
    for x in os.listdir(rollout_set):
        if x != '__init__.py' and x != '__pycache__':
            rollout_set_list = os.listdir(rollout_set+'/'+x)
            for y in rollout_set_list:
                each_frag_name = x+'/'+y
                rollout_set_info.append(each_frag_name)
    fragment_set_info.sort()
    rollout_set_info.sort()
    information = \
    '### Basic Setting ###'+'\n'+\
    '\tcycle: '+str(cycle)+'\n'+\
    '\tmax_branch_core: '+str(max_branch_core)+'\n'+\
    '\tmax_branch_frag: '+str(max_branch_frag)+'\n'+\
    '\n'+\
    '\tresult_conf: '+str(result_conf)+'\n'+\
    '\tmax_conf: '+str(max_conf)+'\n'+\
    '\n'+\
    '\tevaluation_param: '+str(evaluation_param)+'\n'+\
    '\tevaluation_method: '+str(evaluation_method)+'\n'+\
    '\tcutoff_score_init: '+str(cutoff_score_init)+'\n'+\
    '\tcutoff_score_fin: '+str(cutoff_score_fin)+'\n'+\
    '\n'+\
    '\trollout_num: '+str(rollout_num)+'\n'+\
    '\trollout_frag_max: '+str(rollout_frag_max)+'\n'+\
    '\n'+\
    '### Desired Fragment Setting ###'+'\n'+\
    '\tdes_frag_cycle: '+str(des_frag_cycle)+'\n'+\
    '\tdesired_frag: '+str(desired_frag[0])+'\n'+\
    '\tdesired_frag_rollout: '+str(desired_frag_rollout[0])+'\n'+\
    '\n'+\
    '### Protein Information ###'+'\n'
    for x in protein_xyz_dir_list:
        information += ('\t'+str(x)+'\n')
    information += '\n### Core Information ###'+'\n'
    for x in core_dir_list:
        information += ('\t'+str(x)+'\n')
    information += '\n### Fragments Information ###'+'\n'
    information += '\tdirectory: '+str(fragment_set)+'\n'
    for x in fragment_set_info:
        information += ('\t\t'+str(x)+'\n')
    information += '\n### Rollout Fragments Information ###'+'\n'
    information += '\tdirectory: '+str(rollout_set)+'\n'
    for x in rollout_set_info:
        information += ('\t\t'+str(x)+'\n')
    return information

def calc_performance(count_generated,elapsed_time):
    performance_1 = int(86400/(elapsed_time/count_generated))
    performance = 'performance: ' +str(performance_1)+' structures/day'
    return performance

def param_for_logfile(transforming_core_point_int,current_cycle_int,current_branch_core_int,current_branch_frag_int,ref_frag_list_next_int,param_for_next_int,total_charge_int,transform_log_main_int,transform_log_side_int,c_for_next_int,all_frag_int,all_atom_int):
    information = \
    'frag_name: '+str(all_frag_int)+\
    '\nfrag atoms: '+str(all_atom_int)+\
    '\nintermediate coordinate: '
    for i in transforming_core_point_int:
        information += ('\n'+str(i))
    information_2 = \
    '\ncurrent_cycle: '+str(current_cycle_int)+\
    '\ncurrent_branch_core: '+str(current_branch_core_int)+\
    '\ncurrent_branch_frag: '+str(current_branch_frag_int)+\
    '\nref_frag_list_next: '+str(ref_frag_list_next_int)+\
    '\nparam_for_next: '+str(param_for_next_int)+\
    '\ntotal_charge: '+str(total_charge_int)+\
    '\ntransform_log_main: '+str(transform_log_main_int)+\
    '\ntransform_log_side: '+str(transform_log_side_int)+\
    '\nc_for_next: '+str(c_for_next_int)
    information = information + information_2
    return information

def read_result_score(path_result_score):

    with open(path_result_score) as result_score:
        reader = csv.reader(result_score)
        list_row = result_score.read().split('\n')
    result_score = [x.split('"') for x in list_row if x != '']
    original_ID = ['ID_'+str(y[0].split()[0]) for y in result_score][1:]
    strategy = [parse_nested_list(x[-2])[0] for x in result_score[1:]]
    return strategy,original_ID

def parse_nested_list(input_str):
    nested_list = []
    stack = []
    current = nested_list
    i = 0
    n = len(input_str)
    while i < n:
        char = input_str[i]
        if char == '[':
            new_list = []
            current.append(new_list)
            stack.append(current)
            current = new_list
            i += 1
        elif char == ']':
            if not stack:
                raise ValueError("error")
            current = stack.pop()
            i += 1
        elif char == "'":
            i += 1
            start = i
            while i < n and input_str[i] != "'":
                i += 1
            if i >= n:
                raise ValueError("error")
            string_value = input_str[start:i]
            current.append(string_value)
            i += 1
        elif char.isdigit() or (char == '-' and i + 1 < n and input_str[i + 1].isdigit()):
            start = i
            if char == '-':
                i += 1
            while i < n and (input_str[i].isdigit() or input_str[i] == '.'):
                i += 1
            num_str = input_str[start:i]
            if '.' in num_str:
                num = float(num_str)
            else:
                num = int(num_str)
            current.append(num)
        elif char.isspace():
            i += 1
        else:
            i += 1
    if stack:
        raise ValueError("error")
    return nested_list

def remove_vn(cord_list):
    cord_list = [x for x in cord_list if x[0] != 'Vn']
    return cord_list

def import_core(dir_core):
    count = 0
    core_name = dir_core+'.xyz'
    core_file = os.path.basename(core_name)
    initial_remaining_data, S, Q, P = import_structure_2(dir_core,core_file)
    transform_log_main,transform_log_side,c_for_next,param_for_next,ref_frag_list_next = get_initial_x(initial_remaining_data)
    added_frag = ''
    charge_file = dir_core+'/charge.txt'
    charge_core = extract_charge(charge_file)
    total_charge = int(charge_core)

    return initial_remaining_data,S,Q,P,total_charge,transform_log_main,transform_log_side,c_for_next,param_for_next,added_frag,ref_frag_list_next

def get_core_charge(dir_core):
    charge_file = dir_core+'/charge.txt'
    charge_core = extract_charge(charge_file)
    total_charge = int(charge_core)
    return total_charge

def calc_evaluation_v2(LJ_array_sum,coulomb_array_sum,int_E_array_sum,evaluation_param,evaluation_method,heavy_atom_list_init):
    heavy_atom_list = [x if x != 0 else 1 for x in heavy_atom_list_init]
    evaluation = []
    LE_LJ_array_sum = LJ_array_sum/np.array(heavy_atom_list)
    LE_coulomb_array_sum = coulomb_array_sum/np.array(heavy_atom_list)
    LE_int_E_array_sum = int_E_array_sum/np.array(heavy_atom_list)
    LE_LJ_coulomb_sum = LE_LJ_array_sum + coulomb_array_sum
    selectivity_constant = 2
    if evaluation_param == 'LE_LJ+coulomb':
        if evaluation_method == 'average':
            evaluation_array = np.mean(LE_LJ_coulomb_sum,axis=0)
        elif evaluation_method == 'highest':
            evaluation_array = np.amin(LE_LJ_coulomb_sum,axis=0)
        elif evaluation_method == 'selective':
            evaluation_array = np.amin(LE_LJ_coulomb_sum,axis=0) + (selectivity_constant*(np.amin(LE_LJ_coulomb_sum-LE_LJ_coulomb_sum[0],axis=0)))
    if evaluation_param == 'LE_interact_E':
        if evaluation_method == 'average':
            evaluation_array = np.mean(LE_int_E_array_sum,axis=0)
        elif evaluation_method == 'highest':
            evaluation_array = np.amin(LE_int_E_array_sum,axis=0)
        elif evaluation_method == 'selective':
            evaluation_array = np.amin(LE_int_E_array_sum,axis=0) + (selectivity_constant*(np.amin(LE_int_E_array_sum-LE_int_E_array_sum[0],axis=0)))
    if evaluation_param == 'interact_E':
        if evaluation_method == 'average':
            evaluation_array = np.mean(int_E_array_sum,axis=0)
        elif evaluation_method == 'highest':
            evaluation_array = np.amin(int_E_array_sum,axis=0)
        elif evaluation_method == 'selective':
            evaluation_array = np.amin(int_E_array_sum,axis=0) + (selectivity_constant*(np.amin(int_E_array_sum-int_E_array_sum[0],axis=0)))
    if evaluation_param == 'LE_LJ':
        if evaluation_method == 'average':
            evaluation_array = np.mean(LE_LJ_array_sum,axis=0)
        elif evaluation_method == 'highest':
            evaluation_array = np.amin(LE_LJ_array_sum,axis=0)
        elif evaluation_method == 'selective':
            evaluation_array = np.amin(LE_LJ_array_sum,axis=0) + (selectivity_constant*(np.amin(LE_LJ_array_sum-LE_LJ_array_sum[0],axis=0)))
    if evaluation_param == 'LE_LJ':
        if evaluation_method == 'average':
            evaluation_array = np.mean(LJ_array_sum,axis=0)
        elif evaluation_method == 'highest':
            evaluation_array = np.amin(LJ_array_sum,axis=0)
        elif evaluation_method == 'selective':
            evaluation_array = np.amin(LJ_array_sum,axis=0) + (selectivity_constant*(np.amin(LJ_array_sum-LJ_array_sum[0],axis=0)))
    return evaluation_array

def execute_rollout_v2(transforming_core_point,S,Q,P,total_charge_init,transform_log_main_init,transform_log_side_init,c_for_next_init,param_for_next_init,added_frag,ref_frag_list_next_init,current_cycle_init,cycle,current_branch_core_init,current_branch_frag_init,max_branch_core,max_branch_frag,evaluation_param,evaluation_method,desired_list,rollout_set,all_protein_cord,add_cycle,state_for_rollout,rol_des_category,des_frag_cycle,dir_core,cutoff_score_init,rollout_conf_inner,rollout_conf_inner_des,sigma_sum_cut):
    rollout_dict = get_transform_dict_hyd_rollout(rollout_set)


    rol_Frag_F1L = Fragments('F1L',rollout_dict)
    rol_Frag_F2L = Fragments('F2L',rollout_dict)
    rol_Frag_C1L = Fragments('C1L',rollout_dict)
    rol_Frag_C2L = Fragments('C2L',rollout_dict)
    rol_Frag_dict = {'F1L': rol_Frag_F1L, 'F2L': rol_Frag_F2L,'C1L': rol_Frag_C1L, 'C2L': rol_Frag_C2L}

    if 'F3L' in os.listdir(rollout_set):
        rol_Frag_F3L = Fragments('F3L',rollout_dict)
        rol_Frag_dict['F3L'] = rol_Frag_F3L
    if 'F2LB' in os.listdir(rollout_set):
        rol_Frag_F2LB = Fragments('F2LB',rollout_dict)
        rol_Frag_dict['F2LB'] = rol_Frag_F2LB
    if 'C3L' in os.listdir(rollout_set):
        rol_Frag_C3L = Fragments('C3L',rollout_dict)
        rol_Frag_dict['C3L'] = rol_Frag_C3L
    if 'XHD' in os.listdir(rollout_set):
        rol_Frag_XHD = Fragments('XHD',rollout_dict)
        rol_Frag_dict['XHD'] = rol_Frag_XHD
    if 'XFO' in os.listdir(rollout_set):
        rol_Frag_XFO = Fragments('XFO',rollout_dict)
        rol_Frag_dict['XFO'] = rol_Frag_XFO

    rol_Frag_desired = desFragments(desired_list,rol_des_category)



    if des_frag_cycle < 1:
        rollout_conf = rollout_conf_inner
    else:
        rollout_conf = rollout_conf_inner_des
    rollout_eval_list = []


    all_state = []

    depth = 0
    generated_conf = 0

    next_state = []




    while generated_conf < rollout_conf:

        EOF_value = False
        if depth ==0:


            current_cycle = current_cycle_init
            total_charge = total_charge_init
            transform_log_main = transform_log_main_init
            transform_log_side = transform_log_side_init
            c_for_next = c_for_next_init
            param_for_next = param_for_next_init
            ref_frag_list_next = ref_frag_list_next_init
            current_branch_core = current_branch_core_init
            current_branch_frag = current_branch_frag_init

            core = Core_mol(transforming_core_point,current_cycle_init,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next)


            if des_frag_cycle < 1:
                if 100 in rollout_eval_list:
                    generated_conf = rollout_conf
                    break


            next_state = []

            first_list = [sublist[0] for sublist in all_state]

            if not any(dir_core == x for x in first_list):

                temp_state = [dir_core,0,0]
                next_state.append(temp_state)
                all_state.append(temp_state)

            else:
                ind = first_list.index(dir_core)
                temp_state = all_state[ind]
                next_state.append(temp_state)



        to_frons_state = next_state

        temp_state = next_state
        current_state = []
        all_pattern_list = []

        if current_cycle != des_frag_cycle:
            all_atom = [x for i in core.operation for x in rol_Frag_dict[i].atom]
            all_coord = [x for i in core.operation for x in rol_Frag_dict[i].coord]
            all_charge = [x for i in core.operation for x in rol_Frag_dict[i].charge]
            all_each_frag = [x for i in core.operation for x in rol_Frag_dict[i].each_frag]
            all_len_list = [x for i in core.operation for x in rol_Frag_dict[i].len_list]
            all_central_atom = [x for i in core.operation for x in rol_Frag_dict[i].central_atom]
            all_category = [x for i in core.operation for x in rol_Frag_dict[i].category]
            all_param = [x for i in core.operation for x in rol_Frag_dict[i].param]
            all_frag = [x for i in core.operation for x in rol_Frag_dict[i].frag_name]
            all_frag_category_val = [x for i in core.operation for x in rol_Frag_dict[i].frag_category_val]
            all_add_branch_core = [x for i in core.operation for x in rol_Frag_dict[i].add_branch_core]
            all_add_branch_frag = [x for i in core.operation for x in rol_Frag_dict[i].add_branch_frag]
            all_ref_frag_list = [x for i in core.operation for x in rol_Frag_dict[i].ref_frag_list]
            frag_c_for_next = [x for i in core.operation for x in rol_Frag_dict[i].c_for_next]

            split_ind_list_1,split_ind_list_2 = create_split_ind_list(all_len_list)
            possibility = True

        elif current_cycle == des_frag_cycle and rol_Frag_desired.frag_category_val[0] in core.operation:
            all_atom = [x for x in rol_Frag_desired.atom]
            all_coord = [x for x in rol_Frag_desired.coord]
            all_charge = [x for x in rol_Frag_desired.charge]
            all_each_frag = [x for x in rol_Frag_desired.each_frag]
            all_len_list = [x for x in rol_Frag_desired.len_list]
            all_central_atom = [x for x in rol_Frag_desired.central_atom]
            all_category = [x for x in rol_Frag_desired.category]
            all_param = [x for x in rol_Frag_desired.param]
            all_frag = [x for x in rol_Frag_desired.frag_name]
            all_add_branch_core = [x for x in rol_Frag_desired.add_branch_core]
            all_add_branch_frag = [x for x in rol_Frag_desired.add_branch_frag]
            all_ref_frag_list = [x for x in rol_Frag_desired.ref_frag_list]
            all_frag_category_val = [x for x in rol_Frag_desired.frag_category_val]
            split_ind_list_1,split_ind_list_2 = create_split_ind_list(all_len_list)
            frag_c_for_next = [x for x in rol_Frag_desired.c_for_next]
            possibility = True

        else:

            escape_state = copy.deepcopy(to_frons_state)
            first_list = [sublist[0] for sublist in all_state]
            if depth == 0:
                escape_ind= first_list.index(escape_state[0][0])
            else:
                escape_ind= first_list.index(escape_state[0])
            all_state[escape_ind][1] = 100
            all_state[escape_ind][2] = 10000000
            possibility = False

        if possibility:

            if depth == 0:

                temp_state_0 = temp_state[0][0]
                first_elements = [sublist[0] for sublist in all_state if type(sublist[0])==str and sublist[2] != 0]
            else:
                temp_state_0 = temp_state[0]
                first_elements = [sublist[0] for sublist in all_state if len(sublist[0])==(len(temp_state_0)+1) and sublist[0][:len(temp_state_0)]== temp_state_0]
            temp_result = link_result(core,all_atom,all_coord,all_charge,all_central_atom,split_ind_list_1,split_ind_list_2,all_category,all_param,all_len_list,all_frag,all_protein_cord,evaluation_param,evaluation_method,all_ref_frag_list,frag_c_for_next,sigma_sum_cut)

            if first_elements == [] or depth == 0:

                temp_result_list = temp_result.evaluation.tolist()
                if current_cycle <= des_frag_cycle:
                    current_state = [[a,b,1] if c == False else [a,100,1] for a,b,c in zip(all_each_frag,temp_result_list,temp_result.eof_val)]
                else:
                    current_state = [[a,b,1] for a,b in zip(all_each_frag,temp_result_list)]

                each_state = [[[next_state[0][0],sublist[0]]] + sublist[1:] if depth == 0 else [next_state[0]+[sublist[0]]] + sublist[1:] for sublist in current_state]

                all_state.extend(each_state)

            else:

                if depth == 0:
                    current_state = [sublist[0] for sublist in all_state if type(sublist[0])==str]
                else:
                    int_state = [sublist for sublist in all_state if len(sublist[0])==(len(temp_state_0)+1) and sublist[0][:len(temp_state_0)]== temp_state_0]
                    current_state = [[a[0][-1],a[1],a[2]] for a in int_state]

            sorted_current_state = sorted(current_state, key=lambda x: (x[2],x[1]))
            current_score_list = [sublist[1] for sublist in current_state]
            best_index = current_state.index(sorted_current_state[0])

            min_score = current_state[best_index][1]

            if min_score < 100:

                next_frag_atom = temp_result.atom[split_ind_list_1[best_index]:split_ind_list_2[best_index]]
                next_frag_coord = temp_result.coord[split_ind_list_1[best_index]:split_ind_list_2[best_index]]
                EOF_value = temp_result.eof_val[best_index]
                category = temp_result.temp_category[best_index]
                temp_frag = temp_result.temp_frag[best_index]
                ref_frag_list_next = temp_result.ref_frag_list_next_list[best_index]

                first_list = [sublist[0] for sublist in all_state]

                if depth == 0:

                    all_state_ind = first_list.index([next_state[0][0],current_state[best_index][0]])
                else:
                    all_state_ind = first_list.index(next_state[0]+[current_state[best_index][0]])
                if EOF_value == False:
                    all_state[all_state_ind][2] +=1
                    depth += 1
                    next_state = all_state[all_state_ind]
                    new_core_structure = create_new_core_structure(core,next_frag_atom,next_frag_coord,category,temp_frag,ref_frag_list_next)
                    current_cycle = core.add_cycle+core.current_cycle
                    current_branch_core = core.current_branch_core + all_add_branch_core[best_index]
                    current_branch_frag = core.current_branch_frag + all_add_branch_frag[best_index]
                    ref_frag_list_next = temp_result.ref_frag_list_next_list[best_index]
                    param_for_next = temp_result.param_for_next[best_index]
                    total_charge = temp_result.total_charge[best_index]
                    transform_log_main = temp_result.transform_log_main_list[best_index]
                    transform_log_side = temp_result.transform_log_side_list[best_index]

                    if temp_result.temp_frag[best_index] != 'sym24_XHD_0' and temp_result.temp_frag[best_index] != 'sym24_XHD_1' and temp_result.temp_frag[best_index] != 'sym24_XHD_2' and temp_result.temp_frag[best_index] != 'sym24_XHD_3' and temp_result.temp_frag[best_index] != 'sym24_XHD_4':
                        c_for_next = temp_result.c_for_next_list[best_index]
                    else:
                        c_for_next = ['XCA_1']

                    del temp_result
                    del core

                    core = Core_mol(new_core_structure,current_cycle,cycle,current_branch_core,current_branch_frag,max_branch_core,max_branch_frag,ref_frag_list_next,param_for_next,total_charge,transform_log_main,transform_log_side,c_for_next)


                elif EOF_value == True:

                    if current_cycle <= des_frag_cycle:

                        if depth >= 1:
                            escape_ind = first_list.index(next_state[0][0])
                            all_state[escape_ind][1] = 100
                            all_state[escape_ind][2] = 10000
                            if depth != 1:
                                rollout_eval_list.append(100)
                                generated_conf += 1
                            depth = 0
                        else:
                            rollout_eval_list.append(100)
                            generated_conf += 1
                            depth = 0

                            pass
                    else:
                        generated_conf +=1
                        all_state[all_state_ind][2] +=10000
                        depth = 0

                        rollout_eval_list.append(current_state[best_index][1])
                        if current_state[best_index][1] <= cutoff_score_init:
                            generated_conf = rollout_conf

            else:
                if depth >= 1:
                    first_list_to_escape = [sublist[0] for sublist in all_state]

                    escape_ind = first_list_to_escape.index(temp_state_0)

                    all_state[escape_ind][1] = 100
                    all_state[escape_ind][2] = 10000
                    generated_conf += 1
                    rollout_eval_list.append(100)
                    depth = 0

                else:
                    rollout_eval_list.append(100)
                    generated_conf += 1
                    depth = 0

        if possibility == False:
            rollout_eval_list.append(100)
            generated_conf += 1
            depth = 0

    score = min(rollout_eval_list)

    return score


def create_new_core_structure(Core_mol,frag_atom,frad_array,category,temp_frag,ref_frag_list_next):
    core_atom = copy.deepcopy(Core_mol.atom)
    core_coord = copy.deepcopy(Core_mol.coord)
    frad_coord = frad_array.tolist()

    if temp_frag == 'sym24_H1L':
        if category == ['link_f']:

            for n in reversed(range(1,4)):
                for x in ['XC3_','XC3g_','XC2_','XC1_','XCA_']:
                    X_name = x+str(n)

                    if any(X_name in element for element in core_atom):
                        X_name_at_front = X_name
                        n_at_front = copy.deepcopy(n)

                        break
                    else:
                        continue
                else:
                    continue
                break
            for r in ['RH_','RC_','RO_','RN_','RS_','RCl_','RBr_','RI_','RP_','RD_']:
                R_name = r+str(n_at_front)
                if any(R_name in element for element in core_atom):
                    R_name_at_front = R_name
                    break
                else:
                    continue

            for vc in ['VC_3','VC_2','VC_1']:
                if any(vc in element for element in core_atom):
                    VC_name_at_front = vc
                    core_atom,core_coord = remove_VC_after_link_h_f(core_atom,core_coord,VC_name_at_front)
                    break
            for v in ['VLC','VLB','VLA']:
                if any(v in element for element in core_atom):
                    V_name_at_front = v
                    core_atom,core_coord = remove_VC_after_link_h(core_atom,core_coord,V_name_at_front)
                    break
                else:
                    continue


            core_atom = atomname_change_after_link_f(core_atom,X_name_at_front,R_name_at_front,n_at_front)

            core_atom.extend('H')
            core_coord.extend(frad_coord)

        elif category == ['link_c']:

            for v in ['VLC','VLB','VLA']:
                if any(v in element for element in core_atom):
                    V_name_at_front = v
                    core_atom,core_coord = remove_VC_after_link_h(core_atom,core_coord,V_name_at_front)
                    break
                else:
                    continue

            for vc in ['VC_3','VC_2','VC_1']:
                if any(vc in element for element in core_atom):
                    VC_name_at_front = vc
                    core_atom,core_coord = remove_VC_after_link_h_f(core_atom,core_coord,VC_name_at_front)
                    break

            core_atom = atomname_change_after_link_h(core_atom)


            core_atom.extend(frag_atom)
            core_coord.extend(frad_coord)

    else:
        if category == ['link_c']:

            for x in ['XC3_1','XC3g_1','XC2_1','XC1_1','XCA_1']:
                if any(x in element for element in frag_atom):
                    X_name_at_front = x
                    break
                else:
                    continue

            for r in ['RH_1','RC_1','RO_1','RN_1','RS_1','RCl_1','RBr_1','RI_1','RP_1','RD_1']:
                if any(r in element for element in frag_atom):
                    R_name_at_front = r
                    break

            for v in ['VLC','VLB','VLA']:
                if any(v in element for element in core_atom):
                    V_name_at_front = v
                    break
                else:
                    continue

            VC_name_at_front = 'VC_1'

            f_atom_list = ['YCA','YO','YNA','YN3','YN2','YS','YSA']
            for f_atom in f_atom_list:
                if any(f_atom in element for element in frag_atom):
                    core_atom = atomname_change_before_link_n_c(core_atom)
                    core_atom,core_coord = remove_V_only_after_link_c(core_atom,core_coord,V_name_at_front)

                    core_atom.extend(frag_atom)
                    core_coord.extend(frad_coord)

                    if len(ref_frag_list_next) ==1:
                        core_atom = atomname_change_after_link_c_f(core_atom,X_name_at_front,R_name_at_front)
                    break
            else:
                core_atom,core_coord = remove_V_only_after_link_c(core_atom,core_coord,V_name_at_front)
                core_atom.extend(frag_atom)
                core_coord.extend(frad_coord)
                core_atom = atomname_change_after_link_c(core_atom,X_name_at_front,R_name_at_front)
                core_atom = atomname_change_before_link_n(core_atom)

                core_atom,core_coord = remove_VC_only_after_link_c(core_atom,core_coord,VC_name_at_front)

        elif category == ['link_f'] and temp_frag != 'sym24_XHD_0' and temp_frag != 'sym24_XHD_1' and temp_frag != 'sym24_XHD_2' and temp_frag != 'sym24_XHD_3' and temp_frag != 'sym24_XHD_4' and temp_frag != 'XFO':

            for n in reversed(range(1,4)):
                for x in ['XC3_','XC3g_','XC2_','XC1_','XCA_']:
                    X_name = x+str(n)

                    if any(X_name in element for element in core_atom):
                        X_name_at_front = X_name
                        n_at_front = copy.deepcopy(n)

                        break
                    else:
                        continue
                else:
                    continue
                break

            for r in ['RH_','RC_','RO_','RN_','RS_','RCl_','RBr_','RI_','RP_','RD_']:
                R_name = r+str(n_at_front)
                if any(R_name in element for element in core_atom):
                    R_name_at_front = R_name
                    break
                else:
                    continue

            for vc in ['VC_3','VC_2','VC_1']:
                if any(vc in element for element in core_atom):
                    VC_name_at_front = vc
                    break

            for v in ['VLC','VLB','VLA']:
                if any(v in element for element in frag_atom):
                    V_name_at_front = v
                    break

            V_name_at_front = 'VLA'



            core_atom = atomname_change_after_link_f(core_atom,X_name_at_front,R_name_at_front,n_at_front)



            c_atom_list = ['VC_1','VC_2','VC_3']
            for c_atom in c_atom_list:
                if any(c_atom in element for element in frag_atom):
                    core_atom,core_coord = remove_VC_after_link_f(core_atom,core_coord,VC_name_at_front,V_name_at_front)
                    frag_atom,frad_coord = remove_V_name_at_front_of_frag_after_link_f_c(frag_atom,frad_coord,V_name_at_front)
                    core_atom.extend(frag_atom)
                    core_coord.extend(frad_coord)
                    core_atom = atomname_change_after_link_f_to_c(core_atom)
                    core_atom = atomname_change_before_link_n(core_atom)
                    break
            else:
                core_atom.extend(frag_atom)
                core_coord.extend(frad_coord)
                core_atom,core_coord = remove_VC_after_link_f(core_atom,core_coord,VC_name_at_front,V_name_at_front)

        elif temp_frag == 'XFO':

            X_name_at_front = 'XC5'

            for r in ['RH_1','RC_1','RO_1','RN_1','RS_1','RCl_1','RBr_1','RI_1','RP_1','RD_1']:
                if any(r in element for element in core_atom):
                    R_name_at_front = r
                    break

            VC_name_at_front = 'VC_1'

            V_name_at_front = 'VLA'

            core_atom.extend(frag_atom)
            core_coord.extend(frad_coord)

            core_atom,core_coord = remove_VC_after_link_f(core_atom,core_coord,VC_name_at_front,V_name_at_front)
            core_atom = atomname_change_after_link_XFO(core_atom,R_name_at_front)

        elif temp_frag == 'sym24_XHD_0' or temp_frag == 'sym24_XHD_1' or temp_frag == 'sym24_XHD_2' or temp_frag == 'sym24_XHD_3' or temp_frag == 'sym24_XHD_4':
            core_atom.extend(frag_atom)
            core_coord.extend(frad_coord)

            core_atom = atomname_change_after_link_XHD(core_atom)
            core_atom,core_coord = remove_atoms_after_link_XHD(core_atom,core_coord)


    new_core = [[name,float(x),float(y),float(z)]for name,[x,y,z] in zip (core_atom,core_coord)]

    return new_core


def ReLU(x, current_score_list,cutoff_score):
    min_val = min(current_score_list)
    slope = 1/(min_val-cutoff_score)
    y = slope*(x-cutoff_score)
    return y

def calc_UCB(selected_core_state,cutoff_score):
    selected_core_score_list = [sublist[1] for sublist in selected_core_state]
    selected_core_selectnum_list = [sublist[2] for sublist in selected_core_state]
    total_selection = 1
    for iii in selected_core_selectnum_list:
        if 10000000 > iii > 1:
            total_selection += iii
            total_selection -= 1
    ubc_list = [ReLU(v, selected_core_score_list,cutoff_score) + (math.sqrt(2) * math.sqrt(np.log(total_selection)/n)) if v <= cutoff_score and n < 10000 else -100 for v,n in zip(selected_core_score_list,selected_core_selectnum_list)]
    max_ubc = max(ubc_list)
    ubc_ind = ubc_list.index(max_ubc)
    return ubc_ind

def red_print(x):
    if type(x) == str:
        print('\033[31m'+x+'\033[0m')
    else:
        print('\033[31m'+str(x)+'\033[0m')

def blue_print(x):
    if type(x) == str:
        print('\033[34m'+x+'\033[0m')
    else:
        print('\033[34m'+str(x)+'\033[0m')

def green_print(x):
    if type(x) == str:
        print('\033[32m'+x+'\033[0m')
    else:
        print('\033[32m'+str(x)+'\033[0m')

def yellow_print(x):
    if type(x) == str:
        print('\033[33m'+x+'\033[0m')
    else:
        print('\033[33m'+str(x)+'\033[0m')

def separate_atom_and_coord(original_coord):
    coord_list = [[float(xs), float(ys), float(zs)] for names, xs, ys, zs in original_coord]
    atom_list = [names for names, xs, ys, zs in original_coord]
    return atom_list,coord_list

def extract_key_elements(atom_list,category,param,ref_frag):
    if category == ['link_c']:
        for x in ['XC3_1','XC3g_1','XC2_1','XC1_1','XCA_1']:
            if any(x in element for element in atom_list):
                central_atom = x
                break
        for r in ['RH_1','RC_1','RO_1','RN_1','RS_1','RCl_1','RBr_1','RI_1','RP_1','RD_1','RF_1']:
            if any(r in element for element in atom_list):
                reference_atom = r
                break
        virtual_atom = 'VC_1'
    if category == ['link_f']:
        if not 'X4' in atom_list and not 'W1' in atom_list:
            central_atom = param[0]
            reference_atom = ref_frag
            virtual_atom = 'VLA'
        elif 'X4' in atom_list:
            central_atom = 'X4'
            reference_atom = 'A4'
            virtual_atom = 'VLA'
        elif 'W1' in atom_list:
            central_atom = 'W1'
            reference_atom = 'Z1'
            virtual_atom = 'Y1'
    return [central_atom,reference_atom,virtual_atom]

def initialize_fragment(atom_list,coord_list,key_elements,angle):
    coord_array = np.array(coord_list)
    ind_central = atom_list.index(key_elements[0])
    ind_reference = atom_list.index(key_elements[1])
    ind_virtual = atom_list.index(key_elements[2])
    coord_central = np.array(coord_list[ind_central])
    coord_array = coord_array-coord_central
    rho,phi,theta = calculate_polar_coords(coord_array,ind_virtual)
    coord_array = np.dot(create_rotation_matrix_xy(phi), coord_array.T).T
    coord_array = np.dot(create_rotation_matrix_xz(theta), coord_array.T).T
    if float(coord_array[ind_virtual,2]) <= 0:
        coord_array = np.dot(create_rotation_matrix_xz(180), coord_array.T).T
    rho,phi,theta = calculate_polar_coords(coord_array,ind_reference)
    coord_array = np.dot(create_rotation_matrix_xy(phi-angle), coord_array.T).T
    coord_list = coord_array.tolist()
    return coord_list

def initialize_core(atom_list,coord_list,key_elements):
    coord_array = np.array(coord_list)
    ind_central = atom_list.index(key_elements[0])
    ind_reference = atom_list.index(key_elements[1])
    ind_virtual = atom_list.index(key_elements[2])
    coord_central = np.array(coord_list[ind_central])
    coord_array -= coord_central
    rho_1,phi_1,theta_1 = calculate_polar_coords(coord_array,ind_virtual)
    coord_array = np.dot(create_rotation_matrix_xy(phi_1), coord_array.T).T
    coord_array = np.dot(create_rotation_matrix_xz(theta_1), coord_array.T).T
    inv_rotation_matrix_4th = np.linalg.inv(create_rotation_matrix_xy(phi_1))
    inv_rotation_matrix_3nd = np.linalg.inv(create_rotation_matrix_xz(theta_1))
    if float(coord_array[ind_virtual,2]) >= 0:
        coord_array = np.dot(create_rotation_matrix_xz(180), coord_array.T).T
        inv_rotation_matrix_2nd = create_rotation_matrix_xz(180)
    else:
        inv_rotation_matrix_2nd = create_rotation_matrix_xz(0)
    rho_2,phi_2,theta_2 = calculate_polar_coords(coord_array,ind_reference)
    coord_array = np.dot(create_rotation_matrix_xy(phi_2), coord_array.T).T
    coord_list = coord_array.tolist()
    inv_rotation_matrix_1st = np.linalg.inv(create_rotation_matrix_xy(phi_2))
    return inv_rotation_matrix_4th, inv_rotation_matrix_3nd, inv_rotation_matrix_2nd, inv_rotation_matrix_1st, coord_central

def calculate_polar_coords(coord_array,desired_index):
    rho = np.linalg.norm(coord_array[desired_index])
    phi = math.degrees(math.atan2(coord_array[desired_index][1], coord_array[desired_index][0]))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        theta = math.degrees(math.acos(coord_array[desired_index][2] / rho))
    return rho,phi,theta

def create_rotation_matrix_xy(phi):
    minus_phi_rad = np.radians(-phi)
    rotation_matrix = np.array([[np.cos(minus_phi_rad), -np.sin(minus_phi_rad), 0],
                                  [np.sin(minus_phi_rad), np.cos(minus_phi_rad), 0],
                                  [0, 0, 1]])
    return rotation_matrix

def create_rotation_matrix_xz(theta):
    minus_theta_rad = np.radians(-theta)
    rotation_matrix = np.array([[np.cos(minus_theta_rad), 0, np.sin(minus_theta_rad)],
                                  [0, 1, 0],
                                  [-np.sin(minus_theta_rad), 0, np.cos(minus_theta_rad)]])
    return rotation_matrix

def create_structure_debug(Core_mol,coord_array,atom_list,xyz_list):
    pass

def create_structure_from_strategy(each_state,ubc_ind,cycle,max_branch_core,max_branch_frag):
    for frag in each_state[ubc_ind][0]:

        if type(frag) ==  str:
            current_cycle_int = 0
            current_branch_core_int = 0
            current_branch_frag_int = 0

            transforming_core_point_int,S,Q,P,total_charge_int,transform_log_main_int,transform_log_side_int,c_for_next_int,param_for_next_int,added_frag_int,ref_frag_list_next_int = import_core(frag)
            core_int = Core_mol(transforming_core_point_int,current_cycle_int,cycle,current_branch_core_int,current_branch_frag_int,max_branch_core,max_branch_frag,ref_frag_list_next_int,param_for_next_int,total_charge_int,transform_log_main_int,transform_log_side_int,c_for_next_int)
        else:
            temp_fragment = single_Fragment(frag)

            all_atom_int = temp_fragment.atom
            all_coord_int = temp_fragment.coord
            all_charge_int = temp_fragment.charge
            all_each_frag_int = temp_fragment.each_frag
            all_len_list_int = temp_fragment.len_list
            all_central_atom_int = temp_fragment.central_atom
            all_category_int = temp_fragment.category
            all_param_int = temp_fragment.param
            all_frag_int = temp_fragment.frag_name
            all_add_branch_core_int = temp_fragment.add_branch_core
            all_add_branch_frag_int = temp_fragment.add_branch_frag
            all_ref_frag_list_int = temp_fragment.ref_frag_list
            frag_c_for_next = temp_fragment.c_for_next

            param_for_next_int

            result_int = link_result_no_calc_interaction(core_int,all_atom_int,all_coord_int,all_charge_int,all_central_atom_int,None,None,all_category_int,all_param_int,all_len_list_int,all_frag_int,all_ref_frag_list_int,frag_c_for_next)
            next_frag_atom_int = result_int.atom
            next_frag_coord_int = result_int.coord
            EOF_value_int = result_int.eof_val[0]
            category_int = result_int.temp_category[0]
            temp_frag_int = result_int.temp_frag[0]
            ref_frag_list_next_int = result_int.ref_frag_list_next_list[0]
            param_for_next_int = result_int.param_for_next[0]
            total_charge_int = result_int.total_charge[0]
            transform_log_main_int = result_int.transform_log_main_list[0]
            transform_log_side_int = result_int.transform_log_side_list[0]
            c_for_next_int = result_int.c_for_next_list[0]
            transforming_core_point_int = create_new_core_structure(core_int,next_frag_atom_int,next_frag_coord_int,category_int,temp_frag_int,ref_frag_list_next_int)

            if EOF_value_int == False:

                core_int = Core_mol(transforming_core_point_int,current_cycle_int,cycle,current_branch_core_int,current_branch_frag_int,max_branch_core,max_branch_frag,ref_frag_list_next_int,param_for_next_int,total_charge_int,transform_log_main_int,transform_log_side_int,c_for_next_int)

    return transforming_core_point_int

#bonded param

distance_dict = {
                'D,XC3_':1.09,
                'D,XC3g_':1.09,
                'D,XCA_':1.08,
                'D,XC2_':1.08,
                'D,XC1_':1.07,
                'D,YCA':1.09,
                'D,YO':0.97,
                'D,YNA':1.01,
                'D,YN3':1.01,
                'D,YN2':1.01,
                'W1,X4':1.30,
                'X4,W1':1.30,
                'XC5,X4':1.53,
                'X4,XC5':1.53,
                'H,XC3_':1.09,
                'H,XC3g_':1.09,
                'H,XCA_':1.08,
                'H,XC2_':1.08,
                'H,XC1_':1.07,
                'H,YCA':1.09,
                'H,YO':0.97,
                'H,YNA':1.01,
                'H,YN3':1.01,
                'H,YN2':1.01,
                'C,XC3_':1.54,
                'C,XC3g_':1.54,
                'C,XCA_':1.52,
                'C,XC2_':1.52,
                'C,XC1_':1.47,
                'C,YCA':1.53,
                'C,YO':1.42,
                'C,YNA':1.53,
                'C,YN3':1.53,
                'C,YN2':1.53,
                'C,YS':1.82,
                'C,YSA':1.87,
                'XC3_,XC3_':1.54,
                'XC3_,XC3g_':1.54,
                'XC3_,XCA_':1.52,
                'XC3_,XC2_':1.52,
                'XC3_,XC1_':1.47,
                'XC3_,YCA':1.53,
                'XC3_,YCO':1.53,
                'XC3_,YO':1.42,
                'XC3_,YNA':1.47,
                'XC3_,YN3':1.53,
                'XC3_,YN2':1.53,
                'XC3_,YS':1.82,
                'XC3_,YSA':1.87,
                'XC3g_,XC3_':1.54,
                'XC3g_,XC3g_':1.54,
                'XC3g_,XCA_':1.52,
                'XC3g_,XC2_':1.52,
                'XC3g_,XC1_':1.47,
                'XC3g_,YCA':1.53,
                'XC3g_,YCO':1.53,
                'XC3g_,YO':1.42,
                'XC3g_,YNA':1.47,
                'XC3g_,YN3':1.53,
                'XC3g_,YN2':1.53,
                'XC3g_,YS':1.82,
                'XC3g_,YSA':1.87,
                'XCA_,XC3_':1.52,
                'XCA_,XC3g_':1.52,
                'XCA_,XCA_':1.50,
                'XCA_,XC2_':1.50,
                'XCA_,XC1_':1.46,
                'XCA_,YCA':1.53,
                'XCA_,YCO':1.53,
                'XCA_,YO':1.38,
                'XCA_,YNA':1.46,
                'XCA_,YN2':1.40,
                'XCA_,YS':1.76,
                'XCA_,YSA':1.79,
                'XC2_,XC3_':1.52,
                'XC2_,XC3g_':1.52,
                'XC2_,XCA_':1.50,
                'XC2_,XC2_':1.50,
                'XC2_,XC1_':1.46,
                'XC2_,YCA':1.53,
                'XC2_,YCO':1.53,
                'XC2_,YSA':1.79,
                'XC1_,XC3_':1.47,
                'XC1_,XC3g_':1.47,
                'XC1_,XCA_':1.46,
                'XC1_,XC2_':1.46,
                'XC1_,XC1_':1.36,
                'XC1_,YCA':1.45,
                'XC1_,YCO':1.45,
                'XC1_,YSA':1.72,
                'YCA,XC3_':1.53,
                'YCA,XC3g_':1.53,
                'YCA,XCA_':1.50,
                'YCA,XC2_':1.50,
                'YCA,XC1_':1.45,
                'YCO,XC3_':1.53,
                'YCO,XC3g_':1.53,
                'YCO,XCA_':1.50,
                'YCO,XC2_':1.50,
                'YCO,XC1_':1.45,
                'YO,XC3_':1.42,
                'YO,XC3g_':1.42,
                'YO,XCA_':1.43,
                'YNA,XC3_':1.47,
                'YNA,XC3g_':1.47,
                'YNA,XCA_':1.48,
                'YN3,XC3_':1.53,
                'YN3,XC3g_':1.53,
                'YN2,XCA_':1.40,
                'YN2,XC3_':1.53,
                'YN2,XC3g_':1.53,
                'YS,XC3_':1.82,
                'YS,XC3g_':1.82,
                'YS,XCA_':1.76,
                'YSA,XC3_':1.87,
                'YSA,XC3g_':1.87,
                'YSA,XCA_':1.79,
                'YSA,XC2_':1.79,
                'YSA,XC1_':1.72,
                'XC0,XC3_':1.52,
                'XC0,XC3g_':1.52,
                'XC0,XCA_':1.50,
                'XC0,XC2_':1.50,
                'XC0,XC1_':1.46,
                'CAD,XC3_':1.52,
                'CAD,XC3g_':1.52,
                'CAD,XCA_':1.50,
                'CAD,XC2_':1.50,
                'CAD,XC1_':1.46,
                'CAD,YCA':1.53,
                'CAD,YCO':1.53,
                'CAD,YO':1.38,
                'CAD,YNA':1.46,
                'CAD,YN2':1.40,
                'CAD,YS':1.76,
                'CAD,YSA':1.79,
                'CAD,F':1.35,
                'CAD,Cl':1.72,
                'CAD,Br':1.90,
                'CAD,H':1.08,
                'CKD,XC3_':1.52,
                'CKD,XC3g_':1.52,
                'CKD,XCA_':1.50,
                'CKD,XC2_':1.50,
                'CKD,XC1_':1.46,
                'CKD,YCA':1.53,
                'CKD,YCO':1.53,
                'CKD,YSA':1.79,
                'CKD,H':1.08,
                'CRD,XC3_':1.54,
                'CRD,XC3g_':1.54,
                'CRD,XCA_':1.52,
                'CRD,XC2_':1.52,
                'CRD,XC1_':1.47,
                'CRD,YCA':1.53,
                'CRD,YCO':1.53,
                'CRD,YO':1.42,
                'CRD,YNA':1.47,
                'CRD,YN3':1.53,
                'CRD,YS':1.82,
                'CRD,YSA':1.87,
                'CRD,F':1.33,
                'CRD,H':1.09,
                'HET,XC3_':1.47,
                'HET,XC3g_':1.47,
                'HET,XCA_':1.48,
                'HET,D':1.01,
                }

not_permitted_comb_dict = {
                          'YCA,XC3g_,YCA':1,
                          'YCA,XC3g_,YCO':1,
                          'YCO,XC3g_,YCA':1,
                          'YCA,XC3g_,YSA':1,
                          'YNA,XC3g_,YN3':1,
                          'YNA,XC3g_,YN2':1,
                          'YNA,XC3g_,YNA':1,
                          'YNA,XC3g_,YO':1,
                          'YNA,XC3g_,YS':1,
                          'YNA,XC3g_,YSA':1,
                          'YO,XC3g_,YN3':1,
                          'YO,XC3g_,YN2':1,
                          'YO,XC3g_,YNA':1,
                          'YO,XC3g_,YO':1,
                          'YO,XC3g_,YS':1,
                          'YO,XC3g_,YSA':1,
                          'YS,XC3g_,YN3':1,
                          'YS,XC3g_,YN2':1,
                          'YS,XC3g_,YNA':1,
                          'YS,XC3g_,YO':1,
                          'YS,XC3g_,YS':1,
                          'YS,XC3g_,YSA':1,
                          'YSA,XC3g_,YN3':1,
                          'YSA,XC3g_,YN2':1,
                          'YSA,XC3g_,YNA':1,
                          'YSA,XC3g_,YO':1,
                          'YSA,XC3g_,YS':1,
                          'YSA,XC3g_,YSA':1,
                          'YN3,XC3g_,YN3':1,
                          'YN3,XC3g_,YN2':1,
                          'YN3,XC3g_,YNA':1,
                          'YN3,XC3g_,YO':1,
                          'YN3,XC3g_,YS':1,
                          'YN3,XC3g_,YSA':1,
                          'YN2,XC3g_,YN3':1,
                          'YN2,XC3g_,YN2':1,
                          'YN2,XC3g_,YNA':1,
                          'YN2,XC3g_,YO':1,
                          'YN2,XC3g_,YS':1,
                          'YN2,XC3g_,YSA':1,
                          'YSA,XC3g_,YN3':1,
                          'YSA,XC3g_,YN2':1,
                          'YSA,XC3g_,YNA':1,
                          'YSA,XC3g_,YO':1,
                          'YSA,XC3g_,YS':1,
                          'YSA,XC3g_,YSA':1,
                          'XCA_,YO,XCA_':1,
                          'XCA_,YO,XC2_':1,
                          'XCA_,YO,XC1_':1,
                          'XC2_,YO,XCA_':1,
                          'XC1_,YO,XCA_':1,
                          'XC3_,YO,XC2_':1,
                          'XC3_,YO,XC1_':1,
                          'XC3g_,YO,XC2_':1,
                          'XC3g_,YO,XC1_':1,
                          'XC3g_,YO,XC3g_':1,
                          'XC3_,YN2,XC3_':1,
                          'XC3g_,YN2,XC3_':1,
                          'XC3_,YN2,XC3g_':1,
                          'XC3g_,YN2,XC3g_':1,
                          }
