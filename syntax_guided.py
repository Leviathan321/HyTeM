#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import random 
import numpy as np
import copy
import math
import os
import pickle
import sys
import time
import itertools
import shutil

import tree
import SMT
import syntax_eval as se
import create_grammar as gramm
import binary_search

def print_formula(quantifiers, formula):
    
    m = len(quantifiers)
    
    f = ''
    
    for i in range(m):
        if quantifiers[i]==0: f += f'E pi[{i}] '
        elif quantifiers[i]==1: f += f'FA pi[{i}] '
    
    f = f + ': ' + formula
    return f

def parse_gramm(grammar):
    
    '''The function transforms the grammar (either of the quntifiers or of the quantifier-free formula structure)
    given as a string with different levels in the format needed by the algorithm.'''
    
    #Different levels divided by the new line
    levels = grammar.split('\n')
    grammar_structure = []
    
    for level in levels:
        line = []
        
        if '='  in level: #Remove first and last line because they are empty
            
            right = level.split('==')
            if '|' not in right[1]: #right[1] is the right  part of the equation
                line = [right[1]]
            else: 
                item = right[1].split('|')
                for i in range(len(item)):
                    if '|' in item[i]: item[i] = item[i].replace('|', '')
                    if 'forall' in item[i]: item[i] = item[i].replace('forall', 'A')
                    if 'exists' in item[i]: item[i] = item[i].replace('exists', 'E')
                    line.append(item[i])
            grammar_structure.append(line)
    return grammar_structure 


def main(inputs):
    
    ##Initialization
    ##Set all the inputs
    sed = inputs.seed
    traces = inputs.traces
    
    min_nb_quantifiers = inputs.numb_quantifiers[0]
    max_nb_quantifiers = inputs.numb_quantifiers[1]
    
    grammar_quantifiers = parse_gramm(inputs.grammar_quantifiers)
    grammar_structure = parse_gramm(inputs.grammar_structure)
    initial_grammar_predicate = inputs.grammar_predicate
    grammar_predicate = copy.deepcopy(initial_grammar_predicate)
    
    learning_stl = inputs.learning_stl 
    if learning_stl is True:  par_bounds = inputs.par_bounds
    fraction = inputs.fraction
    
    kind_score = inputs.kind_score
    
    bool_interface_aware = inputs.bool_interface_aware
    input_variables = inputs.input_variables
    output_variables = inputs.output_variables
        
    bool_temporal_operators = inputs.bool_temporal_operators
    bool_mutation_quantifiers = inputs.bool_mutation_quantifiers
    bool_different_variables = inputs.bool_different_variables
    bool_same_number = inputs.bool_same_number 
    second_variable = inputs.second_variable
    store_enumeration_formulas = inputs.store_enumeration_formulas
    store_formulas = inputs.store_formulas
    
    name_set_mined_formulas = inputs.name_set_mined_formulas
    
    set_violated_formulas = set()
    set_satisfied_formulas = set()
    set_database_parametric = set()
    
    ##Start
    
    random.seed(sed)
    # numb_traces = len(traces)
    n = len(traces[0]) # number of variables    
    
    
    #Maximum number of Iterations for the same formula
    n_target = inputs.n_target #50 # number of target formulas
    nb_iterations_outer = n_target * 2
    numb_it_sat_max = 500 # trials on the same starting formula to find a satisfied one
    numb_start_cond_max = 50 # trials on the same inizialization (l , m)
    numb_rej_max = 50
    
    if not os.path.exists('grammar'): os.mkdir('grammar')
    if not store_enumeration_formulas: shutil.rmtree('grammar') 
    
    if store_formulas and not os.path.exists('Learned_formulas'): os.mkdir('Learned_formulas')
    path_to_save = f'Learned_formulas/{name_set_mined_formulas}learned_formula_seed{sed}/' 
    if store_formulas and not os.path.exists(path_to_save): os.mkdir(path_to_save)
    count_added_formula = 0
    
    #Storage of formulas:
    database = []
    iter_outer_loop = 0
    
    time_start = time.time()
    
    #Loop to find a set of formulas (LEVEL 1)
    while True:
        iter_outer_loop += 1
        print('Formulas learned so far:', len(database))
        
        #Already collected n_target formulas --> END of the algorithm
        if len(database) >= n_target: break
    
        #Stop after "too many" iterations
        if iter_outer_loop > nb_iterations_outer: break
    
        
        ##INITIALIZATION
        #Sample the quantifiers
        quantifiers , quantifiers_levels,  spec_level_start = tree.sample_quantifiers(grammar_quantifiers, min_nb_quantifiers, max_nb_quantifiers)
        m = len(quantifiers)
        print(f'Number of quantifiers = {m}')
        
        if learning_stl == False:
        #Replace 'alltuples' in the grammar_predicate
            for index_gr in range(len(grammar_predicate)):
                if 'all tuples' in initial_grammar_predicate[index_gr][0][-1]: grammar_predicate[index_gr][0][-1] = [ item for item in itertools.permutations(range(m), min(m,2))]
        elif learning_stl == True:
        #Allows only predicates that do not involve trace variables greater than the number of quantifiers
            grammar_predicate = []
            for item_gr in initial_grammar_predicate:
                if f'pi[{m}]' not in item_gr and f'pi[{m+1}]' not in item_gr : grammar_predicate.append(item_gr)
                
        
        # Tuples of traces to examine in effient_monitoring
        n_subset = None #1/100 * int(math.factorial(numb_traces)/math.factorial(numb_traces - m)) #None #10000
          
        
        #Set the length of formula structure
        if inputs.l is not None: l = inputs.l
        else:   
            max_length = inputs.max_length
            l = random.randint(2*math.ceil(m/2)-1, max_length) # length of the formula
            if kind_score == 'fitness' and l ==1: l=2
        print(f'Target length = {l}')
        
        
        #Counter on the iteration in the following loop (trials with same starting conditions)    
        numb_start_cond = 0
        
        #Loop to find a SINGLE SATISFIED FORMULA (LEVEL 2)
        while True:
            
            numb_start_cond += 1
            if numb_start_cond > numb_start_cond_max:  break
            
            ##REJECTION SAMPLING FOR TRIVIAL OR NOT WELL-FORMED FORMULA
            #not well-formed means that not all trace variables named by the quantifiers appear in the internal formula
            num_rej = 0# number of rejections. After numb_rej_max accepts a double negation
            
            ##Loop to find a valid starting formula (LEVEL 3)
            # print('Start rejection sampling')
            while True:  #REJECTION SAMPLING 
                num_rej += 1
                
                
                #If formulas have not been listed, do so
                if not os.path.exists('grammar'): os.mkdir('grammar')
                if not os.path.exists(f'grammar/spec_level{spec_level_start}_numbquant{m}'): os.mkdir(f'grammar/spec_level{spec_level_start}_numbquant{m}')
                if not os.path.exists(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}'): 
                    os.mkdir(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}')
                    listing_nodes = []
                    listing_nodes.append(tree.BinaryTreeNode('start'))
                    
                    listing_nodes =  gramm.build_formulas(l, grammar_structure, spec_level_start, listing_nodes, 0, 'left')
                    gramm.transform_listing_to_nodes(listing_nodes, l, spec_level_start,m)
                    
                    
                    
                'Create the formula'
                if len(os.listdir(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}')) != 0:
                    list_formulas = os.listdir(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}')
                    list_formulas.sort()
                    sampled_f = random.sample(list_formulas,1)[0]
                    with open(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}/{sampled_f}', 'rb') as f:
                      nodes = pickle.load(f)
                    
                    if learning_stl == True: nodes = tree.sample_predicate_STL(grammar_predicate, nodes, bool_different_variables)
                    else: nodes = tree.sample_predicate(grammar_predicate, nodes, bool_different_variables,  bool_same_number, second_variable)
                    
                    formula = tree.tree_to_formula(nodes)
                
                    if inputs.to_be_present is not None : to_be_present = copy.deepcopy(inputs.to_be_present)
                    else: to_be_present = [ [f'pi[{j}]' for j in range(m) ] , m]
                        
                    # if second_variable == 'true' : to_be_present = [ ['pi[0][0]', 'pi[0][1]'], 2]  
                    # elif second_variable == 'variable': to_be_present = [ [f'pi[{j}]' for j in range(m)], m ] #Check if all pi[j] for j=0,..,m-1 appear in formula
                    # else: to_be_present = [[f'pi[0][{agent}] == {action}' for agent in range(n) for action in second_variable ], 2 ]
                       
                    if sum([to_be_present[0][i] in formula for i in range(len(to_be_present[0]))])>=to_be_present[1]:
                        
                        #Check if there is a double negation
                        if not tree.check_double_negation(nodes) or num_rej > numb_rej_max :
                            
                            #If there are NOT temporal operators we can use z3 to check whether the formula is trivially satisfied
                            if bool_temporal_operators == False:
                                #Exit if the formula is not-trivial
                                if not(SMT.evaluate_triviality(nodes, m ,n)):break
                                
                            elif bool_temporal_operators == True and learning_stl == False:
                                break
                            
                            else: #learning_stl True
                                #Check monotonicity for stl, otherwise rejection sampling
                                mono = binary_search.check_monotonicity(nodes)
                                if mono != False: # ok monotonic
                                    break
                            
                    #Too many rejected -> sample new length for the formula        
                    if num_rej > numb_rej_max + 1:
                        #Change length
                        l = random.randint(2*math.ceil(m/2)-1,max_length) # sample the length of the formula
                        print(f'Need to change the length. New length = {l}')
                        num_rej = 0
            
                elif l < 20: 
                    print(f'It was not possible building a formula with length {l}: I will increment length by one. New length = {l+1}\n')
                    l += 1
                
                else:
                    print('No formula found. Start with a smaller value of length or change the grammar\n')
                    sys.exit()

            # print('End rejection sampling')

            #In case of STL, replace the parameters with the values that are MOST likely to be satisfied
            if learning_stl == True:
                index_mono = 0
                for index_epsilon in range(10):
                    if f'epsilon{index_epsilon}-' in formula:
                        if mono[index_mono] == '+' : concrete_value = par_bounds[index_epsilon][1]#the maximum
                        elif mono[index_mono] == '-': concrete_value = par_bounds[index_epsilon][0] #the minimum
                        formula = formula.replace(f'epsilon{index_epsilon}-', f'{concrete_value}')
                        index_mono += 1
            
            #Avoid recomputing violated formulas
            if set([print_formula(quantifiers, formula)]).issubset(set_violated_formulas) and ( kind_score == 'efficient_monitoring' or kind_score == 'efficient_monitoring_stl'):
                    cost = -1
            #Avoid recomputing satisfied formulas   
            elif set([print_formula(quantifiers, formula)]).issubset(set_satisfied_formulas) and ( kind_score == 'efficient_monitoring' or kind_score == 'efficient_monitoring_stl'):
                    cost = 1
                    
            else: 
                # 'Cost'
                if kind_score == 'fitness': cost = se.compute_cost(traces, formula, quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables)
                elif kind_score == 'stl_robustness' : cost = se.compute_robustness_STL(traces, formula, quantifiers)
                elif kind_score == 'efficient_monitoring': cost = se.efficient_monitoring(traces, formula, quantifiers, bool_temporal_operators, n_subset,  bool_interface_aware, input_variables, output_variables)
                elif kind_score == 'efficient_monitoring_stl': cost = se.monitor_hyperSTL_correctness(traces, formula, quantifiers)
                
                if cost == 'runout' or cost < 0:  set_violated_formulas = set_violated_formulas.union(set([print_formula(quantifiers, formula)]))
                elif cost > 0:  set_satisfied_formulas = set_satisfied_formulas.union(set([print_formula(quantifiers, formula)]))
            
            #Collect number of iterations to exit the following loop with unsuccess
            numb_it_sat = 0
            
            #Loop to find a satisfied formula (LEVEL 3)
            # print('Start mutating')
            while True:
                numb_it_sat += 1
                
                #A satisfied formula is found
                if cost != 'runout' and cost > 0:
                    
                     #In case of parametric STL: --> find the tightest satisfaction
                     if learning_stl == True:
                         formula = tree.tree_to_formula(nodes)#parametric formula
                         
                         if not set([print_formula(quantifiers, formula)]).issubset(set_database_parametric):
                             is_new_formula = True
                             set_database_parametric = set_database_parametric.union(set([print_formula(quantifiers, formula)]))
                             
                             #Pass only the parameter bounds of the parameters in the formula
                             par_bounds_input = []
                             for index_epsilon in range(10):
                               if f'epsilon{index_epsilon}-' in formula: par_bounds_input.append(par_bounds[index_epsilon])
                             
                             area  = math.prod([abs(par_bounds_input[i][1]-par_bounds_input[i][0]) for i in range(len(par_bounds_input)) ])
                             threshold =  area/fraction  
                             
                             #Mine parameters
                             mined_parameter = binary_search.compute_monotonic_parameters(quantifiers, formula, traces, threshold, par_bounds_input, mono)
                       
                             #Replace the formula with concrete values
                             index_mono = 0
                             for index_epsilon in range(10):
                               if f'epsilon{index_epsilon}-' in formula:
                                  formula = formula.replace(f'epsilon{index_epsilon}-', f'{mined_parameter[index_mono]}')
                                  for node in nodes: node.data = node.data.replace(f'epsilon{index_epsilon}-', f'{mined_parameter[index_mono]}')
                                  index_mono += 1
                                  
                         else:  is_new_formula = False
                     
                     elif learning_stl == False:
                         #Check if the formula is already in the database:
                         is_new_formula = True
                         for i_database in range(len(database)):
                             
                             if tree.tree_to_formula(database[i_database].nodes) == tree.tree_to_formula(nodes): 
                                 is_new_formula = False
                                 break
                             
                     
                    
                     
                     if is_new_formula:
                         time_end = time.time()  
                         #Store the new learned formula
                         formula_to_be_stored = gramm.StoreFormula(nodes, quantifiers, cost, time_end-time_start,l)
                         database.append(formula_to_be_stored)
                         print('Found a new satisfied formula!\n')
                         numb_start_cond = numb_start_cond_max + 1
                         print(f'Formula: {print_formula(formula_to_be_stored.quantifiers, tree.NATURAL_tree_to_formula(formula_to_be_stored.nodes))};\ncost = {formula_to_be_stored.cost}, length={formula_to_be_stored.length}, time = {formula_to_be_stored.time}\n\n')
                         count_added_formula += 1
                         
                         #Store only if explicitely asked
                         if store_formulas is True:
                             with open(f'{path_to_save}formula_to_be_stored{count_added_formula}.obj', 'wb') as f:
                                     pickle.dump(formula_to_be_stored, f)
                         
                         time_start = time.time()  
                             
                     break
                
                
                #Stop after numb_it_sat_max unsuccessful iterations    
                if numb_it_sat > numb_it_sat_max:  break
                
                if bool_mutation_quantifiers: rand_number =  random.uniform(0,1)
                else: rand_number = 0.7
                
                
                
                # With probability 0.5, the changes are applied to the formula structure
                if rand_number >= 0.5: 
                    
                    ##REJECTION SAMPLING FOR TRIVIAL OR NOT WELL-FORMED FORMULA
                    num_rej = 0
                    #Loop to find a valid new formula (LEVEL 4)
                    while True: 
                        num_rej += 1
                        #Select the nodes that will be deleted
                        indices_old_node, length_node = tree.select_node(nodes)
                        
                        'Change the formula'
                        spec_level_change = nodes[indices_old_node[0]].level_grammar
                        if not os.path.exists(f'grammar/spec_level{spec_level_change}_numbquant{m}'): os.mkdir(f'grammar/spec_level{spec_level_change}_numbquant{m}')
                        if not os.path.exists(f'grammar/spec_level{spec_level_change}_numbquant{m}/length_{length_node}'): 
                            os.mkdir(f'grammar/spec_level{spec_level_change}_numbquant{m}/length_{length_node}')
                            listing_nodes = []
                            listing_nodes.append(tree.BinaryTreeNode('start'))
                            listing_nodes =  gramm.build_formulas(length_node, grammar_structure, spec_level_change, listing_nodes, 0, 'left')
                            gramm.transform_listing_to_nodes(listing_nodes, length_node , spec_level_change, m)
                        list_changed_formulas = os.listdir(f'grammar/spec_level{spec_level_change}_numbquant{m}/length_{length_node}')
                        list_changed_formulas.sort()
                        sampled_f = random.sample(list_changed_formulas,1)[0]
                        with open(f'grammar/spec_level{spec_level_change}_numbquant{m}/length_{length_node}/{sampled_f}', 'rb') as f:
                            new_node = pickle.load(f)
                        
                        if learning_stl == True: new_node = tree.sample_predicate_STL(grammar_predicate, new_node, bool_different_variables)
                        else: new_node = tree.sample_predicate(grammar_predicate, new_node, bool_different_variables, bool_same_number,second_variable)
                        
                        #Embed the new nodes inside the tree of the original formula
                        new_nodes = tree.change_node(nodes, indices_old_node, new_node)
                        
                        #Transform the new tree into a string
                        new_formula = tree.tree_to_formula(new_nodes)
                        
                        if inputs.to_be_present is not None : to_be_present = copy.deepcopy(inputs.to_be_present)
                        else: to_be_present = [ [f'pi[{j}]' for j in range(m) ] , m]
                        # if second_variable == 'true' : to_be_present = [ ['pi[0][0]' , 'pi[0][1]'], 2]  
                        # elif second_variable == 'variable': to_be_present = [ [f'pi[{j}]' for j in range(m) ] , m]  #Check if all pi[j] for j=0,..,m-1 appear in the new formula
                        # else: to_be_present = [[f'pi[0][{agent}] == {action}' for agent in range(n) for action in second_variable ], 2 ]
                        
                        if sum([to_be_present[0][i] in new_formula for i in range(len(to_be_present[0]))])>=to_be_present[1]: #2<->m
                            # Check if there is a double negation
                            if not tree.check_double_negation(nodes) or num_rej > numb_rej_max :
                               
                                if bool_temporal_operators == False:
                                #Exit if the new formula is not-trivial
                                    if not(SMT.evaluate_triviality(new_nodes, m ,n)):break
                                elif bool_temporal_operators == True and learning_stl == False:
                                    break
                                else:#Before THE EXIT: CHECK MONOTONICITY FOR STL
                                    mono = binary_search.check_monotonicity(new_nodes)
                                    if mono != False: # ok monotonic
                                        break
                               
                    ## Various PRINTS to understand what is happening
                    # formula_to_be_replaced = tree.tree_to_formula([nodes[i] for i in indices_old_node])    
                    # f_to_be_replaced = print_formula(quantifiers, formula_to_be_replaced)
                    # tree.print_tree(new_nodes)
                    new_quantifiers = quantifiers
                    new_quantifiers_levels = quantifiers_levels
                    
                # With probability 0.5, the changes are applied to the formula's quantifiers
                else: 
                    new_quantifiers, new_quantifiers_levels = tree.sample_new_quantifiers(quantifiers, quantifiers_levels, grammar_quantifiers)
                    new_formula = formula
                    new_nodes = copy.deepcopy(nodes)
                 
                    
                ##In case of STL, replace the parameters with the values that are MOST likely to be satisfied
                if learning_stl == True:
                    index_mono = 0
                
                    for index_epsilon in range(10):
                        
                        if f'epsilon{index_epsilon}-' in new_formula:
                            
                            if mono[index_mono] == '+' : concrete_value = par_bounds[index_epsilon][1]#the maximum
                            elif mono[index_mono] == '-': concrete_value = par_bounds[index_epsilon][0] #the minimum
                            new_formula = new_formula.replace(f'epsilon{index_epsilon}-', f'{concrete_value}')
                            index_mono += 1
                
                
                #If the formula has already been evaluated and it is not satisfied do not evaluate it again
                if set([print_formula(new_quantifiers, new_formula)]).issubset(set_violated_formulas) and ( kind_score == 'efficient_monitoring' or kind_score == 'efficient_monitoring_stl'):
                    new_cost = -1
                
                elif set([print_formula(new_quantifiers, new_formula)]).issubset(set_satisfied_formulas) and ( kind_score == 'efficient_monitoring' or kind_score == 'efficient_monitoring_stl'):
                    new_cost = 1
                    
                else:
                    #Compute new cost
                    #Compute cost of the new function 
                    if kind_score == 'fitness': new_cost = se.compute_cost(traces, new_formula, new_quantifiers, bool_temporal_operators,  bool_interface_aware, input_variables, output_variables)
                    elif kind_score == 'stl_robustness': new_cost = se.compute_robustness_STL(traces, new_formula, new_quantifiers)
                    elif kind_score == 'efficient_monitoring': new_cost = se.efficient_monitoring(traces, new_formula, new_quantifiers, bool_temporal_operators, n_subset,  bool_interface_aware, input_variables, output_variables)
                    elif kind_score == 'efficient_monitoring_stl': new_cost = se.monitor_hyperSTL_correctness(traces, new_formula, new_quantifiers)
                    # print(cost, new_cost)
                    if new_cost == 'runout' or new_cost < 0:  set_violated_formulas = set_violated_formulas.union(set([print_formula(new_quantifiers, new_formula)]))
                    elif new_cost > 0: set_satisfied_formulas = set_satisfied_formulas.union(set([print_formula(new_quantifiers, new_formula)]))
                #hastings acceptance ratio --> compute the score from the cost function
                
                
                
                ## Attention because when cost is less than -10 000 , the exp becomes equal to 0
                if new_cost != 'runout':
                    if new_cost > cost: acceptance_ratio = 1
                    elif new_cost == cost: acceptance_ratio = 0.5
                    else: acceptance_ratio = min(1, np.exp(new_cost * 0.5)/np.exp(cost * 0.5))#this cost makes sense when current cost is negative
                    
                    rand_number = random.uniform(0,1)
                    # Apply changes to the formula
                    if rand_number <= acceptance_ratio:
                        formula = new_formula
                        cost = new_cost
                        quantifiers = new_quantifiers
                        nodes = copy.deepcopy(new_nodes)
                        # print(f'\nChange formula into {print_formula(quantifiers, tree.NATURAL_tree_to_formula(nodes))}')
            # print('End mutating')

    #Order database according to the cost value
    # database.sort(key=lambda x: x.cost)
    if not store_enumeration_formulas: shutil.rmtree('grammar')
    #PRINT formula in the learned database:
    for i in range(len(database)): print(f'\nFormula {i+1}: {print_formula(database[i].quantifiers, tree.NATURAL_tree_to_formula(database[i].nodes))};\ncost = {database[i].cost}, length={database[i].length}, time = {database[i].time}\n\n')

    return database

def sample_violated_formulas_stl(inputs):
    
    
        '''The function samples violated formulas for plot stl monitoring '''
        
        grammar_quantifiers = parse_gramm(inputs.grammar_quantifiers)
        grammar_structure = parse_gramm(inputs.grammar_structure)
        
        min_nb_quantifiers = inputs.numb_quantifiers[0]
        max_nb_quantifiers = inputs.numb_quantifiers[1]
    
    
        random.seed(inputs.seed)
        numb_rej_max = 50
        numb_start_cond_max = 50
        database = []
        database_set = set()
        time_start = time.time()
        
        while True:
            #Already collected n_target formulas --> END of the algorithm
            if len(database) >= inputs.n_target: break
            
            ##INITIALIZATION
            #Sample the quantifiers
            quantifiers , quantifiers_levels,  spec_level_start = tree.sample_quantifiers(grammar_quantifiers, min_nb_quantifiers, max_nb_quantifiers)
            m = len(quantifiers)
            print(f'Number of quantifiers = {m}')
            
            initial_grammar_predicate = inputs.grammar_predicate
            grammar_predicate = []
            for item_gr in initial_grammar_predicate:
                if f'pi[{m}]' not in item_gr and f'pi[{m+1}]' not in item_gr : grammar_predicate.append(item_gr)
                    
            max_length = inputs.max_length
            l = random.randint(2*math.ceil(m/2)-1, max_length) # length of the formula
            print(f'Target length = {l}')
            
            numb_start_cond = 0
            #Loop to find a SINGLE VIOLATED FORMULA (LEVEL 2)
            while True:
                
                numb_start_cond += 1
                if numb_start_cond > numb_start_cond_max:  break
                
                ##REJECTION SAMPLING FOR TRIVIAL OR NOT WELL-FORMED FORMULA
                #not well-formed means that not all trace variables named by the quantifiers appear in the internal formula
                num_rej = 0# number of rejections. After numb_rej_max accepts a double negation
                
                ##Loop to find a valid formula (LEVEL 3)
                # print('Start rejection sampling')
                while True:  #REJECTION SAMPLING 
                    num_rej += 1
                    
                    #If formulas have not been listed, do so
                    if not os.path.exists('grammar'): os.mkdir('grammar')
                    if not os.path.exists(f'grammar/spec_level{spec_level_start}_numbquant{m}'): os.mkdir(f'grammar/spec_level{spec_level_start}_numbquant{m}')
                    if not os.path.exists(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}'): 
                        os.mkdir(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}')
                        listing_nodes = []
                        listing_nodes.append(tree.BinaryTreeNode('start'))
                        
                        listing_nodes =  gramm.build_formulas(l, grammar_structure, spec_level_start, listing_nodes, 0, 'left')
                        gramm.transform_listing_to_nodes(listing_nodes, l, spec_level_start,m)
                        
                        
                        
                    'Create the formula'
                    if len(os.listdir(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}')) != 0:
                        list_formulas = os.listdir(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}')
                        list_formulas.sort()
                        sampled_f = random.sample(list_formulas,1)[0]
                        with open(f'grammar/spec_level{spec_level_start}_numbquant{m}/length_{l}/{sampled_f}', 'rb') as f:
                          nodes = pickle.load(f)
                        
                        nodes = tree.sample_predicate_STL(grammar_predicate, nodes, inputs.bool_different_variables)
                        
                        formula = tree.tree_to_formula(nodes)
                        
                        to_be_present = [ [f'pi[{j}]' for j in range(m) ] , m]
                            
                        if sum([to_be_present[0][i] in formula for i in range(len(to_be_present[0]))])>=to_be_present[1]:
                            
                            #Check if there is a double negation
                            if not tree.check_double_negation(nodes) or num_rej > numb_rej_max :
                                mono = binary_search.check_monotonicity(nodes)
                                if mono != False: # ok monotonic
                                        break
                                
                        #Too many rejected -> sample new length for the formula        
                        if num_rej > numb_rej_max + 1:
                            #Change length
                            l = random.randint(2*math.ceil(m/2)-1,max_length) # sample the length of the formula
                            print(f'Need to change the length. New length = {l}')
                            num_rej = 0
                
                    elif l < 20: 
                        print(f'It was not possible building a formula with length {l}: I will increment length by one. New length = {l+1}\n')
                        l += 1
                    
                    else:
                        print('No formula found. Start with a smaller value of length or change the grammar\n')
                        sys.exit()
    
                # print('End rejection sampling')
    
                index_mono = 0
                for index_epsilon in range(10):
                    if f'epsilon{index_epsilon}-' in formula:
                        if mono[index_mono] == '+' : concrete_value = inputs.par_bounds[index_epsilon][1]#the maximum
                        elif mono[index_mono] == '-': concrete_value = inputs.par_bounds[index_epsilon][0] #the minimum
                        formula = formula.replace(f'epsilon{index_epsilon}-', f'{concrete_value}')
                        index_mono += 1
    
                cost = se.efficient_monitoring(inputs.traces, formula, quantifiers, inputs.bool_temporal_operators, None,  inputs.bool_interface_aware, inputs.input_variables, inputs.output_variables)
                time_end = time.time()
               
                
                if cost != 'runout' and cost < 0 and not( set([print_formula(quantifiers, formula)]).issubset(database_set)): 
                         
                         # Instantiate nodes
                         p_formula = tree.tree_to_formula(nodes)#Replace the formula with concrete values
                         index_mono = 0
                         for index_epsilon in range(10):
                             if f'epsilon{index_epsilon}-' in p_formula:
                                if mono[index_mono] == '+' : concrete_value = inputs.par_bounds[index_epsilon][1]#the maximum
                                elif mono[index_mono] == '-': concrete_value = inputs.par_bounds[index_epsilon][0] #the minimum
                                for node in nodes: node.data = node.data.replace(f'epsilon{index_epsilon}-', f'{concrete_value}')
                                index_mono += 1
                         database_set = database_set.union(set([print_formula(quantifiers, formula)]))
                         formula_to_be_stored = gramm.StoreFormula(nodes, quantifiers, cost, time_end-time_start,l)
                         database.append(formula_to_be_stored)
                         time_start = time.time()
                         print(f'Formula: {print_formula(formula_to_be_stored.quantifiers, tree.NATURAL_tree_to_formula(formula_to_be_stored.nodes))};\ncost = {formula_to_be_stored.cost}, length={formula_to_be_stored.length}, time = {formula_to_be_stored.time}\n\n')
                         
                         break
                     
                        
        return database                 
  
            
               
