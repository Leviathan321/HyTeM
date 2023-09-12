#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
sys.path.insert(0, '../../../')
import time
import numpy as np
import matplotlib.pyplot as plt
import random
import tree
import math

import inputs
import data_collection as dc
import syntax_guided
import syntax_eval as se
import pickle
'''
IMPORTANT: 
    in syntax_eval.py, inside the function monitor_hyperSTL_correctness
set the boolean variable bool_return_counter to True.

SET IT BACK TO FALSE WHEN FINISHING USING THIS SCRIPT!
'''

    
file_name = 'Parking_monitor.text'


numb_traces = 185  # number of traces
numb_formulas = 10 # number of formulas to be learned
fraction = 1 # threshold x binary search unknown region
              
start_time_learning_1 = time.time()

seed_learning_1 =   502259#random.randint(0, 1000000)

data_path = '../../../data/DataDemo/'

learning_traces , _, min_length = dc.collect_parking_data(numb_traces, seed_learning_1, path =  data_path)

with open(f'{file_name}','a') as file:
    file.write(f'seed_learning_1 = {seed_learning_1}')

  
###########################    
# FORALL FORALL ###########    
###########################      

with open(f'{file_name}','a') as file:
    file.write('\n\n\n FORALL FORALL\n\n')
seed_learning_FF = 136918#random.randint(0, 1000000)


#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall forall psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']

input_data_FF = inputs.InputData(
    seed=seed_learning_FF,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)
 

# Additional options:
input_data_FF.learning_stl = True
input_data_FF.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_FF.fraction = fraction

input_data_FF.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_FF.store_enumeration_formulas = False
input_data_FF.store_formulas = True
input_data_FF.n_target = numb_formulas
input_data_FF.name_set_mined_formulas = 'FF_'
input_data_FF.return_violated_formulas = True

output_data_FF = syntax_guided.main(input_data_FF)
violated_formulas_FF = syntax_guided.sample_violated_formulas_stl(input_data_FF)
        
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

counter_monitor_FF = []
counter_monitorSTL_FF = []
with open(f'{file_name}','a') as file:
    file.write(f'\nSatisfied formulas:')
    
for index, item in enumerate(output_data_FF):
    outcome, counter = se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers, True, False,  False, [], [])
    counter_monitor_FF.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers)
    counter_monitorSTL_FF.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')

with open(f'{file_name}','a') as file:
    file.write(f'\nViolated formulas:')
    
for index, item in enumerate(violated_formulas_FF):
    outcome, counter = se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes) , item.quantifiers , True, False,  False, [], [])
    counter_monitor_FF.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers )
    counter_monitorSTL_FF.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula( item.quantifiers,tree.NATURAL_tree_to_formula(item.nodes))};\n\n')



with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_FF = {seed_learning_FF}')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
    file.write(f'\nmean(counter_monitor_FF) = {np.mean(counter_monitor_FF)}')
    file.write(f'\nmean(counter_monitorSTL_FF) = {np.mean(counter_monitorSTL_FF)}')



###########################    
# FORALL EXISTS ###########    
###########################      

with open(f'{file_name}','a') as file:
    file.write('\n\n\n FORALL EXISTS\n\n')
seed_learning_FE = 729717#random.randint(0, 1000000)


#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']

input_data_FE = inputs.InputData(
    seed=seed_learning_FE,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)
 

# Additional options:
input_data_FE.learning_stl = True
input_data_FE.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_FE.fraction = fraction

input_data_FE.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_FE.store_enumeration_formulas = False
input_data_FE.store_formulas = True
input_data_FE.n_target = numb_formulas
input_data_FE.name_set_mined_formulas = 'FE_'
input_data_FE.return_violated_formulas = True

output_data_FE = syntax_guided.main(input_data_FE)
violated_formulas_FE = syntax_guided.sample_violated_formulas_stl(input_data_FE)
        
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

counter_monitor_FE = []
counter_monitorSTL_FE = []

with open(f'{file_name}','a') as file:
    file.write(f'\nSatisfied formulas:')

for index, item in enumerate(output_data_FE):
    outcome, counter =se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers, True, False,  False, [], [])
    counter_monitor_FE.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers)
    counter_monitorSTL_FE.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')


with open(f'{file_name}','a') as file:
    file.write(f'\nViolated formulas:')
    
for index, item in enumerate(violated_formulas_FE):
    outcome, counter = se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes) , item.quantifiers , True, False,  False, [], [])
    counter_monitor_FE.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers )
    counter_monitorSTL_FE.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula( item.quantifiers,tree.NATURAL_tree_to_formula(item.nodes))};\n\n')


with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_FE = {seed_learning_FE}')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
    file.write(f'\nmean(counter_monitor_FE) = {np.mean(counter_monitor_FE)}')
    file.write(f'\nmean(counter_monitorSTL_FE) = {np.mean(counter_monitorSTL_FE)}')


###########################    
# EXISTS FORALL ###########    
###########################      

with open(f'{file_name}','a') as file:
    file.write('\n\n\n EXISTS FORALL\n\n')
    
seed_learning_EF = 739984 #random.randint(0, 1000000)


#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == exists forall psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']

input_data_EF = inputs.InputData(
    seed=seed_learning_EF,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)
 

# Additional options:
input_data_EF.learning_stl = True
input_data_EF.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_EF.fraction = fraction

input_data_EF.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_EF.store_enumeration_formulas = False
input_data_EF.store_formulas = True
input_data_EF.n_target = numb_formulas
input_data_EF.name_set_mined_formulas = 'EF_'
input_data_EF.return_violated_formulas = True

output_data_EF = syntax_guided.main(input_data_EF)
violated_formulas_EF = syntax_guided.sample_violated_formulas_stl(input_data_EF)

        
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

counter_monitor_EF = []
counter_monitorSTL_EF = []


with open(f'{file_name}','a') as file:
    file.write(f'\nSatisfied formulas:')
    
for index, item in enumerate(output_data_EF):
    outcome, counter =se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers, True, False,  False, [], [])
    counter_monitor_EF.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers)
    counter_monitorSTL_EF.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')

with open(f'{file_name}','a') as file:
    file.write(f'\nViolated formulas:')
    
for index, item in enumerate(violated_formulas_EF):
    outcome, counter = se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes) , item.quantifiers , True, False,  False, [], [])
    counter_monitor_EF.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers )
    counter_monitorSTL_EF.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula( item.quantifiers,tree.NATURAL_tree_to_formula(item.nodes))};\n\n')



with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_EF = {seed_learning_EF}')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
    file.write(f'\nmean(counter_monitor_EF) = {np.mean(counter_monitor_EF)}')
    file.write(f'\nmean(counter_monitorSTL_EF) = {np.mean(counter_monitorSTL_EF)}')


###########################    
# EXISTS EXISTS ###########    
###########################      

with open(f'{file_name}','a') as file:
    file.write('\n\n\n EXISTS EXISTS\n\n')
seed_learning_EE = 364481#random.randint(0, 1000000)


#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == exists exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']

input_data_EE = inputs.InputData(
    seed=seed_learning_EE,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)
 

# Additional options:
input_data_EE.learning_stl = True
input_data_EE.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_EE.fraction = fraction

input_data_EE.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_EE.store_enumeration_formulas = False
input_data_EE.store_formulas = True
input_data_EE.n_target = numb_formulas
input_data_EE.name_set_mined_formulas = 'EE_'
    
input_data_EE.return_violated_formulas = True

output_data_EE  = syntax_guided.main(input_data_EE)
violated_formulas_EE = syntax_guided.sample_violated_formulas_stl(input_data_EE)

        
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

counter_monitor_EE = []
counter_monitorSTL_EE = []

with open(f'{file_name}','a') as file:
    file.write(f'\nSatisfied formulas:')
    
for index, item in enumerate(output_data_EE):
    outcome, counter =se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers, True, False,  False, [], [])
    counter_monitor_EE.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers)
    counter_monitorSTL_EE.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')

with open(f'{file_name}','a') as file:
    file.write(f'\nViolated formulas:')
    
for index, item in enumerate(violated_formulas_EE):
    outcome, counter = se.efficient_monitoring(learning_traces, tree.tree_to_formula(item.nodes) , item.quantifiers , True, False,  False, [], [])
    counter_monitor_EE.append(counter)
    outcome_stl, counter_stl, skip_stl = se.monitor_hyperSTL_correctness(learning_traces, tree.tree_to_formula(item.nodes), item.quantifiers )
    counter_monitorSTL_EE.append(counter_stl)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula( item.quantifiers,tree.NATURAL_tree_to_formula(item.nodes))};\n\n')



with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_EE = {seed_learning_EE}')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
    file.write(f'\nmean(counter_monitor_EE) = {np.mean(counter_monitor_EE)}')
    file.write(f'\nmean(counter_monitorSTL_EE) = {np.mean(counter_monitorSTL_EE)}')




#####################
##  PLOT ############
#####################



lengths = [r'$ \forall \forall $ ', r'$ \forall \exists $', r'$ \exists \forall $' , r'$ \exists \exists $']
  
X_axis = np.arange(len(lengths))  

fitness = [math.log10(len(learning_traces)*(len(learning_traces)-1))]*4
eff_mon = [math.log10(np.mean(counter_monitor_FF)), math.log10(np.mean(counter_monitor_FE)) , math.log10(np.mean(counter_monitor_EF)), math.log10(np.mean(counter_monitor_EE))]
monitor_stl =[math.log10(np.mean(counter_monitorSTL_FF)), math.log10(np.mean(counter_monitorSTL_FE)), math.log10(np.mean(counter_monitorSTL_EF)), math.log10(np.mean(counter_monitorSTL_EE))]

plt.bar(X_axis - 0.4,  fitness , width = 0.2, label = 'Fitness', color = 'grey')
plt.bar(X_axis - 0.2, eff_mon , width = 0.2, label = 'Effective  monitoring', color = 'darkmagenta')
plt.bar(X_axis - 0,  monitor_stl, width = 0.2, label = 'Effective  monitoring + correctness', color = 'sandybrown')

plt.xticks([-0.2, 0.8, 1.8, 2.8], lengths, fontsize=50) 
   
plt.yticks(fontsize=50) 
plt.legend(fontsize = 30)#,  bbox_to_anchor=(0.3, 0.4, 0.5, 0.5))
plt.savefig("plot_exp_counter.pdf", format='pdf', bbox_inches="tight")  
plt.close()