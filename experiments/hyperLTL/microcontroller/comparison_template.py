#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../')
import time
import numpy as np

import inputs
import random
import data_collection as dc
import syntax_guided
import syntax_eval as se
import tree

#####################
#First batch of data#
#####################

file_name = 'Microcontroller.text'

with open(f'{file_name}','a') as file:
    file.write('TEMPLATE SEARCH \n\n\nExperiment 1:\n\n')


numb_traces = 50  # number of traces
numb_formulas = 50 # number of formulas to be learned


seed_learning_1 =  442331#random.randint(0, 1000000)
seed_learning_template_1 = 465442#random.randint(0, 1000000)

learning_traces_1 , positive_test_traces_1 , min_length = dc.data_microcontroller(numb_traces, seed_learning_1)

start_time_learning_1 = time.time()


n = len(learning_traces_1[0])  # number of variables

# Learn with template

#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall forall psi0 | phi0
'''
#Definition of the grammar for the formula structure        
grammar_structure = '''
phi0 == phi1 implies phi3
phi1 == always phi2
phi2 == always phi2 | eventually phi2 | phi2 until phi2 | s_next phi2 | s_prev phi2 | phi2 and phi2 | phi2 or phi2 | phi2 implies phi2 | phi2 equiv phi2 | not phi2 | predicate0
phi3 == always phi4
phi4 == always phi4 | eventually phi4 | phi4 until phi4 | s_next phi4 | s_prev phi4 | phi4 and phi4 | phi4 or phi4 | phi4 implies phi4 | phi4 equiv phi4 | not phi4 | predicate0
'''
#Definition of the grammar for  predicates
grammar_predicate = [ ]
predicate0 = [['=='], list(range(0,n)), 'all tuples' ] #[ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])

input_data_template_1 = inputs.InputData(
    seed=seed_learning_template_1,
    traces=learning_traces_1,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_template_1.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_template_1.store_enumeration_formulas = False
input_data_template_1.store_formulas = True
input_data_template_1.n_target = numb_formulas
input_data_template_1.name_set_mined_formulas = 'template_'

output_data_template_1 = syntax_guided.main(input_data_template_1)

nb_false_positives_1 = 0
nb_false_negatives_1 = 0
nb_false_positives_1_runout = 0
nb_false_negatives_1_runout = 0

for index, item in enumerate(output_data_template_1):
       
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
      
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

with open(f'{file_name}','a') as file:
    file.write(f'\n|T_1| = {len(output_data_template_1)}')
  
#####################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nComputational time (seconds):\n')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
   
######################
# SAVE LIST OF SEEDS #  
######################

with open(f'{file_name}','a') as file:
    file.write('\n\n\n\n\n\nList of seeds used:\n')
    file.write(f'\n\nseed_learning_1 = {seed_learning_1}\nseed_learning_template_1 ={seed_learning_template_1}')
    


#####################
#     TIME          #  
#####################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nComputational time (seconds):\n')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
   


    
        
    