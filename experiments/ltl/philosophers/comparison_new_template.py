#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
import syntax_guided
import syntax_eval as se
import tree
import data_collection as dc
import log_dining_philosophers


#####################
#First batch of data#
#####################

file_name = 'Philosophers_NEW_template.text'

with open(f'{file_name}','a') as file:
    file.write('TEMPLATE SEARCH \n\n\nExperiment 1:\n\n')

numb_formulas = 50 # number of formulas to be learned

start_time_learning_1 = time.time()

#seed_learning_1 = random.randint(0, 1000000)

seed_learning_template_1 = 365756#random.randint(0, 1000000)

learning_traces = log_dining_philosophers.dining_philosophers()


n = len(learning_traces[0])  # number of variables

# Learn with template

# Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''
grammar_structure = '''
phi0 == always phi1
phi1 == phi2 implies phi3
phi2 ==  predicate0
phi3 == not phi4
phi4 == s_next phi2
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)),  'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])

input_data_template_1 = inputs.InputData(
    seed=seed_learning_template_1,
    traces=learning_traces,
    numb_quantifiers=[1, 1],  # Exactly one quantifier
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_template_1.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_template_1.second_variable = [0,1,2] # 3 actions
input_data_template_1.to_be_present = [[f'pi[0][{agent}] == {action}' for agent in range(n) for action in input_data_template_1.second_variable ], 2 ]

input_data_template_1.l = 6
input_data_template_1.store_enumeration_formulas = False
input_data_template_1.store_formulas = True
input_data_template_1.n_target = numb_formulas
input_data_template_1.name_set_mined_formulas = 'NEW_template_'

output_data_template_1 = syntax_guided.main(input_data_template_1)


for index, item in enumerate(output_data_template_1):
        
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
   
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

  
 
  
##################################################
#     VARIABILITY SAME TRAINING SAMPLES          #
##################################################

with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nVariability (SAME TRAINING SAMPLES):\n')
    file.write(f'\n|T_1|      = {len(output_data_template_1)}')
 
  
#####################
#     TIME          #  
#####################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nComputational time (seconds):\n')
    file.write(f'\ntime_learning_1 = {time_learning_1}')



