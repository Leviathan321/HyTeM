#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, '../../../')
import time
import numpy as np
import matplotlib.pyplot as plt
import inputs
import random
import data_collection as dc
import syntax_guided
import tree

##SETUP 

file_name = 'SerialAdder_bars.text'

numb_traces = 50  # number of traces
length_time = 1000  # number of time steps
numb_formulas = 51 # number of formulas to be learned
free = []
constrained = []
template= []


seed_learning = 420810 #random.randint(0, 1000000)
learning_traces = dc.data_serial_adder(length_time, numb_traces, seed_learning)
n = len(learning_traces[0])  # number of variables

########################
# Learn free: length 3 #
########################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning ={seed_learning}\n\n')
    file.write('FREE \n\n\nLength 3:\n\n')

start_time_learning = time.time()
seed_learning_free_3 = 21716 #random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)),
              'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_free_3 = inputs.InputData(
    seed= seed_learning_free_3 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_free_3.l = 3
input_data_free_3.store_enumeration_formulas = False
input_data_free_3.store_formulas = True
input_data_free_3.n_target = numb_formulas
input_data_free_3.name_set_mined_formulas = 'free3_'
input_data_free_3.bool_mutation_quantifiers = True

output_data_free_3 = syntax_guided.main(input_data_free_3)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_free_3 = {seed_learning_free_3}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_free_3):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
free.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')



########################
# Learn free: length 4 #
########################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning ={seed_learning}\n\n')
    file.write('FREE \n\n\nLength 4:\n\n')

start_time_learning = time.time()
seed_learning_free_4 = 815027#random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)),
              'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_free_4 = inputs.InputData(
    seed= seed_learning_free_4 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_free_4.l = 4
input_data_free_4.store_enumeration_formulas = False
input_data_free_4.store_formulas = True
input_data_free_4.n_target = numb_formulas
input_data_free_4.name_set_mined_formulas = 'free4_'
input_data_free_4.bool_mutation_quantifiers = True

output_data_free_4 = syntax_guided.main(input_data_free_4)

end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_free_4 = {seed_learning_free_4}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_free_4):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
free.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')
    



########################
# Learn free: length 5 #
########################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning ={seed_learning}\n\n')
    file.write('FREE \n\n\nLength 5:\n\n')

start_time_learning = time.time()
seed_learning_free_5 = 518802 #random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)),
              'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_free_5 = inputs.InputData(
    seed= seed_learning_free_5 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_free_5.l = 5
input_data_free_5.store_enumeration_formulas = False
input_data_free_5.store_formulas = True
input_data_free_5.n_target = numb_formulas
input_data_free_5.name_set_mined_formulas = 'free5_'
input_data_free_5.bool_mutation_quantifiers = True

output_data_free_5 = syntax_guided.main(input_data_free_5)

end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_free_5 = {seed_learning_free_5}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_free_5):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
free.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')



########################
# Learn free: length 6 #
########################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning ={seed_learning}\n\n')
    file.write('FREE \n\n\nLength 6:\n\n')

start_time_learning = time.time()
seed_learning_free_6 = 952825#random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)),
              'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_free_6 = inputs.InputData(
    seed= seed_learning_free_6 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_free_6.l = 6
input_data_free_6.store_enumeration_formulas = False
input_data_free_6.store_formulas = True
input_data_free_6.n_target = numb_formulas
input_data_free_6.name_set_mined_formulas = 'free6_'
input_data_free_6.bool_mutation_quantifiers = True

output_data_free_6 = syntax_guided.main(input_data_free_6)

end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_free_6 = {seed_learning_free_6}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_free_6):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
free.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')




########################
# Learn free: length 7 #
########################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning ={seed_learning}\n\n')
    file.write('FREE \n\n\nLength 7:\n\n')

start_time_learning = time.time()
seed_learning_free_7 =  161208 #random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)),
              'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_free_7 = inputs.InputData(
    seed= seed_learning_free_7 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_free_7.l = 7
input_data_free_7.store_enumeration_formulas = False
input_data_free_7.store_formulas = True
input_data_free_7.n_target = numb_formulas
input_data_free_7.name_set_mined_formulas = 'free7_'
input_data_free_7.bool_mutation_quantifiers = True

output_data_free_7 = syntax_guided.main(input_data_free_7)

end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_free_7 = {seed_learning_free_7}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_free_7):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
free.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')

  

###############################
# Learn constrained: length 3 #
###############################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning = {seed_learning}\n\n')
    file.write('CONSTRAINED \n\n\nLength 3:\n\n')

start_time_learning = time.time()
seed_learning_constrained_3 = 114655 #random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_constrained_3 = inputs.InputData(
    seed= seed_learning_constrained_3 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_constrained_3.l = 3
input_data_constrained_3.store_enumeration_formulas = False
input_data_constrained_3.store_formulas = True
input_data_constrained_3.n_target = numb_formulas
input_data_constrained_3.name_set_mined_formulas = 'constrained3_'
input_data_constrained_3.bool_mutation_quantifiers = False

output_data_constrained_3 = syntax_guided.main(input_data_constrained_3)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_constrained_3 = {seed_learning_constrained_3}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_constrained_3):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
constrained.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')


###############################
# Learn constrained: length 4 #
###############################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning = {seed_learning}\n\n')
    file.write('CONSTRAINED \n\n\nLength 4:\n\n')

start_time_learning = time.time()
seed_learning_constrained_4 = 323956 #random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_constrained_4= inputs.InputData(
    seed= seed_learning_constrained_4 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_constrained_4.l = 4
input_data_constrained_4.store_enumeration_formulas = False
input_data_constrained_4.store_formulas = True
input_data_constrained_4.n_target = numb_formulas
input_data_constrained_4.name_set_mined_formulas = 'constrained4_'
input_data_constrained_4.bool_mutation_quantifiers = False

output_data_constrained_4 = syntax_guided.main(input_data_constrained_4)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_constrained_4 = {seed_learning_constrained_4}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_constrained_4):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
constrained.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')



###############################
# Learn constrained: length 5 #
###############################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning = {seed_learning}\n\n')
    file.write('CONSTRAINED \n\n\nLength 5:\n\n')

start_time_learning = time.time()
seed_learning_constrained_5 = 864447#random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_constrained_5 = inputs.InputData(
    seed= seed_learning_constrained_5 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_constrained_5.l = 5
input_data_constrained_5.store_enumeration_formulas = False
input_data_constrained_5.store_formulas = True
input_data_constrained_5.n_target = numb_formulas
input_data_constrained_5.name_set_mined_formulas = 'constrained5_'
input_data_constrained_5.bool_mutation_quantifiers = False

output_data_constrained_5= syntax_guided.main(input_data_constrained_5)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_constrained_5 = {seed_learning_constrained_5}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_constrained_5):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
constrained.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')




###############################
# Learn constrained: length 6 #
###############################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning = {seed_learning}\n\n')
    file.write('CONSTRAINED \n\n\nLength 6:\n\n')

start_time_learning = time.time()
seed_learning_constrained_6 = 721436#random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_constrained_6 = inputs.InputData(
    seed= seed_learning_constrained_6 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_constrained_6.l = 6
input_data_constrained_6.store_enumeration_formulas = False
input_data_constrained_6.store_formulas = True
input_data_constrained_6.n_target = numb_formulas
input_data_constrained_6.name_set_mined_formulas = 'constrained6_'
input_data_constrained_6.bool_mutation_quantifiers = False

output_data_constrained_6= syntax_guided.main(input_data_constrained_6)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_constrained_6 = {seed_learning_constrained_6}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_constrained_6):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
constrained.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')



###############################
# Learn constrained: length 7 #
###############################
with open(f'{file_name}','a') as file:
    file.write(f'seed_learning = {seed_learning}\n\n')
    file.write('CONSTRAINED \n\n\nLength 7:\n\n')

start_time_learning = time.time()
seed_learning_constrained_7= 553214#random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | s_next phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | phi0 equiv phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_constrained_7 = inputs.InputData(
    seed= seed_learning_constrained_7 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_constrained_7.l = 7
input_data_constrained_7.store_enumeration_formulas = False
input_data_constrained_7.store_formulas = True
input_data_constrained_7.n_target = numb_formulas
input_data_constrained_7.name_set_mined_formulas = 'constrained7_'
input_data_constrained_7.bool_mutation_quantifiers = False

output_data_constrained_7 = syntax_guided.main(input_data_constrained_7)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_constrained_7 = {seed_learning_constrained_7}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_constrained_7):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
constrained.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')




###############################
# Learn template: length 5 #
###############################
with open(f'{file_name}','a') as file:
    file.write('TEMPLATE \n\n\nLength 5:\n\n')

start_time_learning = time.time()
seed_learning_template_5= random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''

# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == phi1 implies phi3
phi1 == always phi2
phi2 == phi2 and phi2 | phi2 or phi2 | phi2 implies phi2 | phi2 equiv phi2 | not phi2 | predicate0
phi3 == always phi4
phi4 == phi4 and phi4 | phi4 or phi4 | phi4 implies phi4 | phi4 equiv phi4 | not phi4 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_template_5 = inputs.InputData(
    seed= seed_learning_template_5 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_template_5.l = 5
input_data_template_5.store_enumeration_formulas = False
input_data_template_5.store_formulas = True
input_data_template_5.n_target = numb_formulas
input_data_template_5.name_set_mined_formulas = 'template5_'
input_data_template_5.bool_mutation_quantifiers = False

output_data_template_5 = syntax_guided.main(input_data_template_5)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_template_5 = {seed_learning_template_5}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_template_5):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
template.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')
    
    
    
    
###############################
# Learn template: length 6 #
###############################
with open(f'{file_name}','a') as file:
    file.write('TEMPLATE \n\n\nLength 6:\n\n')

start_time_learning = time.time()
seed_learning_template_6= random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''

# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == phi1 implies phi3
phi1 == always phi2
phi2 == phi2 and phi2 | phi2 or phi2 | phi2 implies phi2 | phi2 equiv phi2 | not phi2 | predicate0
phi3 == always phi4
phi4 == phi4 and phi4 | phi4 or phi4 | phi4 implies phi4 | phi4 equiv phi4 | not phi4 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_template_6 = inputs.InputData(
    seed= seed_learning_template_6 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_template_6.l = 6
input_data_template_6.store_enumeration_formulas = False
input_data_template_6.store_formulas = True
input_data_template_6.n_target = numb_formulas
input_data_template_6.name_set_mined_formulas = 'template6_'
input_data_template_6.bool_mutation_quantifiers = False

output_data_template_6 = syntax_guided.main(input_data_template_6)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_template_6 = {seed_learning_template_6}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_template_6):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
template.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')
    
    


###############################
# Learn template: length 7 #
###############################
with open(f'{file_name}','a') as file:
    file.write('TEMPLATE \n\n\nLength 7:\n\n')

start_time_learning = time.time()
seed_learning_template_7= 580289 #random.randint(0, 1000000)

grammar_quantifiers = '''
psi0 == forall psi0 | phi0
'''

# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == phi1 implies phi3
phi1 == always phi2
phi2 == phi2 and phi2 | phi2 or phi2 | phi2 implies phi2 | phi2 equiv phi2 | not phi2 | predicate0
phi3 == always phi4
phi4 == phi4 and phi4 | phi4 or phi4 | phi4 implies phi4 | phi4 equiv phi4 | not phi4 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = []
predicate0 = [['=='], list(range(0, n)), 'all tuples']  # [ item for item in itertools.permutations(range(number_quantifiers), min(number_quantifiers,2))]
grammar_predicate.append([predicate0])


input_data_template_7 = inputs.InputData(
    seed= seed_learning_template_7 ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_template_7.l = 7
input_data_template_7.store_enumeration_formulas = False
input_data_template_7.store_formulas = True
input_data_template_7.n_target = numb_formulas
input_data_template_7.name_set_mined_formulas = 'template7_'
input_data_template_7.bool_mutation_quantifiers = False

output_data_template_7 = syntax_guided.main(input_data_template_7)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_template_7 = {seed_learning_template_7}')
    file.write(f'\ntime_learning (total) = {time_learning} \n\n')
aux = []
for index, item in enumerate(output_data_template_7):
    if index !=0: 
        aux.append(item.time)
        with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
template.append(np.mean(aux))
with open(f'{file_name}','a') as file:
    file.write(f'mean time = {np.mean(aux)}')
    
    
    
    
###########################
######### PLOTS ############
###########################
lengths = [3, 4, 5 , 6, 7]

X_axis = np.arange(len(lengths))  


plt.bar(X_axis - 0.5, free, width = 0.25, label = 'Free quant.')
plt.bar(X_axis - 0.25, constrained , width = 0.25, label = 'Constrained quant.')

plt.bar( X_axis[2:] , template , width = 0.25, label = 'Template')



plt.xticks([-0.2, 0.8, 1.8, 2.8, 4], lengths, fontsize=40) 
plt.xlabel("Formula length", fontsize = 40)
plt.ylabel("Time", fontsize = 40)   

plt.yticks(fontsize=40) 
plt.legend(fontsize = 30)
plt.savefig("plot_bars.pdf", format='pdf', bbox_inches="tight")



