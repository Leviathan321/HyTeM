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
import syntax_eval as se
import tree


file_name = 'SerialAdder_cumulative.text'

numb_traces = 50  # number of traces
length_time = 1000  # number of time steps
numb_formulas = 20 # number of formulas to be learned
monitoring = []
fitness = []


seed_learning = 707631#random.randint(0, 1000000)
learning_traces = dc.data_serial_adder(length_time, numb_traces, seed_learning)
n = len(learning_traces[0])  # number of variables

#############################################
# Learn constrained: efficient monitoring   #
#############################################

seed_learning_monitoring = 861319#random.randint(0, 1000000)

with open(f'{file_name}','a') as file:
    file.write(f'seed_learning = {seed_learning}\n\n')
    file.write(f'\nseed_learning_monitoring = {seed_learning_monitoring}')
    file.write('Monitoring \n\n')

start_time_learning = time.time()


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


input_data_monitoring= inputs.InputData(
    seed= seed_learning_monitoring ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_monitoring.store_enumeration_formulas = False
input_data_monitoring.store_formulas = True
input_data_monitoring.n_target = numb_formulas
input_data_monitoring.name_set_mined_formulas = 'monitoring_'
input_data_monitoring.bool_mutation_quantifiers = False

output_data_monitoring = syntax_guided.main(input_data_monitoring)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\ntime_learning (total) = {time_learning} seconds\n\n')

for index, item in enumerate(output_data_monitoring):
    monitoring.append(item.time/60)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')




#################################
# Learn constrained:  fitness   #
#################################

seed_learning_fitness = 123402#random.randint(0, 1000000)

with open(f'{file_name}','a') as file:
    file.write(f'\nseed_learning_fitness = {seed_learning_fitness}')
    file.write('Fitness \n\n')

start_time_learning = time.time()


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


input_data_fitness= inputs.InputData(
    seed= seed_learning_fitness ,
    traces=learning_traces,
    numb_quantifiers=[2, 2],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=7)

# Additional options:
input_data_fitness.store_enumeration_formulas = False
input_data_fitness.store_formulas = True
input_data_fitness.n_target = numb_formulas
input_data_fitness.name_set_mined_formulas = 'fitness_'
input_data_fitness.bool_mutation_quantifiers = False
input_data_fitness.kind_score = 'fitness'

output_data_fitness = syntax_guided.main(input_data_fitness)


end_time_learning = time.time()  
time_learning = end_time_learning - start_time_learning

with open(f'{file_name}','a') as file:
    file.write(f'\ntime_learning (total) = {time_learning} seconds\n\n')

for index, item in enumerate(output_data_fitness):
    fitness.append(item.time/60)
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')




############################
##   PLOT ##################
############################



value_x = np.arange(1,21)



plt.plot(value_x,np.cumsum(monitoring), 'o-', label = 'Efficient Monitoring' )
plt.plot(value_x,np.cumsum(fitness), 'o-', label = 'Fitness function')
plt.xlabel('Number of learned formulas', fontsize = 50)    
plt.ylabel('Time', fontsize = 50) 
plt.xticks(fontsize=50)  
plt.xticks([5, 10, 15, 20])  
plt.yticks(fontsize=50) 
plt.legend(fontsize = 40)
plt.savefig(f"time_cumlative_comparison.pdf", format='pdf', bbox_inches="tight")  



