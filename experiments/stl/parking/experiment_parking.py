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

file_name = 'Parking.text'
union_set_different_samples = set()
union_set_same_samples = set()
set_T1 = set()
set_T2 = set()
set_T3 = set()
set_T1_run2 = set()
set_T1_run3 = set()

with open(f'{file_name}','a') as file:
    file.write('TEMPLATE SEARCH \n\n\nExperiment 1:\n\n')


numb_traces = 50  # number of traces
numb_formulas = 50 # number of formulas to be learned
fraction = 10 # threshold x binary search unknown region

start_time_learning_1 = time.time()
# random.seed(1)
seed_learning_1 = 991781#random.randint(0, 1000000)
seed_negative_test_1 = 961339#random.randint(0, 1000000)
seed_learning_template_1 =342134 #random.randint(0, 1000000)

learning_traces_1, positive_test_traces_1, min_length = dc.collect_parking_data(numb_traces, seed_learning_1)
negative_test_traces_1 = dc.data_random_binary_real(min_length, numb_traces, seed_negative_test_1)


n = len(learning_traces_1[0])  # number of variables

# Learn with template

#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall forall psi0 | phi0
'''
#Definition of the grammar for the formula structure        
grammar_structure = '''
phi0 == always phi1
phi1 == phi2 implies phi2
phi2 == eventually phi2 | always phi2 | phi2 until phi2 | phi2 implies phi2 | phi2 and phi2 | phi2 or phi2 | not phi2 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']


input_data_template_1 = inputs.InputData(
    seed=seed_learning_template_1,
    traces=learning_traces_1,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)
 

# Additional options:
input_data_template_1.learning_stl = True
input_data_template_1.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_template_1.fraction = fraction
input_data_template_1.kind_score = 'efficient_monitoring_stl'

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
    positive_result = se.efficient_monitoring(positive_test_traces_1, tree.tree_to_formula(item.nodes), item.quantifiers,
                            True, None, False, [], [])
    negative_result = se.efficient_monitoring(negative_test_traces_1, tree.tree_to_formula(item.nodes),
                                              item.quantifiers,
                                              True, None, False, [], [])
    
    if positive_result == 'runout' : nb_false_positives_1_runout += 1
    elif positive_result == -1: nb_false_negatives_1 = nb_false_negatives_1 + 1
    
    if negative_result == 'runout': nb_false_negatives_1_runout += 1
    elif negative_result == 1: nb_false_positives_1 = nb_false_positives_1 + 1
        
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    
    union_set_different_samples = union_set_different_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    union_set_same_samples = union_set_same_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_T1 = set_T1.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    
ratio_false_negatives_1 = nb_false_negatives_1/(len(output_data_template_1))
ratio_false_positives_1 = nb_false_positives_1/(len(output_data_template_1))
        
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

  

########################      
# Second batch of data #
########################

with open(f'{file_name}','a') as file: file.write('\n\nExperiment 2:\n\n')


start_time_learning_2 = time.time()

seed_learning_2 = 361441 #random.randint(0, 1000000)
seed_negative_test_2 = 817553#random.randint(0, 1000000)
seed_learning_template_2 = 832195#random.randint(0, 1000000)



learning_traces_2, positive_test_traces_2, min_length = dc.collect_parking_data(numb_traces, seed_learning_2)
negative_test_traces_2 = dc.data_random_binary_real(min_length, numb_traces, seed_negative_test_2)

n = len(learning_traces_2[0])  # number of variables



# Learn with template
# Definition of the grammar for the quantifiers
#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall forall psi0 | phi0
'''
#Definition of the grammar for the formula structure        
grammar_structure = '''
phi0 == always phi1
phi1 == phi2 implies phi2
phi2 == eventually phi2 | always phi2 | phi2 until phi2 | phi2 implies phi2 | phi2 and phi2 | phi2 or phi2 | not phi2 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']


input_data_template_2 = inputs.InputData(
    seed=seed_learning_template_2,
    traces=learning_traces_2,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_template_2.learning_stl = True
input_data_template_2.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_template_2.fraction = fraction
input_data_template_2.kind_score = 'efficient_monitoring_stl'

input_data_template_2.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_template_2.store_enumeration_formulas = False
input_data_template_2.store_formulas = True
input_data_template_2.n_target = numb_formulas
input_data_template_2.name_set_mined_formulas = 'template_'



output_data_template_2 = syntax_guided.main(input_data_template_2)

nb_false_positives_2 = 0
nb_false_negatives_2 = 0
nb_false_positives_2_runout = 0
nb_false_negatives_2_runout = 0
for index, item in enumerate(output_data_template_2):
    positive_result = se.efficient_monitoring(positive_test_traces_2, tree.tree_to_formula(item.nodes), item.quantifiers,
                            True, None, False, [], [])
    negative_result = se.efficient_monitoring(negative_test_traces_2, tree.tree_to_formula(item.nodes),
                                              item.quantifiers,
                                              True, None, False, [], [])
    
    if positive_result == 'runout': nb_false_negatives_2_runout +=1
    elif positive_result == -1:  nb_false_negatives_2 = nb_false_negatives_2 + 1
    
    if negative_result == 'runout': nb_false_positives_2_runout +=1 
    elif negative_result == 1: nb_false_positives_2 = nb_false_positives_2 + 1

    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    union_set_different_samples = union_set_different_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_T2 = set_T2.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    
    
ratio_false_negatives_2 = nb_false_negatives_2/(len(output_data_template_2))
ratio_false_positives_2 = nb_false_positives_2/(len(output_data_template_2))
end_time_learning_2 = time.time()  
time_learning_2 = end_time_learning_2 - start_time_learning_2    
   
 
#####################   
#Third batch of data#
#####################

with open(f'{file_name}','a') as file: file.write('\n\nExperiment 3:\n\n')

start_time_learning_3 = time.time()
seed_learning_3 = 744858#random.randint(0, 1000000)
seed_negative_test_3 = 513738#random.randint(0, 1000000)
seed_learning_template_3 = 882559#random.randint(0, 1000000)


learning_traces_3, positive_test_traces_3, min_length = dc.collect_parking_data(numb_traces, seed_learning_3)
negative_test_traces_3 = dc.data_random_binary_real(min_length, numb_traces, seed_negative_test_3)


n = len(learning_traces_3[0])  # number of variables

# Learn with template

#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall forall psi0 | phi0
'''
#Definition of the grammar for the formula structure        
grammar_structure = '''
phi0 == always phi1
phi1 == phi2 implies phi2
phi2 == eventually phi2 | always phi2 | phi2 until phi2 | phi2 implies phi2 | phi2 and phi2 | phi2 or phi2 | not phi2 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']


input_data_template_3 = inputs.InputData(
    seed=seed_learning_template_3,
    traces=learning_traces_3,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_template_3.learning_stl = True
input_data_template_3.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_template_3.fraction = fraction
input_data_template_3.kind_score = 'efficient_monitoring_stl'

input_data_template_3.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_template_3.store_enumeration_formulas = False
input_data_template_3.store_formulas = True
input_data_template_3.n_target = numb_formulas
input_data_template_3.name_set_mined_formulas = 'template_'


output_data_template_3 = syntax_guided.main(input_data_template_3)

nb_false_positives_3 = 0
nb_false_negatives_3 = 0
nb_false_positives_3_runout = 0
nb_false_negatives_3_runout = 0

for index, item in enumerate(output_data_template_3):
    positive_result = se.efficient_monitoring(positive_test_traces_3, tree.tree_to_formula(item.nodes), item.quantifiers,
                            True, None, False, [], [])
    negative_result = se.efficient_monitoring(negative_test_traces_3, tree.tree_to_formula(item.nodes),
                                              item.quantifiers,
                                              True, None, False, [], [])
    
    if positive_result == 'runout': nb_false_negatives_3_runout += 1
    elif positive_result == -1:  nb_false_negatives_3 = nb_false_negatives_3 + 1
    
    if negative_result == 'runout': nb_false_positives_3_runout += 1
    elif negative_result == 1: nb_false_positives_3 = nb_false_positives_3 + 1
    
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    
    union_set_different_samples = union_set_different_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_T3 = set_T3.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    
ratio_false_negatives_3 = nb_false_negatives_3/(len(output_data_template_3))
ratio_false_positives_3 = nb_false_positives_3/(len(output_data_template_3))

end_time_learning_3 = time.time()  
time_learning_3 = end_time_learning_3 - start_time_learning_3   


    
##################################################
#     VARIABILITY DIFFERENT TRAINING SAMPLES     #
##################################################


variability_different_samples = (len(union_set_different_samples.difference(set_T1)) + len(union_set_different_samples.difference(set_T2)) + len(union_set_different_samples.difference(set_T3)) )/( 6 * np.mean([len(set_T1), len(set_T2), len(set_T3)]))
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nVariability (DIFFERENT TRAINING SAMPLES):\n')
    file.write(f'\n|T_1| = {len(output_data_template_1)}')
    file.write(f'\n|T_2| = {len(output_data_template_2)}')
    file.write(f'\n|T_3| = {len(output_data_template_3)}')
    file.write(f'\n|T|   = {len(union_set_different_samples)}')
    file.write(f'\nVariability   = {variability_different_samples}')



################################
#   RATIO OF FALSE POSITIVES   #       
################################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nRatio of false positives:\n')
    file.write(f'\nratio_false_positives_1 = {ratio_false_positives_1}')
    if nb_false_positives_1_runout != 0 : file.write(f'runout_false_positives_1 = {nb_false_positives_1_runout}')
    
    file.write(f'\nratio_false_positives_2 = {ratio_false_positives_2}')
    if nb_false_positives_2_runout != 0 : file.write(f'runout_false_positives_2 = {nb_false_positives_2_runout}')
    
    file.write(f'\nratio_false_positives_3 = {ratio_false_positives_3}')
    if nb_false_positives_3_runout != 0 : file.write(f'runout_false_positives_3 = {nb_false_positives_3_runout}')
    
    file.write(f'\n\nAverage = {np.mean([ratio_false_positives_1, ratio_false_positives_2, ratio_false_positives_3])} +- {np.std([ratio_false_positives_1, ratio_false_positives_2, ratio_false_positives_3])}')
  
    
################################
#   RATIO OF FALSE NEGATIVES   #       
################################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nRatio of false negatives:\n')
    file.write(f'\nratio_false_negatives_1 = {ratio_false_negatives_1}')
    if nb_false_negatives_1_runout != 0 : file.write(f'runout_false_negatives_1 = {nb_false_negatives_1_runout}')
    
    file.write(f'\nratio_false_negatives_2 = {ratio_false_negatives_2}')
    if nb_false_negatives_2_runout != 0 : file.write(f'runout_false_negatives_2 = {nb_false_negatives_2_runout}')
    
    file.write(f'\nratio_false_negatives_3 = {ratio_false_negatives_3}')
    if nb_false_negatives_3_runout != 0 : file.write(f'runout_false_negatives_3 = {nb_false_negatives_3_runout}')
    
    file.write(f'\n\nAverage = {np.mean([ratio_false_negatives_1, ratio_false_negatives_2, ratio_false_negatives_3])} +- {np.std([ratio_false_negatives_1, ratio_false_negatives_2, ratio_false_negatives_3])}')
  

    
#########################################    
# Second run with the same batch of data#    
#########################################

seed_learning_template_1_run2 = 646289 #random.randint(0, 1000000)

with open(f'{file_name}','a') as file: file.write('\n\n\nSecond run with the same traces as Experiment 1:\n\n')


start_time_learning_1_run2 = time.time() 

n = len(learning_traces_1[0])  # number of variables

# Learn with template
#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall forall psi0 | phi0
'''
#Definition of the grammar for the formula structure        
grammar_structure = '''
phi0 == always phi1
phi1 == phi2 implies phi2
phi2 == eventually phi2 | always phi2 | phi2 until phi2 | phi2 implies phi2 | phi2 and phi2 | phi2 or phi2 | not phi2 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']


input_data_template_1_run2 = inputs.InputData(
    seed=seed_learning_template_1_run2,
    traces=learning_traces_1,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_template_1_run2.learning_stl = True
input_data_template_1_run2.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_template_1_run2.fraction = fraction
input_data_template_1_run2.kind_score = 'efficient_monitoring_stl'

input_data_template_1_run2.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_template_1_run2.store_enumeration_formulas = False
input_data_template_1_run2.store_formulas = True
input_data_template_1_run2.n_target = numb_formulas
input_data_template_1_run2.name_set_mined_formulas = 'template_'


output_data_template_1_run2 = syntax_guided.main(input_data_template_1_run2)

for index, item in enumerate(output_data_template_1_run2):
   
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    union_set_same_samples = union_set_same_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_T1_run2 = set_T1_run2.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
      
end_time_learning_1_run2 = time.time()  
time_learning_1_run2 = end_time_learning_1_run2 - start_time_learning_1_run2    
    
 
#########################################    
# Third run with the same batch of data#    
#########################################

seed_learning_template_1_run3 = 71180 #random.randint(0, 1000000)

with open(f'{file_name}','a') as file: file.write('\n\n\nThird run with the same traces as Experiment 1:\n\n')


start_time_learning_1_run3 = time.time() 

n = len(learning_traces_1[0])  # number of variables

# Learn with template
#Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall forall psi0 | phi0
'''
#Definition of the grammar for the formula structure        
grammar_structure = '''
phi0 == always phi1
phi1 == phi2 implies phi2
phi2 == eventually phi2 | always phi2 | phi2 until phi2 | phi2 implies phi2 | phi2 and phi2 | phi2 or phi2 | not phi2 | predicate0
'''
#Definition of the grammar for predicates
grammar_predicate = ['( abs( pi[0][0] - pi[1][0] ) < epsilon0- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon1- )']

input_data_template_1_run3 = inputs.InputData(
    seed=seed_learning_template_1_run3,
    traces=learning_traces_1,
    numb_quantifiers=[2, 2],  # Exactly 2 quantifiers
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_template_1_run3.learning_stl = True
input_data_template_1_run3.par_bounds = [ [0.1, 50], [0.1, 50] ]
input_data_template_1_run3.fraction = fraction
input_data_template_1_run3.kind_score = 'efficient_monitoring_stl'

input_data_template_1_run3.bool_mutation_quantifiers = False  # Quantifiers not modified during the learning process (since one forall is the only option admissible from grammar_predicate)
input_data_template_1_run3.store_enumeration_formulas = False
input_data_template_1_run3.store_formulas = True
input_data_template_1_run3.n_target = numb_formulas
input_data_template_1_run3.name_set_mined_formulas = 'template_'

output_data_template_1_run3 = syntax_guided.main(input_data_template_1_run3)

for index, item in enumerate(output_data_template_1_run3):
   
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    union_set_same_samples = union_set_same_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_T1_run3 = set_T1_run3.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    
end_time_learning_1_run3 = time.time()  
time_learning_1_run3 = end_time_learning_1_run3 - start_time_learning_1_run3    
 
  
##################################################
#     VARIABILITY SAME TRAINING SAMPLES          #
##################################################
variability_same_samples = (len(union_set_same_samples.difference(set_T1)) + len(union_set_same_samples.difference(set_T1_run2)) + len(union_set_same_samples.difference(set_T1_run3)) )/( 6 * np.mean([len(set_T1), len(set_T1_run2), len(set_T1_run3)]))

with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nVariability (SAME TRAINING SAMPLES):\n')
    file.write(f'\n|T_1|      = {len(output_data_template_1)}')
    file.write(f'\n|T_1_run2| = {len(output_data_template_1_run2)}')
    file.write(f'\n|T_1_run3| = {len(output_data_template_1_run3)}')
    file.write(f'\n|T_same_s| = {len(union_set_same_samples)}')
    file.write(f'\nVariability = {variability_same_samples}')

  
#####################
#     TIME          #  
#####################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nComputational time (seconds):\n')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
    file.write(f'\ntime_learning_2 = {time_learning_2}')
    file.write(f'\ntime_learning_3 = {time_learning_3}')
    file.write(f'\n\nAverage = {np.mean([time_learning_1, time_learning_2, time_learning_3])} +- {np.std([time_learning_1, time_learning_2, time_learning_3])}')
    file.write(f'\ntime_learning_1_run2 = {time_learning_1_run2}')
    file.write(f'\ntime_learning_1_run3 = {time_learning_1_run3}')
    
    file.write(f'\n\nTotal_time = {end_time_learning_1_run3 - start_time_learning_1}     (~{(end_time_learning_1_run3 - start_time_learning_1)/3600} hours)')
    
######################
# SAVE LIST OF SEEDS #  
######################

with open(f'{file_name}','a') as file:
    file.write('\n\n\n\n\n\nList of seeds used:\n')
    file.write(f'\n\nseed_learning_1 = {seed_learning_1}\nseed_negative_test_1 = {seed_negative_test_1}\nseed_learning_template_1 ={seed_learning_template_1}')
    file.write(f'\n\nseed_learning_2 = {seed_learning_2}\nseed_negative_test_2 = {seed_negative_test_2}\nseed_learning_template_2 ={seed_learning_template_2}')
    file.write(f'\n\nseed_learning_3 = {seed_learning_3}\nseed_negative_test_3 = {seed_negative_test_3}\nseed_learning_template_3 ={seed_learning_template_3}')
    file.write(f'\n\nseed_learning_template_1_run2 = {seed_learning_template_1_run2}\nseed_learning_template_1_run3 = {seed_learning_template_1_run3}')  
    



##############################################################################################################
###################### UNRESTRICTED SEARCH ########################################################################################
################################################################################################################################################################

#####################
#First batch of data#
#####################

union_set_different_samples = set()
union_set_same_samples = set()
set_U1 = set()
set_U2 = set()
set_U3 = set()
set_U1_run2 = set()
set_U1_run3 = set()

with open(f'{file_name}','a') as file:
    file.write('\n\n\n\n######################################################################################################')
    file.write('\nUNRESTRICTED SEARCH \n\n\nExperiment 1:\n\n')

start_time_learning_1 = time.time()
seed_learning_free_1 = 106954 #random.randint(0, 1000000)

n = len(learning_traces_1[0])  # number of variables

# Learn free
# Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = [ ' ( pi[0][0] < epsilon0- )' , ' ( pi[0][1] > epsilon1- )' , 
                     ' ( pi[1][0] < epsilon2- )' , ' ( pi[1][1] > epsilon3- )' , 
                     ' ( pi[2][0] < epsilon4- )' , ' ( pi[2][1] > epsilon5- )' , 
                     '( abs( pi[0][0] - pi[1][0] ) < epsilon6- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon7- )',
                     ' (( pi[0][0] + pi[1][0] + pi[2][0] ) > epsilon8- )', ' (( pi[0][1] + pi[1][1] + pi[2][1] ) < epsilon9- )' ]


input_data_free_1 = inputs.InputData(
    seed= seed_learning_free_1,
    traces=learning_traces_1,
    numb_quantifiers=[1, 3],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_free_1.learning_stl = True
input_data_free_1.par_bounds = [ [0.1, 50], [0.1, 50] , [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50]]
input_data_free_1.fraction = fraction
input_data_free_1.kind_score = 'efficient_monitoring_stl'
    
input_data_free_1.store_enumeration_formulas = False
input_data_free_1.store_formulas = True
input_data_free_1.n_target = numb_formulas
input_data_free_1.name_set_mined_formulas = 'free_'

output_data_free_1 = syntax_guided.main(input_data_free_1)

nb_false_positives_1 = 0
nb_false_negatives_1 = 0
nb_false_positives_1_runout = 0
nb_false_negatives_1_runout = 0

for index, item in enumerate(output_data_free_1):
    positive_result = se.efficient_monitoring(positive_test_traces_1, tree.tree_to_formula(item.nodes), item.quantifiers,
                            True, None, False, [], [])
    negative_result = se.efficient_monitoring(negative_test_traces_1, tree.tree_to_formula(item.nodes),
                                              item.quantifiers,
                                              True, None, False, [], [])
    
    if positive_result == 'runout' : nb_false_positives_1_runout += 1
    elif positive_result == -1: nb_false_negatives_1 = nb_false_negatives_1 + 1
    
    if negative_result == 'runout': nb_false_negatives_1_runout += 1
    elif negative_result == 1: nb_false_positives_1 = nb_false_positives_1 + 1
        
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    
    union_set_different_samples = union_set_different_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    union_set_same_samples = union_set_same_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_U1 = set_U1.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    
ratio_false_negatives_1 = nb_false_negatives_1/(len(output_data_free_1))
ratio_false_positives_1 = nb_false_positives_1/(len(output_data_free_1))
        
end_time_learning_1 = time.time()  
time_learning_1 = end_time_learning_1 - start_time_learning_1     

  

########################      
# Second batch of data #
########################

with open(f'{file_name}','a') as file: file.write('\n\nExperiment 2:\n\n')

start_time_learning_2 = time.time()

seed_learning_free_2 = 467101#random.randint(0, 1000000)

n = len(learning_traces_2[0])  # number of variables

# Learn free
# Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = [ ' ( pi[0][0] < epsilon0- )' , ' ( pi[0][1] > epsilon1- )' , 
                     ' ( pi[1][0] < epsilon2- )' , ' ( pi[1][1] > epsilon3- )' , 
                     ' ( pi[2][0] < epsilon4- )' , ' ( pi[2][1] > epsilon5- )' , 
                     '( abs( pi[0][0] - pi[1][0] ) < epsilon6- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon7- )',
                     ' (( pi[0][0] + pi[1][0] + pi[2][0] ) > epsilon8- )', ' (( pi[0][1] + pi[1][1] + pi[2][1] ) < epsilon9- )' ]


input_data_free_2 = inputs.InputData(
    seed= seed_learning_free_2 ,
    traces=learning_traces_2,
    numb_quantifiers=[1, 3], 
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_free_2.learning_stl = True
input_data_free_2.par_bounds = [ [0.1, 50], [0.1, 50] , [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50]]
input_data_free_2.fraction = fraction
input_data_free_2.kind_score = 'efficient_monitoring_stl'
    
input_data_free_2.store_enumeration_formulas = True
input_data_free_2.store_formulas = True
input_data_free_2.n_target = numb_formulas
input_data_free_2.name_set_mined_formulas = 'free_'
    

output_data_free_2 = syntax_guided.main(input_data_free_2)

nb_false_positives_2 = 0
nb_false_negatives_2 = 0
nb_false_positives_2_runout = 0
nb_false_negatives_2_runout = 0
for index, item in enumerate(output_data_free_2):
    positive_result = se.efficient_monitoring(positive_test_traces_2, tree.tree_to_formula(item.nodes), item.quantifiers,
                            True, None, False, [], [])
    negative_result = se.efficient_monitoring(negative_test_traces_2, tree.tree_to_formula(item.nodes),
                                              item.quantifiers,
                                              True, None, False, [], [])
    
    if positive_result == 'runout': nb_false_negatives_2_runout +=1
    elif positive_result == -1:  nb_false_negatives_2 = nb_false_negatives_2 + 1
    
    if negative_result == 'runout': nb_false_positives_2_runout +=1 
    elif negative_result == 1: nb_false_positives_2 = nb_false_positives_2 + 1

    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    union_set_different_samples = union_set_different_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_U2 = set_U2.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    
    
ratio_false_negatives_2 = nb_false_negatives_2/(len(output_data_free_2))
ratio_false_positives_2 = nb_false_positives_2/(len(output_data_free_2))
end_time_learning_2 = time.time()  
time_learning_2 = end_time_learning_2 - start_time_learning_2    
   
 
#####################   
#Third batch of data#
#####################


with open(f'{file_name}','a') as file: file.write('\n\nExperiment 3:\n\n')

start_time_learning_3 = time.time()

seed_learning_free_3 = 861702#random.randint(0, 1000000)

n = len(learning_traces_3[0])  # number of variables

# Learn ree

# Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = [ ' ( pi[0][0] < epsilon0- )' , ' ( pi[0][1] > epsilon1- )' , 
                     ' ( pi[1][0] < epsilon2- )' , ' ( pi[1][1] > epsilon3- )' , 
                     ' ( pi[2][0] < epsilon4- )' , ' ( pi[2][1] > epsilon5- )' , 
                     '( abs( pi[0][0] - pi[1][0] ) < epsilon6- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon7- )',
                     ' (( pi[0][0] + pi[1][0] + pi[2][0] ) > epsilon8- )', ' (( pi[0][1] + pi[1][1] + pi[2][1] ) < epsilon9- )' ]


input_data_free_3 = inputs.InputData(
    seed= seed_learning_free_3 ,
    traces=learning_traces_3,
    numb_quantifiers=[1, 3],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_free_3.learning_stl = True
input_data_free_3.par_bounds = [ [0.1, 50], [0.1, 50] , [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50]]
input_data_free_3.fraction = fraction
input_data_free_3.kind_score = 'efficient_monitoring_stl'
    
input_data_free_3.store_enumeration_formulas = True
input_data_free_3.store_formulas = True
input_data_free_3.n_target = numb_formulas
input_data_free_3.name_set_mined_formulas = 'free_'

output_data_free_3 = syntax_guided.main(input_data_free_3)

nb_false_positives_3 = 0
nb_false_negatives_3 = 0
nb_false_positives_3_runout = 0
nb_false_negatives_3_runout = 0

for index, item in enumerate(output_data_free_3):
    positive_result = se.efficient_monitoring(positive_test_traces_3, tree.tree_to_formula(item.nodes), item.quantifiers,
                            True, None, False, [], [])
    negative_result = se.efficient_monitoring(negative_test_traces_3, tree.tree_to_formula(item.nodes),
                                              item.quantifiers,
                                              True, None, False, [], [])
    
    if positive_result == 'runout': nb_false_negatives_3_runout += 1
    elif positive_result == -1:  nb_false_negatives_3 = nb_false_negatives_3 + 1
    
    if negative_result == 'runout': nb_false_positives_3_runout += 1
    elif negative_result == 1: nb_false_positives_3 = nb_false_positives_3 + 1
    
    with open(f'{file_name}','a') as file: file.write(f'\n\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    
    union_set_different_samples = union_set_different_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_U3 = set_U3.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    
ratio_false_negatives_3 = nb_false_negatives_3/(len(output_data_free_3))
ratio_false_positives_3 = nb_false_positives_3/(len(output_data_free_3))

end_time_learning_3 = time.time()  
time_learning_3 = end_time_learning_3 - start_time_learning_3   


    
##################################################
#     VARIABILITY DIFFERENT TRAINING SAMPLES     #
##################################################


variability_different_samples = (len(union_set_different_samples.difference(set_U1)) + len(union_set_different_samples.difference(set_U2)) + len(union_set_different_samples.difference(set_U3)) )/( 6 * np.mean([len(set_U1), len(set_U2), len(set_U3)]))
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nVariability (DIFFERENT TRAINING SAMPLES):\n')
    file.write(f'\n|U_1| = {len(output_data_free_1)}')
    file.write(f'\n|U_2| = {len(output_data_free_2)}')
    file.write(f'\n|U_3| = {len(output_data_free_3)}')
    file.write(f'\n|U|   = {len(union_set_different_samples)}')
    file.write(f'\nVariability   = {variability_different_samples}')



################################
#   RATIO OF FALSE POSITIVES   #       
################################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nRatio of false positives:\n')
    file.write(f'\nratio_false_positives_1 = {ratio_false_positives_1}')
    if nb_false_positives_1_runout != 0 : file.write(f'runout_false_positives_1 = {nb_false_positives_1_runout}')
    
    file.write(f'\nratio_false_positives_2 = {ratio_false_positives_2}')
    if nb_false_positives_2_runout != 0 : file.write(f'runout_false_positives_2 = {nb_false_positives_2_runout}')
    
    file.write(f'\nratio_false_positives_3 = {ratio_false_positives_3}')
    if nb_false_positives_3_runout != 0 : file.write(f'runout_false_positives_3 = {nb_false_positives_3_runout}')
    
    file.write(f'\n\nAverage = {np.mean([ratio_false_positives_1, ratio_false_positives_2, ratio_false_positives_3])} +- {np.std([ratio_false_positives_1, ratio_false_positives_2, ratio_false_positives_3])}')
  
    
################################
#   RATIO OF FALSE NEGATIVES   #       
################################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nRatio of false negatives:\n')
    file.write(f'\nratio_false_negatives_1 = {ratio_false_negatives_1}')
    if nb_false_negatives_1_runout != 0 : file.write(f'runout_false_negatives_1 = {nb_false_negatives_1_runout}')
    
    file.write(f'\nratio_false_negatives_2 = {ratio_false_negatives_2}')
    if nb_false_negatives_2_runout != 0 : file.write(f'runout_false_negatives_2 = {nb_false_negatives_2_runout}')
    
    file.write(f'\nratio_false_negatives_3 = {ratio_false_negatives_3}')
    if nb_false_negatives_3_runout != 0 : file.write(f'runout_false_negatives_3 = {nb_false_negatives_3_runout}')
    
    file.write(f'\n\nAverage = {np.mean([ratio_false_negatives_1, ratio_false_negatives_2, ratio_false_negatives_3])} +- {np.std([ratio_false_negatives_1, ratio_false_negatives_2, ratio_false_negatives_3])}')
  

    
#########################################    
# Second run with the same batch of data#    
#########################################

with open(f'{file_name}','a') as file: file.write('\n\Second run with the same traces as Experiment 1:\n\n')
seed_learning_free_1_run2 = 664015#random.randint(0, 1000000)

start_time_learning_1_run2 = time.time() 

n = len(learning_traces_1[0])  # number of variables

# Learn free
# Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''

# Definition of the grammar for  predicates
grammar_predicate = [ ' ( pi[0][0] < epsilon0- )' , ' ( pi[0][1] > epsilon1- )' , 
                     ' ( pi[1][0] < epsilon2- )' , ' ( pi[1][1] > epsilon3- )' , 
                     ' ( pi[2][0] < epsilon4- )' , ' ( pi[2][1] > epsilon5- )' , 
                     '( abs( pi[0][0] - pi[1][0] ) < epsilon6- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon7- )',
                     ' (( pi[0][0] + pi[1][0] + pi[2][0] ) > epsilon8- )', ' (( pi[0][1] + pi[1][1] + pi[2][1] ) < epsilon9- )' ]

input_data_free_1_run2 = inputs.InputData(
    seed = seed_learning_free_1_run2 ,
    traces=learning_traces_1,
    numb_quantifiers=[1, 3],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_free_1_run2.learning_stl = True
input_data_free_1_run2.par_bounds = [ [0.1, 50], [0.1, 50] , [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50]]
input_data_free_1_run2.fraction = fraction
input_data_free_1_run2.kind_score = 'efficient_monitoring_stl'
    
input_data_free_1_run2.store_enumeration_formulas = True
input_data_free_1_run2.store_formulas = True
input_data_free_1_run2.n_target = numb_formulas
input_data_free_1_run2.name_set_mined_formulas = 'free_'    

output_data_free_1_run2 = syntax_guided.main(input_data_free_1_run2)

for index, item in enumerate(output_data_free_1_run2):
   
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    union_set_same_samples = union_set_same_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_U1_run2 = set_U1_run2.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
      
end_time_learning_1_run2 = time.time()  
time_learning_1_run2 = end_time_learning_1_run2 - start_time_learning_1_run2    
    
 
#########################################    
# Third run with the same batch of data#    
#########################################
seed_learning_free_1_run3 = 556821#random.randint(0, 1000000)
with open(f'{file_name}','a') as file: file.write('\n\Third run with the same traces as Experiment 1:\n\n')


start_time_learning_1_run3 = time.time() 

n = len(learning_traces_1[0])  # number of variables

# Learn free
# Definition of the grammar for the quantifiers
grammar_quantifiers = '''
psi0 == forall psi0 | exists psi0 | phi0
'''
# Definition of the grammar for the formula structure
grammar_structure = '''
phi0 == always phi0 | eventually phi0 | phi0 until phi0 | phi0 and phi0 | phi0 or phi0 | phi0 implies phi0 | not phi0 | predicate0
'''

#Definition of the grammar for  predicates
grammar_predicate = [ ' ( pi[0][0] < epsilon0- )' , ' ( pi[0][1] > epsilon1- )' , 
                     ' ( pi[1][0] < epsilon2- )' , ' ( pi[1][1] > epsilon3- )' , 
                     ' ( pi[2][0] < epsilon4- )' , ' ( pi[2][1] > epsilon5- )' , 
                     '( abs( pi[0][0] - pi[1][0] ) < epsilon6- )', '( abs( pi[0][1] - pi[1][1] ) < epsilon7- )',
                     ' (( pi[0][0] + pi[1][0] + pi[2][0] ) > epsilon8- )', ' (( pi[0][1] + pi[1][1] + pi[2][1] ) < epsilon9- )' ]

input_data_free_1_run3 = inputs.InputData(
    seed = seed_learning_free_1_run3 ,
    traces=learning_traces_1,
    numb_quantifiers=[1, 3],  
    grammar_quantifiers=grammar_quantifiers,
    grammar_structure=grammar_structure,
    grammar_predicate=grammar_predicate,
    max_length=6)

# Additional options:
input_data_free_1_run3.learning_stl = True
input_data_free_1_run3.par_bounds = [ [0.1, 50], [0.1, 50] , [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50], [0.1, 50]]
input_data_free_1_run3.fraction = fraction
input_data_free_1_run3.kind_score = 'efficient_monitoring_stl'
    
input_data_free_1_run3.store_enumeration_formulas = True
input_data_free_1_run3.store_formulas = True
input_data_free_1_run3.n_target = numb_formulas
input_data_free_1_run3.name_set_mined_formulas = 'free_'        
    

output_data_free_1_run3 = syntax_guided.main(input_data_free_1_run3)

for index, item in enumerate(output_data_free_1_run3):
   
    with open(f'{file_name}','a') as file: file.write(f'\nFormula {index+1}: {syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))};\n\n')
    union_set_same_samples = union_set_same_samples.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))]))
    set_U1_run3 = set_U1_run3.union(set([ syntax_guided.print_formula(item.quantifiers, tree.NATURAL_tree_to_formula(item.nodes))])) 
    
end_time_learning_1_run3 = time.time()  
time_learning_1_run3 = end_time_learning_1_run3 - start_time_learning_1_run3    
 
  
##################################################
#     VARIABILITY SAME TRAINING SAMPLES          #
##################################################
variability_same_samples = (len(union_set_same_samples.difference(set_U1)) + len(union_set_same_samples.difference(set_U1_run2)) + len(union_set_same_samples.difference(set_U1_run3)) )/( 6 * np.mean([len(set_U1), len(set_U1_run2), len(set_U1_run3)]))

with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nVariability (SAME TRAINING SAMPLES):\n')
    file.write(f'\n|U_1|      = {len(output_data_free_1)}')
    file.write(f'\n|U_1_run2| = {len(output_data_free_1_run2)}')
    file.write(f'\n|U_1_run3| = {len(output_data_free_1_run3)}')
    file.write(f'\n|U_same_s| = {len(union_set_same_samples)}')
    file.write(f'\nVariability = {variability_same_samples}')

  
#####################
#     TIME          #  
#####################
with open(f'{file_name}','a') as file:
    file.write('\n\n\n\nComputational time (seconds):\n')
    file.write(f'\ntime_learning_1 = {time_learning_1}')
    file.write(f'\ntime_learning_2 = {time_learning_2}')
    file.write(f'\ntime_learning_3 = {time_learning_3}')
    file.write(f'\n\nAverage = {np.mean([time_learning_1, time_learning_2, time_learning_3])} +- {np.std([time_learning_1, time_learning_2, time_learning_3])}')
    file.write(f'\ntime_learning_1_run2 = {time_learning_1_run2}')
    file.write(f'\ntime_learning_1_run3 = {time_learning_1_run3}')
    
    file.write(f'\n\nTotal_time = {end_time_learning_1_run3 - start_time_learning_1}     (~{(end_time_learning_1_run3 - start_time_learning_1)/3600} hours)')


    
######################
# SAVE LIST OF SEEDS #  
######################

with open(f'{file_name}','a') as file:
    file.write('\n\n\n\n\n\nList of seeds used:\n')
    file.write('\nFor the data: same as the template setting.\nFor the learning process:')
    file.write(f'\nseed_learning_free_1 = {seed_learning_free_1}\nseed_learning_free_2 = {seed_learning_free_2}\nseed_learning_free_3 = {seed_learning_free_3}')  
    file.write(f'\nseed_learning_free_1_run2 = {seed_learning_free_1_run2}\nseed_learning_free_1_run3 = {seed_learning_free_1_run3}')  
    

    
        
    