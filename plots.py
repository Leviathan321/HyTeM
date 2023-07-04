#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pickle
import matplotlib.pyplot as plt
import numpy as np
import math
import random

'Plot trigonometric function'

# value = np.arange(-0.999,1, 0.001)
  
# plt.plot(value, value, label = r'$f(x) = x$' )
# plt.plot(value,[ math.cos(i) for i in value], label = r'$f(x) = \cos(x)$' )
# plt.plot(value,[ math.sin(i)for i in value], label = r'$f(x) = \sin(x)$' )
# plt.plot(value,[ math.acos(i)for i in value], label = r'$f(x) = \arccos(x)$' )
# plt.plot(value,[ math.asin(i)for i in value], label = r'$f(x) = \arcsin(x)$' )

# plt.xlabel(r'$x$', fontsize = 50)    
# plt.ylabel(f'$f(x)$', fontsize = 50) 
# plt.xticks(fontsize=30)    
# plt.yticks(fontsize=30) 
# plt.legend(fontsize = 28)
# plt.savefig(f"trig.pdf", format='pdf', bbox_inches="tight") 




# # 'Plot time of cost vs boolean satisfaction '
# ##Constrained (FF, and, always, implies) Serial Adder with 50 traces , 100 times
# path = 'Learned/SerialAdderwithCarryOn/cost vs boolean satisfaction/cost_constrained_50traces_100times/'

# database = []
# for i in range(1,51):
#     data = f'formula_to_be_stored{i}.obj'
#     with open(f'{path}{data}', 'rb') as f:
#               read = pickle.load(f)
#               database.append(read)
             
# cost_time = []
# aux = 0
# for i in range(50):
#     aux += (database[i].time/60)
#     cost_time.append(aux)
    
# path = 'Learned/SerialAdderwithCarryOn/cost vs boolean satisfaction/Satisfaction_constrained_50traces_100times/'


# database = []
# for i in range(1,51):
#     data = f'formula_to_be_stored{i}.obj'
#     with open(f'{path}{data}', 'rb') as f:
#               read = pickle.load(f)
#               database.append(read)
             
# satisfaction_time = []
# aux = 0
# for i in range(50):
#     aux += (database[i].time/60)
#     satisfaction_time.append(aux)

# value_x = np.arange(1,51)



# plt.plot(value_x,satisfaction_time, 'o-', label = 'Effective Monitoring' )
# plt.plot(value_x,cost_time, 'o-', label = 'Fitness function')
# plt.xlabel('Number of learned formulas', fontsize = 50)    
# plt.ylabel('Time', fontsize = 50) 
# plt.xticks(fontsize=50)    
# plt.yticks(fontsize=50) 
# plt.legend(fontsize = 30)
# plt.savefig(f"{path}plot.pdf", format='pdf', bbox_inches="tight")  



    
# 'Plot time using fixed length'
# ## No cost, ONLY BOOLEAN SATISFACTION with early stopping 
# ##Constrained (FF, and, always, implies) Serial Adder with 50 traces , 100 times


# path = 'Learned/SerialAdderwithCarryOn/Length_free/'#constrained/'

# lengths = [3, 4, 5 , 6, 7]
# time_lengths_2quant = []

# for l in range(3,8):
#     data = f'2quant_len{l}_free_sat/' #constrained_sat/
    
#     database = []
#     for i in range(1,51):
#         formula = f'formula_to_be_stored{i}.obj'
#         with open(f'{path}{data}{formula}', 'rb') as f:
#                   read = pickle.load(f)
#                   database.append(read)
                 
#     aux = []
#     for i in range(50): aux.append(database[i].time)
#     time_lengths_2quant.append(np.mean(aux))  
    
# path = 'Learned/SerialAdderwithCarryOn/Length_constrained/'#constrained/'

# time_lengths_3quant = []

# for l in range(3,8):
#     data = f'2quant_len{l}_constrained_sat/'#constrained_sat/
    
#     database = []
#     for i in range(1,51):
#         formula = f'formula_to_be_stored{i}.obj'
#         with open(f'{path}{data}{formula}', 'rb') as f:
#                   read = pickle.load(f)
#                   database.append(read)
                 
#     aux = []
#     for i in range(50): aux.append(database[i].time)
#     time_lengths_3quant.append(np.mean(aux))  
  
# time_template = []
# path = 'Learned/SerialAdderwithCarryOn/'


# data = '2quant_template_sat/'

# database = []
# for i in range(1,51):
#     formula = f'formula_to_be_stored{i}.obj'
#     with open(f'{path}{data}{formula}', 'rb') as f:
#               read = pickle.load(f)
#               database.append(read)
             
# aux = []
# for i in range(50): aux.append(database[i].time)
# time_template.append(np.mean(aux))  
  
  
# X_axis = np.arange(len(lengths))  
  
# # plt.bar(X_axis - 0.2, time_lengths_2quant, width = 0.4, label = '2 quantifiers')
# # plt.bar(X_axis + 0.2, time_lengths_3quant , width = 0.4, label = '3 quantifiers')

# plt.bar(X_axis - 0.4, time_lengths_2quant, width = 0.3, label = 'Free')
# plt.bar(X_axis - 0.1, time_lengths_3quant , width = 0.3, label = 'Constrained')

# plt.bar( [4 + 0.2] , time_template , width = 0.3, label = 'Template')



# plt.xticks([-0.2, 0.8, 1.8, 2.8, 4], lengths, fontsize=50) 
# plt.xlabel("Formula length", fontsize = 50)
# plt.ylabel("Time", fontsize = 50)   
    
# plt.yticks(fontsize=50) 
# plt.legend(fontsize = 30)
# plt.savefig(f"{path}final.pdf", format='pdf', bbox_inches="tight")  



# lengths = [r'$ \forall \forall $ ', r'$ \forall \exists $', r'$ \exists \forall $' , r'$ \exists \exists $']
  
# X_axis = np.arange(len(lengths))  
 
# #Computed using effective monitoring and correctness theorem
# corr = [ 3976,  54,  52, 3404 ]

# #Computed with just effective monitoring (no correctness theorem)
# eff_mon = [27249, 185, 184, 3404   ] 

# # ev= [20441, 117, 119, 203]
# # skip = [13598, 67, 67, 1]
# # plt.bar(X_axis - 0.2, time_lengths_2quant, width = 0.4, label = '2 quantifiers')
# # plt.bar(X_axis + 0.2, time_lengths_3quant , width = 0.4, label = '3 quantifiers')
# plt.bar(X_axis - 0.4, eff_mon , width = 0.3, label = 'Effective monitoring')
# plt.bar(X_axis - 0.1,  corr, width = 0.3, label = 'Effective monitoring and correctness theorem')

# plt.xticks([-0.2, 0.8, 1.8, 2.8], lengths, fontsize=50) 
# # plt.xlabel("Formula length", fontsize = 50)
# # plt.ylabel("Time", fontsize = 50)   
    
# plt.yticks(fontsize=50) 
# plt.legend(fontsize = 30)
# plt.savefig(f"plot_all.pdf", format='pdf', bbox_inches="tight")  
# plt.close()

lengths = [r'$ \forall \forall $ ', r'$ \forall \exists $', r'$ \exists \forall $' , r'$ \exists \exists $']
  
X_axis = np.arange(len(lengths))  
 
#Computed using effective monitoring and correctness theorem
corr = [ 3.59, 1.73, 1.71, 3.53]

#Computed with just effective monitoring (no correctness theorem)
eff_mon = [ 4.43, 2.26, 2.26, 3.53 ] 

al = [ 4.53, 4.53, 4.53, 4.53]

# ev= [20441, 117, 119, 203]
# skip = [13598, 67, 67, 1]
# plt.bar(X_axis - 0.2, time_lengths_2quant, width = 0.4, label = '2 quantifiers')
# plt.bar(X_axis + 0.2, time_lengths_3quant , width = 0.4, label = '3 quantifiers')

plt.bar(X_axis - 0.4,  al , width = 0.2, label = 'Fitness', color = 'grey')
plt.bar(X_axis - 0.2, eff_mon , width = 0.2, label = 'Effective  monitoring', color = 'darkmagenta')
plt.bar(X_axis - 0,  corr, width = 0.2, label = 'Effective  monitoring + correctness', color = 'sandybrown')

plt.xticks([-0.2, 0.8, 1.8, 2.8], lengths, fontsize=50) 
# plt.xlabel("Formula length", fontsize = 50)
# plt.ylabel("Time", fontsize = 50)     
plt.yticks(fontsize=50) 
plt.legend(fontsize = 25,  bbox_to_anchor=(0.3, 0.4, 0.5, 0.5))
plt.savefig(f"plot.pdf", format='pdf', bbox_inches="tight")  
plt.close()

#Computed using effective monitoring and correctness theorem
# corr = [ 8.38,  3.98, 3.95, 8.13 ]

# #Computed with just effective monitoring (no correctness theorem)
# eff_mon = [10.21, 5.22, 5.21, 8.13   ] 

# plt.bar(X_axis - 0.4, eff_mon , width = 0.3, label = 'Effective monitoring')

# plt.bar(X_axis - 0.1,  corr, width = 0.3, label = 'Effective monitoring and correctness theorem')

# plt.xticks([-0.2, 0.8, 1.8, 2.8], lengths, fontsize=50)    
# plt.yticks(fontsize=50) 
# plt.legend(fontsize = 30)
# plt.savefig(f"plot_log.pdf", format='pdf', bbox_inches="tight")  
# plt.close()


# lengths = [r'$ \forall \forall $ ', r'$ \exists \exists $'] 
# X_axis = np.arange(len(lengths))  
# #Computed using effective monitoring and correctness theorem
# corr = [ 3976,  3404 ]
# #Computed with just effective monitoring (no correctness theorem)
# eff_mon = [27249, 3404   ] 
# plt.bar(X_axis - 0.4, eff_mon , width = 0.3, label = 'Effective monitoring')
# plt.bar(X_axis - 0.1,  corr, width = 0.3, label = 'Effective monitoring and correctness theorem')
# plt.xticks([-0.2, 0.8], lengths, fontsize=50) 
# plt.yticks(fontsize=50) 
# plt.legend(fontsize = 30)
# plt.savefig(f"plot_ff_ee.pdf", format='pdf', bbox_inches="tight")  
# plt.close()


# lengths = [ r'$ \forall \exists $', r'$ \exists \forall $' ] 
# X_axis = np.arange(len(lengths))  
# #Computed using effective monitoring and correctness theorem
# corr = [  54,  52 ]
# #Computed with just effective monitoring (no correctness theorem)
# eff_mon = [ 185, 184  ] 
# plt.bar(X_axis - 0.4, eff_mon , width = 0.3, label = 'Effective monitoring')

# plt.bar(X_axis - 0.1,  corr, width = 0.3, label = 'Effective monitoring and correctness theorem')
# plt.xticks([-0.2, 0.8], lengths, fontsize=50) 
# plt.yticks(fontsize=50) 
# plt.legend(fontsize = 30)
# plt.savefig(f"plot_fe_ef.pdf", format='pdf', bbox_inches="tight")  
# plt.close()

