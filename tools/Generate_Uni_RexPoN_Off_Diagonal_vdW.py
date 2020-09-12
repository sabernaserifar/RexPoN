#!/usr/bin/env python
# coding: utf-8

# In[7]:


from math import sqrt 

# Atoms should be provided with the same order as in the ffied.RexPoN
# The parameters are: De, Re, and L
# Use Table II of "J. Chem. Phys. 151, 154111 (2019); doi: 10.1063/1.5113811" for the above parameters
# If not available use UFF values 


atom_dic = {'Hw':[0.0528, 3.2541, 0.5241], 
            'Ow':[0.1498, 3.4249, 0.4349],
            'Cl':[0.3500, 4.0748, 0.5654],
            'Na':[0.0300, 2.9830, 0.5654],
           }





# RexPoN water parameters (never change!)
print('1 1  0.0212 2.4002 10.5460 0.7450 -1.0000 -1.0000 -0.8113 0.0014 -0.2607 -0.5640 0.0 1.6346 239.4803')
print('1 2  0.0906 3.1455 10.0823 0.9587 -1.0000 -1.0000 -0.0007 -1.0000 0.0000 0.0000 0.0 1.5642 389.8714')
print('2 2  0.0272 3.6468 14.5278 0.8100 -1.0000 -1.0000 -0.3235 0.0305 0.0331 -0.5688 0.0 1.4969 634.7066')
# print the vdW section of the RexPoN ffield
i = 1
j = 1
number_of_lines = 0
for atomi in atom_dic:
    j = 1
    for atomj in atom_dic:
        #print(atomi, atomj)
        if (j >= i) and not ( (atomi == 'Ow' or atomi == 'Hw') and (atomj == 'Hw' or atomj == 'Ow')) :
            De = sqrt(atom_dic[atomi][0] * atom_dic[atomj][0] )
            Re = sqrt(atom_dic[atomi][1] * atom_dic[atomj][1] )
            L = sqrt(atom_dic[atomi][2] * atom_dic[atomj][2] )
            print("%i %i %0.5f %0.5f %0.5f %s "%(i, j, De, Re, L, tail))
            number_of_lines += 1
        j += 1
    i += 1
print("Number of lines = ", number_of_lines+3)

