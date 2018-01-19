#!/usr/bin/env python


Andrew W. Brooks
Vanderbilt Genetics Institute
andrew.w.brooks(at)vanderbilt.edu
------------------------------------------------------------------------------
Released under the MIT and Beerware Licenses
Copyright (c) 2018 Andrew W. Brooks
"MIT License" - Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the software or the use or other dealings in the software.
"Beerware License" - Andrew W. Brooks (2017) wrote this code. As long as you retain this notice and attribute my name, you can do whatever you want with this stuff. If we meet someday, and you think this stuff is worth it, you can buy me a beer in return.

### IMPORT PACKAGES ###
import os
import numpy as np
import pandas as pd
from ete3 import Tree
import random
import subprocess
random.seed(191919)  # Set Random Seed for Reproducibility


##### USER INPUT #####################################################################

### INPUT/OUTPUT FOLDER ###
outPath = 'peromyscus'

######################################################################################

### PATH TO THE HOST TREE ###
htPath = outPath+'/ht.newick'

### PATH TO MICROBIOME DENDROGRAM (OR ALTERNATIVE TREE) ###
mtPath = outPath+'/mt.newick'

### NUMBER OF RANDOM TREES ###
numRandom = 10000

##### IMPORT AND CLEAN TREES #########################################################

### IMPORT HOST TREE ###
htIn = Tree(htPath)
### IMPORT MICROBIOME TREE ###
mtIn = Tree(mtPath)

### PRINT THE TREES ###
print('\n --- HOST TREE --- ')
print(htIn)
print('\n --- MICROBIOME TREE --- ')
print(mtIn)

### LIST TO STORE HOST LEAVES ###
htLeafs = []

### CHECK THAT TREES HAVE SAME NUMBER OF LEAVES ###
print('\n --- CHECKING TREES --- \n')
if len(htIn) == len(mtIn):   
    
    ### FOR EACH LEAF IN THE HOST TREE MAKE SURE ALSO IN MICROBIOME TREE ###
    for leaf in mtIn:
        
        ### STORE LEAF NAME ###
        htLeafs.append(leaf.name)
        
        ### CHECK MICROBIOME TREE ###
        if leaf.name not in mtIn:
            print('ERROR - '+leaf.name+' Not Found in Microbiome Tree!')

### # IF NOT SAME NUMBER OF LEAVES PRINT ERROR ###
else: print('ERROR - Host and Microbiome Tree Have Different Numbers of Nodes!')

##### GENERATING RANDOM TREES ########################################################

### PATH TO STORE RANDOM TREES ###
rtPath = outPath+'/random_trees/'

### CREATE OUTPUT FOLDER ###
if not os.path.isdir(rtPath): os.makedirs(rtPath)

### FOR EACH RANDOM TREE ###
for curRandom in np.arange(numRandom):
    
    ### INTITALIZE RANDOM TREE ###
    t = Tree()
    
    ### SHUFFLE LEAVES ###
    random.shuffle(htLeafs)
    
    ### POPULATE TREE ###
    t.populate(len(htIn), names_library=htLeafs)
    
    ### WRITE TREE ###
    t.write(outfile=rtPath+'/tree_'+str(curRandom)+'.newick')

### CONCATENATE THE RANDOM TREES INTO SINGLE FILE ###
subprocess.check_output(str('for i in '+outPath+'/random_trees/*; do  cat $i >> '+outPath+'/random_trees.newick; done'), shell=True)

##### COMARE WITH TREECMP METHODS ####################################################

### FOR EACH TREECMP METHOD ###
methods = {"rc":'Rooted Robinson-Foulds', "mc":'Rooted Matching Cluster', "ms":'Unrooted Matching Split', "rf":'Unrooted Robinson-Foulds'}

for method in ["rc", "mc", "rf", "ms"]:
    
    ### MAKE COMPARISON FOLDER ###
    if not os.path.isdir(outPath+'/compare_'+method+'/'): os.makedirs(outPath+'/compare_'+method+'/')
        
    ### COMPARE HOST AND MICROBIOME TREE ###
    subprocess.check_output(str('java -jar TreeCmp/bin/TreeCmp.jar -r '+htPath+
              ' -d '+method+' -i '+mtPath+' -o '+outPath+'/compare_'+method+'/compare_ht_mt.txt -N'), shell=True)
    
    ### COMPARE HOST AND RANDOM TREES ###
    subprocess.check_output(str('java -jar TreeCmp/bin/TreeCmp.jar -r '+htPath+
              ' -d '+method+' -i '+outPath+'/random_trees.newick -o '+outPath+'/compare_'+method+'/compare_ht_random.txt -N'), shell=True)

    ### READ IN COMPARISON RESULTS ###
    mtCompare = pd.read_csv(outPath+'/compare_'+method+'/compare_ht_mt.txt', sep='\t')
    rtCompare = pd.read_csv(outPath+'/compare_'+method+'/compare_ht_random.txt', sep='\t')

    ### PRINT RESULTS ###
    print(' --- METHOD: '+methods[method]+' --- ')
    print("Host-Microbe Score:   " + str(mtCompare.iloc[0,1]))
    print("Max Stochastic Metric:  " + str(max(rtCompare.iloc[:,1])))
    print('Normalized Score: '+str(mtCompare.iloc[0,1]/max(rtCompare.iloc[:,1])))
    print('Random Trees with Equivalent or More Congruent Score: '+str(len(rtCompare[rtCompare.iloc[:,1] <= mtCompare.iloc[0,1]])))
    print('Total Trees: '+ str(len(rtCompare)))
    print('P-Value: '+str(len(rtCompare[rtCompare.iloc[:,1] <= mtCompare.iloc[0,1]]) / len(rtCompare))+'\n')
