
def gnS(g,surf,index):
# index is the position in theg vector  of the last conductance to be optimized
# entries from index+1 to the end are kinetic parameters that can be optimized by the algotithm


    import math
    import numpy

       
    g_RESCALED=[]
    
    
    for i in range(len(g)):
        if i<=index:
         gg=(g[i]*surf)/(1e-9)
        else: 
         gg=g[i]
         
    g_RESCALED.append(gg)
    
    g_RESCALED=numpy.array(g_RESCALED)
    
    return g_RESCALED