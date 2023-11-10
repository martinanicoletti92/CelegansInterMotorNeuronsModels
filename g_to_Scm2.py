def gScm2(g,surf,index):


    import math
    import numpy
    
    g_RESCALED=[]
    
    
    for i in range(len(g)):
        if i<=index: 
         gg=(g[i]*1e-9)/(surf)
        else:
         gg=g[i]
         
        g_RESCALED.append(gg)

    g_RESCALED=numpy.array(g_RESCALED)
    return g_RESCALED