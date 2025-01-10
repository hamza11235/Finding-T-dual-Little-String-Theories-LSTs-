from sage.all import *
load("utils.sage")
from matplotlib import pyplot as plt


    

def verts(vertices: list[[list]]):
    
    """
    Generating 4D reflexive polytope
    """
    
    p=get_polytope(vertices)
    return p



def pick_fib(p: LatticePolytope, n):
    
    """
    Returns: 

    1. Total number of total 2D subpolytopes, 
    2. A particular 2D subpolytope n, along with projection vectors for the base of the corresponding fibration, and a plot of the        2D subploytope inside the 3D polytope.
    """
    
    K3,num=get_3D_2D_reflexive(p)
    torus,xproj,yproj,plot=pick_fibration(K3, n)
    return num,torus,xproj,yproj,plot


def CY_data(p: LatticePolytope, torus,xproj,yproj):
    
    """
    Returns: 

    1. CY hypersurface equation from the 4D polytope, using Batyrev's construction, with ineffective divisors and non-flat flavor        divisors removed. 
    2. Corrected vertices and coordinates (torus coordinates distinguished).
    3. Base of the fibration.
    """
    
    p,coords, weirs,letvar,dualp,arrcoeff=CY_Hypersurface(p, torus)
    baseindex1=0
    baseindex2=1
    base=get_2Dbase(p, xproj, yproj)
    letvar,model,p,base,coords,nfindex,nbcoords,nbdowns=rem_nonflat_flav(p, coords,  letvar,  base, baseindex1,baseindex2,weirs)
    weirs=sym.expand(batryev(p,dualp,arrcoeff,coords))
    base=get_2Dbase(p, xproj, yproj)
    return p,letvar, coords, baseindex1, baseindex2, weirs, base



def base_dat(base, baseindex1, baseindex2):
    
    """
    Returns base divisors and a plot of the 2D base.
    """
    
    basedivs=getbasedivs(base,baseindex1,baseindex2)
    xc=list()
    yc=list()
    for j in range(len(basedivs)):
        xc.append(basedivs[j][0])
        yc.append(basedivs[j][1])
    for i in range(len(basedivs)):
        x = xc[i]
        y = yc[i]
        plt.scatter(x, y, marker = 'o')

        # V this adds the lines V
        plt.plot([0,xc[i]], [0, yc[i]], color="black")
        # ^ this adds the lines ^

        plt.plot(x,y,0,0, linestyle = '--' )
    return basedivs

    
    
    
    
        
        
        