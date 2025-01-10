import math
import numpy as np
import scipy as sci
import sympy as sym
from numpy.linalg import eig
from scipy.linalg import null_space
from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL



#Functions 


def get_polytope(vertices: list[[list]]):
    
    """
    Generate 4D Reflexive Polytope from a list of 4D vertices.

    Every such polytope is in a one-to-one correspondence with a toric Calabi-Yau (CY) Threefold.
    """
    
    P0 = Polyhedron(vertices)
    if P0.is_reflexive():
        print('All Good');
        Ptest = LatticePolytope(vertices);
        print(Ptest.poly_x("g"));
    else:
        PN = Polyhedron(vertices=P0.polar().integral_points())
        PM = PN.polar()
        print(PM.is_reflexive())
        Ptest=LatticePolytope(P0.polar().integral_points())
        print(Ptest.polar().poly_x("g"))
    return Ptest.polar()



def get_3D_2D_reflexive(p: LatticePolytope):
    
    """
    The class of CY's we are interested in admit both K3 and T2 (torus) fibrations, which means the 4D polytope has reflexive subpolytopes of dimensions 3 and 2. 

    Here we extract the 3D subpolytope, as well as the number of 2D subpolytopes (which gives us the number of T-dual LSTs)
    """
    
    K3=list()
    lemp=p.points()
    for i in range(len(lemp)):
        if lemp[i][3]==0:
            K3.append(lemp[i])
    templo=list()
    for i in range(len(K3)):
        templo.append([K3[i][0],K3[i][1],K3[i][2]])
    shaba=LatticePolytope(templo)
    shaba.points()
    blks=shaba.vertices()
    nk3=np.array(blks).tolist()
    pq = LatticePolytope_PPL(nk3)
    return shaba, len(list(pq.fibration_generator(2)))



def pick_fibration(shaba: LatticePolytope, fibind: int):
    
    """
    Extracting the T2 fibration (and hence the LST), the projection vectors for the 2D base, and a plot of the 2D subpolytope in the 3D polytope.
    """
    
    #Get polytope data
    blks=shaba.vertices()
    nk3=np.array(blks).tolist()
    pq = LatticePolytope_PPL(nk3)
    ct1=LatticePolytope(list(pq.fibration_generator(2))[fibind].integral_points()).plot3d(facet_color='red')
    nm1=LatticePolytope(shaba.points()).plot3d(facet_color='white')
    kaplz=LatticePolytope(list(pq.fibration_generator(2))[fibind].integral_points()).vertices()
    kmn=np.array(kaplz).tolist()
    oldnorm=np.cross(kmn[0], kmn[2])/(math.gcd(np.cross(kmn[0], kmn[2])[0],np.cross(kmn[0], kmn[2])[1],np.cross(kmn[0], kmn[2])[2]))
    comp1=oldnorm.tolist()
    comp1.append(1)
    comp2=[0,0,0,1]
    torusn=LatticePolytope(list(pq.fibration_generator(2))[fibind].integral_points()).vertices()
    torus=np.array(torusn).tolist()
    for i in range (len(torus)):
        torus[i].append(0)
    torus.append([0,0,0,0])
    plot=ct1+nm1
    return torus, comp1,comp2, ct1+nm1



def divgen(n,toruscoords):
    
    """
    Assigning variables to the points on the polytope (divisors), distinguishing the torus vertices.
    """
    
    var = [0] * n
    letvar = [0] * toruscoords
    lettersf=['X', 'Y', 'Z', 'L', 'M', 'N','P', 'Q']
    for i in range(toruscoords):
        var[i]=sym.symbols(lettersf[i])
        letvar[i]=lettersf[i]
    for i in range(toruscoords,n):
        var[i]=sym.expand(sym.symbols('a_' +str(i)))
    return letvar,var



def coeffgen(n):
    
    """
    Generating constants for the CY hypersurface equation.
    """
    
    var = [0] * n
    for i in range(n):
        var[i]=sym.expand(sym.symbols('c_' +str(i)))
    return var



def batryev(polyp,dualp,arrcoeff,arrcoords):
    
    """
    Implementing the Batyrev construction which maps the set of vertices corresponding to a 4D reflexive polytope to the hypersurface equation for a CY threefold.
    """
    
    hyp=0
    for i in range(len(dualp)):
        t=1
        coeff=arrcoeff[i]
        for j in range(len(polyp)):
            t1=arrcoords[j]**(np.dot(polyp[j],dualp[i])+1)
            t=t*t1
        hyp=hyp+coeff*t
    return hyp



def CY_Hypersurface(p: LatticePolytope, torus: list[[list]]):
    
    """
    Generating the CY hypersurface, with the ineffective divisors removed.
    """
    
    #Identify fiber coordinates
    verts=p
    poly_verts = LatticePolytope(verts)
    temp=poly_verts.points()
    polyp=np.array(temp)
    #Rearranging
    countemp=0
    for j in range(len(polyp)):
        temp1=polyp.tolist()
        if temp1[j] in torus and temp1[j]!=[0,0,0,0]:
            polyp[countemp]=temp1[j]
            polyp[j]=temp1[countemp]
            countemp=countemp+1
    #Dual Polytope
    dualpoly=poly_verts.polar()
    dualpoly.is_reflexive()
    tempp=dualpoly.points()
    dualp=np.array(tempp)
    #Label divisors
    letvar,arrcoords=divgen(len(polyp),len(torus)-1)
    #coefficients
    arrcoeff=coeffgen(len(dualp))
    #Hypersurface from Batryev
    q=sym.expand(batryev(polyp,dualp,arrcoeff,arrcoords))
    #Remove divisors that don't intersect the CY (Ineffective)
    ineff=list()
    ineffin=list()
    for j in range(len(polyp)):
        temp5=q.subs({arrcoords[j]:0})
        if len(sym.Add.make_args(temp5))==1:
            ineff.append(arrcoords[j])
            if str(arrcoords[j]) in letvar:
                letvar.remove(str(arrcoords[j]))
            ineffin.append(j)  
    #New list of divisors and labels
    kp=polyp.tolist()  
    kcoords=arrcoords.copy()
    temp6=polyp.tolist()
    for j in range(len(ineffin)):
        kp.remove(temp6[ineffin[j]])
        kcoords.remove(arrcoords[ineffin[j]])
    #CY hypersurface with ineffective divisors removed
    weirs=sym.expand(batryev(kp,dualp,arrcoeff,kcoords))
    return kp,kcoords, weirs,letvar,dualp,arrcoeff



def get_2Dbase(p: LatticePolytope, xproj, yproj):
    
    """
    Getting 2D base of the fibration, which is birational to P1 X C (after decompactifying)
    """
    
    newbase1=list()
    for i in range(len(p)):
        newbase1.append([int(np.dot(p[i],xproj)),int(np.dot(p[i],yproj))])
    return newbase1



def rem_nonflat_flav(p: LatticePolytope, coords,  letvar,  base, baseindex1,baseindex2,weirs):
    
    """
    Removing non-flat flavor (non-compact) divisors
    """
    
    flavs=nfflavs(letvar,p,base,baseindex1,baseindex2)
    nbdowns,nbcoords=idnf(letvar,p,base,coords,flavs,baseindex1,baseindex2)
    nf=singmodel(weirs,coords,p,nbdowns)
    nf2=sym.collect(nf,sym.symbols(letvar))
    letvar,model,p,base,coords,nfindex,nbcoords,nbdowns=remnf(letvar,nf2,nbcoords,nbdowns,p,base,coords)
    return letvar,model,p,base,coords,nfindex,nbcoords,nbdowns



def rem_nonflat_comp(p: LatticePolytope, letvar,  base, basedivs,coords, baseindex1,baseindex2,weirs, checkdivs: list[[list]]):
    
    """
    Removing non-flat compact divisors
    """
    
    nonflatdat=list()
    nonflatdatT=list()
    jk=quiverdata(letvar,weirs,checkdivs,p,base,basedivs,coords,baseindex1,baseindex2)
    for i in range(len(checkdivs)):
        if len(jk[i][4])!=0:
            for j in range(len(jk[i][4])):
                nfvert=jk[i][6][jk[i][4][j]][0]
                basecoordj=jk[i][6][jk[i][4][j]][1]
                nfdegree=numint(letvar,jk[i],j,checkdivs[i],baseindex1,baseindex2)[0]
                nonflatdat.append([nfvert,basecoordj,nfdegree])
        else:
            nonflatdatT.append([])
    nftotdats=nonflatdat
    print(nftotdats)
    modelf,p,base,coords,nfindexf,nbcoordsf,nbdownsf=quiverdata(letvar,weirs,checkdivs,p,base,basedivs,coords,baseindex1,baseindex2)[-1]
    basedivs=getbasedivs(base,baseindex1,baseindex2)
    return p,base,nftotdats,coords,baseindex1,baseindex2,basedivs




def singmodel(model,coords,divs,nbdowndivs):
    
    """
    Blowing down divisors in the hypersurface equation (setting them to 1)
    """
    
    tempbk=list()
    for i in range(len(nbdowndivs)):
        tempbk.append(nbdowndivs[i][0])
    for j in range(len(divs)):
        if divs[j] not in tempbk:
            model=model.subs({coords[j]:1})
    return model




def nfflavs(letvar,divs,divsbase,baseindex1,baseindex2):
    
    """
    Identify flavor divisors
    """
    
    kaps=len(letvar)
    flavs=list()
    for i in range(len(divs)):
        if (divsbase[i][baseindex1]!=0 and divsbase[i][baseindex2]==0) or divsbase[i][baseindex2]<0:
            flavs.append(divs[i])
    return flavs




def idnf(letvar,divs,divsbase,coords,targetdivs,baseindex1,baseindex2):
    
    """
    Identify non-flat divisors
    """
    
    kaps=len(letvar)
    nbdowns=list()
    nbcoords=list()
    for i in range(len(divs)):
        if divs[i] in targetdivs:
            nbdowns.append([divs[i],divsbase[i]])
            nbcoords.append(coords[i])
        elif i in range(kaps):
            nbdowns.append([divs[i],divsbase[i]])
            nbcoords.append(coords[i])            
    return nbdowns, nbcoords




def remnf(letvar,model,nbcoords,nbdowns,divs,divsbase,coords):
    
    """
    Remove non-flat divisors
    """
    
    nfindex=list()
    temp=divs.copy()
    tempb=divsbase.copy()
    tempc=coords.copy()
    for i in range(len(nbcoords)):
        if len(sym.Add.make_args(model.subs({nbcoords[i]:0})))==1:
            print(sym.Add.make_args(model.subs({nbcoords[i]:0})))
            inds=tempc.index(sym.symbols(str(nbcoords[i])))
            print(temp[inds])
            temp.remove(temp[inds])
            tempb.remove(tempb[inds])
            if str(tempc[inds]) in letvar:
                letvar.remove(str(tempc[inds]))
            tempc.remove(tempc[inds])
            nfindex.append(i)
            print(len(nfindex))
    return letvar,model,temp,tempb, tempc,nfindex,nbcoords,nbdowns




def getbasedivs(divs,baseindex1,baseindex2):
    
    """
    Getting normalized base divisors
    """
    
    basedivs=list()
    for j in range(len(divs)):
        if (divs[j][baseindex1]!=0 or divs[j][baseindex2]!=0):
            basedivs.append(([int(divs[j][baseindex1]/(math.gcd(divs[j][baseindex1],divs[j][baseindex2]))),int(divs[j][baseindex2]/(math.gcd(divs[j][baseindex1],divs[j][baseindex2])))]))
    return basedivs




def undivs(divs):
    
    """
    Getting the set of unique base divisors
    """
    
    undivs=list()
    for i in range(len(divs)):
        if divs[i] not in undivs:
            undivs.append(divs[i])
    return undivs




def slcomp(divs):
    
    """
    Slope of compact base divisors (slope of a 2D line)
    """
    
    slopes=list()
    compactdivs=list()
    
    for i in range(len(divs)):
        if divs[i][1]>0:
            if divs[i][0]>0:
                slopes.append(math.degrees(math.atan(divs[i][1]/divs[i][0])))
            if divs[i][0]<0:
                slopes.append(math.degrees(math.atan(divs[i][1]/divs[i][0]))+180)
            if divs[i][0]==0:
                slopes.append(90)
    
    #Arranging
    slopes.sort()
    for i in range(len(slopes)):
        for j in range(len(divs)):
            if divs[j][1]>0:
                if divs[j][1]>0 and divs[j][0]>0 and slopes[i]==math.degrees(math.atan(divs[j][1]/divs[j][0])):
                    compactdivs.append(divs[j]) 
                if divs[j][1]>0 and divs[j][0]<0 and slopes[i]==math.degrees(math.atan(divs[j][1]/divs[j][0]))+180:
                    compactdivs.append(divs[j])
                if divs[j][0]==0 and slopes[i]==90:
                    compactdivs.append(divs[j])
    return compactdivs,slopes




def selfints(divs):
    
    """
    Computing intersection data of base divisors (curves)
    """
    
    nums=list()
    nums.append('flavor curve')
    for i in range(len(divs)):
        if len(divs)==1:
            nums.append(0)
        elif i>0 and i<len(divs)-1:
            templ1=np.dot([0,1],divs[i+1])
            templ2=np.dot([0,1],divs[i-1])
            nums.append((templ1+templ2)/(-np.dot([0,1],divs[i])))
        elif i==0:
            templ1=np.dot([0,1],divs[i+1])
            nums.append((templ1)/(-np.dot([0,1],divs[i])))
        elif i==len(divs)-1:
            templ2=np.dot([0,1],divs[i-1])
            nums.append((templ2)/(-np.dot([0,1],divs[i])))
    nums.append('flavor curve')
    return nums



def neighdiv(basedivs,checkdiv):
    
    """
    Identifying neighboring divisors
    """
    
    temp=slcomp(undivs(basedivs))[0]
    inds=temp.index(checkdiv)
    if inds==0:
        return temp[inds],temp[inds+1]
    elif inds==len(temp)-1:
        return temp[inds-1],temp[inds]
    else:
        return temp[inds-1],temp[inds],temp[inds+1]

    

    
def get4Dvert(divs2D,divs4D,divsbase,baseindex1,baseindex2):
    
    """
    Getting 4D vertices from the 2D base vertices
    """
    
    temp=list()
    temp2=list()
    for i in range(len(divsbase)):
         temp.append([divsbase[i][baseindex1],divsbase[i][baseindex2]])
    for i in range(len(temp)):
        if temp[i][0]!=0 or temp[i][1]!=0:
            if [temp[i][0]/math.gcd(int(temp[i][0]),int(temp[i][1])),temp[i][1]/math.gcd(int(temp[i][0]),int(temp[i][1]))] in divs2D:
                temp2.append(divs4D[i])
    return temp2
        


    
def finalmodel(letvar,model,checkdiv,divs,divsbase,bdivs,coords,baseindex1,baseindex2):
    
    """
    Final model after non-flat divisors are removed
    """
    
    klo=neighdiv(bdivs,checkdiv)
    cap=get4Dvert(klo,divs,divsbase,baseindex1,baseindex2)
    #Add torus coordinates
    for i in range(len(letvar)):
        cap.append(divs[i])
    nbdowns,nbcoords=idnf(letvar,divs,divsbase,coords,cap,baseindex1,baseindex2)
    nf=singmodel(model,coords,divs,nbdowns)
    nf2=sym.collect(nf,sym.symbols(letvar))
    return remnf(letvar,nf2,nbcoords,nbdowns,divs,divsbase,coords)



def quiverdata(letvar,model,checkdivs,divs,divsbase,bdivs,coords,baseindex1,baseindex2):
    
    """
    LST represented as a quiver
    """
    
    nfsdat=list()
    if len(checkdivs)!=0:
        for i in range(len(checkdivs)):
            letvar,modelf,kpf,divsbasef,kcoordsf,nfindexf,nbcoordsf,nbdownsf=finalmodel(letvar,model,checkdivs[i],divs,divsbase,bdivs,coords,baseindex1,baseindex2)
            nfsdat.append([modelf,kpf,divsbasef,kcoordsf,nfindexf,nbcoordsf,nbdownsf])
            if len(nfindexf)!=0:    
                for i in range(len(nfindexf)):
                    model=model.subs({nbcoordsf[nfindexf[i]]:1})
            divs=kpf
            divsbase=divsbasef
            bdivs=getbasedivs(divsbasef,baseindex1,baseindex2)
            coords=kcoordsf
    else:
        nfsdat.append([model,divs,divsbase,coords,[],[],[[],[]]])
    return nfsdat



def numint(letvar,quivs,multipli,checkdivs,baseindex1,baseindex2):   
    
    """
    Check number of non-Kahler intersections of a non-flat divisor
    """
    
    if len(quivs[4])!=0:
        quivs[0]=sym.Add.make_args(quivs[0].subs({quivs[5][quivs[4][multipli]]:0}))[0]
        tempty=list()
        for i in range(len(quivs[-1])):
            tempty.append([quivs[-1][i][1][baseindex1],quivs[-1][i][1][baseindex2]])
        for i in range(len(quivs[5])):
            if i not in range(len(letvar)):
                if [tempty[i][0]/math.gcd(int(tempty[i][0]),int(tempty[i][1])),tempty[i][1]/math.gcd(int(tempty[i][0]),int(tempty[i][1]))]==checkdivs:
                    quivs[0]=sym.Add.make_args(quivs[0].subs({quivs[5][i]:1}))[0]
        kaps=quivs[0].args[1]
        kaps2=sym.collect(kaps,sym.symbols(letvar))
        degree=int(sym.degree(kaps2.args[1].as_poly()))
        return degree,'flat fiber with ' +str(degree)+' solutions'
    else:
        return 'No flat fiber on this divisor'




def basedata(divs,divsbase,nonflatdat,baseindex1,baseindex2):
    
    """
    Generating quiver: including gauge and flavor algebras, their rank, Kacs labels and Coulomb Branch dimension
    """
    
    l1=baseindex1
    l2=baseindex2
    gaugevert=list()
    tempgdivs=list()
    for i in range(len(divs)):
        gaugevert.append([])
        if (divsbase[i][l1]!=0 or divsbase[i][l2]!=0):
            compt=([divsbase[i][l1]/(math.gcd(divsbase[i][l1],divsbase[i][l2])),divsbase[i][l2]/(math.gcd(divsbase[i][l1],divsbase[i][l2]))]) 
            if compt not in tempgdivs:
                tempgdivs.append(compt)
                gaugevert[i].append([divs[i],divsbase[i]])
                for j in range(len(divs)):
                    if (divsbase[j][l1]!=0 or divsbase[j][l2]!=0):
                        compt2=([divsbase[j][l1]/(math.gcd(divsbase[j][l1],divsbase[j][l2])),divsbase[j][l2]/(math.gcd(divsbase[j][l1],divsbase[j][l2]))])
                        if compt2==compt and i!=j:  
                            gaugevert[i].append([divs[j],divsbase[j]])
    fgaugevert=list()
    for i in range(len(gaugevert)):
        if len(gaugevert[i])>0 and gaugevert[i][0][1][l2]>=0:
            fgaugevert.append(gaugevert[i])
    arrtempg=list()
    arrtempg.append(0)
    
    basedivstemp=list()
    for j in range(len(divs)):
        if (divsbase[j][l1]!=0 or divsbase[j][l2]!=0):
            basedivstemp.append(([divsbase[j][l1]/(math.gcd(divsbase[j][l1],divsbase[j][l2])),divsbase[j][l2]/(math.gcd(divsbase[j][l1],divsbase[j][l2]))]))

    for i in range(len(slcomp(undivs(basedivstemp))[1])):
        arrtempg.append(slcomp(undivs(basedivstemp))[1][i])
    intersecs=selfints(slcomp(undivs(basedivstemp))[0])    
    arrtempg.append(180)
    arrvert=list()
    for j in range(len(arrtempg)):
        for i in range(len(fgaugevert)):
            if fgaugevert[i][0][1][l1]<0:
                sema=math.degrees(math.atan(fgaugevert[i][0][1][l2]/fgaugevert[i][0][1][l1]))+180
            if fgaugevert[i][0][1][l1]>0:
                sema=math.degrees(math.atan(fgaugevert[i][0][1][l2]/fgaugevert[i][0][1][l1]))
            if fgaugevert[i][0][1][l1]==0:
                sema=90
            if sema==arrtempg[j]:
                arrvert.append(fgaugevert[i])
    kaclabels=list()
    for i in range(len(arrvert)):
        kaclabels.append([])
        for j in range(len(arrvert[i])):
            kaclabels[i].append(math.gcd(arrvert[i][j][1][l1::l2-l1][0],arrvert[i][j][1][l1::l2-l1][1]))
        for k in range(len(nonflatdat)):
            if len(nonflatdat[k])!=0:
                if nonflatdat[k][0] in arrvert[i]:
                    for m in range(nonflatdat[k][1]-1):
                        kaclabels[i].append(1)
                
    CBdim=len(intersecs)-2
    for i in range(len(kaclabels)):
        if i==0 or i==len(kaclabels)-1:
            kaclabels[i].append(['rank ' +str(len(kaclabels[i])-1)+ ' flavor group'])
        else:
            kaclabels[i].append(['rank ' +str(len(kaclabels[i])-1)+ ' gauge group on '+ str(int(abs(intersecs[i])))+' curve'])
            if len(kaclabels[i])>2:
                CBdim=CBdim+len(kaclabels[i])-2
    return kaclabels,CBdim-1



def lscharge(divs):    
    
    """
    Little String Charges computed from base data
    """
    
    selfs=selfints(slcomp(undivs(divs))[0])
    selfs.pop(0)
    selfs.pop(-1)
    intmat= np.zeros((len(selfs),len(selfs)))
    for i in range(len(selfs)):
        for j in range(len(selfs)):
            if j==i+1 or j==i-1:
                intmat[i,j]=1
            elif i==j:
                intmat[j,j]=selfs[i]
            else:    
                intmat[i,j]=0
    ns = null_space(intmat)
    LScharges=abs(ns/-min(abs(ns)))
    return intmat,LScharges



def Two_Group(p: LatticePolytope,base,basedivs,nftotdats,baseindex1,baseindex2): 
    
    """
    Computing 2-group structure constants (k_{P}, k_{R})
    """
    
    intmat,LScharges=lscharge(basedivs)
    kpoincare=0
    kR=0
    selfs=selfints(slcomp(undivs(basedivs))[0])
    selfs.pop(0)
    selfs.pop(-1)
    cf=basedata(p,base,nftotdats,baseindex1,baseindex2)[0]
    cf.pop(0)
    cf.pop(-1)
    for i in range(len(LScharges)):
        kpoincare=kpoincare-(LScharges[i]*(-selfs[i]-2))
    for i in range(len(LScharges)):
        tempk=cf[i]
        tempk.pop(-1)
        kR=kR+((np.sum(tempk))*LScharges[i])
    print([round(kpoincare[0]),round(kR[0])])

