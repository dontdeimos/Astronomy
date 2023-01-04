#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 14:00:27 2018

@author: sarahdeveny
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from scipy import interpolate
import numpy.polynomial.polynomial as poly
from scipy.interpolate import interp1d
from matplotlib import rc
import unicodedata
from matplotlib import pyplot




name = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3',usecols=[3], skiprows=4,dtype='str')
F300X, F390W = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3',usecols=[4,10], skiprows=4, unpack=True)
q300, q390 = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3',usecols=[6,12], skiprows=4, unpack=True)
o300, o390 = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3',usecols=[7,13], skiprows=4, unpack=True)
g300, g390 = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3',usecols=[9,15], skiprows=4, unpack=True)

#Cut the 2 duplicates
u, indexdupremoval = np.unique(name, return_index=True)

name = name[indexdupremoval]
F300X = F300X[indexdupremoval]
F390W = F390W[indexdupremoval]
q300 = q300[indexdupremoval]
q390 = q390[indexdupremoval]
o300 = o300[indexdupremoval]
o390 = o390[indexdupremoval]
g300 = g300[indexdupremoval]
g390 = g390[indexdupremoval]

#Cut the g's right away
cutg = []
for j in range(0,len(g300)):
    if g300[j] >= 2 and g390[j] >= 2:
        cutg.append(j)
        
name = name[cutg]
F300X = F300X[cutg]
F390W = F390W[cutg]
q300 = q300[cutg]
q390 = q390[cutg]
o300 = o300[cutg]
o390 = o390[cutg]
g300 = g300[cutg]
g390 = g390[cutg]

#Cut the flux that's less than zero
cutflux = []
for j in range(0,len(g300)):
    if F300X[j] > 0 and F390W[j] > 0:
        cutflux.append(j)

name = name[cutflux]
F300X = F300X[cutflux]
F390W = F390W[cutflux]
q300 = q300[cutflux]
q390 = q390[cutflux]
o300 = o300[cutflux]
o390 = o390[cutflux]
g300 = g300[cutflux]
g390 = g390[cutflux]

#Adjust parameters
o300 = o300/F300X
o390 = o390/F390W

#To get instrumental magnitudes
#F300X = -2.5*np.log10(F300X)
#F390W = -2.5*np.log10(F390W)

##To get apparent magnitudes
m300inst =  -2.5*np.log10(F300X)
m390inst =  -2.5*np.log10(F390W)
#                
#F300X = np.add(m300inst,32)    
#F390W = np.add(m390inst,32) 

cal300 = 1.643
cal390 = 0.044

####STILL NEED TO ADD 32 TO MAKE IT COMPLETELY CALIBRATED *DONE AT END*                    
                     
F300X = m300inst# - 1.643
F390W = m390inst# - 0.044

#F300X = list(map(lambda x:x +32, F300X))
#F390W = list(map(lambda x:x +32, F390W))

#F300X = np.array(F300X)
#F390W = np.array(F390W)
#F300X = [x + 32 for x in F300X]
#F390W = [x + 32 for x in F390W]

g300 = np.divide(g300,24)
g390 = np.divide(g390,24)

#print len(name)  #134407 stars

      
#-------------------Test arrays--------------------------
#name = np.array([003,002,001,004])
#F300X = np.array([10.,-6.,-1.,2.])
#q300 = np.array([.2,.9,.1,.8])

#print name
#print F300X
#print q300

######################only want white dwarf sequence###########################
#bv = F300X - F390W
#bvnew = np.where(bv < 2.)
#F300X = F300X[bvnew]
#F390W = F390W[bvnew]
#name = name[bvnew]
#q300 = q300[bvnew]
#q390 = q390[bvnew]
#o300 = o300[bvnew]
#o390 = o390[bvnew]


#########################Defining functions#########################3
#def sort(param,mag,name):
#    sorting = param.argsort()
#    magsorted = mag[sorting[::]]
#    paramsorted = param[sorting[::]]
#    namesorted = name[sorting[::]]
#    return magsorted, paramsorted, namesorted
#
#def binning(xmin,xmax,Nbin,magsorted,paramsorted,namesorted):
#    bins = np.linspace(5,-15,200)
#    whichbin = np.digitize(F300Xsorted,bins)
#
#    mag = [magsorted[whichbin == i] for i in range(1,len(bins)+1)]
#    param = [paramsorted[whichbin == i] for i in range(1,len(bins)+1)]
#    Name = [namesorted[whichbin == i] for i in range(1,len(bins)+1)]
#    
#    return 

################################Sorting by q300################################
#binsort = F300X.argsort()
#F300X = F300X[binsort[::]]
#q300 = q300[binsort[::]]
#name = name[binsort[::]]
#F390W = F390W[binsort[::]]

#------------The argsort gives me an array of the indices that sort the q 
# parameter, applying this sorting array to each parameter will organize
# every array the same way ------------------
sorting = q300.argsort()
F300Xsorted = F300X[sorting[::]]
q300sorted = q300[sorting[::]]
namesorted300 = name[sorting[::]]
#o300sorted = o300[sorting[::]]
#F390Wsorted = F390W[sorting[::]]
#print sorting
#print whichbin
#print namesorted300
#print F300Xsorted
#print q300sorted

#---------------I don't think I need to sort, the linspace breaks up the bins
# and then the digitize organizes the magnitude array into each bin -----------
bins = np.linspace(-15.,5,200)
whichbin = np.digitize(F300Xsorted,bins)
#print whichbin
#print name
#print F300X
#print q300

#-------------This will take the sorted paramters and apply to each bin -------

F300q = [F300Xsorted[whichbin == i] for i in range(1,len(bins)+1)]
Q300q = [q300sorted[whichbin == i] for i in range(1,len(bins)+1)]
Name300q = [namesorted300[whichbin == i] for i in range(1,len(bins)+1)]

#O300q = [o300sorted[whichbin == i] for i in range(1,len(bins)+1)]
#F390 = [F390Wsorted[whichbin == i] for i in range(1,len(bins)+1)]
#print Name300
#print F300
#print F300[5][0]

#-----------This determines the length of each bin------------------------
lengths=[]
for i in range(0,len(F300q)):
    lengths.append(len(F300q[i]))
    
#-----------This defines a new array that matches the shape of the mag and 
# name arrays -------------------------------------    
percentq300 = []
for i in range(len(F300q)):
    percentq300.append(np.zeros(lengths[i]))

#-------------This finally creates the new array that defines each index with 
# a percentage -----------------------------------------
for i in range(0,len(Name300q)):
    for j in range(0,lengths[i]):
        percentq300[i][j]=float(j+1)/float(lengths[i])
        
################################Concatenate####################################
#------------This unravels the binned arrays into a single array---------------
Name300q = np.concatenate(Name300q).ravel()
F300q = np.concatenate(F300q).ravel()
Q300q = np.concatenate(Q300q).ravel()
Percentq300 = np.concatenate(percentq300).ravel()
#idk = np.where(F300q)
#Percentq300 = Percentq300[idk]
#O300 = np.concatenate(O300).ravel()
#F390 = np.concatenate(F390).ravel()
#print len(Name300q)
#print F300
#print Q300
#print F300[5]
#print len(Percentq300)

####Sort back#####
#-----------This sorts the array by name 1-N --------------------
sort = Name300q.argsort()
Name300q = Name300q[sort]
Name300q = np.reshape(Name300q,len(Name300q))
F300q = F300q[sort]
F300q = np.reshape(F300q,len(F300q))
Q300q = Q300q[sort]
Q300q = np.reshape(Q300q,len(Q300q))
Percentq300 = Percentq300[sort] 
Percentq300 = np.reshape(Percentq300,len(Percentq300))
#O300q = O300q[sort]
#F390 = F390[sort]
#print Name300
#print F300

#sortnew = name.argsort()
#F390 = F390W[sortnew]
#best300 = [F300[i][int(len(Q300[i]) * .02):] for i in range(0,len(Q300))]

##############################Sorting by q390##################################

sort390 = q390.argsort()
F390Wsorted = F390W[sort390[::]]
q390sorted = q390[sort390[::]]
namesorted390 = name[sort390[::]]
#o390sorted = o300[sort390[::]]

bins = np.linspace(5,-15,200)
whichbin = np.digitize(F390Wsorted,bins)

F390q = [F390Wsorted[whichbin == i] for i in range(1,len(bins)+1)]
Q390q = [q390sorted[whichbin == i] for i in range(1,len(bins)+1)]
Name390q = [namesorted390[whichbin == i] for i in range(1,len(bins)+1)]
#O390 = [o390sorted[whichbin == i] for i in range(1,len(bins)+1)]
#print F300[5][0]

lengths=[]
for i in range(0,len(F390q)):
    lengths.append(len(F390q[i]))
    
percent390q = []
for i in range(len(F390q)):
    percent390q.append(np.zeros(lengths[i]))

for i in range(0,len(Name390q)):
    for j in range(0,lengths[i]):
        percent390q[i][j]=float(j+1)/float(lengths[i])
        
####Concatenate#####
Name390q = np.concatenate(Name390q).ravel()
F390q = np.concatenate(F390q).ravel()
Q390q = np.concatenate(Q390q).ravel()
#O390 = np.concatenate(O390).ravel()
Percent390q = np.concatenate(percent390q).ravel()
#idk = np.where(F390q)
#Percent390q = Percent390q[idk]
#print F300[5]

####Sort back#####
sort = Name390q.argsort()
Name390q = Name390q[sort]
Name390q = np.reshape(Name390q,len(Name390q))
F390q = F390q[sort]
F390q = np.reshape(F390q,len(F390q))
Q390q = Q390q[sort]
Q390q = np.reshape(Q390q,len(Q390q))
#O390 = O390[sort]
Percent390q = Percent390q[sort] 
Percent390q = np.reshape(Percent390q,len(Percent390q))

############################Sorting by o300####################################

sorto300 = o300.argsort()
F300sorted = F300X[sorto300[::]]
#q300sorted = q390[sorto300[::]]
o300sorted = o300[sorto300[::]]
namesorted300 = name[sorto300[::]]

bins = np.linspace(5,-15,200)
whichbin = np.digitize(F300sorted,bins)

F300o = [F300sorted[whichbin == i] for i in range(1,len(bins)+1)]
#Q300o = [q300sorted[whichbin == i] for i in range(1,len(bins)+1)]
O300o = [o300sorted[whichbin == i] for i in range(1,len(bins)+1)]
Name300o = [namesorted300[whichbin == i] for i in range(1,len(bins)+1)]
#print F300[5][0]

lengths=[]
for i in range(0,len(F300o)):
    lengths.append(len(F300o[i]))
    
percent300o = []
for i in range(len(F300o)):
    percent300o.append(np.zeros(lengths[i]))

for i in range(0,len(Name300o)):
    for j in range(0,lengths[i]):
        percent300o[i][j]=float(j+1)/float(lengths[i])
        
####Concatenate#####
Name300o = np.concatenate(Name300o).ravel()
F300o = np.concatenate(F300o).ravel()
#Q300o = np.concatenate(Q300).ravel()
O300o = np.concatenate(O300o).ravel()
Percent300o = np.concatenate(percent300o).ravel()
#idk = np.where(F300o)
#Percent300o = Percent300o[idk]
#print F300[5]

####Sort back#####
sort = Name300o.argsort()
Name300o = Name300o[sort]
Name300o = np.reshape(Name300o,len(Name300o))
F300o = F300o[sort]
F300o = np.reshape(F300o,len(F300o))
#Q390o = Q300[sort]
O300o = O300o[sort]
O300o = np.reshape(O300o,len(O300o))
Percent300o = Percent300o[sort] 
Percent300o = np.reshape(Percent300o,len(Percent300o))

######################## Sorting o390 #########################################

sort390 = o390.argsort()
F390sorted = F390W[sort390[::]]
#q390sorted = Q390[sort390[::]]
o390sorted = o390[sort390[::]]
namesorted390 = name[sort390[::]]

bins = np.linspace(5,-15,200)
whichbin = np.digitize(F390sorted,bins)

F390o = [F390sorted[whichbin == i] for i in range(1,len(bins)+1)]
O390o = [o390sorted[whichbin == i] for i in range(1,len(bins)+1)]
Name390o = [namesorted390[whichbin == i] for i in range(1,len(bins)+1)]
#Q390 = [q390sorted[whichbin == i] for i in range(1,len(bins)+1)]
#print F300[5][0]

lengths=[]
for i in range(0,len(F390o)):
    lengths.append(len(F390o[i]))
    
percent390o = []
for i in range(len(F390o)):
    percent390o.append(np.zeros(lengths[i]))

for i in range(0,len(Name390o)):
    for j in range(0,lengths[i]):
        percent390o[i][j]=float(j+1)/float(lengths[i])
        
####Concatenate#####
Name390o = np.concatenate(Name390o).ravel()
F390o = np.concatenate(F390o).ravel()
#Q390o = np.concatenate(Q390).ravel()
O390o = np.concatenate(O390o).ravel()
Percent390o = np.concatenate(percent390o).ravel()
#idk = np.where(F390o)
#Percent390o = Percent390o[idk]
#print F300[5]

####Sort back#####
sort = Name390o.argsort()
Name390o = Name390o[sort]
Name390o = np.reshape(Name390o,len(Name390o))
F390o = F390o[sort]
F390o = np.reshape(F390o,len(F390o))
#Q390o = Q390o[sort]
O390o = O390o[sort]
O390o = np.reshape(O390o,len(O390o))
Percent390o = Percent390o[sort] 
Percent390o = np.reshape(Percent390o,len(Percent390o))

#print len(O390o)
#print len(Percent390o)

######################## Sorting g300 #########################################

sortg300 = g300.argsort()
F300gsorted = F300X[sortg300[::]]
#q390sorted = Q390[sort390[::]]
g300sorted = g300[sortg300[::]]
namesorted300g = name[sortg300[::]]

bins = np.linspace(5,-15,200)
whichbin = np.digitize(F300gsorted,bins)

F300g = [F300gsorted[whichbin == i] for i in range(1,len(bins)+1)]
G300g = [g300sorted[whichbin == i] for i in range(1,len(bins)+1)]
Name300g = [namesorted300g[whichbin == i] for i in range(1,len(bins)+1)]
#Q390 = [q390sorted[whichbin == i] for i in range(1,len(bins)+1)]
#print F300[5][0]

lengths=[]
for i in range(0,len(F300g)):
    lengths.append(len(F300g[i]))
    
percent300g = []
for i in range(len(F300g)):
    percent300g.append(np.zeros(lengths[i]))

for i in range(0,len(Name300g)):
    for j in range(0,lengths[i]):
        percent300g[i][j]=float(j+1)/float(lengths[i])
        
####Concatenate#####
Name300g = np.concatenate(Name300g).ravel()
F300g = np.concatenate(F300g).ravel()
#Q390o = np.concatenate(Q390).ravel()
G300g = np.concatenate(G300g).ravel()
Percent300g = np.concatenate(percent300g).ravel()
#idk = np.where(F300g)
#Percent300g = Percent300g[idk]
#print F300[5]

####Sort back#####
sort = Name300g.argsort()
Name300g = Name300g[sort]
Name300g = np.reshape(Name300g,len(Name300g))
F300g = F300g[sort]
F300g = np.reshape(F300g,len(F300g))
#Q390o = Q390o[sort]
G300g = G300g[sort]
G300g = np.multiply(G300g,24)
G300g = np.reshape(G300g,len(G300g))
Percent300g = Percent300g[sort] 
Percent300g = np.reshape(Percent300g,len(Percent300g))

#print len(O390o)
#print len(Percent390o)

######################## Sorting g390 #########################################

sortg390 = g390.argsort()
F390gsorted = F390W[sortg390[::]]
#q390sorted = Q390[sort390[::]]
g390sorted = g390[sortg390[::]]
namesorted390g = name[sortg390[::]]

bins = np.linspace(5,-15,200)
whichbin = np.digitize(F390gsorted,bins)

F390g = [F390gsorted[whichbin == i] for i in range(1,len(bins)+1)]
G390g = [g390sorted[whichbin == i] for i in range(1,len(bins)+1)]
Name390g = [namesorted390g[whichbin == i] for i in range(1,len(bins)+1)]
#Q390 = [q390sorted[whichbin == i] for i in range(1,len(bins)+1)]
#print F300[5][0]

lengths=[]
for i in range(0,len(F390g)):
    lengths.append(len(F390g[i]))
    
percent390g = []
for i in range(len(F390g)):
    percent390g.append(np.zeros(lengths[i]))

for i in range(0,len(Name390g)):
    for j in range(0,lengths[i]):
        percent390g[i][j]=float(j+1)/float(lengths[i])
        
####Concatenate#####
Name390g = np.concatenate(Name390g).ravel()
F390g = np.concatenate(F390g).ravel()
#Q390o = np.concatenate(Q390).ravel()
G390g = np.concatenate(G390g).ravel()
Percent390g = np.concatenate(percent390g).ravel()
#idk = np.where(F390g)
#Percent390g = Percent390g[idk]
#print F300[5]

####Sort back#####
sort = Name390g.argsort()
Name390g = Name390g[sort]
Name390g = np.reshape(Name390g,len(Name390g))
F390g = F390g[sort]
F390g = np.reshape(F390g,len(F390g))
#Q390o = Q390o[sort]
G390g = G390g[sort]
G390g = np.multiply(G390g,24)
G390g = np.reshape(G390g,len(G390g))
Percent390g = Percent390g[sort] 
Percent390g = np.reshape(Percent390g,len(Percent390g))

#print len(O390o)
#print len(Percent390o)


########################Find the matching names################################

#def returnMatches(a,b):
#       return list(set(a) & set(b))    
#
#matches = returnMatches(Name300,Name390)

######################Create final mag lists with matching indices#############

#def indexmatch(a,b):
#    return a.searchsorted(b)

#returns an array of the matching values
#def indexmatch(a,b):
#    return set(a).intersection(b)

#def indexmatch(a,b):
#    return np.in1d(a,b)
    
def indexmatch(a,b):
  b_set = set(b)
  return [i for i, item in enumerate(a) if item in b_set]

#def indexmatch(a,b):
#    return np.where(a==b)

#def indexmatch(a,b):
#    return [x for x in a if x in b]
#-------------------Matching the q's---------------------------------------
####PUT THE LONGER ARRAY FIRST TO MATCH WITH THE SHORTER ARRAY#########
matchindexq300 = indexmatch(Name300q,Name390q)

Name300q = Name300q[matchindexq300]
F300q = F300q[matchindexq300]
Q300q = Q300q[matchindexq300]
Percent300q = Percentq300[matchindexq300]

#######THEN SWITCH TO FIX THE SHORTER ONE##########
matchindexq390 = indexmatch(Name390q,Name300q)

Name390q = Name390q[matchindexq390]
F390q = F390q[matchindexq390]
Q390q = Q390q[matchindexq390]
Percentq390 = Percent390q[matchindexq390]

#----------------------Matching the o's--------------------------------------
matchindexo390 = indexmatch(Name390o,Name300o)

Name390o = Name390o[matchindexo390]
F390o = F390o[matchindexo390]
O390o = O390o[matchindexo390]
Percent390o = Percent390o[matchindexo390]

matchindexo300 = indexmatch(Name300o,Name390o)

Name300o = Name300o[matchindexo300]
F300o = F300o[matchindexo300]
O300o = O300o[matchindexo300]
Percent300o = Percent300o[matchindexo300]

#----------------------Matching the g's--------------------------------------
matchindexg300 = indexmatch(Name300g,Name390g)

Name300g = Name300g[matchindexg300]
F300g = F300g[matchindexg300]
G300g = G300g[matchindexg300]
Percent300g = Percent300g[matchindexg300]

matchindexg390 = indexmatch(Name390g,Name300g)

Name390g = Name390g[matchindexg390]
F390g = F390g[matchindexg390]
G390g = G390g[matchindexg390]
Percent390g = Percent390g[matchindexg390]



#print len(Name390g)
#print len(Name300g)

#------------Match the o's, q's and g's--------------------
matchqo = indexmatch(Name300q,Name300o)

Name300Q = Name300q[matchqo] 
Name390Q = Name390q[matchqo]
F300Q = F300q[matchqo]
F390Q = F390q[matchqo]
Q300Q = Q300q[matchqo]
Q390Q = Q390q[matchqo]
Percent300Q = Percent300q[matchqo]
Percent390Q = Percent390q[matchqo]

matchoq = indexmatch(Name300o,Name300Q)

Name300QO = Name300o[matchoq] 
Name390QO = Name390o[matchoq]
F300QO = F300o[matchoq]
F390QO = F390o[matchoq]
O300QO = O300o[matchoq]
O390QO = O390o[matchoq]
Percent300QO = Percent300o[matchoq]
Percent390QO = Percent390o[matchoq]

matchqg = indexmatch(Name300Q,Name300g)

Name300Qg = Name300Q[matchqg] 
Name390Qg = Name390Q[matchqg]
F300Qg = F300Q[matchqg]
F390Qg = F390Q[matchqg]
Q300Qg = Q300Q[matchqg]
Q390Qg = Q390Q[matchqg]
Percent300Qg = Percent300Q[matchqg]
Percent390Qg = Percent390Q[matchqg]
Name300QOg = Name300QO[matchoq] 
Name390QOg = Name390QO[matchoq]
F300QOg = F300QO[matchoq]
F390QOg = F390QO[matchoq]
O300QOg = O300QO[matchoq]
O390QOg = O390QO[matchoq]
Percent300QOg = Percent300QO[matchoq]
Percent390QOg = Percent390QO[matchoq]

matchgq = indexmatch(Name300g,Name300Qg)

Name300QOG = Name300g[matchgq] 
Name390QOG = Name390g[matchgq]
F300QOG = F300g[matchgq]
F390QOG = F390g[matchgq]
G300QOG = G300g[matchgq]
G390QOG = G390g[matchgq]
Percent300QOG = Percent300g[matchgq]
Percent390QOG = Percent390g[matchgq]

#Needs to be sorted back to alphabetical just in case##
sortq = Name300Q.argsort()
sorto = Name300QO.argsort()
sortg = Name300QOG.argsort()

Name300Qg = Name300Q[sortq] 
Name390Qg = Name390Q[sortq]
F300Qg = F300Q[sortq]
F390Qg = F390Q[sortq]
Q300Qg = Q300Q[sortq]
Q390Qg = Q390Q[sortq]
Percent300Qg = Percent300Q[sortq]
Percent390Qg = Percent390Q[sortq]
Name300QOg = Name300QO[sorto] 
Name390QOg = Name390QO[sorto]
F300QOg = F300QO[sorto]
F390QOg = F390QO[sorto]
O300QOg = O300QO[sorto]
O390QOg = O390QO[sorto]
Percent300QOg = Percent300QO[sorto]
Percent390QOg = Percent390QO[sorto]
Name300QOG = Name300g[sortg] 
Name390QOG = Name390g[sortg]
F300QOG = F300g[sortg]
F390QOG = F390g[sortg]
G300QOG = G300g[sortg]
G390QOG = G390g[sortg]
Percent300QOG = Percent300g[sortg]
Percent390QOG = Percent390g[sortg]

######Final Values before the cuts#######
Name = Name300QOG
F300 = F300QOG
F390 = F390QOG
Q300 = Q300Qg
Q390 = Q390Qg
O300 = O300QOg
O390 = O390QOg
G300 = G300QOG
G390 = G390QOG
Percent300q = Percent300Qg
Percent390q = Percent390Qg
Percent300o = Percent300QOg
Percent390o = Percent390QOg
Percent300g = Percent300QOG
Percent390g = Percent390QOG

###################BEST##################
#best300 = []
#for i in range(0,len(Q300)):
#    if Q300[i] >= .5:
#        best300.append(np.where(Q300[i]))#best300.append(Q300[i])

###Tell me where (gives an index) in the array this condition is met###
            #which index of idxno has a Q that is > .99 apply those indices to idxno 
            #and that should give me an array of indices from the original 
            #that shouldn't be thrown out    
            
##idxq300 = np.where(Percent300q >= 0.1) 
#Tell me (give index) which values in percent are >= 'whatever'
idxq300 = np.array([i for i, x in enumerate(Percent300q) if x >= 0.2])
#Tell me (give index) which values didn't make it into 'whatever'
idxno300 = np.array([i for i, x in enumerate(Percent300q) if x < 0.2])
##Apply the indices that didn't make it to Q
Q300test = Q300[idxno300]
##Which values (give index) that didn't make it are greater than 'whatever'
indices300 = np.array([i for i, x in enumerate(Q300test) if x >= 0.990])
##Apply the indices from shorter array to get values for indices on the longer array
###########To add##########
idxadd300 = idxno300[indices300]
idxq300 = np.append(idxq300,idxadd300)
#Add the values that are > 'whatever' back to the original indices
#We don't want to lose well measured stars
######To delete########
#idxq300 = np.delete(idxq300,indices300)
#----------------
idxq300.sort()
Nameq = Name[idxq300]
F300q = F300[idxq300]
F390q = F390[idxq300]
Q300q = Q300[idxq300]
Q390q = Q390[idxq300]
O300q = O300[idxq300]
O390q = O390[idxq300]
G300q = G300[idxq300]
G390q = G390[idxq300]
##idxq390 = np.where(Percent390q >= 0.1)
idxq390 = np.array([i for i, x in enumerate(Percent390q) if x >= 0.2])
idxno390 = np.array([i for i, x in enumerate(Percent390q) if x < 0.2])
Q390test = Q390[idxno390]
indices390 = np.array([i for i, x in enumerate(Q390test) if x >= 0.990])
########To add#########
idxadd390 = idxno390[indices390]
idxq390 = np.append(idxq390,idxadd390)
########To delete#######
#idxq390 = np.delete(idxq390,indices390)
##-----------------
idxq390.sort()
Nameqq = Name[idxq390]
F300qq = F300[idxq390]
F390qq = F390[idxq390]
Q300qq = Q300[idxq390]
Q390qq = Q390[idxq390]
O300qq = O300[idxq390]
G300qq = G300[idxq390]
G390qq = G390[idxq390]
#####To filter the O's through the new total, we need to change where we get 
#####the initial array
idxo300 = np.where(Percent300o <= 0.8)   #This needs to be "<= .8" for 20% cut
Nameo = Name[idxo300]
F300o = F300[idxo300]
F390o = F390[idxo300]
Q300o = Q300[idxo300]
Q390o = Q390[idxo300]
O300o = O300[idxo300]
O390o = O390[idxo300]
G300o = G300[idxo300]
G390o = G390[idxo300]
idxo390 = np.where(Percent390o <= 0.8)
Nameoo = Name[idxo390]
F300oo = F300[idxo390]
F390oo = F390[idxo390]
Q300oo = Q300[idxo390]
Q390oo = Q390[idxo390]
O300oo = O300[idxo390]
O390oo = O390[idxo390]
G300oo = G300[idxo390]
G390oo = G390[idxo390]
idxg300 = np.where(Percent300g >= 0.0)
Nameg = Name[idxg300]
F300g = F300[idxg300]
F390g = F390[idxg300]
Q300g = Q300[idxg300]
Q390g = Q390[idxg300]
O300g = O300[idxg300]
O390g = O390[idxg300]
G300g = G300[idxg300]
G390g = G390[idxg300]
idxg390 = np.where(Percent390g >= 0.0)
Namegg = Name[idxg390]
F300gg = F300[idxg390]
F390gg = F390[idxg390]
Q300gg = Q300[idxg390]
Q390gg = Q390[idxg390]
O300gg = O300[idxg390]
O390gg = O390[idxg390]
G300gg = G300[idxg390]
G390gg = G390[idxg390]

#############re-match everything#############
#rematch the q's#
matchindexq = indexmatch(Nameq,Nameqq)

NameQ = Nameq[matchindexq]
F300Q = F300q[matchindexq]
Q300Q = Q300q[matchindexq]
Percent300Q = Percent300q[matchindexq]

matchindexqq = indexmatch(Nameqq,Nameq)

NameQQ = Nameqq[matchindexqq]
F390QQ = F390qq[matchindexqq]
Q390QQ = Q390qq[matchindexqq]
Percent390QQ = Percent390q[matchindexqq]

#rematch the o's#
matchindexo = indexmatch(Nameo,Nameoo)

NameO = Nameo[matchindexo]
F300O = F300o[matchindexo]
O300O = O300o[matchindexo]
Percent300O = Percent300o[matchindexo]

matchindexoo = indexmatch(Nameoo,Nameo)

NameOO = Nameoo[matchindexoo]
F390OO = F390oo[matchindexoo]
O390OO = O390oo[matchindexoo]
Percent390OO = Percent390o[matchindexoo]

#rematch the g's#

matchindexg = indexmatch(Nameg,Namegg)

NameG = Nameg[matchindexg]
F300G = F300g[matchindexg]
G300G = G300g[matchindexg]
Percent300G = Percent300g[matchindexg]

matchindexgg = indexmatch(Namegg,Nameg)

NameGG = Namegg[matchindexgg]
F390GG = F390gg[matchindexgg]
G390GG = G390gg[matchindexgg]
Percent390GG = Percent390g[matchindexgg]

#print len(NameG)
#print len(NameGG)

#------------Re-match the o's, q's and g's--------------------
######Match O and Q#########
matchqo = indexmatch(NameQ,NameO)

Name300Qo = NameQ[matchqo] 
Name390Qo = NameQQ[matchqo]
F300Qo = F300Q[matchqo]
F390Qo = F390QQ[matchqo]
Q300Qo = Q300Q[matchqo]
Q390Qo = Q390QQ[matchqo]
Percent300Qo = Percent300Q[matchqo]
Percent390Qo = Percent390QQ[matchqo]

matchoq = indexmatch(NameO,NameQ)

Name300QO = NameO[matchoq] 
Name390QO = NameOO[matchoq]
F300QO = F300O[matchoq]
F390QO = F390OO[matchoq]
O300QO = O300O[matchoq]
O390QO = O390OO[matchoq]
Percent300QO = Percent300O[matchoq]
Percent390QO = Percent390OO[matchoq]

####Match g with Q (and O)#######
matchqg1 = indexmatch(Name300Qo,NameG)

Name300Qg = Name300Qo[matchqg1] 
Name390Qg = Name390Qo[matchqg1]
F300Qg = F300Qo[matchqg1]
F390Qg = F390Qo[matchqg1]
Q300Qg = Q300Qo[matchqg1]
Q390Qg = Q390Qo[matchqg1]
Percent300Qg = Percent300Qo[matchqg1]
Percent390Qg = Percent390Qo[matchqg1]
Name300QOg = Name300QO[matchqg1]
Name390QOg = Name390QO[matchqg1]
F300QOg = F300QO[matchqg1]
F390QOg = F390QO[matchqg1]
O300QOg = O300QO[matchqg1]
O390QOg = O390QO[matchqg1]
Percent300QOg = Percent300QO[matchqg1]
Percent390QOg = Percent390QO[matchqg1]

matchgq2 = indexmatch(NameG,Name300Qg)

Name300Gq = NameG[matchgq2] 
Name390Gq = NameGG[matchgq2]
F300Gq = F300G[matchgq2]
F390Gq = F390GG[matchgq2]
G300Gq = G300G[matchgq2]
G390Gq = G390GG[matchgq2]
Percent300Gq = Percent300G[matchgq2]
Percent390Gq = Percent390GG[matchgq2]

#print len(G300Gq)
#print len(G390Gq)

#############LAST ALPHABETICAL SORT############
sortq = Name300Qg.argsort()
sorto = Name300QOg.argsort()
sortg = Name300Gq.argsort()

Name300Qg = Name300Qg[sortq] 
Name390Qg = Name390Qg[sortq]
F300Qg = F300Qg[sortq]
F390Qg = F390Qg[sortq]
Q300Qg = Q300Qg[sortq]
Q390Qg = Q390Qg[sortq]
Percent300Qg = Percent300Qg[sortq]
Percent390Qg = Percent390Qg[sortq]
Name300QOg = Name300QOg[sorto] 
Name390QOg = Name390QOg[sorto]
F300QOg = F300QOg[sorto]
F390QOg = F390QOg[sorto]
O300QOg = O300QOg[sorto]
O390QOg = O390QOg[sorto]
Percent300QOg = Percent300QOg[sorto]
Percent390QOg = Percent390QOg[sorto]
Name300QOG = Name300Gq[sortg] 
Name390QOG = Name390Gq[sortg]
F300QOG = F300Gq[sortg]
F390QOG = F390Gq[sortg]
G300QOG = G300Gq[sortg]
G390QOG = G390Gq[sortg]
Percent300QOG = Percent300Gq[sortg]
Percent390QOG = Percent390Gq[sortg]

########FINAL ARRAYS AFTER CUTTING#############
####The + 30.357 and 31.956 are the calibration.... they needed to add 32 to both magnitudes then subtract
####u-1.643 and b-0.044, I just added 32 minus the subtraction all in one go.
Name = Name300Qg
F300 = F300Qg
F300 = F300 + 30.357
F390 = F390Qg
F390 = F390 + 31.956
Q300 = Q300Qg
Q390 = Q390Qg
O300 = O300QOg
O390 = O390QOg
G300 = G300Gq
G390 = G390Gq
Percent300q = Percent300Qg
Percent390q = Percent390Qg
Percent300o = Percent300QOg
Percent390o = Percent390QOg
Percent300g = Percent300QOG
Percent390g = Percent390QOG

######################## CUT OUT EVERYTHING BUT THE WD SEQUENCE ###############
#F300 = F300[(F300>18.5)&(F300<26.5)]
#F390 = F390[(F300>18.5)&(F300<26.5)]

BV = F300 - F390

points = np.column_stack([BV, F300])

verts = np.array([[-2.0,-2.0,0.3,0.3,-2.0], [18.0,26.5,26.5,18.0,18.0]]).T
path = mpath.Path(verts)
points_inside = points[path.contains_points(points)]

length = str(len(points_inside))

##plt.plot([-2.0,0.5,0.5,-2.0,-2.0], [-5.0,-5.0,-14.0,-14.0,-5.0], color='r', linewidth=1, linestyle='--',label='MS = '+length)
#plt.plot(points_inside[:,0], points_inside[:,1],'k.', markersize=1)


##########################PLOTTING CMDS########################################

##Plot the CMD
#BV= F300-F390
#plt.figure(figsize=(10,10))
#plt.plot(BV,F300,'k.',markersize=1)
##plt.ylim(0,-14)
#plt.xlim(-1,4)
#plt.xlabel('F300X-F390W', fontsize=14) 
#plt.ylabel('F300X', fontsize=14)
#plt.tight_layout() 
#plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/Thesis/Figures/CMD.eps', format='eps', dpi=800)

######which stars are in this box?#######

#Box = Name[(F300>-11.76) & (F300<-9.59) & (BV>-0.42) & (BV<0.08)]

#----------Qfit vs mag-------------
#plt.figure(figsize=(10,10))
#plt.plot(F300,Q300,'k.',markersize=1)
#plt.xlabel('F300')
#plt.ylabel('Q300')
#plt.tight_layout()
#plt.figure(figsize=(10,10))
#plt.plot(F390,Q390,'k.',markersize=1)
#plt.xlabel('F390')
#plt.ylabel('Q390')
#plt.tight_layout()

#----------Ofit vs mag-------------
#plt.figure(figsize=(10,10))
#plt.plot(F300,O300,'k.',markersize=1)
#plt.xlabel('F300')
#plt.ylabel('O300')
#plt.ylim(0,20)
#plt.tight_layout()
#plt.figure(figsize=(10,10))
#plt.plot(F390,O390,'k.',markersize=1)
#plt.xlabel('F390')
#plt.ylabel('O390')
#plt.ylim(0,20)
#plt.tight_layout()


#--------Creating CMD subplots!-----------

#f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#ax1.plot(BV, y)
#ax1.set_title('Sharing Y axis')
#ax2.scatter(x, y)

#---------COOLING TRACKS--------------------
track172U, track172B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/serenelli_HeWD_models/z001/0172.wfc3x',usecols=[23,16], skiprows=3, unpack=True)
track300U, track300B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/serenelli_HeWD_models/z001/0300.wfc3x',usecols=[23,16], skiprows=3, unpack=True)
track449U, track449B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/serenelli_HeWD_models/z001/0449.wfc3x',usecols=[23,16], skiprows=3, unpack=True)

#--------RIVERA COOLING TRACKS--------------
track160U, track160B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/RiveraTracks/m0160_z003.txt',usecols=[7,8], skiprows=3, unpack=True)
track175U, track175B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/RiveraTracks/m0175_z003.txt',usecols=[7,8], skiprows=3, unpack=True)
track200U, track200B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/RiveraTracks/m0200_z003.txt',usecols=[7,8], skiprows=3, unpack=True)
track250U, track250B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/RiveraTracks/m0250_z003.txt',usecols=[7,8], skiprows=3, unpack=True)
Ntrack300U, Ntrack300B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/RiveraTracks/m0300_z003.txt',usecols=[7,8], skiprows=3, unpack=True)
track450U, track450B = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/47tuc/RiveraTracks/m0450_z003.txt',usecols=[7,8], skiprows=3, unpack=True)

#Convert tracks to apparent magnitude
tr172U = track172U + 0.26 + 13.36
tr300U = track300U + 0.26 + 13.36
tr449U = track449U + 0.26 + 13.36
tr172B = track172B + 0.18 + 13.36
tr300B = track300B + 0.18 + 13.36
tr449B = track449B + 0.18 + 13.36

tr160U = track160U + 0.26 + 13.36
tr175U = track175U + 0.26 + 13.36
tr200U = track200U + 0.26 + 13.36
tr250U = track250U + 0.26 + 13.36
Ntr300U = Ntrack300U + 0.26 + 13.36
tr450U = track450U + 0.26 + 13.36
tr160B = track160B + 0.18 + 13.36
tr175B = track175B + 0.18 + 13.36
tr200B = track200B + 0.18 + 13.36
tr250B = track250B + 0.18 + 13.36
Ntr300B = Ntrack300B + 0.18 + 13.36
tr450B = track450B + 0.18 + 13.36

tBV172 = tr172U - tr172B
tBV300 = tr300U - tr300B
tBV449 = tr449U - tr449B

tBV160 = tr160U - tr160B
tBV175 = tr175U - tr175B
tBV200 = tr200U - tr200B
tBV250 = tr250U - tr250B
NtBV300 = Ntr300U - Ntr300B
tBV450 = tr450U - tr450B
#BV = F300 - F390

plt.plot(BV,F300,'k.',markersize=1)
#plt.plot(points_inside[:,0], points_inside[:,1],'k.', markersize=1)
#plt.plot(tBV172,tr172U,'r')
#plt.plot(tBV300,tr300U,'r')
#plt.plot(tBV449,tr449U,'r')
#plt.xlabel('plot 1')
#
##plt.ylim(17,35)
##plt.xlim(-1.5,2)
##plt.ylim(-14,-5)
plt.xlim(-4,3)
plt.ylim(17,28)
plt.gca().invert_yaxis()  


######TRY STRAIGHTENING THE WD SEQUENCE#####

tBV = tr449U - tr449B

U449 = tr449U[(tr449U>18.0)&(tr449U<26.0)]#&(tBV>-0.7)&(tBV<0.7)]
B449 = tr449B[(tr449U>18.0)&(tr449U<26.0)]#&(tBV>-0.7)&(tBV<0.7)]
tBV = U449 - B449

F300cut = F300[(F300>18.0)&(F300<26.0)]#&(tBV>-0.7)&(tBV<0.7)]
F390cut = F390[(F300>18.0)&(F300<26.0)]#&(tBV>-0.7)&(tBV<0.7)]
FBV = F300cut - F390cut

#plt.figure()
#plt.plot(FBV,F300cut,'k.')
#plt.gca().invert_yaxis() 

#FBV = FBV[(FBV>-0.7)&(FBV<0.7)]
#F300cut = F300cut[(FBV>-0.7)&(FBV<0.7)]
#F390cut = F390cut[(FBV>-0.7)&(FBV<0.7)]

f = interpolate.interp1d(tr449U,tBV449-0.1,fill_value="extrapolate")
tBV_interp = f(F300cut) 

#plt.figure()
#plt.plot(tBV_interp,F300cut,'r.', markersize=2)    # show interpolation...   WORKS :)
# plot up the straightened sequence...
bvdiff = FBV - tBV_interp 
#tbvdiff = tBV449 - tBV_interp
#plt.plot(bvdiff,F300cut,'k.')
###PLOT WILL GIVE FIRST VIEW OF STRAIGHTENED SEQUENCE
#plt.plot(bvdiff, F300cut, 'k.', markersize=2)
#plt.plot([0, 0], [26, 18], 'r-', lw=2)
#plt.xlabel('first straightened sequence')
##plt.plot(bvdiff[abs(bvdiff)<1.0], F300cut[abs(F300cut)<1.0], 'b.', markersize=4) 
#plt.gca().invert_yaxis() 

#THIS CUTS THE CMD TO THE IMPORTANT B-V AND B VALUES
bvdiffcut = bvdiff[abs(bvdiff)<0.7]
F300cut = F300cut[abs(bvdiff)<0.7]
F390cut = F390cut[abs(bvdiff)<0.7]
bvcut = F300cut - F390cut
tBV_interpCut =  tBV_interp[abs(bvdiff)<0.7]     # need this later on to make plots etc...

#THis is making some headway
#plt.figure()
plt.plot(tBV_interpCut,F300cut,'r.', markersize=1)                            
                            
plt.figure()  ####The lines 1105,1212,1228,1231 plot a unstraightened sequence with line through
###THIS PLOTS A ZOOMED IN VERSION OF THE ONE ABOVE
#plt.figure(figsize=(10,10))
plt.plot(bvdiffcut,F300cut,'k.', markersize=1)
#plt.plot([0, 0], [26, 18], 'r-', lw=2)
##plt.plot([-0.7, 0.7], [20, 20 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [20.2, 20.2 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [20.4, 20.4 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [20.6, 20.6 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [20.8, 20.8 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [22.0, 22.0 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [22.2, 22.2 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [22.4, 22.4 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [22.6, 22.6 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [22.8, 22.8 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [24.0, 24.0 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [24.2, 24.2 ], 'b-', lw=1)
##plt.plot([-0.7, 0.7], [24.4, 24.4 ], 'b-', lw=1)
##plt.xlabel('first straightened sequence 2')
##plt.gca().invert_yaxis() 

####################binning and straightening tweak############################
##########Adrienne's code##########

###TOP BIN B is 18-22###
dbrm = bvdiffcut[(F300cut>18)&(F300cut<=19)&(bvdiffcut<0.7)]
rm = F300cut[(F300cut>18)&(F300cut<=19)&(bvdiffcut<0.7)]

flag_2s = np.zeros(len(dbrm))

for j in range(20):    # try up to 20 iterations...
   #This gets the average and standard deviation
   avg = np.mean(dbrm[flag_2s==0], axis=None, dtype=None, out=None)
   stdev = np.std(dbrm[flag_2s==0], axis=None, dtype=None, out=None, ddof=0)

#   print '{0:7.0f}{1:9.3f}{2:9.3f}{3:7.0f}{4:7.0f}'.format(j, avg, stdev, len(flag_2s[flag_2s==1]), len(flag_2s[flag_2s==0]))
   #This flags the stars that are outside of the sigma denoted previously
   for k in range(len(dbrm)):
      if ((dbrm[k]>(avg + 2.5*stdev))|(dbrm[k]<(avg - 2.5*stdev))):
         flag_2s[k] = 1
   
#plt.figure()
#plt.plot(dbrm[dbrm>(avg + stdev)], rm[dbrm>(avg + stdev)], 'c^', markersize=5)
##plt.plot(dbrm[flag_2s==1], rm[flag_2s==1], 'r.', markersize=4)
#plt.plot([0, 0], [22, 18], 'g-', lw=1)
##plt.plot(avg,rm,'m^',markersize=5)
##plt.xlim([-.7,0.7])
##plt.xlabel('stars outside of 2.5 sigma in top bin')
##plt.gca().invert_yaxis() 

#####Rolling bins starts here###
rbin = np.arange(19.0,26.01,0.4)     

avg_save = np.zeros(len(rbin))
stdev_save = np.zeros(len(rbin))
nkept = np.zeros(len(rbin))
nrej = np.zeros(len(rbin))
fkept = np.zeros(len(rbin))

#plt.gca().invert_yaxis()
for i in range(len(rbin)):

# extract set of stars that are in the current bin... (using "m" for "member")
   dbrm = bvdiffcut[(F300cut>(rbin[i]-0.7))&(F300cut<=(rbin[i]+0.7))]
   rm = F300cut[(F300cut>(rbin[i]-0.7))&(F300cut<=(rbin[i]+0.7))]

## special for bright end only:  exclude extreme outliers to the red side and use larger bins
#   if (rbin[i]<23.0):
#      dbrm = bvdiffcut[(F300cut>(rbin[i]-0.35))&(F300cut<=(rbin[i]+0.35))&(bvdiffcut<0.25)]
#      rm = F300cut[(F300cut>(rbin[i]-0.35))&(F300cut<=(rbin[i]+0.35))&(bvdiffcut<0.25)]

# define flag that I'll use to flag stars that don't make the 2sig cut
   flag_2s = np.zeros(len(dbrm))

   for j in range(15):    # try up to 15 iterations...

      avg = np.mean(dbrm[flag_2s==0], axis=None, dtype=None, out=None)
      stdev = np.std(dbrm[flag_2s==0], axis=None, dtype=None, out=None, ddof=0)

      for k in range(len(dbrm)):
         if ((dbrm[k]>(avg + 2.5*stdev))|(dbrm[k]<(avg - 2.5*stdev))):
            flag_2s[k] = 1
                   
#      print '{0:9.0f}{1:9.2f}{2:7.0f}{3:7.0f}{4:9.3f}{5:9.3f}{6:7.0f}{7:7.0f}'.format(i, rbin[i], len(dbrm), j, avg, stdev, len(flag_2s[flag_2s==1]), len(flag_2s[flag_2s==0]))
#      plt.plot(dbrm[dbrm>(avg + 3.5*stdev)], rm[dbrm>(avg + 3.5*stdev)], 'c^', markersize=5)
#      plt.plot(dbrm[flag_2s==1], rm[flag_2s==1], 'r.', markersize=4)
#      plt.plot([0, 0], [26, 20], 'r-', lw=2)
#      plt.xlabel('stars outside of 2.5 sigma in rolling bins')
#      plt.plot(dbrm[flag_2s==1], rm[flag_2s==1], 'm.', markersize=4)

   avg_save[i] = avg
   stdev_save[i] = stdev
   nkept[i] = len(flag_2s[flag_2s==0])
   nrej[i] = len(flag_2s[flag_2s==1])
   fkept[i] = float(nkept[i])/float(len(dbrm))

##plt.figure()   
#plt.plot(avg_save,rbin,'m^',markersize=5)
#plt.plot(avg_save + 2.5*stdev_save,rbin,'m-',markersize=3)
#plt.plot(avg_save - 2.5*stdev_save,rbin,'m-',markersize=3) 
#plt.plot([0, 0], [26, 18], 'r-', lw=2)
#plt.gca().invert_yaxis()
    

#FITTING A CURVE          
coeff = np.polyfit(rbin,avg_save,6)    # where "13" is order of polynomial
Poly = np.poly1d(coeff)                # this defines the polynomial with those coefficients
avg_fit = Poly(F300cut) 
                 # to evaluate the fit at a bunch of (other) x values
#plt.figure()
plt.plot(avg_fit,F300cut,'c.',markersize=2)  # plots the average fit of the sequence
              
dbruf = bvdiffcut

for i in range(len(bvdiffcut)):
   dbruf[i] = bvdiffcut[i] #- 0.0008
   if (F300cut[i]>=20):
      dbruf[i] = bvdiffcut[i] - avg_fit[i]


#plt.figure()   
#plt.plot(avg_save,rbin,'r-', lw=2)
#plt.plot(avg_save,rbin,'b^',markersize=5)
##plt.plot(rbin,avg_fit,'r-', lw=2)
#plt.plot(avg_save + 2.5*stdev_save,rbin,'m-',markersize=3)
#plt.plot(avg_save - 2.5*stdev_save,rbin,'m-',markersize=3) 
plt.plot([0, 0], [26, 20], 'g-', lw=2)
#plt.xlabel('F300X-F390W',fontsize=14)
#plt.ylabel('F300X',fontsize=14)
plt.gca().invert_yaxis()
#plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/Thesis/Thesis/WDsequencecurvefit.eps', format='eps', dpi=800)   

####################################################################################
##########################PLOT THIS FOR THE STRAIGHTENED SEQUENCE!!!!!###################
# let's see how it looks...
##clf()
#plt.figure(figsize=(10,10))
#plt.plot(dbruf, F300cut, 'k.', markersize=2)
#####plt.xlim([-3,3])
#####plt.xlim([-1.5,1.5])
####plt.ylim([27,19])
#plt.xlim([-0.7,0.7])
#plt.ylim([26.5,18.0])
#plt.xlabel('F300X - F390W',fontsize=14)
#plt.ylabel('F300X',fontsize=14)
#plt.plot([0, 0], [27, 18], 'r-', lw=1) 
#plt.tight_layout()  
#plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/Thesis/Thesis/StraightenedWDs.eps', format='eps', dpi=800)   

##NOW ON TO 2ND ITERATION, IN WHICH I WILL COMPUTE SIGMA 
## ...and experiment with what sigma clipping works well...
## TOP BIN: --------  r = 19.5-22   --------------------------------------------------------

# exclude dbr>0.25 (OBVIOUS ouliers; won't coverge otherwise)
# 3.0 sigma clip ==>  0/12 stars excluded ==> sigma = 0.018
# 2.5 sigma clip ==>  1/12 stars excluded ==> sigma = 0.011    
# 2.2 sigma clip ==>  1/12 stars excluded ==> sigma = 0.011    
# going with 2.2sig, since that's what I wound up needing in 22-22.5 range below

########To match Adrienne's variables######

ru = F300cut

###########################################

# extract straightened colors for stars in the current bin... (using "m" for "member")
dbrm = dbruf[(ru>18)&(ru<=20)&(dbruf<0.7)]
rm = ru[(ru>18)&(ru<=20)&(dbruf<0.7)]
# define flag that I'll use to flag stars that don't make the 2sig cut
flag_2s = np.zeros(len(dbrm))

for j in range(15):    # try up to 15 iterations...

   stdev = np.std(dbrm[flag_2s==0], axis=None, dtype=None, out=None, ddof=0)

#   print '{0:7.0f}{1:9.3f}{2:7.0f}{3:7.0f}'.format(j, stdev, len(flag_2s[flag_2s==1]), len(flag_2s[flag_2s==0]))

   for k in range(len(dbrm)):
      if ((dbrm[k]>(3.0*stdev))|(dbrm[k]<(-3.0*stdev))):
         flag_2s[k] = 1

#plt.plot(dbrm[flag_2s==1], rm[flag_2s==1], 'r.', markersize=4)

# show which stars are outside 3 sig (potential He WD candidates):  

#plt.plot(dbrm[dbrm>(3*stdev)], rm[dbrm>(3*stdev)], 'c^', markersize=4)

# setting up for later when I want to join result for this top bin with all the others...
toprbin = np.arange(18.0,20.0,0.1)
topstdev = np.zeros(len(toprbin)) + stdev


# ROLLING BINS: --------  r = 22-26   --------------------------------------------------------

# for r = 22-23, where numbers of stars per bin is very low (~8-30), I'm using 2.2 sig clipping (prob ~ 1/36)
# for r > 23 I'll switch to 3.0 sig clip
#plt.figure()
rbin = np.arange(20.0,26.2,0.1)     

stdev_save = np.zeros(len(rbin))
nkept = np.zeros(len(rbin))
nrej = np.zeros(len(rbin))
fkept = np.zeros(len(rbin))

for i in range(len(rbin)):

# extract set of stars that are in the current bin... (using "m" for "member")
   dbrm = dbruf[(ru>(rbin[i]-0.2))&(ru<=(rbin[i]+0.2))]
   rm = ru[(ru>(rbin[i]-0.2))&(ru<=(rbin[i]+0.2))]

# special for bright end only:  exclude extreme outliers to the red side and use larger bins
#   if (rbin[i]<23.0):
#      dbrm = dbruf[(ru>(rbin[i]-0.1))&(ru<=(rbin[i]+0.1))&(dbruf<0.1)]
#      rm = ru[(ru>(rbin[i]-0.1))&(ru<=(rbin[i]+0.1))&(dbruf<0.1)]

# define flag that I'll use to flag stars that don't make the 2sig cut
   flag_2s = np.zeros(len(dbrm))

   for j in range(15):    # try up to 15 iterations...

      stdev = np.std(dbrm[flag_2s==0], axis=None, dtype=None, out=None, ddof=0)

      for k in range(len(dbrm)):

         if (rbin[i]<23.0):
            if ((dbrm[k]>(3.0*stdev))|(dbrm[k]<(-3.0*stdev))):
               flag_2s[k] = 1

         if (rbin[i]>=23.0):
            if ((dbrm[k]>(3.0*stdev))|(dbrm[k]<(-3.0*stdev))):
               flag_2s[k] = 1

#      print '{0:9.0f}{1:9.2f}{2:7.0f}{3:7.0f}{4:9.3f}{5:7.0f}{6:7.0f}'.format(i, rbin[i], len(dbrm), j, stdev, len(flag_2s[flag_2s==1]), len(flag_2s[flag_2s==0]))
#      plt.plot(dbrm[flag_2s==1], rm[flag_2s==1], 'c^', markersize=4)
#      plt.plot(dbrm[flag_2s==1], rm[flag_2s==1], 'r.', markersize=4)

   stdev_save[i] = stdev
   nkept[i] = len(flag_2s[flag_2s==0])
   nrej[i] = len(flag_2s[flag_2s==1])
   fkept[i] = float(nkept[i])/float(len(dbrm))
#   print '{0:9.0f}{1:9.2f}{2:7.0f}{3:9.3f}{4:7.0f}{5:7.0f}'.format(i, rbin[i], len(dbrm), stdev, len(flag_2s[flag_2s==1]), len(flag_2s[flag_2s==0]))

# plot of up sigma values...

#plt.plot(3*stdev_save,rbin,'m-',lw=1)
#plt.plot(-3*stdev_save,rbin,'m-',lw=1)

# combine these stdev results with value for top bin (r=19.5-22)...
# (looking good!)

allrbin = np.concatenate((toprbin,rbin), axis=0)
allstdev = np.concatenate((topstdev,stdev_save), axis=0)

#plt.plot(2*allstdev,allrbin,'r-',lw=2)
#plt.plot(-2*allstdev,allrbin,'r-',lw=2)

# fit a smooth curve to allrbin, allstdev

space = np.linspace(20.6,26)
coefficient,residual = poly.polyfit(allrbin, allstdev, 6,full=True)
ffit = poly.polyval(space, coefficient)


f_fitting = np.array([0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048,0.048])             
      
newstdev = np.concatenate((f_fitting,ffit))
newspace = np.linspace(18,26,num=72)

#plt.figure()
#plt.plot(allrbin,allstdev,'k.')
#plt.plot(newspace,newstdev,'g-',lw=2)
#plt.plot(allrbin,stdev_fitting,'r-',lw=2)
#plt.xlabel('Mag (F300X)')
#plt.ylim([0,0.35])
#plt.ylabel(unicodedata.lookup("GREEK SMALL LETTER SIGMA"))
#plt.xlim(([17,27]))
##plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/Thesis/Thesis/magsigma.eps', format='eps', dpi=800)   


# now fit a smooth curve to stdev values.....

##################try to unstraighten the data#################

#bvdiffrestraight = bvdiff + tBV_interpCut

#plt.plot(bvdiffrestraight,F300)

########################################################
####This section is old#######
coeff = np.polyfit(allrbin,allstdev,6)    
Poly = np.poly1d(coeff)                # this defines the polynomial with those coefficients
stdev_fit = Poly(ru)                  # to evaluate the fit at a bunch of (other) x values
     
#plt.figure()
#plt.plot(stdev_fit,ru,'r.',markersize=1)  
#plt.plot(2.0*stdev_fit,ru,'b.',markersize=1)  
#plt.plot(3.0*stdev_fit,ru,'k.',markersize=1)  
#plt.plot(3.5*stdev_fit,ru,'m.',markersize=1)  
#plt.plot(4.0*stdev_fit,ru,'c.',markersize=1)  
#plt.plot(-stdev_fit,ru,'r.',markersize=1)  
#plt.plot(-2.0*stdev_fit,ru,'b.',markersize=1)  
#plt.plot(-3.0*stdev_fit,ru,'k.',markersize=1) 
#plt.plot(-3.5*stdev_fit,ru,'m.',markersize=1)
#plt.plot(-4.0*stdev_fit,ru,'c.',markersize=1)
#
#plt.plot(bvdiff,F300cut,'k.')
###############################
###########################################################
####This plots the straightened sequence with whatever sigma line you want
#plt.figure()
#plt.plot(dbruf, F300cut, 'k.', markersize=1)
#plt.plot(newstdev,newspace,'m-',markersize=1,label='1'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))  
#plt.plot(2.0*newstdev,newspace,'b-',markersize=1,label='2'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))  
#plt.plot(3.0*newstdev,newspace,'c-',markersize=1,label='3'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))  
###plt.plot(3.5*newstdev,newspace,'c-',markersize=2)  
#plt.plot(4.0*newstdev,newspace,'g-',markersize=1,label='4'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))  
#plt.plot(-newstdev,newspace,'m-',markersize=1)  
#plt.plot(-2.0*newstdev,newspace,'b-',markersize=1)  
#plt.plot(-3.0*newstdev,newspace,'c-',markersize=1) 
###plt.plot(-3.5*newstdev,newspace,'c-',markersize=2)
#plt.plot(-4.0*newstdev,newspace,'g-',markersize=1) 
#plt.plot([0, 0], [26, 18], 'r-', lw=1)  
#plt.gca().invert_yaxis() 
#plt.xlabel('F300X - F390W',fontsize=14)
#plt.ylabel('F300X',fontsize=14)
#plt.legend()
#plt.tight_layout()
#plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/Thesis/Thesis/CMDgaussian.eps', format='eps', dpi=800)   

## count up HeWD candidates (ignoring if too red for now...)

#inter = np.linspace(18,26,num=2469)
#newfits = np.interp(inter,newspace,newfit)
#inter = np.linspace(18,26,num=89521)
#newfits = np.interp(inter,newspace,newfit)

#newfits = newfits + tBV_interpCut
#plt.plot(newfits,ru,'r.',markersize=2)  

#print len(ru[dbruf>(4*stdev_fit)])
#print len(ru[dbruf>(3*stdev_fit)])      # 86 stars make this cut... NICE :)
#print len(ru[dbruf>(2*stdev_fit)])      # 138 stars make this cut (extra 52 stars in 2-3sig range)
#print len(ru[dbruf>(1*stdev_fit)])
#print len(ru[dbruf<=(1*stdev_fit)])

###########LETS TRY GETTING THE MIDDLE LINE CENTERED ON THE UNSTRAIGHTENED SEQ#######


seq = [tBV_interpCut,avg_fit]
centerline = np.sum(seq,axis=0)
####This plots the centerline moved over
#plt.figure()
#plt.plot(BV,F300,'k.',markersize=1)
#plt.plot(centerline,F300cut,'r.',markersize=1)
#plt.plot(newstdev,newspace,'r-',markersize=2)  
#plt.xlabel('centerline')
#plt.xlim(-4,3)
#plt.ylim(17,28)
#plt.gca().invert_yaxis() 

##### WE NEED TO ADJUST THE SIGMA LINES TO MATCH THE CURVE OF THE CENTERLINE
#####TO DO THIS WE NEED TO INFLICT THE CHANGES FROM THE OTHER TWO ADJUSTMENT LINES
##### ON THE SIGMA LINES
### PERHAPS WE CAN DO AN INTERP LIKE THE FIRST ONE BUT FROM THE STRAIGHT LINE AT
### ZERO TO THE CENTERLINE....

#Need to match the length of tBV160 to newstdev

#smallstdev = allstdev[::3.28]
#smallstdev = np.resize(smallstdev,(25))
#
##BV160_move_1sig = tBV160 + smallstdev
#BV160_move_neg1sig = tBV160 - smallstdev


#lets get newstdev to match the same length as the other two

#1sigma positive side
g = interpolate.interp1d(newspace,newstdev,fill_value="extrapolate")
newstdfit = g(F300cut) 
stdev_move_1sig = newstdfit + centerline

#1sigma negative side
h = interpolate.interp1d(newspace,-newstdev,fill_value="extrapolate")
newstdfitneg = h(F300cut) 
stdev_move_neg1sig = newstdfitneg + centerline

#2sigma positive side
i = interpolate.interp1d(newspace,2*newstdev,fill_value="extrapolate")
newstdfit2pos = i(F300cut) 
stdev_move_2sig = newstdfit2pos + centerline

#2sigma negative side
j = interpolate.interp1d(newspace,-2*newstdev,fill_value="extrapolate")
newstdfit2neg = j(F300cut) 
stdev_move_neg2sig = newstdfit2neg + centerline

#3sigma positive side
k = interpolate.interp1d(newspace,3*newstdev,fill_value="extrapolate")
newstdfit3pos = k(F300cut) 
stdev_move_3sig = newstdfit3pos + centerline

#3sigma negative side
l = interpolate.interp1d(newspace,-3*newstdev,fill_value="extrapolate")
newstdfit3neg = l(F300cut) 
stdev_move_neg3sig = newstdfit3neg + centerline

#4sigma positive side
k = interpolate.interp1d(newspace,4*newstdev,fill_value="extrapolate")
newstdfit4pos = k(F300cut) 
stdev_move_4sig = newstdfit4pos + centerline

#4sigma negative side
l = interpolate.interp1d(newspace,-4*newstdev,fill_value="extrapolate")
newstdfit4neg = l(F300cut) 
stdev_move_neg4sig = newstdfit4neg + centerline

###need to turn the 0.16 cooling track into function to get the stars between
###

coeff = np.polyfit(tr160U,tBV160,8)    
Poly = np.poly1d(coeff)                # this defines the polynomial with those coefficients
trackfunction = Poly(ru) 

BV160_move_1sig = trackfunction + newstdfit
BV160_move_neg1sig = trackfunction - newstdfit


####This plots all of the sigma lines onto the CMD around the WD sequence
plt.figure()
##plt.plot(newstdfit,F300cut,'y.',markersize=1)
##plt.plot(newstdfitneg,F300cut,'y.',markersize=1)
plt.plot(BV,F300,'k.',markersize=1)
plt.plot(stdev_move_1sig,F300cut,'m.',markersize=0.5)
plt.plot(stdev_move_neg1sig,F300cut,'m.',markersize=0.5)
plt.plot(stdev_move_2sig,F300cut,'b.',markersize=0.5)
plt.plot(stdev_move_neg2sig,F300cut,'b.',markersize=0.5)
plt.plot(stdev_move_3sig,F300cut,'c.',markersize=0.5)
plt.plot(stdev_move_neg3sig,F300cut,'c.',markersize=0.5)
plt.plot(stdev_move_4sig,F300cut,'g.',markersize=0.5)
plt.plot(stdev_move_neg4sig,F300cut,'g.',markersize=0.5)
plt.plot(centerline,F300cut,'r.',markersize=0.5)
plt.plot(tBV160,tr160U,'y-',markersize=1)
##plt.plot(tBV450,tr450U,'y-',markersize=1)
plt.plot(BV160_move_1sig,F300cut,'m.',markersize=0.5)
plt.plot(BV160_move_neg1sig,F300cut,'m.',markersize=0.5)
##plt.plot(trackfunction,F300cut,'b.',markersize=1)
plt.xlim(-2.3,1.3)
plt.ylim(18,27)
#plt.xlabel('F300X - F390W',fontsize=14)
#plt.ylabel('F300X',fontsize=14)
plt.gca().invert_yaxis()
#plt.tight_layout()


####MAKE A HISTOGRAM OF STARS IN EACH SIGMA BIN
############################USE THIS ONE##########################
foursigma = len(ru[(dbruf>(4*newstdfit))&(dbruf<=(5*newstdfit))&(ru<=24.1)])
threesigma = len(ru[(dbruf>(3*newstdfit))&(dbruf<=(4*newstdfit))&(ru<=24.1)])      
twosigma = len(ru[(dbruf>(2*newstdfit))&(dbruf<=(3*newstdfit))&(ru<=24.1)])
onesigma = len(ru[(dbruf>(1*newstdfit))&(dbruf<=(2*newstdfit))&(ru<=24.1)])
zerosigma = len(ru[(dbruf<=(1*newstdfit))&(ru<=24.1)])

foursigmaneg = len(ru[(dbruf<(-4*newstdfit))&(dbruf>=(-5*newstdfit))&(ru<=24.1)])
threesigmaneg = len(ru[(dbruf<(-3*newstdfit))&(dbruf>=(-4*newstdfit))&(ru<=24.1)])      
twosigmaneg = len(ru[(dbruf<(-2*newstdfit))&(dbruf>=(-3*newstdfit))&(ru<=24.1)])
onesigmaneg = len(ru[(dbruf<(-1*newstdfit))&(dbruf>=(-2*newstdfit))&(ru<=24.1)])
zerosigmaneg = len(ru[(dbruf>=(-1*newstdfit))&(ru<=24.1)])

positiveside = [zerosigma,onesigma,twosigma,threesigma,foursigma]
negativeside = [zerosigmaneg,onesigmaneg,twosigmaneg,threesigmaneg,foursigmaneg]
bothsides = np.array([foursigmaneg,threesigmaneg,twosigmaneg,onesigmaneg,zerosigmaneg,zerosigma,onesigma,twosigma,threesigma,foursigma])

#bins = np.linspace(-4, 4, 10)
#negbins = np.linspace(-4,0, 5)

###This plot show the same as below, but in a different way
#plt.figure()
#N= 5
#number = np.arange(N)
#p1=plt.bar(number,positiveside)
#p2=plt.bar(number,negativeside,bottom=positiveside)
#plt.legend((p1[0], p2[0]), ('Right', 'Left'))
#pyplot.hist(bothsides, bins, alpha=0.5)
#plt.xlabel(unicodedata.lookup("GREEK SMALL LETTER SIGMA"))

###This plot is for the histogram of stars within each sigma around WD sequence
#plt.figure()
#plt.yscale('log')
#bar_width = 1
#positions = np.arange(-5,5)
#plt.bar(positions, bothsides, bar_width,alpha=0.7,edgecolor='k',align='edge')
#plt.xlabel(unicodedata.lookup("GREEK SMALL LETTER SIGMA"),fontsize=13)
#plt.ylabel('Number of Sources', fontsize=10)
#plt.tight_layout()

#####Give total number of stars within the boundaries
####Bin the y from 23-24.5 to find how many stars are within each bin
#Need to convert these lines to have the same shape as F300 and BV

coeff = np.polyfit(F300cut,stdev_move_4sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_4sig = Poly(F300)

coeff = np.polyfit(F300cut,stdev_move_3sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_3sig = Poly(F300)

coeff = np.polyfit(F300cut,stdev_move_2sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_2sig = Poly(F300)

coeff = np.polyfit(F300cut,stdev_move_1sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_1sig = Poly(F300)

coeff = np.polyfit(F300cut,stdev_move_neg1sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_neg1sig = Poly(F300)

coeff = np.polyfit(F300cut,stdev_move_neg2sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_neg2sig = Poly(F300)

coeff = np.polyfit(F300cut,stdev_move_neg3sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_neg3sig = Poly(F300)

coeff = np.polyfit(F300cut,stdev_move_neg4sig,11)    
Poly = np.poly1d(coeff)                
newstdev_move_neg4sig = Poly(F300)

coeff = np.polyfit(F300cut,centerline,11)    
Poly = np.poly1d(coeff)                
newcenterline = Poly(F300)

##This allows me to get the stars inbetween the curves

coeff = np.polyfit(tr160U,tBV160,8)    
Poly = np.poly1d(coeff)                
newtrackfunction = Poly(F300) 

coeff = np.polyfit(F300cut,BV160_move_1sig,8)
Poly = np.poly1d(coeff)
newtrackfunction1sig = Poly(F300)

#HeWDs outside 4
BVRange = BV[(BV>newstdev_move_4sig)&(BV<newtrackfunction1sig)&(F300<=24.1)&(BV<0.25)&(BV<-0.02)]
F300Range = F300[(BV>newstdev_move_4sig)&(BV<newtrackfunction1sig)&(F300<=24.1)&(BV<0.25)&(BV<-0.02)]
IDs = Name[(BV>newstdev_move_4sig)&(BV<newtrackfunction1sig)&(F300<=24.1)&(BV<0.25)&(BV<-0.02)]

#HeWDs between 3 and 4
BVRange34 = BV[(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)&(F300<=24.5)&(BV<0.25)&(BV<-0.02)]
F300Range34 = F300[(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)&(F300<=24.5)&(BV<0.25)&(BV<-0.02)]
IDs34 = Name[(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)&(F300<=24.5)&(BV<0.25)&(BV<-0.02)]

#COWD
BVRangeCO = BV[(BV>newstdev_move_neg2sig)&(BV<newstdev_move_2sig)&(F300<=24.1)]
F300RangeCO = F300[(BV>newstdev_move_neg2sig)&(BV<newstdev_move_2sig)&(F300<=24.1)]
IDsCO = Name[(BV>newstdev_move_neg2sig)&(BV<newstdev_move_2sig)&(F300<=24.1)]

#Turnoff stars
BVRangeTO = BV[(BV>0.45)&(BV<0.67)&(F300>18.5)&(F300<19.0)]
F300RangeTO = F300[(BV>0.45)&(BV<0.67)&(F300>18.5)&(F300<19.0)]
IDsTO = Name[(BV>0.45)&(BV<0.67)&(F300>18.5)&(F300<19.0)]

#possible BS?
BVRangeBS = BV[(BV<0.45)&(BV>0.25)&(F300>18.0)&(F300<18.62)]
F300RangeBS = F300[(BV<0.45)&(BV>0.25)&(F300>18.0)&(F300<18.62)]
IDsBS = Name[(BV<0.45)&(BV>0.25)&(F300>18.0)&(F300<18.62)]

#MS 20.5-21
BVRangeMS21 = BV[(BV<1.0)&(BV>0.7)&(F300>20.5)&(F300<21.0)]
F300RangeMS21 = F300[(BV<1.0)&(BV>0.7)&(F300>20.5)&(F300<21.0)]
IDsMS21 = Name[(BV<1.0)&(BV>0.7)&(F300>20.5)&(F300<21.0)]

#MS 22.5-23
BVRangeMS23 = BV[(BV<1.7)&(BV>1.2)&(F300>22.5)&(F300<23.0)]
F300RangeMS23 = F300[(BV<1.7)&(BV>1.2)&(F300>22.5)&(F300<23.0)]
IDsMS23 = Name[(BV<1.7)&(BV>1.2)&(F300>22.5)&(F300<23.0)]

#MS 24.5-25
BVRangeMS25 = BV[(BV<2.2)&(BV>1.5)&(F300>24.5)&(F300<25.0)]
F300RangeMS25 = F300[(BV<2.2)&(BV>1.5)&(F300>24.5)&(F300<25.0)]
IDsMS25 = Name[(BV<2.2)&(BV>1.5)&(F300>24.5)&(F300<25.0)]

#MS 26-26.5
BVRangeMS26 = BV[(BV<2.5)&(BV>1.0)&(F300>26.0)&(F300<26.5)]
F300RangeMS26 = F300[(BV<2.5)&(BV>1.0)&(F300>26.0)&(F300<26.5)]
IDsMS26 = Name[(BV<2.5)&(BV>1.0)&(F300>26.0)&(F300<26.5)]

#Test the chosen stars
#plt.figure()
#plt.plot(BV,F300,'k.',markersize=1)
#plt.plot(BVRangeTO,F300RangeTO,'r.',markersize=1)
#plt.plot(BVRangeMS21,F300RangeMS21,'r.',markersize=1)
#plt.xlim(-3,3)
#plt.ylim(17,28)
#plt.gca().invert_yaxis()


######FIGURE OUT HOW TO ISOLATE THE 12 COMPARISON SOURCES AND PLOT THEM
comparisonCVIndex = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/Regions/listCVIDS.txt',usecols=[1])
comparCVBV = []
for i in comparisonCVIndex:
    comparCVBV.append(BVRange)
#
#comparCVB = []
#for i in comparisonCVIndex:
#    comparCVB.append(F300Range[i])
#    
#comparisonMSPIndex = np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/Regions/listMSPIDS.txt',usecols=[1])
#comparMSPBV = []
#for i in comparisonMSPIndex:
#    comparMSPBV.append(BVRange[i])
#
#comparMSPB = []
#for i in comparisonMSPIndex:
#    comparMSPB.append(F300Range[i])
    
# 14 is the index from listotherIDS.txt
#comparotherBV = BVRange[14]
#comparotherB = F300Range[14]

#plt.figure()
#plt.plot(BV,F300,'k.')
plt.plot(BVRange,F300Range,'r*',markersize=3,label='He WDs >4'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))
plt.plot(BVRange34,F300Range34,'m*',markersize=3,label='"maybe" He WDs')
#plt.plot(comparCVBV,comparCVB,'c^',markersize=4,label='CVs')
#plt.plot(comparMSPBV,comparMSPB,'bs',markersize=3,label='MSPs')
#plt.plot(comparotherBV,comparotherB,'gD',markersize=3,label='X-ray source')
#plt.plot(BVRangeCO,F300RangeCO,'c*',markersize=3)
##plt.plot(newstdev_move_4sig,F300,'.',markersize=1)
#plt.legend(loc=4)
plt.xlabel('$m_{300}-m_{390}$',fontsize=14)
plt.ylabel('$m_{300}$',fontsize=14)
#plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/CMDcandidatesCVMSP.png', format='png', dpi=800)   

#plot the cooling tracks on the CMD
#plt.figure()
#plt.plot(BV,F300,'k.',markersize=0.8)
#plt.plot(tBV160,tr160U,'r-',label='0.16$M_{\odot}$')
#plt.plot(tBV175,tr175U,'r-',label='0.17$M_{\odot}$')
#plt.plot(tBV200[(tr200U<=26.5)],tr200U[(tr200U<=26.5)],'r-',label='0.2$M_{\odot}$')
#plt.plot(tBV250[(tr250U<=26.5)],tr250U[(tr250U<=26.5)],'r-')
#plt.plot(NtBV300[(Ntr300U<=26.5)],Ntr300U[(Ntr300U<=26.5)],'r-')
#plt.plot(tBV450[(tr450U<=26.5)],tr450U[(tr450U<=26.5)],'r-')
#plt.xlim(-3,3)
#plt.ylim(17,28)
#plt.gca().invert_yaxis()
#plt.xlabel('F300X - F390W',fontsize=14)
#plt.ylabel('F300X',fontsize=14)
#plt.text(-0.2, 20.2, '0.16$M_{\odot}$')
#plt.text(-0.4,18.6,'0.17$M_{\odot}$')
#plt.text(-1.5,17.6,'0.2$M_{\odot}$')
#plt.text(-2.26,17.6,'0.25$M_{\odot}$')
#plt.text(-2.28,18.1,'0.3$M_{\odot}$')
#plt.text(-2.3,18.5,'0.45$M_{\odot}$')
#plt.tight_layout()
#plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/Thesis/Thesis/coolingtracks.eps', format='eps', dpi=800)   

####Find the stars in the data file with the names from IDs and save as a txtfile
#Sources outside 4 sigma 
#sources=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDs)):
#        if IDs[i] in line: 
#            sources.append(line)
#searchfile.close()
#
##sources between 3 and 4 sigma
#sources34=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDs34)):
#        if IDs34[i] in line: 
#            sources34.append(line)
#searchfile.close()
#
#sourcesCO=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDsCO)):
#        if IDsCO[i] in line: 
#            sourcesCO.append(line)
#searchfile.close()
#
#sourcesHe = sources + sources34
#
#sourcesTO=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDsTO)):
#        if IDsTO[i] in line: 
#            sourcesTO.append(line)
#searchfile.close()
#
#sourcesBS=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDsBS)):
#        if IDsBS[i] in line: 
#            sourcesBS.append(line)
#searchfile.close()

#sourcesMS21=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDsMS21)):
#        if IDsMS21[i] in line: 
#            sourcesMS21.append(line)
#searchfile.close()

#sourcesMS23=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDsMS23)):
#        if IDsMS23[i] in line: 
#            sourcesMS23.append(line)
#searchfile.close()

#sourcesMS25=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDsMS25)):
#        if IDsMS25[i] in line: 
#            sourcesMS25.append(line)
#searchfile.close()
#
#sourcesMS26=[]
#searchfile = open("/Users/sarahdeveny/Research/GraduateThesis/47tuc/heinke_GO-12950/LOGR.XYVIQ3", "r")
#for line in searchfile:
#    for i in range(len(IDsMS26)):
#        if IDsMS26[i] in line: 
#            sourcesMS26.append(line)
#searchfile.close()

#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/candidates_new.txt',sources,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/candidates_maybe.txt',sources34,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/IDsNEW.txt',IDs,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/IDsmaybe.txt',IDs34,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/Allcandidates.txt',sourcesHe,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COcandidates.txt',sourcesCO,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/TOcandidates.txt',sourcesTO,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/BScandidates.txt',sourcesBS,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MScandidates21.txt',sourcesMS21,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MScandidates23.txt',sourcesMS23,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MScandidates25.txt',sourcesMS25,fmt='%s')
#np.savetxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MScandidates26.txt',sourcesMS26,fmt='%s')

#####Get distance from the center
#######Uncomment to calculate distances from the center of different populations
X_o=2964.3838
Y_o=2962.0067

##Candidates outside of 4 sigma
#x4sig=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandys/candidates_x_noCVs.txt',usecols=[0])
#y4sig=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandys/candidates_y_noCVs.txt',usecols=[0])
#
#k=len(x4sig)
#r4sig=np.zeros(k)
#for j in range(0,k):
#    r4sig[j]=(np.sqrt((X_o-x4sig[j])**2+(Y_o-y4sig[j])**2))
#
##Candidates outside 4 sigma with no 3's or 4's  !!!!!! we want this one
#xno344sig=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandyno34/candidates_no34_x_noCVs.txt',usecols=[0])
#yno344sig=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandyno34/candidates_no34_y_noCVs.txt',usecols=[0])
#
#k=len(xno344sig)
#rno344sig=np.zeros(k)
#for j in range(0,k):
#    rno344sig[j]=(np.sqrt((X_o-xno344sig[j])**2+(Y_o-yno344sig[j])**2))
##    print rno344sig
#
##All candidates outside of 3 sigma
#Allx=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandys/Allcandidates_x_noCVs.txt',usecols=[0])
#Ally=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandys/Allcandidates_y_noCVs.txt',usecols=[0])
#
#k=len(Allx)
#Allr=np.zeros(k)
#for j in range(0,k):
#    Allr[j]=(np.sqrt((X_o-Allx[j])**2+(Y_o-Ally[j])**2))
#
##Candidates outside of 3 sigma with no 3's or 4's !!!!!! we want this one
#Allxno34=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandyno34/Allcandidates_no34_x_noCVs.txt',usecols=[0])
#Allyno34=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/xandyno34/Allcandidates_no34_y_noCVs.txt',usecols=[0])
#
#k=len(Allxno34)
#Allrno34=np.zeros(k)
#for j in range(0,k):
#    Allrno34[j]=(np.sqrt((X_o-Allxno34[j])**2+(Y_o-Allyno34[j])**2))
# 
##COWD inside plus and minus 2 sigma    
#COx=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/COcandidates_x.txt',usecols=[0])
#COy=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/COcandidates_y.txt',usecols=[0])
#
#k=len(COx)
#rCO=np.zeros(k)
#for j in range(0,k):
#    rCO[j]=(np.sqrt((X_o-COx[j])**2+(Y_o-COy[j])**2))
#    
##COWD brighter 23
#COx23=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_x23.1.txt',usecols=[0])
#COy23=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_y23.1.txt',usecols=[0])
#
#k=len(COx23)
#rCO23=np.zeros(k)
#for j in range(0,k):
#    rCO23[j]=(np.sqrt((X_o-COx23[j])**2+(Y_o-COy23[j])**2))
#    
##COWD brighter 22
#COx22=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_x22.1.txt',usecols=[0])
#COy22=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_y22.1.txt',usecols=[0])
#
#k=len(COx22)
#rCO22=np.zeros(k)
#for j in range(0,k):
#    rCO22[j]=(np.sqrt((X_o-COx22[j])**2+(Y_o-COy22[j])**2))
#    
##COWD brighter 21
#COx21=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_x21.1.txt',usecols=[0])
#COy21=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_y21.1.txt',usecols=[0])
#
#k=len(COx21)
#rCO21=np.zeros(k)
#for j in range(0,k):
#    rCO21[j]=(np.sqrt((X_o-COx21[j])**2+(Y_o-COy21[j])**2))
#    
##COWD brighter 20
#COx20=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_x20.1.txt',usecols=[0])
#COy20=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/COWD/COcandidates_y20.1.txt',usecols=[0])
#
#k=len(COx20)
#rCO20=np.zeros(k)
#for j in range(0,k):
#    rCO20[j]=(np.sqrt((X_o-COx20[j])**2+(Y_o-COy20[j])**2))
#    
##Turn off stars
#TOx=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/TO/TOfainter/TOcandidates_x.txt',usecols=[0])
#TOy=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/TO/TOfainter/TOcandidates_y.txt',usecols=[0])
#
#k=len(TOx)
#rTO=np.zeros(k)
#for j in range(0,k):
#    rTO[j]=(np.sqrt((X_o-TOx[j])**2+(Y_o-TOy[j])**2))
#    
##BS
#BSx=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/BS/BScandidates_x.txt',usecols=[0])
#BSy=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/BS/BScandidates_y.txt',usecols=[0])
#
#k=len(BSx)
#rBS=np.zeros(k)
#for j in range(0,k):
#    rBS[j]=(np.sqrt((X_o-BSx[j])**2+(Y_o-BSy[j])**2))
#    
##MS 20.5-21
#MSx21=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates21_x.txt',usecols=[0])
#MSy21=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates21_y.txt',usecols=[0])
#
#k=len(MSx21)
#rMS21=np.zeros(k)
#for j in range(0,k):
#    rMS21[j]=(np.sqrt((X_o-MSx21[j])**2+(Y_o-MSy21[j])**2))    
#
##MS 22.5-23
#MSx23=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates23_x.txt',usecols=[0])
#MSy23=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates23_y.txt',usecols=[0])
#
#k=len(MSx23)
#rMS23=np.zeros(k)
#for j in range(0,k):
#    rMS23[j]=(np.sqrt((X_o-MSx23[j])**2+(Y_o-MSy23[j])**2)) 
#    
##MS 24.5-25
#MSx25=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates25_x.txt',usecols=[0])
#MSy25=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates25_y.txt',usecols=[0])
#
#k=len(MSx25)
#rMS25=np.zeros(k)
#for j in range(0,k):
#    rMS25[j]=(np.sqrt((X_o-MSx25[j])**2+(Y_o-MSy25[j])**2)) 
#
##MS 26-26.5
#MSx26=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates26_x.txt',usecols=[0])
#MSy26=np.loadtxt('/Users/sarahdeveny/Research/GraduateThesis/Candidates/TestingHist/MS/MScandidates26_y.txt',usecols=[0])
#
#k=len(MSx26)
#rMS26=np.zeros(k)
#for j in range(0,k):
#    rMS26[j]=(np.sqrt((X_o-MSx26[j])**2+(Y_o-MSy26[j])**2)) 
 
##Convert from pixel distance to arcsecond distance
#r4sig = r4sig*0.04
#rno344sig = rno344sig*0.04
#Allr = Allr*0.04
#Allrno34 = Allrno34*0.04
#rCO = rCO*0.04
#rTO = rTO*0.04
#rBS = rBS*0.04
#rMS21 = rMS21*0.04
#rMS23 = rMS23*0.04
#rMS25 = rMS25*0.04
#rMS26 = rMS26*0.04
#rCO23 = rCO23*0.04
#rCO22 = rCO22*0.04
#rCO21 = rCO21*0.04
#rCO20 = rCO20*0.04
#
#### KS test #########
#
#from scipy.stats import ks_2samp
#ks_2samp(r4sig,rCO)
#
##from scipy.stats import uniform
#
#
##Uniform distribution
#np.random.seed(19680801)
#
#mu = 50
#sigma = 25
#n_bins = 10000
#x = np.random.uniform(0,100,10000)


##This shows a uniform area, but I don't know how it works
#p= 2* rand(3, 1e4)- 1
#p= p[:, sum(p* p, 0)** .5<= 1]
#p.shape
#plt.plot(p[0], p[2], '.')

#length = np.sqrt(np.random.uniform(0, 1))
#angle = np.pi * np.random.uniform(0, 2)
#
#x = length * np.cos(angle)
#y = length * np.sin(angle)
#
#r = x**2 + y**2
#
#plt.figure()
#plt.scatter(x,y)

#y = np.random.uniform(0,100,1000000)
#r = x**2 + y**2

#x = np.random.normal(mu, sigma, size=10000)

#Cumulative histogram
#bins = np.linspace(0,100,10000)
#plt.figure()
#plt.hist(r4sig,bins,cumulative=True,normed=1,histtype='step',label='He WDs >4'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))
#plt.hist(rno344sig,bins,cumulative=True,normed=1,histtype='step',label='He WDs >4'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))#+', no 3 or 4') 
#plt.hist(Allr,bins,cumulative=True,normed=1,histtype='step',label='He WDs >3'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"))
#plt.hist(Allrno34,bins,cumulative=True,normed=1,histtype='step',label='He WDs >3'+unicodedata.lookup("GREEK SMALL LETTER SIGMA"),lw=2)#+', no 3 or 4',lw=2)
#plt.hist(rCO23,bins,cumulative=True,normed=1,histtype='step',label='CO WD $m_{300}$<23.1',lw=1)
#plt.hist(rTO,bins,cumulative=True,normed=1,histtype='step',label='Turn-off',lw=1)
#plt.hist(rMS21,bins,cumulative=True,normed=1,histtype='step',label='20.5<MS<21')
#plt.hist(rMS23,bins,cumulative=True,normed=1,histtype='step',label='22.5<MS<23')
#plt.hist(rMS25,bins,cumulative=True,normed=1,histtype='step',label='24.5<MS<25')
#plt.hist(rMS26,bins,cumulative=True,normed=1,histtype='step',label='26<MS<26.5')
#plt.hist(x,bins,cumulative=True,normed=1,histtype='step',label='Uniform',color='black')
#plt.hist(rBS,bins,cumulative=True,normed=1,histtype='step',label='BS')
#plt.xlabel('Distance from cluster center (arcseconds)')
#plt.ylabel('Fraction of the population')
#plt.xlim(0, 100)
#plt.legend(loc=4)
#plt.tight_layout()
#plt.savefig('/Users/sarahdeveny/Research/GraduateThesis/HistCO<23.eps', format='eps', dpi=800) 

####x and y plot

#fig=plt.figure()
#ax=fig.add_subplot(1,1,1)
#circ=plt.Circle((2964,2962), radius=540, color='g', fill=False)
#ax.add_patch(circ)
#plt.plot(xno344sig,yno344sig,'k.')
#plt.plot(Allxno34,Allyno34,'k.')
#plt.plot(COx,COy,'b.',markersize=1)
#plt.axis('equal')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.tight_layout()

#from astropy.io import ascii  
#ascii.write(sources, '/Users/sarahdeveny/Research/GraduateThesis/candidates.txt') 

    
#with open('/Users/sarahdeveny/Research/GraduateThesis/candidates.txt', 'w') as f:
#    for item in sources:
#        f.write("%s\n" % item)
#    
#    
 #   histtype=u'step'
    
#####Get index and value of list    
#for idx, val in enumerate(F300):
#    print(idx, val)
    

################Uncomment this to plot the histograms of each y bin 
#######################bin 25-25.5
#bin25_fourandthree = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_neg4sig)&(BV<newstdev_move_neg3sig)]
#bin25_threeandtwo = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_neg3sig)s&(BV<newstdev_move_neg2sig)]
#bin25_twoandone = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_neg2sig)&(BV<newstdev_move_neg1sig)]
#bin25_oneandcenter = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_neg1sig)&(BV<newcenterline)]
#bin25_centerandone = BV[(F300>25.0)&(F300<25.5)&(BV>newcenterline)&(BV<newstdev_move_1sig)]
#bin25_oneandtwo = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_1sig)&(BV<newstdev_move_2sig)]
#bin25_twoandthree = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_2sig)&(BV<newstdev_move_3sig)]
#bin25_threeandfour = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)]
##bin25_fourand160 = BV[(F300>25.0)&(F300<25.5)&(BV>newstdev_move_4sig)&(BV<newtrackfunction1sig)]
#
#bin25 = np.array([len(bin25_fourandthree),len(bin25_threeandtwo),len(bin25_twoandone),len(bin25_oneandcenter),len(bin25_centerandone),len(bin25_oneandtwo),len(bin25_twoandthree),len(bin25_threeandfour)])
#
###plot hist
#plt.figure()
#bar_width = 1
#positions = np.arange(-4,4)
#pl4=plt.bar(positions,bin25,bar_width,alpha=0.8,edgecolor='k',align='edge')
#plt.xlabel('Bin 25 - 25.5')
#plt.ylim(0,200)
#plt.tight_layout()
#
#######################bin 23.5-24
#bin245_fourandthree = BV[(F300>24.5)&(F300<25.0)&(BV>newstdev_move_neg4sig)&(BV<newstdev_move_neg3sig)]
#bin245_threeandtwo = BV[(F300>24.5)&(F300<25.0)&(BV>newstdev_move_neg3sig)&(BV<newstdev_move_neg2sig)]
#bin245_twoandone = BV[(F300>24.5)&(F300<25.0)&(BV>newstdev_move_neg2sig)&(BV<newstdev_move_neg1sig)]
#bin245_oneandcenter = BV[(F300>24.5)&(F300<25.0)&(BV>newstdev_move_neg1sig)&(BV<newcenterline)]
#bin245_centerandone = BV[(F300>24.5)&(F300<25.0)&(BV>newcenterline)&(BV<newstdev_move_1sig)]
#bin245_oneandtwo = BV[(F300>24.5)&(F300<25.0)&(BV>newstdev_move_1sig)&(BV<newstdev_move_2sig)]
#bin245_twoandthree = BV[(F300>24.5)&(F300<25.0)&(BV>newstdev_move_2sig)&(BV<newstdev_move_3sig)]
#bin245_threeandfour = BV[(F300>24.5)&(F300<25.0)&(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)]
#
#bin245 = np.array([len(bin245_fourandthree),len(bin245_threeandtwo),len(bin245_twoandone),len(bin245_oneandcenter),len(bin245_centerandone),len(bin245_oneandtwo),len(bin245_twoandthree),len(bin245_threeandfour)])
#
###plot hist
#plt.figure()
#bar_width = 1
#positions = np.arange(-4,4)
#pl4=plt.bar(positions,bin245,bar_width,alpha=0.8,edgecolor='k',align='edge')
#plt.xlabel('Bin 24.5 - 25')
#plt.ylim(0,200)
#plt.tight_layout()
#
######################bin 24-24.5
#bin24_fourandthree = BV[(F300>24.0)&(F300<24.5)&(BV>newstdev_move_neg4sig)&(BV<newstdev_move_neg3sig)]
#bin24_threeandtwo = BV[(F300>24.0)&(F300<24.5)&(BV>newstdev_move_neg3sig)&(BV<newstdev_move_neg2sig)]
#bin24_twoandone = BV[(F300>24.0)&(F300<24.5)&(BV>newstdev_move_neg2sig)&(BV<newstdev_move_neg1sig)]
#bin24_oneandcenter = BV[(F300>24.0)&(F300<24.5)&(BV>newstdev_move_neg1sig)&(BV<newcenterline)]
#bin24_centerandone = BV[(F300>24.0)&(F300<24.5)&(BV>newcenterline)&(BV<newstdev_move_1sig)]
#bin24_oneandtwo = BV[(F300>24.0)&(F300<24.5)&(BV>newstdev_move_1sig)&(BV<newstdev_move_2sig)]
#bin24_twoandthree = BV[(F300>24.0)&(F300<24.5)&(BV>newstdev_move_2sig)&(BV<newstdev_move_3sig)]
#bin24_threeandfour = BV[(F300>24.0)&(F300<24.5)&(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)]
#
#bin24 = np.array([len(bin24_fourandthree),len(bin24_threeandtwo),len(bin24_twoandone),len(bin24_oneandcenter),len(bin24_centerandone),len(bin24_oneandtwo),len(bin24_twoandthree),len(bin24_threeandfour)])
#
###plot hist
#plt.figure()
#bar_width = 1
#positions = np.arange(-4,4)
#pl3=plt.bar(positions,bin24,bar_width,alpha=0.8,edgecolor='k',align='edge')
#plt.xlabel('Bin 24 - 24.5')
#plt.ylim(0,200)
#plt.tight_layout()
#
#######################bin 23.5-24
#bin235_fourandthree = BV[(F300>23.5)&(F300<24.0)&(BV>newstdev_move_neg4sig)&(BV<newstdev_move_neg3sig)]
#bin235_threeandtwo = BV[(F300>23.5)&(F300<24.0)&(BV>newstdev_move_neg3sig)&(BV<newstdev_move_neg2sig)]
#bin235_twoandone = BV[(F300>23.5)&(F300<24.0)&(BV>newstdev_move_neg2sig)&(BV<newstdev_move_neg1sig)]
#bin235_oneandcenter = BV[(F300>23.5)&(F300<24.0)&(BV>newstdev_move_neg1sig)&(BV<newcenterline)]
#bin235_centerandone = BV[(F300>23.5)&(F300<24.0)&(BV>newcenterline)&(BV<newstdev_move_1sig)]
#bin235_oneandtwo = BV[(F300>23.5)&(F300<24.0)&(BV>newstdev_move_1sig)&(BV<newstdev_move_2sig)]
#bin235_twoandthree = BV[(F300>23.5)&(F300<24.0)&(BV>newstdev_move_2sig)&(BV<newstdev_move_3sig)]
#bin235_threeandfour = BV[(F300>23.5)&(F300<24.0)&(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)]
#
#bin235 = np.array([len(bin235_fourandthree),len(bin235_threeandtwo),len(bin235_twoandone),len(bin235_oneandcenter),len(bin235_centerandone),len(bin235_oneandtwo),len(bin235_twoandthree),len(bin235_threeandfour)])
#
###plot hist
#plt.figure()
#bar_width = 1
#positions = np.arange(-4,4)
#pl2=plt.bar(positions,bin235,bar_width,alpha=0.8,edgecolor='k',align='edge')
#plt.xlabel('Bin 23.5 - 24')
#plt.ylim(0,200)
#plt.tight_layout()
#
######################bin 23-23.5
#bin23_fourandthree = BV[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_neg4sig)&(BV<newstdev_move_neg3sig)]
#bin23_threeandtwo = BV[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_neg3sig)&(BV<newstdev_move_neg2sig)]
#bin23_twoandone = BV[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_neg2sig)&(BV<newstdev_move_neg1sig)]
#bin23_oneandcenter = BV[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_neg1sig)&(BV<newcenterline)]
#bin23_centerandone = BV[(F300>23.0)&(F300<23.5)&(BV>newcenterline)&(BV<newstdev_move_1sig)]
#bin23_oneandtwo = BV[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_1sig)&(BV<newstdev_move_2sig)]
#bin23_twoandthree = BV[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_2sig)&(BV<newstdev_move_3sig)]
#bin23_threeandfour = BV[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)]
#
##plt.plot(bin23_threeandtwo,F300[(F300>23.0)&(F300<23.5)&(BV>newstdev_move_neg3sig)&(BV<newstdev_move_neg2sig)],'r*')
#
#bin23 = np.array([len(bin23_fourandthree),len(bin23_threeandtwo),len(bin23_twoandone),len(bin23_oneandcenter),len(bin23_centerandone),len(bin23_oneandtwo),len(bin23_twoandthree),len(bin23_threeandfour)])
##print bin23
#
###plot hist
#plt.figure()
#bar_width = 1
#positions = np.arange(-4,4)
#pl1=plt.bar(positions,bin23,bar_width,alpha=0.8,edgecolor='k',align='edge')
#plt.xlabel('Bin 23 - 23.5')
#plt.ylim(0,200)
#plt.tight_layout()
#
#######################bin 22.5-23
#bin225_fourandthree = BV[(F300>22.5)&(F300<23.0)&(BV>newstdev_move_neg4sig)&(BV<newstdev_move_neg3sig)]
#bin225_threeandtwo = BV[(F300>22.5)&(F300<23.0)&(BV>newstdev_move_neg3sig)&(BV<newstdev_move_neg2sig)]
#bin225_twoandone = BV[(F300>22.5)&(F300<23.0)&(BV>newstdev_move_neg2sig)&(BV<newstdev_move_neg1sig)]
#bin225_oneandcenter = BV[(F300>22.5)&(F300<23.0)&(BV>newstdev_move_neg1sig)&(BV<newcenterline)]
#bin225_centerandone = BV[(F300>22.5)&(F300<23.0)&(BV>newcenterline)&(BV<newstdev_move_1sig)]
#bin225_oneandtwo = BV[(F300>22.5)&(F300<23.0)&(BV>newstdev_move_1sig)&(BV<newstdev_move_2sig)]
#bin225_twoandthree = BV[(F300>22.5)&(F300<23.0)&(BV>newstdev_move_2sig)&(BV<newstdev_move_3sig)]
#bin225_threeandfour = BV[(F300>22.5)&(F300<23.0)&(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)]
#
#bin225 = np.array([len(bin225_fourandthree),len(bin225_threeandtwo),len(bin225_twoandone),len(bin225_oneandcenter),len(bin225_centerandone),len(bin225_oneandtwo),len(bin225_twoandthree),len(bin225_threeandfour)])
#
###plot hist
#plt.figure()
#bar_width = 1
#positions = np.arange(-4,4)
#pl5=plt.bar(positions,bin225,bar_width,alpha=0.8,edgecolor='k',align='edge')
#plt.xlabel('Bin 22.5 - 23')
#plt.ylim(0,200)
#plt.tight_layout()
#
#######################bin 22-23.5
#bin22_fourandthree = BV[(F300>22.0)&(F300<22.5)&(BV>newstdev_move_neg4sig)&(BV<newstdev_move_neg3sig)]
#bin22_threeandtwo = BV[(F300>22.0)&(F300<22.5)&(BV>newstdev_move_neg3sig)&(BV<newstdev_move_neg2sig)]
#bin22_twoandone = BV[(F300>22.0)&(F300<22.5)&(BV>newstdev_move_neg2sig)&(BV<newstdev_move_neg1sig)]
#bin22_oneandcenter = BV[(F300>22.0)&(F300<22.5)&(BV>newstdev_move_neg1sig)&(BV<newcenterline)]
#bin22_centerandone = BV[(F300>22.0)&(F300<22.5)&(BV>newcenterline)&(BV<newstdev_move_1sig)]
#bin22_oneandtwo = BV[(F300>22.0)&(F300<22.5)&(BV>newstdev_move_1sig)&(BV<newstdev_move_2sig)]
#bin22_twoandthree = BV[(F300>22.0)&(F300<22.5)&(BV>newstdev_move_2sig)&(BV<newstdev_move_3sig)]
#bin22_threeandfour = BV[(F300>22.0)&(F300<22.5)&(BV>newstdev_move_3sig)&(BV<newstdev_move_4sig)]
#
#bin22 = np.array([len(bin22_fourandthree),len(bin22_threeandtwo),len(bin22_twoandone),len(bin22_oneandcenter),len(bin22_centerandone),len(bin22_oneandtwo),len(bin22_twoandthree),len(bin22_threeandfour)])
#
###plot hist
#plt.figure()
#bar_width = 1
#positions = np.arange(-4,4)
#pl5=plt.bar(positions,bin22,bar_width,alpha=0.8,edgecolor='k',align='edge')
#plt.xlabel('Bin 22 - 22.5')
#plt.ylim(0,200)
#plt.tight_layout()

#plt.legend((pl4[0],pl3[0],pl2[0],pl1[0]), ('24.5-25', '24-24.5','23.5-24','23-23.5'))

##########################Put everything into a class###########################
##N=len(name)
##Star=np.zeros(N,dtype=object)
##class Stars:
##    def _int_(Self):
##        Self.name='ID'
##        Self.mag=0.0
##        Self.q=0.0
##        Self.per=0.0
##for i in range(0,N):
##    Star[i]=Stars()        
##def create(name, mag, q, per, Class):
##    for i in range(0,len(name)):
##        Star[i].name=name[i]
##        Star[i].mag=mag[i]
##        Star[i].q=q[i]
##        Star[i].per=per[i]
##    return Star
##        
##Star = create(Name,F300,Q300,percent,Stars)
##
###--------------------------Creating data file-----------------------------------
#
##datafile_path = "/Users/sarahdeveny/Research/GraduateThesis/Parameters/Data.txt"
##datafile_id = open(datafile_path, 'w+')
##
##data = np.array([Name300,F300,Q300,Percent300])#,F390,Q390,Percent390])
##data = data.T
##
##np.savetxt(datafile_id,data,fmt='%s',delimiter='\t',
##header='ID\t\t F300\t\t Q300\t Percent')
##
##datafile_id.close()
#
#
####################USING CLASSES##############################3
##N=len(name)
##Star=np.zeros(N,dtype=object)
##class Stars:
##    def _int_(Self):
##        Self.name='ID'
##        Self.mag=0.0
##        Self.q=0.0
##for i in range(0,N):
##    Star[i]=Stars()        
##def create(name, mag, q, Class):
##    for i in range(0,len(name)):
##        Star[i].name=name[i]
##        Star[i].mag=mag[i]
##        Star[i].q=q[i]
##    return Star
##        
##Star = create(name,F300X,q300,Stars)
##
##qqs = []
##for i in range(0,N):
##    qqs.append(Star[i].q)
##
##names = []
##for i in range(0,N):
##    names.append(Star[i].name)
##    
##mags = []
##for i in range(0,N):
##    mags.append(Star[i].mag)
##    
##sorting = np.array(qqs).argsort()
##magsorted = np.array(mags)[sorting[::]]
##qsorted = np.array(qqs)[sorting[::]]
##namesorted = name[sorting[::]]
##
##bins = np.linspace(5,-15,200)
##whichbin = np.digitize(magsorted,bins)
##
##mags = [magsorted[whichbin == i] for i in range(1,len(bins)+1)]
##qqs = [qsorted[whichbin == i] for i in range(1,len(bins)+1)]
##name = [namesorted[whichbin == i] for i in range(1,len(bins)+1)]
#
##Try this.... IDK about this -You can use the for loop with the bin to get the 
##same amount of stars in each bin.
#
##q300.sort()
##bins = []
##for i in range(0, len(q300), 400):
##    bin = q300[i: i+400]
##    bins.append(bin)
#
#################################CHECK FOR DUPLICATES###########################
#
##my_list.sort()
##for i in range(0,len(my_list)-1):
##               if my_list[i] == my_list[i+1]:
##                   print str(my_list[i]) + ' is a duplicate'
#    
#######GET INDICES##################
#
##indices = [i for i, x in enumerate(my_list) if x == "whatever"]
#    
