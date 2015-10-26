# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt 
import math as mat

"""
Created on Sun Aug  9 14:04:59 2015
@author: Guillaume Boudoire (LGSR/OVPF/IPGP)

# TO DO LIST
    1) Créer un vrai fichier pour comparer
    2) Vérifier le calcul des coefficients de diffusion
    3) Vérifier le passage de fO2 à NNO
    4) Description of input file and file ReadMe (abscisses 
    dans ordre croissant et bonnes unités notamment D m2/s et Fo 0.8)
    5) Calculation of the error (coeff & deviation)
    6) Calculation fitting function
    7) Multi-plateaux
    7) Interface Tkinter
    8) Publication "Estimating residence time through crystal zoning analysis:
       #simple diffusion VolcTimePy program in Python 2.7 (Boudoire,
       Boissier, Simon, DiMuro, Michon, Gillet)
"""

###################################################################################################
##################################### VolcTimePy1.0 programm: #####################################
################ TIMESCALES CALCULATION FROM DIFFUSION PROFILES IN ZONED CRYSTALS #################
###################################################################################################

############################## INPUT (CAN BE MODIFY) and INITIALISATION ###########################

Time = 0 # Residence time
Deltat =  3600 # Time pitch (s) CAN BE MODIFIED if the stability test is bad
res_a = 0.3 # Initial plot resolution (µm) CAN BE MODIFIED if the stability test is bad
res_a2 = float(res_a)*0.000001 # Interval between each measurments (meters)
mineral = raw_input("What mineral do you want to study? (olivine/clinopyroxene/plagioclase):")
kind = raw_input("What does your diffusion profile look like? (crenel/echelon):")
path = raw_input("Measurements file directory (.csv):") # File import from Desktop (first column with abscissas and second with values)
coef = raw_input("Coefficients file directory (.csv):") # File import from Desktop (first column with parameters names and second with values)
f = open(path)
f_param = open(coef)
plt.ion()

###################################### FUNCTIONS DEFINITION #######################################

### Routine functions ###


def getIndice(dist, x): # Obtain the point index in a list
    i = -1
    for k in dist:
        if k<= x:
            i = i+1
        else:
            return i 
            
                       
def drange(start, stop, step): # Range (for float values) for running the whole interval in ascending order
    r = start
    while r < stop:
        yield r
        r += step    
   

### Building functions ###


def data_function(files): # Remplacer par f
    tlx = 0 # Profil lenght
    dist = [] #  Measurements distance (X)
    val = []  #  Measurements values for the adopted element (Y)
    for line in files: # Create lists from columns file
        tlx = tlx+1 
        line = line[0:len(line)-1]
        re = line.split(";") # Column separator in the input file
        dist.append(float(re[0]))
        val.append(float(re[1]))    
    return [dist, val, tlx]
     

def initial_crenel(tl,fom,va,dis,resa): # Define the initial function on the appropriate interval whose lenght = tlx with (Normal/Inverse)
    ptsmin = 0
    ptsmax = 0
    xmin = 0
    xmax = 0
    ymin = input("Minimum plateau value:") #min(va) 
    ymax = input("Maximum plateau value:") #max(va) 
    manual = raw_input("Have you enough points to define the crenel interfaces position? (Y/N):")    
    tl = dis[tl-1]
    tly = dis[0]
    a = []
    dist_a = []  
    save_a = [] #  Initial window diffusion function 
    if fom == "normal":
        if manual=="Y":
            while xmin==0:
                if float(va[ptsmin]) <= (ymax+ymin)/2 and float(va[ptsmin+1]) >= (ymax+ymin)/2:
                    xmin=float((float(dis[ptsmin+1])+float(dis[ptsmin]))/2) # Instead of input ("Abscissa of the first interface:")
                    ptsmax=ptsmin+1
                else:
                    ptsmin=ptsmin+1
            while xmax==0:
                if float(va[ptsmax]) >= (ymax+ymin)/2 and float(va[ptsmax+1]) <= (ymax+ymin)/2:
                    xmax=float((float(dis[ptsmax+1])+float(dis[ptsmax]))/2) # Instead of input ("Value of the maximum plateau for the element:")
                else:
                    ptsmax=ptsmax+1
        else:
            xmin = input ("First interface position:")
            xmax = input ("Last interface position:")
        for x in drange(int(tly), int(tl), resa): # Need tlx+1 sometimes
            dist_a.append(x);
            if(xmin<=x and x<=xmax):
                a.append(ymax)
            else:
                a.append(ymin)      
    else:
        if manual=="Y":
            while xmin==0:
                if float(va[ptsmin]) >= (ymax+ymin)/2 and float(va[ptsmin+1]) <= (ymax+ymin)/2:
                    xmin=float((float(dis[ptsmin+1])+float(dis[ptsmin]))/2) # Instead of input ("Abscissa of the first interface:")
                    ptsmax=ptsmin+1
                else:
                    ptsmin=ptsmin+1
            while xmax==0:
                if float(va[ptsmax]) <= (ymax+ymin)/2 and float(va[ptsmax+1]) >= (ymax+ymin)/2:
                    xmax=float((float(dis[ptsmax+1])+float(dis[ptsmax]))/2) # Instead of input ("Value of the maximum plateau for the element:")
                else:
                    ptsmax=ptsmax+1
        else:
            xmin = input ("First interface position:")
            xmax = input ("Last interface position:")            
        for x in drange(int(tly), int(tl), resa): # Need tlx+1 sometimes
            dist_a.append(x);
            if(xmin>=x and x>=xmax):
                a.append(ymax)
            else:
                a.append(ymin)      
    for h in a:
        save_a.append(h) 
    return [dist_a, save_a]
   
    
def initial_echelon(tl,form,va,dis,resa): # Define the initial function on the appropriate interval whose lenght = tlx with (Right/Left)
    ptsmin = 0
    xlim = 0
    ymin = input("Minimum plateau value:") #min(va) 
    ymax = input("Maximum plateau value:") #max(va) 
    manual = raw_input("Have you enough points to define the crenel interfaces position? (Y/N):")    
    tl = dis[tl-1]
    tly = dis[0]
    a = []
    dist_a = []  
    save_a = [] #  Initial window diffusion function 
    if form == "right":
        if manual == "Y":
            while xlim==0:
                if float(va[ptsmin]) <= (ymax+ymin)/2 and float(va[ptsmin+1]) >= (ymax+ymin)/2:
                    xlim=float((float(dis[ptsmin+1])+float(dis[ptsmin]))/2) # Instead of input ("Abscissa of the first interface:")
                else:
                    ptsmin=ptsmin+1
        else:
            xlim = input("Interface position:")
        for x in drange(int(tly), int(tl), resa): # Need tlx+1 sometimes
            dist_a.append(x);
            if(xlim<=x):
                a.append(ymax)
            else:
                a.append(ymin)
    else:
        if manual == "Y":
            while xlim==0:
                if float(va[ptsmin]) >= (ymax+ymin)/2 and float(va[ptsmin+1]) <= (ymax+ymin)/2:
                    xlim=float((float(dis[ptsmin+1])+float(dis[ptsmin]))/2) # Instead of input ("Abscissa of the first interface:")
                else:
                    ptsmin=ptsmin+1
        else:
            xlim = input("Interface position:")
        for x in drange(int(tly), int(tl), resa): # Need tlx+1 sometimes
            dist_a.append(x);
            if(xlim>=x):
                a.append(ymax)
            else:
                a.append(ymin)          
    for h in a:
        save_a.append(h) 
    return [dist_a, save_a]
    
### Diffusion coefficient calculation functions ###

R = 8.3144621 # Gas constant (J/K*mol)

def diff_olivine(Forsterite): # Give a diffusion coefficient from the parameters file for olivine from formula of Girona and al. (2013)
    valeurs = []
    for line in f_param: # Create lists from columns file
        line = line[0:len(line)]
        re = line.split(";") # Column separator in the input file
        valeurs.append(float(re[1]))
    Angle_a = float(valeurs[0])
    Angle_b = float(valeurs[1])
    Angle_c = float(valeurs[2])
    Pressure = float(valeurs[3])
    Temperature = float(valeurs[4])
    NNO = float(valeurs[5])
    fO2_NNO=101325*(NNO) # Pa 
    Dc = (10**(-9.21))*((fO2_NNO/((10**(-7))))**(1/6))*10**(3*(0.9-Forsterite)) \
    *mat.exp((-201+(7*10**(-6))*(Pressure-(10**(5))))/(R*Temperature)) 
    Da = Dc/6
    Db = Da
    D = (Da*(mat.cos(Angle_a))**2)+(Db*(mat.cos(Angle_b))**2)+(Dc*(mat.cos(Angle_c))**2)
    return D
        
        
def diff_others(): # Give a diffusion coefficient from the parameters file for other mineral with Arrhenius' law
    valeurs = []
    for line in f_param: # Create lists from columns file
        line = line[0:len(line)]
        re = line.split(";") # Column separator in the input file
        valeurs.append(float(re[1]))
    Da_init = float(valeurs[0])
    Db_init = float(valeurs[1])
    Dc_init = float(valeurs[2])
    Q = float(valeurs[3])
    Angle_a = float(valeurs[4])
    Angle_b = float(valeurs[5])
    Angle_c = float(valeurs[6])
    Temperature = float(valeurs[7])
    Da = Da_init*mat.exp(-Q/(R*Temperature))
    Db = Db_init*mat.exp(-Q/(R*Temperature))
    Dc = Dc_init*mat.exp(-Q/(R*Temperature))
    D = (Da*(mat.cos(Angle_a))**2)+(Db*(mat.cos(Angle_b))**2)+(Dc*(mat.cos(Angle_c))**2)
    return D
    
############################################ MAIN #################################################

        # Measurements function

[dist,val,tlx]=data_function(f)
# Abscissas list of each measurement points / Values list of each measurement points / Number of analysis

        # Initial function

if kind == "crenel":
    form = raw_input("What is the crenel form? (normal/inverse):")
    [da,sa] = initial_crenel(tlx,form,val,dist,res_a)
    dist_a = [da,sa][0] 
    save_a = [da,sa][1]
# Abscissas list of each crenel function points and Values list of each crenel function points
else:
    form = raw_input("What is the echelon form? (right/left):")
    [da,sa] = initial_echelon(tlx,form,val,dist,res_a)
    dist_a = [da,sa][0] 
    save_a = [da,sa][1]
# Abscissas list of each echelon function points and Values list of each echelon function points

        # Diffusion coefficients & Stability test
    
if mineral == "olivine":
    Testmin=0
    Testmax=0
    Dol = []
    minFo = min(val)
    maxFo = max(val)
    Testmin = float(diff_olivine(minFo))*float(Deltat)/(float(res_a2)*float(res_a2))
    Testmax = float(diff_olivine(maxFo))*float(Deltat)/(float(res_a2)*float(res_a2))
    if max(Testmin,Testmax)<0.5:
        print "Stability test fulfilled (<0.5):", max(Testmin,Testmax)
    else:
        print "No stability: please change the time pitch (Deltat) or the resolution (res_a)"
    for Dfo in val:
        Dol.append(diff_olivine(Dfo))
else:
    Test=0
    D = diff_others()
    Test = float(D)*float(Deltat)/(float(res_a2)*float(res_a2));
    if Test<0.5:
        print "Stability test fulfilled (<0.5):", Test
    else:
        print "No stability: please change the time pitch (Deltat) or the resolution (res_a)"
        # Iteration model based on Fick's Law (finite elements)

timef = (float(Deltat)/((res_a*(10^(-6)))*(res_a*(10^(-6))))) # Time parameter in the equation 
impress = "N"
while impress != "Y":      
    timestep = input ("Iteration number:") 
    timestep = int(timestep)
    a = [] # Final function
    for h in save_a:
        a.append(h)    
    for k in range(timestep):   
        b = [] # Function (step n-1) used to define the final function (step n)
        for h in a:
            b.append(h)       
        for i in range(1, len(a)-1):
            if mineral == "olivine":
                a[i] = b[i]+timef*((Dol[i+1]-Dol[i])*(b[i+1]-b[i])+Dol[i]*(b[i+1]-2*b[i]+b[i-1])) 
            else:
                a[i] = b[i]+D*timef*(b[i+1]-2*b[i]+b[i-1]) 
    plt.plot(dist_a, a, 'go')
    plt.plot(dist, val, 'bx')
    plt.plot(dist_a, save_a, 'rx')
    plt.show(block=True)
    impress = raw_input ("Are you agree with this iteration (Y/N):")
    
        # Final corresponding residence time
    
Time = (timestep*Deltat)/(3600*24)
print "Minimum residence time corresponds to:" , str(Time), "days"




    









