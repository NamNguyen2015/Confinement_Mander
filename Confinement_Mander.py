#!/usr/bin/env python
# coding: utf-8

# # Mander's Confinement (1988)

# This Code looks bad beacuse I tried to use classes for the first time.
# It works but needs a serious review. Also check the results.
# References: 'Mander_1988' and Mander 1983/4 PhD thesis.

# In[14]:


#Just scroll down to see the example
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import seaborn as sns
sns.set_style("whitegrid")
#class Concrete(object):
#    def __init__(self,name,fc,type):
#        self.name=name
#        self.fc=fc
#        self.type=type
        
        
#self.formulation='Mander_1988'        

class Mander(object):
    
    def __init__(self,f_co,f_cc=None,eps_co=0.002,eps_max=0.08):
        self.f_co=f_co
        self.f_cc=f_cc
        self.eps_co=eps_co
        self.eps_max=eps_max
        self.eps_c=np.arange(0,eps_max,0.0005)

    def fc(self,confined=True):
        f_co = self.f_co
        f_cc = self.f_cc
        eps_c = self.eps_c
        eps_co = self.eps_co
        eps_cc = eps_co*(1+5*((f_cc/f_co)-1))  # equation (5) in the paper
        self.eps_cc=eps_cc
        E_sec =(f_cc/eps_cc)  # equation (8)
        self.E_sec=E_sec
        E_c=5000*np.sqrt(f_co) # equation (7)
        self.E_c=E_c
        r=E_c/(E_c-E_sec)  # equation (6)
        x=eps_c/eps_cc  # equation (4)
        f_c=(f_cc*x*r)/(r-1+x**r)  # equation (3)
        df=pd.DataFrame({'$\epsilon_c$':eps_c,'$f_c$':f_c})
        return df
    
    
    def fcc(self,Sect_type,f_lat_eff):
        """
        f_lat_eff: the effective lateral confining stress on the concrete
        """
        f_co=self.f_co
        if Sect_type=='Spiral' or Sect_type=='Circular':
            self.f_cc=f_co*(2.254*np.sqrt(1+(7.94*f_lat_eff/f_co))-(2*f_lat_eff/f_co)-1.254)   # equtaion (29)
        elif Sect_type=='Rectangular' or Sect_type=='Square':
            #For details see Mander 1984 (PhD thesis) section 3.4.5
            S2=np.max(f_lat_eff)*-1
            S1=np.min(f_lat_eff)*-1
            
            #Iterate to find S3=f_cc
            # Initial guess:
            S3=S1
            i=0
            while i<100:
                # Estimate the confined strength S3 ($\sigma_3$) From this the octrahedral normal stress S_oct ($\sigma_{oct}$),
                # the octrahedral shear stress tau_oct ($\tau_{oct}), and the Lode angle theta ($\theta$), are determined:
                S_oct=(S1+S2+S3)/3
                tau_oct=np.sqrt(((S1-S2)**2+(S2-S3)**2+(S3-S1)**2))/3
                theta=np.arccos((S1-S_oct)/(np.sqrt(2)*tau_oct))
                S_oct_ratio=S_oct/f_co
                # these two sees to be polynominals based on Kupfer, Hilsdorf and Rusch ’ and Schickert and Winkler.
                T=0.069232-0.661091*S_oct_ratio-0.049350*S_oct_ratio**2
                C=0.122965-1.150502*S_oct_ratio-0.315545*S_oct_ratio**2
                D=4*(C**2-T**2)*(np.cos(theta)**2)
                tau_oct_bar=C*((0.5*D/np.cos(theta))+(2*T-C)*np.sqrt(D+(5*(T**2))-4*T*C))/(D+(2*T-C)**2)
                tau_oct=tau_oct_bar*f_co
                tau_oct # should be positive
                #print(tau_oct)
                S3_new=((S1+S2)/2)-np.sqrt((4.5*tau_oct**2)-(0.75*(S1-S2)**2))
                error=np.abs((S3_new-S3)/S3_new)
                if error<=0.001:
                    print(' Iteration converged with ',format(error*100,'.2f'),'% relative error')
                    print('S3: ',format(S3,'.3f'),',next estimate: ',format(S3_new,'.3f'))
                    f_cc=np.abs(S3)
                    print('f_cc: ',format(f_cc,'.3f'),'MPa')
                    self.f_cc=f_cc
                    break
                else:
                    S3=S3_new
                    i += 1
            else:
                print('No convergance reached!!!')
                print('Please calculate use the figure 4. of Manders 1988 peper to calculate the confined strenght ratio (fcc/fco) manualy')

                

    def Sect_Spiral(self,D,s,d_s,f_yh=345):
        self.Sect_type='Spiral'
        self.f_yh=f_yh
 
        """
        fyh = yield strength of the transverse reinforcement
        A_eff: area of effectively confined concrete core
        d_s: diameter of spiral between bar centers
        Asp = area of transverse reinforcement bar
        s = center to center spacing or pitch of spiral or circular hoop
        """
        
        #A_core
        self.A_core=(np.pi/4)*(d_s**2) # A_c in the paper
        
        # A_eff:  area of effectively confined concrete core
        s1=s-D # s1=s'
        A_sp=np.pi*(D**2)/4
        self.rho_s=(4*A_sp)/(d_s*s)  # equation (17)
        self.A_eff=(np.pi/4)*(d_s**2)*(1-s1/(2*d_s))  #equation (12)
        
        # f_lat  = lateral pressure from the transverse reinforcement
        self.f_lat=0.5*self.rho_s*f_yh # equation (18)

            
    def Sect_Circular(self,D,s,d_s,f_yh=345):
        #Similar to Spiral, Only update the A_eff:
        self.Sect_type='Circular'
        self.Sect_Spiral(D,s,d_s,f_yh)
        #self.Sect_type='Circular'
        s1=s-D
        self.A_eff=(np.pi/4)*(d_s**2)*(1-s1/(2*d_s))**2  # equation (12) where s1= s'
        
    def Sect_Rectangular(self,bc,dc,A_sx,A_sy,w,D,s,f_yh=345):
        """
        Asx and Asy = the total area of transverse bars running in the x and y directions, respectively
        """
        
        self.Sect_type='Rectangular'
        self.f_yh=f_yh
        self.rho_s=np.array([A_sx/(s*dc),A_sy/(s*bc)]) # comes from rho_x and rho_y equations (23) & (24)
        self.A_ineff=np.sum((np.array(w)**2)/6) # equation (20)
        self.A_core=bc*dc
        s1=s-D
        self.A_eff=(self.A_core-self.A_ineff)*(1-(s1/(2*bc)))*(1-(s1/(2*dc))) # equation (21)
        self.f_lat=self.rho_s*self.f_yh # equation (25) and (26)

    def rho_cc(self,A_long):
        """
        A_long: total Area of longitudinal reinforcement"""
        self.A_long=A_long
        self.rho_cc=A_long/self.A_core # equation (11) rho_cc= A_long/ A_c
    
    def f_lat_eff(self):
        A_cc=self.A_core*(1-self.rho_cc) # equation (11)
        self.A_cc=A_cc
        k_eff=self.A_eff/A_cc # equation (10)
        self.k_eff=k_eff
        self.f_lat_eff=self.f_lat*k_eff # equation (9)





# 

# ## UNIFIED STRESS-STRAIN APPROACH FOR CONFINED CONCRETE WITH MONOTONIC LOADING AT SLOW STRAIN RATES
# ![image.png](attachment:image.png)

# 

# ![image.png](attachment:image.png)

# In[ ]:





# ### Effective Lateral Confining Pressure and the Confinement Effectiveness Coefficient
# The maximum transverse pressure from the confining 
# steel can only be exerted effectively on that part of the concrete core where 
# the confining stress has fully developed due to arching action. Figs. 2 and 
# 3 show the arching action that is assumed to occur between the levels of 
# transverse circular and rectangular hoop reinforcement. Midway between 
# the levels of the transverse reinforcement, the area of ineffectively 
# confined concrete will be largest and the area of effectively confined 
# concrete core $A_e$ will be smallest.
# 
# <div>
# <img src="attachment:image.png" align="left" width="350"/>
# <img src="attachment:image-2.png" align="right" width="350"/>
# </div>

# $f_{l}=f_{lat} k_e$
# 
# $f_{lat}$ : lateral pressure from the transverse reinforcement, assumed to 
# be uniformly distributed over the surface of the concrete core

# $k_e = \frac{A_e}{A_{cc}}$ : confinement effectiveness coefficient;
# 
# $A_e$ : area of effectively confined concrete core;
# 
# $A_{cc}=A_c (1-\rho_{cc})$
# 
# $\rho_{cc}$ : ratio of area of longitudinal reinforcement to area of core of section;
# 
# $A_c$ : area of core of section enclosed by the center lines of the perimeter spiral or hoop

# ### Confinement Effectiveness for Sections Confined by Spirals or Circular Hoops

# In[23]:


# Diameter of transverse reinforcement (D):
D=0.003

# Diameter of Circular/Spiral Column (d_s): (diameter of the loop/hoop see the figure)
d_s=0.14

# spacing of Transverse reinforcements (s):
s=0.1

# Yeilding stress of Transverse reinforcement steel (f_yh):
f_yh=345


# Creat an Instant for Circular section:
Circular=Mander(34)

#Calculate the Section properties: (saves internaly)
Circular.Sect_Circular(D,s,d_s)

# (for spiral use Mander.Sect_Circular(D,s,d_s))

#calculate the longitudinal reinforcement ratio: (saves internaly)
# A_long: Total Long rebar Area
Circular.rho_cc(A_long=0.01)

#calculate the effective lateral pressure(?) based on the section and reinforcement data:
Circular.f_lat_eff()

#Call the calculated data from the object to use in next calculation (for having the short form):
f_effective=Circular.f_lat_eff
sectType=Circular.Sect_type

# calculates the maximum strength of confined concrete: (saves internaly)
Circular.fcc(sectType,f_effective)

# Plot
sns.lineplot(data=Circular.fc(),x='$\epsilon_c$',y='$f_c$')

# Call Mander.fc() to access the X, y data. (returns a panda dataframe)
Circular.fc()


# ### Confinement Effectiveness for Sections Confined by Rectangular or Square sections

# In[21]:


# Diameter of transverse reinforcement (D):
#D=0.015
D=3/1000


# spacing of Transverse reinforcements (s):
s=0.05

# Yeilding stress of Transverse reinforcement steel (f_yh):
f_yh=345

# dimentions of Rect section: (see the figure)
#bc,dc=1,0.2
bc=(320/1000)*0.9
dc=55/1000
# Total Area of Transverse reinforcement in x and y direction for Rectangle/Square section: (see the figure).
A_sx=2*(np.pi*D**2)/4
A_sy=(2*(np.pi*D**2)/4)+(5*(np.pi*D**2)/4)
#A_sx=0 #set this zero for shell element model
# distance between the Longitudinal rebars (define as a list: w=[w1,w2,w3] ):
w=[0.05 for i in range(6)]

# Creat an Instant for Rectangle or Square section:
rect=Mander(36,eps_max=0.03)

#Calculate the Section properties: (saves internaly)
rect.Sect_Rectangular(bc,dc,A_sx,A_sy,w,D,s,f_yh)


#calculate the longitudinal reinforcement ratio: (saves internaly)
# A_long: Total Long rebar Area
#rect.rho_cc(A_long=0.032397674)
A_long=1*(np.pi*0.010**2)/4
rect.rho_cc(A_long=A_long)
#calculate the effective lateral pressure(?) based on the section and reinforcement data:
rect.f_lat_eff()

#Call the calculated data from the object to use in next calculation (for having the short form):
f_effective=rect.f_lat_eff
sectType=rect.Sect_type

# calculates the maximum strength of confined concrete: (saves internaly)
rect.fcc(sectType,f_effective)

# Plot
sns.lineplot(data=rect.fc(),x='$epsilon_c$',y='$f_c$')

# Call Mander.fc() to access the X, y data. (returns a panda dataframe)
rect.fc()
#print('f_l: ',rect.f_lat, 'lateral pressure from the transverse reinforcement')
print('f_co: ',format(rect.f_co,'.3f'),'MPa')
print('f_l_eff: ',rect.f_lat_eff, 'effective lateral pressure from the transverse reinforcement')
print('E_c: ',format(rect.E_c,'.3f'),'MPa')
print('E_sec: ',format(rect.E_sec,'.3f'),'MPa')


# In[6]:


rect.eps_cc


# In[ ]:





# ### Abaqus CDP format

# In[7]:


unconfined=Mander(68,68,eps_max=0.03)
model=unconfined
#model=rect

data=model.fc()
E_c=model.E_c
f_cc=model.f_cc
eps_c=data['eps_c'].to_numpy()
f_c=data['f_c'].to_numpy()

#E=f_c/eps_c
eps_yeild=0.4*f_cc/E_c
f_yeild=0.4*f_cc
## linearize the elastic part

data.loc[data['eps_c']<eps_yeild,'f_c']=data['eps_c']*E_c

#eps_plastic=eps_c-eps_yeild


#data['f_c'].max()
#data['Damage_c_lublinear'] = 0
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-(data['f_c']/f_cc)

#data['Damage_c'] = 0
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-(data['f_c']/f_cc)




#data['f_c'].max()
#data['Damage_T_lublinear'] = 0

#damage formulation using the paper birtel 2006 (also add it for compression) (NOT USED FOR COMPRESSION)
b_t= 0.1 # (between 0 and 1)( 0.1 for tension and 0.7 for compression)
b_c=0.7
#eps_yeild

#assumin element length l_0=1 as indicated in abaqus manual)

data['eps_inelastic']=data['eps_c']-(data['f_c']/E_c)
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'eps_inelastic']=data['eps_c']-(data['f_c']/E_c)
#data.loc[data['eps_c']<eps_yeild,'eps_inelastic']=0
#data['eps_pl_c']=data['eps_inelastic']*b_c
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-((data['f_c']/E_c)/((data['eps_pl_c']*((1/b_c)-1))+(data['f_c']/E_c)))
#data['Damage_c']=1-((data['f_c']/E_c)/((data['eps_pl_c']*((1/b_c)-1))+(data['f_c']/E_c)))

eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-(data['f_c']/f_cc)
data['Damage_c']=data['Damage_c'].fillna(0)
data['eps_inelastic']=data['eps_inelastic'].fillna(0)
data=data.round(decimals=15)

## convert units for abaqus
data['f_c']=data['f_c']*1e6



# In[8]:


#data=rect.fc()
data.head()


# In[9]:


sns.lineplot(data=data,x='eps_c',y='Damage_c')


# In[10]:


data.to_excel("Takahashi2000_confined.xlsx",sheet_name="Compression")


# In[11]:


rect.E_c


# In[12]:


eps_f_cc


# In[13]:


test=Mander(68,68,eps_max=0.03)


# In[14]:


test.fc()


# In[15]:


sns.lineplot(data=test.fc(),x='eps_c',y='f_c')


# In[16]:


test.E_c


# In[17]:


bc


# In[18]:


data=rect.fc()


# In[19]:


data


# In[ ]:





# # For Silvia

# In[25]:


# Diameter of transverse reinforcement (D):
#D=0.015
D=16/1000


# spacing of Transverse reinforcements (s):
s=0.075

# Yeilding stress of Transverse reinforcement steel (f_yh):
f_yh=420

# dimentions of Rect section: (see the figure)
#bc,dc=1,0.2
bc=2.9
dc=2.4
# Total Area of Transverse reinforcement in x and y direction for Rectangle/Square section: (see the figure).
A_sx=10*(np.pi*D**2)/4
A_sy=(11*(np.pi*D**2)/4)+(5*(np.pi*D**2)/4)
#A_sx=0 #set this zero for shell element model
# distance between the Longitudinal rebars (define as a list: w=[w1,w2,w3] ):
w=[0.15 for i in range(60)]

# Creat an Instant for Rectangle or Square section:
Silvia=Mander(40,eps_max=0.05)

#Calculate the Section properties: (saves internaly)
Silvia.Sect_Rectangular(bc,dc,A_sx,A_sy,w,D,s,f_yh)


#calculate the longitudinal reinforcement ratio: (saves internaly)
# A_long: Total Long rebar Area
#rect.rho_cc(A_long=0.032397674)
A_long=120*(np.pi*0.040**2)/4
Silvia.rho_cc(A_long=A_long)
#calculate the effective lateral pressure(?) based on the section and reinforcement data:
Silvia.f_lat_eff()

#Call the calculated data from the object to use in next calculation (for having the short form):
f_effective=Silvia.f_lat_eff
sectType=Silvia.Sect_type

# calculates the maximum strength of confined concrete: (saves internaly)
Silvia.fcc(sectType,f_effective)

# Plot
sns.lineplot(data=Silvia.fc(confined=False),x='eps_c',y='f_c')

# Call Mander.fc() to access the X, y data. (returns a panda dataframe)
Silvia.fc(confined=False)
#print('f_l: ',rect.f_lat, 'lateral pressure from the transverse reinforcement')
print('f_co: ',format(Silvia.f_co,'.3f'),'MPa')
print('f_l_eff: ',Silvia.f_lat_eff, 'effective lateral pressure from the transverse reinforcement')
print('E_c: ',format(Silvia.E_c,'.3f'),'MPa')
print('E_sec: ',format(Silvia.E_sec,'.3f'),'MPa')


# In[19]:


Silvia.fc()


# ## For Harbour

# In[20]:


# Unconfined
# Creat an Instant for Rectangle or Square section:
HB=Mander(68,68,eps_max=0.01)


# Plot
sns.lineplot(data=HB.fc(confined=False),x='eps_c',y='f_c')

# Call Mander.fc() to access the X, y data. (returns a panda dataframe)
HB.fc(confined=False)
#print('f_l: ',rect.f_lat, 'lateral pressure from the transverse reinforcement')
print('f_co: ',format(HB.f_co,'.3f'),'MPa')
print('f_l_eff: ',HB.f_lat_eff, 'effective lateral pressure from the transverse reinforcement')
print('E_c: ',format(HB.E_c,'.3f'),'MPa')
print('E_sec: ',format(HB.E_sec,'.3f'),'MPa')


# In[21]:


HB.fc(confined=False)


# In[22]:


model=HB
#model=rect

data=model.fc()
E_c=model.E_c
f_cc=model.f_cc
eps_c=data['eps_c'].to_numpy()
f_c=data['f_c'].to_numpy()

#E=f_c/eps_c
eps_yeild=0.4*f_cc/E_c
f_yeild=0.4*f_cc
## linearize the elastic part

data.loc[data['eps_c']<eps_yeild,'f_c']=data['eps_c']*E_c

#eps_plastic=eps_c-eps_yeild


#data['f_c'].max()
#data['Damage_c_lublinear'] = 0
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-(data['f_c']/f_cc)

#data['Damage_c'] = 0
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-(data['f_c']/f_cc)




#data['f_c'].max()
#data['Damage_T_lublinear'] = 0

#damage formulation using the paper birtel 2006 (also add it for compression) (NOT USED FOR COMPRESSION)
b_t= 0.1 # (between 0 and 1)( 0.1 for tension and 0.7 for compression)
b_c=0.7
#eps_yeild

#assumin element length l_0=1 as indicated in abaqus manual)

data['eps_inelastic']=data['eps_c']-(data['f_c']/E_c)
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'eps_inelastic']=data['eps_c']-(data['f_c']/E_c)
#data.loc[data['eps_c']<eps_yeild,'eps_inelastic']=0
#data['eps_pl_c']=data['eps_inelastic']*b_c
#eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
#data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-((data['f_c']/E_c)/((data['eps_pl_c']*((1/b_c)-1))+(data['f_c']/E_c)))
#data['Damage_c']=1-((data['f_c']/E_c)/((data['eps_pl_c']*((1/b_c)-1))+(data['f_c']/E_c)))

eps_f_cc=data.loc[data['f_c']==data['f_c'].max(),'eps_c'].values[0]
data.loc[data['eps_c']>eps_f_cc,'Damage_c']=1-(data['f_c']/f_cc)
data['Damage_c']=data['Damage_c'].fillna(0)
data['eps_inelastic']=data['eps_inelastic'].fillna(0)
data=data.round(decimals=15)

## convert units for abaqus
data['f_c']=data['f_c']*1e6



# In[23]:


data.to_excel("HB_Unconfined.xlsx",sheet_name="Compression")


# In[24]:


data


# In[25]:


E_c


# In[ ]:




