#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 15:19:38 2023

@author: namnguyen
"""

import streamlit as st
import pandas as pd
import numpy as np
import math
import Confinement_Mander as CM



options=['Circular', 'Spiral', 'Rectangular']
selected_option=st.sidebar.selectbox('Select one option:', options)

if selected_option=='Circular':
    st.image('circular.jpeg')
    st.write('Parameters:')
    D=st.number_input('Diameter of transverse reinforcement D [m]:', 0.003 )
    d_s=st.number_input('Diameter of Circular/Spiral Column d_s [m]:', 0.14)
    s=st.number_input('spacing of Transverse reinforcements s [m]:', 0.1)

    f_yh=st.number_input('Yeilding stress of Transverse reinforcement steel f_yh [MPa]:', 345)
    
    f_co=st.number_input('f_co [MPa]', 34)
    A_long=st.number_input('Total Long rebar Area [$m^2$]', 0.01)
    
   #Create an Instant for Circular section:
    Circular= CM.Mander(f_co)
    
    #Calculate the Section properties: (saves internaly)
    Circular.Sect_Circular(D,s,d_s)
    
    #calculate the longitudinal reinforcement ratio: (saves internaly)
    A_long=st.number_input('A_long [$m^2$]:', 0.01)
    Circular.rho_cc(A_long)
    #calculate the effective lateral pressure(?) based on the section and reinforcement data:
    Circular.f_lat_eff()
    #Call the calculated data from the object to use in next calculation (for having the short form):
    f_effective=Circular.f_lat_eff
    sectType=Circular.Sect_type
    # calculates the maximum strength of confined concrete: (saves internaly)
    Circular.fcc(sectType,f_effective)
    # Call Mander.fc() to access the X, y data. (returns a panda dataframe)
    df=Circular.fc()
    st.write(df)
    
    # Plot
    st.line_chart(df, x='eps_c', y= 'f_c')
    # Set the x and y axis labels
    chart.set_axis_labels('Compressive Strain', 'Compressive Stress')
   

   
    
    
elif selected_option== 'Spiral':
    st.image('circular.jpeg')
    st.write('Parameters:')
    D=st.number_input('Diameter of transverse reinforcement (D):', 0.003 )
    d_s=st.number_input('Diameter of Circular/Spiral Column (d_s):', 0.14)
    s=st.number_input('spacing of Transverse reinforcements (s):', 0.1)

    f_yh=st.number_input('Yeilding stress of Transverse reinforcement steel (f_yh):', 345)
    
    f_co=st.number_input('f_co', 34)
    A_long=st.number_input('Total Long rebar Area', 0.01)
    
   #Create an Instant for Circular section:
    Spiral= CM.Mander(f_co)
    
    #Calculate the Section properties: (saves internaly)
    Spiral.Sect_Spiral(D,s,d_s)
    
    #calculate the longitudinal reinforcement ratio: (saves internaly)
    Spiral.rho_cc(A_long=0.01)
    #calculate the effective lateral pressure(?) based on the section and reinforcement data:
    Spiral.f_lat_eff()
    #Call the calculated data from the object to use in next calculation (for having the short form):
    f_effective=Spiral.f_lat_eff
    sectType=Spiral.Sect_type
    # calculates the maximum strength of confined concrete: (saves internaly)
    Spiral.fcc(sectType,f_effective)
    # Call Mander.fc() to access the X, y data. (returns a panda dataframe)
    df=Spiral.fc()
    st.write(df)
    # Plot
    st.line_chart(df, x='eps_c', y='f_c')
    
else:
    st.image('rectangular.jpeg')
    D=st.number_input('Diameter of transverse reinforcement (D):',3/1000)
   

    
    
    
    s=st.number_input('spacing of Transverse reinforcements (s):', 0.05)
    

    f_yh=st.number_input('Yeilding stress of Transverse reinforcement steel (f_yh):',345)
    
    #  (see the figure)
    #bc,dc=1,0.2
    bc=st.number_input('dimentions of Rect section (bc):',(320/1000)*0.9)
    dc=st.number_input('dimentions of Rect section (dc):',55/1000)
    # Total Area of Transverse reinforcement in x and y direction for Rectangle/Square section: (see the figure).
    f_co=st.number_input('f_co', 36)
    eps_max=st.number_input('eps_max', 0.03)
    
    st.write("Calculate A_sx, Asy :")
    A_sx=2*(np.pi*D**2)/4
    st.write('A_sx=',A_sx)
    A_sy=(2*(np.pi*D**2)/4)+(5*(np.pi*D**2)/4)
    st.write('A_sy=',A_sy)
    
    #A_sx=0 #set this zero for shell element model
    # distance between the Longitudinal rebars (define as a list: w=[w1,w2,w3] ):
    w=[0.05 for i in range(6)]
    

    # Creat an Instant for Rectangle or Square section:
    rect=CM.Mander(f_co,eps_max)
    
    #Calculate the Section properties: (saves internaly)
    rect.Sect_Rectangular(bc,dc,A_sx,A_sy,w,D,s,f_yh)
    
    
    #calculate the longitudinal reinforcement ratio: (saves internaly)
    # A_long: Total Long rebar Area
    #rect.rho_cc(A_long=0.032397674)
    st.write('Calculate A_long:')
    A_long=1*(np.pi*0.010**2)/4
    st.write('A_long=',A_long)
    rect.rho_cc(A_long=A_long)
    #calculate the effective lateral pressure(?) based on the section and reinforcement data:
    rect.f_lat_eff()
    
    #Call the calculated data from the object to use in next calculation (for having the short form):
    f_effective=rect.f_lat_eff
    sectType=rect.Sect_type
    
    # calculates the maximum strength of confined concrete: (saves internaly)
    rect.fcc(sectType,f_effective)
    # Call Mander.fc() to access the X, y data. (returns a panda dataframe)
    rect.fc()
    #print('f_l: ',rect.f_lat, 'lateral pressure from the transverse reinforcement')
    st.write('f_co: ',format(rect.f_co,'.3f'),'MPa')
    st.write('f_l_eff: ',rect.f_lat_eff, 'effective lateral pressure from the transverse reinforcement')
    st.write('E_c: ',format(rect.E_c,'.3f'),'MPa')
    st.write('E_sec: ',format(rect.E_sec,'.3f'),'MPa')
    
    df=rect.fc()
    
    st.write(df)
    # Plot
    st.line_chart(df, x='eps_c', y='f_c')
















