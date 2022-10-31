import pandas as pd
import numpy as np
import potentials_complex
from potentials_complex.potentials_complex_calc import Vectors_normal_well_complex
from utilities import read_aquifer_xlsx, read_wells_xlsx
import os

def __init__(self):
        pass

class Vectors_normal_well_complex_without_sepration:
    def __init__(self):  
        pass    

    def region_boundaries_normal_well_complex_without_sepration(self): # Region Coordinates, boundaries and meshgrids.
        Xmin = 0
        Xmax = 200
        Ymin = 0
        Ymax = 200
        xmesh = np.linspace(Xmin, Xmax, 201)
        ymesh = np.linspace(Ymin, Ymax, 201)
        [X_axis, Y_axis] = np.meshgrid(xmesh, ymesh)
        Z = (X_axis + Y_axis* 1j) 
        # print('These are complex numbers or complex coordinates',Z)
        # print('These Real Value of Z array',np.real(Z))
        # print('These are Imaginary values of Z array',np.imag(Z))

        return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z


class Potentials_normal_well_complex_without_sepration:
    def potential_data_complex_without_sepration(self):

        df1 = read_aquifer_xlsx()  # This is for reading values from aquifer file
        baseFlowX = df1['Base Flow in X direction']
        h0 = df1['Reference Head']
        H = df1['Aquifer Thickness'] # H is a variable storing value from H_thickness
        k = df1["Hydraulic Conductivity"] # k is a variable storing value from hydraulic conductivity 
        
        df2 = read_wells_xlsx()    # This is for reading values from well file
        x_coords = df2['x-coordinates'].tolist()
        y_coords = df2['y-coordinates'] .tolist()
        pumping = df2['pumping'].tolist()      
        x_array = np.array(x_coords)
        y_array = np.array(y_coords)
        pumping_array = np.array(pumping)
        Zw = (x_array + y_array * 1j)
        
        alpha = 0
        Zref = 0 + 0j
        
        return baseFlowX[0], x_array, h0[0], H[0], k[0], y_array, pumping_array, Zw, alpha, Zref
    
    def phi_0_complex_without_sepration(self):  # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        _ , _ , h0, H, k, _, _, Zw, alpha,Zref = self.potential_data_complex_without_sepration()
        if h0 > H:          # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location
            phi_0 = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0 = 0.5 *k *h0 * h0

        # print('This is phi_0, at the reference head ', phi_0)
        return phi_0        

    def phi_well_plus_baseflowphi_psi_well_plus_baseflowpsi_complex(self):
        vectors_complex = Vectors_normal_well_complex_without_sepration()
        _, _, _, _, _, _, Z = vectors_complex.region_boundaries_normal_well_complex_without_sepration()
         
        baseFlowX , _ , h0, H, k, _, pumping_array, Zw, alpha, Zref = self.potential_data_complex_without_sepration()
        phi_0_the_constant = self.phi_0_complex_without_sepration()

        modulus_r = []
        Zref_modulus = []
        Theta = []
        Omega_well = 0

        # Calculation of Reference Modulus
        for i in range(len(Zw)):
            ref_modulus = np.absolute(Zref - (Zw[i]))
            Zref_modulus.append(ref_modulus)

        # Calculating Angle/Argument at every mesh 
        for i in range(len(Zw)):
            theta = np.angle(Z-Zw[i])
            Theta.append(theta)

        # calculation for r 
        for i in range(len(Zw)):
            r = np.absolute(Z-(Zw[i]))
            modulus_r.append(r)
            
        #calculation for phi of well at some reference point
        for i in range(len(pumping_array)):
            Potential_Complex_Omega = ((pumping_array[i]/ (2*np.pi)) * (np.log((modulus_r[i])/(Zref_modulus[i])))) + (((pumping_array[i])* 1j * (Theta[i] )) / (2*np.pi))
            Omega_well += Potential_Complex_Omega
        
        
        phi_uniform_flow_field = ((-1* baseFlowX) * (Z * (np.exp(-1j * alpha))))

        Omega_plus_constant = Omega_well + phi_0_the_constant + phi_uniform_flow_field
        
        data1 = pd.DataFrame(Omega_plus_constant.real)
        path1 = os.path.join("potentials_complex", "phi_well_without_river.xlsx")
        data1.to_excel(path1, sheet_name='sheet1', index=False, )
        
        data2 = pd.DataFrame(Omega_plus_constant.imag)
        path2 = os.path.join("potentials_complex", "psi_well_without_river.xlsx")
        data2.to_excel(path2, sheet_name='sheet1', index=False, )
        
        return Omega_plus_constant

    def calculation_for_head_complex_without_sepration(self):
        h_list = [] # append all head values in this list
        Omega_plus_constant = self.phi_well_plus_baseflowphi_psi_well_plus_baseflowpsi_complex()
        _ , _ , _, H, k, _, _, Zw, alpha, Zref = self.potential_data_complex_without_sepration()
        phicrit = 0.5 * k * H**2 

        omega_real = Omega_plus_constant.real

        for phi_item in omega_real:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Book page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            
            h_list.append(h)
        
        head= np.array(h_list)
        data3 = pd.DataFrame(head.real)
        path3 = os.path.join("potentials_complex", "head_well_without_river.xlsx")
        data3.to_excel(path3, sheet_name='sheet1', index=False, )
        

        # print('baba chal jaaaaa')
        # print('these are head values from complex without sepration', head)
        
        return head
    


