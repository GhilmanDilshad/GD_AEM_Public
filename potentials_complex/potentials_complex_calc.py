import pandas as pd
import numpy as np
from utilities import read_aquifer_xlsx, read_wells_xlsx


def __init__(self):
        pass


class Vectors_normal_well_complex:
    def __init__(self):  
        pass    

    def region_boundaries_normal_well_complex(self): # Region Coordinates, boundaries and meshgrids.
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


class Potentials_normal_well_complex:
    def potential_data_complex(self):

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
        # print('These are Zw values of Wells ',Zw)
        return baseFlowX[0], x_array, h0[0], H[0], k[0], y_array, pumping_array, Zw, alpha,Zref
    
    def phi_0_complex(self):  # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        _ , _ , h0, H, k, _, _, Zw, alpha,Zref = self.potential_data_complex()
        if h0 > H:          # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location
            phi_0 = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0 = 0.5 *k *h0 * h0

        # print('This is phi_0, at the reference head ', phi_0)
        return phi_0        


    def phi_well_plus_baseflow_complex(self):
        vectors_complex = Vectors_normal_well_complex()
        _, _, _, _, _, _, Z = vectors_complex.region_boundaries_normal_well_complex()
         
        baseFlowX , _ , h0, H, k, _, pumping_array, Zw, alpha, Zref = self.potential_data_complex()
        
        phi_well_complex = 0 
        
        
        # print('These are Z-Zw Values', Z-Zw)
        modulus_r = []

        # calculation for r 
        for i in range(len(Zw)):
            r = np.absolute(Z-Zw[i])
            modulus_r.append(r)
        # print('This is modulus of r ', r)
        #calculation for phi of well at some reference point
        for i in range(len(pumping_array)):
            omega_well = (pumping_array[i]/ (2*np.pi)) * (np.log((modulus_r[i])/(np.absolute(Zref - Zw[i]))))
            phi_well_complex += omega_well
        # print('These are r Values : ',r,)       
        # print('These are Omega_of_well_Complex: ', phi_well_complex)

        # Uniform Flow Field QoZeialpha
        phi_base_flow_complex = (-1 * baseFlowX) * (((Z.real)*(np.cos(alpha)))  + (Z.imag * (np.sin(alpha))))
        
       
        return phi_well_complex + phi_base_flow_complex
        

    def Stream_Fxn_complex(self):
        vectors_complex = Vectors_normal_well_complex()
        _, _, _, _, _, _, Z = vectors_complex.region_boundaries_normal_well_complex() 
        baseFlowX , _ , h0, H, k, _, pumping_array, Zw, alpha, Zref = self.potential_data_complex()

        psi_well_complex = 0
        
        #Stream_FXN_FOR_WELL_interms_of_complex
        for i in range(len(pumping_array)):
            Theta = np.angle(Z-Zw[i])
            omega_stream_fxn_well = (pumping_array[i] / (2*np.pi)) * (Theta)
            psi_well_complex += omega_stream_fxn_well


        #Stream_FXN_for_Base_Flow
        psi_base_flow_complex = (-1 * baseFlowX) * (((Z.imag) * (np.cos(alpha))) - ((Z.real) * (np.sin(alpha))))

        return psi_well_complex + psi_base_flow_complex
        
    def complex_potential_phi_total_region(self):

        phi_well_plus_base = self.phi_well_plus_baseflow_complex()
        phi_0 = self.phi_0_complex()
        phi_total_complex = phi_0 + phi_well_plus_base

        
        return phi_total_complex
    
    def complex_stream_total_region(self):
        psi_total_complex = self.Stream_Fxn_complex()
    
        return psi_total_complex

    def calculation_for_head_complex(self):
        h_list = [] # append all head values in this list

        phi_total_complex = self.complex_potential_phi_total_region()
        _ , _ , _, H, k, _, _, Zw, alpha, Zref = self.potential_data_complex()
      
        phicrit = 0.5 * k * H**2 


        for phi_item in phi_total_complex:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Book page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            
            h_list.append(h)
        
        head= np.array(h_list)
     
        # print('these are head values from complex with sepration', head)
        
        return head
    


