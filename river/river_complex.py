import pandas as pd
import numpy as np
from sympy import Q
from utilities import read_aquifer_xlsx, read_wells_xlsx
import os

#Defining Vector Class Again Due to Seprate or Bigger Regio Required Changing the vector Parameters from the vector file. So a new Vector Class for river
class Vectors_river_complex:
    def __init__(self):  
        pass    

    def region_boundaries_river_complex(self): # Region Coordinates, boundaries and meshgrids.
        Xmin = -200
        Xmax = 200
        Ymin = 0
        Ymax = 200
        xmesh = np.linspace(Xmin, Xmax, 201)
        ymesh = np.linspace(Ymin, Ymax, 201)
        [X_axis, Y_axis] = np.meshgrid(xmesh, ymesh)
        Z = (X_axis + Y_axis*1j)
        # print('these are z values complex', Z)
        return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z

    def vectors_data_river_complex(self):    # getting coordinates of well(s) and converting them to list 
        df3 = read_wells_xlsx()   # This is for reading values from well file
       
        x_values = df3['x-coordinates'].tolist()   
        y_values = df3['y-coordinates'].tolist()

        x_array = np.array(x_values) # converting list to array
        y_array = np.array(y_values) # converting list to array
        return x_values, y_values, x_array, y_array


class Potentials_River_complex:
    def __init__(self):
        pass

    def potentials_data_river_complex(self): #df1 is data frame for the aquifer
        
        df1 = read_aquifer_xlsx()  # This is for reading values from aquifer file
        baseFlowX = df1['Base Flow in X direction']
        h0 = df1['Reference Head']
        H = df1['Aquifer Thickness'] # H is a variable storing value from H_thickness
        k = df1["Hydraulic Conductivity"] # k is a variable storing value from hydraulic conductivity 
        por = df1["Porosity"]
    
        df2 = read_wells_xlsx()    # This is for reading values from well file
        x_coords = df2['x-coordinates']
        y_coords = df2['y-coordinates'] 
        pumping = df2['pumping']        
        x_array = np.array(x_coords)
        y_array = np.array(y_coords)
        pumping_array = np.array(pumping)
        Zw = (x_array + y_array * 1j)
        # print('These are Z of wells', Zw)
        alpha = 0

        return baseFlowX[0], x_coords[0], h0[0], H[0], k[0], y_coords[0], pumping_array, Zw , alpha, por[0]


    def phi_base_fxn_river_complex(self):   # calculation for Qox.X (base FlowX)
        vector = Vectors_river_complex()
        _, _, _, _, _, _, Z = vector.region_boundaries_river_complex()
        baseFlowX, x_coords, _, _, _, _, _, _, alpha, _ = self.potentials_data_river_complex()
        phi_base = -1 * baseFlowX* Z * np.exp(-1j * alpha )
        
        # print('These are phi_base_river_complex', phi_base.real)
        return phi_base


    def phi_0_river_complex(self):  # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        _ , _ , h0, H, k, _, _, _, _, _ = self.potentials_data_river_complex()
        
        
        if h0 > H:          # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location
            phi_0 = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0 = 0.5 *k *h0 * h0
        
        # print('this is phi_0 with river', phi_0)
        return phi_0 


    def calculation_phi_well_river_complex(self):
        vector = Vectors_river_complex()
        _, _, _, _, _, _, Z = vector.region_boundaries_river_complex()
        _, _, _, _, _, _, pumping_array, Zw, alpha, _ = self.potentials_data_river_complex()

        # Defining paramaters again due to the Method Images creating an aritificial well in opposite mirror direction (at -x = d)
        
        Theta_pumping_well = []
        Theta_recharge_well = []
        Modulus_pumping_well = []
        Modulus_recahrge_well = []
    
        phi_well_and_mirror_well_complex = 0

        # Calculation for Pumping Well Modulus
        for i in range(len(Zw)):
            modulus_pumping_well = np.absolute(Z-Zw[i])
            Modulus_pumping_well.append(modulus_pumping_well)
        # print('These are modulus pumping well', Modulus_pumping_well)

        # Well Images Creation 
        well_images_behind_river = Zw * -1
        corrected_well_image = np.conjugate(well_images_behind_river)
        # print('These are well Images for river comples', corrected_well_image)

        # Calculation for Recharge Well Modulus
        for i in range(len(Zw)):
            modulus_recharge_well = np.absolute(Z - (corrected_well_image[i]))
            Modulus_recahrge_well.append(modulus_recharge_well)

        # print('These are modulus of image recharge well behind river', Modulus_recahrge_well)    
        
        # Calculation for Angle/Argument at every mesh for pumping well
        for i in range(len(Zw)):
            theta1 = np.angle(Z-Zw[i])
            Theta_pumping_well.append(theta1)
        # print('These are theta values for pumping recharge well', Theta_pumping_well)
        # Calculation for Angle/Argument at every mesh for recharge well    
        for i in range(len(Zw)):
            theta2 = np.angle(Z-corrected_well_image[i])
            Theta_recharge_well.append(theta2)
        # print('These are theta values for image recharge well', Theta_recharge_well)

        # Here y will remain real as only the x distance is getting changed in mirror image not the y thus providing constant phi value along y axis
        # for i in range(len(Zw)):
        #     phi_wells = (pumping_array[i]/(2*np.pi))  * (np.log(Modulus_pumping_well[i]) + (1j * Theta_pumping_well[i]) - np.log(Modulus_recahrge_well[i]) - (1j * Theta_recharge_well[i]))
        #     # phi_wells = ( pumping_array[i] / (2 * np.pi)) * np.log(((Modulus_pumping_well[i])*(np.exp(1j * Theta_pumping_well[i]))) / ((Modulus_recahrge_well[i]))* np.exp(-1j*Theta_pumping_well[i]))
        #     # phi_R_I = ((Qr[i] / (2 * np.pi)) * (np.log(np.sqrt((X - (xreal[i]))**2 + (Y - (yreal[i]))**2) ))) + ((Qi[i] / (2 * np.pi)) * (np.log(np.sqrt((X - (ximagin[i]))**2 + (Y - (yreal[i]))**2) )))
        #     phi_well_and_mirror_well_complex += phi_wells
        for i in range(len(Zw)):
            phi_wells = (pumping_array[i]/(2*np.pi)) * (np.log(Z-Zw[i]) - np.log(Z-(corrected_well_image[i])))
            phi_well_and_mirror_well_complex += phi_wells
        # print (' this is phi_well of complex river', phi_well_and_mirror_well_complex)
        return phi_well_and_mirror_well_complex

    def dischargepotential_total_of_region_river_complex(self):
        phi_base = self.phi_base_fxn_river_complex()
        phi_0 = self.phi_0_river_complex()
        phi_well_and_mirro_well = self.calculation_phi_well_river_complex() 
        
        phi_river = phi_base + phi_0 + phi_well_and_mirro_well

        data1 = pd.DataFrame(phi_river.real)
        path1 = os.path.join("river", "phi_wells_with_river.xlsx")
        data1.to_excel(path1, sheet_name='sheet1', index=False, )

        data2 = pd.DataFrame(phi_river.imag)
        path2 = os.path.join("river", "psi_wells_with_river.xlsx")
        data2.to_excel(path2, sheet_name='sheet1', index=False, )
        # print(phi_river, 'these are all phi values with river')

        return phi_river


    def calculation_for_head_river_complex(self):
        h_list = [] # append all head values in this list


        phi_river1 = self.dischargepotential_total_of_region_river_complex()
        phi_river= phi_river1.real
        _ , _ , _, H, k, _, _, _, _, _ = self.potentials_data_river_complex()
        

        phicrit = 0.5 * k * H**2 

        for phi_item in phi_river:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Book page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            h_list.append(h)
        
        head= np.array(h_list)

        data3 = pd.DataFrame(head.real)
        path3 = os.path.join("river", "head_wells_with_river.xlsx")
        data3.to_excel(path3, sheet_name='sheet1', index=False, )
        # print('these are head values with river complex', head)
        return head

