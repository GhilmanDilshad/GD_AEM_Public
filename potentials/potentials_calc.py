import pandas as pd
import numpy as np
from vectorcalculations.vector_cords import Vectors
from utilities import read_aquifer_xlsx, read_wells_xlsx

class Potentials:


    def __init__(self):
        pass

    def potentials_data(self): #df1 is data frame for the aquifer
        
        #df1 = pd.read_excel(file_loc, index_col=None, na_values=['NA'], usecols = 'A,B:AA')     
        df1 = read_aquifer_xlsx()  # This is for reading values from aquifer file
        baseFlowX = df1['Base Flow in X direction']
        h0 = df1['Reference Head']
        H = df1['Aquifer Thickness'] # H is a variable storing value from H_thickness
        k = df1["Hydraulic Conductivity"] # k is a variable storing value from hydraulic conductivity 
        
        df2 = read_wells_xlsx()    # This is for reading values from well file
        x_coords = df2['x-coordinates']
        y_coords = df2['y-coordinates'] 
        pumping = df2['pumping']        
    
        # print(y_coords)
        return baseFlowX[0], x_coords[0], h0[0], H[0], k[0], y_coords[0], pumping[0]


    def phi_base_fxn(self):   # calculation for Qox.X (base FlowX)
        
        vector = Vectors()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries()

        baseFlowX, x_coords, _, _, _, _, _ = self.potentials_data()
        phi_base = -1*baseFlowX*X 
        
        # print(phi_base, 'This is Discharge Qox at base')

        return phi_base


    def phi_0(self):  # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        _ , _ , h0, H, k, _, _ = self.potentials_data()
        
        
        if h0 > H:          # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location
            phi_0 = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0 = 0.5 *k *h0 * h0
        return phi_0 


    def calculation_phi_well(self):
        phi_well = 0
        a_ref = 0 # Imaginary Reference at x grid
        b_ref = 0 # Imaginary Reference at y grid
        vector = Vectors()

        _, _, _, _, _, _, [X, Y] = vector.region_boundaries()
        _, _, x_array, y_array, pumping_array = vector.vectors_data()
        phi_base = self.phi_base_fxn() # This is baseflow -Qox

        # print(x_array, 'This is x coordinates array values of wells')
        for i in range(len(x_array)):
            #phi at all meshgrids
            phi_q = 0 + ((pumping_array[i] / (4 * np.pi))) * np.log(((X - x_array[i])**2 + (Y - y_array[i])**2) / ((a_ref - x_array[i])**2 + (b_ref - y_array[i])**2))
                
            phi_well += phi_q
        # print('These are values for non complex phi_well: ',phi_well)    
        return phi_well

    def dischargepotential_total_of_region(self):
        phi_base = self.phi_base_fxn()
        phi_0 = self.phi_0()
        phi_well = self.calculation_phi_well() 
        
        phi = phi_base + phi_0 + phi_well
        # print('These are phi values non complex withotu river',phi)
        return phi


    def calculation_for_head(self):
        h_list = [] # append all head values in this list

        phi = self.dischargepotential_total_of_region()
        _ , _ , _, H, k, _, _ = self.potentials_data()
        
        phicrit = 0.5 * k * H**2 
        # if phi >= phicrit:
        #     h = (phi_item + (0.5 * k * H**2)) / (k * H)
        # else:
        #     h = np.sqrt((2 * phi_item ) / k )
        # return h
        for phi_item in phi:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Book page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            
            h_list.append(h)
        
        head= np.array(h_list)
        # print('baba chal jaaaaa')
        # print('these are head values non complex ', head)
        
        return head
    

    def calculation_stream_fxn(self): # page 55 of book
        vector = Vectors()
        psi_well = 0

        baseFlowX, _, _, _, _, y_coords, _ = self.potentials_data()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries()
        _, _, x_array, y_array, pumping_array = vector.vectors_data()
              
        psi_base = -baseFlowX * Y

        for i in range(len(x_array)):
            psi_q =  (pumping_array[i] / (2 * np.pi)) * (np.arctan2((Y - y_array[i]), (X - x_array[i])) )
            psi_well += psi_q       
        
        psi = psi_base + psi_well
        # print('These are values for psi non complex for well', psi_well)
        # print(' this is psi for non complex', psi)
        return psi
        