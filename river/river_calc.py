import pandas as pd
import numpy as np
from sympy import Q
from utilities import read_aquifer_xlsx, read_wells_xlsx

#Defining Vector Class Again Due to Seprate or Bigger Regio Required Changing the vector Parameters from the vector file. So a new Vector Class for river
class Vectors_river:
    def __init__(self):  
        pass    

    def region_boundaries_river(self): # Region Coordinates, boundaries and meshgrids.
        Xmin = -200
        Xmax = 200
        Ymin = 0
        Ymax = 200
        xmesh = np.linspace(Xmin, Xmax, 201)
        ymesh = np.linspace(Ymin, Ymax, 201)
        [X, Y] = np.meshgrid(xmesh, ymesh)

        return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, [X, Y]

    def vectors_data_river(self):    # getting coordinates of well(s) and converting them to list 
        df3 = read_wells_xlsx()   # This is for reading values from well file
       
        x_values = df3['x-coordinates'].tolist()   
        y_values = df3['y-coordinates'].tolist()
       
        pumping_rate = df3['pumping'].tolist()
        pumping_array = np.array(pumping_rate)

        x_array = np.array(x_values) # converting list to array
        y_array = np.array(y_values) # converting list to array
        return x_values, y_values, x_array, y_array, pumping_array


class Potentials_River:
    def __init__(self):
        pass

    def potentials_data_river(self): #df1 is data frame for the aquifer
        
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
    
        return baseFlowX[0], x_coords[0], h0[0], H[0], k[0], y_coords[0], pumping[0], por[0]


    def phi_base_fxn_river(self):   # calculation for Qox.X (base FlowX)
        
        vector = Vectors_river()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_river()


        baseFlowX, x_coords, _, _, _, _, _, _ = self.potentials_data_river()
        phi_base = -1 * baseFlowX*X 
        
        # print(phi_base, 'This is Discharge Qox at base with river non complex')
        return phi_base


    def phi_0_river(self):  # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        _ , _ , h0, H, k, _, _, _ = self.potentials_data_river()
        
        
        if h0 > H:          # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location
            phi_0 = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0 = 0.5 *k *h0 * h0
        
        # print('this is phi_0 with river', phi_0)
        return phi_0 


    def calculation_phi_well_river(self):
        vector = Vectors_river()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_river()
        _, _, x_array, y_array, pumping_array = vector.vectors_data_river()

        # Defining paramaters again due to the Method Images creating an aritificial well in opposite mirror direction (at -x = d)
        Qr = pumping_array
        Qi = -1 * pumping_array
        xreal= x_array
        ximagin=-1*x_array
        yreal= y_array
    
        phi_well_and_mirro_well = 0


        # Here y will remain real as only the x distance is getting changed in mirror image not the y thus providing constant phi value along y axis
        for i in range(len(Qr)):
            
            phi_wells = ( pumping_array[i] / (4 * np.pi)) * np.log(((X - (xreal[i]))**2 + (Y - (yreal[i]))**2) / ((X - (ximagin[i]))**2 + (Y - (yreal[i]))**2))
            # phi_R_I = ((Qr[i] / (2 * np.pi)) * (np.log(np.sqrt((X - (xreal[i]))**2 + (Y - (yreal[i]))**2) ))) + ((Qi[i] / (2 * np.pi)) * (np.log(np.sqrt((X - (ximagin[i]))**2 + (Y - (yreal[i]))**2) )))
            phi_well_and_mirro_well += phi_wells


        # print (' this is phi_well non complex', phi_well_and_mirro_well)
        return phi_well_and_mirro_well

    def dischargepotential_total_of_region_river(self):
        phi_base = self.phi_base_fxn_river()
        phi_0 = self.phi_0_river()
        phi_well_and_mirro_well = self.calculation_phi_well_river() 
        
        phi_river = phi_base + phi_0 + phi_well_and_mirro_well

        data = pd.DataFrame(phi_river)
        data.to_excel('phi_well_with_river.xlsx', sheet_name='sheet1', index=False, )
        # print(phi_river, 'these are all phi values with river')
        return phi_river


    def calculation_for_head_river(self):
        h_list = [] # append all head values in this list


        phi_river = self.dischargepotential_total_of_region_river()
        _ , _ , _, H, k, _, _, _ = self.potentials_data_river()
        

        phicrit = 0.5 * k * H**2 

        for phi_item in phi_river:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Book page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            h_list.append(h)
        
        head= np.array(h_list)
        # print('these are head values with river non complex', head)
        return head
    

    def calculation_stream_fxn_river(self): # page 60 of book
        vector = Vectors_river()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_river()
        _, _, x_array, y_array, pumping_array = vector.vectors_data_river()   
        baseFlowX, _, _, _, _, _, _, _ = self.potentials_data_river()

        Q = pumping_array
        xreal= x_array
        ximagin=-1*x_array
        yreal= y_array

        psi_well = 0    
        psi_base = -baseFlowX * Y

        for i in range(len(x_array)):
            psi_q =  (Q[i] / (2 * np.pi)) * (np.arctan2((Y - (y_array[i])), (X - (xreal[i]))) - np.arctan2((Y - (y_array[i])), (X - (ximagin[i]))))
            psi_well += psi_q       
        
        psi = psi_base + psi_well
        # print(' this is psi with river', psi)s
        return psi

