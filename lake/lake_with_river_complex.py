from cmath import exp
import pandas as pd
import numpy as np
from utilities import read_aquifer_xlsx, read_wells_xlsx, read_lake_xlsx
import os

class Vectors_Lake_River_Complex:
    def __init__(self):  
        pass    

    def region_boundaries_lake_river_complex(self): # Region Coordinates, boundaries and meshgrids.
        Xmin = -600
        Xmax = 600
        Ymin = 0
        Ymax = 200
        xmesh = np.linspace(Xmin, Xmax, 601)
        ymesh = np.linspace(Ymin, Ymax, 601)
        [X_axis, Y_axis] = np.meshgrid(xmesh, ymesh)
        Z = (X_axis + Y_axis* 1j) 

        return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis

    def vectors_data_lake_river_complex(self):  
        df1 = read_aquifer_xlsx()  # This is for reading values from aquifer file
        baseFlowX = df1['Base Flow in X direction']
        h0 = df1['Reference Head']
        H = df1['Aquifer Thickness'] # H is a variable storing value from H_thickness
        k = df1["Hydraulic Conductivity"] # k is a variable storing value from hydraulic conductivity 
        por = df1["Porosity"]

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
        
        return baseFlowX[0], x_array, h0[0], H[0], k[0], y_array, pumping_array, Zw, alpha, Zref, por[0]


class Potentials_Lake_River_Complex:


    def __init__(self):
        pass
       

    def potentials_data_lake_river_complex(self): #df1 is data frame for the aquifer
        dataframe = read_lake_xlsx()
        radius = dataframe['radius']
        x_coordinates = dataframe['x-coordinates']
        y_coordinates = dataframe['y-coordinates']
        inflow_outflow = dataframe['inflow_outflow']
            
        # Lake Radius and Location
        R_Lake = radius.to_list()           #Radius of Lake
        x_Lake = x_coordinates.to_list()    #lake_location_at_x_axis
        y_Lake = y_coordinates.to_list()    #lake_location_at_y_axi
        Q_lake = inflow_outflow.to_list()   #Aassuming Discharge rate of LakeW

        
        x_lake_array = np.array(x_Lake)
        y_lake_array = np.array(y_Lake)    
        Z_Lake = (x_lake_array + y_lake_array * 1j)
        

        return  R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake

    def phi_0_lake_river_complex(self):     # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        vector = Vectors_Lake_River_Complex()
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vector.vectors_data_lake_river_complex()
                

        # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location # This Needs to be changed according to strack 
        if h0 > H:          
            phi_0_lake_complex = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0_lake_complex = 0.5 *k *h0 * h0
        
        # print('this is phi_0 with Lake', phi_0_lake_complex)
        return phi_0_lake_complex


    def the_constant_C_river(self):
        phi_0 = self.phi_0_lake_river_complex()
        vector = Vectors_Lake_River_Complex()
        _, _, _, _, _, _, pumping_array, Zw, _, _, _ = vector.vectors_data_lake_river_complex()
        R_Lake, _, _, _, Z_Lake = self.potentials_data_lake_river_complex()
        _, _, _, _, _, _, Z, _, _ = vector.region_boundaries_lake_river_complex()
        # phi_infinity = self.phi_far_field_R_infinity_river()
        C = 0 

        for i in range(len(Zw)):
            C_constant = (pumping_array[i]/(2*np.pi)) * (np.log((R_Lake)/ (np.absolute(Z - Zw[i]))))
            C += C_constant

        C_final_constant =   phi_0 
        return C_final_constant


    def phi_uniform_flow_field_river_complex(self):   # calculation for Qox.X (base FlowX)
        
        vector = Vectors_Lake_River_Complex()
        _, _, _, _, _, _, Z,_, _ = vector.region_boundaries_lake_river_complex()
        baseFlowX, _, _, _, _, _, _, _, alpha, _, _ = vector.vectors_data_lake_river_complex()
        R_Lake, _, _, _, Z_Lake = self.potentials_data_lake_river_complex()
        R_lake_array = np.array(R_Lake)
        

        phi_uniform_flow_field_complex =  (-1 * baseFlowX ) * ((Z - Z_Lake) - ((R_lake_array**2) / (Z - Z_Lake))) * (np.exp(-1j * alpha))

        return phi_uniform_flow_field_complex

    def calculation_inflow_outflow_lake_river_complex(self):
        vector = Vectors_Lake_River_Complex()
        _, _, _, _, _, _, Z, _, _ = vector.region_boundaries_lake_river_complex()
        baseFlowX, _, _, _, _, _, _, _, alpha, _, por = vector.vectors_data_lake_river_complex()
        R_Lake, _, _, Q_lake, Z_Lake = self.potentials_data_lake_river_complex()
        R_infinity =  25**2 / 0.4 

        phi_inflow_outflow_complex = 0
        modulus_Z_and_lake_list = []
        Theta = []
        R_lake_array = np.array(R_Lake)
        Q_lake_array = np.array(Q_lake)
        
        # calculation for r 
        for i in range(len(Z_Lake)):
            r = np.absolute(Z-(Z_Lake[i]))
            modulus_Z_and_lake_list.append(r)
            # print('this are complexed r values', r)
        for i in range(len(Z_Lake)):
            theta = np.angle(Z-Z_Lake[i])
            Theta.append(theta)
            # print('these are theta values complexed', Theta)

        for i in range(len(Q_lake_array)):
           phi_inflow_outflow_complex_q =  ((Q_lake_array / (2 * np.pi) ) * (np.log((Z-Z_Lake)/R_infinity)))  
           phi_inflow_outflow_complex += phi_inflow_outflow_complex_q   

        return phi_inflow_outflow_complex
   
    def calculation_phi_well_lake_river_complex(self):
        vector = Vectors_Lake_River_Complex()
        _, _, _, _, _, _, Z, _, _ = vector.region_boundaries_lake_river_complex()
        baseFlowX, x_array, _, _, _, y_array, pumping_array, Zw, alpha, _, por = vector.vectors_data_lake_river_complex()
        R_Lake, _, _, Q_lake, Z_Lake = self.potentials_data_lake_river_complex()
        
        phi_well_lake_complex = 0
        R_lake_array = np.array(R_Lake)
       
        well_image_list = []
        trying_well_images = Zw * -1
        actual_well_images = np.conjugate(trying_well_images)
            
        for i in range(len(pumping_array)):
            phi_well_lake_q = (pumping_array[i] / (2*np.pi)) * np.log((((Z)-Zw[i]) / ((Z) - (actual_well_images[i])) * (R_lake_array/np.absolute(Zw[i]))))
            
            phi_well_lake_complex += phi_well_lake_q

        return phi_well_lake_complex

    def dischargepotential_total_of_region_lake_river_complex(self):
        phi_uniform_flow_field = self.phi_uniform_flow_field_river_complex()
        phi_0_complex = self.phi_0_lake_river_complex()
        phi_well_lake_complex = self.calculation_phi_well_lake_river_complex() 
        phi_inflow_outflow_complex = self.calculation_inflow_outflow_lake_river_complex()
        C_the_constant= self.the_constant_C_river()


        phi_lake_total_complex =    (phi_uniform_flow_field) +  (phi_well_lake_complex) + (phi_inflow_outflow_complex) + C_the_constant
        
        
        data1 = pd.DataFrame(phi_lake_total_complex.real)
        path1 = os.path.join("lake", "phi_well_lake_with_river_complex.xlsx")
        data1.to_excel(path1, sheet_name='sheet1', index=False, )
        
        data2 = pd.DataFrame(phi_lake_total_complex.imag)
        path2 = os.path.join("lake", "psi_well_lake_with_river_complex.xlsx")
        data2.to_excel(path2, sheet_name='sheet1', index=False, )


        return phi_lake_total_complex


    def calculation_for_head_lake_river_complex(self):
        vector = Vectors_Lake_River_Complex()
        _, _, _, H, k, _, _, _, _, _, _ = vector.vectors_data_lake_river_complex()
        phi = self.dischargepotential_total_of_region_lake_river_complex()
        

        h_list = [] # append all head values in this list       
        phicrit = 0.5 * k * H**2 

        for phi_item in phi:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Strack page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            h_list.append(h)
        
        head_complex = np.array(h_list)

        data3 = pd.DataFrame(head_complex.real)
        path3 = os.path.join("lake", "head_well_lake_with_river_complex.xlsx")
        data3.to_excel(path3, sheet_name='sheet1', index=False, )

        return head_complex