from cmath import exp
from email.mime import image
import pandas as pd
import numpy as np
from utilities import read_aquifer_xlsx, read_wells_xlsx
import os

class Vectors_Line_Complex:
    def __init__(self):  
        pass    

    def region_boundaries_line_complex(self): # Region Coordinates, boundaries and meshgrids.
        Xmin = -400
        Xmax = 400
        Ymin = 0
        Ymax = 200
        xmesh = np.linspace(Xmin, Xmax, 401)
        ymesh = np.linspace(Ymin, Ymax, 401)
        [X_axis, Y_axis] = np.meshgrid(xmesh, ymesh)
        Z = (X_axis + Y_axis* 1j) 

        return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis

    def vectors_data_line_complex(self):  
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


class Potentials_Line_Complex:


    def __init__(self):
        pass
       

    def potentials_data_line_complex(self): #df1 is data frame for the aquifer
        
        # Line Starting and End Point
        x1_Line = [200]     #Starting Point of Line
        x2_Line = [200]     #Ending_Point_Of_Line
        y1_Line = [50]
        y2_Line = [150]
        # Extraction (Positive) or Inflitration (Negative) Rate of Line-Sink
        Sigma_line = [-5] 
               
        x1 = np.array(x1_Line)
        x2 = np.array(x2_Line)
        y1 = np.array(y1_Line)    
        y2 = np.array(y2_Line)
        z1_Line = (x1 + y1 * 1j)
        z2_Line = (x2 + y2 * 1j)
        # print(z1_Line)
        # print(z2_Line)
        # # Height of Creek Bed Above Base Page 316 Strack and it Should Always Be less then H 
        # H0 = [18]
        # # Resistance of Creek Bed of Unit (m/d) page 316 Strack has mistake in unit of c
        # c = [5]
        # # width of Creek Unit (m)
        # b = [25]
        # phi_creek_0 = [20]
        # # Calculation for Sigma 
        # Sigma_line = 

        return  z1_Line, z2_Line, Sigma_line

    def phi_0_line_complex(self):  # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        vector = Vectors_Line_Complex()
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vector.vectors_data_line_complex()
        
        # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location # This Needs to be changed according to strack 
        if h0 > H:          
            phi_0_line_complex = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0_line_complex = 0.5 *k *h0 * h0
        
        # print('this is phi_0 with Lake', phi_0)
        return phi_0_line_complex


    def phi_uniform_flow_field_line_complex(self):   # calculation for Qox.X (base FlowX)
        
        vector = Vectors_Line_Complex()
        _, _, _, _, _, _, Z,_, _ = vector.region_boundaries_line_complex()
        baseFlowX, _, _, _, _, _, _, _, alpha, _, _ = vector.vectors_data_line_complex()
        z1_Line, z2_Line, Sigma_line = self.potentials_data_line_complex()

        phi_uniform_flow_field_line_complex =  (-1 * baseFlowX ) * Z * (np.exp(-1j * alpha))

        return phi_uniform_flow_field_line_complex


    def calculation_inflow_outflow_line_complex(self):
        vector = Vectors_Line_Complex()
        _, _, _, _, _, _, Z, _, _ = vector.region_boundaries_line_complex()
        baseFlowX, _, _, _, _, _, _, _, alpha, _, por = vector.vectors_data_line_complex()
        z1_Line, z2_Line, Sigma_line = self.potentials_data_line_complex()

        phi_inflow_outflow_line_complex = 0
        Capital_Z = []
        Lentgh_of_Creek = []
        Sigma_line_array = np.array(Sigma_line)

        # For Calculation of lentgh of Creek
        for i in range(len(Sigma_line_array)):
            lentgh_of_creek = np.absolute(z2_Line[0]-z1_Line[0]) 
            Lentgh_of_Creek.append(lentgh_of_creek)
        # print('This is Lentgh of Creek',Lentgh_of_Creek)

        # calculation for capital Z
        for i in range(len(Sigma_line_array)):
            capital_Z = (Z - 0.5 * (z1_Line[0] + z2_Line[0])) / (0.5 * (z2_Line[0] - z1_Line[0]))
            Capital_Z.append(capital_Z)
        
        for i in range(len(Sigma_line_array)):
           phi_inflow_outflow_line_complex_q =  ((Sigma_line_array[0] * Lentgh_of_Creek[0]) / (4 * np.pi)) * (((Capital_Z[0] + 1)*np.log(Capital_Z[0] + 1)) - ((Capital_Z[0] - 1)*np.log(Capital_Z[0] - 1)) - 2 )
           
           phi_inflow_outflow_line_complex += phi_inflow_outflow_line_complex_q   

        return phi_inflow_outflow_line_complex
   
   
    def calculation_phi_well_line_complex(self):
        vector = Vectors_Line_Complex()
        _, _, _, _, _, _, Z, _, _ = vector.region_boundaries_line_complex()
        baseFlowX, _, _, _, _, _, pumping_array, Zw, alpha, _, por = vector.vectors_data_line_complex()
        z1_Line, z2_Line, Sigma_line = self.potentials_data_line_complex()
        Zref = 50 + 100j
        phi_well_line_complex = 0
        Coordinates_for_Image_creation_location= 0
        Images_of_wells = []
  

        coordinates_for_Image_creation_location = ((z1_Line[0])+(z2_Line[0]))/2

        Coordinates_for_Image_creation_location+=coordinates_for_Image_creation_location
        
        # print('These are Coordinates for midpoint of line sink', Coordinates_for_Image_creation_location)
        for i in range(len(Zw)):
            images_of_wells = Coordinates_for_Image_creation_location * 1
            Images_of_wells.append(images_of_wells)
        # Creation of Images at particular center of line


        # for i in range(len(Zw)):
        #     images_of_wells = 1 * Coordinates_for_Image_creation_location
        #     Images_of_wells.append(images_of_wells)
        # print('These are Coordinates for images of wells with line sink', Images_of_wells)
        
        for i in range(len(Zw)):
            phi_wells = (pumping_array[i]/(2*np.pi)) * (np.log(Z-Zw[i]) - np.log(Z-(Images_of_wells[i])))
            phi_well_line_complex += phi_wells
        # print (' this is phi_well of complex line source', phi_well_line_complex)
        return phi_well_line_complex


    def dischargepotential_total_of_region_line_complex(self):
        phi_uniform_flow_field = self.phi_uniform_flow_field_line_complex()
        phi_0_complex = self.phi_0_line_complex()
        phi_well_line_complex = self.calculation_phi_well_line_complex() 
        phi_inflow_outflow_line_complex = self.calculation_inflow_outflow_line_complex()
        # C_the_constant= self.the_constant_C()
        # phi_infinity = self.phi_far_field_R_infinity


        phi_line_total_complex =    (phi_uniform_flow_field) +  (phi_well_line_complex) + (phi_inflow_outflow_line_complex) + [phi_0_complex]         
        # phi_line_total_complex =    (phi_uniform_flow_field) +  (phi_well_lake_complex) + (phi_inflow_outflow_complex) + C_the_constant       
        
        data1 = pd.DataFrame(phi_line_total_complex.real)
        path1 = os.path.join("linesink", "phi_well_line_without_river_complex.xlsx")
        data1.to_excel(path1, sheet_name='sheet1', index=False, )
        
        data2 = pd.DataFrame(phi_line_total_complex.imag)
        path2 = os.path.join("linesink", "psi_well_line_without_river_complex.xlsx")
        data2.to_excel(path2, sheet_name='sheet1', index=False, )             
        # print(phi_lake_total, 'these are all phi values with river')
        return phi_line_total_complex


    def calculation_for_head_line_complex(self):
        vector = Vectors_Line_Complex()
        _, _, _, H, k, _, _, _, _, _, _ = vector.vectors_data_line_complex()
        phi = self.dischargepotential_total_of_region_line_complex()

        h_list = [] # append all head values in this list       
        phicrit = 0.5 * k * H**2 

        for phi_item in phi:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Book page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            h_list.append(h)
        
        head_complex = np.array(h_list)
        
        data3 = pd.DataFrame(head_complex.real)
        path3 = os.path.join("linesink", "head_well_line_without_river_complex.xlsx")
        data3.to_excel(path3, sheet_name='sheet1', index=False, )       
        # print('these are head values of complex lake', head_complex)
        return head_complex