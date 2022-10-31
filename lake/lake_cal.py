import pandas as pd
import numpy as np
from utilities import read_aquifer_xlsx, read_wells_xlsx

class Vectors_Lake:
    def __init__(self):  
        pass    

    def region_boundaries_lake(self): # Region Coordinates, boundaries and meshgrids.
        Xmin = -200
        Xmax = 200
        Ymin = 0
        Ymax = 200
        xmesh = np.linspace(Xmin, Xmax, 351)
        ymesh = np.linspace(Ymin, Ymax, 351)
        [X, Y] = np.meshgrid(xmesh, ymesh)

        return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, [X, Y]

    def vectors_data_lake(self):  
        df3 = read_wells_xlsx()  
       
        x_values = df3['x-coordinates'].tolist()   
        y_values = df3['y-coordinates'].tolist()
        pumping_rate = df3['pumping'].tolist()
        pumping_array = np.array(pumping_rate)

        x_array = np.array(x_values) # converting list to array
        y_array = np.array(y_values) # converting list to array
        return x_values, y_values, x_array, y_array, pumping_array



class Potentials_Lake:

    def __init__(self):
        pass
       
    
    def potentials_data_lake(self): #df1 is data frame for the aquifer
        
   
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
             
        # Lake Radius and Location
        R_Lake = 25     #Radius of Lake
        x_Lake = 100  #lake_location_at_x_axis
        y_Lake = [100]    #lake_location_at_y_axi
        Q_lake = 250 #Aassuming Discharge rate of Lake

        # take_a_user_input_of_reference_head_lake_or_use_the_first_reference_head_as_the_reference=0
    
        return baseFlowX[0], x_coords[0], h0[0], H[0], k[0], y_coords[0], pumping[0], por[0], R_Lake, x_Lake, y_Lake, Q_lake


    def phi_base_fxn_lake(self):   # calculation for Qox.X (base FlowX)
        
        vector = Vectors_Lake()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_lake()
        baseFlowX, _, _, _, _, _, _, _, R_Lake, x_Lake, y_Lake, _ = self.potentials_data_lake()
        
        
        phi_base =  -1 * baseFlowX * ((X - x_Lake) - (((X - x_Lake)*(R_Lake)**2)) / ((X - x_Lake)**2 + (Y - y_Lake )**2) ) # Here -1 will stay it is part of formula 
        

        # print( 'This is Discharge uniform flow field Lake Non Complex', phi_base)
        return phi_base


    def phi_0_lake(self):  # by using condition from book at page number39 #this is discharge potential for confined and unconfined 
        
        _ , _ , h0, H, k, _, _, _, _, _, _, Q_lake = self.potentials_data_lake()

        
        if h0 > H:          # this is discharge potential phi_0 at the reference head, defining the aquifer is confined or unconfined at the refernce head location
            phi_0 = k * H * h0 - 0.5 *k *H * H
        else:
            phi_0 = 0.5 *k *h0 * h0
        

        # print('this is phi_0 with Lake', phi_0)
        return phi_0 


    def calculation_phi_well_lake(self):
        vector = Vectors_Lake()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_lake()
        _, _, x_array, y_array, pumping_array = vector.vectors_data_lake()
        _, _, _, _, _, _, _, _, R_Lake, x_Lake, y_Lake, _ = self.potentials_data_lake()

        # Defining paramaters again due to the Method Images creating an aritificial well in opposite mirror direction (at -x = d)
        Qr = pumping_array
        Qi = -1 * pumping_array
        xreal= x_array
        ximagin=-1*x_array
        yreal= y_array
    
        phi_well_lake = 0
        
        # Here y will remain real as only the x distance is getting changed in mirror image not the y thus providing constant phi value along y axis
        # for i in range(len(Qr)):
            
        #     phi_wells_lake_q = (Qr[i] / (4 * np.pi)) * np.log ((((X-(R_Lake**2/(ximagin[i])))**2 + (Y - (yreal[i]))**2) / ((X - (ximagin[i]))**2 + (Y - (yreal[i]))**2 )) * ((ximagin[i])**2 / (R_Lake**2)))

        #     phi_well_lake += phi_wells_lake_q
            # break # only to make it run one time 
        x_Lake_for_image_creation = []

        for i in range(len(x_array)):
            image_creation_of_wells = (1 * x_Lake)+0.4
            x_Lake_for_image_creation.append(image_creation_of_wells)
        
        
        for i in range(len(Qr)):
            
            phi_wells_lake_q = (Qr[i] / (4 * np.pi)) * np.log ((((X-(x_array[i]))**2 + (Y - (yreal[i]))**2) / ((X - (x_Lake))**2 + (Y - (y_Lake))**2 )) * ((x_Lake_for_image_creation[i])**2 / (R_Lake**2)))

            phi_well_lake += phi_wells_lake_q
            # break # only to make it run one time 
        # print('this is non complex phi_of_well_with lake', phi_well_lake ) 
        return phi_well_lake


    def the_constant_C(self):
        phi_0 = self.phi_0_lake()
        vector = Vectors_Lake()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_lake()
        _, _, x_array, y_array, pumping_array = vector.vectors_data_lake()
        _, _, _, _, _, _, _, _, R_Lake, x_Lake, y_Lake, _ = self.potentials_data_lake()
        
        C = 0 

        for i in range(len(x_array)):
            C_constant = (pumping_array[i]/(2*np.pi)) * (np.log((R_Lake)/ (np.sqrt((X-x_array[i])**2+(Y-y_array[i])**2))))

            C += C_constant

        C_final_constant =  C + phi_0
        # print('This is C the final constant',C_final_constant)
        return C_final_constant


    def calculation_inflow_outflow_lake(self):
        vector = Vectors_Lake()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_lake()
        _, _, _, _, _, _, _, _, R_Lake, x_Lake, y_Lake,Q_lake = self.potentials_data_lake()

        phi_inflow_outflow =  ((Q_lake / (4 * np.pi) ) * np.log (((X - x_Lake )**2 + (Y - y_Lake)**2) / (R_Lake)**2))
        
        # print('this is phi inflowoutflow non complex', phi_inflow_outflow)
        return phi_inflow_outflow



    def dischargepotential_total_of_region_lake(self):
        phi_base = self.phi_base_fxn_lake()
        phi_0 = self.phi_0_lake()
        phi_well_lake = self.calculation_phi_well_lake() 
        phi_inflow_outflow = self.calculation_inflow_outflow_lake()
        C_final_constant = self.the_constant_C()

        phi_lake_total = phi_well_lake + (phi_base) + (phi_inflow_outflow) + C_final_constant 
        # print(phi_lake_total, 'these are all phi values with river')
        return phi_lake_total


    def calculation_for_head_lake(self):
        h_list = [] # append all head values in this list


        phi = self.dischargepotential_total_of_region_lake()
        _ , _ , _, H, k, _, _, _, _, _, _, _ = self.potentials_data_lake()
        

        phicrit = 0.5 * k * H**2 

        for phi_item in phi:
            for sub_phi_item in phi_item:
                if sub_phi_item >= phicrit:
                    h = (phi_item + (0.5 * k * H**2)) / (k * H)   # Book page number 38
                else:
                    h = np.sqrt((2 * phi_item ) / k )
            h_list.append(h)
        
        head = np.array(h_list)
        # print('these are head values with river', head)
        return head
    


    def calculation_stream_fxn_lake(self): 
        vector = Vectors_Lake()
        _, _, _, _, _, _, [X, Y] = vector.region_boundaries_lake()
        _, _, x_array, y_array, pumping_array = vector.vectors_data_lake()   

        baseFlowX, _, _, _, _, _, _, _, R_Lake, x_Lake, y_Lake, Q_lake = self.potentials_data_lake()

        
        Qr = pumping_array
        xreal= x_array
        ximagin=-1*x_array
        yreal= y_array

        psi_base01 = 0
        psi_base02 = 0
        psi_base03 = 0    
        

        #psi base # lake in uniform flow field
        for i in range(len(x_array)):
            psi_base_q = -1 * baseFlowX * ((Y - y_Lake) + (((Y - y_Lake)*(R_Lake**2)) /((X-x_Lake)**2 + (Y - y_Lake)**2 )) )
            psi_base01 += psi_base_q
            # break
        # print('this is psi_base01', psi_base01)
            #net inflow out flow # page 76 haitjema, page 260 strack
        for i in range(len(x_array)):
            psi_base_w = (Q_lake / (2 * np.pi)) * (np.arctan2( (Y - y_Lake) , (X - x_Lake) ))
            psi_base02 += psi_base_w
            # break
        # print('this is psi_base02', psi_base02)
        # psi due to well 3.137 page 72                                        #(X - ((R_Lake**2) / (ximagin[i])))) - np.arctan2((Y - (y_array[i])), (X - (ximagin[i]))))
        for i in range(len(x_array)):
            psi_well =  (Qr[i] / (2 * np.pi)) * (np.arctan2((Y - (y_array[i])), (X - (xreal[i]))) - np.arctan2((Y - (y_array[i])), (X - (ximagin[i]))))
            psi_base03 += psi_well       
        #     # break           


        # I am trying here something lets see it works or know

        # for i in range(len(x_array)):
        #     psi_well =  (Qr[i] / (2 * np.pi)) * (np.arctan2((Y - (y_array[i])), (X - (xreal[i]))) - np.arctan2((y_Lake - (y_array[i])), (x_Lake - (ximagin[i]))))
        #     psi_base03 += psi_well       
        #     # break           




        # print('this is psi_base03', psi_base03)
        psi_total_lake_plus_river = (psi_base01) + (psi_base02) + (psi_base03)

        return psi_total_lake_plus_river