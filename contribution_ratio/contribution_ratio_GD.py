from email.mime import base
from symtable import Symbol
from xml.etree.ElementInclude import DEFAULT_MAX_INCLUSION_DEPTH
import pandas as pd
import numpy as np
import sympy as sp
from utilities import read_aquifer_xlsx, read_wells_xlsx
import os
from river.river_complex import Potentials_River_complex
from river.river_complex import Vectors_river_complex

from lake.lake_complex import Vectors_Lake_Complex
from lake.lake_complex import Potentials_Lake_Complex
from lake.lake_with_river_complex import Vectors_Lake_River_Complex
from lake.lake_with_river_complex import Potentials_Lake_River_Complex

class contribution_ratio_river_and_well_strack_GD:
    def __init__(self): 
        pass  

    def limit_of_stagnations_points_by_strack(self):

        data_reading = Potentials_River_complex()
        baseFlowX, _, h0, H, k, _, pumping_array, Zw , alpha, por = data_reading.potentials_data_river_complex()

        condition_for_stagnation_point_to_exist = []

        limit = ( pumping_array[0] ) / (np.pi * (Zw.real[0]) * baseFlowX)
        condition_for_stagnation_point_to_exist.append(limit)

        # print('Method By Strack This is condition for stagnation point to exist GD ',condition_for_stagnation_point_to_exist)
        return condition_for_stagnation_point_to_exist

    
    def calculation_for_stagnation_points_formula_by_holzbecher_and_strack(self): 
        data_reading = Potentials_River_complex()
        baseFlowX, _, h0, H, k, _, pumping_array, Zw , alpha, por = data_reading.potentials_data_river_complex()
   
        
        stagnation_points = []
        limit = self.limit_of_stagnations_points_by_strack()
        x = np.array(limit)
        y = sp.symbols('y')
        well_distance = np.array(Zw.imag[0])

        Ys = sp.Symbol('Ys')

        well_images_behind_river = Zw * -1
        corrected_well_image1 = np.conjugate(well_images_behind_river)
        print(corrected_well_image1)
        print(Zw)
        gs=sp.Symbol('gs')
        new_image = (50+100j)
        
        
        # eq_new = sp.Eq(( baseFlowX - (((pumping_array[0]) / (2*np.pi)) * ((1/(gs-new_image)) - (1/(gs-Zw[0]))) ) ),0 )   
# 
        # eq_new = sp.Eq(( baseFlowX - (((pumping_array[0]) / (2*np.pi)) * ((1/(gs-corrected_well_image1[0])) - (1/(gs-Zw[0]))) ) ),0 )   
        # eq_new_1 = -baseFlowX
        # eq_new_2 = (pumping_array[0]/(2*np.pi)) * (1/((gs)+(-1*(Zw[0])))) 
        # eq_new_3 = -1* (pumping_array[0]/(2*np.pi)) * (1/((gs)+(-1*(corrected_well_image1[0])))) 
        # eq_new   =  sp.Eq((eq_new_1+eq_new_2+eq_new_3),0)
        # solve = sp.solve(eq_new)
        # solve1 = np.fromiter(solve,dtype=complex)
        # print('sssssssssssssssssssssssssssssssssssssssssssssssssssssss',solve1)

        # calculating_discharge_by_this_approach_1 =   (-baseFlowX * solve1[0]) - (pumping_array[0]/(2*np.pi)) * ((np.log(solve1[0]-Zw[0])) - (np.log(solve1[0]-corrected_well_image1[0])))

       
        # print('1stttttttttttttttttttt',calculating_discharge_by_this_approach_1)
        # calculating_discharge_by_this_approach_2 =   (-baseFlowX * solve1[1]) - (pumping_array[0]/(2*np.pi)) * ((np.log(solve1[1]-Zw[0])) - (np.log(solve1[1]-corrected_well_image1[0])))      
        # print('2ndddddddddddddddddddd',calculating_discharge_by_this_approach_2)

       # Staaaack
        if x >= 1:
            print("Stagnation Points Exists Ghilman")

            stagnation_point_by_strack = sp.Eq(Ys**2, ((Zw.real[0])**2) +  (((Zw.real[0])*pumping_array[0])/(np.pi*baseFlowX)))
        
            stagnation_point_by_strack_solution = sp.solve(stagnation_point_by_strack)
            # print('These are two stagnation Points Not Corrected By Strack', stagnation_point_by_strack_solution)

            corrected_stagnation_points_by_strack = (stagnation_point_by_strack_solution + well_distance )
            print('Stracccck',corrected_stagnation_points_by_strack)
        else: 
            print('No Stagnation Points Exist')
        
        #       HOLZBECHER
        if x >= 1:
            print("Stagnation Points Exists Ghilman")

            Equation_for_calc_stag_points = sp.Eq(y**2, -((Zw.real[0])**2) +  (((Zw.real[0])*pumping_array[0])/(np.pi*baseFlowX)))
            two_stagnation_points = sp.solve(Equation_for_calc_stag_points)
            # print('These are two stagnation points by formula of Holzbecher',two_stagnation_points)
            corrected_stagnation_points_by_holzbecher = (two_stagnation_points +well_distance)
            print('Holzbecher', corrected_stagnation_points_by_holzbecher)

        else:
            print("No Stagnation Points Exist")
        
        # here I am trying stagnation points by whole formula 
        if x >= 1:
            print("Stagnation Points Exist ")

            Zs = sp.symbols('Zs')
            unseprated_equation_strack = sp.Eq((baseFlowX - ( (pumping_array[0]/(2*np.pi)) * (1/(Zs-Zw[0])) ) + ( (pumping_array[0]/(2*np.pi)) * (1/(Zs-corrected_well_image1[0])) ) ),0)
            Unserpated_solution_combine = sp.solve(unseprated_equation_strack)
            print('by whole formula',Unserpated_solution_combine)
            corrected_unseprated = Unserpated_solution_combine + Zw[0].imag
            print('correction stagnaton points ',corrected_unseprated)
        else:
            print("No Stagnation Points Exist")
            corrected_stagnation_points_by_holzbecher = []
            corrected_stagnation_points_by_strack = []
            corrected_unseprated = []
        return corrected_stagnation_points_by_holzbecher, corrected_stagnation_points_by_strack, corrected_unseprated

    def calculation_for_discharge_at_stagnation_points(self):
        corrected_stagnation_points_by_holzbecher, corrected_stagnation_points_by_strack, corrected_unseprated = self.calculation_for_stagnation_points_formula_by_holzbecher_and_strack()
        if len(corrected_stagnation_points_by_holzbecher) !=0 and len(corrected_stagnation_points_by_strack) !=0 and len(corrected_unseprated) !=0 :
            Stagnation_points_holzbecher = np.fromiter(corrected_stagnation_points_by_holzbecher, dtype=float)
            Stagnation_points_strack = np.fromiter(corrected_stagnation_points_by_strack, dtype=float)
            S_P_Combine = np.fromiter(corrected_unseprated, dtype=complex)

            S_P_by_me_1 = S_P_Combine.real[0] *1j
            S_P_by_me_2 = S_P_Combine.real[1] *1j
            print('sp by me 1',S_P_by_me_1)
            print('sp by me 2',S_P_by_me_2)

            y1_holzbecher = (Stagnation_points_holzbecher[0]) *1j
            y2_holzbecher = (Stagnation_points_holzbecher[1]) *1j
            print('y1 holzebecher',y1_holzbecher)
            print('y2 holzbecher', y2_holzbecher)

            y1_strack = (Stagnation_points_strack[0]) *1j
            y2_strack = (Stagnation_points_strack[1]) *1j
            print('y1strack',y1_strack)
            vectors_data_river = Vectors_river_complex()
            xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z = vectors_data_river.region_boundaries_river_complex()

            data_reading = Potentials_River_complex()
            baseFlowX, _, h0, H, k, _, pumping_array, Zw , alpha, por = data_reading.potentials_data_river_complex()

            #To Correct the Image in Complex Plane
            well_images_behind_river = Zw * -1
            corrected_well_image = (np.conjugate(well_images_behind_river))
                
            capture_lentgh_holzbecher = np.abs(Stagnation_points_holzbecher[1]) - (Stagnation_points_holzbecher[0])
            capture_lentgh_strack = np.abs(Stagnation_points_strack[1]) - (Stagnation_points_strack[0])
            combine_capture_length_by_unseprated_formula = np.abs(S_P_Combine[1] - S_P_Combine[0])
            print('This is Capture Lentgh By Holzbecher',capture_lentgh_holzbecher)
            print('This is Capture Lentgh By Strack',capture_lentgh_strack)
            print('This is Combine Capture Lentgh By Unseprated Strack Formula',combine_capture_length_by_unseprated_formula)

            Holzbecher_discharge_passing_through_stagnation_point_1 = (-baseFlowX * y1_holzbecher) - (pumping_array[0]/(2*np.pi)) * ((np.log(y1_holzbecher-Zw[0])) - (np.log(y1_holzbecher-corrected_well_image[0])))

            Holzbecher_discharge_passing_through_stagnation_point_2 = (-baseFlowX * y2_holzbecher) - (pumping_array[0]/(2*np.pi)) * ((np.log(y2_holzbecher-Zw[0])) - (np.log(y2_holzbecher-corrected_well_image[0])))        

            Strack_discharge_passing_through_stagnation_point_1 = (-baseFlowX * y1_strack) + (pumping_array[0]/(2*np.pi)) * ((np.log(y1_strack-Zw[0])) - (np.log(y1_strack-corrected_well_image[0])))

            Strack_discharge_passing_through_stagnation_point_2 = (-baseFlowX * y2_strack) + (pumping_array[0]/(2*np.pi)) * ((np.log(y2_strack-Zw[0])) - (np.log(y2_strack-corrected_well_image[0])))     


            # This is Discharge at stagnation points which we got from combine formula calculating only at the stagnation points of Y axis where x= 0
            Ghilman_discharge_passing_through_stagnation_point_1 = (-baseFlowX * (S_P_by_me_1)) + (pumping_array[0]/(2*np.pi)) * ((np.log((S_P_by_me_1)-Zw[0])) - (np.log((S_P_by_me_1)-corrected_well_image[0])))

            Ghilman_discharge_passing_through_stagnation_point_2 = (-baseFlowX * (S_P_by_me_2)) + (pumping_array[0]/(2*np.pi)) * ((np.log((S_P_by_me_2)-Zw[0])) - (np.log((S_P_by_me_2)-corrected_well_image[0])))     

            print('This is discharge passing through stagnation point 1 by Holzbecher',Holzbecher_discharge_passing_through_stagnation_point_1)
            print('This is discharge passing through stagnation point 2 by Holzbecher',Holzbecher_discharge_passing_through_stagnation_point_2)

            print('This is discharge passing through stagnation point 1 by Strack',Strack_discharge_passing_through_stagnation_point_1)
            print('This is discharge passing through stagnation point 2 by Strack',Strack_discharge_passing_through_stagnation_point_2)

            print('This is Discharge passing through at Stag 1 by actual formula', Ghilman_discharge_passing_through_stagnation_point_1)
            print('This is Discharge passing through at Stag 2 by actual formula', Ghilman_discharge_passing_through_stagnation_point_2)
            bank_filtrate_holzbecher = Holzbecher_discharge_passing_through_stagnation_point_2.imag - Holzbecher_discharge_passing_through_stagnation_point_1.imag + pumping_array[0]
            bank_filtrate_strack =  Strack_discharge_passing_through_stagnation_point_2.imag - Strack_discharge_passing_through_stagnation_point_1.imag + pumping_array[0]
            bank_filtrate_by_me = Ghilman_discharge_passing_through_stagnation_point_2.imag - Ghilman_discharge_passing_through_stagnation_point_1.imag + pumping_array[0]

            print('This is Bank Filtrate Portion By Holzbecher',bank_filtrate_holzbecher)
            print('This is Bank Filtrate Portion By Strack',bank_filtrate_strack)
            print('This is Bank Filtrate Portion By Ghilman',bank_filtrate_by_me)
            contribution_ratio_holzbecher = np.round((bank_filtrate_holzbecher / pumping_array[0]),decimals=3)
            contribution_ratio_strack = np.round((bank_filtrate_strack / pumping_array[0]),decimals=3)
            contribution_ratio_ghilman = np.round((bank_filtrate_by_me / pumping_array[0]),decimals=3)


            print('This is Contribution Ratio By Holzbecher',contribution_ratio_holzbecher)
            print('This is Contribution Ratio By Strack',contribution_ratio_strack)
            print('This is Contribution Ratio By Ghilman',contribution_ratio_ghilman)

        # Here I will Try Stagnation Points in Combination Together without sepration
        else:
            exit


        return None


        
# class working_for_lake_plus_well_No_River:
#     def __init__(self): 
#         pass  
        

#     def stagnation_points_lake_plus_well(self):

        
#         vectors_data_lake = Vectors_Lake_Complex()
#         baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vectors_data_lake.vectors_data_lake_complex()

#         potential_data_lake = Potentials_Lake_Complex()
#         R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake = potential_data_lake.potentials_data_lake_complex()
#         d_well_distance = np.abs((Z_Lake.real[0] - Zw.real[0]) / 2 )    
#         R_infinity =  np.array(R_Lake[0]**2 / 0.4 )
#         # Single Well Image at center of Lake    
#         image_of_wells = (1 * Z_Lake[0]) + 0.4
        
#         # Stagnation Point For Lake Symbol Creation
#         ZsL= sp.symbols('ZsL')
#         ZsW= sp.symbols('ZsW')

#         Eq_lake_stag_part_01 = baseFlowX * (1 + ((R_Lake[0]**2) / (ZsL-Z_Lake[0])**2))

#         Eq_lake_stag_part_02 = - (Q_lake[0]/(2*np.pi)) * (1 / ((ZsL-Z_Lake[0])*1))

#         Eq_lake_stag_part_03 = ((((pumping_array[0]) / (2*np.pi)) * ((1/(ZsL-image_of_wells)) - (1/(ZsL-Zw[0]))) * (R_Lake[0]/np.absolute(Zw[0]))) )

#         Eq_lake_plus_well04_new = ( ( pumping_array[0] * R_Lake[0] ) / (2*np.pi * np.abs(Zw[0])) ) * ((1+(R_Lake[0]**2)/(np.conjugate(Zw[0])**2)) - (1- (ZsL-Zw[0])))  
#         EQUATION_FOR_LAKE_AND_UNIFORM_FLOW_FIELD = sp.Eq((Eq_lake_stag_part_01 + Eq_lake_stag_part_02 + 0 ), 0 )
#         SP_Lake = sp.solve(EQUATION_FOR_LAKE_AND_UNIFORM_FLOW_FIELD)
        
#         S_P_UF_LAKE = np.fromiter(SP_Lake, dtype=complex)

#         S_P_UF_LAKE_UPDATED= S_P_UF_LAKE*1

#         print('These are stagnation points by series of equations',S_P_UF_LAKE_UPDATED)

#         # Calculating Discharge At The Stagation Point of Lakes

#         discharge_at_lake_stag_01 = -baseFlowX *((S_P_UF_LAKE_UPDATED[0]-Z_Lake[0])-((R_Lake[0]**2)/(S_P_UF_LAKE_UPDATED[0]-Z_Lake[0]))) + ((Q_lake[0]/(2*np.pi)) * (np.log((S_P_UF_LAKE_UPDATED[0]-Z_Lake[0])/R_infinity))) #+ ((pumping_array[0]/(2*np.pi)) * ((np.log((S_P_UF_LAKE_UPDATED[0] -Zw[0])/(S_P_UF_LAKE_UPDATED[0]-image_of_wells)*(R_Lake[0]/np.absolute(Zw[0]))))))
#         discharge_at_lake_stag_02 = -baseFlowX *((S_P_UF_LAKE_UPDATED[1]-Z_Lake[0])-((R_Lake[0]**2)/(S_P_UF_LAKE_UPDATED[1]-Z_Lake[0]))) + ((Q_lake[0]/(2*np.pi)) * (np.log((S_P_UF_LAKE_UPDATED[1]-Z_Lake[0])/R_infinity))) #+ ((pumping_array[0]/(2*np.pi)) * ((np.log((S_P_UF_LAKE_UPDATED[1] -Zw[0])/(S_P_UF_LAKE_UPDATED[1]-image_of_wells)*(R_Lake[0]/np.absolute(Zw[0]))))))
#         # discharge_at_lake_stag_03 = -baseFlowX *((S_P_UF_LAKE_UPDATED[2]-Z_Lake[0])-((R_Lake[0]**2)/(S_P_UF_LAKE_UPDATED[2]-Z_Lake[0]))) + ((Q_lake[0]/(2*np.pi)) * (np.log((S_P_UF_LAKE_UPDATED[2]-Z_Lake[0])/R_infinity))) + ((pumping_array[0]/(2*np.pi)) * ((np.log((S_P_UF_LAKE_UPDATED[2] -Zw[0])/(S_P_UF_LAKE_UPDATED[2]-image_of_wells)*(R_Lake[0]/np.absolute(Zw[0]))))))

#         print('This is discharge at stagnation point 1',discharge_at_lake_stag_01)
#         print('This is discharge at stagnation point 2',discharge_at_lake_stag_02)
#         # print('This is discharge at stagnation point 3',(disc3harge_at_lake_stag_03))
#         if x_array[0] < x_Lake[0]:
#             # EQ For Well and Uniform Flow Field
#             Eq_lake_stag_part_03 = sp.Eq(( baseFlowX + (((pumping_array[0]) / (2*np.pi)) * ((1/(ZsW-image_of_wells)) - (1/(ZsW-Zw[0]))) * (1)) ),0 )
#             # Eq_lake_plus_well04_new =     sp.Eq(( baseFlowX + ( ( pumping_array[0] * R_Lake[0] ) / (2*np.pi * np.abs(Zw[0])) ) * ((1+(R_Lake[0]**2)/(np.conjugate(Zw[0])**2)) - (1- (ZsL-Zw[0]))) ),0)
#             solution_well_and_unform_flow = sp.solve(Eq_lake_stag_part_03)
#             S_P_UF_WELL = np.fromiter(solution_well_and_unform_flow, dtype=complex)
#             S_P_UF_WELL_UPDATED = S_P_UF_WELL + 0
#             print('Thsese are Stagnation Points Due To Well',S_P_UF_WELL)
#             print('Thsese are Stagnation Points updated Due To Well',S_P_UF_WELL_UPDATED)


#             # Discharge Passing at Stagnation Points Due To Well
#             discharge_at_lake_stag_01_due_to_well=((pumping_array[0]/(2*np.pi)) * ((np.log((S_P_UF_LAKE_UPDATED[0] -Zw[0])/(S_P_UF_LAKE_UPDATED[0]-image_of_wells)*(R_Lake[0]/np.absolute(Zw[0]))))))         #R_Lake[0]/np.absolute(Zw[0])
#             discharge_at_lake_stag_02_due_to_well=((pumping_array[0]/(2*np.pi)) * ((np.log((S_P_UF_LAKE_UPDATED[1] -Zw[0])/(S_P_UF_LAKE_UPDATED[1]-image_of_wells)*(R_Lake[0]/np.absolute(Zw[0]))))))         #R_Lake[0]/np.absolute(Zw[0]
#             print('Discharge at Stag 1 Because of Well Ghilman',discharge_at_lake_stag_01_due_to_well)
#             print('Discharge at Stag 2 Because of Well Ghilman',discharge_at_lake_stag_02_due_to_well)

#             lake_filtrate =  discharge_at_lake_stag_02_due_to_well.imag - discharge_at_lake_stag_01_due_to_well.imag + pumping_array[0] + discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag - Q_lake[0] #+ discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag
#             print(lake_filtrate,'m\u00b3/d')

#             lake_contriubtion_ratio = np.round(lake_filtrate/(pumping_array[0]),decimals=3)
#             print(lake_contriubtion_ratio)
#         else:
#             # EQ For Well and Uniform Flow Field
#             Eq_lake_stag_part_03 = sp.Eq(( baseFlowX + (((pumping_array[0]) / (2*np.pi)) * ((1/(ZsW-Zw[0])) - (1/(ZsW-image_of_wells))) * (1)) ),0 )
#             # Eq_lake_plus_well04_new =     sp.Eq(( baseFlowX + ( ( pumping_array[0] * R_Lake[0] ) / (2*np.pi * np.abs(Zw[0])) ) * ((1+(R_Lake[0]**2)/(np.conjugate(Zw[0])**2)) - (1- (ZsL-Zw[0]))) ),0)
#             solution_well_and_unform_flow = sp.solve(Eq_lake_stag_part_03)
#             S_P_UF_WELL = np.fromiter(solution_well_and_unform_flow, dtype=complex)
#             S_P_UF_WELL_UPDATED = S_P_UF_WELL + 0
#             print('Thsese are Stagnation Points Due To Well',S_P_UF_WELL)
#             print('Thsese are Stagnation Points updated Due To Well',S_P_UF_WELL_UPDATED)


#             # Discharge Passing at Stagnation Points Due To Well
#             discharge_at_lake_stag_01_due_to_well=((pumping_array[0]/(2*np.pi)) * ((np.log((S_P_UF_WELL_UPDATED[0] -Zw[0])/(S_P_UF_WELL_UPDATED[0]-image_of_wells)*(R_Lake[0]/np.absolute(Zw[0]))))))       #R_Lake[0]/np.absolute(Zw[0])
#             discharge_at_lake_stag_02_due_to_well=((pumping_array[0]/(2*np.pi)) * ((np.log((S_P_UF_WELL_UPDATED[1] -Zw[0])/(S_P_UF_WELL_UPDATED[1]-image_of_wells)*(R_Lake[0]/np.absolute(Zw[0]))))))        #R_Lake[0]/np.absolute(Zw[0]
#             print('Discharge at Stag 1 Because of Well',discharge_at_lake_stag_01_due_to_well)
#             print('Discharge at Stag 2 Because of Well',discharge_at_lake_stag_02_due_to_well)

#             lake_filtrate =  discharge_at_lake_stag_01_due_to_well.imag - discharge_at_lake_stag_02_due_to_well.imag + pumping_array[0] + discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag - Q_lake[0]   #+ discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag
#             print('This is Lake Filtrate',lake_filtrate,'m\u00b3/d')

#             lake_contriubtion_ratio = np.round(lake_filtrate/(pumping_array[0]),decimals=3)
#             print(lake_contriubtion_ratio)
            
#         return None

        
class working_for_lake_plus_well_with_River:
    def __init__(self): 
        pass  
        

    def stagnation_points_lake_and_well_plus_river(self):

        
        vectors_data_lake = Vectors_Lake_River_Complex()
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vectors_data_lake.vectors_data_lake_river_complex()

        potential_data_lake = Potentials_Lake_River_Complex()
        R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake = potential_data_lake.potentials_data_lake_river_complex()

        R_infinity =  np.array(R_Lake[0]**2 / 0.4 )
        # Single Well Image behind river
        trying_well_images = Zw[0] * -1
        actual_well_images = np.conjugate(trying_well_images)
        
        # Stagnation Point For Lake Symbol Creation
        ZsL= sp.symbols('ZsL')
        ZsW= sp.symbols('ZsW')
        q_o = (Q_lake[0])/(2*np.pi*baseFlowX*R_Lake[0]) 
        q_o_by_2 = q_o/2
        print('this is q0by2',q_o_by_2)


        
        if q_o_by_2 >= -1 and q_o_by_2 <= 1:
            Eq_lake_stag_part_01 =   baseFlowX * (1 + ((R_Lake[0]**2) / (ZsL-Z_Lake[0])**2))
            Eq_lake_stag_part_02 =  -(Q_lake[0]/(2*np.pi)) * (1 / ((ZsL-Z_Lake[0])*1))

            # Eq_lake_stag_part_03_R = -(pumping_array[0]) / ((2*np.pi)) * ((1/(ZsL-actual_well_images)) - (1/(ZsL-Zw[0]))) * (R_Lake[0]/np.absolute(Zw[0]))
        
            EQUATION_FOR_LALE_AND_UNIFORM_FLOW_FIELD = sp.Eq((Eq_lake_stag_part_01 + Eq_lake_stag_part_02 + 0), 0 )
            SP_Lake = sp.solve(EQUATION_FOR_LALE_AND_UNIFORM_FLOW_FIELD)
            S_P_OF_LAKE = np.fromiter(SP_Lake, dtype=complex)
            print('Stagnation Points of Lake are(River is present):',S_P_OF_LAKE)
            if x_array[0] < x_Lake[0]:
                # Calculating Discharge At The Stagation Point Which is Always I guess at +- iR
                discharge_at_lake_stag_01 = -baseFlowX* S_P_OF_LAKE[0]  - ( (pumping_array[0]/(2*np.pi)) * np.log(( S_P_OF_LAKE[0]-Zw[0])/(S_P_OF_LAKE[0] - actual_well_images )) )
                discharge_at_lake_stag_02 = -baseFlowX* S_P_OF_LAKE[1]  - ( (pumping_array[0]/(2*np.pi)) * np.log(( S_P_OF_LAKE[1]-Zw[0])/(S_P_OF_LAKE[1] - actual_well_images )) )
                
                lake_filtrate = discharge_at_lake_stag_01.imag - discharge_at_lake_stag_02.imag   + np.abs(Q_lake[0])
                lake_contribution_ratio = np.round(lake_filtrate/(pumping_array[0]), decimals=3) 

            else:

                lake_filtrate = -1
                        
            if lake_filtrate < 0:
                print('Can Not Calculate Lake Filtrate')
            else:
                print('Discharge Passing at First Stagnation Point of Lake :',discharge_at_lake_stag_01.imag,' m\u00b3/d')
                print('Discharge Passing at Second Stagnation Point of Lake :',discharge_at_lake_stag_02.imag,' m\u00b3/d')

                print('This is lake Filtrate whenn River is Present',lake_filtrate, 'm\u00b3/d')
                print('This is Lake Contribution Ratio',lake_contribution_ratio )

        else:
            print("Can Not Calculate Stagnation Point for Lake, Please Enter Realistic Values")
        

        # EQ For Well and Uniform Flow Field
        Eq_lake_stag_part_03 = sp.Eq( ( baseFlowX - ((pumping_array[0]/(2*np.pi)*(1)) * (1/(ZsW-Zw.real[0]))) + ((pumping_array[0]/(2*np.pi)) * (1/(ZsW-actual_well_images.real)))),0)
        solution_well_and_unform_flow = sp.solve(Eq_lake_stag_part_03)
        S_P_UF_WELL = np.fromiter(solution_well_and_unform_flow, dtype=complex)
        print('This is sp of well combine without correction',S_P_UF_WELL)
        SP_OF_WELL_1 = (S_P_UF_WELL[0].real + Zw[0].imag) *1j 
        SP_OF_WELL_2 = (S_P_UF_WELL[1].real + Zw[0].imag) *1j
        print('Thsese are Stagnation Points Due To Well and River',SP_OF_WELL_1,SP_OF_WELL_2)

        # Discharge Passing at Stagnation Points Due To Well

        discharge_at_stag_01_due_to_well= -baseFlowX * ((SP_OF_WELL_1-Z_Lake[0])-((R_Lake[0]**2)/(SP_OF_WELL_1-Z_Lake[0]))) + ((Q_lake[0]/(2*np.pi)) * (np.log((SP_OF_WELL_1-Z_Lake[0])/R_infinity))) - ( (pumping_array[0]/(2*np.pi)) * np.log((SP_OF_WELL_1-Zw[0])/(SP_OF_WELL_1 -actual_well_images )) ) 
        discharge_at_stag_02_due_to_well= -baseFlowX * ((SP_OF_WELL_2-Z_Lake[0])-((R_Lake[0]**2)/(SP_OF_WELL_2-Z_Lake[0]))) + ((Q_lake[0]/(2*np.pi)) * (np.log((SP_OF_WELL_2-Z_Lake[0])/R_infinity))) - ( (pumping_array[0]/(2*np.pi)) * np.log((SP_OF_WELL_2-Zw[0])/(SP_OF_WELL_2 -actual_well_images )) )
        
        print('Discharge at Stag 1 Because of Well at River',discharge_at_stag_01_due_to_well)
        print('Discharge at Stag 2 Because of Well at River',discharge_at_stag_02_due_to_well)

        river_filtrate =  discharge_at_stag_02_due_to_well.imag - discharge_at_stag_01_due_to_well.imag + pumping_array[0] 

        river_contriubtion_ratio = np.round((river_filtrate/(pumping_array[0])), decimals=3)
        river_capture_lentgh = np.round(np.abs(SP_OF_WELL_2.imag - SP_OF_WELL_1.imag),decimals=3 )
        
        if river_filtrate >0 :
            print('This is River Filtrate',river_filtrate,'m\u00b3/d')
            print('This is River Capture Lentgh', river_capture_lentgh, 'm')
            print('This is River Contributution Ratio',river_contriubtion_ratio)
        else:
            print('Can Not Calculate River Filtrate or Contribution Ratio ! Please Increase Pumping Rate or Move Well Towards River')    

        return None



class working_for_lake_plus_well_No_River_GD:
    def __init__(self): 
        pass  
        

    def stagnation_points_lake_and_well_noriver_GD(self):

        
        vectors_data_lake = Vectors_Lake_Complex()
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vectors_data_lake.vectors_data_lake_complex()

        potential_data_lake = Potentials_Lake_Complex()
        R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake = potential_data_lake.potentials_data_lake_complex()
        d_well_distance = np.abs((Z_Lake.real[0] - Zw.real[0]) / 2 )    
        R_infinity =  np.array(R_Lake[0]**2 / 0.4 )
        # Single Well Image at center of Lake    
        image_of_wells = (1 * Z_Lake[0]) + 0.4
        
        # Stagnation Point For Lake Symbol Creation
        ZsL= sp.symbols('ZsL')
        ZsW= sp.symbols('ZsW')
        q_o = (Q_lake[0])/(2*np.pi*baseFlowX*R_Lake[0]) 
        q_o_by_2 = q_o/2
        print('This is qoby2',q_o_by_2)


        Eq_lake_stag_part_01 = baseFlowX * (1 + ((R_Lake[0]**2) / (ZsL-Z_Lake[0])**2))

        # Eq_lake_stag_part_01 = 0 
        Eq_lake_stag_part_02 = - (Q_lake[0]/(2*np.pi)) * (1 / ((ZsL-Z_Lake[0])*1))
        # Eq_lake_stag_part_03 = -((((pumping_array[0]) / (2*np.pi)) * ((1/(ZsL-image_of_wells)) - (1/(ZsL-Zw[0]))) * (R_Lake[0]/np.absolute(Zw[0]))) )
        EQUATION_FOR_LAKE_AND_UNIFORM_FLOW_FIELD = sp.Eq((Eq_lake_stag_part_01 + Eq_lake_stag_part_02 + 0 ), 0 )
        SP_Lake = sp.solve(EQUATION_FOR_LAKE_AND_UNIFORM_FLOW_FIELD)
    
        S_P_UF_LAKE = np.fromiter(SP_Lake, dtype=complex)
        print('These are Stagnation Points Due To Lake (No River)', S_P_UF_LAKE)
 

            # Calculating Discharge At The Stagation Point of Lakes
        if q_o_by_2 >= -1 and q_o_by_2 <= 1:
            # print(x_array)
            # print(x_Lake)
            # print('dsafadsfsadfasdf')
            if x_array[0] < x_Lake[0]:
                                               
                discharge_at_stag_01_due_to_well=  (pumping_array[0]/(2*np.pi)) * ((np.log(((S_P_UF_LAKE[0] -Zw[0])/(S_P_UF_LAKE[0]-image_of_wells))*(R_Lake[0]/np.absolute(Zw[0]))))) -((Q_lake[0]/(2*np.pi)) * (np.log((S_P_UF_LAKE[0]-Z_Lake[0])/R_infinity)))      -baseFlowX * ((S_P_UF_LAKE[0]-Z_Lake[0])-((R_Lake[0]**2)/(S_P_UF_LAKE[0]-Z_Lake[0])))                 #-baseFlowX*S_P_UF_LAKE[0]
                discharge_at_stag_02_due_to_well=  (pumping_array[0]/(2*np.pi)) * ((np.log(((S_P_UF_LAKE[1] -Zw[0])/(S_P_UF_LAKE[1]-image_of_wells))*(R_Lake[0]/np.absolute(Zw[0]))))) -((Q_lake[0]/(2*np.pi)) * (np.log((S_P_UF_LAKE[1]-Z_Lake[0])/R_infinity)))      -baseFlowX * ((S_P_UF_LAKE[1]-Z_Lake[0])-((R_Lake[0]**2)/(S_P_UF_LAKE[1]-Z_Lake[0])))                 #-baseFlowX*S_P_UF_LAKE[1]
                print('This is discharge at first stag point of lake', discharge_at_stag_01_due_to_well)
                print('This is discharge at second stag point of lake', discharge_at_stag_02_due_to_well)
                Lake_filtrate = discharge_at_stag_02_due_to_well.imag - discharge_at_stag_01_due_to_well.imag  + pumping_array[0]
            else:
        
                # discharge_at_stag_01_due_to_well=  (pumping_array[0]/(2*np.pi)) * ((np.log(((S_P_UF_LAKE[0] -image_of_wells)/(S_P_UF_LAKE[0]-Zw[0]))*(R_Lake[0]/np.absolute(Zw[0]))))) +((Q_lake[0]/(2*np.pi)) * (np.log((S_P_UF_LAKE[0]-Z_Lake[0])/R_infinity)))  -baseFlowX * ((S_P_UF_LAKE[0]-Z_Lake[0])-((R_Lake[0]**2)/(S_P_UF_LAKE[0]-Z_Lake[0])))#-baseFlowX*S_P_UF_LAKE[0]
                # discharge_at_stag_02_due_to_well=  (pumping_array[0]/(2*np.pi)) * ((np.log(((S_P_UF_LAKE[1] -image_of_wells)/(S_P_UF_LAKE[1]-Zw[0]))*(R_Lake[0]/np.absolute(Zw[0]))))) +((Q_lake[0]/(2*np.pi)) * (np.log((S_P_UF_LAKE[1]-Z_Lake[0])/R_infinity)))  -baseFlowX * ((S_P_UF_LAKE[1]-Z_Lake[0])-((R_Lake[0]**2)/(S_P_UF_LAKE[1]-Z_Lake[0])))   #-baseFlowX*S_P_UF_LAKE[1]        
                print('Warning! x coordinate of well should be less than  x coordinates of Lake')
                Lake_filtrate = -1
            # if x_array < x_Lake and y_array >= y_Lake:
            #     Lake_filtrate = discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag + discharge_at_stag_02_due_to_well.imag - discharge_at_stag_01_due_to_well.imag + pumping_array[0] - Q_lake[0]
            #     lake_CR= Lake_filtrate/pumping_array[0]

            # elif x_array < x_Lake and y_array < y_Lake:
            #     Lake_filtrate = discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag + discharge_at_stag_01_due_to_well.imag - discharge_at_stag_02_due_to_well.imag + pumping_array[0] - Q_lake[0]
            #     lake_CR= Lake_filtrate/pumping_array[0]
            
            # elif x_array > x_Lake and y_array > y_Lake:
            #     Lake_filtrate = discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag + discharge_at_stag_01_due_to_well.imag - discharge_at_stag_02_due_to_well.imag + pumping_array[0] - Q_lake[0]
            #     lake_CR= Lake_filtrate/pumping_array[0]

            # # This is Fourth Case x_array > x_Lake and y_array <= y_Lake
            # elif x_array > x_Lake and y_array <= y_Lake:
            #     Lake_filtrate = discharge_at_lake_stag_02.imag - discharge_at_lake_stag_01.imag + discharge_at_stag_02_due_to_well.imag - discharge_at_stag_01_due_to_well.imag + pumping_array[0] - Q_lake[0]
            #     lake_CR= Lake_filtrate/pumping_array[0]

            
            if Lake_filtrate < 0 :
                print("Can Not Calculate Stagnation Point for Lake, Please Enter Realistic Values")

            elif Lake_filtrate > 1:
                lake_CR = (((Lake_filtrate/pumping_array[0])))    
                print('This is Contribution Ration Of Lake',lake_CR*100,'%')
            # print('Discharge Passing at Stag 1 Because of Well Pumping hehehe',discharge_at_stag_01_due_to_well)
            # print('Discharge Passing at Stag 2 Because of Well Pumping hehehe',discharge_at_stag_02_due_to_well)
                    

            # if lake_CR > 1:
            #     lake_CR =1
            # else:
                
        else:
            print("Can Not Calculate Stagnation Point for Lake, Please Enter Realistic Values")

            
        return None


