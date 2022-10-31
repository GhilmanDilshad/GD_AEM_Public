from email.mime import base
from symtable import Symbol
import pandas as pd
import numpy as np
import sympy as sp
from utilities import read_aquifer_xlsx, read_wells_xlsx
import os
from river.river_complex import Potentials_River_complex
from river.river_complex import Vectors_river_complex

from lake.lake_with_river_complex import Vectors_Lake_River_Complex
from lake.lake_with_river_complex import Potentials_Lake_River_Complex


class contribution_ratio_river_and_well:
    def __init__(self): 
        pass  

    def limit_of_stagnations_points(self):

        data_reading = Potentials_River_complex()
        baseFlowX, _, h0, H, k, _, pumping_array, Zw , alpha, por = data_reading.potentials_data_river_complex()

        condition_for_stagnation_point_to_exist = []
        
        limit = ( pumping_array[0] ) / (np.pi * (Zw.real[0]) * baseFlowX)
        condition_for_stagnation_point_to_exist.append(limit)

        print('This is condition for stagnation point to exist ',condition_for_stagnation_point_to_exist)
        # print('This is type of Zw',type(Zw))
        # print('This is Zw imag value',Zw[0])
        return condition_for_stagnation_point_to_exist

    
    def calculation_for_stagnation_points(self): 
        data_reading = Potentials_River_complex()
        baseFlowX, _, h0, H, k, _, pumping_array, Zw , alpha, por = data_reading.potentials_data_river_complex()
        data_reading_vector = Vectors_river_complex()
        
        stagnation_points = []
        limit = self.limit_of_stagnations_points()
        x = np.array(limit)
        y = sp.symbols('y')
        well_distance = np.array(Zw.imag[0])

        Zs = sp.Symbol('Zs')

        well_images_behind_river = Zw * -1
        corrected_well_image1 = np.conjugate(well_images_behind_river)
        stagnation_point_by_complex = baseFlowX + (pumping_array[0]/(2*np.pi)) * ((1/(Zs-corrected_well_image1)) - (1/(Zs-Zw)))
        print('This is stagnaton point by complex just checking', sp.solve(stagnation_point_by_complex))

        if x >= 1:
            print("Stagnation Points Exists")

            Equation_for_calc_stag_points = sp.Eq(y**2, -((Zw.real[0])**2) +  (((Zw.real[0])*pumping_array[0])/(np.pi*baseFlowX)))
            two_stagnation_points = sp.solve(Equation_for_calc_stag_points)

            corrected_stagnation_points = (two_stagnation_points + well_distance)
            print('Corrected Stagnation Points', corrected_stagnation_points)

        else:
            corrected_stagnation_points = []
            print("No Stagnation Points Exist")

        return corrected_stagnation_points

    def calculation_capture_lentgh(self):
        
        corrected_stagnation_points1 = self.calculation_for_stagnation_points()

        if len(corrected_stagnation_points1) != 0:
            first_point = np.round(float((corrected_stagnation_points1[0])), decimals=1)
            second_point = np.round(float((corrected_stagnation_points1[1])), decimals=1)

            capture_length = np.abs((second_point)-(first_point))
            print('This is Capture Length', capture_length)
            print('These are two points', first_point,second_point)

        else:
            capture_length=0
            first_point = 0 
            second_point = 0 
        return capture_length, first_point, second_point
        
#is say uper tak ka sab sahi hai ab agay ka dekhty hai 

    def calculation_for_contribution_by_stream_fxn(self):
        vectors_data_river = Vectors_river_complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z = vectors_data_river.region_boundaries_river_complex()

        data_reading = Potentials_River_complex()
        baseFlowX, _, h0, H, k, _, pumping_array, Zw , alpha, por = data_reading.potentials_data_river_complex()
        corrected_stagnation_points1= self.calculation_for_stagnation_points()
        capture_length, first_point, second_point = self.calculation_capture_lentgh()

        if capture_length != 0 and first_point != 0 and second_point != 0:

            #Conversion of Ystag to imaginary part
            imaginary_first_point = first_point*(+1j)
            imaginary_second_point = second_point*(+1j)

            #To Correct the Image in Complex Plane
            well_images_behind_river = Zw * -1
            corrected_well_image = (np.conjugate(well_images_behind_river))

            Z_stag_1 = imaginary_first_point
            Z_stag_2 = imaginary_second_point

            complex_stag_1 = (-baseFlowX * Z_stag_1) - (pumping_array/(2*np.pi)) * ((np.log(Z_stag_1-Zw[0])) - (np.log(Z_stag_1-corrected_well_image[0])))

            complex_stag_2 = (-baseFlowX * Z_stag_2) - (pumping_array/(2*np.pi)) * ((np.log(Z_stag_2-Zw[0])) - (np.log(Z_stag_2-corrected_well_image[0])))        

            Q_bf = complex_stag_2.imag - complex_stag_1.imag + pumping_array
            print(complex_stag_1)
            print(complex_stag_2)
            print('Bank Filtrate',Q_bf,'m\u00b3/d ')

            ratio = np.round((Q_bf / pumping_array),decimals=3)

            print('This is Contribution Ratio', ratio)
            print('This is Contribution Ratio Percentage',ratio*100,'%') 


        # def tracing_streamlines(self, delta_s=0.1,minimum_distance_travelled=0.1):

        #     data_reading = Potentials_River_complex()
        #     baseFlowX, _, h0, H, k, _, pumping_array, Zw , alpha, por = data_reading.potentials_data_river_complex()

        #     Z_stag_1, Z_stag_2, Q_bf, ratio = self.calculation_for_stagnation_points()



        # def Qx_discharge_vector():
        return None


# class contribution_ratio_lake_and_well_GD:
#     def __init__(self): 
#         pass  


    def stagnation_points_w_r_t_lake(self):

        lake_data = Potentials_Lake_River_Complex()
        R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake = lake_data.potentials_data_lake_river_complex()
        
        uniform_flow_data = Vectors_Lake_River_Complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis = uniform_flow_data.region_boundaries_lake_river_complex()
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = uniform_flow_data.vectors_data_lake_river_complex()
        
        trying_well_images = Zw * -1
        actual_well_images = np.conjugate(trying_well_images)
        # R_infinity =  np.array(R_Lake)**2 / (0.4 )
        Zs = sp.Symbol('Zs')


        part1_uniform_flow =  baseFlowX  * ( 1 + (( R_Lake[0] **2 ) / (Zs - Z_Lake[0]) ** 2 ))
        part2_due_to_lake  =  Q_lake[0] /  ((2*np.pi) * ( Zs - Z_Lake[0]) )
        part3_due_to_well  =  (pumping_array[0] / (2*np.pi)) * ( (1/ ((Zs - actual_well_images[0]))) - (1/ (Zs - Zw[0]) ))


        Equation_for_stagnation_points_lake_alone = part1_uniform_flow  - part2_due_to_lake
        
        Stagnation_Points_Lake = sp.solve(Equation_for_stagnation_points_lake_alone)
        # print('These are stagnation points for lake', Stagnation_Points_Lake)
        Stagnation_Points_Lake_array = np.fromiter(Stagnation_Points_Lake, dtype=complex)

        print('These are stagnation points for lake', Stagnation_Points_Lake_array)

        
        discharge_at_first_stagnation_point_lake =  -1 * baseFlowX * ((Stagnation_Points_Lake_array[0]-Z_Lake[0])-((R_Lake[0]**2)/(Stagnation_Points_Lake_array[0]-Z_Lake[0]))) + (Q_lake[0]/(2*np.pi)) * np.log((Stagnation_Points_Lake_array[0]-Z_Lake[0])/R_Lake[0])
        discharge_at_second_stagnation_point_lake = -1 * baseFlowX * ((Stagnation_Points_Lake_array[1]-Z_Lake[0])-((R_Lake[0]**2)/(Stagnation_Points_Lake_array[1]-Z_Lake[0]))) + (Q_lake[0]/(2*np.pi)) * np.log((Stagnation_Points_Lake_array[1]-Z_Lake[0])/R_Lake[0])
        discharge_flowing_through_lake = discharge_at_first_stagnation_point_lake - discharge_at_second_stagnation_point_lake

        print('This is discharge potential at first stagnation point', discharge_at_first_stagnation_point_lake)
        print('This is discharge potential at second stagnation point', discharge_at_second_stagnation_point_lake)
        print('This is discharge Flowing Through Lake', discharge_flowing_through_lake.imag,'m\u00b3/d' )


        return None


    # def stagnation_point_lake_river_well_combine_new_way(self):
    #     lake_data = Potentials_Lake_River_Complex()
    #     R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake = lake_data.potentials_data_lake_river_complex()
        
    #     uniform_flow_data = Vectors_Lake_River_Complex()
    #     xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis = uniform_flow_data.region_boundaries_lake_river_complex()
    #     baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = uniform_flow_data.vectors_data_lake_river_complex()
        
    #     trying_well_images = Zw * -1
    #     actual_well_images = np.conjugate(trying_well_images)
    #     R_infinity =  np.array(R_Lake[0])**2 / (0.4 )
    #     Ys = sp.Symbol('Ys')

    #     limit_of_well_stagnation_point_to_exist = pumping_array[0] / (np.pi*Zw[0]*baseFlowX)
    #     q_o =  (Q_lake[0]) / (2*np.pi*baseFlowX*R_Lake[0])


        
    #     # if q_o == 0 and limit_of_well_stagnation_point_to_exist < 1 :
    #     #     print('The Stagnation Points Exists at +-iR')
    #     #     part1_uniform_flow =  baseFlowX  * ( 1 + (( R_Lake[0] **2 ) / (Ys - Z_Lake[0]) ** 2 )) 
    #     #     part2_due_to_lake = 0
    #     #     print('The Lake is Stagnant in State of Equilibrium and The Wells Are not Pumping Strong Enough To Induce Stagnattion Points')

    #     if q_o == 0 and limit_of_well_stagnation_point_to_exist > 1 :
    #         part1_uniform_flow =  baseFlowX  * ( 1 + (( R_Lake[0] **2 ) / (Ys - Z_Lake[0]) ** 2 )) + (pumping_array[0] / (2*np.pi)) * ( (R_Lake[0]/ (np.absolute(Zw[0])*(Ys - actual_well_images[0]))) - (R_Lake[0]/ (np.absolute(Zw[0])*(Ys - Zw[0]) )))
    #         part2_due_to_lake = 0
    #         print('The Lake is Stagnant in State of Equilibrium and The Wells Are Pumping Strong Enough To Induce Stagnattion Points')

    #     if q_o/2 < 1 and limit_of_well_stagnation_point_to_exist < 1:
    #         print('There are two Stagnation Points and The Wells are not pumping enough To induce Stagnation Points')
    #         part1_uniform_flow =  baseFlowX  * ( 1 + (( R_Lake[0] **2 ) / (Ys - Z_Lake[0]) ** 2 ))
    #         part2_due_to_lake  =  -Q_lake[0] /  ((2*np.pi*R_infinity) * ( Ys - Z_Lake[0]) )
    #         print('this is part 1',part1_uniform_flow)
    #         print('this is part 2',part2_due_to_lake)
    #     if q_o/2 < 1 and limit_of_well_stagnation_point_to_exist > 1:
    #         ('There are two Stagnation Points and The Wells are pumping enough To induce Stagnation Points')
    #         part1_uniform_flow =  baseFlowX  * ( 1 + (( R_Lake[0] **2 ) / (Ys - Z_Lake[0]) ** 2 )) + (pumping_array[0] / (2*np.pi)) * ((R_Lake[0]/ (np.absolute(Zw[0])*(Ys - actual_well_images[0]))) -(R_Lake[0]/ (np.absolute(Zw[0])*(Ys - Zw[0]) )))
    #         part2_due_to_lake  =  -Q_lake[0] /  ((2*np.pi*R_infinity) * ( Ys - Z_Lake[0]) )

    #         print('This is Part 1',part1_uniform_flow)
    #         print('This is part 2',part2_due_to_lake)
    #     # else:
    #     #     print('Unable To Calculate Stagnationts With Current Scenario')
    #     #     exit

        
    #     ZS=sp.symbols('ZS')
    #     ZS_imaginary_of_well = ZS*1j
    #     Zs_of_well = sp.Eq((baseFlowX + ((pumping_array[0]/(2*np.pi)) * ((1/(ZS_imaginary_of_well-actual_well_images[0]))-(1/(ZS_imaginary_of_well-Zw[0])))) ),0)
    #     solutionofZS = sp.solve(Zs_of_well)
    #     print('this is wholesome solution in terms of complex with X=0 Holzbecher', solutionofZS)


    #     Equation_for_stagnation_points_lake_alone = part1_uniform_flow  + part2_due_to_lake
    #     Stagnation_Points_Lake = sp.solve(Equation_for_stagnation_points_lake_alone)
    #     Stagnation_Points_Lake_array = np.fromiter(Stagnation_Points_Lake, dtype=complex)       
    #     print('These are stagnation points for' , Stagnation_Points_Lake_array)