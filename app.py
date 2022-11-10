from wells.save_well_data import SaveWellData
from aquifer.save_aquifer_data import SaveAquiferData
from potentials.potentials_calc import Potentials
from vectorcalculations.vector_cords import Vectors
from plots.plots_work import Plotting
from river.river_calc import Potentials_River

from river.river_complex import Potentials_River_complex
from lake.lake_cal import Potentials_Lake
from lake.lake_complex import Potentials_Lake_Complex
from lake.lake_with_river_complex import Potentials_Lake_River_Complex
from potentials_complex.potentials_complex_calc import Potentials_normal_well_complex
from potentials_complex.complex_without_seperation import Potentials_normal_well_complex_without_sepration

from linesink.line_river_complex import Potentials_Line_River_Complex
from linesink.line_complex import Potentials_Line_Complex

from contribution_ratio.contribution_ratio_GD import contribution_ratio_river_and_well_strack_GD
from contribution_ratio.contribution_ratio_GD import working_for_lake_plus_well_with_River
from contribution_ratio.contribution_ratio_GD import working_for_lake_plus_well_No_River_GD

# 
import sys
from lake.read_and_save_data import ReadAndSaveLakeData

####################---------------------Single Or Multiple Well Without River Boundary-----------------------##########################

# potentialcalculations = Potentials()
# dischargepotentialrefrerence = potentialcalculations.phi_0()
# phibaseQox = potentialcalculations.phi_base_fxn()
# phiwell = potentialcalculations.calculation_phi_well()
# totalpotential =  potentialcalculations.dischargepotential_total_of_region()
# psi = potentialcalculations.calculation_stream_fxn()
# head = potentialcalculations.calculation_for_head()


####################---------------------To RUN Complex Numbers
# complex_normal_wells = Potentials_normal_well_complex()
# # complex_phi_total = complex_normal_wells.complex_potential_phi_total_region()
# # complex_psi_total = complex_normal_wells.complex_stream_total_region()
# # head_complex = complex_normal_wells.calculation_for_head_complex()
# potentials_data1 = complex_normal_wells.potential_data_complex()


####################---------------------To RUN With River or Method of Iamges With phi and head being Constant at Y Axis-----------------------##########################

# river = Potentials_River()
# phibaseQoxriver = river.phi_base_fxn_river()
# phi0river = river.phi_0_river()
# phiwellriver = river.calculation_phi_well_river()
# totalphiriver = river.dischargepotential_total_of_region_river()
# headwithriver = river.calculation_for_head_river()
# psiwithriver = river.calculation_stream_fxn_river()

# ####################---------------------To RUN With River_Complex or Method of Iamges With phi and head being Constant at Y Axis-----------------------##########################
# river_complex = Potentials_River_complex()
# potentials_data = river_complex.potentials_data_river_complex()
# phi_base = river_complex.phi_base_fxn_river_complex()
# phi_well = river_complex.calculation_phi_well_river_complex()

####################---------------------To RUN With Lake
# Lake_Data = Potentials_Lake()
# phibaselake = Lake_Data.phi_base_fxn_lake() 
# phi0 = Lake_Data.phi_0_lake()
# phi_of_well = Lake_Data.calculation_phi_well_lake()
# inflow_outflow = Lake_Data.calculation_inflow_outflow_lake()
# phi_total_lake = Lake_Data.dischargepotential_total_of_region_lake()
# head = Lake_Data.calculation_for_head_lake()
# psi_base_lake = Lake_Data.calculation_stream_fxn_lake()


#Complex Lake--------------------------

# Lake_Complex = Potentials_Lake_Complex()
# potentialsdatalakecomplex = Lake_Complex.potentials_data_lake_complex()
# phi_uniform_flow_field_complex = Lake_Complex.phi_uniform_flow_field_complex()
# phi_well_lake_complex = Lake_Complex.calculation_phi_well_lake_complex()
# phi_inflow_outflow = Lake_Complex.calculation_inflow_outflow_lake_complex()
# head_complex_lake = Lake_Complex.calculation_for_head_lake_complex()
# # phi_inifinity = Lake_Complex.phi_far_field_R_infinity()
# C_constant = Lake_Complex.the_constant_C()

#Lake with River Complex--------------------------

# Lake_River_Complex = Potentials_Lake_River_Complex()
# potentialsdatalakecomplex = Lake_River_Complex.potentials_data_lake_river_complex()
# phi_uniform_flow_field_complex = Lake_River_Complex.phi_uniform_flow_field_river_complex()
# phi_well_lake_complex = Lake_River_Complex.calculation_phi_well_lake_river_complex()
# phi_inflow_outflow = Lake_River_Complex.calculation_inflow_outflow_lake_river_complex()
# head_complex_lake = Lake_River_Complex.calculation_for_head_lake_river_complex()
# # phi_inifinity = Lake_Complex.phi_far_field_R_infinity_river()
# C_constant = Lake_River_Complex.the_constant_C_river()





# complex_without_sepration_try = Potentials_normal_well_complex_without_sepration()

# potentialdata = complex_without_sepration_try.potential_data_complex_without_sepration()
# phi_complex_withot_sepration = complex_without_sepration_try.phi_well_plus_baseflowphi_psi_well_plus_baseflowpsi_complex()
# head_complex_without_speration = complex_without_sepration_try.calculation_for_head_complex_without_sepration()

# line_sink
# potentialsdata = Potentials_Line_River_Complex()
# phi_well_data = potentialsdata.calculation_phi_well_line_river_complex()
# phi_inflow_outflow_line = potentialsdata.calculation_inflow_outflow_line_river_complex()

# potentiansdataline = Potentials_Line_Complex()
# phi_well_data_line = potentiansdataline.calculation_phi_well_line_complex()


# Stagnation Points

# GD_trying_stagnation_points_by_strack = contribution_ratio_river_and_well_strack_GD()
# limit_of_stagnation_Points_GD = GD_trying_stagnation_points_by_strack.limit_of_stagnations_points_by_strack()
# stagnation_points_GD = GD_trying_stagnation_points_by_strack.calculation_for_stagnation_points_formula_by_holzbecher_and_strack()
# discharge_at_stagnation_points = GD_trying_stagnation_points_by_strack.calculation_for_discharge_at_stagnation_points()

# Stagnation_points_lake_plus_wells_no_river = working_for_lake_plus_well_No_River()
# stagnation_points_Lake_and_river = Stagnation_points_lake_plus_wells_no_river.stagnation_points_lake_plus_well()

# lake_plus_river = working_for_lake_plus_well_with_River()
# stagnations_Points_and_filtrate = lake_plus_river.stagnation_points_lake_and_well_plus_river()

# lake_no_river_GD = working_for_lake_plus_well_No_River_GD()
# stagnation_points = lake_no_river_GD.stagnation_points_lake_and_well_noriver_GD()


# plots1 = Plotting()
# colorgrid= plots1.river_color_grid()


plots = Plotting()

while True:
     
    print("1. Do You Want To Solve Well(s), Without River")
    print("2. Do You Want To Solve Well(s), With River")
    print("3. Do You Want To Solve Well(s), With Lake, Withour River")
    print("4. Do You Want To Solve Well(s), With Lake, With River")
    print("5. Do You Want To Solve Well(s), With Line Source / Creek, Without River")
    print("6. Do You Want To Solve Well(s), With Line Source, With River")

    choose_answer = input("Which program do you want to choose from 1-6: ")

    '''aquifer_input = SaveAquiferData()
    aquifer_input.aquifier_save_To_Excel()

    wells_input = SaveWellData()
    wells_input.wells_save_To_Excel()'''  

    if choose_answer == "1":
        complex_plot_unseprated = plots.plotting_complex_without_sepration()
        break
    elif choose_answer == "2":
        Plotting_river_Complex = plots.river_complex()
        break
    elif choose_answer == "3":
        lake_data = ReadAndSaveLakeData()
        #lake_data.save_lake_data()
        plotting_lake_complex = plots.lake_complex()
        break
    elif choose_answer == "4":
        lake_data = ReadAndSaveLakeData()
        #lake_data.save_lake_data()
        plotting_lake_river_complex= plots.lake_river_complex()
        break
    elif choose_answer == "5":
        plotting_complex = plots.line_complex()
        break
    elif choose_answer == "6":
        plotting_line_complex_river= plots.line_complex_river()
        break  
