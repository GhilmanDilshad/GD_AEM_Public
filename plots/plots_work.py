from cProfile import label
from email.base64mime import header_length
from email.quoprimime import header_decode
from re import X
from xml.etree.ElementInclude import XINCLUDE_INCLUDE
import matplotlib
from matplotlib import colors
from matplotlib import artist
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from sympy import N
from vectorcalculations.vector_cords import Vectors
from potentials.potentials_calc import Potentials
from river.river_calc import Potentials_River
from river.river_calc import Vectors_river
from river.river_complex import Vectors_river_complex
from river.river_complex import Potentials_River_complex
from lake.lake_cal import Potentials_Lake
from lake.lake_cal import Vectors_Lake
from lake.lake_complex import Vectors_Lake_Complex
from lake.lake_complex import Potentials_Lake_Complex
from lake.lake_with_river_complex import Vectors_Lake_River_Complex
from lake.lake_with_river_complex import Potentials_Lake_River_Complex

from potentials_complex.potentials_complex_calc import Potentials_normal_well_complex
from potentials_complex.potentials_complex_calc import Vectors_normal_well_complex

from potentials_complex.complex_without_seperation import Vectors_normal_well_complex_without_sepration
from potentials_complex.complex_without_seperation import Potentials_normal_well_complex_without_sepration
from linesink.line_river_complex import Vectors_Line_River_Complex
from linesink.line_river_complex import Potentials_Line_River_Complex
from linesink.line_complex import Potentials_Line_Complex
from linesink.line_complex import Vectors_Line_Complex

from contribution_ratio.contrib import contribution_ratio_river_and_well
import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import collections
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import cmasher as cmr

import seaborn as sns


class Plotting:
    def __init__(self): 
        pass    
        
    #Plotting for head Without River/Constant Head Boundary
    def plot_normal_well_without_river(self):
        abc = 0
        potentials = Potentials()
        head = potentials.calculation_for_head()
        phi_total = potentials.dischargepotential_total_of_region()
        vectors = Vectors()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, [X, Y]= vectors.region_boundaries()

        stream_fxn = potentials.calculation_stream_fxn()
    
        figsize=plt.figure(figsize=(8, 5))
        # fig2, ax = plt.subplots()
        # plt.subplot(111, aspect=1)
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 200)
        plt.ylim(0, 200)
        # levels = H
        
#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        contour1 = plt.contourf (X, Y, head, 38, cmap=cm.winter, alpha=0.75)
        contour2 = plt.contour (X, Y, phi_total, levels= 40, colors='black', linestyles='dashed' )
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        contour3 = plt.contour(X, Y, stream_fxn, levels=38, 
            colors='black',
            linestyles='solid')
        contour2.changed()
        labels = ['Streamline', 'Potentialline']
        contour3.collections[8].set_label(labels[0])
        contour2.collections[9].set_label(labels[1])
        plt.legend(loc='upper right')

#-------Quiver
        # slice_interval = 6
        # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
        # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
            # units = 'height',    
            # angles = 'xy',
            # scale = 50,
            # headwidth= 2    )
            # linewidths=1.0, alpha=0.25, width=1, zorder=3)
       
                # rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
        # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 
#-------Streamplot
        # plt.streamplot(X, Y, v , u,
            # color='white',
           
            # linewidth=1.6, density=1.0, arrowsize=1.2, zorder=3)

        # plt.plot(-d, 0, 'ko')
        plt.show()

        return abc
    
    
    #Plotting for head Without River/Constant Head Boundary

    def plotting_complex(self):
        abcomplex = 0
        potentials_complex = Potentials_normal_well_complex()
        head = potentials_complex.calculation_for_head_complex()
        phi_total_complex = potentials_complex.complex_potential_phi_total_region()
        stream_fxn_complex = potentials_complex.complex_stream_total_region()
        
        vectors_complex = Vectors_normal_well_complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z = vectors_complex.region_boundaries_normal_well_complex()
        

        Z_real = Z.real
        Z_imag = Z.imag
        figsize=plt.figure(figsize=(8, 5))
        # fig2, ax = plt.subplots()
        # plt.subplot(111, aspect=1)
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 200)
        plt.ylim(0, 200)
        # levels = H
        
#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        contour1 = plt.contourf (Z_real, Z_imag, head, 38, cmap=cm.plasma, alpha=0.75)
        contour2 = plt.contour (Z_real, Z_imag, phi_total_complex, levels= 40, colors='black', linestyles='dashed' )
        
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        contour3 = plt.contour(Z_real, Z_imag, stream_fxn_complex, levels=38, 
            colors='red',
            linestyles='solid')
        contour2.changed()
        labels = ['Streamline', 'Potentialline']
        contour3.collections[8].set_label(labels[0])
        contour2.collections[9].set_label(labels[1])
        plt.legend(loc='upper right')

#-------Quiver
        # slice_interval = 6
        # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
        # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
            # units = 'height',    
            # angles = 'xy',
            # scale = 50,
            # headwidth= 2    )
            # linewidths=1.0, alpha=0.25, width=1, zorder=3)
       
                # rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
        # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 
#-------Streamplot
        # plt.streamplot(X, Y, v , u,
            # color='white',
           
            # linewidth=1.6, density=1.0, arrowsize=1.2, zorder=3)

        # plt.plot(-d, 0, 'ko')
        plt.show()

        return abcomplex
   
    def plotting_complex_without_sepration(self):
        abc2omplex = 0
        potentials_unseprated_Complex = Potentials_normal_well_complex_without_sepration()
        head = potentials_unseprated_Complex.calculation_for_head_complex_without_sepration()
        phi_total_complex = potentials_unseprated_Complex.phi_well_plus_baseflowphi_psi_well_plus_baseflowpsi_complex()
        
        
        vectors_complex_unseprated = Vectors_normal_well_complex_without_sepration()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z = vectors_complex_unseprated.region_boundaries_normal_well_complex_without_sepration()

        Z_real = Z.real
        Z_imag = Z.imag
        discharge_potential = phi_total_complex.real
        stream_fxn = phi_total_complex.imag

        figsize=plt.figure(figsize=(16, 9))
        # fig2, ax = plt.subplots()
        # plt.subplot(111, aspect=1)
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 200)
        plt.ylim(0, 200)
        # levels = H
        
#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        contour1 = plt.contourf (Z_real, Z_imag, head, 38, cmap=cm.ocean, alpha=0.70)
        contour2 = plt.contour (Z_real, Z_imag, discharge_potential, levels= 41, colors='white', linestyles='dashed' )
        
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        contour3 = plt.contour(Z_real, Z_imag, stream_fxn, levels=38, 
            colors='white',
            linestyles='solid')
        contour2.changed()
        labels = ['Streamline', 'Potentialline']
        contour3.collections[8].set_label(labels[0])
        contour2.collections[9].set_label(labels[1])
        plt.legend(loc='upper right')

#-------Quiver
        # slice_interval = 6
        # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
        # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
            # units = 'height',    
            # angles = 'xy',
            # scale = 50,
            # headwidth= 2    )
            # linewidths=1.0, alpha=0.25, width=1, zorder=3)
       
                # rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
        # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 
#-------Streamplot
        # plt.streamplot(X, Y, v , u,
            # color='white',         
            # linewidth=1.6, density=1.0, arrowsize=1.2, zorder=3)

        saving_folder = os.path.join("potentials_complex","Wells_Without_River.png")
        plt.savefig(saving_folder,figsize=(1920,1080), dpi = 300)
        plt.show()

        return abc2omplex        
        
#     #Not being in used 
#     def river(self):
#         Priver = Potentials_River()
#         head = Priver.calculation_for_head_river()
#         _, x_coords, h0, H, _, _, _, por = Priver.potentials_data_river()
#         psi =  Priver.calculation_stream_fxn_river()
#         phi_river = Priver.dischargepotential_total_of_region_river()
    
#         vectorsriver = Vectors_river()
#         xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, [X, Y]= vectorsriver.region_boundaries_river()    
#         d=x_coords

#         [u, v] = np.gradient(-head)
#         u = u / head / (xmesh[2]- xmesh[1]) /por
#         v = v / head / (ymesh[2] - ymesh[1]) / por
#         e= 1
        
#         figsize=plt.figure(figsize=(8, 5))
#         # fig2, ax = plt.subplots()
#         # plt.subplot(111, aspect=1)
#         plt.xlabel('X-Distance')
#         plt.ylabel('Y-Distance')  
#         plt.xlim(-200, 200)
#         plt.ylim(0, 200)
#         # levels = H
# #--------River Line Only For Graphical Representation of River
#         # plt.axvline(0, color= '#00FFFF', lw=1.5)

# #--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
#         contour1 = plt.contourf (X, Y, head, 38, cmap=cm.winter, alpha=0.75, zorder=1)
#         contour2 = plt.contour (X, Y, phi_river, levels= 40, colors='white', linestyles='dashed', zorder=2 )
#         plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
#         contour3 = plt.contour(X, Y, psi, levels=38, colors='white', linestyles='solid', zorder=3)
#         contour2.changed()
#         labels = ['Streamline', 'Potentialline']
#         contour3.collections[8].set_label(labels[0])
#         contour2.collections[9].set_label(labels[1])
#         plt.legend(loc='upper right')

# #-------Quiver
#         # slice_interval = 6
#         # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
#         # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
#             # units = 'height',    
#             # angles = 'xy',
#             # scale = 50,
#             # headwidth= 2    )
#             # linewidths=1.0, alpha=0.25, width=1, zorder=3)
       
#                 # rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
#         # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 
# #-------Streamplot
#         # plt.streamplot(X, Y, v , u,
#             # color='white',
           
#             # linewidth=1.6, density=1.0, arrowsize=1.2, zorder=3)
#         # plt.text(0, 0, 'river', rotation=-90, ha='center', va='center')
#         # plt.plot(-d, 0, 'ko')
#         plt.savefig('Plot01.png')
#         plt.show()

#         return d

    def river_complex(self):
        Ghilman = 0
        Rivercomplex = Potentials_River_complex()
        
        head_complex_real = Rivercomplex.calculation_for_head_river_complex()
        phi_river_total_complex = Rivercomplex.dischargepotential_total_of_region_river_complex()
        _, _, _, _, _, _, pumping_array, Zw, alpha, por = Rivercomplex.potentials_data_river_complex()

        # head_real = head.real

        vectorsrivercomplex = Vectors_river_complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z= vectorsrivercomplex.region_boundaries_river_complex()    
        
        Z_real = Z.real
        Z_imag = Z.imag
        discharge_potential = phi_river_total_complex.real
        stream_fxn = phi_river_total_complex.imag

        stream_stagnation_data = contribution_ratio_river_and_well()
        capture_length, first_point, second_point = stream_stagnation_data.calculation_capture_lentgh()
        complex_potential_stagnation_point1 = stream_stagnation_data.calculation_for_contribution_by_stream_fxn()
        figsize1=plt.figure(figsize=(16, 9)) 
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 200)
        plt.ylim(0, 200)
        # levels = H
#--------River Line Only For Graphical Representation of River
        plt.axvline(0, color= '#00FFFF', lw=25, zorder= 6)
        plt.text(1, 100, 'River', rotation=-90, ha='center', va='center', zorder=7)
        plt.plot(0, 0, 'ko')
        #Plotting Capture Lentgh
        # plt.plot((0,0),(first_point,second_point), color='maroon',linestyle="--", linewidth=15, zorder=8)
        # contour_stagnation = plt.contour(Z_real,Z_imag, complex_potential_stagnation_point1.imag,color='maroon',linestyle="--", linewidth=10, zorder=9  )
        # circle = plt.Circle((x_Lake,y_Lake), 25, color='#66b6ff', alpha=1, zorder=5)
        # fig = plt.gcf()
        # ax = fig.gca()
        # ax.add_patch(circle)
#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        cmr_corlormap= cmr.ocean
        contour1 = plt.contourf (Z_real, Z_imag, head_complex_real, 201, cmap=cm.winter, alpha=0.75, zorder = 1) #cm.Blues
        contour2 = plt.contour (Z_real, Z_imag, discharge_potential, levels= 51, colors = 'White', 
        linestyles='dashed', zorder=2, ) # correct phi base with phi  gist_rainbow
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        # contour3 = plt.contour(Z_real, Z_imag, stream_fxn, levels=35, 
        #     cmap=cm.Wistia,
        #     linestyles='solid' ,zorder=3)
        contour2.changed()
#-------Streamplot
        [v, u] = np.gradient(-head_complex_real)

        u = u / head_complex_real / (xmesh[2] - xmesh[1]) / por
        v = v / head_complex_real / (ymesh[2] - ymesh[1]) / por
        e= 1
        
        contour4 = plt.streamplot(Z.real, Z.imag , u , v,
        
            linewidth=1.6, density=(1, 1.2), arrowsize=1.75, color = 'chartreuse' ,zorder=4)   #  cmap=cm.PRGn

        # plt.quiver(X_axis, Y_axis, v ,u , color= 'purple', linewidths=0.1, alpha=0.5, width=0.001)

        # label1= contour2.clabel(contour2, inline=1, fontsize=10)
        # lines=[contour2.collection[0]]
        
        # labels = ['Potentialline','Streamline']
        # contour2.collections[8].set_label(labels[0])
        # contour3.collections[9].set_label(labels[1])
        
        # # contour4.set_label(labels[2])
        # plt.legend(loc='upper right')

        saving_folder = os.path.join("river","Wells_With_River.png")
        plt.savefig(saving_folder, dpi = 300)
        plt.show()


   #Not Being in Used 
#     def lake(self):
#         Plake = Potentials_Lake()
#         _, x_coords, h0, H, _, _, _, por, _, _, _, _  = Plake.potentials_data_lake()
        
#         head = Plake.calculation_for_head_lake()
#         phi_base = Plake.phi_base_fxn_lake()  # just for checking lake ka circluar equipotential ban raha hai ya nahi 
#         phi_0_lake = Plake.phi_0_lake()
#         phi_well_lake = Plake.calculation_phi_well_lake()
#         phi_lake_total = Plake.dischargepotential_total_of_region_lake()
#         psi =  Plake.calculation_stream_fxn_lake()
        

#         vectorslake = Vectors_Lake()
#         xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, [X, Y]= vectorslake.region_boundaries_lake()    
    
#         d=x_coords
#         [u, v] = np.gradient(-head)
#         u = u / head / (xmesh[2]- xmesh[1]) /por
#         v = v / head / (ymesh[2] - ymesh[1]) / por
#         e= 1
        
        
        
#         figsize=plt.figure(figsize=(5, 5))
#         plt.xlabel('X-Distance')
#         plt.ylabel('Y-Distance')  
#         plt.xlim(-200, 200)
#         plt.ylim(0, 200)
#         # levels = H
# #--------River Line Only For Graphical Representation of River
#         plt.axvline(0, color= '#00FFFF', lw=0.5)
#         circle = plt.Circle((100,100), 25, color='#66b6ff', alpha=1, zorder=4)
#         fig = plt.gcf()
#         ax = fig.gca()
#         ax.add_patch(circle)
# #--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
#         contour1 = plt.contourf (X, Y, head, 40, cmap=cm.YlGnBu, alpha=1)
#         contour2 = plt.contour (X, Y, (phi_lake_total), levels= 200, colors='white', linestyles='dashed' ) # correct phi base with phi
#         plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
#         contour3 = plt.contour(X, Y, psi, levels=250, 
#             colors='black',
#             linestyles='solid')
#         contour2.changed()
#         labels = ['Streamline', 'Potentialline']
#         contour3.collections[8].set_label(labels[0])
#         contour2.collections[9].set_label(labels[1])
#         plt.legend(loc='upper right')

# #-------Quiver
#         # slice_interval = 6
#         # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
#         # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
#             # units = 'height',    
#             # angles = 'xy',
#             # scale = 50,
#             # headwidth= 2    )
#             # linewidths=1.0, alpha=0.25, width=1, zorder=3)
       
#                 # rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
#         # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 
# #-------Streamplot
#         plt.streamplot(X, Y, v , u,
#             color='black',
#             linewidth=1.6, density=1.0, arrowsize=1.2, zorder=3)
#         # plt.text(0, 0, 'river', rotation=-90, ha='center', va='center')
#         # plt.plot(-d, 0, 'ko')
#         plt.savefig('Normal_Lake.png')
#         plt.show()

#         return d
    
    
    def lake_complex(self):
        Ghilman = 0
        Plakecomplex = Potentials_Lake_Complex()
        
        head = Plakecomplex.calculation_for_head_lake_complex()
        phi_lake_total_complex = Plakecomplex.dischargepotential_total_of_region_lake_complex()
        R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake = Plakecomplex.potentials_data_lake_complex()
        
        head_real = head.real
        vectorslakecomplex = Vectors_Lake_Complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis  = vectorslakecomplex.region_boundaries_lake_complex()    
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vectorslakecomplex.vectors_data_lake_complex()

        Z_real = Z.real
        Z_imag = Z.imag
        discharge_potential = phi_lake_total_complex.real
        stream_fxn = phi_lake_total_complex.imag

        
        figsize=plt.figure(figsize=(16, 9)) 
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 400)
        plt.ylim(0, 200)
        # levels = H
#--------River Line Only For Graphical Representation of River
        # plt.axvline(0, color= '#00FFFF', lw=0.5)
        # plt.text(0, 0, 'river', rotation=-90, ha='center', va='center')
        # plt.plot(-d, 0, 'ko')
        circle = plt.Circle((x_Lake,y_Lake), 25, color='#66b6ff', alpha=1, zorder=5)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_patch(circle)
#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        cmr_corlormap= cmr.ocean
        contour1 = plt.contourf (Z_real, Z_imag, head_real, 80, cmap=cm.YlGnBu, alpha=1, zorder = 1) #cm.Blues
        contour2 = plt.contour (Z_real, Z_imag, discharge_potential, levels= 65, colors = '#13F4EF', linestyles='dashed', zorder=2 ) # correct phi base with phi  gist_rainbow
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        contour3 = plt.contour(Z_real, Z_imag, stream_fxn, levels=200, 
            cmap=cm.Wistia,
            linestyles='solid' ,zorder=3)
        contour2.changed()

# #-------Streamplot
        # [v, u] = np.gradient(-head_real)

        # u = u / head_real / (xmesh[2] - xmesh[1]) / por
        # v = v / head_real / (ymesh[2] - ymesh[1]) / por
        # e= 1
        
        # contour4 = plt.streamplot(Z.real, Z.imag, u , v, linewidth=1.6, density=(2.5, 1.2), arrowsize=1.75, color = 'chartreuse' ,zorder=4)
        

        # plt.quiver(X_axis, Y_axis, v ,u , color= 'purple', linewidths=0.1, alpha=0.5, width=0.001)
        # labels = ['Streamline', 'Potentialline', 'Flow Lines']
        # contour3.collections[8].set_label(labels[0])
        # contour2.collections[9].set_label(labels[1])
        # # contour4.set_label(labels[2])
        # plt.legend(loc='upper right')

        saving_folder = os.path.join("lake","Lake_Complex.png")
        plt.savefig(saving_folder, dpi = 300)
        plt.show()


#-------Quiver

        # slice_interval = 6
        # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
        # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
        #     units = 'height',    
        #     angles = 'xy',
        #     scale = 50,
        #     headwidth= 2    
        #     linewidths=1.0, alpha=0.25, width=1, zorder=3)
        #     rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
        # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 

        return Ghilman

    def lake_river_complex(self):
        Ghilman = 0
        Plakecomplex = Potentials_Lake_River_Complex()
        
        head = Plakecomplex.calculation_for_head_lake_river_complex()
        phi_lake_total_complex = Plakecomplex.dischargepotential_total_of_region_lake_river_complex()
        R_Lake, x_Lake, y_Lake, Q_lake, Z_Lake = Plakecomplex.potentials_data_lake_river_complex()
        head_real = head.real

        vectorslakecomplex = Vectors_Lake_River_Complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis  = vectorslakecomplex.region_boundaries_lake_river_complex()    
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vectorslakecomplex.vectors_data_lake_river_complex()

        Z_real = Z.real
        Z_imag = Z.imag

        discharge_potential = phi_lake_total_complex.real
        stream_fxn = phi_lake_total_complex.imag

        
        figsize=plt.figure(figsize=(16, 9)) 
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 400)
        plt.ylim(0, 200)
        # levels = H
#--------River Line Only For Graphical Representation of River
        plt.axvline(0, color= '#00FFFF', lw=25, zorder= 6)
        plt.text(1, 100, 'River', rotation=-90, ha='center', va='center', zorder=7)
        plt.plot(0, 0, 'ko')
        circle = plt.Circle((x_Lake,y_Lake), 25, color='#66b6ff', alpha=1, zorder=5)
        fig = plt.gcf()
        ax = fig.gca()
        ax.add_patch(circle)
#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        cmr_corlormap= cmr.ocean
        contour1 = plt.contourf (Z_real, Z_imag, head_real, 100, cmap=cm.YlGnBu, alpha=1, zorder = 1) #cm.Blues, cmap=cmr_corlormap
        contour2 = plt.contour (Z_real, Z_imag, discharge_potential, levels= 75, colors = '#13F4EF', linewidth=10, linestyles='dashed', zorder=2) # correct phi base with phi  gist_rainbow
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        contour3 = plt.contour(Z_real, Z_imag, stream_fxn, levels=80, cmap=cm.Wistia, linestyles='solid' ,zorder=3)
        contour2.changed()
        # labels = ['Streamline', 'Potentialline']

        # contour2.collections[15].set_label(labels[1])
        # contour3.collections[15].set_label(labels[0])
        # plt.legend(loc='upper right')
#-------Streamplot
        # [v, u] = np.gradient(-head_real)

        # u = u / head_real / (xmesh[2] - xmesh[1]) / por
        # v = v / head_real / (ymesh[2] - ymesh[1]) / por
        # e= 1
        
        # contour4 = plt.streamplot(Z.real, Z.imag, u , v, linewidth=1.6, density=(7, 5), arrowsize=1.75, color = 'chartreuse' ,zorder=4)   #  cmap=cm.PRGn
        
        # plt.quiver(X_axis, Y_axis, v ,u , color= 'purple', linewidths=0.1, alpha=0.5, width=0.001)

        # # contour4.set_label(labels[2])

        
        saving_folder = os.path.join("lake","Lake_River_Complex.png")
        plt.savefig(saving_folder, dpi = 300)
        plt.show()


#-------Quiver

        # slice_interval = 6
        # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
        # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
        #     units = 'height',    
        #     angles = 'xy',
        #     scale = 50,
        #     headwidth= 2    
        #     linewidths=1.0, alpha=0.25, width=1, zorder=3)
        #     rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
        # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 

    def line_complex(self):
        Ghilman = 0
        Plinecomplex = Potentials_Line_Complex()
        
        head = Plinecomplex.calculation_for_head_line_complex()
        phi_line_total_complex = Plinecomplex.dischargepotential_total_of_region_line_complex()
        z1_Line, z2_Line, Sigma_line =  Plinecomplex.potentials_data_line_complex()
        head_real = head.real

        vectorslakecomplex = Vectors_Line_Complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis  = vectorslakecomplex.region_boundaries_line_complex()    
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vectorslakecomplex.vectors_data_line_complex()

        Z_real = Z.real
        Z_imag = Z.imag

        discharge_potential = phi_line_total_complex.real
        stream_fxn = phi_line_total_complex.imag

        
        figsize=plt.figure(figsize=(16, 9)) 
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 400)
        plt.ylim(0, 200)
        # levels = H
#--------Canal Line Only For Graphical Representation of Line
        x_line=z1_Line[0].real,z2_Line[0].real
        y_line=z1_Line[0].imag,z2_Line[0].imag
        plt.plot(x_line,y_line, color='#66b6ff', linewidth=10, zorder=7)
#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        cmr_corlormap= cmr.ocean
        contour1 = plt.contourf (Z_real, Z_imag, head_real, 40, cmap=cmr_corlormap, alpha=0.8, zorder = 1) 
        contour2 = plt.contour (Z_real, Z_imag, discharge_potential, levels= 40, colors = '#13F4EF', linewidth=10, linestyles='dashed', zorder=2) 
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        contour3 = plt.contour(Z_real, Z_imag, stream_fxn, levels=40, cmap=cm.Wistia, linestyles='solid' ,zorder=3)
        # contour2.changed()

#-------Streamplot
        # [v, u] = np.gradient(-head_real)

        # u = u / head_real / (xmesh[2] - xmesh[1]) / por
        # v = v / head_real / (ymesh[2] - ymesh[1]) / por
        # e= 1
        
        # contour4 = plt.streamplot(Z.real, Z.imag, u , v, linewidth=1.6, density=(5, 1.2), arrowsize=1.75, color = 'chartreuse' ,zorder=4)   #  cmap=cm.PRGn

        # labels = ['Streamline', 'Potentialline', 'Flow Lines']
        # contour3.collections[8].set_label(labels[0])
        # contour2.collections[9].set_label(labels[1])
        # # contour4.set_label(labels[2])
        # plt.legend(loc='upper right')

        saving_folder = os.path.join("linesink","Line_Complex.png")
        plt.savefig(saving_folder, dpi = 300)
        plt.show()


#-------Quiver

        # slice_interval = 6
        # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
        # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
        #     units = 'height',    
        #     angles = 'xy',
        #     scale = 50,
        #     headwidth= 2    
        #     linewidths=1.0, alpha=0.25, width=1, zorder=3)
        #     rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)


        return Ghilman






    def line_complex_river(self):
        Ghilman = 0
        Plinecomplex = Potentials_Line_River_Complex()
        
        head = Plinecomplex.calculation_for_head_line_river_complex()
        phi_line_total_complex = Plinecomplex.dischargepotential_total_of_region_line_river_complex()
        z1_Line, z2_Line, Sigma_line =  Plinecomplex.potentials_data_line_river_complex()
        head_real = head.real

        vectorslakecomplex = Vectors_Line_River_Complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z, X_axis, Y_axis  = vectorslakecomplex.region_boundaries_line_river_complex()    
        baseFlowX, x_array, h0, H, k, y_array, pumping_array, Zw, alpha, Zref, por = vectorslakecomplex.vectors_data_line_river_complex()

        Z_real = Z.real
        Z_imag = Z.imag

        discharge_potential = phi_line_total_complex.real
        stream_fxn = phi_line_total_complex.imag

        
        figsize=plt.figure(figsize=(16, 9)) 
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, 400)
        plt.ylim(0, 200)
        # levels = H
#--------River Line Only For Graphical Representation of River
        plt.axvline(0, color= '#00FFFF', lw=25, zorder= 6)
        plt.text(1, 100, 'River', rotation=-90, ha='center', va='center', zorder=7)
        plt.plot(0, 0, 'ko')

#--------Canal Line Only For Graphical Representation of Line
        x_line=z1_Line.real,z2_Line.real
        y_line=z1_Line.imag,z2_Line.imag
        print('These are x_values of creek',x_line)
        print('These are y_values of Creek', y_line)
        plt.plot(x_line,y_line, color='#66b6ff', linewidth=10, zorder=7)

#--------- Head Contours Color, Line Contours for phi, Dashed iso lines for psi    
        cmr_corlormap= cmr.nuclear
        contour1 = plt.contourf (Z_real, Z_imag, head_real, 100, cmap=cm.GnBu, alpha=1, zorder = 1) 
        contour2 = plt.contour (Z_real, Z_imag, discharge_potential, levels= 75, colors = 'lime', linewidth=10, linestyles='dashed', zorder=2) 
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        
        contour3 = plt.contour(Z_real, Z_imag, stream_fxn, levels=80, cmap=cm.autumn, linestyles='solid' ,zorder=3)
        # contour2.changed()

# #-------Streamplot
        # [v, u] = np.gradient(-head_real)

        # u = u / head_real / (xmesh[2] - xmesh[1]) / por
        # v = v / head_real / (ymesh[2] - ymesh[1]) / por
        # e= 1
        
        # contour4 = plt.streamplot(Z.real, Z.imag, u , v, linewidth=1.6, density=(3.5, 0.8), arrowsize=1.75, color = 'chartreuse' ,zorder=4)   #  cmap=cm.PRGn

        # plt.quiver(X_axis, Y_axis, v ,u , color= 'purple', linewidths=0.1, alpha=0.5, width=0.001)
        # labels = ['Streamline', 'Potentialline', 'Flow Lines']
        # contour3.collections[8].set_label(labels[0])
        # contour2.collections[9].set_label(labels[1])
        # # # contour4.set_label(labels[2])
        # plt.legend(loc='upper right')

        saving_folder = os.path.join("linesink","Line_River_Complex.png")
        plt.savefig(saving_folder, dpi = 300)
        plt.show()


#-------Quiver

        # slice_interval = 6
        # skip = slice(None, None, slice_interval), slice(None, None, slice_interval)
        # plt.quiver(X[skip], Y[skip], v[skip], u[skip],
        #     units = 'height',    
        #     angles = 'xy',
        #     scale = 50,
        #     headwidth= 2    
        #     linewidths=1.0, alpha=0.25, width=1, zorder=3)
        #     rect = Rectangle((500, 200),2 ,400, linewidth=1, edgecolor='b', facecolor='b', zorder=1)
        # plt.contour(X, Y, psi, [-160, 160], color='k')is ki zrort nahi yeh Q/2 sa capture zone batata hai 

        return Ghilman
      

    def river_color_grid(self):

        Rivercomplex = Potentials_River_complex()
        
        head_complex_real = Rivercomplex.calculation_for_head_river_complex()
        phi_river_total_complex = Rivercomplex.dischargepotential_total_of_region_river_complex()
        _, _, _, _, _, _, pumping_array, Zw, alpha, por = Rivercomplex.potentials_data_river_complex()

        # head_real = head.real

        vectorsrivercomplex = Vectors_river_complex()
        xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, Z= vectorsrivercomplex.region_boundaries_river_complex()    
        
        Z_real = Z.real
        Z_imag = Z.imag
        discharge_potential = phi_river_total_complex.real
        stream_fxn = phi_river_total_complex.imag
        [x_axis, y_axis]= np.meshgrid(xmesh,ymesh)
        stream_stagnation_data = contribution_ratio_river_and_well()
        capture_length, first_point, second_point = stream_stagnation_data.calculation_capture_lentgh()
        complex_potential_stagnation_point1 = stream_stagnation_data.calculation_for_contribution_by_stream_fxn()

        fig = plt.figure(figsize=(8, 8),facecolor='w')
        ax = fig.add_subplot(111)

        # Major ticks every 20, minor ticks every 5
        major_ticks = np.arange(0, (Xmax+25), 25)
        minor_ticks = np.arange(0, (Ymax+25), 5)
        print(major_ticks)
        print(minor_ticks)

        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)

        # And a corresponding grid
        # ax.grid(which='both',linewidth=1.1, color='k')

        
        ax.grid(b=True, which='major', color='k', linestyle='-', zorder=5)
        ax.grid(which='minor', alpha=1, color='w')
        
        plt.xlabel('X-Distance')
        plt.ylabel('Y-Distance')  
        plt.xlim(0, Xmax)
        plt.ylim(0, Ymax)

        contour1 = plt.pcolor(x_axis,y_axis,head_complex_real, cmap='rainbow')
        plt.colorbar(contour1, shrink = 0.7, label= ' Head (m) ')
        plt.grid(color='w', linewidth=1.2)


        plt.show()
        # ax =sns.heatmap(head_complex_real, cmap='hot', square=True,linecolor='black')
        # plt.show()

        return None