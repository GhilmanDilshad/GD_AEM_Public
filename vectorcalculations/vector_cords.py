import pandas as pd
import numpy as np
from utilities import read_aquifer_xlsx, read_wells_xlsx


class Vectors:
    def __init__(self):  
        pass    

    def region_boundaries(self): # Region Coordinates, boundaries and meshgrids.
        Xmin = 0
        Xmax = 200
        Ymin = 0
        Ymax = 200
        xmesh = np.linspace(Xmin, Xmax, 201)
        ymesh = np.linspace(Ymin, Ymax, 201)
        [X, Y] = np.meshgrid(xmesh, ymesh)

        return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, [X, Y]

    def vectors_data(self):    # getting coordinates of well(s) and converting them to list 
        df3 = read_wells_xlsx()   # This is for reading values from well file
       
        x_values = df3['x-coordinates'].tolist()   
        y_values = df3['y-coordinates'].tolist()
       
        pumping_rate = df3['pumping'].tolist()
        pumping_array = np.array(pumping_rate)

        x_array = np.array(x_values) # converting list to array
        y_array = np.array(y_values) # converting list to array
        return x_values, y_values, x_array, y_array, pumping_array


    # def region_boundaries_river(self): # Region Coordinates, boundaries and meshgrids.
    #     Xmin = -200
    #     Xmax = 200
    #     Ymin = -200
    #     Ymax = 200
    #     xmesh = np.linspace(Xmin, Xmax, 100)
    #     ymesh = np.linspace(Ymin, Ymax, 100)
    #     [X, Y] = np.meshgrid(xmesh, ymesh)

    #     return xmesh, ymesh, Xmin, Xmax, Ymin, Ymax, [X, Y]

    # def vectors_data_river(self):    # getting coordinates of well(s) and converting them to list 
    #     df3 = read_wells_xlsx()   # This is for reading values from well file
       
    #     x_values = df3['x-coordinates'].tolist()   
    #     y_values = df3['y-coordinates'].tolist()
       
    #     pumping_rate = df3['pumping'].tolist()
    #     pumping_array = np.array(pumping_rate)

    #     x_array = np.array(x_values) # converting list to array
    #     y_array = np.array(y_values) # converting list to array
    #     return x_values, y_values, x_array, y_array, pumping_array