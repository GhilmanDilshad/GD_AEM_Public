import pandas as pd
import numpy as np


def read_aquifer_xlsx():
    dataframe = pd.read_excel("aquifer/aquifer.xlsx", index_col=None, na_values=['NA'], usecols = 'A,B:AA') 
    return dataframe

def read_wells_xlsx():
    dataframe = pd.read_excel("wells/well_sheet.xlsx", index_col=None, na_values=['NA'], usecols = 'A,B:AA') 
    return dataframe

def read_lake_xlsx():
    dataframe = pd.read_excel("lake/lake_data.xlsx", index_col=None, na_values=['NA'], usecols = 'A,B:AA') 
    return dataframe

