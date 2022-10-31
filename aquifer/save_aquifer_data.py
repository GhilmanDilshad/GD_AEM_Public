from aquifer.aquifer_input import AquiferInput
import pandas as pd

class SaveAquiferData:

	def __init__(self):
		pass

	def aquifier_save_To_Excel(self):
		inputData = AquiferInput()
		aq_id, k_value, H_Thickness, baseFlowX, h0, por = inputData.input_id()
		data = pd.DataFrame({"Aquifer ID":aq_id, "Hydraulic Conductivity":k_value, "Aquifer Thickness": H_Thickness, "Base Flow in X direction": baseFlowX, "Reference Head": h0, "Porosity": por}, index=[0])
		data.to_excel('aquifer/aquifer.xlsx', sheet_name='sheet1', index=False)
		return None


