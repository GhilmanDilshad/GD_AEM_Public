from wells.well_input import WellInput
import pandas as pd

class SaveWellData:

	def __init__(self):
		pass

	def wells_save_To_Excel(self):
		inputData = WellInput()
		well_id, x_coords, y_coords, pumping = inputData.number_of_wells()
		data = pd.DataFrame({"Well ID":well_id, "x-coordinates":x_coords, "y-coordinates": y_coords, "pumping": pumping})
		data.to_excel('wells/well_sheet.xlsx', sheet_name='sheet1', index=False)
		return None