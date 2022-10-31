import pandas as pd
	
class ReadAndSaveLakeData():
	
	def __init__(self) -> None:
		pass

	def save_lake_data(self):
		lake_id, radius,x_coordinates, y_coordinates, inflow_outflow = self.number_of_lakes()
		data = pd.DataFrame({"Lake ID":lake_id, "radius":radius, "x-coordinates":x_coordinates, "y-coordinates": y_coordinates, "inflow_outflow": inflow_outflow})
		data.to_excel('lake/lake_data.xlsx', sheet_name='sheet1', index=False)
		return None

	def number_of_lakes(self):
		# enter wells upto 5
		l_id = []
		r = []
		x_coords = []
		y_coords = []
		iNO = []
		while True:
			lake_id = int(input("Please Enter number of lake between 1 and 5: "))
			if lake_id == 1:
				for i in range(lake_id):
					print("Please Enter parameters for well id: {}".format(i+1))
					radius,x_coordinates, y_coordinates, inflow_outflow = self.input_data()
					l_id.append(i + 1)
					x_coords.append(x_coordinates)
					y_coords.append(y_coordinates)
					r.append(radius)
					iNO.append(inflow_outflow)
				break
			print("Please Enter number of wells between 1 and 5")
		return l_id, r, x_coords, y_coords, iNO
	

	def input_data(self):

		# get Base Flow in X direction from user
		radius = float(input("Please Enter Radius of lake: "))
		x_coordinates = float(input("Please Enter X-Coordinate of lake: "))
		y_coordinates = float(input("Please Enter Y-Coordinate of lake: : "))
		inflow_outflow = float(input("Please Enter Inflow/Outflow in 'm\u00b3/d': "))

		return radius,x_coordinates, y_coordinates, inflow_outflow

		