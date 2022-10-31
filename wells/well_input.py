class WellInput:

	def __init__(self):
		pass

	def number_of_wells(self):
		# enter wells upto 5
		w_id = []
		x_coords = []
		y_coords = []
		pumping = []
		while True:
			well_id = int(input("Please Enter number of wells between 1 and 5: "))
			if well_id >= 1 and well_id <= 5:
				for i in range(well_id):
					print("Please Enter parameters for well id: {}".format(i+1))
					x_coordinates, y_coordinates, pumping_rate = self.input_parameters()
					w_id.append(i + 1)
					x_coords.append(x_coordinates)
					y_coords.append(y_coordinates)
					pumping.append(pumping_rate)
				break
			print("Please Enter number of wells between 1 and 5")
		return w_id, x_coords, y_coords, pumping
	
	def input_parameters(self):
		# radius of well
		r = 0.2
				
		while True:
			# get x-coordinate from user
			x_coordinates = float(input("Please Enter x-coordinates: "))
			if x_coordinates >= r:
				x_coordinates = x_coordinates
				break
			print("X- coordiantes must be greater than radius of well (0.2)")

		while True:
			# get y-coordinate from user
			y_coordinates = float(input("Please Enter y-coordinates: "))
			if y_coordinates >= r:
				y_coordinates = y_coordinates
				break
			print("X- coordiantes must be greater than radius of well (0.2)")
		

		# get Base Flow in X direction from user
		pumping_rate = float(input("Please Enter pumping rate: "))

		return x_coordinates, y_coordinates, pumping_rate

		