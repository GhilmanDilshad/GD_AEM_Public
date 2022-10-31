class AquiferInput:

	def __init__(self):
		pass

	def input_id(self):
		# set aquifer id is 1
		aq_id = 1
		# get parameters for aquifer id is 1
		k_value, H_Thickness, baseFlowX, h0, por = self.input_parameters()
		""" while True:
			aid = int(input("Please Enter 1 and 2 for an Acquifer ID: "))
			if aid == 1:
				k_value, H_Thickness, baseFlowX, phi_RH, por = self.input_parameters()
				break
			elif aid == 2:
				k_value, H_Thickness, baseFlowX, phi_RH, por = self.input_parameters()
				break
			print("Please Enter Acquifer ID 1 or 2") """
		return aq_id, k_value, H_Thickness, baseFlowX, h0, por
	
	def input_parameters(self):
		# get hydraulic conductivity from user
		k_value = float(input("Please Enter hydraulic conductivity [m/d]: "))

		# get Aquifer Thickness from user
		H_Thickness = float(input("Please Enter Aquifer Thickness [m]: "))

		# get Base Flow in X direction from user
		baseFlowX = float(input("Please Enter Base Flow in X direction [m\u00B2/d]: "))

		# get Reference Head from user
		h0 = float(input("Please Enter Reference Head [m]: "))

		# get Porosity from User
		por = float(input("Please Enter Porosity: "))

		return k_value, H_Thickness, baseFlowX, h0, por



		