from TREND_4.conversion_table import convert2trend_name, convert2pcsaft_name

class Property:	
	name2unitname = {
					"MW"   : "kg/mol",
					"Ttrip": "K",
					"ptrip": "MPa",
					"Tcrit": "K",
					"pcrit": "MPa",
					"Dcrit": "mol/m^3",
					"AF"   : "#",
					"Tmin" : "K",
					"Tmax" : "K",
					"pmax" : "MPa",
					"Dmax" : "mol/m^3",}

	def __init__ (self, name, value):
		self.name = str(name).strip()
		# print(name, len(value))
		if (len(value) == 0):
			self.value = None
		else:
			self.value = float(value)
		if (self.name in self.name2unitname.keys()):
			self.unitname = self.name2unitname[self.name]
		else:
			self.unitname = "?"

	def __str__(self):
		if (self.value is None):
			return f"{self.name} = {'None'} {self.unitname}"
		else:
			return f"{self.name} = {self.value} {self.unitname}"
# data structures 
class Substances:
	def __init__(self, names, ratios):
		tmp_sum = 0
		self.count = 0
		self.name = []
		self.ratio = []

		for n,r in zip(names,ratios):
			self.name.append(n)
			self.ratio.append(r)
			tmp_sum += r
			if tmp_sum > 1 :
				raise(Error("Molar ratio sum exeeds 1.0"))
			self.count += 1 

	def __eq__(self, other):
		if (self.count != other.count):
			return False
		if (self.name != other.name):
			return False
		if (self.count != other.count):
			return False
		return True

	def __ne__(self, other):
		return not(self == other)

	def get_pcsaft_str(self):
		ret = ""
		for ind,el in enumerate(self.name):
			el = convert2pcsaft_name(el)
			ret += el
			if (ind < self.count - 1):
				ret += ","
		return ret

	def get_trend_name_str(self):
		ret = ""
		# print(self.name)
		for el in self.name:
			el = convert2trend_name(el)
			ret += str.ljust(el,30)
			ret += "\n"
		for i in range(30-self.count):
			ret += " "*30 + "\n"
		return ret

	def get_trend_ratio_str(self):
		ret = ""
		for el in self.ratio:
			ret += str.ljust('{:18.16f}'.format(el),26)
		for i in range(30-self.count):
			ret += str.ljust('{:18.16f}'.format(0.0),26)

		ret += "\n"
		return ret
