from abc import ABC, abstractmethod
from typing import List,Tuple
from numpy import arange

class Eos(ABC):
	@abstractmethod
	def __init__(self):
		pass

	@abstractmethod
	def __eq__(self):
		pass

	def __ne__(self, other):
		return not(self == other)

	@abstractmethod
	def __str__(self):
		pass

	@abstractmethod
	def get_p(self, temperature, density, composition, substance, is_silent: bool= True) ->Tuple[float]:
		pass