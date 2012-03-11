import pyhalo

C = pyhalo.Cosmology()
z = 1
print pyhalo.growth(1,C) / pyhalo.growth(0,C)
