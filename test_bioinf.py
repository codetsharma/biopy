import numpy as np


# region 1 climate
# []

kanto = np.array([73,86,64])

climate2 = [0.5,0.2,0.3]
sabby = np.array(climate2)

print ("\nkanto: ", type(kanto))
print ("\nsabby: ", type(sabby))

apple_yield = np.dot(kanto, sabby)
print ("\napple yield: \n", apple_yield)