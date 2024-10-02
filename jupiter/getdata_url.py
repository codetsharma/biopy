import urllib.request
import numpy as np

urllib.request.urlretrieve(
    'https://hub.jovian.ml/wp-content/uploads/2020/08/climate.csv',
    'climate.txt'
)

climate_data = np.genfromtxt('climate.txt', delimiter=',', skip_header=1)
print (climate_data)
print (climate_data.shape)