"""
GOAL: Exploring how to use python's pandas library:
pandas library is widely used for data manipulation and analysis
"""

import pandas as pd

# ######################################################################
# Creating a Dataframe:
# ######################################################################
# Note: Pandas uses "DataFrames" as the Primary datastructure 
# ~~~~~ df:          DataFrame
# Python variables like list, dictionary can be converted into Dataframe
# 
# 1. DataFrame - From a Dictionary:
data_Dict = {
    'Name': ['Alice',       'Bob',      'Charlie', 'David'],  
    'Age' : [25,            30,             35,       40],
    'City': ['New York', 'Los Angeles', 'Chicago', 'Houston'],              
}

df_Dict = pd.DataFrame(data_Dict)
print ("\nDataframe from DICT:\n")
print (df_Dict)


# 2. DataFrame - From a List of Lists:
data_List = [
    ['Alice',   25, 'New York'],
    ['Bob'  ,   30, 'Los Angeles'],
    ['Charlie', 35, 'Chicago'],
    ['David',   40, 'Houston']
]

df_List = pd.DataFrame(data_List)
print ("\nDataframe from LIST:\n")
print (df_List)

# 3. DataFrame - From a CSV File:
# with pd.read_csv (csv_file, header=0, index_col=0, usecols=['name','type','start_time','duration','mosi','miso'], dtype={'name': str, 'type':str, 'start_time': float, 'duration': float, 'mosi': int, 'miso': int}) as df:
df = pd.read_csv ('csv2.csv')
print ("\nDataframe from file: \n")
print (df.head()) #display 1st 5 rows.


