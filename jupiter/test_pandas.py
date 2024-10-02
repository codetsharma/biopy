# test_pandas
import numpy as np
import pandas as pd

#
# Lists of data
# Columns of data for Revenue, Employees, Sector:
#
data = {'Revenue': [274515,200734,182527,181945,143015,129184,92224,85965,84893,
                    82345,77867,73620,69864,63191],
        'Employees': [147000,267937,135301,878429,163000,197000,158000,58604,
                      109700,350864,110600,364800,85858,243540],
        'Sector': ['Consumer Electronics','Consumer Electronics','Software Services',
                   'Chip Manufacturing','Software Services','Consumer Electronics',
                   'Consumer Electronics','Software Services','Consumer Electronics',
                   'Consumer Electronics','Chip Manufacturing','Software Services',
                   'Software Services','Consumer Electronics'],
        'Founding Date':['01-04-1976','13-01-1969','04-09-1998','20-02-1974',
                         '04-04-1975','15-09-1987','01-02-1984','04-02-2004',
                         '07-04-1946','01-01-1910','18-07-1968','16-06-1911',
                         '11-11-1998','07-03-1918'],
        'Country':['USA','South Korea','USA','Taiwan','USA','China','USA','USA',
                   'Japan','Japan','USA','USA','China','Japan']} 
index = ['Apple','Samsung','Alphabet','Foxconn','Microsoft','Huawei',
         'Dell Technologies','Meta','Sony','Hitachi','Intel','IBM',
         'Tencent','Panasonic']

df = pd.DataFrame(data, index)
print ("\ndf: \n", df)
# First 5 rows by default
print ("\ndf.head(): \n", df.head())
print ("\ndf.head(13): \n", df.head(13))
print ("\ndf.tail(): \n", df.tail())
print ("\ndf.shape: \n", df.shape)
print ("\ndf.describe():\n", df.describe())
print ("\ndf.nunique(): \n", df.nunique()) # No. of Unique Values. 
print ("\ndf.info(): \n", df.info())

rt = np.array([[1,2],[3,4]], dtype=np.int8)
print (rt)
print ("\nrt.dtype: \n", rt.dtype)
print (df[['Revenue', 'Employees','Country']])

df_subset = df[['Revenue', 'Employees']]
print ("\ndf: \n", df)
print ("\ndf_subset: \n", df_subset)
print ("\ndf subset mean: \n", df_subset.mean())

print ("\nOutput a Series with the column Employees: \n")
employees_s = df['Employees']
print (employees_s)
print ("\nNOTE: A Series from a dataframe has the name: <col name> and index = all rows.\n")

employees_median = employees_s.median()
employees_median