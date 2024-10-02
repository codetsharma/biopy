import sys
import numpy as np


a = np.array([1,2,3,4,4,5,6,7,8])
b = np.array([0, .5, 1, 1.5, 2])
print (a[0], a[1])
print ("\na[1] to last, exclude last!! ", a[1:-1])
print ("\na[0:]: ", a[0:], "a[1:3]: ", a[1:3], "a[1:-1]: ", a[1:-1], "a[::2]: ", a[::2])
print ("\nb[0]: ", b[0], "b[2]: ", b[2], "Last element b[-1]: ", b[-1])
print ("\nList of indices passed in Multi-Index Array: b[[0,2,-1]]: ", b[[0,2,-1]])

print (a.dtype)
c = np.array([1,2,3,4], dtype=np.int8)

print (c.dtype)
if (a.dtype != c.dtype):
    print ("Ooh La La!!!")

d = np.array(['a', 'b', 'c'])
A = np.array([[1,2],
             [3,4]])
B = np. array ([
    [1,3,6],
    [4,5,9],
    [2,3,4]
    ])

print ("\nShape of B: ", B.shape)
print ("\nSize of B: ",  B.size)
print ("\nB[1,1]: ",     B[1,1])
print ("All of B :) \n", B)
print ("\nB[0:3,2]: \n", B[0:3,2].transpose())
print ("\nB[0:3,:2]: \n", B[0:3,:2])

B[2] = 99
print ("\nB with row 3 as all 99: \n", B)
print ("\nB[:2, 2:]: \n", B[:2, 2:])

# Summary Statistics:
print ("\nSummary Statistics: \n")
h = np.array ([6,7,8,9,10])
print ("\nh.sum(): \n", h.sum())
print ("\nh.mean(): \n", h.mean())
print ("\nh.var()iance: \n", h.var())

i = np.array([
    [1,2,3],
    [4,5,6],
    [7,8,9]
], dtype=np.int8)
print (i)
print("\ni.sum(): \n", i.sum())
print("\ni.mean(): \n", i.mean())
print("\ni.std(): \n",                      i.std())
print("\ni.sum(axis=0):  column-wise: \n", i.sum(axis=0))
print("\ni.sum(axis=1):  row-wise \n", i.sum(axis=1))
print("\ni.mean(axis=0): column-wise: \n", i.mean(axis=0))
print("\ni.mean(axis=1): row-wise: \n", i.mean(axis=1))

print ("\nBroadcasting , Vectorized Operations\n")
print ("\ni*5\n", i*5)
print ("\ni+5\n", i+5)
print ("\ni-5\n", i-5)
print ("\ni/5\n", i/5)

print ("\nBoolean Arrays\n")
print ("\ni >= 5\n", i>=5)
print ("\ni[i>=5]\n", i[i>=5])
print ("\ni[0:,-1]:\n", i[0:,-1])
print ("i [even numbers]:", i[i % 2 == 0]) 
Y = np.random.randint(150, size=(3,4))
print ("\nnp.random.randint(150, size=(3,4))\n", (Y))

print ("\nIt takes ", sys.getsizeof(1), " bytes to hold 1 in python!!!!!")
print ("\nsys.getsizeof(10*1000): \n",sys.getsizeof(10*1000))
print ("\nnp.dtype(int).itemsize: \n", np.dtype(np.int8).itemsize)
print ("\nnp.dtype(float).itemsize: \n", np.dtype(float).itemsize)
print ("\nsys.getsizeof([1]): \n", sys.getsizeof([1]))
print ("\nnp.array([1]).nbytes: \n", np.array([1]).nbytes)
l = list(range(1000))
s = np.arange(1000)
print ( np.sum(s ** 2))
print ( sum([x ** 2 for x in l]))