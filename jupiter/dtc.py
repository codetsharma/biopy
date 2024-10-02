#
# dataTypeComprehensions
#

# ~~~~~~~~~~~~~~~~~~~~~ #
# List Comprehension:   #
# ~~~~~~~~~~~~~~~~~~~~~ #
# Subst of another 'static' or user-driven list:
#
# Possible Application: 
# - Filter out 1 aspect of data, 
# "[ func (varName) for <varName> in <Iterable> if <condition> "]"
print ("\nList Comprehensions: \n")
print ("\nprint ([x**5 for x in range(19)]): \n", [x**5 for x in range(19)])
print ("\nset ([pow(x,19) in range(-10,11)]): \n", set ([pow(x, 3) for x in range(-10,11) if x%3==0]))

# ~~~~~~~~~~~~~~~~~~~~~ #
# Dict Comprehension:   #
# ~~~~~~~~~~~~~~~~~~~~~ #
print ("\nDict Comprehensions: \n")
print ("\nSquare:\n", {x:x**2 for x in range(10)})
print ("\nCube\n", {x:x**3 for x in range(10)})

# ~~~~~~~~~~~~~~~~~~~~~ #
# Set Comprehension:    #
# ~~~~~~~~~~~~~~~~~~~~~ #
a = set ([x**2 for x in range(-10, 11)])
print ("\na: ", a)

# Assignment expressions
x = 30 + 40
y = x * 5
print(x, y)

y = (x:=30+40)*5
print(x, y)

# Object oriented programming, an expression that makes a class
C = type('C', (), {'x': 10})
print (C.x)

# Function expressions
cube = lambda x: x ** 3
print (cube(5))

#
# JSON parsing:
#
import json
with open("books.json", "r") as catalog_file:
    catalog = json.load(catalog_file)

print ("\nTask 1: How many books in catalog? \n")
print (len(catalog))
for book in catalog.values():
    print (book['description'])


