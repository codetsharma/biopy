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

print ([x**5 for x in range(19)])
print (set ([x**2 for x in range(-10,11)]))