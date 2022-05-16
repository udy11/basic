# Last updated: 07-Jun-2013
# Udaya Maurya (udaya_cbscients@yahoo.com, telegram: https://t.me/udy11)
# Source: https://github.com/udy11, https://gitlab.com/udy11
# Function to binary search nearest number in a sorted array

# ALL YOU NEED TO DO:
# Call function bsa(a, n1, n2, b) with:
#   a = numerical array in which nearest number to be searched
#   n1 = first index from which to begin search (inclusive)
#   n2 = last index upto which search is to be made (not inclusive)
#   b = number whose nearest to be searched
# NOTE that n2 is not included in the search

# OUTPUT:
# Output of the functions is the index of number to be searched

# For ascendingly sorted array
def bsan_asc(a, n1, n2, b):
    if b > a[n2 - 1]:
        return n2 - 1
    if b < a[n1]:
        return n1
    nd = n2 - n1
    nm = n1 + nd // 2
    while nd > 1:
        if a[nm] > b:
            n2 = nm
        else:
            n1 = nm
        nd = n2 - n1
        nm = n1 + nd // 2
    d1 = b - a[n1]
    d2 = a[n2] - b
    if d1 < d2:
        return n1
    else:
        return n2

# For descendingly sorted array
def bsan_dec(a, n1, n2, b):
    if b < a[n2 - 1]:
        return n2 - 1
    if b > a[n1]:
        return n1
    nd = n2 - n1
    nm = n1 + nd // 2
    while nd > 1:
        if a[nm] < b:
            n2 = nm
        else:
            n1 = nm
        nd = n2 - n1
        nm = n1 + nd // 2
    d1 = a[n1] - b
    d2 = b - a[n2]
    if d1 < d2:
        return n1
    else:
        return n2

a = [11 + i * 0.1 for i in  range(1000)]
b = 101.78

print(bsan_asc(a, 0, 1000, b))
a.reverse()
print(bsan_dec(a, 0, 1000, b))
