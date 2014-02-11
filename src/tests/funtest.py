'''
Created on Jan 9, 2014

@author: jwe

def return_a_whole_new_string(the_string):
    new_string = something_to_do_with_the_old_string(the_string)
    return new_string

# then you could call it like
my_string = return_a_whole_new_string(my_string)
'''
def fun(a, b=None):
    a = 1
    b = 2
    return 3

x = 10
y = 20

z = fun(x, b=y)
print x, y, z