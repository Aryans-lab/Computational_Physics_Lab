#Problem: Sum of first 20 even numbers and the factorial of 8 using for loop or any other loop
#Name: Aryan Bandyopadhyay, Roll number: 2411014
def sum_of_first_20_even_numbers():
    sum_even = 0
    for i in range(1,21):
        sum_even = sum_even + 2*i
    return sum_even

def factorial_of_8():
    fact = 1
    for i in range(1,9):
        fact =fact*i
    return fact
#Now I want output for both the functions
print("Sum of first 20 even numbers:", sum_of_first_20_even_numbers())
print("Factorial of 8:", factorial_of_8())

# End of my code


#Output:
########################################################
# Sum of first 20 even numbers: 420
# Factorial of 8: 40320
########################################################