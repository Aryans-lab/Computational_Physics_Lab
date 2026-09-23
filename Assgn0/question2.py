#Problem: Calculate the sum of N = 15 terms of a GP and HP series for common difference
#1.5 and common ratio 0.5 starting from t0 = 1.25. Use of analytical formulae is
#not allowed.
#Name: Aryan Bandyopadhyay, Roll number: 2411014
n=15
t0=1.25 
r=0.5
d=1.5

def sum_of_gp_series(n, t0, r):
    sum_gp = 0
    term = t0
    for i in range(n):
        sum_gp += term
        term *= r
    return sum_gp

def sum_of_hp_series(n, t0, d):
    sum_hp = 0
    term = t0
    for i in range(n):
        sum_hp += 1/term
        term += d
    return sum_hp

print("Sum of GP series:", sum_of_gp_series(n, t0, r))
print("Sum of HP series:", sum_of_hp_series(n, t0, d))

# End of my code


#Output:
########################################################
# Sum of GP series: 2.4999237060546875
# Sum of HP series: 2.4139570733659186
########################################################