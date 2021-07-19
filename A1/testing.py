'''
x1 = 1
y1 = 1
print(x1,y1)
for i in range(2,100):
    x2 = (3*x1 + 4)/7.0
    y2 = (1+1/i)*x1 + 3*y1
    print(x2*y2)
    x1 = x2
    y1 = y2


'''

x0 = 0

for i in range(1,100):
    x1 = x0**2 + 0.2
    print(x1, x1**2)
    x0 = x1


'''
# a5
sgn = -1
summ = 1

for i in range(1,6):
    summ += sgn * (1.0/(i+1))
    sgn = -1*sgn
    print(summ)

a5 =summ
print()
# a10
sgn = -1
summ = 1

for i in range(1,11):
    summ += sgn * (1.0/(i+1))
    sgn = -1*sgn
    print(summ)
a10 = summ

print("asn:", a5-a10)
'''

'''
frac = 2/3.0

mult  = 10000/2

for i in range(1,20):
    curr = mult*(frac**i)/i
    print(i, curr)
    if curr < 1:
        break
'''