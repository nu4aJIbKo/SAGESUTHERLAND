#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def HEC_random_point(C):
    f = C.hyperelliptic_polynomials()[0]
    while True:
        x_r = f.base_ring().random_element()
        y2 = f(x_r)
        if y2.is_square():
            return [x_r, y2.sqrt()]
def HEC_random_points_uniq(C, n):
    f = C.hyperelliptic_polynomials()[0]
    res = []
    for i in range(1,n+1):
        tries = 100
        found = False
        while tries > 0 and (not found):
            P = HEC_random_point(C)
            if not (P in res):
                res.append(P)
                found = True
            tries = tries - 1
    return res
def JC_random_element(C):
    f = C.hyperelliptic_polynomials()[0]
    R.<x> = f.base_ring()['x']
    J = C.jacobian()
    points = HEC_random_points_uniq(C, C.genus())
    u = 1
    for point in points:
        u = u * (x-point[0])
    v = f.parent().lagrange_polynomial(points)
    print('u,v= ', [u,v])
    return J([u, v])
q=433
K = GF(q)
R.<x> = K[]
f = x**7+x**6+x**5+x**3+x**2+x
H = HyperellipticCurve(f); H
C = HyperellipticCurve(f)
D=JC_random_element(C)
print('D= ', D)
B_j={}
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
          73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 
          151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
          229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
          311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389,
          397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
          479, 487, 491, 499, 503, 509, 521, 523, 541]
g=3
M = floor(((sqrt(q) + 1)**2)**g)
w=1
Pw=2
phi_pw=1
pw=2
def primorials(primes, w):
    primes = primes[:w]
    pw = primes[-1]
    Pw = prod(primes)
    phi_pw = prod(p - 1 for p in primes)
    return pw, Pw, phi_pw
primorials(primes, 1)
while Pw <= (sqrt(M)):
    w=w+1
    primorials(primes, w)
    pw, Pw, phi_pw = primorials(primes, w)
w=w-1
primorials(primes, w)
pw, Pw, phi_pw = primorials(primes, w)
m=floor(sqrt(M//(Pw*phi_pw)))
b=m*Pw
print('b- ', b)
primes=primes[:w]
print('primes- ', primes)
E = prod([p^floor(log(M, p)) for p in primes])
beta=E*D
print('beta- ', beta)
def A(ND,E,primes,w):
    primes=primes[:w]
    P=1
    j=1
    for p in primes:
        D_i=(E//prod([p^floor(log(M, p))]))*E*D
        print('E*D - ', E*D)
        print('D_i - ', D_i)
        while D_i != 0:
            D_i=p*D_i
            j=j+1
        P=P*p**j
    return P
#Шаг младенца
for j in range(0,b):
    if gcd(j,Pw):
        B=j*beta
    B_j[j]=B
    if B==0:
        N=j
        print('N- ', N)
        h=A(N*D,E,primes,w)
        print('h- ', h)
        print('h*N-', h*N)
J=0
I=1
i=2
while I !=J:
    B_i=i*b*D
    for j in range(0,b):
        if B_i==B_j:
            N=i*b-j
            h=A(ND,E,primes,w)
            print('h*N-', h*N)
            J=I
            j=b
        break
    i=i+1
print(factor(E))
print(pw)
print(Pw)
print(phi_pw)
print(m)


# In[ ]:




