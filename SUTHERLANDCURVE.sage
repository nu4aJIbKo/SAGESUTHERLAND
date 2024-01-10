import time
start = time.time()
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
q = Primes().next(2^(25)+24)
#q = Primes().next(2^(20)+30)
q = Primes().next(2^(30)+24)
q =33554467
print('q= ' ,q)
K = GF(q)
R.<x> = K[]
f = x**7+x**4+x
C = HyperellipticCurve(f)
print(C)
D=JC_random_element(C)
J = C.jacobian()
print('D= ', D)
primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
          73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 
          151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
          229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
          311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389,
          397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
          479, 487, 491, 499, 503, 509, 521, 523, 541]
g=3
M1=floor(((sqrt(q) + 1)**2)**g)
M=2^31
w=1
Pw=2
phi_pw=1
pw=2
#Функция праймориал
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
print('Pw-',Pw)
b=m*Pw
print('b- ', b)
#print('Pw-' ,Pw)
primes=primes[:w]
#print('primes- ', primes)
E = prod([p^floor(log(M, p)) for p in primes])
print("E=", E,factor(E))
beta=E*D
print('beta- ', beta)
#Классический алгоритм A(AN,E)
def A(ND,E,primes,w):
    primes=primes[:w]
    S_i=1
    S=1
    for p in primes:
        D_i=(E//prod([p^floor(log(M, p))]))*ND
        #print('D_i-',D_i)
        #print(factor(E//prod([p^floor(log(M, p))])))
        while list(D_i) != list(C.jacobian()(0)):
            D_i=p*D_i
            #print('D_i',D_i)
            S_i=S_i*p
        S=S*S_i
        S_i=1
    print("S=",S)    
    return(S)
#Шаг младенца
print("Шаг_младенца")
d = {}
for i in range(0,b):
    d[tuple(D)]=0
for j in range(0,b):
    if gcd(j,Pw)==1:
        X=j*beta
        #print('X=',X)
        d[tuple(X)]=j
        #print(d)

if tuple(J(0)) in d:
    j = d[tuple(J(0))]
    N = j
    h=A(N*D,E,primes,w)
    print('h- ', h)
    print('N- ',N)
    print('порядок- ', h*N)
    print('Найден в Шаг_младенца за',time.time() - start)
#Шаг гиганта   
print("Шаг_Гиганта")
G=0
I=1
i=2
j=1
while I !=G:
    B_i=i*b*beta
    #print("B_i=",B_i)
    if tuple(B_i) in d:
        j = d[tuple(B_i)]
        N = i*b-j
        #print("j=",j)
        #print("i*b=",i*b)
        print('N- ',N)
        h=A(N*D,E,primes,w)
        #print('h- ', h)
        print('порядок- ', h*N)
        print('Найден в Шаг_гиганта за',time.time() - start)
        G=I
        break
    i=i+1
    #print(i)
print(time.time() - start)
if list(J(0))==list(D*h*N):
    print("true")
print("end__end")















