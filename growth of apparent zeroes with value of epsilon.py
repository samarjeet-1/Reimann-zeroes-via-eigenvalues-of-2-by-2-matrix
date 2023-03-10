import numpy as np
import numpy
import matplotlib.pyplot as plt




####################################################################
def primesfrom2to(n):
    """ Input n>=6, Returns a array of primes, 2 <= p < n """
    sieve = numpy.ones(n//3 + (n%6==2), dtype=bool)
    for i in range(1,int(n**0.5)//3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[       k*k//3     ::2*k] = False
            sieve[k*(k-2*(i&1)+4)//3::2*k] = False
    return numpy.r_[2,3,((3*numpy.nonzero(sieve)[0][1:]+1)|1)]

######################################################################

x=int(input('enter upper limit '))
primes=primesfrom2to(x)



peval=[]
zzarray=np.loadtxt("D:\studies\python programs\zero",dtype=float)

#n=smp.symbols('n')
#########################################
#matrix with prime entries that takes in 4 numbers and gives eigenvalues
for i in range(1,len(primes)-5):
    pmatrix=np.array([
      [primes[i],primes[i+1]],
       [primes[i+2],primes[i+3]]
        ])
    en,ev= np.linalg.eig(pmatrix)
    peval.append(en[1])
    
    
peval=np.array(peval)

    
    
result=[]


result_len=[]
diff_boundval=[]


def growth (array1,array2,diff_bound):
    while diff_bound<=0.5:
     i = 0
     j = 0
     result=[]
     while i < len(array1) and j < len(array2):
         if abs(array2[j] - array1[i]) < diff_bound:
            result.append(array2[j])
            i += 1
            j += 1
         elif array2[j] < array1[i]:
            j += 1
         else:
            i += 1
     result_len.append(len(result))
     diff_boundval.append(diff_bound)
     diff_bound=diff_bound+0.01
    
            
growth(zzarray, peval,0)

plt.xlabel("ε")
plt.ylabel("Fz(ϵ)")
plt.plot(diff_boundval,result_len)