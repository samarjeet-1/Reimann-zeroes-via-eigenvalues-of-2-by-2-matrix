import numpy as np
import seaborn as sns
import numpy
import matplotlib.pyplot as plt
import mpmath as mp



#################################################################### T
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

x=int(input('enter upper limit ')) #upper limit of prime numbers you want, Maximum value used during testing
                                   #was 32,000,000 as from 33,000,000 we would run out of riemann zeroes
primes=primesfrom2to(x) #generate prime numbers


peval=[] #stores largest eigenvalue of each matrix
zzarray=np.loadtxt("D:\studies\python programs\zero",dtype=float) # loads pre-downloaded list of imaginary part of Riemann zeroes
detarray=[] # Stores determinant of all matrices
matrixarray=[] # stores the matrices
reqzzarray=[]

#########################################
for i in range(1,len(primes)-5): #Creates our matrices
    pmatrix=np.array([
      [primes[i],primes[i+1]],
       [primes[i+2],primes[i+3]]
        ])
    en,ev= np.linalg.eig(pmatrix) #calculating eigenvalues
    peval.append(en[1])
    reqzzarray.append(zzarray[i-1])
    detarray.append(np.linalg.det(pmatrix))
    matrixarray.append(pmatrix)
    
    
peval=np.array(peval) #converting list into numpy array
reqzzarray=np.array(reqzzarray)
    
    
diff_bound=0.29  # set the value of ε (epsilon)

rzdetarray=[] #determinants of matrices whose eigenvalues are riemann zeroes
rzeigarray=[] #eigenvalues w
result=[] #eigenvalues which are Riemann zeroes
compzzarray=[] #stores the true riemann zero corresponding to each eigenvalue that is also a riemann zero
zmatrices=[] #matrices whose eigenvalues are Riemann zeroes
rzposarray=[]  #stores value of i/position of each eigenvalue that is a riemann zero


def find_pairs(array1, array2): #Program that searches for eigenvalues corresponding to zeroes from the entire array of eigevalues

    i = 0
    j = 0
    while i < len(array1) and j < len(array2):
        if abs(array2[j] - array1[i]) < diff_bound:
            result.append(array2[j])
            compzzarray.append(array1[i])
            rzdetarray.append(detarray[j])
            rzposarray.append(j+1)
            zmatrices.append(matrixarray[j])
   
          
            
            i += 1
            j += 1
        elif array2[j] < array1[i]:
            j += 1
        else:
            i += 1
        
            
find_pairs(zzarray, peval) #zzarray is our reference array with pre downloaded zeroes and peval is the array from which we want to find zeroes
result=np.array(result)
compzzarray=np.array(compzzarray)
rztrarray=np.array(rztrarray)
rzposarray=np.array(rzposarray)



        

og_vs_eig_diff=np.abs(result-compzzarray) #Computes the absolute difference between eigenvalues which are considered as zeroes and true riemann zeroes
plt.hist(og_vs_eig_diff, bins=int(len(og_vs_eig_diff)**0.5), alpha=0.5, label="true vs eigenvalue") #plots their histogram
plt.legend()
plt.show()

plt.figure()
trcontriarray=np.abs(result-rztrarray) #calculates  𝛿r
plt.hist(trcontriarray, bins=int(len(og_vs_eig_diff)**0.5), alpha=0.5, label="non trace contribution") #plotting 𝛿r histogram
plt.legend()
plt.show()
plt.figure()




 


############################################
eigenvalue_diffs=np.diff(result) #calculates difference between consecutive eigenvalues that are Riemann zeroes.
plt.hist(eigenvalue_diffs, bins=int(len(eigenvalue_diffs)**0.5), alpha=0.5, label="Eigenvalues")
plt.legend()
plt.show()

ediff=np.diff(peval)
riemann_zero_diff=np.diff(reqzzarray)
compzzarray_diff=np.diff(compzzarray)
####################################################################
plt.hist(eigenvalue_diffs, bins=int(len(eigenvalue_diffs)**0.5), alpha=0.5, label="Eigenvalues")
#plt.hist(riemann_zero_diff, bins=int(len(riemann_zero_diff)**0.5), alpha=0.5, label="Riemann zeroes")
#plt.hist(ediff, bins=int(len(ediff)**0.5), alpha=0.5, label="eigenvalue diffs")
plt.hist(compzzarray_diff, bins=int(len(compzzarray_diff)**0.5), alpha=0.5, label="similar zeroes")
plt.legend()
plt.show()


################################################################
plt.figure()
# Create the KDE plot for eigenvalues
sns.kdeplot(result, label='Eigenvalues zeroes',bw=0.2,shade=True)

# Create the KDE plot for Riemann zeros
sns.kdeplot(reqzzarray, label='Riemann Zeros', bw=0.2, shade=True) 

# Set plot title and axis labels
plt.title('KDE Plot of Eigenvalues and Riemann Zeros')
plt.xlabel('Imaginary Part')
plt.ylabel('Density')

# Display the plot
################################################################
plt.figure()
num_bins = int(len(result)**0.5)
hist_range = (0, max(result))

# Calculate the histogram of the eigenvalues
hist_eig, bins_eig = np.histogram(result, bins=num_bins, range=hist_range, density=True)

# Calculate the spectral density
spectral_density = hist_eig / (2*np.pi*bins_eig[1:]*np.diff(bins_eig))

# Plot the spectral density
plt.plot(bins_eig[1:], spectral_density)
plt.xlabel('Eigenvalues')
plt.ylabel('Spectral Density')



num_bins = int(len(reqzzarray)**0.5)
hist_range = (0, max(reqzzarray))

# Calculate the histogram of the eigenvalues
hist_zero, bins_zero = np.histogram(reqzzarray, bins=num_bins, range=hist_range, density=True)

# Calculate the spectral density
spectral_density = hist_zero / (2*np.pi*bins_zero[1:]*np.diff(bins_zero))

# Plot the spectral density
plt.plot(bins_zero[1:], spectral_density)
plt.xlabel('zeroes')
plt.ylabel('Spectral Density')
plt.show()
plt.show()


plt.figure()
num_bins = int(len(eigenvalue_diffs)**0.5)
hist_range = (0, max(eigenvalue_diffs))

# Calculate the histogram of the eigenvalues
hist_eig, bins_eig = np.histogram(eigenvalue_diffs, bins=num_bins, range=hist_range, density=True)

# Calculate the spectral density
spectral_density = hist_eig / (2*np.pi*bins_eig[1:]*np.diff(bins_eig))

# Plot the spectral density
plt.plot(bins_eig[1:], spectral_density)
plt.xlabel('Eigenvalues difference')
plt.ylabel('Spectral Density')

plt.figure()

num_bins = int(len(riemann_zero_diff)**0.5)
hist_range = (0, max(riemann_zero_diff))

# Calculate the histogram of the eigenvalues
hist_zero, bins_zero = np.histogram(riemann_zero_diff, bins=num_bins, range=hist_range, density=True)

# Calculate the spectral density
spectral_density = hist_zero / (2*np.pi*bins_zero[1:]*np.diff(bins_zero))

# Plot the spectral density
plt.plot(bins_zero[1:], spectral_density)
plt.xlabel('zeroes difference')
plt.ylabel('Spectral Density')
plt.show()

################################################################
plt.figure()
import numpy as np
import matplotlib.pyplot as plt

# Generate example data
data = result

# Apply Hann window to data
window = np.hanning(len(data))
data_windowed = data * window

# Compute Fourier transform
fft = np.fft.fft(data_windowed)

# Compute frequency axis
freq = np.fft.fftfreq(len(data), d=1)

# Shift DC component to center of spectrum
fft_shifted = np.fft.fftshift(fft)
freq_shifted = np.fft.fftshift(freq)

# Plot results
plt.plot(freq_shifted, np.abs(fft_shifted))
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Generate example data
data = reqzzarray

# Apply Hann window to data
window = np.hanning(len(data))
data_windowed = data * window

# Compute Fourier transform
fft = np.fft.fft(data_windowed)

# Compute frequency axis
freq = np.fft.fftfreq(len(data), d=1)

# Shift DC component to center of spectrum
fft_shifted = np.fft.fftshift(fft)
freq_shifted = np.fft.fftshift(freq)

# Plot results
plt.plot(freq_shifted, np.abs(fft_shifted))
plt.xlabel('Frequency')
plt.ylabel('Magnitude')
plt.show()


'''plt.figure()
x=np.arange(1,len(rzposarray)+1)

y=[]

for i in x:
    func=((1/(diff_bound))*(np.log(i)**2))+(np.e*i)+((np.euler_gamma/(diff_bound))*mp.li(i))+((1/diff_bound)*(i**0.5)*np.log(i))-(3*(i**0.5))-(diff_bound*(np.log(i)**0.5))
    y.append(func)
    
plt.plot(x,rzposarray)
plt.plot(x,y)
plt.xlabel('n')
plt.ylabel('ith position of nth zero')
plt.show()
y=np.array(y)
error=np.abs(y-rzposarray)
plt.figure()
plt.plot(error)'''

plt.figure()
plt.scatter(np.arange(0,len(rzdetarray)),rzdetarray)
plt.figure()
plt.scatter(np.arange(0,len(riemann_zero_diff)),riemann_zero_diff)
plt.figure()
plt.scatter(np.arange(0,len(eigenvalue_diffs)),eigenvalue_diffs)
plt.figure()
plt.scatter(np.arange(0,len(og_vs_eig_diff)),og_vs_eig_diff)

zero=[]
for k in result:
    zero.append(mp.zeta(0.5+1j*k))
    print(mp.zeta(0.5+1j*k))



        
    





