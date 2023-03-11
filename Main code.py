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
trarray=[] #trace of marices
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
    trarray.append(np.trace(pmatrix))
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
rztrarray=[] #trace of matrices whose eigenvalues are riemann zeroes


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
            rztrarray.append(trarray[j])
   
          
            
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
rztrarray=np.array(rztrarray)



        

og_vs_eig_diff=np.abs(result-compzzarray) #Computes the absolute difference between eigenvalues which are considered as zeroes and true riemann zeroes
plt.xlabel("Difference")
plt.ylabel("frequency")
plt.grid()
plt.hist(og_vs_eig_diff, bins=int(len(og_vs_eig_diff)**0.5), alpha=0.5, label="Difference between eigenvlues and true zeroes") #plots their histogram
plt.legend()
plt.show()

plt.figure()
trcontriarray=np.abs(result-rztrarray) #calculates  𝛿r
plt.grid()
plt.xlabel("delta r")
plt.ylabel("frequency")
plt.hist(trcontriarray, bins=int(len(og_vs_eig_diff)**0.5), alpha=0.5, label="Non-Trace contribution to eigenvalues") #plotting 𝛿r histogram
plt.legend()
plt.show()
plt.figure()




 


############################################
eigenvalue_diffs=np.diff(result) #calculates difference between consecutive eigenvalues that are Riemann zeroes.
compzzarray_diff=np.diff(compzzarray) ##calculates difference between corresponding riemann zeroes of eigenvalues that are considered riemann zeroes
plt.xlabel("Spacing between Consecutive eigenvalues corresponding to Reimann zeroes")
plt.ylabel("frequency")
plt.grid()
plt.hist(eigenvalue_diffs, bins=int(len(eigenvalue_diffs)**0.5), alpha=0.5, label="Spacings between Eigenvalues that are Riemann Zeroes")
#plt.hist(compzzarray_diff, bins=int(len(compzzarray_diff)**0.5), alpha=0.5, label="Corresponding True zeroes")
plt.legend()
plt.show()
plt.figure()

ediff=np.diff(peval) # Calculates Difference between all the consecutive eigenvalues
plt.xlabel("Spacing between Consecutive eigenvalues")
plt.ylabel("frequency")
plt.grid()
plt.hist(ediff, bins=int(len(ediff)**0.5), alpha=0.5, label="Spacings between all the eigenvalues")
plt.legend()
plt.show()
plt.figure()

riemann_zero_diff=np.diff(zzarray) #Calculates difference between consective Riemann zeroes
plt.xlabel("Non-normalized spacings between Riemann zeroes ")
plt.ylabel("frequency")
plt.grid()
plt.hist(riemann_zero_diff, bins=int(len(riemann_zero_diff)**0.5), alpha=0.5, label="Non-normalized spacings between Riemann zeroes")
plt.legend()
plt.show()
plt.figure()


################################################################
# Create the KDE plot for eigenvalues
plt.grid()
sns.kdeplot(result, label='Eigenvalue zeroes',bw=0.2,shade=True)

# Create the KDE plot for Riemann zeros
sns.kdeplot(zzarray, label='Riemann Zeros', bw=0.2, shade=True) 

plt.title('KDE Plot of Eigenvalues and Riemann Zeros')
plt.xlabel('Zeroes')
plt.ylabel('Density')
plt.legend()
plt.show()
################################################################
plt.figure()
num_bins = int(len(result)**0.5)
hist_range = (0, max(result))

# Calculate the histogram of the eigenvalues which are zeroes
hist_eig, bins_eig = np.histogram(result, bins=num_bins, range=hist_range, density=True)

# Calculate the spectral density
spectral_density = hist_eig / (2*np.pi*bins_eig[1:]*np.diff(bins_eig))

# Plot the spectral density
plt.plot(bins_eig[1:], spectral_density,label="Eigenvalues that are Riemann zeroes")




num_bins = int(len(zzarray)**0.5)
hist_range = (0, max(zzarray))

# Calculate the histogram of reimann zeroes
hist_zero, bins_zero = np.histogram(zzarray, bins=num_bins, range=hist_range, density=True)

# Calculate the spectral density
spectral_density = hist_zero / (2*np.pi*bins_zero[1:]*np.diff(bins_zero))

# Plot the spectral density
plt.plot(bins_zero[1:], spectral_density,label="Riemann Zeroes")
plt.xlabel('Frequency')
plt.ylabel('Spectral Density')
plt.legend()
plt.grid()
plt.show()
plt.figure()


################################################################ Fourier transform
plt.figure()
import numpy as np
import matplotlib.pyplot as plt

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


data =zzarray

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
plt.grid()
plt.legend()
plt.show()
plt.figure()


############################################ Equations for positions i in the sequence of eigenvalues for the nth eigenvalue corresponding to a riemann zero
x=np.arange(1,len(rzposarray)+1)
y=[]

# below is the equation to plot Predicted i position of the nth zero for epsilon=0.29
for i in x:
    func=(np.log(i)**2)+mp.li(i)-((diff_bound*i*np.log(i))**0.5)+(2*(i**0.5)*np.log(i))+i
    y.append(func)

    
# below is the equation to plot Predicted i position of the nth zero for epsilon=0.2, comment out the above code for 0.29 and then remove comment of this code
#for i in x:
#    func=(1/diff_bound)*((np.log(i)**2)+mp.li(i)-((i*np.log(i))**0.5))+(2*(i**0.5)*np.log(i))+i
#    y.append(func)
    
    
plt.plot(x,rzposarray,label="True i position of the nth zero")
plt.plot(x,y, label="Predicted i position of the nth zero")
plt.xlabel('n')
plt.ylabel('i position of nth zero')
plt.grid()
plt.legend()
plt.show()
y=np.array(y)
error=np.abs(y-rzposarray)
plt.figure()
plt.xlabel("n")
plt.ylabel("Absolute error")
plt.grid()
plt.plot(error)
plt.legend()
plt.show()
plt.figure()

###########################################
plt.xlabel("i")
plt.ylabel("Determinant of ith matrix")
plt.grid()
plt.hist(detarray, bins=int(len(detarray)**0.5), alpha=0.5, label="determinants of all the matrices")
plt.legend()
plt.show()
