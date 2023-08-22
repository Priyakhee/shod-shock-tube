'''          
             
     To solve Eular Equations for shod shock tube with Riemann problem based initial condition using the Steger Warming flux vector scheme'''

import numpy as np
import matplotlib.pyplot as plt
gamma=1.4

def funcP(u1,u2,u3):  #fuction to find pressure
    p=(gamma-1)*(u3-((u2**2)/(2*u1)))
    return p

fi=np.zeros(3) 
def fluxf(u1,u2,u3,a,l1,l2,l3):  #function to find F+ or F- (flux)of a cell
    fi[0]=u1*(l1+2*(gamma-1)*l2+l3)/(2*gamma)
    fi[1]=u1*(l1*((u2/u1)-a)+(2*(gamma-1)*u2*l2/u1)+l3*((u2/u1)+a))/(gamma*2)
    fi[2]=(u1*(l1*(((u3+funcP(u1,u2,u3))/u1)-u2*a/u1)+(l2*(gamma-1)*(u2/u1)**2)+l3*(((u3+funcP(u1,u2,u3))/u1)+u2*a/u1)))/(2*gamma)
    return fi

#initializing variables and taking input
m=int(input("enter the no. equations:"))  #taking input for no.of equations
numcell=int(input("enter the no. cells:")) #taking input for number of cells
L=float(input("enter the length:"))        #taking input for domain size
t=float(input("enter the simulation time:"))#taking input for the time for which velocity profile is required
CFL=0.8         #CFL value
del_x=L/numcell #delx
time=0 #for time iteration
S=np.zeros(numcell)#to store max wave speed 
W=np.zeros((m+1,numcell)) #rho,p,u,ie 
U=np.zeros((m,numcell)) #rho,rho*u,rho*e  at nth time step
Unew=np.zeros((m,numcell))#rho,rho*u,rho*e  at (n+1)th time step
a=np.zeros(numcell)  #speed of sound 
fluxP=np.zeros((m,numcell)) #F+ flux 
fluxM=np.zeros((m,numcell)) #F- flux 
fp=np.zeros(m) #F+ flux in one cell
fm=np.zeros(m) #F- flux in one cell
f=np.zeros((m,numcell+1)) #flux at interface 
eigvalue=np.zeros(m) #eigenvalues in one cell
lambdaP=np.zeros(m) #+ve eigenvalues in one cell
lambdaM=np.zeros(m) #-ve eigenvalues in one cell
x=np.zeros(numcell) #x values from 0 to L
ana=np.zeros((m+2,numcell)) # to store analytical results for comparison

# reading analytical data
filename = 'analytical.csv'
data = np.loadtxt(filename, delimiter=',', dtype=float)
for j in range(numcell): 
    for i in range(m+2):
           ana[i][j]=data[j][i]

#calculating x
x[0]=del_x/2
for i in range (1,numcell):
    x[i]=x[i-1]+del_x
    
#initial value input     
WoL= [] #left side values
print("Enter the intial WoL entries:")
for i in range(m):       
    WoL.append(float(input()))
WoR= [] #right side values
print("Enter the intial WoR entries:")
for i in range(m):          
    WoR.append(float(input()))
    
#writing inital values
for i in range(m):
    for j in range(numcell):
        if (j<(L/2)/del_x):
            W[i][j]=WoL[i]
        else:
            W[i][j]=WoR[i]            
for i in range (numcell):   #values of U at 0th time step
    U[0][i]=W[0][i]
    U[1][i]=W[0][i]*W[2][i]
    U[2][i]=((W[0][i]*(W[2][i])**2)/2)+W[1][i]/(gamma-1)
 
while (time<t): #to iterate time step
    for j in range(numcell):
         a[j]=(gamma*funcP(U[0][j],U[1][j],U[2][j])/U[0][j])**.5 #calculating a in cells
         for i in range (-1,2,1):
             eigvalue[i+1]=(U[1][j]/U[0][j])+i*a[j] #calculating eiganvalues in a cell
           
         for i in range (m): #finding lambdaP &lambdaM
             lambdaP[i]=(eigvalue[i]+abs(eigvalue[i]))/2 
             lambdaM[i]=(eigvalue[i]-abs(eigvalue[i]))/2
        
         S[j]=max(abs(eigvalue)) #storing maximum eigenvalues of each cells
         
         fp=fluxf(U[0][j],U[1][j],U[2][j],a[j],lambdaP[0],lambdaP[1],lambdaP[2]) #finding F+ in a cell
         for i in range (m):
            fluxP[i][j]=fp[i]  #storing F+ values 
            
         fm=fluxf(U[0][j],U[1][j],U[2][j],a[j],lambdaM[0],lambdaM[1],lambdaM[2]) #finding F- in a cell
         for i in range (m):
            fluxM[i][j]=fm[i]  #storing F- values

    del_t=CFL*del_x/max(S)  #calculating delta t
    
    #calculating flux at interface
    for j in range (numcell-1):
         for i in range (m):
            f[i][j+1]=fluxP[i][j]+fluxM[i][j+1] 
    for i in range (m):
         f[i][0]=f[i][1]
         f[i][numcell]=f[i][numcell-1]

    #calculating u for n+1 th time step in a cell   
    for j in range(numcell):
         for i in range (m): 
            Unew[i][j]=U[i][j]-del_t*(f[i][j+1]-f[i][j])/del_x  
       
    for j in range(numcell): #updating u values
        for i in range(m):
            U[i][j]=Unew[i][j]  
            
    time=time+del_t #increasing time step

for i in range (numcell): #converting U to W(rho,P,U,ie)
    W[0][i]=Unew[0][i]
    W[2][i]=Unew[1][i]/Unew[0][i]
    W[1][i]=(gamma-1)*(Unew[2][i]-(W[0][i]*(W[2][i]**2)/2))
    W[3][i]=W[1][i]/((gamma-1)*W[0][i])

#plots
#plotting final velocity profiles along x 
plt.plot(x[range(numcell)],W[2][range(numcell)],label='U steger warming',linestyle='',marker='o' ,markersize=2,color='black') 
plt.plot(x[range(numcell)],ana[2][range(numcell)],label='analytical')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('velocity')
plt.title('velocity Vs X ')
plt.show()
#plotting final density profiles along x
plt.plot(x[range(numcell)],W[0][range(numcell)],label='rho steger warming',linestyle='',marker='o' ,markersize=1.5,color='black' ) 
plt.plot(x[range(numcell)],ana[1][range(numcell)],label='analytical')
plt.legend()
plt.xlabel('x')
plt.ylabel('density')
plt.title('density Vs X ')
plt.show()
#plotting final pressure profiles along x
plt.plot(x[range(numcell)],W[1][range(numcell)],label='P steger warming',linestyle='',marker='o' ,markersize=1.5,color='black' ) 
plt.plot(x[range(numcell)],ana[3][range(numcell)],label='analytical ')
plt.legend()
plt.xlabel('x')
plt.ylabel('pressure')
plt.title('pressure Vs X ')
plt.show()
#plotting final internal energy profiles  along x
plt.plot(x[range(numcell)],W[1][range(numcell)]/(.4*W[0][range(numcell)]),label='ie steger warming',linestyle='',marker='o' ,markersize=2,color='red' )       
plt.plot(x[range(numcell)],ana[4][range(numcell)],label='analytical')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('internal energy')
plt.title('internal energy Vs X ')
plt.show()
                                           
             
                 
                 
                     
