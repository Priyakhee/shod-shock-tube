'''              shod shock tube problem
             
     To solve Eular Equations with Riemann problem based initial condition using the Roe solver with entropy fix'''

import numpy as np
import matplotlib.pyplot as plt
gamma=1.4

fi=np.zeros(3)
def fluxf(u1,u2,u3):#function to calculate flux in a cell
    fi[0]=u1*u3
    fi[1]=u2+u1*u3**2
    fi[2]=(u1*(u3**3)/2)+(gamma*u2*u3)/(gamma-1)
    return fi

def funcH(rho,u,p): #function to calculate H
    h=((u**2)/2)+(p*gamma)/(rho*(gamma-1))
    return h

#initializing variables and taking input
m=int(input("enter the no. equations:"))  #taking input for no.of equations
numcell=int(input("enter the no. cells:")) #taking input for number of cells
L=float(input("enter the length:"))        #taking input for domain size
t=float(input("enter the simulation time:"))#taking input for the time for which velocity profile is required
CFL=0.8         #CFL value
del_x=L/numcell #size of one cell
time=0.0  #for time iteration
e=0.000001 #entropy fix criteria
W=np.zeros((m,numcell)) #rho,P,u
U=np.zeros((m,numcell)) #rho,rho*u,rho*e  at nth time step
Unew=np.zeros((m,numcell)) #rho,rho*u,rho*e  at (n+1)th time step
flux=np.zeros(m)  #F flux in one cell
f=np.zeros((m,numcell+1)) #flux at inter face
eigvec=np.zeros((m,m)) #eigenvectors in a cell
eigval=np.zeros(m) #eigenvalues in one cell
delta=np.zeros(m)
mxm=np.zeros(numcell) #to store max wave speed 
x=np.zeros(numcell)   #x values from 0 to L
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
    
#inital values
for i in range(m):
    for j in range(numcell):
        if (j<(L/2)/del_x):
            W[i][j]=WoL[i]
        else:
            W[i][j]=WoR[i]
for i in range (numcell):  #values of U at 0th time step 
    U[0][i]=W[0][i]
    U[1][i]=W[0][i]*W[2][i]
    U[2][i]=((W[0][i]*(W[2][i])**2)/2)+W[1][i]/(gamma-1)
 
while (time<t):  # to iterate time step
    for j in range(numcell-1):
         #calculating u_tilda,H_tilda,a_tilda,du1,du2,du3,del1,del2,del3
         u_tild=((W[0][j+1]**.5)*W[2][j+1] +(W[0][j]**.5)*W[2][j])/((W[0][j+1]**.5)+(W[0][j]**.5))
         H_tild=((W[0][j+1]**.5)*funcH(W[0][j+1],W[2][j+1],W[1][j+1]) + (W[0][j]**.5)*funcH(W[0][j],W[2][j],W[1][j]))/((W[0][j+1]**.5)+(W[0][j]**.5))
         a_tild=((gamma-1)*(H_tild-(u_tild**2)/2))
         du1=W[0][j+1]-W[0][j]
         du2=W[0][j+1]*W[2][j+1]-W[0][j]*W[2][j]
         du3=((W[0][j+1]*(W[2][j+1])**2)/2)+(W[1][j+1]/(gamma-1))-((W[0][j]*(W[2][j])**2)/2)-(W[1][j]/(gamma-1))
         delta[1]=(gamma-1)*(du1*(H_tild-u_tild**2)+u_tild*du2-du3)/(a_tild**2)
         delta[0]=(du1*(u_tild+a_tild)-du2-a_tild*delta[1])/(2*a_tild)
         delta[2]=du1-delta[0]-delta[1]

         #calculating eiganvalues and eigen vectors in a cell
         for i in range (-1,2,1): 
             eigval[i+1]=u_tild+i*a_tild
             eigvec[i+1][0]=1
             eigvec[i+1][1]=u_tild+i*a_tild
             eigvec[i+1][2]=H_tild+i*u_tild*a_tild
         eigvec[1][2]=(u_tild**2)/2

         #applying entropy fix
         for i in range (-1,2,2):
             if(eigval[i]<e):
                 lam=((eigval[i]**2/e)+e)/2
             elif(eigval[i]>e):
                 lam=eigval[i]
             eigval[i]=lam
         mxm[j]=max(abs(eigval))  

         #calculating flux 
         flux=fluxf(W[0][j],W[1][j],W[2][j])
         for i in range (m):
             summ=0
             for k in range (m):
                 if(eigval[k]<0):
                    summ=summ+ eigval[k]*delta[k]*eigvec[k][i]
             f[i][j+1]=flux[i]+summ
             
    # to find maximum eigen value at last cell 
    u_tild=((W[0][numcell-1]**.5)*W[2][numcell-1] +(W[0][numcell-2]**.5)*W[2][numcell-2])/((W[0][numcell-1]**.5)+(W[0][numcell-2]**.5))
    H_tild=((W[0][numcell-1]**.5)*funcH(W[0][numcell-1],W[2][numcell-1],W[1][numcell-1]) + (W[0][numcell-2]**.5)*funcH(W[0][numcell-2],W[2][numcell-2],W[1][numcell-2]))/((W[0][numcell-1]**.5)+(W[0][numcell-2]**.5))
    a_tild=((gamma-1)*(H_tild-(u_tild**2)/2))
    for i in range (-1,2,1):
             eigval[i+1]=u_tild+i*a_tild
    for i in range (-1,2,2):
             if(eigval[i]<e):
                 lam=((eigval[i]**2/e)+e)/2
             elif(eigval[i]>e):
                 lam=eigval[i]
             eigval[i]=lam
    mxm[numcell-1]=max(abs(eigval)) #storing maximum eigenvalue of last cell

    for i in range (m): #Boundary Condition update
        f[i][0]=f[i][1]
        f[i][numcell]=f[i][numcell-1]

    del_t=CFL*del_x/max(mxm) #calculating delta t
    
    for j in range(numcell):
         for i in range (m):
            Unew[i][j]=U[i][j]-del_t*(f[i][j+1]-f[i][j])/del_x  #calculating u for n+1 th time step in a cell
  
    for j in range(numcell):
        for i in range(m):
            U[i][j]=Unew[i][j] #updating u values
    
    for i in range (numcell): #converting U to W(rho,P,U)
        W[0][i]=Unew[0][i]
        W[2][i]=Unew[1][i]/Unew[0][i]
        W[1][i]=(gamma-1)*(Unew[2][i]-(W[0][i]*(W[2][i]**2)/2))
 
    time=time+del_t #increasing time step

#plots        
#final velocity Vs axial direction plot
plt.plot(x[range(numcell)],W[2][range(numcell)],label='U roes entropy fix',linestyle='',marker='o' ,markersize=2,color='black') 
plt.plot(x[range(numcell)],ana[2][range(numcell)],label='analytical')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('velocity')
plt.title('velocity Vs X ')
plt.show()
#final dansity Vs axial direction plot
plt.plot(x[range(numcell)],W[0][range(numcell)],label='rho roes entropy fix',linestyle='',marker='o' ,markersize=1.5,color='black' ) 
plt.plot(x[range(numcell)],ana[1][range(numcell)],label='analytical')
plt.legend()
plt.xlabel('x')
plt.ylabel('density')
plt.title('density Vs X ')
plt.show()
#final pressure Vs axial direction plot
plt.plot(x[range(numcell)],W[1][range(numcell)],label='P roes entropy fix',linestyle='',marker='o' ,markersize=1.5,color='black' )
plt.plot(x[range(numcell)],ana[3][range(numcell)],label='analytical')
plt.legend()
plt.xlabel('x')
plt.ylabel('pressure')
plt.title('pressure Vs X ')
plt.show()
#final internal energy Vs axial direction plot
plt.plot(x[range(numcell)],W[1][range(numcell)]/(.4*W[0][range(numcell)]),label='ie roes entropy fix',linestyle='',marker='o' ,markersize=2,color='black' )         
plt.plot(x[range(numcell)],ana[4][range(numcell)],label='analytical')
plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('internal energy')
plt.title('internal energy Vs X ')
plt.show()
                       
                                           
             
                 
                 
                     
                
             
             
         
