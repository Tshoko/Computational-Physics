# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 07:09:56 2020

@author: boitshoko
"""
import numpy as np
import matplotlib.pyplot as plt

import json

""""Quetion 1

       a
"""
def boundry(u,k=0):
    u[::,0]=k
    u[0,::]=k
    u[len(u)-1,::]=k
    u[::,len(u)-1]=k
    u[20:80,40:41]=k
    
    return u


def relax(u,p,e,h):
    niter = 0
    diff = 1000
    change=0
    for i in range(0,10000):
      
        Vold=u
        #establish V matrices shifted in each direction needed for vectorization
        #print(u,"new u")
        #print(u,"down")
        uup=np.roll(u, -1, axis=0) # up
        udown=np.roll(u, 1, axis=0) # down
        
        #print(u,"up")
        uright=np.roll(u, 1, axis=1) # right
        
        #print(u,"right")
        uleft=np.roll(u, -1, axis=1) # left
        #print(u,"left")
        #This averaging step solves Laplace's equation
        u=(udown+uup+uright+uleft+p)/4
        #print(u,"full")
        
        u=boundry(u,0)
        
        
        #print(u)
        
        #maintain the boundary condition
        #track the change in the matrix at each step to know when to quit
        deltaV=u-Vold
        change=abs(np.sum(deltaV)/(len(u)*len(u)))
        #conditional to quit the method when changes are within tolerance
        if (i>100000 and change<0.1):
            
            break
    print(change)
    return u

L=100

# Set array size and set the interior value with Tguess
u = np.empty((L, L))
uguess=0

u.fill(uguess)
# let's put a source just off centre
R = np.empty((L, L))
R.fill(uguess)
e=8.85*10**(-12)
R[50,60]=1
       
h =1
e=1
u=relax(u,R,e, h)
boundry(u,0)


print(u[10,90])

print(u[50,50])

print(u[50,10])

fig,ax=plt.subplots(1,1)
ylist = np.linspace(0, 10, 100)
ylist=ylist[::-1]
xlist = np.linspace(0, 10, 100)
colorinterpolation = 50
colourMap = plt.cm.jet
X, Y = np.meshgrid(xlist, ylist)
# Configure the contour
plt.title("Potential")
plt.contourf(X, Y, u, colorinterpolation, cmap=colourMap)
plt.xlabel('x (cm)')
plt.ylabel('y (cm)')
# Set Colorbar
plt.colorbar()
# Show the result in the plot window
plt.show()
plt.imshow(u)
plt.xlabel('x (cm)')
plt.ylabel('y (cm)')


""""Quetion 1

       B
"""




def relaxe_E(u):
    Ex=-(np.roll(u, -1, axis=1) -u)
    Ey=-(np.roll(u, 1, axis=0)-u)
    for i in range(0,10000):
        Vold=u
        
        
        #print(u,"up")
        Exright=np.roll(Ex, -1, axis=1) # right
        
        #print(u,"right")
        Eydown=np.roll(u, 1, axis=0) # left
        #print(u,"left")
        #This averaging step solves Laplace's equation
        Ex1=(Exright-Ex)/2
        Ey1=(Eydown-Ey)/2
        
    return (Ex1,Ey1)
Ey,Ex=np.gradient(np.array(-u),2/100)


fig, ax = plt.subplots(figsize=(9,9))
        

ax.streamplot(X,Y,Ex,Ey,density=3)

ax.set_aspect('equal')
ax.plot(-1,0,'-or')
ax.plot(1,0,'-ob')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_title('Stream Plot of a grounded Plate and a Point Charges')
plt.show()
#xv, yv = np.meshgrid(X, Y, indexing='ij', sparse=False)

fig = plt.figure(11)
ax = fig.gca()
ax.quiver(X,Y,Ex,Ey)
plt.axis('equal')
plt.show()
Ex,Ey=relaxe_E(u)
fig = plt.figure(11)
ax = fig.gca()
ax.quiver(X,Y,Ex,Ey)
plt.axis('equal')
plt.xlabel('x (cm)')
plt.ylabel('y (cm)')
plt.title('Stream Plot of a grounded Plate and a Point Charges')
plt.show()
dx =10
Ey,Ex=np.gradient(np.array(-u),dx)

mag = np.sqrt(Ex.dot(Ey))
dot_product = np. dot(Ex/mag, Ey/mag)
angle = np. arccos(Ey/mag)
angle2=np.arctan(Ey/Ex)
print(mag[20,80],np.degrees(angle2[10,90]))

Ey,Ex=np.gradient(np.array(-u),dx)


""""Quetion 1

       C
"""

def density(y):
    q=1
    d=20
    return q*d/(2*np.pi*(d**2+y**2)**(3/2))
den=[]
Y=np.zeros(60)
k=20
j=0
for i in range(-30,30):
    den.append(density(i))
    Y[j]=k
    k+=1
    j+=1
    
den=np.array(den)
print(np.sum(den))
plt.plot(Y,den)
plt.xlabel("y(mm)")
plt.ylabel("Charge densisty IN C/mm^3")
plt.show()

p=mag[21:80,42]-mag[21:80,41]
Y=np.zeros(59)
k=20
j=0
for i in range(-29,30):
    Y[j]=k
    k+=1
    j+=1
print(np.sum(p))
plt.plot(Y,p)
plt.xlabel("y(mm)")
plt.ylabel("Charge densisty IN C/mm^3")
plt.plot()
plt.show()



""" question 2 """

L=6
V=np.zeros(L)
V[3]=7.50
A = np.zeros([L,L])
A=A.tolist()

def Initial(A,rx):
    """Make the matrix for the curcuit"""
    r1=1.50
    r2=63
    r3=rx
    r4=3
    r5=240
    r6=190
    
    A[0][0]=1
    A[0][1]=-1
    A[0][2]=-1
    A[1][1]=1
    A[1][3]=1
    A[1][4]=-1
    A[2][2]=1
    A[2][3]=-1
    
    A[2][5]=-1
    A[3][0]=r1
    A[3][1]=r2
    A[3][4]=r5
    A[4][1]=-r2
    A[4][2]=r3
    A[4][3]=r4
    A[5][3]=r4
    A[5][4]=r5
    A[5][5]=-r6
    return A
def Current(A,b,L):
    """ Solve for different values of Rx"""
    Irx=[]
    for i in range(0,L):
        A=np.array(Initial(A,i))
        I = np.linalg.solve(A, V)
        Irx.append(I[b])
    return Irx


A=np.array(Initial(A,501))
I = np.linalg.solve(A, V)    
A=np.array(Initial(A,500))
I = np.linalg.solve(A, V)
print("For Rx=500,IA is")
print(I[3])
plt.plot(Current(A,3,500))
plt.xlabel("Rx(Ohms)")
plt.ylabel("I(Ampres)")

Rx=np.linspace(45, 55, 100000)
I=Current(A,3,500)
k=0
b=50
for i in I:
    if(abs(i-0.00)<0.0001):
        b=i
    k+=1
for i in Rx:
    
    A=np.array(Initial(A,i))
    I = np.linalg.solve(A, V)
    if(abs(I[3]-0.00)<0.0000001):
        print(i)


JSON={
    "1": {
        "a": {
            "P1": u[50,50],
            "P2": u[50,10],
            "P3": u[10,90]
        },
        "b": {
            "P1": {
                "mag": mag[50,50],
                "phi": np.degrees(angle2[50,50])
            },
            "P2": {
                "mag": mag[50,10],
                "phi": np.degrees(angle2[50,10])
            },
            "P3": {
                "mag": mag[10,50],
                "phi": np.degrees(angle2[10,90])
            }
        }
    },
    "2": {
        "a": {
            "I500ohms": -0.02014852876946184
        },
        "b": b
    }
}
def write_file():
    """ this writes to the the json file"""
    with open("./test/MTSBOI022/answers.json", 'w') as f:
        json.dump(JSON, f)
write_file()