# -*- coding: utf-8 -*-
"""
Created on Wed May 20 19:27:54 2020

@author: 
"""

from __future__ import division #no idea, thsoki stuff
import numpy
import numpy as np
import matplotlib.pyplot as plt
import timeit #possible use to time code at end, dont see the point, maybe few times for writeup
from mpl_toolkits.mplot3d import Axes3D

h = 40e-3 #height
diam = 40e-3 # diameter
d = 25e-3 # distance from source
thicc = 2e-3 # casing/housing thickness
al_d = 0.3e-3 # aluminium thickness of entrance


#geometry in cylindrical
z1 = d
z2 = z1 + al_d
z3 = z2 + h
p1 = diam/2
p2 = p1 + thicc


def cyl_to_sph_theta(p, z):
    phi = np.arctan(p/z)
    return phi

def cyl_to_sph_radius(p,z):
    r = np.sqrt(p**2 + z**2)
    return r

def car_to_sph_radius(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    return r

def car_to_sph_theta(x, y, z):
    theta = np.arctan(np.sqrt(x**2 +y**2)/z)
    return theta

def car_to_sph_phi(x, y):
    phi = np.arctan(y/x)
    return phi

def cart_to_cyl_rho(x, y):
    rho = np.sqrt(x**2 + y**2)
    return rho


# boundary conditions of entrance;

r1 = cyl_to_sph_radius(p1, z1)



#def ProbOfDec(dt):
#    1#Calculating the decay probabilities in the time interval #related to half life. not that important right now

#    t_half=5
#    return 1-numpy.exp(-dt / t_half_rad * numpy.log(2))


def sim_direct(N0,p,steps):
    """This method counts the number of nucleurs that decay per unite 
       time given the probability of decay 
       Input
            N0   :Number of nucleurs
            p    :probability of decay
            steps:Number of time steps, seconds 
            
       Returns
              Number of decayed nuclies per unite time
       """ 
    #t1=1
    #dt = t1 / steps #Calculating the interval between each time division

    N = np.zeros(steps,dtype=int)
    N[0]=N0 # initial number of nuclei
    #     p = 0.1 # decay probability in one step

    for i in range(len(N)-1):
        # generate n random numbers
        r = np.random.uniform(low=0.0, high=1.0, size=(N[i],))

        # count the number of surviving nuclei (r>p) and store in the variable
        N[i+1] =np.count_nonzero(r>p)
    count_1=0      #number of 1.17 Mev photons
    count_2=0      #Number of 1.33 Mev photons
    Nd=N0-N[len(N)-1] #total number of co-60 that has decayed
    global E
    E = []
    print(Nd)
    for i in range(Nd):
        b = np.random.random()
        if b < 0.9988: #probility of this decay path
            count_1 += 1
            count_2 += 1
            E.append(1.17) #every 1.77 MeV photon emitted means a 1.33 photon gets emitted right after
            E.append(1.33)
        else:
            count_2 += 1 #other decay path, only 2 decay paths 
            E.append(1.33)
        
            
        
        
        # #Generate photons from decayed Co 60 atoms

        # if (random.random()<0.9988):     
            
        #     if(random.random()<0.9985):
        #         #generate 1.17 Mev decay
        #         count_1+=1
        #         if(random.random()<0.9998):
        #             #Generate 1.33 Mev Photon
        #             count_2+=1
            
        # else:
        #     if(random.random()<0.9998):
        #         #Generate 1.33 Photon
        #             count_2+=1
          
    print("Number of 1.17 Mev produced:",count_1)
    print("Number of 1.33Mev photons produced:",count_2)
    print(count_1/count_2)
    print('Number of photons released',len(E),'from', N0, 'nuclei')
    return N
plt.figure(0)
plt.plot(sim_direct(10000,0.003,10000)) #easy to calculate decay probability of Co-60, related to half life, can ramp up the number of particles, low for testing
plt.xlabel("t (A.U.)")
plt.ylabel("N")
#we now have our photons in E with their energies

m_e = 0.511 #mass of electron in mev/c**2
c = 3e8 #speed of light in m/s
r0 = 7.9406e-30 #reduced compton wavelength of electrons, google to check value 
#dont bother with pair production only really relevant when E>1.5MeV, clearly seen in spectrum

def k0(E0): #look in notes at end
    return E0/(m_e*c**2)

def k(E): #same as above, look in notes
    return E/(m_e*c**2)

def sig1(k,k0): #cross section update 
    a = (2*np.pi*r0)*( ((1+k)/(k**2))*(  ((2+2*k)/(1+2*k)) - np.log(1+2*k)/k ) +np.log(1+2*k)/(2*k) - ((1+3*k)/((1+2*k)**2))       )
    return a

def scatter(E0,theta): #scatter formula, 
    a = (E0)/(1 + (E0/(m_e))*(1 - np.cos(theta)))
    return a

def compton(E0, theta):
    a = E0 - (E0)/(1+(E0/(m_e))*(1-np.cos(theta)))
    return a

def transport_and_plot():
    #Mean free paths
    lambda_a = 45 #mfp absorbtion
    lambda_s = 0.3 # mfp scattering

    #cross setions
    sigma_a = 1/lambda_a #absroption cross ssection
    sigma_s = 1/lambda_s #scattering cross section
    #Total cross section
    sigma_t = sigma_a + sigma_s
    
    #Total mean free path
    lambda_t = 1/sigma_t
    
    
    hl = np.zeros(len(E)) # particle history, number of scatters
    x = np.zeros(len(E)) #
    y = np.zeros(len(E)) ##starting position (0,0,0) of photons
    z = np.zeros(len(E)) #
    for j in range(len(E)):
        is_absorbed = 0
        A = 1
        while is_absorbed == 0:
            #s = -lambda_t*np.log(np.random.uniform(low = 0, high =1))
            
            theta = np.arcsin(-1+2*np.random.uniform(0,1)) #-pi/2 to pi/2  ##
            phi = 2*np.pi*np.random.uniform(0,1) #0-2pi                    ##Inverse CDF as well
            
            
            if A == 1:
                E1 = scatter(E[j], theta) #first scatter energy
                A = 0
            else:
                E1 = scatter(E1, theta) #subsequent scatter energies
            #print(E1) #might be broken, needs to be fixed for scattering in geometry
            
            s = -lambda_t*np.log(np.random.uniform(low = 0, high =1)) #distance travelled by photons every time step, relies on inverse CDF
            
            dz = np.sin(theta)                #
            dx = s*np.cos(theta)*np.cos(phi)  ##how position will change
            dy = s*np.cos(theta)*np.sin(phi)  #
        
            x[j] = x[j] + dx #
            y[j] = y[j] + dy ##position udate
            z[j] = z[j] + dz #
            hl[j] = hl[j] + 1 # not importnat, for testing, can be removed if you wish
        
            if np.random.uniform(0, 1) < sigma_a/sigma_t: # absorption criterion
                is_absorbed = 1
    
    print(np.mean(hl))  # should be around 150
    plt.figure(1)  # lol funny this might not work, but doesn't matter at all,still shows up
    ax = plt.axes(projection='3d')
    ax.scatter3D(x, y, z)



#transport_and_plot()




#Just if statements, most likely rejection method to determine if inside geometry and simple detection?
def transport_in_detector():
    # Mean free paths
    lambda_a = 45  # mfp absorbtion
    lambda_s = 0.3  # mfp scattering

    # cross setions
    sigma_a = 1 / lambda_a  # absroption cross ssection
    sigma_s = 1 / lambda_s  # scattering cross section
    # Total cross section
    sigma_t = sigma_a + sigma_s

    print(sigma_a/sigma_t)
    print("Absorption Criterion")

    # Total mean free path
    lambda_t = 1 / sigma_t

    hl = np.zeros(len(E))  # particle history, number of scatters
    x = np.zeros(len(E))  #
    y = np.zeros(len(E))  ##starting position (0,0,0) of photons
    z = np.zeros(len(E))
    rho = np.zeros(len(E))

    New_Energies = []
    x_paths = []
    y_paths = []
    z_paths = []
    comptoncount = 0
    photocount = 0
    for j in range(len(E)):
        # intialising angles of trajectory from point source
        #theta = 1 * np.pi * np.random.uniform(0, 1)
        theta = np.arcsin(-1 + 2 * np.random.uniform(0, 1))  # -pi/2 to pi/2  ##
        phi = 2 * np.pi * np.random.uniform(0, 1)  # 0-2pi
        s = -lambda_t*np.log(np.random.uniform(low=0, high=1))*1e-3
        inside_engaged = 0
        x_path = []
        y_path = []
        z_path = []
        for t in range(100):  # 3000  time steps
            ds = 1e-3
            dt = ds/c
            rho[j] = cart_to_cyl_rho(x[j], y[j])
            relaxation_time = s/c
            electron_interaction = relaxation_time/dt
            Absorbed  = False
            if rho[j] < p1 and z[j] > z2 and z[j] < z3 and not Absorbed:  # inside the detector
                if inside_engaged >= electron_interaction: #interaction occurs/ elecron vs
                    if np.random.uniform(0, 1) < 1/5 and (E[j] == 1.33 or E[j] == 1.17):  # absorption criterion
                        photocount = photocount + 1
                        E[j] = np.random.normal(E[j], E[j]*0.03, size=None)
                        New_Energies.append(E[j])
                        dz = 0  #
                        dx = 0  ##how position will change
                        dy = 0  #
                        x[j] = x[j] + dx  #
                        y[j] = y[j] + dy  ##position uPdate
                        z[j] = z[j] + dz  #
                        hl[j] = hl[j] + 1  # not importnat, for testing, can be removed if you wish
                        Absorbed = True
                    else: #COMPTON Scattering
                        comptoncount = comptoncount + 1
                        #theta = np.arcsin(-1 + 2 * np.random.uniform(0, 1))
                        theta = np.pi*np.random.normal(np.pi/2, np.pi/2, size=None)   #trying a uniform theta
                        #theta = 0    #trying a constant theta
                        E_compton = compton(E[j], theta)
                        E[j] = scatter(E[j], theta)
                        E[j] = np.random.normal(E[j], E[j] * 0.03, size=None)
                        E_compton = np.random.normal(E_compton, E_compton*0.03, size=None)
                        New_Energies.append(E_compton)
                        inside_engaged = 0

                        dz = ds*np.sin(theta)  #
                        dx = ds * np.cos(theta) * np.cos(phi)  ##how position will change
                        dy = ds * np.cos(theta) * np.sin(phi)  #
                        x[j] = x[j] + dx  #
                        y[j] = y[j] + dy  ##position uPdate
                        z[j] = z[j] + dz  #
                        hl[j] = hl[j] + 1  # not importnat, for testing, can be removed if you wish

                else: # NO INTERACRION YET but inside so must count some steps
                    dz = ds*np.sin(theta)  #
                    dx = ds * np.cos(theta) * np.cos(phi)  ##how position will change
                    dy = ds * np.cos(theta) * np.sin(phi)  #
                    x[j] = x[j] + dx  #
                    y[j] = y[j] + dy  ##position uPdate
                    z[j] = z[j] + dz  #
                    hl[j] = hl[j] + 1  # not importnat, for testing, can be removed if you wish
                    inside_engaged = inside_engaged + 1

            elif rho[j] < p1 and z[j] > z2 and z[j] < z3 and Absorbed: #PHOTOELCTRIC HAS OCCURRED
                dz = 0  #
                dx = 0  ##how position will change
                dy = 0  #
                x[j] = x[j] + dx  #
                y[j] = y[j] + dy  ##position uPdate
                z[j] = z[j] + dz  #
                hl[j] = hl[j] + 1  # not importnat, for testing, can be removed if you wish

            else: # outside the detector
                dz = ds*np.sin(theta)  #
                dx = ds * np.cos(theta) * np.cos(phi)  ##how position will change
                dy = ds * np.cos(theta) * np.sin(phi)  #
                x[j] = x[j] + dx  #
                y[j] = y[j] + dy  ##position uPdate
                z[j] = z[j] + dz  #
                hl[j] = hl[j] + 1  # not importnat, for testing, can be removed if you wish

            x_path.append(x[j])
            y_path.append(y[j])
            z_path.append(z[j])
        x_paths.append(x_path)
        y_paths.append(y_path)
        z_paths.append(z_path)


        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Cylinder
        x_cyl = np.linspace(-p1, p1, 100)
        X_case = np.linspace(-p2, p2, 100)
        z_cyl = np.linspace(z2, z3, 100)
        Xc, Zc = np.meshgrid(x_cyl, z_cyl)
        Xc_case, Zc_case = np.meshgrid(X_case, z_cyl)
        Yc = np.sqrt(p1**2 - Xc**2)
        Yc_case = np.sqrt(p2**2 - Xc_case**2)

        # Draw parameters
        rstride = 20
        cstride = 10
        ax.plot_surface(Xc, Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xc, -Yc, Zc, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xc_case, Yc_case, Zc_case, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xc_case, -Yc_case, Zc_case, alpha=0.2, rstride=rstride, cstride=cstride)
        for k in range(len(x_paths)):
            ax.scatter3D(x_paths[k], y_paths[k], z_paths[k], s=1)
        plt.show()


    print("Photoelectric Count")
    print(photocount)
    print("Compton Count")
    print(comptoncount)
    print("New Energies")
    print(New_Energies)
    bins = np.linspace(0, 1.50, 100)
    digitized = np.digitize(New_Energies, bins)
    histogram, bins = np.histogram(digitized, bins)
    plt.figure()
    plt.hist(New_Energies, bins, color="skyblue")
    plt.xlabel("Energy (MeV)")
    plt.ylabel("Counts")
    plt.show()

transport_in_detector()





