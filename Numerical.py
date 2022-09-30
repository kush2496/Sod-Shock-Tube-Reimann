import numpy as np
import matplotlib.pyplot as plt
gamma=1.4


## Test Cases left and right values



unl=np.zeros(3)           #[density  velocity  pressure]
unr=np.zeros(3)

unl[0]=float(input("Input left density :"))
unl[1]=float(input("Input left velocity :"))
unl[2]=float(input("Input left pressure :"))
unr[0]=float(input("Input right density :"))
unr[1]=float(input("Input right velocity :"))
unr[2]=float(input("Input right pressure :"))
t = float(input("Input the time step :"))
Ul=np.zeros(3)                  #[rho    rhoU    rhoE]
Ur=np.zeros(3)
epsilon = 10**(-6)

Ul[0]=unl[0]
Ur[0]=unr[0]
Ul[1]=unl[0]*unl[1]
Ur[1]=unr[0]*unr[1]
Ul[2]=unl[0]*(unl[1]**2)*0.5 + unl[2]/(gamma-1)
Ur[2]=unr[0]*(unr[1]**2)*0.5 + unr[2]/(gamma-1)


# Formulation of grid point
xl=0
xr=1
ncell=int(input("Enter No. of cell within domain:"))
cfl = float(input("Enter CFL No.:"))
xnodes=np.linspace(xl,xr,ncell+1)

x=np.zeros(ncell)
for i in range(ncell):
    x[i]=(xnodes[i]+xnodes[i+1])/2
delx=xnodes[2]-xnodes[1]


U = np.zeros((3,ncell))
# Initializing the domain at zero time step
for i in range(ncell):
    if x[i]<0.5:
        for j in range(3):
            U[j][i]=Ul[j]           ##[rho    rhoU    rhoE]
    else:
        for j in range(3):
            U[j][i] = Ur[j]         ##[rho    rhoU    rhoE]
eigenK = np.zeros((3,3,ncell+1))


time = 0.0
while (time<t):
    delU1 = np.zeros(ncell)
    delU2 = np.zeros(ncell)
    delU3 = np.zeros(ncell)
    delTelda = np.zeros((3, ncell+1))
    lamb =np.zeros((3,ncell+1))
    wr=0
    wl=0
    Hr=0
    Hl=0
    aTelda=0
    uTelda=0
    HTelda=0
    x=np.zeros(ncell)
    for i in range(ncell - 1):
        wr = U[0][i + 1] ** 0.5
        wl = U[0][i] ** 0.5
        Hr = ((gamma*U[2][i + 1]) - 0.5 * (gamma-1)*((U[1][i + 1] **2 )/U[0][i + 1]))/ U[0][i + 1]
        Hl = ((gamma*U[2][i]) - 0.5 * (gamma-1)*((U[1][i] **2)/U[0][i]))/ U[0][i]

        uTelda = ((wr * U[1][i + 1] / U[0][i + 1]) + (wl * U[1][i] / U[0][i])) / (wr + wl)
        HTelda = ((wr * Hr) + (wl * Hl)) / (wr + wl)
        aTelda = ((gamma - 1) * (HTelda - 0.5 * (uTelda ** 2))) ** 0.5
        x[i]  = (HTelda - 0.5 * (uTelda ** 2))

        delU1[i+1] = U[0][i + 1] - U[0][i]
        delU2[i+1] = U[1][i + 1] - U[1][i]
        delU3[i+1] = U[2][i + 1] - U[2][i]

        delTelda[1][i+1] = ((delU1[i+1] * (HTelda - (uTelda ** 2))) + (uTelda * delU2[i+1]) - delU3[i+1]) * (gamma - 1) / (aTelda ** 2)
        delTelda[0][i+1] = ((delU1[i+1] * (uTelda + aTelda)) - delU2[i+1] - (aTelda * delTelda[1][i+1])) / (2 * aTelda)
        delTelda[2][i+1] = delU1[i+1] - delTelda[0][i+1] - delTelda[1][i+1]

        lamb[0][i+1] = uTelda - aTelda
        lamb[1][i+1] = uTelda
        lamb[2][i+1] = uTelda + aTelda

        eigenK[:, :, i+1] = [[1, 1, 1], [(uTelda - aTelda), uTelda, (uTelda + aTelda)],
                           [(HTelda - (uTelda * aTelda)), 0.5 * (uTelda ** 2), (HTelda + (uTelda * aTelda))]]

    # print(lamb)
    maxLamb = np.abs(lamb).max()
    # print(maxLamb)
    # print(eigenK[:, :, :])

    flux = np.zeros((3, ncell))
    fluxFace = np.zeros((3, ncell + 1))
    sum = np.zeros((3, ncell + 1))

    # Entropy Fix
    for i in range(ncell):
        for j in range(3):
            if (abs(lamb[j][i])<epsilon):
                if(lamb[j][i] !=0):
                    lamb[j][i] = 0.5*((lamb[j][i]/epsilon) + epsilon)

    # for i in range(0,ncell-1):
    #     #left expansion wave
    #     ustar1 = U[1][i] + delTelda[0][i+1]*lamb[0][i+1]/(U[0][i] + delTelda[0][i+1])
    #     rhostarL = U[0][i] + delTelda[0][i+1]
    #     pstar1 = (gamma -1) * (U[2][i] + delTelda[0][i+1]*eigenK[2][0][i+1] - 0.5*rhostarL*(ustar1**2))
    #     astarL = (gamma*rhostarL/pstar1)**0.5
    #     p1 = (gamma-1)*(U[2][i] + (U[1][i]**2)/(2*U[0][i]))
    #     aL = (gamma*p1/U[0][i])**0.5
    #     lamb1L = U[0][i] - aL
    #     lamb1R = ustar1 - astarL
    #     if(lamb1L<0 and lamb1R>0):
    #         lamb[0][i+1] = lamb1L*(lamb1R - lamb[0][i+1])/(lamb1R - lamb1L)
    #
    #     #right wxpansion wave
    #
    #     ustar3 = U[1][i+1] -delTelda[2][i+1]*lamb[2][i+1]/(U[0][i+1] - delTelda[2][i+1])
    #     rhostarR = U[0][i+1] - delTelda[2][i+1]
    #     pstar3 = (gamma -1)*(U[2][i+1] -delTelda[2][i+1]*eigenK[2][2][i+1] -0.5*rhostarR*(ustar3**2))
    #     astarR = (gamma*pstar3/rhostarR)**0.5
    #     p3 = (gamma-1)*(U[2][i+1] + (U[1][i+1]**2)/(2*U[0][i+1]))
    #     aR = (gamma*p3/U[0][i+1])**0.5
    #     lamb3L = ustar3 + astarR
    #     lamb3R = U[0][i+1] + aR
    #     if(lamb3L<0 and lamb3R>0):
    #         lamb[2][i+1] = lamb3R*(lamb[2][i+1] -lamb3L)/(lamb3R - lamb3L)


    for i in range(ncell):
        flux[0][i] = U[1][i]
        flux[1][i] = ((3 - gamma)*0.5*(U[1][i]**2)/U[0][i]) + (gamma-1)*U[2][i]
        flux[2][i] = (gamma * U[1][i] * U[2][i] / U[0][i]) - ((gamma - 1) *0.5* (U[1][i] ** 3) / (U[0][i] ** 2))


    # Computing the F i+(1/2) = F i + sum
    # Computing the F i-(1/2) = F i-1 + sum
    for i in range(1,ncell):
        for j in range(3):
            sum[:, i] += abs(lamb[j][i]) * delTelda[j][i] * eigenK[:, j, i]
        # 0th index flux face is the left face of 0th cell
        # 1st index flux face is the right face of 0th cell
        # F i+1/2 = 0.5*(Fl + Fr) - 0.5*summation(delTelda*lambda*eigenVect)
        fluxFace[:, i] = 0.5*(flux[:, i-1] + flux[:, i]) - 0.5*sum[:, i]
    fluxFace[:,0]=fluxFace[:,1]
    fluxFace[:,ncell] = fluxFace[:,ncell-1]

    # U new = U old + delT/delx * ( faceFlux i-1/2  - faceFlux i+1/2 )
    for i in range(ncell):
        U[:, i] = U[:, i] + (cfl / maxLamb) * (fluxFace[:, i] - fluxFace[:, i+1])

    time += cfl * delx / maxLamb
print("\n\n\nDensity :\n\n\n\n")
for i in range(ncell):
    print(i,U[0][i])

#
# plt.plot(U[0],color='purple',linestyle='',marker='o',markersize=1.2,label='density')
# plt.title("density")
# plt.show()
# plt.close()
# plt.plot(U[1,:]/U[0,:],color='green',linestyle='',marker='o',markersize=1.2,label='velocity')
# plt.title("velocity")
# plt.show()
# plt.close()
# plt.plot((gamma-1)*(U[2,:]-(0.5*(U[1,:]**2)/U[0,:])),color='blue',linestyle='',marker='o',markersize=1.2,label='pressure')
# plt.title("pressure")
# plt.show()
# plt.close()
# plt.plot(U[2],color='red',linestyle='',marker='o',markersize=1.2,label='energy')
# plt.title("energy")
# plt.show()
# plt.close()