import numpy as np
import matplotlib.pyplot as plt

gamma =1.4
Al=np.zeros(3)
Ar=np.zeros(3)
Al[0]=float(input("Input left density :"))
Al[1]=float(input("Input left velocity :"))
Al[2]=float(input("Input left pressure :"))
Ar[0]=float(input("Input right density :"))
Ar[1]=float(input("Input right velocity :"))
Ar[2]=float(input("Input right pressure :"))

# Formulation of grid point
xl=0
xr=1
ncell=int(input("Enter No. of cell within domain:"))
xnodes=np.linspace(xl,xr,ncell+1)

x=np.zeros(ncell)
for i in range(ncell):
    x[i]=(xnodes[i]+xnodes[i+1])/2
delx=xnodes[2]-xnodes[1]

# Initilization of grid values
W=np.zeros((3,ncell))
for i in range(3):
    for j in range(ncell):
        if(x[j]<0.5):
            W[i][j]=Al[i]
        else:
            W[i][j]=Ar[i]

# Speed of sound
al= (gamma*Al[2]/Al[0])**(0.5)
ar= (gamma*Ar[2]/Ar[0])**(0.5)

# Computing pressure value on left and right hand side of contact discontinuity
pstar = 0.000001
epsilon=1
iteration=0
Fl=0.0
Fr=0.0
Flprime=0.0
Frprime=0.0
ustar1=0.0
ustar2=0.0
phostarR=0.0
phostarL=0.0
# print("Iteration         pstar                              Fr                                  Fl")

# Using Newton Raphson technique to iteratively calculate the vaue
# of pressure left and right hand side of contact discontinuity
while(iteration<100):
    iteration=iteration+1
    if(pstar<=Al[2]):
        Fl= (2*al/(gamma -1))*((pstar/Al[2])**((gamma-1)*0.5/gamma)-1)
        Flprime = (1/(Al[0]*al))*(pstar/Al[2])**(-(gamma+1)*0.5/gamma)
    else:
        Fl= (pstar - Al[2])*((2/((gamma+1)*Al[0]))/(pstar +((gamma-1)*Al[2]/(gamma+1))))**0.5
        Flprime = (1-((pstar-Al[2])/(2*(pstar+((gamma-1)*Al[2]/(gamma+1))))))*((2/((gamma+1)*Al[0]))/(pstar +((gamma-1)*Al[2]/(gamma+1))))**0.5
    if(pstar<=Ar[2]):
        Fr= (2*ar/(gamma -1))*((pstar/Ar[2])**((gamma-1)*0.5/gamma)-1)
        Frprime = (1/(Ar[0]*ar))*(pstar/Ar[2])**(-(gamma+1)*0.5/gamma)
    else:
        Fr= (pstar - Ar[2])*((2/((gamma+1)*Ar[0]))/(pstar +((gamma-1)*Ar[2]/(gamma+1))))**0.5
        Frprime = (1-((pstar-Ar[2])/(2*(pstar+((gamma-1)*Ar[2]/(gamma+1))))))*((2/((gamma+1)*Ar[0]))/(pstar +((gamma-1)*Ar[2]/(gamma+1))))**0.5

    F= Fr +Fl +Ar[1] - Al[1]
    Fprime = Frprime + Flprime
    pstarnew = pstar - F/Fprime
    epsilon = (pstarnew -pstar)/pstar
    pstar=pstarnew
    # print(iteration,"          ",pstar,"          ",Fr,"          ",Fl)

ustar= Al[1] - Fl

if(pstar>Al[2]):

    print("Shock on left hand of contact discontinuity")
    rhostarL = Al[0]*((pstar/Al[2] + (gamma-1)/(gamma+1))/((pstar/Al[2])*(gamma-1)/(gamma+1) + 1))
    SL = Al[1] - al*(((gamma+1)*(pstar/(2*gamma*Al[2]))+(gamma-1)/(2*gamma))**0.5)

else:
    print("Rarefaction on left hand of contact discontinuity")
    rhostarL = Al[0]*(pstar/Al[2])**(1/gamma)
    astarL = al*(pstar/Al[2])**((gamma-1)/(2*gamma))
    SHL = Al[1] - al
    STL = ustar -astarL

if(pstar>Ar[2]):

    print("Shock on right hand of contact discontinuity")
    rhostarR = Ar[0]*((pstar/Ar[2] + (gamma-1)/(gamma+1))/((pstar/Ar[2])*(gamma-1)/(gamma+1) + 1))
    SR = Ar[1] + ar * (((gamma + 1) * (pstar / (2 * gamma * Ar[2])) + (gamma - 1) / (2 * gamma)) ** 0.5)

else:
    print("Rarefaction on right hand of contact discontinuity")
    rhostarR = Ar[0]*(pstar/Ar[2])**(1/gamma)
    astarR = ar*(pstar/Ar[2])**((gamma-1)/(2*gamma))
    SHR = Ar[1] + ar
    STR = ustar + astarR


WstarL = np.zeros(3)
WstarR = np.zeros(3)

WstarL[0] = rhostarL
WstarL[1] = ustar
WstarL[2] = pstar

WstarR[0] = rhostarR
WstarR[1] = ustar
WstarR[2] = pstar
time = float(input("Enter the time at which we need to evaluate :"))
for i in range(ncell):
    slope = (x[i] - x[int(ncell/2)-1])/time
    if(slope<ustar):
        if(pstar>Al[2]):
            #left shock
            if(slope<SL):
                for j in range(3):
                    W[j][i] = Al[j]
            else:
                for j in range(3):
                    W[j][i] = WstarL[j]
        else:
            #left rarefaction
            if(slope<SHL):
                for j in range(3):
                    W[j][i] = Al[j]
            else:
                if(slope>STL):
                    for j in range(3):
                        W[j][i] = WstarL[j]
                else:
                    W[0][i] = Al[0] * ((2/(gamma+1))+((gamma-1)/(al*(gamma+1)))*(Al[1] - slope))**(2/(gamma-1))
                    W[1][i] = (2/(gamma+1))*(al + (gamma -1)*0.5*Al[1] +slope )
                    W[2][i] = Al[2] * ((2/(gamma+1))+((gamma-1)/(al*(gamma+1)))*(Al[1] - slope))**(2*gamma/(gamma-1))
    else:
        if(pstar>Ar[2]):
            #right shock
            if(slope<SR):
                for j in range(3):
                    W[j][i] = WstarR[j]
            else:
                for j in range(3):
                    W[j][i] = Ar[j]
        else:
            #right rarefaction
            if(slope>SHR):
                for j in range(3):
                    W[j][i] = Ar[j]
            else:
                if(slope<STR):
                    for j in range(3):
                        W[j][i] = WstarR[j]
                else:
                    W[0][i] = Ar[0] * ((2/(gamma+1))-((gamma-1)/(ar*(gamma+1)))*(Ar[1] - slope))**(2/(gamma-1))
                    W[1][i] = (2/(gamma+1))*(-ar + (gamma -1)*0.5*Ar[1] +slope )
                    W[2][i] = Ar[2] * ((2/(gamma+1))-((gamma-1)/(ar*(gamma+1)))*(Ar[1] - slope))**(2*gamma/(gamma-1))

print("******************************")
print("**********ustar***************")
print("******************************")
print(ustar)
print("******************************")
print("**********Rhostar**************")
print("******************************")
print("Rho Star Left: ",rhostarL)
print("Rho Star Right:",rhostarR)
print("\n\n\n\n")
#print(STR,STL)
for i in range(ncell):
    print(i,W[0][i])
    # print(0.5*W[0,i]*(W[1,i]**2) + W[2,i]/(gamma-1))
# plt.plot(x,W[0],color='black',linestyle='-',label='density')
# plt.plot(x,W[1],color='green',linestyle='-',label='velocity')
# plt.plot(x,W[2],color='blue',linestyle='-',label='pressure')
# plt.plot(x,0.5*W[0,:]*(W[1,:]**2) + W[2,:]/(gamma-1),color='red',linestyle='-',label='energy')
# #plt.title("U-final")
# plt.legend()
# # plt.savefig('Analytical.jpeg')
# plt.show()
# plt.close()