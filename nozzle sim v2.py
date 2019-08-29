#hT1ps://www.grc.nasa.gov/WWW/K-12/airplane/rkT1hsum.html
#hT1ps://www.grc.nasa.gov/WWW/K-12/airplane/isentrop.html
#hT1ps://www.grc.nasa.gov/www/k-12/airplane/mflchk.html

import math
import sys


p1=6894757.29318 #1 is chamber pressure
k = 1.2 #property of your working fluid
T1=3633.3333 #Fahrenheit_to_Kelvin(3539.93) #temperature in the chamber in kelvin
p3=67913.3593 #free stream pressure outside nozzle -> Pa        1 atm?
p2=p3 #True for best performance
AverageMolecularWeight=0.02267 #this is in kg. Not typical for molecular weight. Rs=Ru/MolarMass

"""
print('Number of arguments:'+ len(sys.argv))
print('Argument List:'+ str(sys.argv))
if len(sys.argv) == 6:
    p1 = sys.argv[0]
    k = sys.argv[1]
    T1 = sys.argv[2]
    p3 = sys.argv[3]
    p2 = sys.argv[4]
    AverageMolecularWeight = sys.argv[5]
elif len(sys.argv) == 5:
    p1 = sys.argv[0]
    k = sys.argv[1]
    T1 = sys.argv[2]
    p3 = sys.argv[3]
    p2 = p3
    AverageMolecularWeight = sys.argv[4]
elif len(sys.argv) == 1:
    if(lower(sys.argv[0]) == "help"):
        print("p1, k, T1, p3, p2(optional), MolecWeight")
"""
def Fahrenheit_to_Kelvin(F):
    return 273.5 + ((F - 32.0) * (5.0/9.0))

pi=3.14159265
R=8.314462 #gas constant
Rspecific=R/AverageMolecularWeight

pt=((2/(k+1))**(k/(k-1)))*p1 #57
Tt=(2*T1)/(k+1) #page 57
vt=math.sqrt(((2*k)/(k+1))*Rspecific*T1) #page 57


v2=math.sqrt(((2*k)/(k-1))*Rspecific*T1*(1-(p2/p1)**((k-1)/k))) #page 54
T2=T1*(p2/p1)**((k-1)/k)
M2=v2/math.sqrt(k*Rspecific*T2)

V1=Rspecific*T1/p1 #page 55
Vt=V1*((k+1)/2)**(1/(k-1)) #page 57
while True:
    print("""Select input parameter by number
    1. Mass Flow
    2. Throat radius
    3. Thrust""")
    mode=input()
    if mode=="1":
        mdot=float(input("Input Desired Mass Flow in kg/s"))
        At=(mdot*Vt)/vt
        #print("At 1: " + (mdot/pt)*math.sqrt((Rspecific*Tt)/AverageMolecularWeight*k)) #wrong
        #print("At 2: " + str((mdot*Vt)/vt)) #derived from 3-24 on page 59 #correct
        #print("At 3: " + str((mdot*math.sqrt(T1)*((1/2)*(k-1)*(M2**2)+1)**((-1*k)/(2-2*k)-1/(2-2*k)))/(M2*p1*math.sqrt(k/Rspecific)))) #wrong
        #A2=At*(1/M2)*((1+((k-1)/2))/((k+1)/2))**((k+1)/(2*(k-1)))
        A2overAt=((1+M2**2*(k-1)/2)**((k+1)/(k-1)/2))*(((k+1)/2)**(-1*((k+1)/(k-1)/2)))/M2
        A2=A2overAt*At
        break
    elif mode=="2":
        r=float(input("Input Desired throat radius in meters"))
        At=pi*(r**2) #area of a circle from radius
        A2=At/((((k+1)/2)**(1/(k-1)))*((p2/p1)**(1/k))*math.sqrt(((k+1)/(k-1))*(1-((p2/p1)**((k-1)/k)))))
        mdot=At*p1*k*((2/(k+1))**((k+1)/(k-1)))/math.sqrt(k*Rspecific*T1)
        break
    else:
        print("Input not recognised. Please try again.")


F = mdot * v2 + (p2 - p3) * A2

print("Force: "+str(round(F, 4))+"\n")

print("Mass Flow Rate: "+str(round(mdot, 4)))
print("Mach at exit: "+str(round(M2, 4))+"\n")

#print("p1: "+str(round(p1, 4))) #given
#print("T1: "+str(round(T1, 4))+"\n") #given

print("vt: "+str(round(vt, 4))) #HH p10
print("pt: "+str(round(pt, 4))) #HH p9 close nuf
print("At: "+str(round(At, 4))) #HH p9
print("Tt: "+str(round(Tt, 4))+"\n") #HH pg 9

print("v2: "+str(round(v2, 4))) #HH pg 10 close nuf
#print("p2: "+str(round(p2, 4))) #given
print("A2: "+str(round(A2, 4))) #HH p10
print("T2: "+str(round(T2, 4))+"\n") #HH p9

#print("p3: "+str(round(p3, 4))) #given

print("p2/p1: "+str(round(p2/p1, 4))) #HH
print("A2/At: "+str(round(A2/At, 4))) ##HH pg8
print("Tt/T2: "+str(round(T2/T1, 4))) #investigate further

print("Specific Gas constant: " + str(Rspecific))

