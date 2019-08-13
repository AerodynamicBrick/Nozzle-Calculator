#https://www.grc.nasa.gov/WWW/K-12/airplane/rktthsum.html
#https://www.grc.nasa.gov/WWW/K-12/airplane/isentrop.html
#https://www.grc.nasa.gov/www/k-12/airplane/mflchk.html

import math


def Fahrenheit_to_Kelvin(F): 
    return 273.5 + ((F - 32.0) * (5.0/9.0))


Pt=2068000 #1500psi --> Pa #Pt is chamber pressure
gam = 1.3 #property of your working fluid
Tt=Fahrenheit_to_Kelvin(3539.93) #temperature in the chamber
p0=101352.932 #free stream pressure outside nozzle -> Pa        1 atm?
AverageMolecularWeight=0.024049 #this is in kg. Not typical for molecular weight

R=8.314462 #gas constant

while True:
    print("""Select input parameter by number
    1. Mass Flow
    2. Throat radius""")
    mode=input()
    if mode=="1":
        userMassFlow=float(input("Input Desired Mass Flow in kg/s"))
        r = (2**((-1*(gam+1))/(4*(gam-1)))*(gam+1)**(-(gam+1)/(2*(2-(2*gam))))*math.sqrt(userMassFlow)*Tt**(1/4))/(math.sqrt(3.14159)*math.sqrt(Pt)*(gam/R)**(1/4))
        break
    elif mode=="2":
        r=float(input("Input Desired throat radius in meters"))
        break
    else:
        print("Input not recognised. Please try again.")

Rspecific=R/AverageMolecularWeight

#machExit=math.sqrt((2/(gam-1))*(Pt/p0)**((gam-1)/gam)-1) #WROOOOOOONNNG?
Texit=((p0/Pt)**((gam - 1)/gam))*Tt
Vexit=math.sqrt(((2*gam)/(gam-1))*Rspecific*Tt*(1-(p0/Pt)**((gam-1)/gam)))
machExit=Vexit/math.sqrt(gam*R*Texit) #derived. more below

#AoverAstar=9.297 #get from calculator
#to get area at exit for use in base equations
Athroat=3.1415*(r**2.0)
AoverAstar=(((gam+1)/2)**(-1*((gam+1)/(2*(gam-1)))))*(((1+((gam-1)/2)*((machExit)**2))**((gam+1)/(2*(gam-1))))/machExit)
Aexit=Athroat*AoverAstar

#base equasions
mdot=((Athroat*Pt)/math.sqrt(Tt))*math.sqrt(gam/R)*((gam+1)/2)**(-1*((gam+1)/(2*(gam-1))))


#ToverTt=((1 + (machExit**2) * ((gam-1)/2))**(-1)) updated below with RPE eq above
#Texit=ToverTt*Tt #calculated earlier from Tx/Ty=(px/py)^((k-1)/k) on page 48 of RPE
ToverTt=Texit/Tt


PeOverPt=((1+(machExit**2)*((gam-1)/2))**(-1*(gam/(gam-1))))
Pexit=PeOverPt*Pt

#VexitFromMach = machExit * math.sqrt(gam * R * Texit) #Already solved for without using mach@exit. See RPE page 52
F = mdot * Vexit + (Pexit - p0) * Aexit

print("Force: "+str(round(F, 4))+"\n")
print("Mass Flow Rate: "+str(round(mdot, 4)))

print("Choke Flow diameter: "+str(round(2*r, 4)))
print("Exit diameter: "+str(round(2*r*AoverAstar, 4)))
print("A/A*: "+str(round(AoverAstar, 4))+"\n")

print("Mach at exit: "+str(round(machExit, 4)))
print("Temperature at exit: "+str(round(Texit, 4)))
print("Pressure at exit: "+str(round(Pexit, 4)))
print("Velocity at exit: "+str(round(Vexit, 4))+"\n")
#print("Velocity at exit from Mach: "+str(round(VexitFromMach, 4))+"\n") #equal to above line

print("Specific Gas constant: " + str(Rspecific))
print("T/Tt: " + str(ToverTt))
print("Pe/Pt: " + str(PeOverPt))
