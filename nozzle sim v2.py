
#huzel and huang
#Rocket Propulsion Elements
#www.grc.nasa.gov/
#http://www.aspirespace.org.uk
#Rao

import math
import sys
import openpyxl


p1=6894.76*350 #p1 is chamber pressure
k = 1.23 #property of your working fluid
T1=2900 #Fahrenheit_to_Kelvin(3539.93) #temperature in the chamber in kelvin
p3=6894.76*14.65 #free stream pressure outside nozzle -> Pa
p2=p3 #True for best performance
AverageMolecularWeight=0.025 #this is in kg. Not typical for molecular weight. Rs=Ru/MolarMass

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
def Meters_to_Inches(inch):
    return 39.370*inch

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
    1. Mass Flow""")
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
    else:
        print("Input not recognised. Please try again.")
    """
    elif mode=="2":
        r=float(input("Input Desired throat radius in meters"))
        At=pi*(r**2) #area of a circle from radius
        A2=At/((((k+1)/2)**(1/(k-1)))*((p2/p1)**(1/k))*math.sqrt(((k+1)/(k-1))*(1-((p2/p1)**((k-1)/k)))))
        mdot=At*p1*k*((2/(k+1))**((k+1)/(k-1)))/math.sqrt(k*Rspecific*T1)
        break"""



F = mdot * v2 + (p2 - p3) * A2

print("Force: "+str(round(F, 4))+"\n")

print("Mass Flow Rate: "+str(round(mdot, 4)))
print("Mach at exit: "+str(round(M2, 4))+"\n")

print("p1: "+str(round(p1, 4))) #given
print("T1: "+str(round(T1, 4))+"\n") #given

print("vt: "+str(round(vt, 6))) #HH p10
print("pt: "+str(round(pt, 6))) #HH p9 close nuf
print("At: "+str(round(At, 8))) #HH p9
print("Tt: "+str(round(Tt, 6))+"\n") #HH pg 9

print("v2: "+str(round(v2, 6))) #HH pg 10 close nuf
print("p2: "+str(round(p2, 4))) #given
print("A2: "+str(round(A2, 8))) #HH p10
print("T2: "+str(round(T2, 6))+"\n") #HH p9

print("p3: "+str(round(p3, 4))) #given

print("p2/p1: "+str(round(p2/p1, 6))) #HH
print("A2/At: "+str(round(A2/At, 8))) ##HH pg8
print("Tt/T2: "+str(round(Tt/T2, 6))) #fixed probably?

#print("Pc-injector/Pc-stagnation")

print("Specific Gas constant: " + str(round(Rspecific,6)))

genfileq=input("Generate a file? (yes/no)").lower()
if(genfileq=="yes"):
    print("look up desired Rao exit and throat angles. https://i.imgur.com/gq6j8UO.png")

    wb = openpyxl.load_workbook('partgen.xlsx') #name of the included file
    ws = wb.active #default sheet?

    ws['A1']="Parameter Name"
    ws['B1']="Value"
    ws['C1']="Units"

    ws['B8'] = str(float(input("Exit angle in deg:")))
    ws['c8'] = "deg"

    ws['B3'] = str(float(input("Throat angle in deg:")))
    ws['c3'] = "deg"

    ws['B4'] = str(str(A2)) #area exit
    ws['c4'] = "m^2"

    ws['B5'] = str(pi*(float(input("Combustion Chamber radius in INCHES:"))**2)) #area inlet
    ws['C5'] = "in^2"

    ws['B6'] = str(At)
    ws['C6'] = "m^2" #area throat

    ws['B7'] = str(float(input("Converging angle in degrees:")))
    ws['C7'] = "deg"

    #ws['B8'] = str(str(float(input("Rao Throat radius in INCHES:")))+" in")
    #ws['C8'] = "in"

    ws['B2'] = str(math.sqrt(At/3.14))
    ws['C2'] = "m"

    wb.save("partgenout.xlsx")
elif(genfileq.lower()=="no"):
    print("Goodbye!")
else:
    print("Input not recognized")

#heat transfer and chamber
print("\nHeat Transfer:")
LOption=input("Attempt L* calculation?")
if(LOption==lower("yes")):
    print("Still in progress")#stub
elif(LOption==lower("no")):
    LStar=1.5 #meters
Vc=LStar*At #chamber volume #rpe 287
print("Vc: "+str(round(Vc, 4)))
ts=Vc/(mdot*V1) #stay time #rpe 287
print("ts: "+str(round(ts, 4)))

Dc=input("Input chamber diameter: ")
if(Dc==""):
    Dc=0.05 #~2 in default chamber diameter

A1=pi*(Dc/2)^2 #circles dawg

ConvergenceAngle=input("Input Convergence half angle in deg (default 30 deg): ")
if(ConvergenceAngle==""):
    ConvergenceAngle=math.radians(30) #default convergence angle

Lc=Dc*math.tan(math.rad(ConvergenceAngle)) #triangles dawg
L1=(Vc-A1*Lc*(1+math.sqrt(At/A1)+At/A1))/A1 #modified from RPE pg 285

DeltaPQuestion=input(">INPUT< Pressure drop or >CALC< pressure drop")
if(DeltaPQuestion==lower("input")):
    DeltaP_regen=input("Input Pressure Drop due to cooling passages")
    if(DeltaP_regen==""):
        DeltaP_regen=p1*.1 #default of 10% of chamber pressure
elif(DeltaPQuestion==lower("calc")):
    LenCoolantPipe=input("Input Length of coolant passage")
    if(LenCoolantPipe==""):
        LenCoolantPipe=0.4572 #~18in to m default passage length
    eqivD=input("Input equivilent chamber diameter")
    if(eqivD==""):
        eqivD=0.0095 #~.375 in to m
    f_regen=input("Input friction loss coefficient")
    if(f_regen==""):
        f_regen=.03 #default friction loss coefficient #RPE pg 296
    v_regen = input("Input friction loss coefficient")
    if (v_regen == ""):
        v_regen = 1  # default flow velocity
    DeltaP_regen=.5*f_regen*v_regen^2*(LenCoolantPipe/eqivD)
    print("Pressure drop: " + DeltaP_regen)

ExhaustDensity=input("Input exhaust density: ")
if(ExhaustDensity==""):
    ExhaustDensity=0 #stub
Prandtl=input("Input Prandtl number: ")
if(Prandtl==""):
    Prandtl=0 #stub
GasConductivity=input("Input Gas Conductivity: ")
if(GasConductivity==""):
    GasConductivity=0 #stub
AbsoluteGasViscosity=input("Input Absolute Gas Viscosity: ")
if(AbsoluteGasViscosity==""):
    AbsoluteGasViscosity=0 #stub
GasFilmCoefficient = 0.023*(((ExhaustDensity)**.8)/(Dc**.2))*(Prandtl**.4)*GasConductivity/(AbsoluteGasViscosity**.8)

Cbar=input("Input Average Specific Heat of Coolant: ") #average specific heat
if(Cbar==""):
    Cbar=0  #stub
CoolantVelocity=input("Input Coolant Velocity: ")
if(CoolantVelocity==""):
    CoolantVelocity=0 #stub
Rho_coolant=input("Input Coolant Density: ")
if(Rho_coolant==""):
    Rho_coolant=0 #stub
AbsoluteCoolantViscosity=input("Input Absolute Coolant Viscosity: ")
if(AbsoluteCoolantViscosity==""):
    AbsoluteCoolantViscosity=0 #stub
CoolantConductivity=input("input Coolant Conductivity: ")
if(CoolantConductivity==""):
    CoolantConductivity=0 #stub
ChamberSurfaceArea=2*pi*(Dc/2)*L1
LiquidFilmCoefficient = 0.023*Cbar*(mdot/ChamberSurfaceArea)*((Dc*CoolantVelocity*Rho_coolant)/AbsoluteCoolantViscosity)**-2*(AbsoluteCoolantViscosity*Cbar/CoolantConductivity)**(-2/3) #rpe pg 318

CoolantTemperature=input("Input Coolant Temperature: ")
if(CoolantTemperature==""):
    CoolantTemperature=0 #stub
ChamberWallThickness=input("Input Chamber Wall Thickness ")
if(ChamberWallThickness==""):
    ChamberWallThickness=0.0254
ChamberWallConductivity=input("Input Chamber Wall Conductivity: ")
if(ChamberWallConductivity==""):
    ChamberWallConductivity=0 #stub
ConvectionHeatTransfer = (T1-CoolantTemperature)/(1/GasFilmCoefficient+ChamberWallThickness/ChamberWallConductivity+1/LiquidFilmCoefficient)

