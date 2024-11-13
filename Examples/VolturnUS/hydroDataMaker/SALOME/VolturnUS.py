def writeHydroStaticProperties(Solid_Nemoh,WaterPlaneAreaSolid,Dt,OUT):
        
    props=geompy.BasicProperties(Solid_Nemoh)
    SubmergedVolume = props[2]
    COBpoint = geompy.MakeCDG(Solid_Nemoh) 
    COB = geompy.PointCoordinates(COBpoint)
    props=geompy.BasicProperties(WaterPlaneAreaSolid)
    Area = props[2]/Dt
    Inertia = geompy.Inertia(WaterPlaneAreaSolid)
    COG_WaterPlaneAreaSolid_point = geompy.MakeCDG(WaterPlaneAreaSolid)
    COG_WaterPlaneAreaSolid = geompy.PointCoordinates(COG_WaterPlaneAreaSolid_point) 
    Sx_WaterPlaneAreaSolid = COG_WaterPlaneAreaSolid[1]/Area
    Sy_WaterPlaneAreaSolid = COG_WaterPlaneAreaSolid[0]/Area
    Inertia = np.array(Inertia)/Dt
    
    if (OUT!=''):
        OUT.write("%s %6f\n" %("SubmergedVolume",SubmergedVolume))        
        OUT.write("%s %6f\n" %("COBx", COB[0])) 
        OUT.write("%s %6f\n" %("COBy", COB[1])) 
        OUT.write("%s %6f\n" %("COBz", COB[2]))      
        OUT.write("%s %6f\n" %("WaterPlaneArea", Area))    
        OUT.write("%s %6f\n" %("WaterPlaneInertia_Ixx", Inertia[0])) 
        OUT.write("%s %6f\n" %("WaterPlaneInertia_Iyy", Inertia[4])) 
        OUT.write("%s %6f\n" %("WaterPlaneInertia_Ixy", Inertia[1]))    
        OUT.write("%s %6f\n" %("WaterPlaneStaticMoment_Sx", Sx_WaterPlaneAreaSolid)) 
        OUT.write("%s %6f\n" %("WaterPlaneStaticMoment_Sy", Sy_WaterPlaneAreaSolid))  

    return SubmergedVolume

def writeMassInertia(Shape,Name,OUT,density,t):
    #if Shape is a shell, put thickness in t
    #if Shape is solid, put t to 0
    
    if t==0: #Solid
        props=geompy.BasicProperties(Shape)
        Volume = props[2]
        Mass=Volume*density
        COGpoint = geompy.MakeCDG(Shape) 
        COG = geompy.PointCoordinates(COGpoint) 
        Inertia = geompy.Inertia(Shape)
        Inertia = np.array(Inertia[:9])*density   
        Inertia = Inertia.reshape((3, 3))        
    else: #Surface
        faces=[]
        faces=geompy.ExtractShapes(Shape, geompy.ShapeType["FACE"], True)
        if faces==[]:
            faces=[Shape]
    
        Mass_old=0
        Mass_new=0
        COG_old=np.zeros(3)
        COG_new=np.zeros(3)
        Inertia_old=np.zeros([3,3])
        Inertia_new=np.zeros([3,3])
        for face in faces:
            props_i=geompy.BasicProperties(face)
            Surface_i = props_i[1]
            Mass_i=Surface_i*t*density  
            COGpoint_i = geompy.MakeCDG(face) 
            COG_i = np.array(geompy.PointCoordinates(COGpoint_i)) 
            Inertia_i = geompy.Inertia(face) 
            Inertia_i = np.array(Inertia_i[:9])*t*density
            Inertia_i = Inertia_i.reshape((3, 3))
            
            Mass_new=Mass_old+Mass_i
            COG_new=(Mass_old*COG_old+Mass_i*COG_i)/Mass_new
            Delta_COG_i=COG_i-COG_new
            Delta_COG_old=COG_old-COG_new
            Inertia_new=Inertia_i+Mass_i*(Delta_COG_i.T@Delta_COG_i*np.eye(3)-np.outer(Delta_COG_i,Delta_COG_i))+Inertia_old+Mass_old*(Delta_COG_old.T@Delta_COG_old*np.eye(3)-np.outer(Delta_COG_old,Delta_COG_old)) 
            
            Mass_old=Mass_new
            COG_old=COG_new
            Inertia_old=Inertia_new
            
        Mass=Mass_new
        COG=COG_new
        Inertia=Inertia_new
        
        
        
    if(OUT!=''):    
        OUT.write("%s %6f\n" %(Name+"_Mass", Mass))
        OUT.write("%s %6f\n" %(Name+"_COGx", COG[0])) 
        OUT.write("%s %6f\n" %(Name+"_COGy", COG[1])) 
        OUT.write("%s %6f\n" %(Name+"_COGz", COG[2])) 
        OUT.write("%s %6f\n" %(Name+"_Ixx", Inertia[0,0])) 
        OUT.write("%s %6f\n" %(Name+"_Ixy", Inertia[0,1])) 
        OUT.write("%s %6f\n" %(Name+"_Ixz", Inertia[0,2])) 
        OUT.write("%s %6f\n" %(Name+"_Iyx", Inertia[1,0])) 
        OUT.write("%s %6f\n" %(Name+"_Iyy", Inertia[1,1])) 
        OUT.write("%s %6f\n" %(Name+"_Iyz", Inertia[1,2])) 
        OUT.write("%s %6f\n" %(Name+"_Izx", Inertia[2,0])) 
        OUT.write("%s %6f\n" %(Name+"_Izy", Inertia[2,1])) 
        OUT.write("%s %6f\n" %(Name+"_Izz", Inertia[2,2])) 
    
    return Mass, COG

def computecomtot(M, com):
    # M is a vector of size n (bodies)
    # com is a matrix 3*n
    Mtot = np.sum(M)
    comtot = np.zeros(3)
    for gdl in range(3):
        comtot[gdl] = np.sum(M * com[:, gdl]) / Mtot
    return Mtot, comtot

def import_PlatformData(file_path):

    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    PlatformData = {}
    
    for line in lines:
        line = line.strip()
        parts = line.split('\t')
        field_name = parts[0]
        field_value = parts[1]
        try:
            field_value = float(field_value)
        except ValueError:
             
            pass
        
        PlatformData[field_name] = field_value
       
        
    return PlatformData

def get_script_directory():
    return os.path.dirname(os.path.realpath(__file__))    
    

###
### Script Start
###
import os
import numpy as np
import sys
import salome
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS  


###
### Initialization
###
GUI=False

if(GUI!=True):
    PlatformData = import_PlatformData(os.path.join(get_script_directory(), 'PlatformData.txt'))
    os.makedirs(os.path.join(get_script_directory(), 'Mesh'), exist_ok=True)
    os.makedirs(os.path.join(get_script_directory(), 'Mass_Inertia'), exist_ok=True)
    os.makedirs(os.path.join(get_script_directory(), 'CAD'), exist_ok=True)
          
    D_MC=PlatformData['D_MC']                          #10
    D_EC=PlatformData['D_EC']                          #12.5
    Draft=PlatformData['Draft']                        #20
    FB=PlatformData['FB']                              #15
    L_EXT=PlatformData['L_EXT']                        #51.75
    W_ARM=PlatformData['W_ARM']                        #12.5
    H_ARM=PlatformData['H_ARM']                        #7
    W_BR=PlatformData['W_BR']                          #4
    H_BR=PlatformData['H_BR']                          #4
    t_MC=PlatformData['t_MC']                          #3/100
    t_EC=PlatformData['t_EC']                          #2.75/100
    t_AR=PlatformData['t_AR']                          #2.75/100
    t_BR=PlatformData['t_BR']                          #2.5/100
    Maxsize_Nemoh=PlatformData['Maxsize_Nemoh']        #8
    Minsize_Nemoh=PlatformData['Minsize_Nemoh']        #3.5
    Platform_Name=PlatformData['Platform_Name']        #'VolturnUS'
    Rotation_z=PlatformData['Rotation_z']              #0
    Weight=PlatformData['Weight']                      #2.4168e+06   
    
    
    OUT = open(os.path.join(get_script_directory(), 'Mass_Inertia', Platform_Name +'.dat'), "w")

else:
    D_MC=10
    D_EC=12.5
    Draft=20
    FB=15
    L_EXT=51.75
    W_ARM=12.5
    H_ARM=7
    W_BR=4
    H_BR=4
    t_MC=3/100
    t_EC=2.75/100
    t_AR=2.75/100
    t_BR=2.5/100
    Maxsize_Nemoh=8
    Minsize_Nemoh=3.5
    Platform_Name='VolturnUS'
    Rotation_z=0
    Weight=2.4168e+06     
    
    OUT = ''




water_density=1025
steel_density=7850
solid_ballast_density=2450
Dt=0.1 #WaterPlaneAreaSolid thickness (arbitrary value)
W_ARM=W_ARM-0.0001 #Correction for topological reasons regarding solid nemoh shell and related mesh




###
### Geometry
###
geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

#Solids
CentralCol = geompy.MakeCylinderRH(D_MC/2, FB+Draft)
geompy.TranslateDXDYDZ(CentralCol, 0, 0, -Draft)

ExtCol1 = geompy.MakeCylinderRH(D_EC/2, FB+Draft)
geompy.TranslateDXDYDZ(ExtCol1,-L_EXT, 0, -Draft)
ExtCols = geompy.MultiRotate1DNbTimes(ExtCol1, OZ, 3)

Arm1_tmp1 = geompy.MakeBoxDXDYDZ(L_EXT, W_ARM, H_ARM)
geompy.TranslateDXDYDZ(Arm1_tmp1, -L_EXT, -W_ARM/2, -Draft)
Arm1_tmp2 = geompy.MakeCutList(Arm1_tmp1, [CentralCol, ExtCol1], True)
Arm2_tmp2=geompy.MakeRotation(Arm1_tmp2, OZ, 120*math.pi/180.0)
Arm3_tmp2=geompy.MakeRotation(Arm1_tmp2, OZ, 240*math.pi/180.0)
Arms=geompy.MakeFuseList([Arm1_tmp2, Arm2_tmp2, Arm3_tmp2],True, True)

Brace1_tmp1 = geompy.MakeBoxDXDYDZ(L_EXT, W_BR, H_BR)
geompy.TranslateDXDYDZ(Brace1_tmp1, -L_EXT, -W_BR/2, FB-H_BR)
Brace1_tmp2 = geompy.MakeCutList(Brace1_tmp1, [CentralCol, ExtCol1], True)
Braces = geompy.MultiRotate1DNbTimes(Brace1_tmp2, OZ, 3)

SWL = geompy.MakeCylinderRH(200, 300)

Solid = geompy.MakeFuseList([CentralCol, ExtCols, Arms, Braces],False, False)

Solid_Nemoh = geompy.MakeCutList(Solid, [SWL], False)


geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

geompy.addToStudy(CentralCol,'CentralCol' )

geompy.addToStudy( ExtCol1, 'ExtCol1' )
geompy.addToStudy( ExtCols, 'ExtCols' )

geompy.addToStudy( Arm1_tmp1, 'Arm1_tmp1')
geompy.addToStudy( Arm1_tmp2, 'Arm1_tmp2')
geompy.addToStudy( Arm2_tmp2, 'Arm2_tmp2')
geompy.addToStudy( Arm3_tmp2, 'Arm3_tmp2')
geompy.addToStudy( Arms, 'Arms' )

geompy.addToStudy( Brace1_tmp1, 'Brace1_tmp1' )
geompy.addToStudy( Brace1_tmp2, 'Brace1_tmp2' )
geompy.addToStudy( Braces, 'Braces' )

geompy.addToStudy( SWL, 'SWL' )
geompy.addToStudy( Solid, 'Solid' )

geompy.addToStudy( Solid_Nemoh, 'Solid_Nemoh' )


# Shells
[Shell_CentralCol] = geompy.ExtractShapes(CentralCol, geompy.ShapeType["SHELL"], True)

[Shell_ExtCols_tmp] = geompy.ExtractShapes(ExtCol1, geompy.ShapeType["SHELL"], True)
Shell_ExtCols = geompy.MultiRotate1DNbTimes(Shell_ExtCols_tmp, OZ, 3)

[Shell_Arms_tmp] = geompy.ExtractShapes(Arms, geompy.ShapeType["SHELL"], True)
Shell_Arms = geompy.MakeCutList(Shell_Arms_tmp, [Shell_CentralCol, Shell_ExtCols], True)

[Shell_Braces_tmp] = geompy.ExtractShapes(Brace1_tmp2, geompy.ShapeType["SHELL"], True)
Shell_Braces = geompy.MultiRotate1DNbTimes(Shell_Braces_tmp, OZ, 3)
Shell_Braces = geompy.MakeCutList(Shell_Braces, [Shell_CentralCol, Shell_ExtCols], True)

[Shell_Nemoh] = geompy.ExtractShapes(Solid_Nemoh, geompy.ShapeType["SHELL"], True)


geompy.addToStudy( Shell_CentralCol, 'Shell_CentralCol')

geompy.addToStudyInFather( ExtCol1, Shell_ExtCols_tmp, 'Shell_ExtCols_tmp' )
geompy.addToStudy(Shell_ExtCols, 'Shell_ExtCols')

geompy.addToStudyInFather( Arms, Shell_Arms_tmp, 'Shell_Arms_tmp' )
geompy.addToStudy(Shell_Arms, 'Shell_Arms')

geompy.addToStudyInFather( Braces, Shell_Braces_tmp, 'Shell_Braces_tmp' )
geompy.addToStudy(Shell_Braces, 'Shell_Braces' )

geompy.addToStudyInFather( Solid_Nemoh, Shell_Nemoh, 'Shell_Nemoh' )


# Waterplane
WaterPlaneArea_tmp = geompy.MakeSection(Solid_Nemoh, SWL, True)
WaterPlaneArea = geompy.MakeFaceWires([WaterPlaneArea_tmp], True)
WaterPlaneAreaSolid = geompy.MakePrismVecH2Ways(WaterPlaneArea, OZ, Dt/2)

geompy.addToStudy( WaterPlaneArea_tmp, 'WaterPlaneArea_tmp' )
geompy.addToStudy( WaterPlaneArea, 'WaterPlaneArea' )
geompy.addToStudy( WaterPlaneAreaSolid, 'WaterPlaneAreaSolid' )


###
### Hydrostatic properties
###
SubmergedVolume = writeHydroStaticProperties(Solid_Nemoh, WaterPlaneAreaSolid , Dt, OUT)


###
### Ballast eq and evaluation of mass and inertia properties
### 
MassCentralColumn,COGCentralColumn=writeMassInertia(Shell_CentralCol,"ShellCentralColumn",OUT,steel_density,t_MC)
MassExternalColumns,COGExternalColumns=writeMassInertia(Shell_ExtCols ,"ShellExternalColumns",OUT,steel_density,t_EC)
MassArms,COGArms=writeMassInertia(Shell_Arms ,"ShellArms",OUT,steel_density,t_AR)
MassBraces,COGBraces=writeMassInertia(Shell_Braces ,"ShellBraces",OUT,steel_density,t_BR)

WaterBallastMass,COGWaterBallast=writeMassInertia(Arms ,"WaterBallast",OUT,water_density ,0)

SteelMass=MassCentralColumn+MassExternalColumns+MassArms+MassBraces
SolidBallastMass=SubmergedVolume*water_density-SteelMass-WaterBallastMass-Weight
if(SolidBallastMass<0):
    Negative_BallastMass=1
else:
    Negative_BallastMass=0
    
H_ballast=SolidBallastMass/solid_ballast_density/np.pi/(D_EC/2)**2 /3
SolidBallast = geompy.MakeCylinderRH(D_EC/2, H_ballast)
geompy.TranslateDXDYDZ(SolidBallast, -L_EXT, 0, -Draft)
SolidBallast = geompy.MultiRotate1DNbTimes(SolidBallast, OZ, 3)
SolidBallastMass,COGSolidBallast=writeMassInertia(SolidBallast ,"SolidBallast",OUT,solid_ballast_density,0)


M=np.array([MassCentralColumn, MassExternalColumns, MassArms, MassBraces, WaterBallastMass, SolidBallastMass])
com=np.vstack([COGCentralColumn,COGExternalColumns, COGArms, COGBraces, COGWaterBallast, COGSolidBallast])
Mtot,comtot = computecomtot(M,com)


Solid_MOST=geompy.MakeTranslation(Solid, -comtot[0], -comtot[1], -comtot[2])
geompy.Rotate(Solid_MOST, OZ, Rotation_z)



print("Mtot:")
print(Mtot)
print("comtot: ")
print(comtot)
print("Negative BallastMass: ")
print(Negative_BallastMass)


geompy.addToStudy(SolidBallast,'SolidBallast')
geompy.addToStudy( Solid_MOST, 'Solid_MOST' )


###
### SMESH component
###
import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New()

Mesh_nemoh = smesh.Mesh(Shell_Nemoh,'Mesh_nemoh')
NETGEN_1D_2D = Mesh_nemoh.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( Maxsize_Nemoh )
NETGEN_2D_Parameters_1.SetMinSize( Minsize_Nemoh )
isDone = Mesh_nemoh.Compute()

###
### Export
###

if(GUI!=True):
    OUT.write("%s %6f\n" %("Negative_BallastMass",Negative_BallastMass))
    OUT.close()
    Mesh_nemoh.ExportDAT(os.path.join(get_script_directory(), 'Mesh', Platform_Name + '.dat'))
    geompy.ExportSTEP(Solid_MOST,os.path.join(get_script_directory(), 'CAD', Platform_Name + '.STEP'))
    geompy.ExportSTL(Solid_MOST,os.path.join(get_script_directory(), 'CAD', Platform_Name + '.STL'))
