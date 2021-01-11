import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import os
from mpi4py import MPI
from lammps import lammps
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D 
me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

abs1=pd.read_csv('silicon_1photon_extinction.txt', sep='\t') #cm-1
abs2=pd.read_csv('silicon_2photon_absorbance.txt', sep=' ') #cm/GW
heatcapacity=0.7 #(J/(g*K))
density=2.3 #(g/cm3)
def get1phabsorbance(wavelength, data):
    return np.interp(wavelength, data['wavelength(nm)'], data['a(/cm)'])
def get2phabsorbance(wavelength, data):
    # in cm/W
    return np.interp(wavelength, data['wavelength(nm)'], data['b(cm/GW)'])/1000000000
def calculate_tempdistr(wavelength, fluence, pulselength, porosity, zmax):
    T=[]
    I0=fluence/(pulselength/10**15)
    a=get1phabsorbance(wavelength, abs1)*(1-porosity)
    b=get2phabsorbance(wavelength, abs2)*(1-porosity)
    #H=-I0*a/((I0*b+a)*np.exp(-ax)-I0*b)
    #C=b/a*np.log(I0/(I0+a/b))
    
    for i in range(zmax):
        #I=-a/b*np.exp(a/b*C)/(np.exp(a/b*C)-np.exp(a*i*10**-7))
        H=I0*a/((I0*b+a)*np.exp(a*i*10**-7)-I0*b)
        Iabs=a*H+b*H*H
        T.append(Iabs*(pulselength/10**15)*density/heatcapacity)
        #T.append(Iabs*(pulselength/10**15))
        #T.append(H)
    return T
def calculate_tempdistrstring(wavelength, fluence, pulselength, porosity, zmax):
    I0=fluence/(pulselength/10**15)
    a=get1phabsorbance(wavelength, abs1)*(1-porosity)
    b=get2phabsorbance(wavelength, abs2)*(1-porosity)
    I=str(I0*a)+'/('+str(I0*b+a)+'*exp({})'.format(str(a*10**(-7))+'*(5.43*{Size}-z))').format(Size=zmax)+'-{})'.format(I0*b)
    #C=b/a*np.log(I0/(I0+a/b))
    #I='('+str(-a/b*np.exp(a/b*C))+'/'+'('+str(np.exp(a/b*C))+'-'+ 'exp({})'.format(str(a*10**(-7))+'*(5.43*{Size}-z)').format(Size=zmax)+')'+')'
    Iabs='('+str(a)+'*'+I+'+'+str(b)+'*'+I+'*'+I+')'+'*'+str(pulselength/10**15)
    return Iabs
def run(A):
    return 'run ' + str(A)

def createinputelectrontemp(filename ,nx,ny,nz,T0):
    x_ = np.linspace(0., nx-1., nx)
    y_ = np.linspace(0., ny-1., ny)
    z_ = np.linspace(0., nz-1., nz)
    x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')
    x=x.astype(int)
    y=y.astype(int)
    z=z.astype(int)
    hdl2=open('inputstemp/'+filename,'w+')
    for obj in zip(x.flatten(),y.flatten(),z.flatten()):
        T=float(T0*np.exp(-obj[2]))
        hdl2.write(str(obj[0])+' '+str(obj[1])+' '+str(obj[2])+' '+str(T)+'\n')
    hdl2.close()
    return 0
def calculate_heatstring(wavelength, fluence, pulselength, porosity, zmax):
    I0=fluence/(pulselength/10**15)
    a=get1phabsorbance(wavelength, abs1)*(1-porosity)
    b=get2phabsorbance(wavelength, abs2)*(1-porosity)
    C=b/a*np.log(I0/(I0+a/b))
    I=str(I0*a)+'/('+str(I0*b+a)+'*exp({})'.format(str(a*10**(-7))+'*(5.43*{Size}-z)').format(Size=zmax)+'-{})'.format(I0*b)
    #I='('+str(-a/b*np.exp(a/b*C))+'/'+'('+str(np.exp(a/b*C))+'-'+ 'exp({})'.format(str(a*10**(-7))+'*(5.43*{Size}-z)').format(Size=zmax)+')'+')'
    Heatsource='('+str(a)+'*'+I+'+'+str(b)+'*'+I+'*'+I+')'
    return Heatsource


C_e=0.00144
rho_e = 0.781
kappa_e = 0
gamma_p = 0.277 
gamma_s = 0.2191
v_0 = 0.0590
Nx = 10
Ny = 10
Nz = 10
N = 50

def getscript( wavelength, fluence, pulselength,porosity,zmax, size, center, radius,inputfilename=''):
    Cylindertop = 50
    Cylinderbottom = 0
    filename= inputfilename
    fluxfilename='NPs_flux_output{wavelength}nm{fluence}Jcm2{pulselength}fs{porosity}por{size}nmsize'.format(wavelength=wavelength,
    fluence=fluence, pulselength=pulselength, porosity=porosity,size=size)
    moviename ='3d{wavelength}nm{fluence}Jcm2{pulselength}fs{porosity}por{size}nmsize.mpg'.format(wavelength=wavelength,
    fluence=fluence, pulselength=pulselength, porosity=porosity,size=size)
    flatmoviename='flat{wavelength}nm{fluence}Jcm2{pulselength}fs{porosity}por{size}nmsize.mpg'.format(wavelength=wavelength,
    fluence=fluence, pulselength=pulselength, porosity=porosity,size=size)
    slicemoviename='slice{wavelength}nm{fluence}Jcm2{pulselength}fs{porosity}por{size}nmsize.mpg'.format(wavelength=wavelength,
    fluence=fluence, pulselength=pulselength, porosity=porosity,size=size)
    flatslicemoviename='flatslice{wavelength}nm{fluence}Jcm2{pulselength}fs{porosity}por{size}nmsize.mpg'.format(wavelength=wavelength,
    fluence=fluence, pulselength=pulselength, porosity=porosity,size=size)
    #INIT
    Initialize='''
    units metal
    atom_style atomic
    atom_modify map array
    dimension       3
    boundary        p p f
    '''
    #DEFINE REGIONS
    Regions=''' 
    lattice  diamond 5.43
    region 1 block 0 {Size} 0 {Size} 1 50 units lattice
    
    #region 3 cylinder z {Center} {Center} {Radius} {Cylinderbottom} {Cylindertop} units lattice 
    region 3 block 0 {Size} 0 {Size} 1 50 units lattice
    region centerslice block {slicex1} {slicex2} 0 {Size} 1 100 units lattice
    region todump block 0 {Size} 0 {Size} 52 100 units lattice
    '''.format(Size=size, Center=center, Radius=radius, 
               Cylinderbottom=Cylinderbottom, Cylindertop=Cylindertop,slicex1=size//2-2, slicex2=size//2)
    #CREATE ATOMS
    Atomcreation=''' 
        create_atoms 1 region 1
        '''
    #OR LOAD THEM FROM FILE
    Atomreading=''' 
        read_data {filename}
        '''.format(filename=filename)
    #ATOM GROUPS AND PARAMS
    Atomconfig=''' 
    group allatoms type 1
    group centeratoms1 region 3
    variable Nn equal count(all,todump)
    mass 1 28.0
    pair_style tersoff 
    pair_coeff * * Si.tersoff  Si
    neighbor 0.5 bin
    neigh_modify every 1 delay 1 check yes 
    '''
    #CALC CONFIG
    Calcconf=''' 
    thermo 10
    thermo_modify lost ignore flush yes
    velocity allatoms create 300 12345
    fix 2 all nve
    fix 3 all ave/time 1 1 1 v_Nn file {fluxfilename}
    compute coord all coord/atom cutoff 2.6
    '''.format(fluxfilename='outputs/'+fluxfilename)
    #DUMP CONFIG
    Dumpoptions=''' 
    dump  8  all  movie  10  {moviename}  c_coord  type  zoom  1.0  adiam  2.2  axes  yes  0.8  0.02  view  60  -30  size  2048  2048
    dump  9  all  movie  10  {flatmoviename}  c_coord  type  zoom  1.6  adiam  2.2  axes  yes  0.8  0.02  view  90  0  size  2048  2048
    dump  12  all  movie  10 {slicemoviename}  c_coord  type  zoom  1.0  adiam  2.2  axes  yes  0.8  0.02  view  60  -30  size  2048  2048
    dump  13  all  movie  10 {flatslicemoviename}  c_coord  type  zoom  1.6  adiam  2.2  axes  yes  0.8  0.02  view  90  0  size  2048  2048

    dump_modify 8 amap 1 6 cf 0.0 3 min blue 0.6 yellow max green
    dump_modify 9 amap 1 6 cf 0.0 3 min blue 0.6 yellow max green
    dump_modify 12 amap 1 6 cf 0.0 3 min blue 0.6 yellow max green region centerslice
    dump_modify 13 amap 1 6 cf 0.0 3 min blue 0.6 yellow max green region centerslice
    '''.format(moviename='outputs/'+moviename,flatmoviename='outputs/'+flatmoviename,
               slicemoviename='outputs/'+slicemoviename,flatslicemoviename='outputs/'+flatslicemoviename)
    #CREATE PULSEVELOCITYDIST
    #ev/ps
    createpulse=''' 
    variable heatsource atom ({heatsource})*1.1*10e-16
    fix 4 allatoms heat 1 v_heatsource region 3
    '''.format(heatsource=calculate_heatstring(wavelength, fluence, pulselength, porosity, zmax))
    #Removeheat
    Removeheat='''
    unfix 4
    '''
    #CREATE SCRIPT
    script=Initialize+Regions+Atomreading+ Atomconfig +Calcconf +Dumpoptions +run(100)+createpulse+run(pulselength)+Removeheat+run(1100-pulselength)
    return script


inputfilename='microporSi50percent_1010100.csv'
wavelength=300
fluence=1
pulselength=1000000 #fs
porosity=0.5
zmax=50
Size = 10
Center = 5.0
Radius = 2.0

for fluence in [1,2,5,10,20]:
    inputfilename='microporSi50percent_1010100.csv'
    porosity=0.5
    for wavelength in [300,400,500,600,700,800,900,1000]:
        for pulselength in [10,20,40,100,200,400,1000]:
            #hdl=open('torun/script_{}nm_{}Jcm2_{}fs_{}_por'.format(wavelength, fluence, pulselength, porosity),'w')
            lmp = lammps()
            print(wavelength,pulselength)
            script=getscript(wavelength, fluence, pulselength, porosity,zmax, Size, Center, Radius,inputfilename)
            for line in script.split('\n'):
            #hdl.write(script)
            #hdl.close()
                lmp.command(line)
    inputfilename='Si_1010100.csv'
    porosity=0
    for wavelength in [300,400,500,600,700,800,900,1000]:
        for pulselength in [10,20,40,100,200,400,1000]:
            #hdl=open('torun/script_{}nm_{}Jcm2_{}fs_{}_por'.format(wavelength, fluence, pulselength, porosity),'w')
            lmp = lammps()
            print(wavelength,pulselength)
            script=getscript(wavelength, fluence, pulselength, porosity,zmax, Size, Center, Radius,inputfilename)
            for line in script.split('\n'):
            #hdl.write(script)
            #hdl.close()
                lmp.command(line)
print("Proc %d out of %d procs has" % (me,nprocs), lmp)
