# coding=utf-8
import netCDF4 as nc
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import time as timelibrary
import os,sys,stat,csv,math,shutil
import mpi4py.MPI as MPI   
comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

ave_file = nc.Dataset("/BIGDATA1/iocas_mmu_7/support/PSLClm.nc")
SPLA = ave_file.variables['PSLClm'][0]

nao_file = nc.Dataset("/BIGDATA1/iocas_mmu_7/support/NAOPSL.nc")
NAO = nao_file.variables['NAOPSL'][:]

start_date = '0052-11-10-00000'
output_date = '0052-11-11-00000'
base_dir = '/BIGDATA1/iocas_mmu_7/cesm/1.2.2/ice/14/parallel/' + str(comm_rank) + '/mycase/'
run_dir = '/BIGDATA1/iocas_mmu_7/cesm/1.2.2/ice/14/parallel/' + str(comm_rank) + '/mycase/run/'
fin_dir = '/BIGDATA1/iocas_mmu_7/cesm/1.2.2/ice/14/parallel/' + str(comm_rank) + '/mycase/finished'
file_name = '/BIGDATA1/iocas_mmu_7/cesm/1.2.2/ice/14/parallel/result_min_1110.txt'
duration = 15
basic = 1.58523
init_value = 1e10
error_value = 10

pre_file = nc.Dataset(run_dir + 'mycase.cam.r.' + start_date + '.nc', mode='a')
U = pre_file.variables['U'][:]
PT = pre_file.variables['PT'][:]
V = pre_file.variables['V'][:]
Q = pre_file.variables['Q'][:]
PS = pre_file.variables['PS'][:]
PHIS = pre_file.variables['PHIS'][:]
pre_file.close()

cur_file = nc.Dataset(run_dir + 'mycase.cam.h0.0052-11-11-00000' + '.nc', mode='r')
lat = cur_file.variables['lat'][:]
lev = cur_file.variables['lev'][:][25]
lon = cur_file.variables['lon'][:]
cur_file.close()

class PSO():
    def __init__(self, pN, dim, max_iter):
        self.w = 0.8
        self.c1 = 2
        self.c2 = 2
        self.r1 = 0.6
        self.r2 = 0.3
        self.pN = pN
        self.dim = dim
        self.max_iter = max_iter
        self.X = np.zeros((self.pN, self.dim))
        self.V = np.zeros((self.pN, self.dim))
        self.pbest = np.zeros((self.pN, self.dim))
        self.gbest = np.zeros((1, self.dim))
        self.p_fit = np.zeros(self.pN)
        self.fit = init_value

    def write_to_file(self, val):
        with open(file_name, 'a') as f:
            f.write(str(val))
            f.write('\n')

    def Edry(self, u, v, t, ps):
        Ra = 287.04
        Tr = 270
        Pir = 1000
        Cp = 1005.7
        u_value = 0
        v_value = 0
        t_value = 0
        ps_value = 0
        len_u = len(u)
        len_ps = len(ps)
        for i in range(0, len_u):
            u_value += u[i] * u[i]
            v_value += v[i] * v[i]
            t_value += Cp / Tr * t[i] * t[i]
        for j in range(0, len_ps):           
                ps_value +=(ps[j] / Pir) * (ps[j] / Pir)

        value = u_value + v_value + Cp / Tr *t_value + Ra * Tr *ps_value
        return value * 0.1

    def funT(self, t):
        t_value=0
        for i in range(0, 32):
            for j in range(0, 288):
                t_value += t[i*288+j] * t[i*288+j]*math.cos((60.78534+i*0.94241)/180*math.pi)
        return t_value


    def norm(self, ini):
        return np.linalg.norm(ini, ord=2)

    def write_per(self, ini, value):
        level = 1
        row = 32
        column = 288
        writer = nc.Dataset(run_dir + 'ini_' + str(value) + '.nc', 'w', format='NETCDF4')
        writer.createDimension('lev', level)
        writer.createDimension('lat', row)
        writer.createDimension('lon', column)
        writer.createVariable("lev", 'd', ("lev"))
        writer.createVariable("lon", 'd', ("lon"))
        writer.createVariable("lat", 'd', ("lat"))
        writer.variables['lat'][:] = lat[160: 192]
        writer.variables['lon'][:] = lon
        writer.variables['lev'][:] = lev
        writer.createVariable("T_per", "d", ("lev", "lat", "lon"))
        expT = [[[0 for i in range(0, column)] for j in range(0, row)] for k in range(0, level)]
        for i in range(0, row):
            for j in range(0, column):
                expT[0][i][j] = ini[i * column + j]
        writer.variables['T_per'][:] = expT

    def write_output(self, value):
        out_file = nc.Dataset(run_dir + 'mycase.cam.h0.' + output_date + '.nc', mode='r')
        PSL = out_file.variables['PSL'][:]
        row = 192
        column = 288
        writer = nc.Dataset(run_dir + 'NAOI_' + str(value) + '.nc', 'w', format='NETCDF4')
        writer.createDimension('lat', row)
        writer.createDimension('lon', column)
        writer.createDimension('time', duration)
        writer.createVariable("lon", 'd', ("lon"))
        writer.createVariable("lat", 'd', ("lat"))
        writer.createVariable("time", 'd', ("time"))
        writer.variables['lat'][:] = lat
        writer.variables['lon'][:] = lon
        writer.variables['time'][:] = out_file.variables['time'][:]
        writer.createVariable("PSL", "d", ("time", "lat", "lon"))
        psl = [[[0 for i in range(0, column)] for j in range(0, row)]for k in range(0, duration)]
        writer.variables['PSL'][:] = PSL


    def function(self, x):
        level = 26
        row = 192
        column = 288
        pre_file = nc.Dataset(run_dir + 'mycase.cam.r.' + start_date + '.nc', mode='a')
        ini = np.dot(np.array(x),np.array(reduction))
        print('restore finish.')
        ini_T = ini[0 : 32 * 288]
        funt = self.funT(ini_T)
        # self.write_to_file('edry: ' + str(edry))
        if (funt > 100):
            ini = 100 / funt * ini
        ini_expT = [[[0 for i in range(0, column)] for j in range(0, row)] for k in range(0, level)]
        for k in range(0, level):
            for i in range(0, row):
                for j in range(0, column):
                    if k==25:
                        if i > 159:
                            ini_expT[k][i][j] = ini[(i-160) * column + j] + PT[k][i][j]
                        else:
                            ini_expT[k][i][j] = PT[k][i][j]
                    else:
                        ini_expT[k][i][j] = PT[k][i][j]
        pre_file.variables['PT'][:] = ini_expT[:] 
        pre_file.close()
        if comm_rank is 1:
            timelibrary.sleep(30)
        elif comm_rank is 2:
            timelibrary.sleep(60)
        os.system('yhbatch -n 48 ' + base_dir +  'mycase.run')
        check = True
        while check:
            timelibrary.sleep(30)
            if os.path.exists(fin_dir):
                check = None
                os.rmdir(fin_dir)
        print("break ok")
        if(os.path.exists(run_dir + 'rpointer.rof')):
            os.remove(run_dir + 'rpointer.rof')
        if(os.path.exists(run_dir + 'rpointer.ocn')):
            os.remove(run_dir + 'rpointer.ocn')
        if(os.path.exists(run_dir + 'rpointer.lnd')):
            os.remove(run_dir + 'rpointer.lnd')
        if(os.path.exists(run_dir + 'rpointer.ice')):
            os.remove(run_dir + 'rpointer.ice')
        if(os.path.exists(run_dir + 'rpointer.drv')):
            os.remove(run_dir + 'rpointer.drv')
        if(os.path.exists(run_dir + 'rpointer.atm')):
            os.remove(run_dir + 'rpointer.atm')

        with open(run_dir + 'rpointer.rof','w') as f:
                 f.write('./mycase.rtm.r.' + start_date + '.nc')
        with open(run_dir + 'rpointer.ocn','w') as f:
                 f.write('mycase.docn.r.' + start_date + '.nc\n')
                 f.write('mycase.docn.rs1.' + start_date + '.bin')
        with open(run_dir + 'rpointer.lnd','w') as f:
                 f.write('./mycase.clm2.r.' + start_date + '.nc')
        with open(run_dir + 'rpointer.ice','w') as f:
                 f.write('mycase.cice.r.' + start_date + '.nc')
        with open(run_dir + 'rpointer.drv','w') as f:
                 f.write('mycase.cpl.r.' + start_date + '.nc')
        with open(run_dir + 'rpointer.atm','w') as f:
                 f.write('mycase.cam.r.' + start_date + '.nc')
        if(not os.path.exists(run_dir + 'mycase.cam.h0.' + output_date + '.nc')):
            return error_value
        fun_file = nc.Dataset(run_dir + 'mycase.cam.h0.' + output_date + '.nc')
        try:
            fun_PSL2 = fun_file.variables['PSL'][duration - 1][:]
        except IndexError:
            self.write_to_file('Index Error')
            return error_value
        fun_PSL = [[0 for i in range(0, 288)] for j in range(0, 192)]
        for i in range(0, 192):
            for j in range(0, 288):
                if j < 144:
                    fun_PSL[i][j]=fun_PSL2[i][j+144]
                else:
                    fun_PSL[i][j]=fun_PSL2[i][j-144]
        fun_PSL = fun_PSL - SPLA
        value_a = 0
        value_b = 0
        SPLA_e = [[0 for i in range(0, 105)] for j in range(0, 65)]
        for i in range(0, 65):
            for j in range(0, 105):
                SPLA_e[i][j] = fun_PSL[i+117][j+72]           
        for i in range(0, 65):
            for j in range(0, 105):
                value_a += SPLA_e[i][j]*NAO[i][j]*math.sqrt(math.cos((20.26178+i*0.94241)/180*math.pi))
                value_b += NAO[i][j]*NAO[i][j]
        value_fun = (value_a/31914.543)/value_b
        value = value_fun - basic
        self.write_to_file(str(comm_rank) + ': ' + str(value))
        if value < self.fit:
            self.write_per(ini, value)
            self.write_output(value)
        return value

    def init_Population(self):
        for i in range(0, 10 * comm_rank + self.pN):
            for j in range(0, self.dim):
                self.X[i][j] = random.uniform(1,10)
                self.V[i][j] = random.uniform(1,10)
            self.pbest[i] = self.X[i]
            tmp = self.function(self.X[i])
            self.p_fit[i] = tmp
            if tmp < self.fit:
                self.fit = tmp
                self.gbest = self.X[i]
        if comm_rank is 0:
            all_fit = comm.gather(self.fit, root=0)
            all_gbest = comm.gather(self.gbest, root=0)
            self.write_to_file(all_fit)
            for k in range(0, len(all_fit)):
                if all_fit[k] < self.fit:
                    self.fit = all_fit[k]
                    self.gbest = all_gbest[k]
        self.fit = comm.bcast(self.fit if rank == 0 else None, root=0)
        self.gbest = comm.bcast(self.gbest if rank == 0 else None, root=0)
        if comm_rank is 0:
            self.write_to_file(self.fit)


    def iterator(self):
        fitness = []
        for t in range(0, self.max_iter):
            for i in range(0, self.pN):
                temp = self.function(self.X[i])
                if (temp < self.p_fit[i]):
                    self.p_fit[i] = temp
                    self.pbest[i] = self.X[i]
                    if (self.p_fit[i] < self.fit):
                        self.gbest = self.X[i]
                        self.fit = self.p_fit[i]
            if comm_rank is 0:
                all_fit = comm.gather(self.fit, root=0)
                all_gbest = comm.gather(self.gbest, root=0)
                self.write_to_file(all_fit)
                for k in range(0, len(all_fit)):
                    if all_fit[k] < self.fit:
                        self.fit = all_fit[k]
                        self.gbest = all_gbest[k]
            self.fit = comm.bcast(self.fit if rank == 0 else None, root=0)
            self.gbest = comm.bcast(self.gbest if rank == 0 else None, root=0)
            for i in range(0, self.pN):
                self.V[i] = self.w * self.V[i] + self.c1 * self.r1 * (self.pbest[i] - self.X[i]) + self.c2 * self.r2 * (
                self.gbest - self.X[i])
                self.X[i] = self.X[i] + self.V[i]
            fitness.append(self.fit)
            if comm_rank is 0:
                self.write_to_file('\n')
                self.write_to_file(self.fit)
                self.write_to_file('\n')
        return fitness


reduction = []
reader = csv.reader(open('/BIGDATA1/iocas_mmu_7/cesm/1.2.2/ice/sample2/generate/reductions_T.csv', "r"))
for line in reader:
    lineArr = [float(n.strip()) for n in line]
    reduction.append(lineArr)
my_pso = PSO(pN=10, dim=50, max_iter=100)
my_pso.init_Population()
fitness = my_pso.iterator()
print(my_pso.fit)