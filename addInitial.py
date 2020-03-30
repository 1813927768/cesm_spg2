import netCDF4 as nc
import numpy as np
import csv

numprocs = 5

base_dir = ['/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/' + str(i) + '/mycase' for i in range(numprocs)]
run_dir = [route + '/run/' for route in base_dir]
start_date = '0052-11-01-00000'

reduction = []
reader = csv.reader(open('/BIGDATA1/iocas_mmu_3/support/reductions_T_w.csv', "r"))
for line in reader:
    lineArr = [float(n.strip()) for n in line]
    reduction.append(lineArr)

pre_file = nc.Dataset(run_dir[0] + 'mycase.cam.r.' + start_date + '.nc', mode='a')
# U = pre_file.variables['U'][:]
PT = pre_file.variables['PT'][:]
# V = pre_file.variables['V'][:]
# Q = pre_file.variables['Q'][:]
# PS = pre_file.variables['PS'][:]
# PHIS = pre_file.variables['PHIS'][:]
pre_file.close()

def addInitial(index,x):
    level = 26
    row = 192
    column = 288
    pre_file = nc.Dataset(run_dir[index] + 'mycase.cam.r.' + start_date + '.nc', mode='a')
    ini = np.dot(np.array(x),np.array(reduction))
    ini_T = ini[0 : 32 * 288]
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