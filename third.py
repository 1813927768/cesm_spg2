# coding=utf-8
# /BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/third.py

import netCDF4 as nc
import numpy as np
import random
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import time as timelibrary
import os,sys,stat,csv,math,shutil,json

#伪并行数
numprocs = 5


ave_file = nc.Dataset("/BIGDATA1/iocas_mmu_3/support/PSLClm.nc")
SPLA = ave_file.variables['PSLClm'][0]

nao_file = nc.Dataset("/BIGDATA1/iocas_mmu_3/support/NAOPSL.nc")
NAO = nao_file.variables['NAOPSL'][:]

start_date = '0052-11-15-00000'
output_date = '0052-11-16-00000'
base_dir = ['/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/' + str(i) + '/mycase' for i in range(numprocs)]
run_dir = [route + '/run/' for route in base_dir]
fin_dir = [route + '/finished' for route in base_dir]

file_name = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/result_min_1110.txt'
process_file = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/process.txt'
restart_file = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/restart.txt'
duration = 15
basic = -0.17638014536610233     # 不加扰动模式输出-0.8047237110735179
init_value = 1e10
error_value = -1

# 所有进程初始PT是一样的
pre_file = nc.Dataset('/BIGDATA1/iocas_mmu_3/support/input_15/mycase.cam.r.%s.nc'%(start_date), mode='a')
# U = pre_file.variables['U'][:]
PT = pre_file.variables['PT'][:]
# V = pre_file.variables['V'][:]
# Q = pre_file.variables['Q'][:]
# PS = pre_file.variables['PS'][:]
# PHIS = pre_file.variables['PHIS'][:]
pre_file.close()

cur_file = nc.Dataset(run_dir[0] + 'mycase.cam.h0.0052-11-01-00000' + '.nc', mode='r')
lat = cur_file.variables['lat'][:]
lev = cur_file.variables['lev'][:][25]
lon = cur_file.variables['lon'][:]
cur_file.close()


#SPG2全局变量设置
M = 5    #lineSearch考虑的之前的特征方程值个数
RMAX = 1e15
RMIN = 1e-15
MAXIT = 50
MAXIFCNT = 480
eps = 1e-6 #参照fortran代码

DELTA = 1e-15

class SPG2():
    def __init__(self, dim, numprocs, ini = None):

        self.dim = dim
        self.numprocs = numprocs
        
        # itern ifcnt igcnt 迭代次数 特征方程计算次数 梯度计算次数
        self.ifcnt = self.itern = self.igcnt = 0
        self.g = [0 for i in range(self.dim)]
        self.f = init_value
        self.fvalues = [-1.0e+99 for i in range(M)]
        self.fvalue = [0 for i in range(MAXIT)]
        self.fcnorm = [0 for i in range(MAXIT)]
        # 用来控制伪并行的进行（-1表示对应case空闲,i表示此时正在计算某个维度的偏微分值）
        self.ctl = [-1 for i in range(self.numprocs)]
        # 用来暂时存储delta值
        self.delta = [0 for i in range(self.numprocs)]
        if not ini:
            # 随机初始化x位置
            self.x =  [0 for i in range(self.dim)]
            for i in range(self.dim):
                self.x[i] = random.uniform(1,10)
        else:
            self.x = ini
        self.DELTA = DELTA


    # 重入
    def restart(self):
        # 读取重入文件
        with open(restart_file, 'r') as f:
            for line in f.readlines():
                parts = line.split("=")
                head = parts[0].strip()
                body = json.loads(parts[1])
                if head == "iter":
                    self.ifcnt = self.itern = self.igcnt = body
                elif head == "x":
                    self.x = body
                elif head == "g":
                    self.g = body
                elif head == "f":
                    self.f = body
                elif head =="fvalue":
                    self.fvalue = body
                elif head == "fvalues":
                    self.fvalues = body
                elif head == "fcnorm":
                    self.fcnorm = body
                elif head == "fbest":
                    self.fbest = body
                elif head == "xbest":
                    self.xbest = body
                elif head == "rambda":
                    self.rambda = body
                
    # 保存iter信息
    def write_to_rfile(self):
        with open(file_name, 'a') as f:
            f.write('\n'+"iter = "+ str(self.itern)+'\n')
            f.write("fvalue = "+ str(self.fvalue)+'\n')
            f.write("fvalues = "+ str(self.fvalues)+'\n')
            f.write("x = "+ str(self.x)+'\n')
            f.write("f = "+ str(self.f)+'\n')
            f.write("g = "+ str(self.g)+'\n')
            f.write("fcnorm = "+ str(self.fcnorm)+'\n')
            f.write("fbest = "+ str(self.fbest)+'\n')
            f.write("xbest = "+ str(self.xbest)+'\n')
            f.write("rambda = "+ str(self.rambda)+'\n')

    # 将扰动x投影到特征空间（修改）
    def project2(self,x):
        ini = np.dot(np.array(x),np.array(reduction))
        self.write_to_pfile('--proj start (%d)--'%(self.itern))
        self.write_to_pfile(x)
        ini_T = ini[0 : 32 * 288]
        funt = self.funT(ini_T)
        self.write_to_pfile('funt = ' + str(funt))
        # 约束
        if (funt > 100):
            # x = 100 / funt * x
            x = [100 / funt * i for i in x ]     
            self.write_to_pfile('need constraint')  
        self.write_to_pfile(x)
        self.write_to_pfile('--proj done (%d)--'%(self.itern))
        return x
    
    def project(self,x):
        self.write_to_pfile(x)
        return x

    # 通过lineSearch更新扰动位置信息
    def lineSearch(self, d, gtd):
        self.write_to_pfile("--line search start (%d)--"%(self.itern))
        self.write_to_pfile("   old f = "+str(self.f))
        self.write_to_pfile("   old x = "+str(self.x))
        self.write_to_pfile("   gtd = "+str(gtd))
        self.write_to_pfile("   d = "+str(d))
        
        GAMMA =  1.0e-04
        dnew = d
        x = self.x
        f = self.f
        fvalues = self.fvalues

        # 计算fmax
        fmax = fvalues[0]
        for i in range(1,M):
            fmax = max(fmax,fvalues[i])
        self.write_to_pfile("   fmax = "+str(fmax))

        iter = alpha = 1.0

        while True:
            # 计算xnew,新的扰动
            xnew = [x[i]+dnew[i] for i in range(len(x))]
            fcomp = fmax + alpha * GAMMA * gtd  #gtd对应δ，alpha对应λ

            self.write_to_pfile('   alpha = '+ str(alpha))
            self.write_to_pfile("   xnew = "+str(xnew))
            self.write_to_pfile('   fmax + alpha * GAMMA * gtd = '+ str(fcomp))

            # 计算fnew
            if (np.array(xnew) == np.zeros((50,))).all(): 
                fnew = 0          #如果xnew是原点，直接给出结果
            elif fcomp < -3:
                fnew = 0          # 因为f不可能小于-3所以，当fcomp<-3,直接跳过这一循环
            else:
                iter = iter + 1
                fnew = self.evalfg2(x=xnew,f=None,choice=1)
            
            self.write_to_pfile('   second call evalfg')
            self.write_to_pfile('   fnew = '+ str(fnew))

            # 判断是否满足跳出条件
            if fnew <= fcomp:
                d = [xnew[i]-x[i] for i in range(len(x))]
                self.write_to_pfile("d*alpha = "+str(d))
                self.write_to_pfile('final iter = '+ str(iter))
                break

            # 计算新的alpha
            if(alpha <= 0.1):
                alpha = alpha/3         # speed up the line-search process
            else:
                atemp = (-gtd*alpha**2)/(2.0*(fnew-f-alpha*gtd))
                if(atemp < 0.1 or atemp > 0.9*alpha):
                    atemp = alpha/2
                alpha = atemp         
            
            # d变化超出精度范围    
            dnew = [alpha*d[i] for i in range(len(d))]
            dnorm = np.linalg.norm(np.array(dnew),np.inf)   
            self.write_to_pfile("   dnorm = "+str(dnorm)) 
            if(dnorm <= 1e-17):
                raise Exception("d beyond accuracy")

        #更新扰动位置信息
        # self.x = xnew
        self.f = fnew
        self.fvalues[(self.itern % M)] = fnew
        self.fvalue[self.itern] = fnew
        #如果找到更优值
        if(self.f < self.fbest):
            self.fbest = fnew
            self.xbest = xnew
            self.write_per(self.x2ini(xnew),fnew)
            self.write_output(fnew)
            self.write_to_pfile("better f found")

        self.write_to_pfile("new f = "+str(self.f))
        self.write_to_pfile("--line search done (%d)--"%(self.itern))
        
        return xnew

    # choice
    # choice = 1 计算特征方程值
    # choice = 2 计算both
    # choice = 3 计算梯度
    # x 扰动
    # f 目标方程值
    # g 梯度
    def evalfg2(self,x,f,choice):   
        if choice == 1:
            fnew = self.function(x,0)
            self.ifcnt= self.ifcnt + 1
            return fnew
        elif choice == 2:
            fnew = self.function(x,0)
            self.ifcnt= self.ifcnt + 1
            self.igcnt= self.igcnt + 1
            gnew = self.gradient(x,fnew)
            return fnew,gnew
        elif choice == 3:
            fnew = f
            self.igcnt= self.igcnt + 1
            gnew = self.gradient(x,fnew)
            return gnew
        else:
            self.write_to_pfile("unknown choice")

    # 找到空闲的caseID
    def findFree(self):
        free = []
        for i in range(self.numprocs):
            if(self.ctl[i] == -1):
                free.append(i)
        return free

    def findBusy(self):
        busy = []
        for i in range(self.numprocs):
            if(self.ctl[i] != -1):
                busy.append(i)
        return busy

    # 求梯度用，返回某个维度加delta后的x值
    def createX(self,xold,dim):
        # 防止浅拷贝影响xold
        x = xold.copy()
        deltall = abs(x[dim]*self.DELTA)
        if(deltall == 0.0):
            deltall = self.DELTA
        x[dim] += deltall
        return x,deltall
    
    def createX2(self,xold,dim):
        x = xold.copy()
        deltall = abs(x[1]*self.DELTA*(dim))
        x[1] += deltall
        return x,deltall


    # 求梯度（伪并行）
    def gradient(self,xold,f):

        self.write_to_pfile("--gradient start (%d)--"%(self.itern))
        gnew = [0 for i in range(self.dim)]

        # 当前维度
        current = 0
        while(True):

            # 检查是否有计算完成的case返回结果
            for busyCase in self.findBusy():
                val = self.completeTask(busyCase)
                if not val is None:
                    dim = self.ctl[busyCase]
                    gnew[dim] = (val-f)/self.delta[busyCase]
                    self.ctl[busyCase] = -1
                    self.write_to_pfile("   Done! dim :" + str(dim) + " caseID :" +str(busyCase) + " val :" +str(val))


            # 检查是否有空闲case可供计算
            freeCases = self.findFree()
            # 结束条件，所有case空闲且所有梯度都已计算完成
            if(len(freeCases) == self.numprocs and current >= self.dim):
                self.write_to_pfile("all cases are free and all dims are done")
                break
            if(len(freeCases) > 0):
                for freeCase in freeCases:
                    # 防止溢出（去跑第51个case）
                    if current >= 50:
                        break
                    x,delta = self.createX(xold,current)
                    self.submitTask(x,freeCase)
                    self.ctl[freeCase] = current
                    self.delta[freeCase] = delta
                    self.write_to_pfile("   Start! dim :" + str(current) + " caseID :" +str(freeCase) + " delta :" +str(delta))
                    current += 1
                    
            
            timelibrary.sleep(30)

        self.write_to_pfile(gnew)
        self.write_to_pfile("--gradient done (%d)--"%(self.itern))

        # 将所有进程重置为空闲状态
        for i in range(self.numprocs):
            self.ctl[i] = -1

        return gnew

    def spg2(self,restartIndex = 0):


        if restartIndex > 0:
            self.write_to_pfile('restarting now!!')

            self.restart()           
            cgnorm = self.fcnorm[self.itern]
            cg = [0 for _ in range(self.dim)]

        # 第一次开始
        else:
            #将初始扰动映射到特征空间内
            self.x = self.project(self.x)
            self.xbest = self.x
            self.write_to_pfile("project done, start evalfg2")

            #求梯度和特征方程值
            self.f,self.g = self.evalfg2(choice=2,x=self.x,f=self.f)
            self.write_to_pfile('first call evalfg')

            self.fbest = self.fvalue[0] = self.fvalues[0] = self.f

            cg = [self.x[i] - self.g[i] for i in range(self.dim)]
            self.write_to_pfile('proj cg')
            cg = self.project(cg)

            cgnorm = 0.0  #cgnorm用于判断是否收敛
            for i in range(self.dim):
                cgnorm = max(cgnorm,abs(cg[i]-self.x[i]))
            #如果不收敛
            if (cgnorm != 0):
                self.rambda= 1/cgnorm    #设置α的初始值

            self.fcnorm[self.itern] = cgnorm
            self.write_to_pfile(str(self.itern)+" iter: cgnm="+str(cgnorm))
            
            self.write_to_rfile()

        # 主循环开始
        while(cgnorm > eps and self.itern < MAXIT-1 and self.ifcnt < MAXIFCNT):
            
            #控制信息更新
            self.itern += 1

            #准备lineSearch所需参数: d,gtd
            gtd = 0
            d = [self.x[i]-self.rambda*self.g[i] for i in range(self.dim)]
            self.write_to_pfile('proj d')
            d = self.project(d)
            for i in range(self.dim):
                d[i] = d[i]-self.x[i]
                gtd += self.g[i]*d[i]

            #通过lineSearch更新扰动x和特征方程值f
            xnew = self.lineSearch(d,gtd)

            #计算更新后的梯度
            gnew = self.evalfg2(choice=3,x=xnew,f=self.f) 
            self.write_to_pfile('fourth call evalfg')

            s = [0 for i in range(self.dim)]
            y = [0 for i in range(self.dim)]
            sts = sty = 0.0
            for i in range(self.dim):
                s[i] = xnew[i]-self.x[i]
                y[i]= gnew[i]-self.g[i]
                sts = sts + s[i]*s[i]
                sty = sty + s[i]*y[i]
                # 更新x,g
                self.x[i] = xnew[i]
                self.g[i] = gnew[i]
                cg[i] = self.x[i]-self.g[i]
            
            self.write_to_pfile('proj cg2')
            cg = self.project(cg)

            cgnorm = 0
            for i in range(self.dim):
                cgnorm = max(cgnorm,abs(cg[i]-self.x[i]))
            
            if(sty <= 0):
                self.rambda = RMAX
            else:
                self.rambda = min(RMAX,max(RMIN,sts/sty))

            self.fcnorm[self.itern] = cgnorm
            self.write_to_pfile(str(self.itern)+" iter: cgnm="+str(cgnorm))

            self.write_to_rfile()
            #主循环结束

        self.f = self.fbest
        self.x = self.xbest

        self.write_to_pfile("fbest = "+ str(self.f))
        self.write_to_pfile("xbest = "+ str(self.x))

        if(cgnorm <= eps):
            self.write_to_pfile('convergence')
        elif(self.itern > MAXIT):
            self.write_to_pfile('too many iterations')
        elif(self.ifcnt>MAXIFCNT):
            self.write_to_pfile('too many function evaluations')
        else:
            self.write_to_pfile('unknown stop')

        return
        
    def write_to_file(self, val):
        with open(file_name, 'a') as f:
            f.write(str(val))
            f.write('\n')

    # 记录过程信息
    def write_to_pfile(self,val):
        with open(process_file, 'a') as f:
            f.write(str(val))
            f.write('\n')

    # 原空间里x的约束函数（本质是二范数的非负加权，是凸函数）
    def funT(self, t):
        t_value=0
        for i in range(0, 32):
            for j in range(0, 288):
                t_value += t[i*288+j] * t[i*288+j]*math.cos((60.78534+i*0.94241)/180*math.pi)
        return t_value

    def write_per(self, ini, value):
        level = 1
        row = 32
        column = 288
        writer = nc.Dataset(run_dir[0] + 'ini_' + str(value) + '.nc', 'w', format='NETCDF4')
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
        out_file = nc.Dataset(run_dir[0] + 'mycase.cam.h0.' + output_date + '.nc', mode='r')
        PSL = out_file.variables['PSL'][:]
        row = 192
        column = 288
        writer = nc.Dataset(run_dir[0] + 'NAOI_' + str(value) + '.nc', 'w', format='NETCDF4')
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
        # psl = [[[0 for i in range(0, column)] for j in range(0, row)]for k in range(0, duration)]
        writer.variables['PSL'][:] = PSL
   
    # 还原
    def x2ini(self,x):
        ini = np.dot(np.array(x),np.array(reduction))
        # self.write_to_pfile('restore finish.')
        ini_T = ini[0 : 32 * 288]
        funt = self.funT(ini_T)
        # 约束
        if (funt > 100):
            ini = 100 / funt * ini
        return ini

    def submitTask(self, x, index):
        level = 26
        row = 192
        column = 288
        pre_file = nc.Dataset(run_dir[index] + 'mycase.cam.r.' + start_date + '.nc', mode='a')
        ini = np.dot(np.array(x),np.array(reduction))
        ini_T = ini[0 : 32 * 288]
        funt = self.funT(ini_T)
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
        # if rank is 1:
        #     timelibrary.sleep(30)
        # elif rank is 2:
        #     timelibrary.sleep(60)
        os.system('yhbatch -n 48 ' + base_dir[index] +  '/mycase.run')

    def completeTask(self,index):
        if not os.path.exists(fin_dir[index]):
            return None
        os.rmdir(fin_dir[index])
        if(os.path.exists(run_dir[index] + 'rpointer.rof')):
            os.remove(run_dir[index] + 'rpointer.rof')
        if(os.path.exists(run_dir[index] + 'rpointer.ocn')):
            os.remove(run_dir[index] + 'rpointer.ocn')
        if(os.path.exists(run_dir[index] + 'rpointer.lnd')):
            os.remove(run_dir[index] + 'rpointer.lnd')
        if(os.path.exists(run_dir[index] + 'rpointer.ice')):
            os.remove(run_dir[index] + 'rpointer.ice')
        if(os.path.exists(run_dir[index] + 'rpointer.drv')):
            os.remove(run_dir[index] + 'rpointer.drv')
        if(os.path.exists(run_dir[index] + 'rpointer.atm')):
            os.remove(run_dir[index] + 'rpointer.atm')

        with open(run_dir[index] + 'rpointer.rof','w') as f:
                 f.write('./mycase.rtm.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.ocn','w') as f:
                 f.write('mycase.docn.r.' + start_date + '.nc\n')
                 f.write('mycase.docn.rs1.' + start_date + '.bin')
        with open(run_dir[index] + 'rpointer.lnd','w') as f:
                 f.write('./mycase.clm2.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.ice','w') as f:
                 f.write('mycase.cice.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.drv','w') as f:
                 f.write('mycase.cpl.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.atm','w') as f:
                 f.write('mycase.cam.r.' + start_date + '.nc')
        if(not os.path.exists(run_dir[index] + 'mycase.cam.h0.' + output_date + '.nc')):
            self.write_to_pfile("CompleteTask: output not found in case$" + str(index))
            raise Exception("output not found")
        fun_file = nc.Dataset(run_dir[index] + 'mycase.cam.h0.' + output_date + '.nc')
        try:
            fun_PSL2 = fun_file.variables['PSL'][duration - 1][:]
        except IndexError:
            self.write_to_pfile("CompleteTask: ErrorValue in case$" + str(index))
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
        value = basic -value_fun
        return value

    # 求目标函数值
    def function(self, x, index):
        level = 26
        row = 192
        column = 288
        pre_file = nc.Dataset(run_dir[index] + 'mycase.cam.r.' + start_date + '.nc', mode='a')
        ini = np.dot(np.array(x),np.array(reduction))
        self.write_to_pfile('--func start (%d)--'%(self.itern))
        ini_T = ini[0 : 32 * 288]
        funt = self.funT(ini_T)
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
        # if rank is 1:
        #     timelibrary.sleep(30)
        # elif rank is 2:
        #     timelibrary.sleep(60)
        os.system('yhbatch -n 48 ' + base_dir[index] +  '/mycase.run')
        check = True
        while check:
            timelibrary.sleep(30)
            if os.path.exists(fin_dir[index]):
                check = None
                os.rmdir(fin_dir[index])

        if(os.path.exists(run_dir[index] + 'rpointer.rof')):
            os.remove(run_dir[index] + 'rpointer.rof')
        if(os.path.exists(run_dir[index] + 'rpointer.ocn')):
            os.remove(run_dir[index] + 'rpointer.ocn')
        if(os.path.exists(run_dir[index] + 'rpointer.lnd')):
            os.remove(run_dir[index] + 'rpointer.lnd')
        if(os.path.exists(run_dir[index] + 'rpointer.ice')):
            os.remove(run_dir[index] + 'rpointer.ice')
        if(os.path.exists(run_dir[index] + 'rpointer.drv')):
            os.remove(run_dir[index] + 'rpointer.drv')
        if(os.path.exists(run_dir[index] + 'rpointer.atm')):
            os.remove(run_dir[index] + 'rpointer.atm')

        with open(run_dir[index] + 'rpointer.rof','w') as f:
                 f.write('./mycase.rtm.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.ocn','w') as f:
                 f.write('mycase.docn.r.' + start_date + '.nc\n')
                 f.write('mycase.docn.rs1.' + start_date + '.bin')
        with open(run_dir[index] + 'rpointer.lnd','w') as f:
                 f.write('./mycase.clm2.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.ice','w') as f:
                 f.write('mycase.cice.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.drv','w') as f:
                 f.write('mycase.cpl.r.' + start_date + '.nc')
        with open(run_dir[index] + 'rpointer.atm','w') as f:
                 f.write('mycase.cam.r.' + start_date + '.nc')
        if(not os.path.exists(run_dir[index] + 'mycase.cam.h0.' + output_date + '.nc')):
            self.write_to_pfile("output not found")
            raise Exception("output not found")
        fun_file = nc.Dataset(run_dir[index] + 'mycase.cam.h0.' + output_date + '.nc')
        try:
            fun_PSL2 = fun_file.variables['PSL'][duration - 1][:]
        except IndexError:
            self.write_to_pfile("error_value")
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
        value = basic -value_fun
        self.write_to_pfile("value = "+ str(value))
        self.write_to_pfile('--func done (%d)--'%(self.itern))
        # self.write_per(ini, value)
        # self.write_output(value)
        return value

    def testDelta(self):
        self.f,self.g = self.evalfg2(choice=2,x=self.x,f=self.f)
        for _ in range(10):
            self.DELTA = self.DELTA/10
            self.write_to_pfile("--DELTA = %f --"%(self.DELTA))
            self.g = self.evalfg2(choice=3,x=self.x,f=self.f)
    
    def testButterfly(self):

        current = 0
        xold = self.x
        while(True):

            # 检查是否有计算完成的case返回结果
            for busyCase in self.findBusy():
                val = self.completeTask(busyCase)
                if not val is None:
                    dim = self.ctl[busyCase]
                    self.ctl[busyCase] = -1
                    self.write_to_pfile("   Done! dim :" + str(dim) + " caseID :" +str(busyCase) + " val :" +str(val))


            # 检查是否有空闲case可供计算
            freeCases = self.findFree()
            # 结束条件，所有case空闲且所有梯度都已计算完成
            if(len(freeCases) == self.numprocs and current >= 20):
                self.write_to_pfile("all cases are free and all dims are done")
                break
            if(len(freeCases) > 0):
                for freeCase in freeCases:
                    # 防止溢出（去跑第51个case）
                    if current >= 20:
                        break
                    x,delta = self.createX2(xold,current)
                    self.submitTask(x,freeCase)
                    self.ctl[freeCase] = current
                    self.delta[freeCase] = delta
                    self.write_to_pfile("   Start! dim :" + str(current) + " caseID :" +str(freeCase) + " delta :" +str(delta))
                    current += 1
                    
            
            timelibrary.sleep(30)

        self.write_to_pfile("--test done (%d)--"%(self.itern))

        # 将所有进程重置为空闲状态
        for i in range(self.numprocs):
            self.ctl[i] = -1



reduction = []
reader = csv.reader(open('/BIGDATA1/iocas_mmu_3/support/reductions_T_w.csv', "r"))
for line in reader:
    lineArr = [float(n.strip()) for n in line]
    reduction.append(lineArr)

neg = [3.694162002055574, 2.3007709012494426, 4.8564804561266595, 4.971143796030497, 6.6417859699045785, -2.7814456395462903, 0.7354560896390208, 7.584109812333855, 4.999347166525484, 9.516322163733989, 0.9918733730045819, 1.899068376285662, 4.263410578352719, 2.67338268145171, 0.5891255224536921, -0.5273946133362423, 1.3420659347804285, 8.605162827867703, -1.9399323246716704, 7.635651247648104, 3.305485212205186, 0.4343931366877543, 4.263861838797711, 3.8129153563787903, 6.05191528727672, 0.7570705608916601, 0.5339783017709662, 3.633647504154749, 3.2836359036036127, 3.8366305560815626, -1.5335418930844278, 0.16657571615455743, 3.592546671730668, -1.8181997971891057, 5.942634891232685, 1.6339432339362756, 7.410652098422764, 3.907538962682855, -1.6389248398426584, 2.8762962228812405, 3.34335792203602, 5.109376191755419, 2.398497731454465, -0.0691665703460284, 4.15884432736265, -2.4337141032173233, 0.9693897374041158, 1.996891324207389, 6.971879735451466, 5.3292478661414435]
pos = [-1.7957734205160336, -2.4832257799548296, -6.518950693905748, -2.8716890588652095, -5.208429912253714, -1.8939349224180435, -2.211837703126254, -0.9219129507694348, -0.4203642737983837, -6.917585049477449, 1.0114212234853013, 4.309530142287708, -1.4154000343104505, 5.8887225960226735, -1.728333231354218, -0.3953315654001579, -1.311986426192005, -3.52117319655655, -2.768678539554144, 2.831811575392007, -0.9685200193952151, 2.91744554786807, -1.5589255333178458, -5.225812527578701, 0.7003509197313242, -2.5443162996837554, -5.090484528995176, -1.0200368279530765, -3.766734768839968, 3.8542476194380786, 1.8663605582645508, 3.201635529775394, 1.0106769168551193, 1.809287177408816, -4.481481330099257, -1.5731893849669127, -0.9298522171938908, -5.822129061587142, -2.7964681304752212, -2.0521914210352996, 3.0540597054958, -6.697304296331807, 3.4313218962524914, 2.473482385396798, 1.3423745178058568, 3.305570016478261, -5.14159557513977, -5.43097336446301, -2.709339432188396, -3.206242847697322]
my_spg2 = SPG2(dim = 50,numprocs = numprocs,ini=pos)
my_spg2.testButterfly()