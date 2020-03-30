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

start_date = '0052-11-01-00000'
output_date = '0052-11-02-00000'
base_dir = ['/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/' + str(i) + '/mycase' for i in range(numprocs)]
run_dir = [route + '/run/' for route in base_dir]
fin_dir = [route + '/finished' for route in base_dir]
# base_dir = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/' + str(rank) + '/mycase/'
# run_dir = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/' + str(rank) + '/mycase/run/'
# fin_dir = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/' + str(rank) + '/mycase/finished'
file_name = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/result_min_1110.txt'
process_file = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/plot.txt'
restart_file = '/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/restart.txt'
duration = 15
basic = 1.58523
init_value = 1e10
error_value = 10

# 所有进程初始PT是一样的
pre_file = nc.Dataset(run_dir[0] + 'mycase.cam.r.' + start_date + '.nc', mode='a')
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
normp = 3 #暂不考虑约束问题
M = 10    #lineSearch考虑的之前的特征方程值个数
RMAX = 1e30
RMIN = 1e-30
MAXIT = 50
MAXIFCNT = 480
eps = 1e-6 #参照fortran代码

class SPG2():
    def __init__(self, dim, numprocs):

        self.dim = dim
        self.numprocs = numprocs
        
        # itern ifcnt igcnt 迭代次数 特征方程计算次数 梯度计算次数
        self.ifcnt = self.itern = self.igcnt = 0
        self.x =  [0 for i in range(self.dim)]
        self.g = [0 for i in range(self.dim)]
        self.f = init_value
        self.fvalues = [-1.0e+99 for i in range(M)]
        self.fvalue = [0 for i in range(MAXIT)]
        self.fcnorm = [0 for i in range(MAXIT)]
        # 用来控制伪并行的进行（-1表示对应case空闲,i表示此时正在计算某个维度的偏微分值）
        self.ctl = [-1 for i in range(self.numprocs)]
        # 用来暂时存储delta值
        self.delta = [0 for i in range(self.numprocs)]
        # 随机初始化x位置
        for i in range(self.dim):
            self.x[i] = random.uniform(1,10)

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

    # 重入2
    def restart2(self):
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
                elif head == "xf":
                    self.f = body
                elif head == "gf":
                    self.f = body

                
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

    # normp = 1:  box  constraints
    # normp = 2:  ball constraints
    # delta > 0:  parameter in the constraint
    # 将扰动x投影到特征空间
    def project(self,x):
        delta = deltm = 0.0
        # box constraint
        if(normp == 1):
            x = [min(delta,max(i,deltm)) for i in x]
        # ball constraint
        elif(normp == 2):
            sum = 0
            for i in x:
                sum += i*i
            sqrtsum = math.sqrt(sum)
            if(sqrtsum <= delta):
                self.write_to_pfile('proj done')
                return x
            else:
                c = delta/sqrtsum
            x = [c*i for i in x]
        self.write_to_pfile('proj done')
        return x

    # 将扰动x投影到特征空间（修改）
    def project2(self,x):
        ini = np.dot(np.array(x),np.array(reduction))
        self.write_to_pfile('--proj start--')
        self.write_to_pfile(x)
        ini_T = ini[0 : 32 * 288]
        funt = self.funT(ini_T)
        self.write_to_pfile('f = ' + str(funt))
        # 约束
        if (funt > 100):
            # x = 100 / funt * x
            x = [100 / funt * i for i in x ]     
            self.write_to_pfile('funt = '+str(funt))  
        self.write_to_pfile('--proj done--')
        self.write_to_pfile(x)
        return x

    # 通过lineSearch更新扰动位置信息
    def lineSearch(self, d, gtd):
        self.write_to_pfile("--line search start--")
        self.write_to_pfile("old f = "+str(self.f))
        self.write_to_pfile("gtd = "+str(gtd))
        self.write_to_pfile("d = "+str(d))
        
        GAMMA =  1.0e-04
        x = self.x
        f = self.f
        fvalues = self.fvalues

        self.write_to_pfile("x = "+str(x))
        # 计算fmax
        fmax = fvalues[0]
        for i in range(1,M):
            fmax = max(fmax,fvalues[i])
        
        # 计算xnew新的扰动
        xnew = [x[i]+d[i] for i in range(len(x))]
        self.write_to_pfile("xnew = "+str(xnew))

        # 计算特征方程值
        fnew = self.evalfg2(x=xnew,f=None,choice=1)
        self.write_to_pfile('second call evalfg')

        alpha = 1.0

        #主循环
        while(fnew > fmax + alpha * GAMMA * gtd): #gtd对应δ，alpha对应λ
            if(alpha <= 0.1):
                alpha = alpha/2
            else:
                atemp = (-gtd*alpha**2)/(2.0*(fnew-f-alpha*gtd))
                if(atemp < 0.1 or atemp > 0.9*alpha):
                    atemp = alpha/2
                alpha = atemp
            
            xnew = [x[i]+alpha*d[i] for i in range(self.dim)]

            self.write_to_pfile('alpha = '+ str(alpha))
            self.write_to_pfile('xnew = '+ str(xnew))
            self.write_to_pfile('fmax + alpha * GAMMA * gtd = '+ str(fmax + alpha * GAMMA * gtd))

            # 计算特征方程值
            fnew = self.evalfg2(x=xnew,f=None,choice=1)


            self.write_to_pfile('fnew = '+ str(fnew))
            self.write_to_pfile('third call evalfg')

       # 更新扰动位置信息
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

        self.write_to_pfile("--line search end--")
        self.write_to_pfile("new f = "+str(self.f))

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
        deltall = abs(x[dim]*0.01)
        if(deltall == 0.0):
            deltall = 0.01
        x[dim] += deltall
        return x,deltall


    # 求梯度（伪并行）
    def gradient(self,xold,f):

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
                    self.write_to_pfile("gradient done current dim :" + str(dim) + " current caseID :" +str(busyCase))


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
                    self.write_to_pfile("gradient start current dim :" + str(current) + " current caseID :" +str(freeCase))
                    current += 1
                    
            
            timelibrary.sleep(30)

        self.write_to_pfile("gradient done")

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


            # with open(restart_file, 'r') as f:
            #     for line in f.readlines():
            #         parts = line.split("=")
            #         head = parts[0].strip()
            #         body = json.loads(parts[1])
            #         if head == "iter":
            #             self.ifcnt = self.itern = self.igcnt = body
            #         elif head == "x":
            #             self.x = body
            #         elif head == "g":
            #             self.g = body
            #         elif head == "xf":
            #             xf = body
            #         elif head == "gf":
            #             gf = body
            #         elif head =="fvalue":
            #             for i in range(len(body)):
            #                 self.fvalue[i] = body[i]
            #         elif head == "fvalues":
            #             for i in range(len(body)):
            #                 self.fvalues[i] = body[i]
            #         elif head == "fcnorm":
            #             for i in range(len(body)):
            #                 self.fcnorm[i] = body[i]
            #         elif head == "fbest":
            #             self.f = self.fbest = body
            #         elif head == "xbest":
            #             self.xbest = body
            
            # s = [0 for i in range(self.dim)]
            # y = [0 for i in range(self.dim)]
            # cg = [0 for i in range(self.dim)]
            # sts = sty = 0.0
            # for i in range(self.dim):
            #     s[i] = self.x[i]-xf[i]
            #     y[i]= self.g[i]-gf[i]
            #     sts = sts + s[i]*s[i]
            #     sty = sty + s[i]*y[i]
            #     cg[i] = self.x[i]-self.g[i]
            

            # self.write_to_pfile('proj cg2')
            # cg = self.project(cg)
            
            # cgnorm=0

            # for i in range(self.dim):
            #     cgnorm = max(cgnorm,abs(cg[i]-self.x[i]))
            
            # if(sty <= 0):
            #     self.rambda = RMAX
            # else:
            #     self.rambda = min(RMAX,max(RMIN,sts/sty))

        # 第一次开始
        else:
            #将初始扰动映射到特征空间内
            self.x = self.project2(self.x)
            self.xbest = self.x
            self.write_to_pfile("project done, start evalfg2")

            #求梯度和特征方程值
            self.f,self.g = self.evalfg2(choice=2,x=self.x,f=self.f)
            self.write_to_pfile('first call evalfg')

            self.fbest = self.fvalue[0] = self.fvalues[0] = self.f

            cg = [self.x[i] - self.g[i] for i in range(self.dim)]
            self.write_to_pfile('proj cg')
            cg = self.project2(cg)

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

            gtd = 0
            d = [self.x[i]-self.rambda*self.g[i] for i in range(self.dim)]
            self.write_to_pfile('proj d')
            d = self.project2(d)
            for i in range(self.dim):
                d[i] = d[i]-self.x[i]
                gtd += self.g[i]*d[i]

            #通过lineSearch更新扰动x和特征方程值f
            xnew = self.lineSearch(d,gtd)

            #计算更新后的梯度
            gnew = self.evalfg2(choice=3,x=xnew,f=self.f) 
            self.write_to_pfile('fouth call evalfg')

            s = [0 for i in range(self.dim)]
            y = [0 for i in range(self.dim)]
            sts = sty = 0.0
            for i in range(self.dim):
                s[i] = xnew[i]-self.x[i]
                y[i]= gnew[i]-self.g[i]
                sts = sts + s[i]*s[i]
                sty = sty + s[i]*y[i]
                self.x[i] = xnew[i]
                self.g[i] = gnew[i]
                cg[i] = self.x[i]-self.g[i]
            
            self.write_to_pfile('proj cg2')
            cg = self.project2(cg)

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
        self.write_to_pfile('restore finish.')
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
        self.write_to_pfile('start ok')
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
        self.write_to_pfile("break ok")
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
            return error_value
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
        value = value_fun - basic
        self.write_to_pfile("value = "+ str(value))
        return value

    # 求目标函数值
    def function(self, x, index):
        level = 26
        row = 192
        column = 288
        pre_file = nc.Dataset(run_dir[index] + 'mycase.cam.r.' + start_date + '.nc', mode='a')
        ini = np.dot(np.array(x),np.array(reduction))
        self.write_to_pfile('start ok')
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
        self.write_to_pfile("break ok")
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
            return error_value
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
        value = value_fun - basic
        self.write_to_pfile("value = "+ str(value))
        self.write_per(ini, value)
        self.write_output(value)
        return value



reduction = []
reader = csv.reader(open('/BIGDATA1/iocas_mmu_3/support/reductions_T_w.csv', "r"))
for line in reader:
    lineArr = [float(n.strip()) for n in line]
    reduction.append(lineArr)
my_spg2 = SPG2(dim = 50,numprocs = numprocs)
my_spg2.spg2()
# my_spg2.function([0 for _ in range(50)],0)
