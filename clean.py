# clean the case directory for a restart

import os

nums = 5

base_dir = ['/BIGDATA1/iocas_mmu_3/cesm/1.2.2/ice/14/parallel/' + str(i) + '/mycase' for i in range(nums)]
run_dir = [route + '/run/' for route in base_dir]
fin_dir = [route + '/finished' for route in base_dir]

start_date = '0052-11-01-00000'

for index in range(nums):
    if os.path.exists(fin_dir[index]):
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
    
    os.system("find %s -name '*log*' -delete"%(run_dir[index]))