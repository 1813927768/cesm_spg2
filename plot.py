# make plots to valiadate final results for comparison

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as p
import random
from mpl_toolkits.basemap import Basemap, cm

direc ="D:\\repos\\cesm_spg2\\result\\"

def plot_nao(fileName,imgType):
    meteo_file = direc + fileName
    basic_file = direc +"NAOI_basic.nc"
    mean_file = direc + "PSLClm.nc"
    fh = Dataset(meteo_file, mode='r')
    basic = Dataset(basic_file, mode='r')
    mean = Dataset(mean_file, mode='r')

    #lons = fh.variables['lon'][:]
    lats = fh.variables['lat'][:]
    output = fh.variables['PSL'][14][:]
    if imgType == "mean":
        basic_data = mean.variables['PSLClm'][0][:]
    elif imgType == "basic":
        basic_data = basic.variables['PSL'][14][:]
    else:
        return
    # mean = mean.variables['PSLClm'][0][:]
    tlml = output - basic_data
    lons = np.arange(-180., 181.25, 1.25)

    fun_PSL = [[0.0 for i in range(0, 289)] for j in range(0, 192)]
    fun_PSL = np.array(fun_PSL)
    for i in range(0, 192):
        for j in range(0, 289):
            if j<144:
                fun_PSL[i][j] = tlml[i][j+144]
            else:
                fun_PSL[i][j] = tlml[i][j-144]
    # tlml_units = basic.variables['PSL'].units

    lon_0 = lons.mean()
    lat_0 = lats.mean()

    #m = Basemap(projection='ortho', lat_0=90, lon_0=0)
    #20N-80N, 90W-40E square region
    m = Basemap(lat_0=lat_0, lon_0=lon_0, llcrnrlat=20, urcrnrlat=80, llcrnrlon=-90, urcrnrlon=40, resolution='l')
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    tlml_0 = fun_PSL[::, ::]

    #m.drawparallels(np.arange(-90., 91., 20.), labels=[0,0,0,0], fontsize=10)
    #m.drawmeridians(np.arange(-180., 181.25, 40.), labels=[0,0,0,0], fontsize=10)

    #20N-80N, 90W-40E square region
    m.drawparallels(np.arange(20., 80.5, 20.), labels=[1,0,0,1], fontsize=12)
    m.drawmeridians(np.arange(-90., 40.5, 32.5), labels=[1,0,0,1], fontsize=12)

    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()

    #cm = plt.cm.get_cmap('coolwarm')
    #cs = m.pcolor(xi, yi, np.squeeze(tlml_0), cmap=cm, vmin=np.min(tlml_0), vmax=np.max(tlml_0))

    a = np.arange(-4000, 0, 500)
    b = np.arange(0, 4000, 500)

    clevs = np.append(a, b)
    cs = m.contourf(xi, yi, np.squeeze(tlml_0), clevs, cmap=p.cm.coolwarm)
    cbar = m.colorbar(cs, location='bottom', size="8%", pad="10%")
    # cbar.set_label(tlml_units, fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    plt.title("positive_unconstrained")
    plt.show()
    # plt.savefig(direc+imgType+"_"+fileName+".png")
    fh.close()


def plot_ini(fileName):

    meteo_file = direc + fileName
    basic_file = direc +"ini_basic.nc"
    mean_file = direc + "PSLClm.nc"
    fh = Dataset(meteo_file, mode='r')
    basic = Dataset(basic_file, mode='r')
    mean = Dataset(mean_file, mode='r')

    lons = fh.variables['lon'][:]
    lats = fh.variables['lat'][:]
    output = fh.variables['T_per'][0][:]
    basic_data = basic.variables['T_per'][0][:]
    mean = mean.variables['PSLClm'][0][:]
    tlml = output - basic_data
    # lons = np.arange(-180., 181.25, 1.25)

    # fun_PSL = [[0 for i in range(0, 289)] for j in range(0, 192)]
    # fun_PSL = np.array(fun_PSL)
    # for i in range(0, 192):
    #     for j in range(0, 289):
    #         if j<144:
    #             fun_PSL[i][j] = tlml[i][j+144]
    #         else:
    #             fun_PSL[i][j] = tlml[i][j-144]
    # tlml_units = basic.variables['PSL'].units

    lon_0 = lons.mean()
    lat_0 = lats.mean()

    m = Basemap(projection='ortho', lat_0=90, lon_0=0)
    #20N-80N, 90W-40E square region
    # m = Basemap(lat_0=lat_0, lon_0=lon_0, llcrnrlat=20, urcrnrlat=80, llcrnrlon=-90, urcrnrlon=40, resolution='l')
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    tlml_0 = tlml[::, ::]

    m.drawparallels(np.arange(-90., 91., 20.), labels=[0,0,0,0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181.25, 40.), labels=[0,0,0,0], fontsize=10)

    #20N-80N, 90W-40E square region
    # m.drawparallels(np.arange(20., 80.5, 20.), labels=[1,0,0,1], fontsize=12)
    # m.drawmeridians(np.arange(-90., 40.5, 32.5), labels=[1,0,0,1], fontsize=12)

    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()

    cm = plt.cm.get_cmap('coolwarm')
    cs = m.pcolor(xi, yi, np.squeeze(tlml_0), cmap=cm, vmin=np.min(tlml_0), vmax=np.max(tlml_0))

    # a = np.arange(-3500, 0, 500)
    # b = np.arange(0, 4000, 500)
    # clevs = np.append(a, b)
    # cs = m.contourf(xi, yi, np.squeeze(tlml_0), clevs, cmap=p.cm.coolwarm)

    cbar = m.colorbar(cs, location='bottom', size="8%", pad="10%")
    cbar.ax.tick_params(labelsize=12)
    plt.show()
    fh.close()


if __name__ == '__main__':
    # plot_nao(fileName = "NAOI_0.9377883543892195.nc",imgType="mean")
    plot_nao(fileName = "NAOI_-0.21755857522757985.nc",imgType="basic")
    # plot_ini(fileName="ini_-1.2563443051793428.nc")