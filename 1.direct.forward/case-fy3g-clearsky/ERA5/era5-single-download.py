#%%
import sys  
import os
import cdsapi
if len(sys.argv) <= 1:
    print("请传递参数！")
    exit()
else:
    print("download- YYYYMMDD", sys.argv[1]," UTC: ", sys.argv[2])

yyyymmdd= sys.argv[1]
UTC= sys.argv[2]

if not os.path.exists("ERA5/"+yyyymmdd):
    os.mkdir("ERA5/"+yyyymmdd)
c = cdsapi.Client()
c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', #'2m_dewpoint_temperature',
            '2m_temperature', 'sea_surface_temperature', 'surface_pressure',# 'total_cloud_cover',           
        ],
        'year': yyyymmdd[0:4],
        'month': yyyymmdd[4:6],
        'day': yyyymmdd[6:8],
        'time': UTC+':00',
        'area': [65, -180, -65, 180,],
        # 'area': [30, 120, 20, 130,],
        'grid': ['0.25','0.25'],
        'format': 'netcdf',
    },
    'ERA5/'+yyyymmdd+'/ERA5-Single-GBL-'+yyyymmdd+'-'+UTC+'00.nc')
