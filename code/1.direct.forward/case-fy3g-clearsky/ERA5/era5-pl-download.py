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
    'reanalysis-era5-pressure-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            'fraction_of_cloud_cover', 'geopotential', 'relative_humidity',
            'specific_cloud_ice_water_content', 'specific_cloud_liquid_water_content', 'specific_humidity',
            'specific_rain_water_content', 'specific_snow_water_content', 'temperature',
        ],
        'pressure_level': [
            '50', '70', '100',
            '125', '150', '175',
            '200', '225', '250',
            '300', '350', '400',
            '450', '500', '550',
            '600', '650', '700',
            '750', '775', '800',
            '825', '850', '875',
            '900', '925', '950',
            '975', '1000',
        ],
        'year': yyyymmdd[0:4],
        'month': yyyymmdd[4:6],
        'day': yyyymmdd[6:8],
        'time': UTC+':00',
        'area': [65, -180, -65, 180,],
        # 'area': [30, 120, 20, 130,],
        'grid': ['0.25','0.25']
    },
    'ERA5/'+yyyymmdd+'/ERA5-PL-GBL-'+yyyymmdd+'-'+UTC+'00.nc')
    # 'ERA5-PL-GBL-'+year+month+day+time+'.nc')

# %%
