# to get input for MEOWSS TMRCA put into file: active_zars.txt

stubs = [('/ZARR_NEW/0.5/100/', '_13_1000_0.00025_2.5e-05_0.5_100.vcz'),
         ('/ZARR_NEW/0.7/100/', '_13_1000_0.00025_2.5e-05_0.7_100.vcz'),
         ('/ZARR_NEW/0.5/100/', '_14_1000_0.00025_0.00025_0.5_100.vcz'),
         ('/ZARR_NEW/0.7/100/', '_14_1000_0.00025_0.00025_0.7_100.vcz'),
         ('/ZARR_NEW/0.5/100/', '_15_1000_0.00025_0.0025_0.5_100.vcz'),
         ('/ZARR_NEW/0.7/100/', '_15_1000_0.00025_0.0025_0.7_100.vcz'),
         ('/ZARR_NEW/0.5/100/', '_16_1000_0.0025_2.5e-05_0.5_100.vcz'),
         ('/ZARR_NEW/0.7/100/', '_16_1000_0.0025_2.5e-05_0.7_100.vcz'),
         ('/ZARR_NEW/0.5/100/', '_17_1000_0.0025_0.00025_0.5_100.vcz'),
         ('/ZARR_NEW/0.7/100/', '_17_1000_0.0025_0.00025_0.7_100.vcz')
         ]
seeds = range(1,11)
points = range(10, 110, 10)
thresholds = [0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97]
window=100
with open("active_zars.txt", "w") as file:
    file.write('filename,threshold,points,window\n')
    
    for seed in seeds:
        for folder_stub, filename_stub in stubs:
            for point in points:
                for threshold in thresholds:
                    file.write(f'{seed}{folder_stub}{seed}{filename_stub},{threshold},{point},{window}\n')

