# SCAN OVERLAP MANAGEMENT 
## dependencies
DGtal version 1.4 or later, see [DGtal installation] (https://github.com/DGtal-team/DGtal).
## compilation instructions
```
mkdir build
cd build
cmake ..  -DDGtal_DIR=/path/to/DGtal
make
```
## example of use
You can use the code on a sample provided in the data directory. For example, for the WSPS4 sample, the command to type in a terminal is : 
```
./logRec -i ../data/RGB_R_WSPS4_Cloud.xyz -c ../data/WSPS4_centerline.xyz -p ../data/WSPS4_clean.id  -o test
```
This command generates two files: test_rec.id and test_rec_colorByScan.xyz

- **test_rec.id** : Contains the identifiers of the points after applying the scan overlap management algorithm. 
- **test_rec_colorByScan.xyz** : Contains the resulting point cloud with a color associated with each scan.