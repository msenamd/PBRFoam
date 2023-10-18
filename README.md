# PBRFOAM

A coupled solver of an in-house Particle Burning Rate (PBR) model and OpenFOAM for wildland fire spread modeling at flame scale. 

The solver is built on OpenFOAM version 4.x build number: dev-8dafde6048ab. 
See this repository for help installing this version: https://github.com/mmahmed15/OpenFOAM_UMD


# Installation instructions 

1. Download the source files from this repository in the home directory
```
cd ~
git clone https://github.com/mmahmed15/PBRFoam.git
```


2. Navigate to the compile directory and compile the code
```
cd ~/PBRFOAM
chmod -R +x *
./Allwmake -j 6 >& log.Allwmake &
```

check for errors in the log.Allwmake file 
