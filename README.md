# PBRFOAM

A coupled solver of an in-house Particle Burning Rate (PBR) model and OpenFOAM for wildland fire spread modeling at flame scale. 

The solver is built on OpenFOAM version 4.x build number: dev-8dafde6048ab. 
See this repository for help installing this version: https://github.com/mmahmed15/OpenFOAM_UMD


# Installation instructions 

1. Download the source files from this repository in the home directory
```
cd ~
git clone https://github.com/msenamd/PBRFoam.git
```


2. Navigate to the compile directory and compile the code
```
cd ~/PBRFOAM
chmod -R +x *
./Allwmake -j 6 >& log.Allwmake &
```

check for errors in the log.Allwmake file 

# Reference Paper
Mohamed Mohsen Ahmed , Arnaud Trouvé, Jason Forthofer, and Mark Finney. Simulations of Flaming Combustion and Flaming-to-Smoldering Transition in Wildland Fire Spread at Flame Scale. Combustion and Flame 262 (2024) 113370. https://doi.org/10.1016/j.combustflame.2024.113370


