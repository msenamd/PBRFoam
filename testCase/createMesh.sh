#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:0:0           # Run time (hh:mm:ss)
#SBATCH --partition=parallel
#SBATCH --account=atrouve
#SBATCH --mail-user=mmahmed@umd.edu
#SBATCH --mail-type=begin     # email when the job starts
#SBATCH --mail-type=end       # email when the job finishes

module load gcc/5.5.0   
module load openmpi/3.1

rm -rf constant/polyMesh
rm -rf system/topoSetDict
rm -rf *.obj

blockMesh

snappyHexMesh -overwrite

###extrudeMesh
###transformPoints -translate "(0 0 5.05)"

\cp -rf system/topoSetDict_grass system/topoSetDict
\cp -rf system/createPatchDict_grass system/createPatchDict  
topoSet
createPatch -overwrite

\cp -rf system/topoSetDict_burner system/topoSetDict
\cp -rf system/createPatchDict_burner system/createPatchDict  
topoSet
createPatch -overwrite

rm -rf *.obj

checkMesh | tee log.checkMesh

setParticlesSeparation
