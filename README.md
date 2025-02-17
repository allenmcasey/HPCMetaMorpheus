HPCMetaMorpheus

This is a C++ version of the MetaMorpheus library developed for HPC environments. 

The original C# version was developed by the L.M. Smith Group https://github.com/smith-chem-wisc/
The code has been translated to C++ by the PSTL group at the University of Houston, 
and parallelized using OpenMP and MPI. The software is distributed under a MIT license.

The repository contains 3 branches, the main (master) sequential branch with the C++ code, a multi-threaded OpenMP based source (topic/OpenMP), 
and an MPI version for clusters (topic/MPI). As of right now (July 2021), only XLSeaarchTask has received significant testing and considered reasonably stable, although bugs and issues most likely still exist.  SearchTask and CalibrationTasks are work in progress should be available in the near future.

The syntax used by HPCMetaMorpheus follows closely the syntax of MetaMorpheus on Linux systems. Please refer to https://github.com/smith-chem-wisc/MetaMorpheus/wiki for details. An example for a Crosslink SearchTask would be:

```
./HPCMetaMorpheus -t Task1XLSearchTask.toml -s Input.mgf -d database.fasta
```

HPCMetaMorpheus supports at the moment mgf and mzml Spectra files.  
