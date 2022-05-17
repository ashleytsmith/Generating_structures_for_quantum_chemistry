# Generating structures for quantum chemistry

## Summary


Python module for generating different zeolite structures as the input for quantum chemistry calculations. Generates all possible neighboring acid sites given the first one is already fixed. This gif shows a few of the generated structures with the neighbor distance getting larger and larger.

<p align="center">
<img src="https://github.com/ashleytsmith/Generating_structures_for_quantum_chemistry/blob/main/Images_for_GitHub/Neighbours.gif" width="400" alt="movie of some of the genrated structures"> 
</p>

## Background chemistry

*	Zeolites are aluminosilicate materials made up of tetrahedra of Si/Al and oxygen. 
*	Swapping a Si for and Al atom and balancing the charge with a proton produces an acid site (see the first picture below).

<p align="center">
<img src="https://github.com/ashleytsmith/Generating_structures_for_quantum_chemistry/blob/main/Images_for_GitHub/First_acid_site.png" width="400" alt="First acid site">

  
 <img src="https://github.com/ashleytsmith/Generating_structures_for_quantum_chemistry/blob/main/Images_for_GitHub/Add_next_nearest_neighbour_site.png" width="400" alt = "Illustration of two sites"> 
 
  </p>

*	The structure CHA has 36 tetrahedra which means 35 places to put a second Al atom and 35 by 4 ways to add a second site (see the second picture above for an example of a two-site structure).
*	The module generate_sites.py recursively walks through starting from the first site and categorises the neighbors into nearest neighbor, next nearest neighbors and so on which are labelled neighbours_0, neighbours_1 in the file structure respectively.

# Requirements


This module works with atomic simulation environment and a minimum you will need to install 

* Python 3.6 or above
* atomic simulation environment
* numpy
* scipy

There is a good guide on how to do this in the ase documentation. If you prefer using anaconda and virtual environments like me it is also possible to install ase via anaconda which also has great documentation.


```
conda install -c conda-forge ase
```



