# Generating structures for quantum chemistry

## Summary


Python module for generating different zeolite structures as the input for quantum chemistry calculations. Generates all possible neighbouring acid sites given the first one is already fixed. This gif shows a few of the generated structures with the neighbour distance getting larger and larger.

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
*	The module generate_sites.py recursively walks through starting from the first site and categorises the neighbours into nearest neighbour, next nearest neighbours and so on which are labelled neighbours_0, neighbours_1 in the file structure respectively.

## Requirements


This module works with atomic simulation environment and a minimum you will need to install 

* Python 3.6 or above
* atomic simulation environment
* numpy
* scipy

There is a good guide on how to do this in the ase documentation. If you prefer using anaconda and virtual environments like me it is also possible to install ase via anaconda which also has great documentation.


```
conda install -c conda-forge ase
```

## Summary of the code

**recursive_neighbours()**	Walks through the structure starting from the first site and categorises a given Si/Al atom in terms of how many neighbours it is away e.g. nearest neighbours are connected by one oxygen i.e. Site1-O-Site2. Second nearest neighbours are separated by another intermediate tetrahedra i.e. Site1-O-Si-O-Site2 and so on. 

**find_neighbours()**		Finds the nearest neighbours of a given atom and checks 
the expected result is attained.

**expected_neighbours()**	Returns expected number of neighbours for a given atom type.

**replace_Si_with_Al()**	Replaces a Si atom with Al keeping original the indices intact.

**add_hydrogen()** Adds the hydrogen on the outer edge of the triangle made by two Si and an oxygen.


**make_files_by_Si_to_Si_distance()** Make the top layer of folders.

**make_subfolders_for_neighbours()**	Make folders for each silicon and another subfolder for its oxygen neighbours. Folders are named with the indices from the input file.

**change_folder()**	Change the folder and create a folder if the path doesn’t exist already.

**check_for_fails()**   Returns the indices for which an incorrect number of neighbours are found.


**Input parameters and required files**:

number_of_atoms_in_structure = 108 <br>
oxygen_of_the_first_acid_site =  107 *index of the first acid site that we find the neighbours of*. <br>
Al_of_the_first_acid_site = 7 <br>

**Geometry and distance criteria**

cut_off_distance = 2<br>
silicon_expected_neighbours = 4<br>
oxygen_expected_neighbours = 2<br>
hydrogen_expected_neighbours = 1<br>



