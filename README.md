# isingQuench

*Author:* Adam Iaizzi  
*email:* iaizzi@bu.edu  
*website:* https://www.iaizzi.me  
*gitlab repo:* https://gitlab.com/adazi/afm-ising2d  
Copyright Adam Iaizzi 2020

This code is a Metropolis algorithm Monte Carlo simulation of the 2D Ising antiferromagnet with an external magnetic field. This code is a modified version of `ising2d.f90` originally written by Prof. Anders W. Sandvik [source](http://physics.bu.edu/~py502/lectures5/examples/index.html), used with permission. 

## Quick Start

Clone the repository to your machine:  
```bash
git clone https://gitlab.com/adazi/afm-ising2d.git
```

Navigate to the source code directory:  
```bash
cd afm-ising2d/
```

To compile:  
```bash
gfortran -O3 ising.f90
```

To run:  
```bash
./a.out
```

Example output to stdout: 
```
 h=   1.0000000000000000       t =    1.0000000000000000E-002
 Bin            1  of           10
 Bin            2  of           10
 Bin            3  of           10
 Bin            4  of           10
 Bin            5  of           10
 Bin            6  of           10
 Bin            7  of           10
 Bin            8  of           10
 Bin            9  of           10
 Bin           10  of           10
 ```

## Prerequisites

This program is pretty simple and self contained. The only prerequisite is a fortran compiler. I have used gfortran. 

## Program Files

`ising.f90` --- Calls routines to run Monte Carlo and coordinates MPI communication between replicas, also reads and writes to files.  
`read.in` --- Options for almost all possible user-defined behavior, described later.  
`seed.in` --- Random seed input file. Random seed should be a large integer. Should be automatically updated when the program finishes. **Note:** may not be properly updated if program does not exit normally.  


## Input files: read.in

The main user input file is `read.in`. All parameters are required and must be provided in this order, separated by spaces. 

```
ll tt hh
steps1 steps2 bins
```

`ll` -- Integer  
Linear size of system (will be ll x ll)

`tt` -- Real  
Temperature in units of T/J

`hh` -- Real  
Magnetic field strength, in units of h/J. Negative values not guaranteed to work. 

`steps1` -- Integer  
Number of initial equilibration Monte Carlo sweeps to be performed before beginning measurements.  
Suggested value: 10,000

`steps2` -- Integer  
Number of Monte Carlo sweeps per bin.  
Suggested value: 10,000

`bins` -- Integer  
Number of MC bins to be preformed.  
Suggested value: >10

### An example read.in file

As an example, this read.in file will do a simulation of a 16x16 system at T=1.5 and h=3 with 10,000 equilibration steps followed by 10 bins with 10,000 MC sweeps each. 

```
16 1.5d0 3.d0
10000 10000 10
```

Note that the format `1.5d0` gives 1.5x10^0. This format makes it explicit (to Fortran) that the number is a real number. 


## Output Files

This program produces a number of output files. All data files end in `.txt` and are formatted as simple ASCII tab-separated rectangular arrays for easy importation into a post processing program.

The first column of most output files is temperature followed by one or more measured quantities. *There are no labels.* One line for each bin. 

The program appends to existing files when run again. The user must manually delete old data if they want to start fresh. 

### Details of specific files:

Below I will describe what will appear on each line of the output files. 

```
enrg.txt
--> energy
--> Real, extensive
--> 3 columns
T	<E> 	<E^2>

amag.txt
--> uniform magnetization
--> 4 col
T		<m>		<m^2>		<m^4>

smag.txt
--> staggered magnetization
--> 4 col
T		<sm>		<sm^2>		<sm^4>

spins.txt
---> Full spin configuration, tab-separated +1 -1...
---> One configuration per line 
---> No averaging
---> By default: written once per bin. 
---> Caution: Will generate large output files for large systems. 
s_0		s_1		s_2		s_3 .... s_{n-1}

census.txt
---> Number of spins of each operator type C^y_x
---> Also includes the sum of probability of flipping each spin
---> By default, only one measurement per bin
---> If y=sum of neighboring spins, x=center spin
---> The state has index y if x=-1 and y+1 if x=+1
---> sumFlips is the sum of the flip probabilities of all spins in the system. When this number is <<1, the system is frozen. 
T P[-4] P[-3] P[-2] P[-1] P[0] P[+1] P[+2] P[+3] P[+4] P[+5] sumFlips

```


### Preventing outputs from using up too much disk space

Writing the full spin state out to disk can be time consuming, and for large systems, require large amounts of storage. As such, users should be cautious about running the program with large sizes and especially with large numbers of bins. The number of entries on each line is `L^2`, and there are `bins`, so as a rough estimate, the space required will be `16*L*L*bins` bytes. 

To turn off the writing out to `spins.txt` the user can just comment out lines 158-160: 

```fortran
!write spin configuration out to file at end of bin
open(10,file='spins.txt',position='append')
write(10,*) spin(:)
close(10)
``` 

**Tip:** Since this ACSII text format is extremely in efficient, compressing the `spins.txt` file using `gzip` can reduce the space required by more than a factor of 10. 

## Performing a quench

Reproducing the freezing process in [my paper](https://arxiv.org/abs/2001.09268) with this code is quite simple. If we set the program up with L=32, T=0.01, h=1 and run it, it will quickly reach a frozen state. The problem is that only one such state will be reached, and then the program is stuck there forever. 

We will need to run the program many times to average over many final states. That can be done with a simple bash script: 
```bash
for x in {1..200}
do
    cat seed.in #check that random seed is different each time
    echo "Running " $x " of " 2000
    time ./a.out > log.dat
done
```

Fortunately, since we are only interested in the final state of this individual simulations can be very simple. We need only a little bit of equilibration time and then we can make one measurement and move on. 

Example `read.in` file:
```
32 0.01d0 1.0d0
20000 20000 1
``` 

Once can then do the averaging in a separate program. 


## References

This code was used for a paper on quenches in the 2D Ising antiferromagnetic. Available on [arXiv](https://arxiv.org/abs/2001.09268)

The [original program](http://physics.bu.edu/~py502/lectures5/examples/index.html) on which this code was based was part of Anders Sandvik's [PY502 course](http://physics.bu.edu/~py502/) at Boston University. 

[My dissertation](https://www.springer.com/us/book/9783030018023)  

[Author's website](https://www.iaizzi.me)  



## Known Issues

N/A
