# Razorback
Astrodyanmics Research Code at The University of Kansas
Author: Daniel Owen (daniel.owen@ku.edu)

The main purpose of this code is to collect functions and data structures as well as driver and analysis scripts for modeling and solving the multiple gravity assist problem. This research is focused on comparing both optimization algorithms and problem formulation to gain further insight into the problem. 

For further information of this type of work and where the idea originated please see my paper: 
    
    Owen, D., Rhodes, Z., Kaplinger, B., "A Uranus Mission Design Demonstrating a Simulated Annealing Algorithm", AAS/AIAA Spaceflight Mechanics Meetings, AAS 23-371

# Current Development
This is a brand-new repository (May 14, 2024) with the current goal to set up the repo and development environments and begin porting code over from previous work

# Conda
- python 3.12
- plotly
- numpy
- scipy

# Future Development
The near-term goals are to set up an early version of the first test problem (outer planet transfer) and hook it up to an early version of the first optimizer (Particle Swarm Optimizer). An interface also needs to be developed between the problem definition (and cost function) and optimization algorithm. 

From there the test problems as well as optimizers will be built out. Eventually problems and optimizers will be run in parallel and analysis scripts will be developed to perform comparisons and dive into the results of each combination

# Name
The Razorback is the name of a ship in the book series by James S.A. Corey The Expanse. The ship is also featured in the television adaptation from sci-fi (and later Amazon Prime). 

The ship is a racing yacht owned by billionaire Something Mao and piloted by his daughter Julie Mao. In the world of The Expanse these races are not on oceans or lake but in space, racing between planets or around the moons of the gas giants. As such the Razorback is fitted by a large rocket motor (an Epstein drive in the lore of the novels) for propelling around the solar system. The races however are done ballistically, by planning gravity assist maneuvers using cunning astrodynamics, navigation, and optimization. As such the name razorback felt apt for this repository as the tools and analysis developed here may one day be used by daring pilots and astronauts in the future to race similar vessels. 