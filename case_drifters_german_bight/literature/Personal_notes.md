22.11.03
Random walk is composed of a sequence of random steps.
Research question: work towards a 3D particle simulation of salt in the Rhine Meuse delta (with the RMM model)
Challenges: 
	- adapt particle-models to large-scale parallel computation, i.e. with hundreds of domains(processors); 
	- model the surpressed vertical mixing;
	- application of the model to understand the physics of the intrusion phenomena
Recent work on literature review:
	- application of particle model to track drifters:Ruiter(2020) & Meyerjürgens(2018)
	- general knowledge about a particle model: Sebille et. al.(2018) and several other reports from previous student
	- components of the particle model: Particles.jl 
My words on random walk model: velocity = drift + diffusion coefficient* brownian motion velocity


Question in bachelor thesis of N. Scheijen
- equivalence between the random walk model and the diffusion-advection equation; 
	- equivalence obtained when sampled enough + relate the diffusion parameter in both equations
	- RW: mass conserving; a large amount of particles needed, accumulated error or instability for long time or large space
	- DA: suitable for large scale; not mass conserving or even negative concentration can occur
	- 
- section 2.2.4: relation between particle density(unit mass?) and PDF of particles? --> enough particles?
- section 2.1.2: near wall turbulence - eddy-viscosity model - eddy viscocity parabolic in vertical direction --> for exploration only
- equation(2.7): why only two directions instead of three? --> only 2D is considered
- equation(2.13): where does it come from? --> evolution of Ito process?
- corrected (2.22): https://www.wolframalpha.com/input?i=dv%2F%28a-bv%5E2%29%3D+dt%2Fc
- added mass in section 2.3
 
22.11.10
What is a NetCDF file?
A format for data storage for scientific computing.
https://rwoconne.github.io/rwoclass/astr511/IDLresources/idl_5.1_html/dftoc.htm#pgfId=57020
What's contain in a NetCDF file?
Variables, Attributes and Dimensions
Variables that share the same name with dimensions, e.g. time, space e.t.c.
Variables can have attributes.
Load .nc file
_______
using NetCDF
ncinfo(filename)
x = ncread(filename, variablename)
_______

_______
using NetCDF
x = nccreate(filename, variablename, "dim1name",...,"dim2name",...,atts=global_attributes)
ncwrite(var, filename, variablename)
_______
Note that multiple variables with their own dimension can be created at the same time.
When "dim1name" is followed by an integer, only the named integer(dimension) is created.
When "dim1name" is followed by an 1D array, it does three thing. Set the named dimension as the length of the array, create a variable that stores the 1D array and also create a attribute with missing value.

What's included in the source data file?
In the German Bight directory,
julia> filename = "data\\dcsm-fm_201703\\DCSM-FM_05nm_0001_map.nc"
"data\\dcsm-fm_201703\\DCSM-FM_05nm_0001_map.nc"

julia> filepath = joinpath(@__DIR__, filename)
"I:\\Master_Thesis\\particles.jl\\case_drifters_german_bight\\data\\dcsm-fm_201703\\DCSM-FM_05nm_0001_map.nc"

julia> ncinfo(filename)