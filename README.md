# Particles.jl
Particle modelling tools for coastal flow

# package development on windows
- add your julia development folder to the LOAD_PATH using `push!(LOAD_PATH, "some_drive:\\some_folder\\julia"`. This can be done on the prompt or in your startup.jl file in folder c:\Users\username\.julia\config. Create the config folder if it does not exist.  
- activate Particles  
cd to folder eg `cd("d:\\some_folder\\julia\\Particles")`  
`using Pkg`  
`Pkg.activate(".")`  
or activate the julia package management environment with the ']'-key on the REPL and use commands like `activate .`  
- update packages
If you expect missing packages on you machine then run  
`Pkg.update()`  
- running unit tests
The unit tests reside in tests. The main script is runtests.jl. To run them type:  
`Pkg.test("Particles")`  

## example of using the package
`cd("d:\\some_folder\\julia\\Particles")`  
`Pkg.activate(".")`  
`using Particles`  
`cd("examples")`  
`include("example_loop.jl")`  
`run_simulation(d)`  

# package development on linux
- activate Particles
cd to folder on bash prompt e.g. `cd $HOME/src/julia/Particles`
start julia with `julia --project`
or activate the julia package management environment with the ']'-key on the REPL and use commands like "activate ."
- update packages
If you expect missing packages on you machine then run  
`Pkg.update()`
- running unit tests
The unit tests reside in tests. The main script is runtests.jl. To run them type:  
`Pkg.test("Particles")`  
## example of using the package
`cd $HOME/src/julia/Particles`  
`julia --project`  
`using Particles`  
`cd("examples")`  
`include("example_loop.jl")`  
`run_simulation(d)`  

