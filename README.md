# OpenSystem-tDMRG
In this repository we share the program used to generate the data for the paper entitled "Efficient quantum information probes of non-equilibrium quantum criticality". The paper can be found at https://doi.org/10.1038/s41534-022-00671-8

This program uses Matrix Product State techniques to obtain the steady-states of Markovian open quantum systems. The MPS implementation is based on the ITensor library.


# Installing Package

- Clone this repository
- Install the ITensor directory and change point 2 of the Makefile to link to its location
- Run on the terminal `./setup.sh`. This will download the Eigen library and create the directories used by this package.


# Bechmark

Here some instruction to run the code are provided

This program was implemented specifically for a fermionic chain and a fermionic ladder model. Alter the Makefile (point 3 and 5) to chose between the two.


## Running code for chain

To simulate the fermionic chain model run on the terminal `./XXZ_sim L D delta mu Nt dt step restore D_old delta_old mu_old`, where we have the following variables
- `L` is the number of sites in the chain
- `D` is the bond dimension
- `delta` is a rescaled interaction, such that the fermionic repulsive interaction in the paper is given by `V = 4 delta`
- `mu` is the reservoir's bias
- `Nt` is the number of time evolution steps
- `dt` is the time-step
- `step` is the number of time-steps between each measurement
- `restore` if `0` the simulation starts from an initial infinite temperature state, if `1` it continues from a previously obtained state
- `D_old`, `delta_old` and `mu_old` are the parameters of this previous state

The output files of a simulation are labeled by a variable `ID = Hash(Physical Parameters)`. To help the user read the data we included a Mathematica file `analysis_chain.nb`. There is an example in it that can be replicated by the user with the following inputs `./XXZ_sim 4 16 0.25 0.1 10000 0.05 50 0 16 0.25 0.1`


## Running code for ladder

To simulate the fermionic ladder model run on the terminal `./FermionLadder_sim L D V mu Nt dt step restore D_old mode`, where we have the following variables
- `L` is the number of rungs in the ladder
- `D` is the bond dimension
- `V` is the repulsive interaction
- `mu` is the reservoir's bias
- `Nt` is the number of time evolution steps
- `dt` is the time-step
- `step` is the number of time-steps between each measurement
- `restore` if `0` the simulation starts from an initial infinite temperature state, if `1` it continues from a previously obtained state
- `D_old` is the parameter of this previous state
- `mode` is variable that controls the reservoir's configuration. The one used in the paper corresponds to `0`

There is an example in the file `analysis_ladder.nb` that can be replicated by the user with the following inputs `./FermionLadder_sim 3 16 1.5 0.1 10000 0.05 50 0 16 0`
