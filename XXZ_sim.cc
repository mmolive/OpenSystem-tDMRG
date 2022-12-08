//
//Copyright 2022 Miguel M. Oliveira
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.
//
#include "XXZ_open.h"


int main(int argc, char* argv[]){

    int L, d, D, B_L, Nt, step, restore, D_old;
    double J, h, delta, delta_old, gamma, Gamma, mu, mu_old, dt;

    //Physical parameters
    sscanf(argv[1],"%d",&L);
    sscanf(argv[2],"%d",&D);
    sscanf(argv[3],"%lf",&delta);
    sscanf(argv[4],"%lf",&mu);
    sscanf(argv[5],"%d",&Nt);
    sscanf(argv[6],"%lf",&dt);
    sscanf(argv[7],"%d",&step);
    sscanf(argv[8],"%d",&restore);
    sscanf(argv[9],"%d",&D_old);
    sscanf(argv[10],"%lf",&delta_old);
    sscanf(argv[11],"%lf",&mu_old);


    d=2;
    B_L= 8;
    J= -0.5;
    h= 0.0;
    gamma= 0.0;
    Gamma= 1.0;


    //Initializations
    XXZ_open *XXZ = new XXZ_open(L, d*d, D, B_L, J, h, delta, gamma, Gamma, mu, Nt/step/10, restore, D_old, delta_old, mu_old);
    //XXZ->Lindbladian_MPO();
    //XXZ->time_evolution_operator(dt, 2);
    XXZ->time_evolution_operator_4th(dt);


    //Simulation
    XXZ->tDMRG_4th_order(Nt, step, dt);
    //XXZ->tDMRG_2nd_order(Nt, step, dt);


    delete XXZ;

    return 0;
}
