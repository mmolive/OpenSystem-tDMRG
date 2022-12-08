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
#include "FermionLadder_open.h"


int main(int argc, char* argv[]){

    int L, d, D, B_L, Nt, step, restore, D_res, mode;
    double tl, tr, V, Vl, Vr, mu, muin, gamma, dt;
    double muL1, muL2, muR1, muR2, GammaL1, GammaL2, GammaR1, GammaR2;

    //Physical parameters
    sscanf(argv[1],"%d",&L);
    sscanf(argv[2],"%d",&D);
    sscanf(argv[3],"%lf",&V);
    sscanf(argv[4],"%lf",&muin);
    sscanf(argv[5],"%d",&Nt);
    sscanf(argv[6],"%lf",&dt);
    sscanf(argv[7],"%d",&step);
    sscanf(argv[8],"%d",&restore);
    sscanf(argv[9],"%d",&D_res);
    sscanf(argv[10],"%d",&mode);

    d=2;
    B_L= 14;
    tl= -1.0;
    tr= -1.0;
    Vl= V;
    Vr= V;
    mu= 0.0;
    gamma= 0.0;
    GammaL1= 1.0;
    GammaL2= 1.0;
    GammaR1= 1.0;
    GammaR2= 1.0;

    if(mode==0){
        muL1= muin;
        muL2= muin;
        muR1= muin;
        muR2= muin;
    }
    if(mode==1){
        muL1= muin;
        muL2= -muin;
        muR1= muin;
        muR2= -muin;
    }
    if(mode==2){
        muL1= muin;
        muL2= -muin;
        muR1= -muin;
        muR2= muin;
    }


    //Initializations
    Ladder_open *Ladder = new Ladder_open(L, d*d*d*d, D, B_L, tl, tr, Vl, Vr, mu, gamma, muL1, muL2, muR1, muR2, GammaL1, GammaL2, GammaR1, GammaR2, Nt/step/10, restore, D_res);
    //Ladder->time_evolution_operator(dt, 2);
    Ladder->time_evolution_operator_4th(dt);
    //Ladder->init_domain();


    //Simulation
    Ladder->tDMRG_4th_order(Nt, step, dt);

    //Ladder->tDMRG_2nd_order(Nt, step, dt);



    delete Ladder;

    return 0;
}
