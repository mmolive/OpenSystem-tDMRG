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

    int L, d, D, B_L;
    double J, h, delta, gamma, Gamma, mu;

    //Physical parameters
    sscanf(argv[1],"%d",&L);
    sscanf(argv[2],"%d",&D);
    sscanf(argv[3],"%lf",&delta);
    sscanf(argv[4],"%lf",&mu);


    d=2;
    B_L= 8;
    J= -0.5;
    h= 0.0;
    gamma= 0.0;
    Gamma= 1.0;


    XXZ_open *XXZ = new XXZ_open(L, d*d, D, B_L, J, h, delta, gamma, Gamma, mu, 10, 1, D, delta, mu);
    XXZ->Measure_final();
    XXZ->Singular_values();


    delete XXZ;

    return 0;
}
