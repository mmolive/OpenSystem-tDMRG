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
#include "open.h"


#ifndef XXZ_H
#define XXZ_H


class XXZ_open: public open {
private:
    int B;
    double J, delta, h, gamma, Gamma, mu, *monitor;
    ITensor Id, Sx, Sy, Sz, Su, Sd, Szm, SdSz;
    ofstream outM, outJ, out_corr, out_Jsus, out_cov, outMonitor, outSpectrum;

public:
    XXZ_open (int, int, int, int, double, double, double, double, double, double, int, int, int, double, double);
    ~XXZ_open ();

    void init_all_down ();

    uint64_t ID_XXZ (int, double, double);

    void time_evolution_operator (double, int);

    void time_evolution_operator_4th (double);

    void Lindbladian_MPO ();

    void Measure (int, int, double);

    void Monitor (int);

    void Measure_final ();

    void Singular_values ();
};


#endif
