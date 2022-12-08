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


#ifndef LADDER_H
#define LADDER_H


class Ladder_open: public open {
private:
    int B;
    double tl, tr, Vl, Vr, mu, gamma, *monitor;
    double muL1, muL2, muR1, muR2, GammaL1, GammaL2, GammaR1, GammaR2;
    ITensor Id, Sx_up_even, Sx_up_odd, Sy_up_even, Sy_up_odd, Sx_down_even, Sx_down_odd, Sy_down_even, Sy_down_odd, nf_up, nf_down, J_rung;
    ITensor Su_up_even, Su_up_odd, Sd_up_even, Sd_up_odd, Su_down_even, Su_down_odd, Sd_down_even, Sd_down_odd, Szm, ST;
    ofstream outM_up, outM_down, out_corr_up, out_corr_down, outJ_up, outJ_down, outJ_rung, out_Jsus, out_cov, out_covF, outMonitor, outSpectrum;

public:
    Ladder_open (int, int, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, int, int, int);
    ~Ladder_open ();

    void init_all_down ();

    void init_infT ();

    void init_domain ();

    void init_PreState (int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

    uint64_t ID_SpinLadder (int);

    uint64_t ID_SpinLadder (int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double);

    void time_evolution_operator (double, int);

    void time_evolution_operator_4th (double);

    void Lindbladian_MPO ();

    void Measure (int, int, double);

    void Monitor (int);

    void Measure_final ();

    void Singular_values ();
};


#endif
