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
#include <iostream>
#include <fstream>
//#include <mkl.h>
#include "itensor/all.h"
#include "../Eigen/Eigen/Eigen"
#include "ManyBodyOp.h"

using namespace std;
using namespace itensor;


#ifndef MPS_H
#define MPS_H


class myMPS {
public:
    int L, d, D;
    uint64_t ID;
    double Time;
    MPS M;
    ITensor *U, *U1, *U2;

    myMPS (int, int, int);
    ~myMPS ();

    template <class T>
    uint64_t FNVHash(T a, uint64_t hash){

        uint64_t prime;
        prime = 1099511628211;

        unsigned char* data = reinterpret_cast<unsigned char*>(&a);
        for(uint64_t i=0; i<sizeof(T); ++i)
            hash = (data[i] ^ hash) * prime;

        return hash;
    }

    void tDMRG_1st_order(int, int, double, double&);

    void tDMRG_2nd_order(int, int, double);

    void tDMRG_4th_order(int, int, double);

    void tDMRG_step(int);

    void tDMRG_step(int, ITensor *);

    void PrintTime();

    int CheckStop();

    //void tDMRG_step_old(int);

    virtual void Measure (int, int, double) =0;

    virtual void Monitor (int) =0;

    virtual void Measure_final () =0;

    virtual void Singular_values () =0;
};



#endif
