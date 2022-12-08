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
#include "mps.h"

clock_t t_ref, t_position= 0, t_measure= 0, t_normalize= 0, t_svd= 0, t_full;


myMPS::myMPS (int Lin, int din, int Din) {

    L = Lin;
    d = din;
    D = Din;

    auto sites = SiteSet(L, d);
    M = MPS(sites);

    U = new ITensor [L-1];
    U1 = new ITensor [L-1];
    U2 = new ITensor [L-1];

    for(int i=0; i<L-1; ++i){
        U[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
        U1[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
        U2[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
    }
}

myMPS::~myMPS(){
    delete[] U;
    delete[] U1;
    delete[] U2;
}


void myMPS::tDMRG_1st_order(int Nt, int step, double dt, double &Time){

    char name[100];
    ofstream out;
    sprintf(name, "time/XXZ_comp_%d_%d.txt", L, D);
    out.open(name, std::ios_base::app);
    sprintf(name, "output/open/Final_state_%lu.bin", ID);

    int l, t;

    M.orthogonalize({"MaxDim",D});

    for(t=1; t<=Nt; ++t){

        if(t%(Nt/10)==0)
            cout << Time + t*dt << endl;

        //apply time evolution operator to odd bonds
        for(l=1; l<L; l+=2)
            tDMRG_step(l);

        //write state in right canonical form
        t_ref = time(NULL);
        M.position(1);
        t_position += time(NULL) - t_ref;

        t_ref = time(NULL);
        M.normalize();
        t_normalize += time(NULL) - t_ref;

        //apply time evolution operator to even bonds
        for(l=2; l<L; l+=2)
            tDMRG_step(l);

        //write state in right canonical form again
        t_ref = time(NULL);
        M.position(1);
        t_position += time(NULL) - t_ref;

        t_ref = time(NULL);
        M.normalize();
        t_normalize += time(NULL) - t_ref;

        //measurements
        if(t%step==0){
            t_ref = time(NULL);
            Measure(t/step-1, Nt/step, Time + t*dt);
            t_measure += time(NULL) - t_ref;
        }

        if(t%(Nt/10)==0)
            writeToFile(name, M);
    }
    Time = Time + Nt*dt;

    Measure_final();
    Singular_values();


    out << t_svd << " " << t_position << " " << t_measure << " " << t_normalize << endl;
    out.close();
}


void myMPS::tDMRG_2nd_order(int Nt, int step, double dt){

    int l, t, Time0;
    char name[100];
    ofstream out, tTotal;
    sprintf(name, "time/ltime_%d_%d_%lu.txt", L, D, ID);
    out.open(name, std::ios_base::app);
    sprintf(name, "time/Ttime_%d_%d_%lu.txt", L, D, ID);
    tTotal.open(name, std::ios_base::app);
    sprintf(name, "output/open/Final_state_%lu.bin", ID);

    M.orthogonalize({"MaxDim",D});

    Time0 = Time;
    t_full = clock();
    for(t=1; t<=Nt; ++t){

        if(t%(Nt/10)==0)
            cout << Time0 + t*dt << endl;

        //apply time evolution operator to even bonds
        for(l=2; l<L; l+=2)
            tDMRG_step(l);

        //write state in right canonical form again
        /*t_ref = clock();
        M.position(1);
        t_position += clock() - t_ref;*/


        //apply time evolution operator to odd bonds
        for(l=1; l<L; l+=2)
            tDMRG_step(l);

        //write state in right canonical form
        /*t_ref = clock();
        M.position(1);
        t_position += clock() - t_ref;*/


        //apply time evolution operator to even bonds again
        for(l=2; l<L; l+=2)
            tDMRG_step(l);

        //Normalize
        t_ref = clock();
        M.position(1);
        t_position += clock() - t_ref;

        t_ref = clock();
        M.normalize();
        t_normalize += clock() - t_ref;


        //measurements
        if(t%step==0){
            t_ref = clock();
            Measure(t/step-1, Nt/step/10, Time0 + t*dt);
            t_measure += clock() - t_ref;

            if(CheckStop()==1)
                break;
        }

        if(t%(Nt/10)==0){
            writeToFile(name, M);
            Monitor(Nt/step/10);
            Time = Time + Nt/10*dt;
            PrintTime();
        }
    }
    //Time = Time + (t-1)*dt;
    t_full = clock() - t_full;


    Measure_final();
    Singular_values();

    tTotal << (double)t_full/CLOCKS_PER_SEC << endl;
    out << (double)t_svd/CLOCKS_PER_SEC << " " << (double)t_position/CLOCKS_PER_SEC << " " << (double)t_measure/CLOCKS_PER_SEC << " " << (double)t_normalize/CLOCKS_PER_SEC << endl;
    out.close();
    tTotal.close();
}


void myMPS::tDMRG_4th_order(int Nt, int step, double dt){

    int l, t, Time0;
    char name[100];
    sprintf(name, "output/open/Final_state_%lu.bin", ID);

    M.orthogonalize({"MaxDim",D});


    Time0 = Time;
    for(t=1; t<=Nt; ++t){

        if(t%(Nt/10)==0)
            cout << Time0 + t*dt << endl;


        //apply time evolution operator to odd bonds 1st time
        for(l=1; l<L; l+=2)
            tDMRG_step(l, U);

        //apply time evolution operator to even bonds 1st time
        for(l=2; l<L; l+=2)
            tDMRG_step(l, U);

        //apply time evolution operator to odd bonds 2nd time
        for(l=1; l<L; l+=2)
            tDMRG_step(l, U1);

        //apply time evolution operator to even bonds 2nd time
        for(l=2; l<L; l+=2)
            tDMRG_step(l, U1);

        //apply time evolution operator to odd bonds 3rd time
        for(l=1; l<L; l+=2)
            tDMRG_step(l, U2);

        //renormalize
        M.position(1);
        M.normalize();

        //measurements
        if(t%step==0){
            Measure(t/step-1, Nt/step/10, Time0 + t*dt);
            if(CheckStop()==1)
                break;
        }

        if(t%(Nt/10)==0){
            writeToFile(name, M);
            Monitor(Nt/step/10);
            Time = Time + Nt/10*dt;
            PrintTime();
        }
    }

    Measure_final();
    Singular_values();
}


void myMPS::tDMRG_step(int loc){

    double err;

    t_ref = clock();
    M.position(loc);
    t_position += clock() - t_ref;

    t_ref = clock();

    auto phi = M(loc) * M(loc+1);
    phi *= U[loc-1];
    phi.noPrime();

    auto [u,S,V] = svd(phi, inds(M(loc)),{"MaxDim",D});
    M.set(loc, u);
    M.set(loc+1, S*V);

    t_svd += clock() - t_ref;

    //err = norm(u*S*V-phi)/norm(phi);
    //error += err*err;
}


void myMPS::tDMRG_step(int loc, ITensor *U){

    double err;

    M.position(loc);
    //M.normalize();

    auto phi = M(loc) * M(loc+1);
    phi *= U[loc-1];
    phi.noPrime();

    auto [u,S,V] = svd(phi, inds(M(loc)),{"MaxDim",D});
    M.set(loc, u);
    M.set(loc+1, S*V);
}


int myMPS::CheckStop(){

    char name[100];
    ifstream stop;
    sprintf(name, "output/stop/stop_%lu.txt", ID);
    stop.open(name);

    return stop.is_open();
}


void myMPS::PrintTime(){

    char name[100];
    ofstream outTime;
    sprintf(name, "output/open/Time_%lu.txt", ID);
    outTime.open(name);
    outTime << Time << endl;
    outTime.close();
}


/*
void myMPS::tDMRG_step_alt(int loc){

    t_ref = time(NULL);

    auto phi = M(loc) * M(loc+1);
    phi *= U[loc-1];
    phi.noPrime();

    auto [u,S,V] = svd(phi, inds(M(loc)),{"MaxDim",D});
    M.set(loc, u);
    M.set(loc+1, S*V);

    t_svd += time(NULL) - t_ref;


    if(loc+1 != L && loc+1 != L-1){

        t_ref = time(NULL);

        phi = M(loc+1) * M(loc+2);

        auto [u,S,V] = svd(phi, inds(M(loc+1)),{"MaxDim",D});
        M.set(loc+1, u);
        M.set(loc+2, S*V);

        t_position += time(NULL) - t_ref;
    }
}*/
