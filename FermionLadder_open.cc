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


Ladder_open::Ladder_open(int Lin, int din, int Din, int Bin, double tlin, double trin, double Vlin, double Vrin, double muin, double gammain, double muL1in, double muL2in, double muR1in, double muR2in, double GammaL1in, double GammaL2in, double GammaR1in, double GammaR2in, int Nt, int restore, int D_res) : open(Lin, din, Din) {

    B = Bin;
    tl= tlin;
    tr= trin;
    Vl= Vlin;
    Vr= Vrin;
    mu= muin;
    gamma= gammain;
    muL1= muL1in;
    muL2= muL2in;
    muR1= muR1in;
    muR2= muR2in;
    GammaL1= GammaL1in;
    GammaL2= GammaL2in;
    GammaR1= GammaR1in;
    GammaR2= GammaR2in;
    monitor = new double[3*Nt];

    auto i = Index(d);
    Id = ITensor(i);

    Sx_up_even = ITensor(i);
    Sx_up_odd = ITensor(i);
    Sy_up_even = ITensor(i);
    Sy_up_odd = ITensor(i);
    Sx_down_even = ITensor(i);
    Sx_down_odd = ITensor(i);
    Sy_down_even = ITensor(i);
    Sy_down_odd = ITensor(i);
    Su_up_even = ITensor(i);
    Su_up_odd = ITensor(i);
    Sd_up_even = ITensor(i);
    Sd_up_odd = ITensor(i);
    Su_down_even = ITensor(i);
    Su_down_odd = ITensor(i);
    Sd_down_even = ITensor(i);
    Sd_down_odd = ITensor(i);
    Szm = ITensor(i);
    ST = ITensor(i);
    nf_up = ITensor(i);
    nf_down = ITensor(i);
    J_rung = ITensor(i);
    for(int i=1; i<=d; ++i){
        Id.set(i, 0.0);
        Sx_up_even.set(i, 0.0);
        Sx_up_odd.set(i, 0.0);
        Sy_up_even.set(i, 0.0);
        Sy_up_odd.set(i, 0.0);
        Sx_down_even.set(i, 0.0);
        Sx_down_odd.set(i, 0.0);
        Sy_down_even.set(i, 0.0);
        Sy_down_odd.set(i, 0.0);
        Su_up_even.set(i, 0.0);
        Su_up_odd.set(i, 0.0);
        Sd_up_even.set(i, 0.0);
        Sd_up_odd.set(i, 0.0);
        Su_down_even.set(i, 0.0);
        Su_down_odd.set(i, 0.0);
        Sd_down_even.set(i, 0.0);
        Sd_down_odd.set(i, 0.0);
        Szm.set(i, 0.0);
        ST.set(i, 0.0);
        nf_up.set(i, 0.0);
        nf_down.set(i, 0.0);
        J_rung.set(i, 0.0);
    }
    Id.set(1, 1.0); Id.set(6, 1.0); Id.set(11, 1.0); Id.set(16, 1.0);
    Sx_up_even.set(3, 1.0); Sx_up_even.set(8, 1.0); Sx_up_even.set(9, 1.0); Sx_up_even.set(14, 1.0);
    Sx_up_odd.set(3, 1.0); Sx_up_odd.set(8, -1.0); Sx_up_odd.set(9, 1.0); Sx_up_odd.set(14, -1.0);
    Sy_up_even.set(3, 1.0_i); Sy_up_even.set(8, 1.0_i); Sy_up_even.set(9, -1.0_i); Sy_up_even.set(14, -1.0_i);
    Sy_up_odd.set(3, 1.0_i); Sy_up_odd.set(8, -1.0_i); Sy_up_odd.set(9, -1.0_i); Sy_up_odd.set(14, 1.0_i);
    Sx_down_even.set(2, 1.0); Sx_down_even.set(5, 1.0); Sx_down_even.set(12, -1.0); Sx_down_even.set(15, -1.0);
    Sx_down_odd.set(2, 1.0); Sx_down_odd.set(5, 1.0); Sx_down_odd.set(12, 1.0); Sx_down_odd.set(15, 1.0);
    Sy_down_even.set(2, 1.0_i); Sy_down_even.set(5, -1.0_i); Sy_down_even.set(12, -1.0_i); Sy_down_even.set(15, 1.0_i);
    Sy_down_odd.set(2, 1.0_i); Sy_down_odd.set(5, -1.0_i); Sy_down_odd.set(12, 1.0_i); Sy_down_odd.set(15, -1.0_i);
    Szm.set(1, 1.0); Szm.set(6, -1.0); Szm.set(11, -1.0); Szm.set(16, 1.0);
    Su_up_even.set(9, -1.0); Su_up_even.set(14, 1.0);
    Su_up_odd.set(9, 1.0); Su_up_odd.set(14, 1.0);
    Sd_up_even.set(3, 1.0); Sd_up_even.set(8, 1.0);
    Sd_up_odd.set(3, -1.0); Sd_up_odd.set(8, 1.0);
    Su_down_even.set(5, 1.0); Su_down_even.set(15, 1.0);
    Su_down_odd.set(5, -1.0); Su_down_odd.set(15, 1.0);
    Sd_down_even.set(2, -1.0); Sd_down_even.set(12, 1.0);
    Sd_down_odd.set(2, 1.0); Sd_down_odd.set(12, 1.0);
    nf_up.set(1, 1.0); nf_up.set(6, 1.0);
    nf_down.set(1, 1.0); nf_down.set(11, 1.0);
    J_rung.set(7, -1.0_i); J_rung.set(10, 1.0_i);
    ST.set(7, 1.0);

    ID = ID_SpinLadder(D);
    char name[100];

    //input older state
    if(restore==1){

        uint64_t ID_res = ID_SpinLadder(D_res);
        sprintf(name, "output/open/Final_state_%lu.bin", ID_res);
        readFromFile(name, M);

        if(ID==ID_res){
            ifstream inTime;
            sprintf(name, "output/open/Time_%lu.txt", ID);
            inTime.open(name);
            inTime >> Time;
            inTime.close();
        }
        else
            Time=0;

        for(int i=0; i<L-1; ++i){
            U[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
            U1[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
            U2[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
        }
    }
    else{
        Time=0;
        init_infT();
        sprintf(name, "output/open/Monitor_%lu.bin", ID);
        outMonitor.open(name, ios::out | ios::binary);
        outMonitor.close();
    }


    Lindblad = MPO(siteInds(M));
    auto sites = siteInds(Lindblad,1);
    auto links = IndexSet(Index(B,"l=1,Link"), Index(B,"l=2,Link"));
    for(int l=2; l<L-1; ++l){
        sprintf(name, "l=%d,Link", l+1);
        links = unionInds(links, Index(B,name));
        sites = unionInds(sites, siteInds(Lindblad,l));
    }
    sites = unionInds(sites, siteInds(Lindblad,L-1));
    sites = unionInds(sites, siteInds(Lindblad,L));

    Lindblad.set(1, ITensor( sites(1), sites(2), links(1)));
    for(int l=2; l<L; ++l)
        Lindblad.set(l, ITensor( sites(2*l-1), sites(2*l), links(l-1), links(l)));
    Lindblad.set(L, ITensor( sites(2*L-1), sites(2*L), links(L-1)));
}

Ladder_open::~Ladder_open (){
    delete[] monitor;
}


//Sets initial state with all spins pointing down
void Ladder_open::init_all_down (){

    for(int j=0; j<d-1; ++j)
        M.ref(1).set(j+1,1, 0.0);
    M.ref(1).set(d,1, 1.0);

    for(int i=2; i<L; ++i){
        for(int j=0; j<d-1; ++j)
            M.ref(i).set(1,j+1,1, 0.0);
        M.ref(i).set(1,d,1, 1.0);
    }

    for(int j=0; j<d-1; ++j)
        M.ref(L).set(1,j+1, 0.0);
    M.ref(L).set(1,d, 1.0);
}

//Sets initial infinite temperature state
void Ladder_open::init_infT (){

    for(int j=0; j<d; ++j)
        M.ref(1).set(j+1,1, 1.0);

    for(int i=2; i<L; ++i){
        for(int j=0; j<d; ++j)
            M.ref(i).set(1,j+1,1, 1.0);
    }

    for(int j=0; j<d; ++j)
        M.ref(L).set(1,j+1, 1.0);
}

//Sets initial state with domain wall
void Ladder_open::init_domain (){

    for(int j=1; j<d; ++j)
        M.ref(1).set(j+1,1, 0.0);
    M.ref(1).set(1,1, 1.0);

    for(int i=2; i<=L/2; ++i){
        for(int j=1; j<d; ++j)
            M.ref(i).set(1,j+1,1, 0.0);
        M.ref(i).set(1,1,1, 1.0);
    }

    for(int i=L/2+1; i<L; ++i){
        for(int j=0; j<d-1; ++j)
            M.ref(i).set(1,j+1,1, 0.0);
        M.ref(i).set(1,d,1, 1.0);
    }

    for(int j=0; j<d-1; ++j)
        M.ref(L).set(1,j+1, 0.0);
    M.ref(L).set(1,d, 1.0);
}

void Ladder_open::init_PreState (int Lin, int Din, double tlin, double trin, double Vlin, double Vrin, double muin, double gammain, double muL1in, double muL2in, double muR1in, double muR2in, double GammaL1in, double GammaL2in, double GammaR1in, double GammaR2in){

  char name[100];
  uint64_t ID_res = ID_SpinLadder(Lin, Din, tlin, trin, Vlin, Vrin, muin, gammain, muL1in, muL2in, muR1in, muR2in, GammaL1in, GammaL2in, GammaR1in, GammaR2in);

  sprintf(name, "output/open/Final_state_%lu.bin", ID_res);
  readFromFile(name, M);

  for(int i=0; i<L-1; ++i){
    U[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
    U1[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
    U2[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
  }
}


//gets unique ID computed from simulation parameters
uint64_t Ladder_open::ID_SpinLadder (int Din){

    uint64_t hash;
    hash = 14695981039346656037;

    hash = FNVHash<int>(L, hash);
    hash = FNVHash<int>(Din, hash);
    hash = FNVHash<double>(tl, hash);
    hash = FNVHash<double>(tr, hash);
    hash = FNVHash<double>(Vl, hash);
    hash = FNVHash<double>(Vr, hash);
    hash = FNVHash<double>(mu, hash);
    hash = FNVHash<double>(gamma, hash);
    hash = FNVHash<double>(muL1, hash);
    hash = FNVHash<double>(muL2, hash);
    hash = FNVHash<double>(muR1, hash);
    hash = FNVHash<double>(muR2, hash);
    hash = FNVHash<double>(GammaL1, hash);
    hash = FNVHash<double>(GammaL2, hash);
    hash = FNVHash<double>(GammaR1, hash);
    hash = FNVHash<double>(GammaR2, hash);

    return hash;
}

uint64_t Ladder_open::ID_SpinLadder (int Lin, int Din, double tlin, double trin, double Vlin, double Vrin, double muin, double gammain, double muL1in, double muL2in, double muR1in, double muR2in, double GammaL1in, double GammaL2in, double GammaR1in, double GammaR2in){

    uint64_t hash;
    hash = 14695981039346656037;

    hash = FNVHash<int>(Lin, hash);
    hash = FNVHash<int>(Din, hash);
    hash = FNVHash<double>(tlin, hash);
    hash = FNVHash<double>(trin, hash);
    hash = FNVHash<double>(Vlin, hash);
    hash = FNVHash<double>(Vrin, hash);
    hash = FNVHash<double>(muin, hash);
    hash = FNVHash<double>(gammain, hash);
    hash = FNVHash<double>(muL1in, hash);
    hash = FNVHash<double>(muL2in, hash);
    hash = FNVHash<double>(muR1in, hash);
    hash = FNVHash<double>(muR2in, hash);
    hash = FNVHash<double>(GammaL1in, hash);
    hash = FNVHash<double>(GammaL2in, hash);
    hash = FNVHash<double>(GammaR1in, hash);
    hash = FNVHash<double>(GammaR2in, hash);

    return hash;
}


//Function that writes the time-evolution operators
void Ladder_open::time_evolution_operator(double dt, int order){

    using namespace Eigen;

    int dim=2;
    MatrixXcd LindBulk, LindLeft, LindRight, expo;
    MatrixXcd Sx(dim,dim), Sy(dim,dim), Sz(dim,dim), Su(dim,dim), Sd(dim,dim), nf(dim,dim), Id(dim,dim);
    Pauli_Matrices(dim, Sx, Sy, Sz, Su, Sd, Id);
    nf= Sz/2;


    //Bulk even bonds
    LindBulk = Hamiltonian_bond(dim, 2, 0, 2, tl, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 2, tl, Sd, Su, Id);
    LindBulk += Hamiltonian_fermion_bond(2, 1, 3, tl, Su, Sd, Id, Sz, 1);
    LindBulk += Hamiltonian_fermion_bond(2, 1, 3, tl, Sd, Su, Id, Sz, 1);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Sd, Su, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 2, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 1, 3, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, Vr, nf, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 0, mu, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 1, mu, nf, Id);
    LindBulk += Jump_site(dim, 2, 0, gamma, nf, Id);
    LindBulk += Jump_site(dim, 2, 1, gamma, nf, Id);
    if(L%2==1)
        LindRight = LindBulk;
    if(order==2)
        LindBulk = LindBulk/2;

    expo = matrix_complex_exponential(d*d, dt, LindBulk);
    for(int i=1; i<L-2; i+=2)
        TEvol_Eigen_itensor(d, expo, &U[i]);


    //Bulk odd bonds
    LindBulk = Hamiltonian_bond(dim, 2, 1, 3, tl, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 1, 3, tl, Sd, Su, Id);
    LindBulk += Hamiltonian_fermion_bond(2, 0, 2, tl, Su, Sd, Id, Sz, 1);
    LindBulk += Hamiltonian_fermion_bond(2, 0, 2, tl, Sd, Su, Id, Sz, 1);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Sd, Su, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 2, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 1, 3, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, Vr, nf, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 0, mu, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 1, mu, nf, Id);
    LindBulk += Jump_site(dim, 2, 0, gamma, nf, Id);
    LindBulk += Jump_site(dim, 2, 1, gamma, nf, Id);
    LindLeft = LindBulk;
    if(L%2==0)
        LindRight = LindBulk;

    expo = matrix_complex_exponential(d*d, dt, LindBulk);
    for(int i=2; i<L-2; i+=2)
        TEvol_Eigen_itensor(d, expo, &U[i]);


    //Left bond
    LindLeft += Jump_site(dim, 2, 0, 0.5*GammaL1*(1+muL1), Su, Id);
    LindLeft += Jump_site(dim, 2, 0, 0.5*GammaL1*(1-muL1), Sd, Id);
    LindLeft += Jump_fermion_site(2, 1, 0, 0.5*GammaL2*(1+muL2), Su, Id, Sz);
    LindLeft += Jump_fermion_site(2, 1, 0, 0.5*GammaL2*(1-muL2), Sd, Id, Sz);

    expo = matrix_complex_exponential(d*d, dt, LindLeft);
    TEvol_Eigen_itensor(d, expo, &U[0]);


    //Right bond
    LindRight += Hamiltonian_bond(dim, 2, 2, 3, tr, Su, Sd, Id);
    LindRight += Hamiltonian_bond(dim, 2, 2, 3, tr, Sd, Su, Id);
    LindRight += Hamiltonian_bond(dim, 2, 2, 3, Vr, nf, nf, Id);
    LindRight += Hamiltonian_site(dim, 2, 2, mu, nf, Id);
    LindRight += Hamiltonian_site(dim, 2, 3, mu, nf, Id);
    LindRight += Jump_site(dim, 2, 2, gamma, nf, Id);
    LindRight += Jump_site(dim, 2, 3, gamma, nf, Id);
    if(L%2==0){
        LindRight += Jump_site(dim, 2, 2, 0.5*GammaR1*(1-muR1), Su, Id);
        LindRight += Jump_site(dim, 2, 2, 0.5*GammaR1*(1+muR1), Sd, Id);

        LindRight += Jump_fermion_site(2, 3, 2, 0.5*GammaR2*(1-muR2), Su, Id, Sz);
        LindRight += Jump_fermion_site(2, 3, 2, 0.5*GammaR2*(1+muR2), Sd, Id, Sz);
    }
    else{
        LindRight += Jump_fermion_site(2, 2, 3, 0.5*GammaR1*(1-muR1), Su, Id, Sz);
        LindRight += Jump_fermion_site(2, 2, 3, 0.5*GammaR1*(1+muR1), Sd, Id, Sz);

        LindRight += Jump_site(dim, 2, 3, 0.5*GammaR2*(1-muR2), Su, Id);
        LindRight += Jump_site(dim, 2, 3, 0.5*GammaR2*(1+muR2), Sd, Id);
        if(order==2)
            LindRight = LindRight/2;
    }

    expo = matrix_complex_exponential(d*d, dt, LindRight);
    TEvol_Eigen_itensor(d, expo, &U[L-2]);
}


//Function that writes the time-evolution operators in 4th order
void Ladder_open::time_evolution_operator_4th(double dt){

    using namespace Eigen;

    int dim=2;
    MatrixXcd LindBulk, LindLeft, LindRight, expo;
    MatrixXcd Sx(dim,dim), Sy(dim,dim), Sz(dim,dim), Su(dim,dim), Sd(dim,dim), nf(dim,dim), Id(dim,dim);
    Pauli_Matrices(dim, Sx, Sy, Sz, Su, Sd, Id);
    nf= Sz/2;


    //Bulk even bonds
    LindBulk = Hamiltonian_bond(dim, 2, 0, 2, tl, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 2, tl, Sd, Su, Id);
    LindBulk += Hamiltonian_fermion_bond(2, 1, 3, tl, Su, Sd, Id, Sz, 1);
    LindBulk += Hamiltonian_fermion_bond(2, 1, 3, tl, Sd, Su, Id, Sz, 1);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Sd, Su, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 2, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 1, 3, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, Vr, nf, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 0, mu, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 1, mu, nf, Id);
    LindBulk += Jump_site(dim, 2, 0, gamma, nf, Id);
    LindBulk += Jump_site(dim, 2, 1, gamma, nf, Id);
    if(L%2==1)
        LindRight = LindBulk;

    expo = LindBulk * complex<double>(0.5, -sqrt(3)/6.0);
    expo = matrix_complex_exponential(d*d, dt, expo);
    for(int i=1; i<L-2; i+=2)
        TEvol_Eigen_itensor(d, expo, &U[i]);

    expo = LindBulk * complex<double>(0.5, sqrt(3)/6.0);
    expo = matrix_complex_exponential(d*d, dt, expo);
    for(int i=1; i<L-2; i+=2)
        TEvol_Eigen_itensor(d, expo, &U1[i]);


    //Bulk odd bonds
    LindBulk = Hamiltonian_bond(dim, 2, 1, 3, tl, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 1, 3, tl, Sd, Su, Id);
    LindBulk += Hamiltonian_fermion_bond(2, 0, 2, tl, Su, Sd, Id, Sz, 1);
    LindBulk += Hamiltonian_fermion_bond(2, 0, 2, tl, Sd, Su, Id, Sz, 1);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Su, Sd, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, tr, Sd, Su, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 2, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 1, 3, Vl, nf, nf, Id);
    LindBulk += Hamiltonian_bond(dim, 2, 0, 1, Vr, nf, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 0, mu, nf, Id);
    LindBulk += Hamiltonian_site(dim, 2, 1, mu, nf, Id);
    LindBulk += Jump_site(dim, 2, 0, gamma, nf, Id);
    LindBulk += Jump_site(dim, 2, 1, gamma, nf, Id);
    LindLeft = LindBulk;
    if(L%2==0)
        LindRight = LindBulk;

    expo = LindBulk * complex<double>(0.25, -sqrt(3)/12.0);
    expo = matrix_complex_exponential(d*d, dt, expo);
    for(int i=2; i<L-2; i+=2)
        TEvol_Eigen_itensor(d, expo, &U[i]);

    expo = LindBulk / 2;
    expo = matrix_complex_exponential(d*d, dt, expo);
    for(int i=2; i<L-2; i+=2)
        TEvol_Eigen_itensor(d, expo, &U1[i]);

    expo = LindBulk * complex<double>(0.25, sqrt(3)/12.0);
    expo = matrix_complex_exponential(d*d, dt, expo);
    for(int i=2; i<L-2; i+=2)
        TEvol_Eigen_itensor(d, expo, &U2[i]);


    //Left bond
    LindLeft += Jump_site(dim, 2, 0, 0.5*GammaL1*(1+muL1), Su, Id);
    LindLeft += Jump_site(dim, 2, 0, 0.5*GammaL1*(1-muL1), Sd, Id);
    LindLeft += Jump_fermion_site(2, 1, 0, 0.5*GammaL2*(1+muL2), Su, Id, Sz);
    LindLeft += Jump_fermion_site(2, 1, 0, 0.5*GammaL2*(1-muL2), Sd, Id, Sz);

    expo = LindLeft * complex<double>(0.25, -sqrt(3)/12.0);
    expo = matrix_complex_exponential(d*d, dt, expo);
    TEvol_Eigen_itensor(d, expo, &U[0]);

    expo = LindLeft / 2;
    expo = matrix_complex_exponential(d*d, dt, expo);
    TEvol_Eigen_itensor(d, expo, &U1[0]);

    expo = LindLeft * complex<double>(0.25, sqrt(3)/12.0);
    expo = matrix_complex_exponential(d*d, dt, expo);
    TEvol_Eigen_itensor(d, expo, &U2[0]);


    //Right bond
    LindRight += Hamiltonian_bond(dim, 2, 2, 3, tr, Su, Sd, Id);
    LindRight += Hamiltonian_bond(dim, 2, 2, 3, tr, Sd, Su, Id);
    LindRight += Hamiltonian_bond(dim, 2, 2, 3, Vr, nf, nf, Id);
    LindRight += Hamiltonian_site(dim, 2, 2, mu, nf, Id);
    LindRight += Hamiltonian_site(dim, 2, 3, mu, nf, Id);
    LindRight += Jump_site(dim, 2, 2, gamma, nf, Id);
    LindRight += Jump_site(dim, 2, 3, gamma, nf, Id);
    if(L%2==0){
        LindRight += Jump_site(dim, 2, 2, 0.5*GammaR1*(1-muR1), Su, Id);
        LindRight += Jump_site(dim, 2, 2, 0.5*GammaR1*(1+muR1), Sd, Id);

        LindRight += Jump_fermion_site(2, 3, 2, 0.5*GammaR2*(1-muR2), Su, Id, Sz);
        LindRight += Jump_fermion_site(2, 3, 2, 0.5*GammaR2*(1+muR2), Sd, Id, Sz);

        expo = LindRight * complex<double>(0.25, -sqrt(3)/12.0);
        expo = matrix_complex_exponential(d*d, dt, expo);
        TEvol_Eigen_itensor(d, expo, &U[L-2]);

        expo = LindRight / 2;
        expo = matrix_complex_exponential(d*d, dt, expo);
        TEvol_Eigen_itensor(d, expo, &U1[L-2]);

        expo = LindRight * complex<double>(0.25, sqrt(3)/12.0);
        expo = matrix_complex_exponential(d*d, dt, expo);
        TEvol_Eigen_itensor(d, expo, &U2[L-2]);
    }
    else{
        LindRight += Jump_fermion_site(2, 2, 3, 0.5*GammaR1*(1-muR1), Su, Id, Sz);
        LindRight += Jump_fermion_site(2, 2, 3, 0.5*GammaR1*(1+muR1), Sd, Id, Sz);

        LindRight += Jump_site(dim, 2, 3, 0.5*GammaR2*(1-muR2), Su, Id);
        LindRight += Jump_site(dim, 2, 3, 0.5*GammaR2*(1+muR2), Sd, Id);

        expo = LindRight * complex<double>(0.5, -sqrt(3)/6.0);
        expo = matrix_complex_exponential(d*d, dt, expo);
        TEvol_Eigen_itensor(d, expo, &U[L-2]);

        expo = LindRight * complex<double>(0.5, sqrt(3)/6.0);
        expo = matrix_complex_exponential(d*d, dt, expo);
        TEvol_Eigen_itensor(d, expo, &U1[L-2]);
    }
}


//Function that writes the Lindblad operator as a MPO
void Ladder_open::Lindbladian_MPO (){

    using namespace Eigen;

 /*   int dim=2;
    MatrixXcd LBulk = MatrixXcd::Zero(d*B, d*B), LLeft = MatrixXcd::Zero(d, d*B), LRight = MatrixXcd::Zero(d*B, d);
    MatrixXcd Sx(dim,dim), Sy(dim,dim), Sz(dim,dim), Su(dim,dim), Sd(dim,dim), Id(dim,dim);
    Pauli_Matrices(dim, Sx, Sy, Sz, Su, Sd, Id);

    //Initialize
    LBulk.block(0, 0, d, d) = MatrixXd::Identity(d, d);
    LBulk.block((B-1)*d, (B-1)*d, d, d) = MatrixXd::Identity(d, d);
    LLeft.block(0, (B-1)*d, d, d) = MatrixXd::Identity(d, d);
    LRight.block(0, 0, d, d) = MatrixXd::Identity(d, d);


    //Bonds
    Lindblad_bond(dim, B, 2, 0, 1, LBulk, LLeft, LRight, Jl, Sx, Id);
    Lindblad_bond(dim, B, 2, 1, 3, LBulk, LLeft, LRight, Jl, Sx, Id);
    Lindblad_bond(dim, B, 2, 0, 5, LBulk, LLeft, LRight, Jl, Sy, Id);
    Lindblad_bond(dim, B, 2, 1, 7, LBulk, LLeft, LRight, Jl, Sy, Id);
    Lindblad_bond(dim, B, 2, 0, 9, LBulk, LLeft, LRight, deltal, Sz, Id);
    Lindblad_bond(dim, B, 2, 1, 11, LBulk, LLeft, LRight, deltal, Sz, Id);
    Lindblad_rung(dim, B, 2, 0, 1, LBulk, LLeft, LRight, Jr, Sx, Sx, Id);
    Lindblad_rung(dim, B, 2, 0, 1, LBulk, LLeft, LRight, Jr, Sy, Sy, Id);
    Lindblad_rung(dim, B, 2, 0, 1, LBulk, LLeft, LRight, deltar, Sz, Sz, Id);

    //Sites
    Lindblad_site(dim, B, 2, 0, LBulk, LLeft, LRight, h, Sz, Id);
    Lindblad_site(dim, B, 2, 1, LBulk, LLeft, LRight, h, Sz, Id);
    Lindblad_jump(dim, B, 2, 0, LBulk, LLeft, LRight, gamma, Sz, Id, 0);
    Lindblad_jump(dim, B, 2, 1, LBulk, LLeft, LRight, gamma, Sz, Id, 0);

    //Reservoirs
    Lindblad_jump(dim, B, 2, 0, LBulk, LLeft, LRight, 0.5*GammaL1*(1+muL1), Su, Id, 1);
    Lindblad_jump(dim, B, 2, 0, LBulk, LLeft, LRight, 0.5*GammaL1*(1-muL1), Sd, Id, 1);
    Lindblad_jump(dim, B, 2, 1, LBulk, LLeft, LRight, 0.5*GammaL2*(1+muL2), Su, Id, 1);
    Lindblad_jump(dim, B, 2, 1, LBulk, LLeft, LRight, 0.5*GammaL2*(1-muL2), Sd, Id, 1);
    Lindblad_jump(dim, B, 2, 0, LBulk, LLeft, LRight, 0.5*GammaR1*(1-muR1), Su, Id, 2);
    Lindblad_jump(dim, B, 2, 0, LBulk, LLeft, LRight, 0.5*GammaR1*(1+muR1), Sd, Id, 2);
    Lindblad_jump(dim, B, 2, 1, LBulk, LLeft, LRight, 0.5*GammaR2*(1-muR2), Su, Id, 2);
    Lindblad_jump(dim, B, 2, 1, LBulk, LLeft, LRight, 0.5*GammaR2*(1+muR2), Sd, Id, 2);


    LLeft_Eigen_itensor(d, B, LLeft, &Lindblad.ref(1));
    for(int i=2;i<=L-1;++i)
        LBulk_Eigen_itensor(d, B, LBulk, &Lindblad.ref(i));
    LRight_Eigen_itensor(d, B, LRight, &Lindblad.ref(L));*/
}


void Ladder_open::Measure(int t, int Nt, double Time){

    int bond;
    double Lind, Lind2, N, J, Mag;
    N = Norm(Id);

    bond = maxLinkDim(M);
    //Lind = abs(innerC(M, Lindblad, prime(M)));
    //Lind2 = abs(innerC( prime(Lindblad, "Link"), prime(M), Lindblad, prime(M)));
    //Mag = local_observable(L/4+1, nf_up, Id)/N;
    if((L/4+1)%2==1)
        J= tl*(two_point_correlation(L/4+1, L/4+2, Sx_up_odd, Sy_up_odd, Id) - two_point_correlation(L/4+1, L/4+2, Sy_up_odd, Sx_up_odd, Id))/2/N;
    else
        J= tl*(two_point_correlation(L/4+1, L/4+2, Sx_up_even, Sy_up_even, Id) - two_point_correlation(L/4+1, L/4+2, Sy_up_even, Sx_up_even, Id))/2/N;

    monitor[3*(t%Nt)] = Time;
    monitor[3*(t%Nt) + 1] = J;
    monitor[3*(t%Nt) + 2] = bond;
}


void Ladder_open::Monitor(int Nt){

    char* buffer;
    char name[100];

    sprintf(name, "output/open/Monitor_%lu.bin", ID);
    outMonitor.open(name, ios::out | ios::binary | ios::app);

    buffer = reinterpret_cast<char*>(monitor);
    outMonitor.write(buffer, 3*Nt*sizeof(double));

    outMonitor.close();
}


void Ladder_open::Measure_final(){

    char* buffer;
    char name[100];
    ITensor O1, O2, O3, O4;
    double N, Mag_up[L], Mag_down[L], J_up[L-1], J_down[L-1], J[L], corr_up[L][L], corr_down[L][L], chi[L-1];
    complex<double> cov[2*L][2*L], cov_Full[4*L][4*L];
    N = Norm(Id);

    sprintf(name, "output/open/Mup_%lu.bin", ID);
    outM_up.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Mdown_%lu.bin", ID);
    outM_down.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Jup_%lu.bin", ID);
    outJ_up.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Jdown_%lu.bin", ID);
    outJ_down.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Jrung_%lu.bin", ID);
    outJ_rung.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Corrup_%lu.bin", ID);
    out_corr_up.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Corrdown_%lu.bin", ID);
    out_corr_down.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Cov_%lu.bin", ID);
    out_cov.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Cov_Full_%lu.bin", ID);
    out_covF.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Jsus_%lu.bin", ID);
    out_Jsus.open(name, ios::out | ios::binary);


    for(int loc=1; loc<=L-1; ++loc){

        Mag_up[loc-1] = local_observable(loc, nf_up, Id)/N;
        Mag_down[loc-1] = local_observable(loc, nf_down, Id)/N;
        chi[loc-1] = 2.0*tl*tl*two_point_correlation(loc, loc+1, nf_up, nf_up, Id)/N;

        if(loc%2==1){
            J_up[loc-1] = tl*(two_point_correlation(loc, loc+1, Sx_up_odd, Sy_up_odd, Id) - two_point_correlation(loc, loc+1, Sy_up_odd, Sx_up_odd, Id))/2/N;
            J_down[loc-1] = tl*(two_point_correlation(loc, loc+1, Sx_down_odd, Sy_down_odd, Id) - two_point_correlation(loc, loc+1, Sy_down_odd, Sx_down_odd, Id))/2/N;
        }
        else{
            J_up[loc-1] = tl*(two_point_correlation(loc, loc+1, Sx_up_even, Sy_up_even, Id) - two_point_correlation(loc, loc+1, Sy_up_even, Sx_up_even, Id))/2/N;
            J_down[loc-1] = tl*(two_point_correlation(loc, loc+1, Sx_down_even, Sy_down_even, Id) - two_point_correlation(loc, loc+1, Sy_down_even, Sx_down_even, Id))/2/N;
        }
        J[loc-1] = tr*local_observable(loc, J_rung, Id)/N;

        for(int loc2=0; loc2<L; ++loc2){

            if(loc != loc2+1){
              corr_up[loc-1][loc2] = two_point_correlation(loc, loc2+1, nf_up, nf_up, Id)/N;
              corr_down[loc-1][loc2] = two_point_correlation(loc, loc2+1, nf_down, nf_down, Id)/N;
            }
            else{
              corr_up[loc-1][loc2] = 1.0;
              corr_down[loc-1][loc2] = 1.0;
            }


            if(loc < loc2+1){

                if( loc%2==1){
                    O1 = Sd_up_odd;
                    O3 = Sd_down_odd;
                }
                else{
                    O1 = Sd_up_even;
                    O3 = Sd_down_even;
                }
                if( (loc2+1)%2==1){
                    O2 = Su_up_odd;
                    O4 = Su_down_odd;
                }
                else{
                    O2 = Su_up_even;
                    O4 = Su_down_even;
                }

                cov[loc-1][loc2] = -two_point_correlation_string(loc, loc2+1, O1, O2, Id, Szm)/N;
                cov[loc2][loc-1] = conj(cov[loc-1][loc2]);
                cov[L + loc-1][L + loc2] = -cov[loc2][loc-1];
                cov[L + loc2][L + loc-1] = -cov[loc-1][loc2];

                cov_Full[loc-1][loc2] = cov[loc-1][loc2];
                cov_Full[loc2][loc-1] = cov[loc2][loc-1];
                cov_Full[2*L + loc-1][2*L + loc2] = cov[L + loc-1][L + loc2];
                cov_Full[2*L + loc2][2*L + loc-1] = cov[L + loc2][L + loc-1];
                cov_Full[L + loc-1][L + loc2] = -two_point_correlation_string(loc, loc2+1, O3, O4, Id, Szm)/N;
                cov_Full[L + loc2][L + loc-1] = conj(cov_Full[L + loc-1][L + loc2]);
                cov_Full[3*L + loc-1][3*L + loc2] = -cov_Full[L + loc2][L + loc-1];
                cov_Full[3*L + loc2][3*L + loc-1] = -cov_Full[L + loc-1][L + loc2];

                cov_Full[loc-1][L + loc2] = -two_point_correlation_string(loc, loc2+1, O1, O4, Id, Szm)/N;
                cov_Full[L + loc-1][loc2] = -two_point_correlation_string(loc, loc2+1, O3, O2, Id, Szm)/N;
                cov_Full[L + loc2][loc-1] = conj(cov_Full[loc-1][L + loc2]);
                cov_Full[loc2][L + loc-1] = conj(cov_Full[L + loc-1][loc2]);
                cov_Full[2*L + loc-1][3*L + loc2] = -cov_Full[L + loc2][loc-1];
                cov_Full[2*L + loc2][3*L + loc-1] = -cov_Full[L + loc-1][loc2];
                cov_Full[3*L + loc-1][2*L + loc2] = -cov_Full[loc2][L + loc-1];
                cov_Full[3*L + loc2][2*L + loc-1] = -cov_Full[loc-1][L + loc2];


                if( (loc2+1)%2==1)
                    O2 = Sd_up_even;
                else
                    O2 = Sd_up_odd;

                cov[loc-1][L + loc2] = -two_point_correlation_string(loc, loc2+1, O1, O2, Id, Szm)/N;
                cov[loc2][L + loc-1] = -cov[loc-1][L + loc2];
                cov[L + loc-1][loc2] = conj(cov[loc2][L + loc-1]);
                cov[L + loc2][loc-1] = conj(cov[loc-1][L + loc2]);
            }
            if(loc == loc2+1){
                cov[loc-1][loc-1] = 1-local_observable(loc, nf_up, Id)/N;
                cov[L + loc-1][L + loc-1] = 1 - cov[loc-1][loc-1];
                cov_Full[loc-1][loc-1] = cov[loc-1][loc-1];
                cov_Full[2*L + loc-1][2*L + loc-1] = cov[L + loc-1][L + loc-1];
                cov_Full[L + loc-1][L + loc-1] = 1-local_observable(loc, nf_down, Id)/N;
                cov_Full[3*L + loc-1][3*L + loc-1] = 1 - cov_Full[L + loc-1][L + loc-1];

                cov_Full[loc-1][L + loc-1] = -local_observable(loc, ST, Id)/N;
                cov_Full[L + loc-1][loc-1] = conj(cov_Full[loc-1][L + loc-1]);
                cov_Full[2*L + loc-1][3*L + loc-1] = -cov_Full[L + loc-1][loc-1];
                cov_Full[3*L + loc-1][2*L + loc-1] = -cov_Full[loc-1][L + loc-1];
            }
        }
    }

    Mag_up[L-1] = local_observable(L, nf_up, Id)/N;
    Mag_down[L-1] = local_observable(L, nf_down, Id)/N;
    J[L-1] = tr*local_observable(L, J_rung, Id)/N;
    for(int loc2=0; loc2<L; ++loc2){

        if(L != loc2+1){
          corr_up[L-1][loc2] = two_point_correlation(L, loc2+1, nf_up, nf_up, Id)/N;
          corr_down[L-1][loc2] = two_point_correlation(L, loc2+1, nf_down, nf_down, Id)/N;
        }
        else{
          corr_up[L-1][loc2] = 1.0;
          corr_down[L-1][loc2] = 1.0;
        }

        if(loc2+1 == L){
            cov[L-1][L-1] = 1-local_observable(L, nf_up, Id)/N;
            cov[2*L-1][2*L-1] = 1 - cov[L-1][L-1];
            cov_Full[L-1][L-1] = cov[L-1][L-1];
            cov_Full[3*L-1][3*L-1] = cov[2*L-1][2*L-1];
            cov_Full[2*L-1][2*L-1] = 1-local_observable(L, nf_down, Id)/N;
            cov_Full[4*L-1][4*L-1] = 1 - cov_Full[2*L-1][2*L-1];

            cov_Full[L-1][2*L -1] = -local_observable(L, ST, Id)/N;
            cov_Full[2*L-1][L -1] = conj(cov_Full[L-1][2*L -1]);
            cov_Full[3*L -1][4*L -1] = -cov_Full[2*L -1][L-1];
            cov_Full[4*L -1][3*L -1] = -cov_Full[L-1][2*L -1];
        }
    }


    buffer = reinterpret_cast<char*>(Mag_up);
    outM_up.write(buffer, L*sizeof(double));
    buffer = reinterpret_cast<char*>(Mag_down);
    outM_down.write(buffer, L*sizeof(double));
    buffer = reinterpret_cast<char*>(J);
    outJ_rung.write(buffer, L*sizeof(double));
    buffer = reinterpret_cast<char*>(J_up);
    outJ_up.write(buffer, (L-1)*sizeof(double));
    buffer = reinterpret_cast<char*>(J_down);
    outJ_down.write(buffer, (L-1)*sizeof(double));
    buffer = reinterpret_cast<char*>(corr_up);
    out_corr_up.write(buffer, L*L*sizeof(double));
    buffer = reinterpret_cast<char*>(corr_down);
    out_corr_down.write(buffer, L*L*sizeof(double));
    buffer = reinterpret_cast<char*>(cov);
    out_cov.write(buffer, 4*L*L*sizeof(complex<double>));
    buffer = reinterpret_cast<char*>(cov_Full);
    out_covF.write(buffer, 16*L*L*sizeof(complex<double>));
    buffer = reinterpret_cast<char*>(chi);
    out_Jsus.write(buffer, (L-1)*sizeof(double));


    outM_up.close();
    outM_down.close();
    outJ_up.close();
    outJ_down.close();
    outJ_rung.close();
    out_corr_up.close();
    out_corr_down.close();
    out_cov.close();
    out_covF.close();
    out_Jsus.close();
}


void Ladder_open::Singular_values (){

    char* buffer;
    char name[100];
    M.position(L/2);
    auto phi = M(L/2) * M(L/2+1);
    ITensor U(inds(M(L/2))), D, V;

    Spectrum spec = svd(phi, U, D, V);

    sprintf(name, "output/open/Spectrum_%lu.bin", ID);
    outSpectrum.open(name, ios::out | ios::binary);

    for(auto eig : spec.eigsKept()){
        buffer = reinterpret_cast<char*>(&eig);
        outSpectrum.write(buffer, sizeof(double));
    }
    outSpectrum.close();
}
