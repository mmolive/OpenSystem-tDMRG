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


XXZ_open::XXZ_open (int Lin, int din, int Din, int Bin, double Jin, double hin, double deltain, double gammain, double Gammain, double muin, int Nt, int restore, int D_old, double delta_old, double mu_old) : open(Lin, din, Din) {

    B = Bin;
    J = Jin;
    h = hin;
    delta = deltain;
    gamma = gammain;
    Gamma = Gammain;
    mu = muin;
    monitor = new double[3*Nt];

    auto i = Index(d);
    Id = ITensor(i);
    Sx = ITensor(i);
    Sy = ITensor(i);
    Sz = ITensor(i);
    Szm = ITensor(i);
    Su = ITensor(i);
    Sd = ITensor(i);
    SdSz = ITensor(i);

    Id.set(1, 1.0); Id.set(2, 0.0); Id.set(3, 0.0); Id.set(4, 1.0);
    Sx.set(1, 0.0); Sx.set(2, 1.0); Sx.set(3, 1.0); Sx.set(4, 0.0);
    Sy.set(1, 0.0); Sy.set(2, 1.0_i); Sy.set(3, -1.0_i); Sy.set(4, 0.0);
    Sz.set(1, 1.0); Sz.set(2, 0.0); Sz.set(3, 0.0); Sz.set(4, -1.0);
    Szm.set(1, -1.0); Szm.set(2, 0.0); Szm.set(3, 0.0); Szm.set(4, 1.0);
    Su.set(1, 0.0); Su.set(2, 0.0); Su.set(3, 1.0); Su.set(4, 0.0);
    Sd.set(1, 0.0); Sd.set(2, 1.0); Sd.set(3, 0.0); Sd.set(4, 0.0);
    SdSz.set(1, 0.0); SdSz.set(2, -1.0); SdSz.set(3, 0.0); SdSz.set(4, 0.0);


    ID = ID_XXZ(D, delta, mu);
    char name[100];

    //input older state
    if(restore==1){

        uint64_t ID_res = ID_XXZ(D_old, delta_old, mu_old);
        sprintf(name, "output/open/Final_state_%lu.bin", ID_res);
        readFromFile(name, M);

        if(ID==ID_res){
            ifstream inTime;
            sprintf(name, "output/open/Time_%lu.txt", ID);
            inTime.open(name);
            inTime >> Time;
            inTime.close();
        }
        else{
            Time=0;
            sprintf(name, "output/open/Monitor_%lu.bin", ID);
            outMonitor.open(name, ios::out | ios::binary);
            outMonitor.close();
        }

        for(int i=0; i<L-1; ++i){
            U[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
            U1[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
            U2[i] = ITensor(siteIndex(M, i+1).prime(), siteIndex(M, i+2).prime(), siteIndex(M, i+1), siteIndex(M, i+2));
          }
    }
    else{
        Time=0;
        init_all_down();
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

XXZ_open::~XXZ_open (){
    delete[] monitor;
}


//Sets initial state with all spins pointing down
void XXZ_open::init_all_down (){

    M.ref(1).set(1,1, 0.0);
    M.ref(1).set(2,1, 0.0);
    M.ref(1).set(3,1, 0.0);
    M.ref(1).set(4,1, 1.0);

    for(int i=2; i<L; ++i){

        M.ref(i).set(1,1,1, 0.0);
        M.ref(i).set(1,2,1, 0.0);
        M.ref(i).set(1,3,1, 0.0);
        M.ref(i).set(1,4,1, 1.0);
    }

    M.ref(L).set(1,1, 0.0);
    M.ref(L).set(1,2, 0.0);
    M.ref(L).set(1,3, 0.0);
    M.ref(L).set(1,4, 1.0);
}


//gets unique ID computed from simulation parameters
uint64_t XXZ_open::ID_XXZ (int D_in, double delta_in, double mu_in){

    uint64_t hash;
    hash = 14695981039346656037;

    hash = FNVHash<int>(L, hash);
    hash = FNVHash<int>(D_in, hash);
    hash = FNVHash<double>(J, hash);
    hash = FNVHash<double>(delta_in, hash);
    hash = FNVHash<double>(h, hash);
    hash = FNVHash<double>(gamma, hash);
    hash = FNVHash<double>(mu_in, hash);
    hash = FNVHash<double>(Gamma, hash);

    return hash;
}


//Function that writes the time-evolution operators
void XXZ_open::time_evolution_operator(double dt, int order){

    using namespace Eigen;

    int dim=2;
    MatrixXcd LindBulk, LindLeft, LindRight, LindAux, expo[2];
    MatrixXcd Sx(dim,dim), Sy(dim,dim), Sz(dim,dim), Su(dim,dim), Sd(dim,dim), Id(dim,dim);
    Pauli_Matrices(d, Sx, Sy, Sz, Su, Sd, Id);


    //Bulk
    LindBulk = Hamiltonian_bond(dim, 1, 0, 1, J, Sx, Sx, Id);
    LindBulk += Hamiltonian_bond(dim, 1, 0, 1, J, Sy, Sy, Id);
    LindBulk += Hamiltonian_bond(dim, 1, 0, 1, delta, Sz, Sz, Id);
    LindBulk += Hamiltonian_site(dim, 1, 0, h, Sz, Id);
    LindBulk += Jump_site(dim, 1, 0, gamma, Sz, Id);
    LindAux = LindBulk/2;

    expo[0] = matrix_complex_exponential(d*d, dt, LindBulk);
    if(order==2)
        expo[1] = matrix_complex_exponential(d*d, dt, LindAux);
    else
        expo[1] = matrix_complex_exponential(d*d, dt, LindBulk);

    for(int i=1; i<L-2; ++i)
        TEvol_Eigen_itensor(d, expo[i%2], &U[i]);

    LindLeft = LindBulk;
    LindRight = LindBulk;


    //Left
    LindLeft += Jump_site(dim, 1, 0, 0.5*Gamma*(1+mu), Su, Id);
    LindLeft += Jump_site(dim, 1, 0, 0.5*Gamma*(1-mu), Sd, Id);

    expo[0] = matrix_complex_exponential(d*d, dt, LindLeft);
    TEvol_Eigen_itensor(d, expo[0], &U[0]);


    //Right
    LindRight += Hamiltonian_site(dim, 1, 1, h, Sz, Id);
    LindRight += Jump_site(dim, 1, 1, gamma, Sz, Id);
    LindRight += Jump_site(dim, 1, 1, 0.5*Gamma*(1-mu), Su, Id);
    LindRight += Jump_site(dim, 1, 1, 0.5*Gamma*(1+mu), Sd, Id);
    LindAux = LindRight/2;

    if(order==2 && L%2==1)
        expo[0] = matrix_complex_exponential(d*d, dt, LindAux);
    else
        expo[0] = matrix_complex_exponential(d*d, dt, LindRight);

    TEvol_Eigen_itensor(d, expo[0], &U[L-2]);

   /*ofstream out;
    out.open("input_MPO/test_real.txt");
    for(int i=0; i<d; ++i){
        for(int j=0; j<d; ++j){
            for(int k=0; k<d; ++k){
                for(int l=0; l<d; ++l)
                    out << elt(realPart(U[2]), i+1, j+1, k+1, l+1) << " ";
            }
            out << endl;
        }
    }
    out.close();

    out.open("input_MPO/test_imag.txt");
    for(int i=0; i<d; ++i){
        for(int j=0; j<d; ++j){
            for(int k=0; k<d; ++k){
                for(int l=0; l<d; ++l)
                    out << elt(imagPart(U[2]), i+1, j+1, k+1, l+1) << " ";
            }
            out << endl;
        }
    }
    out.close();*/

}


//Function that writes the time-evolution operators in 4th order
void XXZ_open::time_evolution_operator_4th(double dt){

  using namespace Eigen;

  int dim=2;
  MatrixXcd LindBulk, LindLeft, LindRight, expo;
  MatrixXcd Sx(dim,dim), Sy(dim,dim), Sz(dim,dim), Su(dim,dim), Sd(dim,dim), Id(dim,dim);
  Pauli_Matrices(d, Sx, Sy, Sz, Su, Sd, Id);


  //Bulk even bonds
  LindBulk = Hamiltonian_bond(dim, 1, 0, 1, J, Sx, Sx, Id);
  LindBulk += Hamiltonian_bond(dim, 1, 0, 1, J, Sy, Sy, Id);
  LindBulk += Hamiltonian_bond(dim, 1, 0, 1, delta, Sz, Sz, Id);
  LindBulk += Hamiltonian_site(dim, 1, 0, h, Sz, Id);
  LindBulk += Jump_site(dim, 1, 0, gamma, Sz, Id);

  expo = LindBulk * complex<double>(0.5, -sqrt(3)/6.0);
  expo = matrix_complex_exponential(d*d, dt, expo);
  for(int i=1; i<L-2; i+=2)
      TEvol_Eigen_itensor(d, expo, &U[i]);

  expo = LindBulk * complex<double>(0.5, sqrt(3)/6.0);
  expo = matrix_complex_exponential(d*d, dt, expo);
  for(int i=1; i<L-2; i+=2)
      TEvol_Eigen_itensor(d, expo, &U1[i]);


  //Bulk odd bonds
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

  LindLeft = LindBulk;
  LindRight = LindBulk;


  //Left bond
  LindLeft += Jump_site(dim, 1, 0, 0.5*Gamma*(1+mu), Su, Id);
  LindLeft += Jump_site(dim, 1, 0, 0.5*Gamma*(1-mu), Sd, Id);

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
  LindRight += Hamiltonian_site(dim, 1, 1, h, Sz, Id);
  LindRight += Jump_site(dim, 1, 1, gamma, Sz, Id);
  LindRight += Jump_site(dim, 1, 1, 0.5*Gamma*(1-mu), Su, Id);
  LindRight += Jump_site(dim, 1, 1, 0.5*Gamma*(1+mu), Sd, Id);
  if(L%2==0){
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
      expo = LindRight * complex<double>(0.5, -sqrt(3)/6.0);
      expo = matrix_complex_exponential(d*d, dt, expo);
      TEvol_Eigen_itensor(d, expo, &U[L-2]);

      expo = LindRight * complex<double>(0.5, sqrt(3)/6.0);
      expo = matrix_complex_exponential(d*d, dt, expo);
      TEvol_Eigen_itensor(d, expo, &U1[L-2]);
  }
}


void XXZ_open::Lindbladian_MPO (){

    using namespace Eigen;

    int dim=2;
    MatrixXcd LBulk = MatrixXcd::Zero(d*B, d*B), LLeft = MatrixXcd::Zero(d, d*B), LRight = MatrixXcd::Zero(d*B, d);
    MatrixXcd Sx(dim,dim), Sy(dim,dim), Sz(dim,dim), Su(dim,dim), Sd(dim,dim), Id(dim,dim);
    Pauli_Matrices(dim, Sx, Sy, Sz, Su, Sd, Id);

    //Initialize
    LBulk.block(0, 0, d, d) = MatrixXd::Identity(d, d);
    LBulk.block((B-1)*d, (B-1)*d, d, d) = MatrixXd::Identity(d, d);
    LLeft.block(0, (B-1)*d, d, d) = MatrixXd::Identity(d, d);
    LRight.block(0, 0, d, d) = MatrixXd::Identity(d, d);


    //Bonds
    Lindblad_bond(dim, B, 1, 0, 1, LBulk, LLeft, LRight, J, Sx, Id);
    Lindblad_bond(dim, B, 1, 0, 3, LBulk, LLeft, LRight, J, Sy, Id);
    Lindblad_bond(dim, B, 1, 0, 5, LBulk, LLeft, LRight, delta, Sz, Id);

    //Sites
    Lindblad_site(dim, B, 1, 0, LBulk, LLeft, LRight, h, Sz, Id);
    Lindblad_jump(dim, B, 1, 0, LBulk, LLeft, LRight, gamma, Sz, Id, 0);

    //Reservoirs
    Lindblad_jump(dim, B, 1, 0, LBulk, LLeft, LRight, 0.5*Gamma*(1+mu), Su, Id, 1);
    Lindblad_jump(dim, B, 1, 0, LBulk, LLeft, LRight, 0.5*Gamma*(1-mu), Sd, Id, 1);
    Lindblad_jump(dim, B, 1, 0, LBulk, LLeft, LRight, 0.5*Gamma*(1-mu), Su, Id, 2);
    Lindblad_jump(dim, B, 1, 0, LBulk, LLeft, LRight, 0.5*Gamma*(1+mu), Sd, Id, 2);


    LLeft_Eigen_itensor(d, B, LLeft, &Lindblad.ref(1));
    for(int i=2;i<=L-1;++i)
        LBulk_Eigen_itensor(d, B, LBulk, &Lindblad.ref(i));
    LRight_Eigen_itensor(d, B, LRight, &Lindblad.ref(L));



    /*ofstream out;
    out.open("input_MPO/test_real.txt");
    for(int i=0; i<B; ++i){
        for(int j=0; j<d; ++j){
            for(int k=0; k<B; ++k){
                for(int l=0; l<d; ++l)
                    out << elt(realPart( Lindblad(2) ), j+1, l+1, i+1, k+1) << " ";
            }
            out << endl;
        }
    }
    out.close();

    out.open("input_MPO/test_imag.txt");
    for(int i=0; i<B; ++i){
        for(int j=0; j<d; ++j){
            for(int k=0; k<B; ++k){
                for(int l=0; l<d; ++l)
                    out << elt(imagPart( Lindblad(2) ), j+1, l+1, i+1, k+1) << " ";
            }
            out << endl;
        }
    }
    out.close();


   ofstream out;
   out.open("input_MPO/test_real.txt");
   for(int i=0; i<d; ++i){
       for(int j=0; j<B; ++j){
           for(int k=0; k<d; ++k)
                   out << elt(realPart( Lindblad(1) ), i+1, k+1, j+1) << " ";
       }
       out << endl;
   }
   out.close();

   out.open("input_MPO/test_imag.txt");
   for(int i=0; i<d; ++i){
       for(int j=0; j<B; ++j){
           for(int k=0; k<d; ++k)
                   out << elt(imagPart( Lindblad(1) ), i+1, k+1, j+1) << " ";
       }
       out << endl;
   }
   out.close();


   ofstream out;
   out.open("input_MPO/test_real.txt");
   for(int i=0; i<B; ++i){
       for(int j=0; j<d; ++j){
           for(int k=0; k<d; ++k)
                   out << elt(realPart( Lindblad(L) ), j+1, k+1, i+1) << " ";
           out << endl;
       }
   }
   out.close();

   out.open("input_MPO/test_imag.txt");
   for(int i=0; i<B; ++i){
       for(int j=0; j<d; ++j){
           for(int k=0; k<d; ++k)
                   out << elt(imagPart( Lindblad(L) ), j+1, k+1, i+1) << " ";
           out << endl;
       }
   }
   out.close(); */




}


void XXZ_open::Measure(int t, int Nt, double Time){

    int bond;
    double Lind, Lind2, N, Ji;
    N = Norm(Id);

    bond = maxLinkDim(M);
    //Lind = abs(innerC(M, Lindblad, prime(M)));
    //Lind2 = abs(innerC( prime(Lindblad, "Link"), prime(M), Lindblad, prime(M)));
    Ji= J*(two_point_correlation(L/2, L/2+1, Sx, Sy, Id) - two_point_correlation(L/2, L/2+1, Sy, Sx, Id))/N;

    monitor[3*(t%Nt)] = Time;
    monitor[3*(t%Nt) + 1] = Ji;
    monitor[3*(t%Nt) + 2] = bond;


    /*char* buffer = reinterpret_cast<char*>(&t);
    outJ.write(buffer, sizeof(double));
    outBond.write(buffer, sizeof(double));
    outL.write(buffer, sizeof(double));

    J= (2*two_point_correlation(L/2, L/2+1, Sx, Sy, Id) - 2*two_point_correlation(L/2, L/2+1, Sy, Sx, Id))/N;
    buffer = reinterpret_cast<char*>(&J);
    outJ.write(buffer, sizeof(double));

    bond = maxLinkDim(M);
    buffer = reinterpret_cast<char*>(&bond);
    outBond.write(buffer, sizeof(int));

    Lind = abs(innerC(M, Lindblad, M));
    buffer = reinterpret_cast<char*>(&Lind);
    outL.write(buffer, sizeof(double));*/
}


void XXZ_open::Monitor(int Nt){

    char* buffer;
    char name[100];

    sprintf(name, "output/open/Monitor_%lu.bin", ID);
    outMonitor.open(name, ios::out | ios::binary | ios::app);

    buffer = reinterpret_cast<char*>(monitor);
    outMonitor.write(buffer, 3*Nt*sizeof(double));

    outMonitor.close();
}


void XXZ_open::Measure_final(){

    char* buffer;
    char name[100];
    double N, Mag[L], Ji[L-1], chi[L-1], corr[L][L];
    complex<double> cov[2*L][2*L];
    N = Norm(Id);

    sprintf(name, "output/open/M_%lu.bin", ID);
    outM.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/J_%lu.bin", ID);
    outJ.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Corr_%lu.bin", ID);
    out_corr.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Cov_%lu.bin", ID);
    out_cov.open(name, ios::out | ios::binary);
    sprintf(name, "output/open/Jsus_%lu.bin", ID);
    out_Jsus.open(name, ios::out | ios::binary);


    for(int loc=1; loc<=L-1; ++loc){
        Mag[loc-1] = local_observable(loc, Sz, Id)/N;
        Ji[loc-1] = J*(two_point_correlation(loc, loc+1, Sx, Sy, Id) - two_point_correlation(loc, loc+1, Sy, Sx, Id))/N;
        chi[loc-1] = 2.0*J*J*(1 - two_point_correlation(loc, loc+1, Sz, Sz, Id)/N);

        for(int loc2=0; loc2<L; ++loc2){

            if(loc != loc2+1)
              corr[loc-1][loc2] = two_point_correlation(loc, loc2+1, Sz, Sz, Id)/N;
            else
              corr[loc-1][loc2] = 1.0;

            if(loc < loc2+1){
                cov[loc-1][loc2] = two_point_correlation_string(loc, loc2+1, SdSz, Su, Id, Szm)/N;
                cov[loc2][loc-1] = conj(cov[loc-1][loc2]);
                cov[L + loc-1][L + loc2] = -cov[loc2][loc-1];
                cov[L + loc2][L + loc-1] = -cov[loc-1][loc2];
                cov[loc-1][L + loc2] = two_point_correlation_string(loc, loc2+1, SdSz, Sd, Id, Szm)/N;
                cov[loc2][L + loc-1] = -cov[loc-1][L + loc2];
                cov[L + loc-1][loc2] = conj(cov[loc2][L + loc-1]);
                cov[L + loc2][loc-1] = conj(cov[loc-1][L + loc2]);
            }
            if(loc == loc2+1){
                cov[loc-1][loc-1] = (1-local_observable(loc, Sz, Id)/N)/2.0;
                cov[L + loc-1][L + loc-1] = 1 - cov[loc-1][loc-1];
            }
        }
    }
    Mag[L-1] = local_observable(L, Sz, Id)/N;
    for(int loc2=0; loc2<L; ++loc2){

        if(L != loc2+1)
          corr[L-1][loc2] = two_point_correlation(L, loc2+1, Sz, Sz, Id)/N;
        else
          corr[L-1][loc2] = 1.0;

        if(loc2+1 == L){
            cov[L-1][L-1] = (1-local_observable(L, Sz, Id)/N)/2.0;
            cov[2*L-1][2*L-1] = 1 - cov[L-1][L-1];
        }
    }


    buffer = reinterpret_cast<char*>(Mag);
    outM.write(buffer, L*sizeof(double));
    buffer = reinterpret_cast<char*>(Ji);
    outJ.write(buffer, (L-1)*sizeof(double));
    buffer = reinterpret_cast<char*>(corr);
    out_corr.write(buffer, L*L*sizeof(double));
    buffer = reinterpret_cast<char*>(cov);
    out_cov.write(buffer, 4*L*L*sizeof(complex<double>));
    buffer = reinterpret_cast<char*>(chi);
    out_Jsus.write(buffer, (L-1)*sizeof(double));

    outM.close();
    outJ.close();
    out_corr.close();
    out_cov.close();
    out_Jsus.close();


    //double test1;
    //test1 = two_point_correlation(1, 3, Sz, Sz, Id);


    //char* buffer = reinterpret_cast<char*>(monitor);
    //outMonitor.write(monitor, sizeof(monitor));

    //cout << 4*Nt*sizeof(monitor) << endl;

    /*int bond;
    double N, Mag, J;
    N = Norm(Id);

    char* buffer = reinterpret_cast<char*>(&t);
    outM.write(buffer, sizeof(double));
    outJ.write(buffer, sizeof(double));
    outBond.write(buffer, sizeof(double));
    outL.write(buffer, sizeof(double));
    outL2.write(buffer, sizeof(double));

    for(int loc=1; loc<=L; ++loc){

        Mag= local_observable(loc, Sz, Id)/N;

        buffer = reinterpret_cast<char*>(&Mag);
        outM.write(buffer, sizeof(double));
    }

    for(int loc=1;loc<=L-1;++loc){

        J= (2*two_point_correlation(loc, loc+1, Sx, Sy, Id) - 2*two_point_correlation(loc, loc+1, Sy, Sx, Id))/N;

        buffer = reinterpret_cast<char*>(&J);
        outJ.write(buffer, sizeof(double));
    }

    bond = maxLinkDim(M);
    buffer = reinterpret_cast<char*>(&bond);
    outBond.write(buffer, sizeof(int));

    double Lind = abs(innerC(M, Lindblad, M));
    buffer = reinterpret_cast<char*>(&Lind);
    outL.write(buffer, sizeof(double));

    double Lind2 = innerC(Lindblad, M, Lindblad, M).real(); //- Lind*Lind;
    buffer = reinterpret_cast<char*>(&Lind2);
    outL2.write(buffer, sizeof(double));*/


    //auto Lind =innerC(M, Lindblad, M);
    //outL << t << " " << abs(Lind) << endl;
}


void XXZ_open::Singular_values (){

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
