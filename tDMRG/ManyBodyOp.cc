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
#include "ManyBodyOp.h"

using namespace std;
using namespace Eigen;
using namespace itensor;


//Functions to give dimension of matrix
int power (int a, int b){

    int res=1;
    for(int i=0;i<b;++i)
        res= res*a;

    return res;
}

int min (int a, int b){

    if(a<b)
        return a;
    else
        return b;
}

int dim(int i, int L, int d, int D){

    int res= min(i,L-i);

    if( res <= log(D)/log(d) )
        res= power(d,res);
    else
        res= D;

    return res;
}


//Function to multiply complex numbers
/*lapack_complex_double Complex_mul (lapack_complex_double a, lapack_complex_double b) {

    lapack_complex_double res;

    res.real= a.real * b.real - a.imag * b.imag;
    res.imag= a.real * b.imag + a.imag * b.real;

    return res;
}*/


void Pauli_Matrices(int d, MatrixXcd &Sx, MatrixXcd &Sy, MatrixXcd &Sz, MatrixXcd &Su, MatrixXcd &Sd, MatrixXcd &Id){

    Sx(0,0) = complex<double>(0,0); Sx(0,1) = complex<double>(1,0); Sx(1,0) = complex<double>(1,0); Sx(1,1) = complex<double>(0,0);
    Sy(0,0) = complex<double>(0,0); Sy(0,1) = complex<double>(0,-1); Sy(1,0) = complex<double>(0,1); Sy(1,1) = complex<double>(0,0);
    Sz(0,0) = complex<double>(1,0); Sz(0,1) = complex<double>(0,0); Sz(1,0) = complex<double>(0,0); Sz(1,1) = complex<double>(-1,0);
    Su(0,0) = complex<double>(0,0); Su(0,1) = complex<double>(1,0); Su(1,0) = complex<double>(0,0); Su(1,1) = complex<double>(0,0);
    Sd(0,0) = complex<double>(0,0); Sd(0,1) = complex<double>(0,0); Sd(1,0) = complex<double>(1,0); Sd(1,1) = complex<double>(0,0);
    Id(0,0) = complex<double>(1,0); Id(0,1) = complex<double>(0,0); Id(1,0) = complex<double>(0,0); Id(1,1) = complex<double>(1,0);
}


MatrixXcd KroneckerProduct(int m, int n, int p, int q, const Ref<const MatrixXcd>&O1, const Ref<const MatrixXcd>&O2){

    MatrixXcd res(m*p, n*q);

    for(int i=0; i<m*p; ++i){
        for(int j=0; j<q*n; ++j)
            res(i,j) = O1(i/p, j/q) * O2(i%p, j%q);
    }

    return res;
}


MatrixXcd matrix_complex_exponential( int d, double dt,  MatrixXcd &Lind){

    MatrixXcd expo(d,1);
    ComplexEigenSolver<MatrixXcd> dia;
    dia.compute(Lind);

    for(int i=0;i<d;++i)
        expo(i,0) = exp(dia.eigenvalues()(i,0) * dt);

    return dia.eigenvectors() * expo.asDiagonal() * dia.eigenvectors().inverse();
}


//Converts Eigen matrices to itensor
void TEvol_Eigen_itensor(int d, MatrixXcd &Lind, ITensor *U){

    for(int i=0; i<d; ++i){
        for(int j=0; j<d; ++j){
            for(int k=0; k<d; ++k){
                for(int l=0; l<d; ++l)
                    U->set(i+1, j+1, k+1, l+1, Lind(i*d + j, k*d + l).real() + Lind(i*d + j, k*d + l).imag() *1_i );
            }
        }
    }
}

void LBulk_Eigen_itensor(int d, int B, MatrixXcd &Lind, ITensor *U){

    for(int i=0; i<B; ++i){
        for(int j=0; j<d; ++j){
            for(int k=0; k<B; ++k){
                for(int l=0; l<d; ++l)
                    U->set(j+1, l+1, i+1, k+1, Lind(i*d + j, k*d + l).real() + Lind(i*d + j, k*d + l).imag() *1_i );
            }
        }
    }
}

void LLeft_Eigen_itensor(int d, int B, MatrixXcd &Lind, ITensor *U){

    for(int i=0; i<d; ++i){
        for(int j=0; j<B; ++j){
            for(int k=0; k<d; ++k)
                U->set(i+1, k+1, j+1, Lind(i, j*d + k).real() + Lind(i, j*d + k).imag() *1_i );
        }
    }
}

void LRight_Eigen_itensor(int d, int B, MatrixXcd &Lind, ITensor *U){

    for(int i=0; i<B; ++i){
        for(int j=0; j<d; ++j){
            for(int k=0; k<d; ++k)
                U->set(j+1, k+1, i+1, Lind(i*d + j, k).real() + Lind(i*d + j, k).imag() *1_i );
        }
    }
}
