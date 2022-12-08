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
//#include <mkl.h>
#include "../Eigen/Eigen/Eigen"
#include "itensor/all.h"


#ifndef MANYBODYOP_H
#define MANYBODYOP_H


int dim(int i, int L, int d, int D);

//lapack_complex_double Complex_mul (lapack_complex_double a, lapack_complex_double b);

void Pauli_Matrices(int, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &);

Eigen::MatrixXcd KroneckerProduct(int, int, int, int, const Eigen::Ref<const Eigen::MatrixXcd>&, const Eigen::Ref<const Eigen::MatrixXcd>&);

Eigen::MatrixXcd matrix_complex_exponential( int, double,  Eigen::MatrixXcd &);

void TEvol_Eigen_itensor(int, Eigen::MatrixXcd &, itensor::ITensor *);


void LBulk_Eigen_itensor(int, int, Eigen::MatrixXcd &, itensor::ITensor *);


void LLeft_Eigen_itensor(int, int, Eigen::MatrixXcd &, itensor::ITensor *);


void LRight_Eigen_itensor(int, int, Eigen::MatrixXcd &, itensor::ITensor *);



/*void Update_Lind(int, double, lapack_complex_double*, std::ifstream*);

void Copy_temp(int, lapack_complex_double*, Eigen::MatrixXcd &);

void Copy_MPO_Bulk(int, int, lapack_complex_double*, itensor::ITensor*);

void Copy_MPO_Left(int, int, lapack_complex_double*, itensor::ITensor*);

void Copy_MPO_Right(int, int, lapack_complex_double*, itensor::ITensor*);
*/

#endif
