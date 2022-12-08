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


#ifndef OPEN_H
#define OPEN_H


class open: public myMPS {

public:
    MPO Lindblad;

    open (int Lin, int din, int Din) : myMPS(Lin, din, Din) {}

    double Norm(ITensor &);

    double local_observable(int, ITensor &, ITensor &);

    double two_point_correlation(int, int, ITensor &, ITensor &, ITensor &);

    complex<double> two_point_correlation_string(int, int, ITensor &, ITensor &, ITensor &, ITensor &);
};


Eigen::MatrixXcd Hamiltonian_site(int, int, int, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &Id);

Eigen::MatrixXcd Hamiltonian_bond(int, int, int, int, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &);

Eigen::MatrixXcd Hamiltonian_fermion_bond(int, int, int, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, int);

Eigen::MatrixXcd Jump_site(int, int, int, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &);

Eigen::MatrixXcd Jump_fermion_site(int, int, int, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &);

void Lindblad_site(int, int, int, int, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &);

void Lindblad_bond(int, int, int, int, int, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &);

void Lindblad_rung(int, int, int, int, int, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &);

void Lindblad_jump(int, int, int, int, Eigen::MatrixXcd &, Eigen::MatrixXcd &, Eigen::MatrixXcd &, double, Eigen::MatrixXcd &, Eigen::MatrixXcd &, int);


#endif
