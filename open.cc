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


double open::Norm(ITensor &Id){

   auto Aux = M(1);

    for(int loc=1; loc<=L; ++loc){

        Id *= delta( inds(Id)(1), findIndex(M(loc),"Site") );

        Aux *= Id;
        if(loc != L)
            Aux *= M(loc+1);
    }

    return elt(realPart(Aux));
}


double open::local_observable(int loc, ITensor &O, ITensor &Id){

    auto Aux = M(1);

    for(int l=1; l<=L; ++l){

        if( l != loc ){
            Id *= delta( inds(Id)(1), findIndex(M(l),"Site") );
            Aux *= Id;
        }
        else{
            O *= delta( inds(O)(1), findIndex(M(l),"Site") );
            Aux *= O;
        }

        if(l != L)
            Aux *= M(l+1);
    }

    return elt(realPart(Aux));
}


double open::two_point_correlation(int loc1, int loc2, ITensor &O1, ITensor &O2, ITensor &Id){

    auto Aux = M(1);

    for(int l=1; l<=L; ++l){

        if( l == loc1 ){
            O1 *= delta( inds(O1)(1), findIndex(M(l),"Site") );
            Aux *= O1;
        }
        else if( l == loc2 ){
            O2 *= delta( inds(O2)(1), findIndex(M(l),"Site") );
            Aux *= O2;
        }
        else{
            Id *= delta( inds(Id)(1), findIndex(M(l),"Site") );
            Aux *= Id;
        }

        if(l != L)
            Aux *= M(l+1);
    }

    return elt(realPart(Aux));
}


complex<double> open::two_point_correlation_string(int loc1, int loc2, ITensor &O1, ITensor &O2, ITensor &Id, ITensor &Sz){

    auto Aux = M(1);

    for(int l=1; l<=L; ++l){

        if( l == loc1 ){
            O1 *= delta( inds(O1)(1), findIndex(M(l),"Site") );
            Aux *= O1;
        }
        else if( l == loc2 ){
            O2 *= delta( inds(O2)(1), findIndex(M(l),"Site") );
            Aux *= O2;
        }
        else if( l > loc1 && l < loc2 ){
            Sz *= delta( inds(Sz)(1), findIndex(M(l),"Site") );
            Aux *= Sz;
        }
        else{
            Id *= delta( inds(Id)(1), findIndex(M(l),"Site") );
            Aux *= Id;
        }

        if(l != L)
            Aux *= M(l+1);
    }

    return eltC(Aux);
}


using namespace Eigen;


MatrixXcd Hamiltonian_site(int d, int size, int loc, double parameter, MatrixXcd &O, MatrixXcd &Id){

    MatrixXcd Aux, Iblock, res(1,1);
    res(0,0)=1;

    for(int i=0; i<loc%size; ++i)
        res = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, res, Id);
    Aux = res;
    Iblock = res;

    res = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, res, O);
    Aux = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, Aux, O.transpose());
    Iblock = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, Iblock, Id);

    for(int i=loc%size+1; i< size; ++i){
        res = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, res, Id);
        Aux = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, Aux, Id);
        Iblock = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, Iblock, Id);
    }

    res = KroneckerProduct(pow(2, size), pow(2, size), pow(2, size), pow(2, size), res, Iblock);
    Aux = KroneckerProduct(pow(2, size), pow(2, size), pow(2, size), pow(2, size), Iblock, Aux);
    Iblock = KroneckerProduct(pow(2, size), pow(2, size), pow(2, size), pow(2, size), Iblock, Iblock);

    if(loc/size==0)
        res = -KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), res, Iblock) + KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), Aux, Iblock);
    else
        res = -KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), Iblock, res) + KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), Iblock, Aux);

    return complex<double>(0,1) * parameter * res;
}


void Hamiltonian_bond_step(int size, int loc, int &dim, int &dima, MatrixXcd &O, MatrixXcd &res, MatrixXcd &Aux, MatrixXcd &Iblock){

    res = KroneckerProduct(dim, dim, 2, 2, res, O);
    dim = 2 * dim;

    if((loc+1)%size==0){
        res = KroneckerProduct(dim, dim, pow(2, size), pow(2, size), res, Iblock);
        dim = pow(2, size) * dim;
    }

    if(loc%size==0){
        Aux = KroneckerProduct(dima, dima, pow(2, size), pow(2, size), Aux, Iblock);
        dima = pow(2, size) * dima;
    }

    Aux = KroneckerProduct(dima, dima, 2, 2, Aux, O.transpose());
    dima = 2 * dima;
}


MatrixXcd Hamiltonian_bond(int d, int size, int loc1, int loc2, double parameter, MatrixXcd &O1, MatrixXcd &O2, MatrixXcd &Id){

    int dim=1, dima=1;
    MatrixXcd Aux(1,1), Iblock(1,1), res(1,1);
    res(0,0)=1, Aux(0,0)=1, Iblock(0,0)=1;

    for(int i=0; i<=size; ++i)
        Iblock = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, Iblock, Id);

    for(int i=0; i<loc1; ++i)
        Hamiltonian_bond_step(size, i, dim, dima, Id, res, Aux, Iblock);

    Hamiltonian_bond_step(size, loc1, dim, dima, O1, res, Aux, Iblock);

    for(int i=loc1+1; i<loc2; ++i)
        Hamiltonian_bond_step(size, i, dim, dima, Id, res, Aux, Iblock);

    Hamiltonian_bond_step(size, loc2, dim, dima, O2, res, Aux, Iblock);

    for(int i=loc2+1; i<2*size; ++i)
        Hamiltonian_bond_step(size, i, dim, dima, Id, res, Aux, Iblock);

    return complex<double>(0,1) * parameter * (-res + Aux);
}


MatrixXcd Hamiltonian_fermion_bond(int size, int loc1, int loc2, double parameter, MatrixXcd &O1, MatrixXcd &O2, MatrixXcd &Id, MatrixXcd &Sz, int string){

    int dim=1, dima=1;
    MatrixXcd Block(1,1), Iblock(1,1), res, mirror, aux1, aux2;
    Block(0,0)=1, Iblock(0,0)=1;
    MatrixXcd Op[2*size];

    if(string==0){
        for(int i=0; i<2*size; ++i)
            Op[i]= Id;
        Op[loc1]= O1;
        Op[loc2]= O2;
    }
    else{
        for(int i=0; i<2*size; ++i)
            Op[i]= Sz;
        Op[loc1]= O1;
        Op[loc2]= O2;
    }

    for(int i=0; i<size; ++i)
        Iblock = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, Iblock, Id);

    aux1= KroneckerProduct(2, 2, 2, 2, Op[0], Op[1]);
    aux2= KroneckerProduct(2, 2, 2, 2, Op[2], Op[3]);

    res= KroneckerProduct(4, 4, 4, 4, aux1, Iblock);
    res= KroneckerProduct(16, 16, 4, 4, res, aux2);
    res= KroneckerProduct(64, 64, 4, 4, res, Iblock);

    mirror= KroneckerProduct(4, 4, 4, 4, Iblock, aux1.transpose());
    mirror= KroneckerProduct(16, 16, 4, 4, mirror, Iblock);
    mirror= KroneckerProduct(64, 64, 4, 4, mirror, aux2.transpose());

    return complex<double>(0,1) * parameter * (-res + mirror);
}


MatrixXcd Jump_site(int d, int size, int loc, double parameter, MatrixXcd &O, MatrixXcd &Id){

    MatrixXcd Aux, Iblock, left, right, res(1,1);
    res(0,0)=1;

    for(int i=0; i<loc%size; ++i)
        res = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, res, Id);
    Aux = res;
    left = res;
    right = res;
    Iblock = res;

    res = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, res, O);
    Aux = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, Aux, O.adjoint().transpose());
    left = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, left, O.adjoint()*O);
    right = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, right, (O.adjoint()*O).transpose());
    Iblock = KroneckerProduct(pow(2, loc%size), pow(2, loc%size), 2, 2, Iblock, Id);

    for(int i=loc%size+1; i< size; ++i){
        res = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, res, Id);
        Aux = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, Aux, Id);
        left = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, left, Id);
        right = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, right, Id);
        Iblock = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, Iblock, Id);
    }

    res = KroneckerProduct(pow(2, size), pow(2, size), pow(2, size), pow(2, size), res, Aux);
    left = KroneckerProduct(pow(2, size), pow(2, size), pow(2, size), pow(2, size), left, Iblock);
    right = KroneckerProduct(pow(2, size), pow(2, size), pow(2, size), pow(2, size), Iblock, right);
    Iblock = KroneckerProduct(pow(2, size), pow(2, size), pow(2, size), pow(2, size), Iblock, Iblock);

    if(loc/size==0){
        res = KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), res, Iblock);
        left = KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), left, Iblock);
        right = KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), right, Iblock);
    }
    else{
        res = KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), Iblock, res);
        left = KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), Iblock, left);
        right = KroneckerProduct(pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), pow(2, 2*size), Iblock, right);
    }

    return parameter * (res -0.5*left -0.5*right);
}


MatrixXcd Jump_fermion_site(int size, int loc1, int loc2, double parameter, MatrixXcd &O, MatrixXcd &Id, MatrixXcd &Sz){

    MatrixXcd Aux, Iblock(1,1), left, right, res(1,1);
    MatrixXcd Op[2*size];
    Iblock(0,0)=1;

    for(int i=0; i<2*size; ++i)
        Op[i]= Id;
    Op[loc1%2]= O;
    Op[loc2%2]= Sz;

    for(int i=0; i<2*size; ++i)
        Iblock = KroneckerProduct(pow(2, i), pow(2, i), 2, 2, Iblock, Id);

    res = KroneckerProduct(2, 2, 2, 2, Op[0], Op[1]);
    res = KroneckerProduct(4, 4, 2, 2, res, Op[0].adjoint().transpose());
    res = KroneckerProduct(8, 8, 2, 2, res, Op[1].adjoint().transpose());

    left = KroneckerProduct(2, 2, 2, 2, Op[0].adjoint()*Op[0], Op[1].adjoint()*Op[1]);
    left = KroneckerProduct(4, 4, 2, 2, left, Id);
    left = KroneckerProduct(8, 8, 2, 2, left, Id);

    right = KroneckerProduct(2, 2, 2, 2, Id, Id);
    right = KroneckerProduct(4, 4, 2, 2, right, (Op[0].adjoint()*Op[0]).transpose());
    right = KroneckerProduct(8, 8, 2, 2, right, (Op[1].adjoint()*Op[1]).transpose());

    if(loc1/2==0){
        res = KroneckerProduct(16, 16, 16, 16, res, Iblock);
        left = KroneckerProduct(16, 16, 16, 16, left, Iblock);
        right = KroneckerProduct(16, 16, 16, 16, right, Iblock);
    }
    else{
        res = KroneckerProduct(16, 16, 16, 16, Iblock, res);
        left = KroneckerProduct(16, 16, 16, 16, Iblock, left);
        right = KroneckerProduct(16, 16, 16, 16, Iblock, right);
    }

    return parameter * (res -0.5*left -0.5*right);
}


void Lindblad_site(int d, int B, int size, int loc, MatrixXcd &LBulk, MatrixXcd &LLeft, MatrixXcd &LRight, double parameter, MatrixXcd &O, MatrixXcd &Id){

    int dim = pow(d, 2*size);
    MatrixXcd Aux, res(1,1), Iblock = MatrixXd::Identity(pow(d,size), pow(d,size));;
    res(0,0)=1;

    for(int i=0; i<loc; ++i)
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
    Aux = res;

    res = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, res, O);
    Aux = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, Aux, O.transpose());

    for(int i=loc+1; i<size; ++i){
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
        Aux = KroneckerProduct(pow(d, i), pow(d, i), d, d, Aux, Id);
    }
    res = complex<double>(0,1) * parameter * (-KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), res, Iblock) + KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), Iblock, Aux));

    LBulk.block((B-1)*dim, 0, dim, dim) += res;
    LLeft.block(0, 0, dim, dim) += res;
    LRight.block((B-1)*dim, 0, dim, dim) += res;
}


void Lindblad_bond(int d, int B, int size, int loc, int pos, MatrixXcd &LBulk, MatrixXcd &LLeft, MatrixXcd &LRight, double parameter, MatrixXcd &O, MatrixXcd &Id){

    int dim = pow(d, 2*size);
    MatrixXcd Aux, res(1,1), Iblock = MatrixXd::Identity(pow(d,size), pow(d,size));;
    res(0,0)=1;

    for(int i=0; i<loc; ++i)
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
    Aux = res;

    res = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, res, O);
    Aux = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, Aux, O.transpose());

    for(int i=loc+1; i<size; ++i){
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
        Aux = KroneckerProduct(pow(d, i), pow(d, i), d, d, Aux, Id);
    }
    res = KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), res, Iblock);
    Aux = KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), Iblock, Aux);

    //Bulk
    LBulk.block(pos*dim, 0, dim, dim) += res;
    LBulk.block((pos+1)*dim, 0, dim, dim) += Aux;
    LBulk.block((B-1)*dim, pos*dim, dim, dim) += -complex<double>(0,1) * parameter * res;
    LBulk.block((B-1)*dim, (pos+1)*dim, dim, dim) += complex<double>(0,1) * parameter * Aux;

    //Left
    LLeft.block(0, pos*dim, dim, dim) += -complex<double>(0,1) * parameter * res;
    LLeft.block(0, (pos+1)*dim, dim, dim) += complex<double>(0,1) * parameter * Aux;

    //Right
    LRight.block(pos*dim, 0, dim, dim) += res;
    LRight.block((pos+1)*dim, 0, dim, dim) += Aux;
}


void Lindblad_rung(int d, int B, int size, int loc1, int loc2, MatrixXcd &LBulk, MatrixXcd &LLeft, MatrixXcd &LRight, double parameter, MatrixXcd &O1, MatrixXcd &O2, MatrixXcd &Id){

    int dim = pow(d, 2*size);
    MatrixXcd Aux, res(1,1), Iblock = MatrixXd::Identity(pow(d,size), pow(d,size));;
    res(0,0)=1;

    for(int i=0; i<loc1; ++i)
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
    Aux = res;

    res = KroneckerProduct(pow(d, loc1), pow(d, loc1), d, d, res, O1);
    Aux = KroneckerProduct(pow(d, loc1), pow(d, loc1), d, d, Aux, O1.transpose());

    for(int i=loc1+1; i<loc2; ++i){
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
        Aux = KroneckerProduct(pow(d, i), pow(d, i), d, d, Aux, Id);
    }

    res = KroneckerProduct(pow(d, loc2), pow(d, loc2), d, d, res, O2);
    Aux = KroneckerProduct(pow(d, loc2), pow(d, loc2), d, d, Aux, O2.transpose());

    for(int i=loc2+1; i<size; ++i){
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
        Aux = KroneckerProduct(pow(d, i), pow(d, i), d, d, Aux, Id);
    }
    res = complex<double>(0,1) * parameter * (-KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), res, Iblock) + KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), Iblock, Aux));

    LBulk.block((B-1)*dim, 0, dim, dim) += res;
    LLeft.block(0, 0, dim, dim) += res;
    LRight.block((B-1)*dim, 0, dim, dim) += res;
}


void Lindblad_jump(int d, int B, int size, int loc, MatrixXcd &LBulk, MatrixXcd &LLeft, MatrixXcd &LRight, double parameter, MatrixXcd &O, MatrixXcd &Id, int ctrl){

    int dim = pow(d, 2*size);
    MatrixXcd Aux, left, right, res(1,1), Iblock = MatrixXd::Identity(pow(d,size), pow(d,size));;
    res(0,0)=1;

    for(int i=0; i<loc; ++i)
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
    Aux = res;
    left = res;
    right = res;

    res = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, res, O);
    Aux = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, Aux, O.adjoint().transpose());
    left = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, left, O.adjoint()*O);
    right = KroneckerProduct(pow(d, loc), pow(d, loc), d, d, right, (O.adjoint()*O).transpose());

    for(int i=loc+1; i<size; ++i){
        res = KroneckerProduct(pow(d, i), pow(d, i), d, d, res, Id);
        Aux = KroneckerProduct(pow(d, i), pow(d, i), d, d, Aux, Id);
        left = KroneckerProduct(pow(d, i), pow(d, i), d, d, left, Id);
        right = KroneckerProduct(pow(d, i), pow(d, i), d, d, right, Id);
    }

    res = KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), res, Aux);
    left = KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), left, Iblock);
    right = KroneckerProduct(pow(d, size), pow(d, size), pow(d, size), pow(d, size), Iblock, right);
    res = parameter * (res -0.5*left -0.5*right);

    if(ctrl==0){
        LBulk.block((B-1)*dim, 0, dim, dim) += res;
        LLeft.block(0, 0, dim, dim) += res;
        LRight.block((B-1)*dim, 0, dim, dim) += res;
    }
    else if(ctrl==1)
        LLeft.block(0, 0, dim, dim) += res;
    else
        LRight.block((B-1)*dim, 0, dim, dim) += res;
}
