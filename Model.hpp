#pragma once

#include "Vect.hpp"

enum ModelType
{
	D2Q9,
	D3Q19
};

template<ModelType M> struct Model;

template<> struct Model<D2Q9> 
{
	static const int Dim = 2;
	static const int Q = 9;
	static const Vect<int, 2> Es[9];
	static const int Qneg[9];
	static const double Ws[9];
};

template<> struct Model<D3Q19>
{
	static const int Dim = 3;
	static const int Q = 19;
	static const Vect<int, 3> Es[19];
	static const int Qneg[19];
	static const double Ws[19];
};