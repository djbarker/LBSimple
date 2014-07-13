#pragma once

#include "Vect.hpp"

enum ModelType
{
	D2Q9
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