#include "Model.hpp"

const Vect<int,2> Model<D2Q9>::Es[9]   = { Vect<int,2>{ 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } };
const int         Model<D2Q9>::Qneg[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
const double      Model<D2Q9>::Ws[9]   = { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };