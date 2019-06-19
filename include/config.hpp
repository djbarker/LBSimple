#pragma once

// compile time config

#include <include/Vect.hpp>
#include <include/Model.hpp>

// select the model
static constexpr ModelType the_model_type = ModelType::D2Q9;
using the_model = Model<the_model_type>;
static constexpr size_t Dims = the_model::Dim;
static constexpr size_t Q = the_model::Q;

// useful typedefs
typedef unsigned char RawType;

typedef Vect<double, Dims> vect_t;
typedef Vect<int, Dims> sub_t;