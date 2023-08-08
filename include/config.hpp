#pragma once

// compile time config

#include <include/Model.hpp>
#include <include/Vect.hpp>

// select the model
static constexpr ModelType the_model_type = ModelType::D2Q9;
using the_model = Model<the_model_type>;
static constexpr size_t Dims = the_model::Dim;
static constexpr size_t Q = the_model::Q;

// useful typedefs
typedef unsigned char RawType; // TODO: RawType -> cell_t

typedef double scalar_t;
typedef Vect<scalar_t, Dims> vect_t;
typedef Vect<int, Dims> sub_t; // NOTE: int not size_t since we can have negative numbers before periodic wrapping
