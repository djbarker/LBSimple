#include "utils.hpp"

using namespace std;

template<>
int sub2idx<2>(const Vect<int, 2>& sub, const Vect<int, 2>& extent)
{
	return sub[0] + sub[1]*extent[0];
}

template<>
int sub2idx<3>(const Vect<int, 3>& sub, const Vect<int, 3>& extent)
{
	return sub[0] + sub[1] * extent[1] + sub[2] * extent[0] * extent[1];
}


template<>
Vect<int, 2> idx2sub<2>(size_t idx, const Vect<int, 2>& extent)
{
	Vect<int, 2> out;
	out[0] = idx % extent[0];
	out[1] = (idx / extent[0]) % extent[1];
	return out;
}

template<>
Vect<int, 3> idx2sub<3>(size_t idx, const Vect<int, 3>& extent)
{
	Vect<int, 3> out;
	out[0] = idx % extent[0];
	out[1] = (idx / extent[0]) % extent[1];
	out[2] = ((idx / extent[0]) / extent[1]) % extent[2];
	return out;
}

template<>
Vect<int, 2> calc_num_domains<2>(size_t nproc)
{
	int root = 1;
	for (int i = 2; i <= sqrt(nproc); ++i)
	if (nproc%i == 0) root = i;

	return Vect<int, 2>{root, (int)nproc / root};
}

template<>
Vect<int, 3> calc_num_domains<3>(size_t nproc)
{
	int root1 = 1;
	for (int i = 2; i <= pow(nproc, 1. / 3.); ++i)
	if (nproc%i == 0) root1 = i;

	int root2 = 1;
	for (int i = 2; i <= sqrt(nproc / root1); ++i)
	if ((nproc / root1) % i == 0) root2 = i;

	return Vect<int, 3>{root1, root2, (int)nproc / (root1*root2)};
}