#include "utils.h"

using namespace std;

template<>
void raster(Vect<int, 2>& sub, const Vect<int, 2>& extent)
{
	sub[0] += 1;
	if (sub[0] == extent[0])
	{
		sub[0] = 0;
		sub[1] += 1;
	}
}

template<>
Vect<int,2> raster_end(const Vect<int, 2>& extent)
{
	Vect<int, 2> end = { extent[0] - 1, extent[1] - 1 };
	raster(end, extent);
	return end;
	//return Vect<int, 2>{ 0, extent[1] };
}

template<>
int sub2idx<2>(const Vect<int, 2>& sub, const Vect<int, 2>& extent)
{
	return sub[0] + sub[1]*extent[0];
}
