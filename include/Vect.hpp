#pragma once

#include <array>
#include <cassert>
#include <initializer_list>
#include <iostream>

/** 
 * Thin wrapper for std::array with vector addition, etc.
 */
template<class T, size_t D>
class Vect
{
public:

	template<class S, size_t E> friend class Vect;

	/*Vect() = default;
	Vect(const Vect<T, D>&) = default; // This confuses MSCV since we have a non const copy constructor below.
	Vect(Vect<T, D>&&) = default;
	~Vect() = default;
	Vect<T, D>& operator=(const Vect<T, D>&) = default;
	Vect<T, D>& operator=(Vect<T, D>&&) = default;*/

	Vect()
	{
		for (int i = 0; i < D; ++i)
			components[i] = (T)0;
	}

	Vect(const Vect<T, D>& v)
	{
		components = v.components;
	}

	Vect(const std::array<T, D>& v) {
		for (int i = 0; i < D; ++i) {
			components[i] = v[i];
		}
	}

	~Vect()
	{
	}

	Vect<T, D>& operator=(const Vect<T, D>& v)
	{
		if (&v != this)
		{
			components = v.components;
		}
		return *this;
	}
	
	Vect(std::initializer_list<T> list)
	{
		assert(list.size() == D && "Incorrectly sized initializer_list provided to Vect<T,D>!");
		int i = 0;
		for (auto item : list)
		{
			components[i] = item;
			++i;
		}
	}

	/* 
	// stop the variadic template constructor hiding copy constructor for non-const objects
	Vect(Vect<T, D>& v)
		:Vect<T, D>(const_cast<const Vect<T,D>&>(v))
	{
	}

	// init each component to a different value
	template<typename U, typename... Us>
	Vect(U&& u, Us&&... us) : components{ { std::forward<U>(u), std::forward<Us>(us)... } } {
		static_assert(sizeof...(Us) == D - 1, "Not enough args supplied!");
	}
	*/

	template<class S>
	Vect<S, D> as() const
	{
		Vect<S, D> out;
		for (size_t i = 0; i < D; ++i)
			out.components[i] = (S)components[i];
		return out;
	}

	T& operator[] (size_t i)
	{
		return components[i];
	}

	const T& operator[] (size_t i) const
	{
		return components[i];
	}

	bool operator==(const Vect<T, D>& v) const
	{
		return components == v.components;
	}

	bool operator!=(const Vect<T, D>& v) const
	{
		return !(*this == v);
	}

	Vect<T, D> operator+(const Vect<T, D>& v) const
	{
		Vect<T, D> out;
		for (size_t i = 0; i < D; ++i)
			out.components[i] = components[i] + v.components[i];
		return out;
	}

	Vect<T, D> operator-(const Vect<T, D>& v) const
	{
		Vect<T, D> out;
		for (size_t i = 0; i < D; ++i)
			out.components[i] = components[i] - v.components[i];
		return out;
	}

	Vect<T, D>& operator+=(const Vect<T, D>& v)
	{
		for (size_t i = 0; i < D; ++i)
			components[i] += v.components[i];
		return *this;
	}

	Vect<T, D>& operator-=(const Vect<T, D>& v)
	{
		for (size_t i = 0; i < D; ++i)
			components[i] -= v.components[i];
		return *this;
	}

	Vect<T, D> operator*(T a) const
	{
		Vect<T, D> out;
		for (size_t i = 0; i < D; ++i)
			out.components[i] = a*components[i];
		return out;
	}

	Vect<T, D> operator/(T a) const
	{
		Vect<T, D> out;
		for (size_t i = 0; i < D; ++i)
			out.components[i] = components[i]/a;
		return out;
	}

	Vect<T, D>& operator *= (T a)
	{
		for (size_t i = 0; i < D; ++i)
			components[i] *= a;
		return *this;
	}

	Vect<T, D>& operator /= (T a)
	{
		for (size_t i = 0; i < D; ++i)
			components[i] /= a;
		return *this;
	}

	static Vect<T, D> zero()
	{
		return Vect<T, D>();
	}

	static Vect<T, D> ones()
	{
		Vect<T, D> out;
		for (int i = 0; i < D; ++i)
			out[i] = (T)1;
		return out;
	}

private:
	std::array<T, D> components;
};

template<class T, size_t D>
Vect<T, D> operator*(T a, const Vect<T, D>& v)
{
	return v*a;
}

template<class T, size_t D>
T dot(const Vect<T, D>& a, const Vect<T, D>& b)
{
	T out = (T)0;
	for (size_t i = 0; i < D; ++i)
		out += a[i] * b[i];
	return out;
}

template<class T, size_t D>
T trace(const Vect<T, D>& v)
{
	T out = v[0];
	for (size_t i = 1; i < D; ++i)
		out *= v[i];
	return out;
}

template<class T, size_t D>
std::ostream& operator<<(std::ostream& out, const Vect<T, D>& v)
{
	out << "(";
	for (int i = 0; i < D - 1; ++i)
		out << v[i] << ", ";
	out << v[D - 1] << ")";
	return out;
}