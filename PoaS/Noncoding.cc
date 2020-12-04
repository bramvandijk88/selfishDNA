#include "Noncoding.hh"

Noncoding::Noncoding() : Pearl()
{
	type=0;
}

Pearl* Noncoding::clone() const
{
	return new Noncoding(*this);
}

Noncoding::~Noncoding()	// No default constructor
{
}

Noncoding::Noncoding(int t) : Pearl()
{
	type = t;
}

Noncoding::Noncoding(const Noncoding &nc) : Pearl(nc)
{
	type = nc.type;
	mobility = nc.mobility;
}
