#include "Antitoxin.hh"

Antitoxin::Antitoxin() : Pearl() 	// No default constructor DO NOT USE
{
	type=0;
}

Pearl* Antitoxin::clone() const
{
	return new Antitoxin(*this);
}

Antitoxin::~Antitoxin()
{
}

Antitoxin::Antitoxin(int t) : Pearl()
{
	type = t;
	gene_age_=0;
	num_vertical_transfers_=0;
	num_horizontal_transfers_=0;
}

Antitoxin::Antitoxin(const Antitoxin &antitoxin) : Pearl(antitoxin)
{
	type = antitoxin.type;
	mobility = antitoxin.mobility;
	gene_age_ = antitoxin.gene_age_;
	num_vertical_transfers_ = antitoxin.num_vertical_transfers_;
	num_horizontal_transfers_ = antitoxin.num_horizontal_transfers_;
}
