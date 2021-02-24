#include "Toxin.hh"

Toxin::Toxin() : Pearl()
{
	type=0;
}

Pearl* Toxin::clone() const
{
	return new Toxin(*this);
}

Toxin::~Toxin()	// No default constructor
{
}

Toxin::Toxin(int t) : Pearl()
{
	type = t;
	gene_age_=0;
	num_vertical_transfers_=0;
	num_horizontal_transfers_=0;
}

Toxin::Toxin(const Toxin &toxin) : Pearl(toxin)
{
	type = toxin.type;
	mobility = toxin.mobility;
	gene_age_ = toxin.gene_age_;
	num_vertical_transfers_ = toxin.num_vertical_transfers_;
	num_horizontal_transfers_ = toxin.num_horizontal_transfers_;
}
