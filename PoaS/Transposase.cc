#include "Transposase.hh"

Transposase::Transposase() : Pearl()
{
	type=0;
}

Pearl* Transposase::clone() const
{
	return new Transposase(*this);
}

Transposase::~Transposase()	// No default constructor
{
}

Transposase::Transposase(int t) : Pearl()
{
	type = t;
	gene_age_=0;
	num_vertical_transfers_=0;
	num_horizontal_transfers_=0;
}

Transposase::Transposase(const Transposase &transposase) : Pearl(transposase)
{
	type = transposase.type;
	mobility = transposase.mobility;
	gene_age_ = transposase.gene_age_;
	num_vertical_transfers_ = transposase.num_vertical_transfers_;
	num_horizontal_transfers_ = transposase.num_horizontal_transfers_;
}
