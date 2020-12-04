#include "ARG.hh"

ARG::ARG() : Pearl() 	// No default constructor DO NOT USE
{
	type=0;
}

Pearl* ARG::clone() const
{
	return new ARG(*this);
}

ARG::~ARG()
{
}

ARG::ARG(int t) : Pearl()
{
	type = t;
	gene_age_=0;
	num_vertical_transfers_=0;
	num_horizontal_transfers_=0;
}

ARG::ARG(const ARG &ARG) : Pearl(ARG)
{
	type = ARG.type;
	mobility = ARG.mobility;
	gene_age_ = ARG.gene_age_;
	num_vertical_transfers_ = ARG.num_vertical_transfers_;
	num_horizontal_transfers_ = ARG.num_horizontal_transfers_;
}
