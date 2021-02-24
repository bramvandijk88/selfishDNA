#include "NonEss.hh"
// NonEss = house keeping genes
NonEss::NonEss() : Pearl()
{
	type=0;
}

Pearl* NonEss::clone() const
{
	return new NonEss(*this);
}

NonEss::~NonEss()	// No default constructor
{
}

NonEss::NonEss(int t) : Pearl()
{
	type = t;
	gene_age_=0;
	num_vertical_transfers_=0;
	num_horizontal_transfers_=0;
}

NonEss::NonEss(const NonEss &noness) : Pearl(noness)
{
	type = noness.type;
	mobility = noness.mobility;
	gene_age_ = noness.gene_age_;
	num_vertical_transfers_ = noness.num_vertical_transfers_;
	num_horizontal_transfers_ = noness.num_horizontal_transfers_;
}
