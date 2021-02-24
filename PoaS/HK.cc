#include "HK.hh"
// HK = house keeping genes
HK::HK() : Pearl()
{
	type=0;
}

Pearl* HK::clone() const
{
	return new HK(*this);
}

HK::~HK()	// No default constructor
{
}

HK::HK(int t) : Pearl()
{
	type = t;
	gene_age_=0;
	num_vertical_transfers_=0;
	num_horizontal_transfers_=0;
}

HK::HK(const HK &hk) : Pearl(hk)
{
	type = hk.type;
	mobility = hk.mobility;
	gene_age_ = hk.gene_age_;
	num_vertical_transfers_ = hk.num_vertical_transfers_;
	num_horizontal_transfers_ = hk.num_horizontal_transfers_;
}
