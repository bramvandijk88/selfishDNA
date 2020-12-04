#include "Core.hh"

Core::Core() : Pearl()
{
	type=0;
}

Pearl* Core::clone() const
{
	return new Core(*this);
}

Core::~Core()	// No default constructor
{
}

Core::Core(int t) : Pearl()
{
	type = t;
	gene_age_=0;
	num_vertical_transfers_=0;
	num_horizontal_transfers_=0;
}

Core::Core(const Core &core) : Pearl(core)
{
	type = core.type;
	mobility = core.mobility;
	gene_age_ = core.gene_age_;
	num_vertical_transfers_ = core.num_vertical_transfers_;
	num_horizontal_transfers_ = core.num_horizontal_transfers_;
}
