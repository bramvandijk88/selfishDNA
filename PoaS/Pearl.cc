#include "Pearl.hh"

/*
All constructors and deconstructors are empty because this is just a baseclass for all kinds of Pearls such 
as genes, transcription factors, transposons, or anything you can think of to add to the Pearls-on-a-string
representation of a genome

Currently, it includes:
Toxin_Antitoxin.hh
*/
Pearl::Pearl()
{
	gene_age_=0;
	num_vertical_transfers_=0;
    num_jumps_=0;
	num_horizontal_transfers_=0;
}

Pearl::~Pearl()
{
}

Pearl::Pearl(const Pearl &bd)
{
    type = bd.type;
	gene_age_=bd.gene_age_;
    num_jumps_=bd.num_jumps_;
	num_vertical_transfers_=bd.num_vertical_transfers_;
	num_horizontal_transfers_=bd.num_horizontal_transfers_;
}
