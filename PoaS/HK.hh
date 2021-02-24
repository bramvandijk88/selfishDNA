#ifndef HKHeader
#define HKHeader

#include "Pearl.hh"

class HK : public Pearl
{
	public:
		int type;	// Number for HKgene.. you need them ALLLLLLL >:)
		double mobility;

		// int gene_age_;						// Time since this gene was born from a non-coding gene		
		// int num_horizontal_transfer_;		// HGT counter (transposon jumps between cells)
		// int num_vertical_transfers_;		// VGT counter
		// int num_jumps_;						// Transposon jump counter (transposon jumps in a cell)

		HK();  // constructor
		virtual Pearl* clone() const;
		~HK(); // deconstructor
		HK(int t); // Constructor that gets toxin type as an argument
		explicit HK(const HK &hk); // Copy constructor
};
#endif
