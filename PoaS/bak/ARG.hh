#ifndef ARGHeader
#define ARGHeader

#include "Pearl.hh"

class ARG : public Pearl
{
	public:
		int type;	                       // Number ranging from 0 to total number of existing resistance genes
		double mobility;

		// int gene_age_;						// Time since this gene was born from a non-coding gene		
		// int num_horizontal_transfer_;		// HGT counter (transposon jumps between cells)
		// int num_vertical_transfers_;		// VGT counter
		// int num_jumps_;						// Transposon jump counter (transposon jumps in a cell)

		ARG();  // constructor
		virtual Pearl* clone() const;
		~ARG(); // deconstructor
		ARG(int t); // Constructor that gets ARG type as an argument
		explicit ARG(const ARG &ARG); // Copy constructor
};

#endif
