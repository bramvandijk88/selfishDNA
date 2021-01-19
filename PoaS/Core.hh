#ifndef CoreHeader
#define CoreHeader

#include "Pearl.hh"

class Core : public Pearl
{
	public:
		int type;	// Number for coregene.. you need them ALLLLLLL >:)
		double mobility;

		// int gene_age_;						// Time since this gene was born from a non-coding gene		
		// int num_horizontal_transfer_;		// HGT counter (transposon jumps between cells)
		// int num_vertical_transfers_;		// VGT counter
		// int num_jumps_;						// Transposon jump counter (transposon jumps in a cell)

		Core();  // constructor
		virtual Pearl* clone() const;
		~Core(); // deconstructor
		Core(int t); // Constructor that gets toxin type as an argument
		explicit Core(const Core &core); // Copy constructor
};
#endif
