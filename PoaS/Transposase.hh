#ifndef TransposaseHeader
#define TransposaseHeader

#include "Pearl.hh"

class Transposase : public Pearl
{
	public:
		int type;	
		double mobility;

		int gene_age_;						// Telt bijv. 1000 tijdstappen
		int num_horizontal_transfers_;		// Elke keer als HGT (reset after x time steps)
		int num_vertical_transfers_;		// Elke keer als VGT

		Transposase();  // constructor
		virtual Pearl* clone() const;
		~Transposase(); // deconstructor
		Transposase(int t); // Constructor that gets toxin type as an argument
		explicit Transposase(const Transposase &transposon); // Copy constructor
};
#endif