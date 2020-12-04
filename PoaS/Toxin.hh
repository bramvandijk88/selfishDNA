#ifndef ToxinHeader
#define ToxinHeader

#include "Pearl.hh"

class Toxin : public Pearl
{
	public:
		int type;	// Number ranging from 0 to NR_POSSIBLE_TOXINS which defines the produced toxin
		double mobility;

		int gene_age_;						// Telt bijv. 1000 tijdstappen
		int num_horizontal_transfers_;		// Elke keer als HGT (reset after x time steps)
		int num_vertical_transfers_;		// Elke keer als VGT

		Toxin();  // constructor
		virtual Pearl* clone() const;
		~Toxin(); // deconstructor
		Toxin(int t); // Constructor that gets toxin type as an argument
		explicit Toxin(const Toxin &toxin); // Copy constructor

};
#endif
