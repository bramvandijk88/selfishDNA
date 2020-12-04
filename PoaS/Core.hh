#ifndef CoreHeader
#define CoreHeader

#include "Pearl.hh"

class Core : public Pearl
{
	public:
		int type;	// Number for coregene.. you need them ALLLLLLL >:)
		double mobility;

		int gene_age_;						// Telt bijv. 1000 tijdstappen
		int num_horizontal_transfers_;		// Elke keer als HGT (reset after x time steps)
		int num_vertical_transfers_;		// Elke keer als VGT

		Core();  // constructor
		virtual Pearl* clone() const;
		~Core(); // deconstructor
		Core(int t); // Constructor that gets toxin type as an argument
		explicit Core(const Core &core); // Copy constructor
};
#endif
