#ifndef ARGHeader
#define ARGHeader

#include "Pearl.hh"

class ARG : public Pearl
{
	public:
		int type;	                       // Number ranging from 0 to total number of existing resistance genes
		double mobility;

		int gene_age_;						
		int num_horizontal_transfers_;		// HGT counter 
		int num_vertical_transfers_;		// VGT counter

		ARG();  // constructor
		virtual Pearl* clone() const;
		~ARG(); // deconstructor
		ARG(int t); // Constructor that gets ARG type as an argument
		explicit ARG(const ARG &ARG); // Copy constructor
};

#endif
