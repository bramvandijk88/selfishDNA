#ifndef AntitoxinHeader
#define AntitoxinHeader

#include "Pearl.hh"

class Antitoxin : public Pearl
{
	public:
		int type;	// Number ranging from 0 to NR_POSSIBLE_TOXINS which defines the produced antitoxin
		double mobility;

///////////////////////////////

		Antitoxin();  // constructor
		virtual Pearl* clone() const;
		~Antitoxin(); // deconstructor
		Antitoxin(int t); // Constructor that gets antitoxin type as an argument
		explicit Antitoxin(const Antitoxin &antitoxin); // Copy constructor
};

#endif
