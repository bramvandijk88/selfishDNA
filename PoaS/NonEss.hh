#ifndef NonEssHeader
#define NonEssHeader

#include "Pearl.hh"

class NonEss : public Pearl
{
	public:
		int type;	// Number for NonEssgene.. you need them ALLLLLLL >:)
		double mobility;


		NonEss();  // constructor
		virtual Pearl* clone() const;
		~NonEss(); // deconstructor
		NonEss(int t); // Constructor that gets toxin type as an argument
		explicit NonEss(const NonEss &noness); // Copy constructor
};
#endif
