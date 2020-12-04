#ifndef NoncodingHeader
#define NoncodingHeader

#include "Pearl.hh"

class Noncoding : public Pearl
{
	public:
		int type;	// Type is maintained when it turns into noncoding
        double mobility;
		Noncoding();  // constructor
		virtual Pearl* clone() const;
		~Noncoding(); // deconstructor
		Noncoding(int t); // Constructor that gets type as an argument
		explicit Noncoding(const Noncoding &nc); // Copy constructor
};
#endif
