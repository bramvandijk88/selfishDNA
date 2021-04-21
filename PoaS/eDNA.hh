#ifndef eDNAHeader
#define eDNAHeader

#include <typeinfo>
#include "Header.hh"
#include "Genome.hh"
#include "Transposase.hh"

class eDNA
{
	public:
		std::vector< list<Pearl*>  > *Fragments;	// eDNA is a vector of Genome pointers (but then, instead of complete genomes, they are fragments)
		typedef std::list<Pearl*>::iterator iter;
		typedef std::vector< std::list<Pearl*> >::iterator fr_iter;
		eDNA();
		~eDNA();
		void FragmentiseGenome(Genome * G, float fraction);
		void Degradate(float);
		void InfluxDNA(float);

		eDNA::fr_iter DeleteFragment(fr_iter fri);
		bool IsTransposase(Pearl *Pearl) const;
		bool HasTransposase(list<Pearl*>* frag);
};

#endif
