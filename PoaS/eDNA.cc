#include "eDNA.hh"
#include "Genome.hh"

/**********
		DEFAULT CONSTRUCTOR AND DECONSTRUCTOR STUFF
																							*********/
eDNA::eDNA()					// Default constructor doesn't do anything :)
{
	Fragments=new vector< list<Pearl*> >();			// This holds a list of fragments of DNA (i.e. small pieces of pearl-on-a-string)
}

eDNA::~eDNA()					// Deconstructor deletes the genome and all Pearls it's pointing to
{
	bool V = FALSE;
	iter i;
	if(V) cout << "Deconstructing eDNA pool" << endl;
	/* NOG MAKEN!! */
}

void eDNA::FragmentiseGenome(Genome * G, float fraction)																		// When the cell dies, the DNA is chopped up before it's dropped on the grid
{
	bool V = false;
	if(G->StringOfPearls->size() == 0)
		return;
	if(V)cout << "Before new genome fragmentised into environment: " << endl;
	if(V)
	{
		for (vector< list<Pearl*> >::iterator it = Fragments->begin(); it != Fragments->end(); ++it)
		{
			cout << G->ListContent(&(*it)) << ", ";
		}
		cout << endl;
		//	delete G;			// Because the simulation is synchronous, this happens AFTER everyone has updated or it goes horribly wrong :)
	}
	Pearl *pearl;
	if(V)cout << "Genome that will be fragmented into environment:" << endl;
	if(V)cout << G->ListContent(NULL) << endl;
	iter i_start = G->StringOfPearls->begin();
	iter i_end = G->StringOfPearls->begin();
	int stringsize = G->StringOfPearls->size();	// To keep std::advance from going over the edge
	int fragsize;	// Generate a piece of X pearls, where X should not be bigger than the remaining stringsize

	do  																			// Do while is like a normal while loop, but ensures the thing is executed at least once
	{
		i_start = i_end;	
		fragsize = genrand_int(min(3,stringsize),(min(8,stringsize)));			// NOTE HARD CODED, ALSO CHANGE IN JUMP FUNCTION
		
		if(V)cout << stringsize << endl;
		if(V)cout << fragsize << endl;

		list<Pearl*> DNAfragment; 	//(i_start, i_end);
		G->CopyPartOfGenomeToTempList(i_start,fragsize,DNAfragment);
		if(genrand_real1() < fraction) Fragments->push_back(DNAfragment);
		stringsize-=fragsize;
		DNAfragment.clear();
		advance(i_end,fragsize);
	}	while(i_end!=G->StringOfPearls->end());

	if(V)cout << "After new genome fragmentised into environment: " << endl;
	
	if(V)
	{
		for (vector< list<Pearl*> >::iterator it = Fragments->begin(); it != Fragments->end(); ++it)
		{
			cout << G->ListContent(&(*it)) << ", ";
		}
		cout << endl;
	}
}

void eDNA::Degradate(float degr)																		// When the cell dies, the DNA is chopped up before it's dropped on the grid
{
	bool V = FALSE;

	if(V)
	{
		cout << "Before degr: " << endl;
		cout << Fragments->size() << endl;
		cout << endl;
	}

	int pos = 0;
	for (fr_iter frit = Fragments->begin(); frit != Fragments->end(); frit++)
	{
		pos ++;

		if(genrand_real1() < degr)								// default 0.1 LS: shouldn't these chances better be parameterised into the par file?
		{
			if(V) cout << "Deleting eDNAfragment at pos " << pos << endl;
			frit = DeleteFragment(frit);
			frit--;
		}
	}
	if(V)
	{
		cout << "After degr: " << endl;
		cout << Fragments->size() << endl;
		cout << endl;
	}
}

eDNA::fr_iter eDNA::DeleteFragment(fr_iter frit)
{
	bool V = FALSE;
	if(V) cout << "Deleting fragment of size " << frit->size() << endl;
	iter it;

	it = frit->begin();
	while(it!=frit->end())
	{
		if(V) cout << "Deleting bead" << endl;
		delete *it;
		it = frit->erase(it);
	}

	frit = Fragments->erase(frit);
	return frit;
}
