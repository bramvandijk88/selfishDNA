#include "Genome.hh"
#include "eDNA.hh"

/**********
		DEFAULT CONSTRUCTOR AND DECONSTRUCTOR STUFF
																							*********/
Genome::Genome()					// Default constructor doesn't do anything :)
{
	int i;	// Doet volgens mij nix :D
	StringOfPearls=NULL;

}

Genome::~Genome()					// Deconstructor deletes the genome and all Pearls it's pointing to
{
	bool V = FALSE;
	iter i;
	if(V) cout << "Deconstructing genome of generation " << generation_ << endl;
	if(V) cout << ListContent(StringOfPearls) << endl;
	if(StringOfPearls!=NULL)
	{																				// E.g. StringOfPearls   1   2   3   4   5
		i=StringOfPearls->begin();						//											 ^
		while(i!=StringOfPearls->end())				//                           ^ etc.
		{
			delete(*i);
			i++;
		}
		i=StringOfPearls->erase(StringOfPearls->begin(), StringOfPearls->end());


		delete StringOfPearls;
		StringOfPearls=NULL;
		if(V) cout << "Done cleaning up genome" << endl;
	}
	
}

/**********
		GENOME GENERATION / COPYING
			GenerateGenome = random genome generation (i.o.w. only at start of simulation)
			CloneGenome = clone genome from parent (or part of genome)
			CopyPartOfGenome = After making a new (empty genome), this copies part of a GIVEN genome to this genome
																				*********/
void Genome::GenerateGenome(int init_nr_hkgenes, int init_nr_noness, int init_nr_noncoding, int init_nr_tra, 
							float gene_cost, float transp_cost, float genome_size_cost, float init_mob, float effect_noness)
{
	gene_cost_ = gene_cost;
	transp_cost_ = transp_cost;
	genome_size_cost_ = genome_size_cost;
	string rowcoltime = "";
	parent = NULL;

	bool V = false;
	int i,j;

	fitness_effect_noness = effect_noness;
	Transposase* tra;
	HK* hk;
	Noncoding* nc;
	NonEss* noness; 

	transformant=0;

	//printf("%p\n", genomesize_);	// With this you can check if pointer segfaults, which means that it doesnt exist
	// Set all copy-number and genome size thingies to 0
	genomesize_= 0;
	HKgenes_=0;	

	if(V) cout << "Making new Pearl list" << endl;
	StringOfPearls=new list<Pearl*>(); 		// Make a new (empty) list of Pearl pointers	

	for(i=0;i<init_nr_hkgenes;i++)
	{																																			//LS: everyone gets (some?) 		?
		hk = new HK(i);					// Let hk point to this HK gene
		hk->mobility = 0.0;
		StringOfPearls->push_back(hk);  	    	// Add pointer to hk to list
		genomesize_++;
		HKgenes_++;																													//LS: when making the HKgenes also the minimum necessity is set? (somewhere else you say HKgenes_ is the minimum of HKgenes necessary)
	}
	for(i=0;i<init_nr_noness;i++)
	{																																			//LS: everyone gets (some?) 		?
		noness = new NonEss(i);					// Let hk point to this HK gene
		noness->mobility = 0.0;
		StringOfPearls->push_back(noness);  	    	// Add pointer to hk to list
		genomesize_++;		
	}
	//int nr_nc = genrand_real1()*init_nr_noncoding*2;
	int nr_nc = init_nr_noncoding;
	for(i=0;i<nr_nc;i++)
	{
		nc = new Noncoding(-1);					// Default the "type" of a non-coding part corresponds to a random type of toxin/antitoxin
		nc->mobility = 0.00;
		StringOfPearls->push_back(nc);  	    	// Add pointer to hk to list
		genomesize_++;
	}
	

	if(V) cout << "VOOR:" << endl;
	if(V) cout << ListContent(NULL) << endl;

	// Shuffling for lists isn't very efficient, but the 3 lines below do the trick
	   vector< Pearl* > Vec( StringOfPearls->begin(), StringOfPearls->end() );
	   random_shuffle( Vec.begin(), Vec.end() );
	   StringOfPearls->assign( Vec.begin(), Vec.end() );

	
		for(i=0;i<init_nr_tra;i++)
		{
			iter ii = StringOfPearls->begin();
			int randpos = (int)(genrand_real1()*genomesize_);
			advance(ii,randpos);			
			tra = new Transposase(-1);					// Default the "type" of a non-coding part corresponds to a random type of toxin/antitoxin
			tra->mobility = init_mob;
			StringOfPearls->insert(ii,tra);  	    	// Add pointer to hk to list
			genomesize_++;
		}
	
	if(V) cout << "NA:" << endl;
    if(V) cout << ListContent(NULL) << endl;
	if(V)cout << endl << endl;
	viable = TRUE;
	Create_Gene_Lists();			// This creates a look-up structure to quickly figure our resistance and stuff. :)
	CalculateCompStrength();
	GenomeAtBirth = ListContentShort(NULL);
//	sleep(1);
	generation_ = 1;
	time_infected_ = 0;
}

void Genome::CloneGenome(Genome *parent_genome)
{
	//double mob_start = GetPearlMobility(parent->StringOfPearls->front());
	//double mob_end = GetPearlMobility(parent->StringOfPearls->back());

	bool V = FALSE;
	StringOfPearls=new list<Pearl*>();
	if(V) cout << "Trying to clone genome of size " << parent_genome->StringOfPearls->size() << endl;
	if(V) cout << ListContent(parent_genome->StringOfPearls) << endl;
	CopyPartOfGenome(parent_genome->StringOfPearls->begin(),parent_genome->StringOfPearls->end());		// Copies full genome to this genome
	
	string rowcoltime = "";
	genomesize_=parent_genome->genomesize_;
	parent=parent_genome;
	generation_=parent_genome->generation_;	
	HKgenes_=parent_genome->HKgenes_;
	genome_size_cost_=parent_genome->genome_size_cost_;
	gene_cost_=parent_genome->gene_cost_;
	transp_cost_=parent_genome->transp_cost_;
	fitness_effect_noness=parent_genome->fitness_effect_noness;
	extra_death=parent_genome->extra_death;
	compstrength=parent_genome->compstrength;
	time_infected_=parent_genome->time_infected_;
	transformant=0;	
	generation_++;
}

void Genome::RecombineGenomes(Genome *parent_genome, Genome *parent_genome2)
{
	bool V = FALSE;
	StringOfPearls=new list<Pearl*>();
	if(V) cout << "Trying to clone genome of size " << parent_genome->StringOfPearls->size() << endl;
	if(V) cout << ListContent(parent_genome->StringOfPearls) << endl;
	Recombine(parent_genome,parent_genome2);			// Copies full genome to this genome	
	string rowcoltime = "";
	genomesize_=parent_genome->genomesize_;
	parent=parent_genome;
	generation_=parent_genome->generation_;	
	HKgenes_=parent_genome->HKgenes_;
	genome_size_cost_=parent_genome->genome_size_cost_;
	transp_cost_=parent_genome->transp_cost_;
	gene_cost_=parent_genome->gene_cost_;
	fitness_effect_noness=parent_genome->fitness_effect_noness;
	extra_death=parent_genome->extra_death;
	compstrength=parent_genome->compstrength;
	time_infected_=parent_genome->time_infected_;
	transformant=0;	
	generation_++;
}

void Genome::CopyPartOfGenome(iter begin,iter end)
{
	iter i;
	Pearl *pearl;
	i=begin;
	while(i!=end)
	{
		pearl=(*i)->clone();
		pearl->num_vertical_transfers_++;
		StringOfPearls->push_back(pearl);
		i++;
	}
}

void Genome::Recombine(Genome *parent_genome, Genome *parent_genome2)
{
	iter i;
	iter e;
	bool V = false;
	if(V) cout << "Parent 1: " << ListContent(parent_genome->StringOfPearls) << endl;
	if(V) cout << "Parent 2: " << ListContent(parent_genome2->StringOfPearls) << endl;
	int halfway = parent_genome->StringOfPearls->size()/2;
	int halfway_2 = parent_genome2->StringOfPearls->size()/2;
	if(V) cout << "Halfway point 1: " << halfway << endl;
	if(V) cout << "Halfway point 2: " << halfway_2 << endl;

	Pearl *pearl;
	i=parent_genome->StringOfPearls->begin();
	e=parent_genome->StringOfPearls->begin();
	advance(e,halfway);

	while(i!=e)
	{
		pearl=(*i)->clone();
		pearl->num_vertical_transfers_++;
		StringOfPearls->push_back(pearl);
		i++;
	}

	i=parent_genome2->StringOfPearls->begin();
	advance(i,halfway_2);
	e=parent_genome2->StringOfPearls->end();

	while(i!=e)
	{
		pearl=(*i)->clone();
		pearl->num_vertical_transfers_++;
		StringOfPearls->push_back(pearl);
		i++;
	}
	if(V) cout << "Result: " << ListContent(NULL) << endl;
	if(V) cout << "Done with recombining" << endl << endl;

}


void Genome::CopyPartOfGenomeToTempList(iter begin,int size,list<Pearl*> &StringOfPearlsTemp)			// For inversions and duplications of multiple genes
{
	iter i;
	Pearl *pearl;
	i=begin;
	while(i!=StringOfPearls->end())
	{
		pearl=(*i)->clone();																					// IS THIS CLONE NECESSARY AND EVEN DESIRED?? CHECK LATER!!
		StringOfPearlsTemp.push_back(pearl);
		if(StringOfPearlsTemp.size()==size) return;
		i++;
	}
	return;
}



bool Genome::Viable()
{
	bool V = false;
	// cout << ListContent(NULL) << endl;
	// cout << compstrength << endl;
	if(compstrength<=0.0) { return false; }
	Create_Gene_Lists();
	list<int> uniquehks;
	if(HKgenes_ > 0)
	{
		uniquehks = HKgenes;
		uniquehks.sort();
		uniquehks.unique();
		if(V)
		{
			cout << "Found " << uniquehks.size() << " hk genes. Need " << HKgenes_ << endl;		// For clarity: HKgenes_ is just a number to check the number of genes the species should minimally have, while without the hyphen (HKgenes) is actually the list of all HKgenes in that genome
		}
	}
	return (HKgenes_ == uniquehks.size());
}

void Genome::CalculateCompStrength()
{
			bool V=false;
			if(StringOfPearls->size() == 0)
			{
				compstrength = 0.0;
				extra_death = 1.0;
			}
			else 			// A) Check if all hk genes exist in genome LS:!! nu mag je niet meer unieke HKgenes hebben dan het minimum noodzakelijk?!?
			{
				if(V) cout << "Determine fitness because all hk genes are present" << endl;					

				list<int> noness;
				if(nonessential.size() > 0)
				{
					noness = nonessential;
					noness.sort();
					noness.unique();				
				}

				compstrength = 1+noness.size()*fitness_effect_noness;

				if(V) cout << "Fitness basis set to " << compstrength  << endl;
				int HKgenes_redundant = 0;
				if(HKgenes.size() > HKgenes_) HKgenes_redundant = HKgenes.size() - HKgenes_;				
				//HKgenes_redundant = 0;
				if(V) cout << "genes:" << HKgenes_redundant+nonessential.size() << endl;
				compstrength -= (HKgenes_redundant+nonessential.size())*gene_cost_;
				if(V) cout << "Transposases:" << transposases.size() << endl;
				if(V) cout << "Cost: " << transp_cost_ << endl;
				compstrength -= transposases.size()*transp_cost_;
				if(V) cout << "Comstrength after cost:" << compstrength << endl;
				// compstrength = 1.0;
				compstrength -= (StringOfPearls->size()*genome_size_cost_);
				if(V) cout << "Fitness l furthermore " << (StringOfPearls->size()*genome_size_cost_)  << endl;
				if(StringOfPearls->size()>1000) compstrength = 0.0;
				compstrength = max(compstrength,0.0);
				
			}

			if(V)cout << "Fitness set to " << compstrength << endl;

}



/**********

		GENOME MUTATIONS (In the examples, ^ = iterator)
			MutateGenome applies all mutations to a genome
												i++ moves the iterator by 1 step																						A  B  C  D  E
																											   	   							   										       ^
												GeneDeletions = Gene at position of iterator is deleted										  A  C  D  E
																											   	   							   													 ^
				GeneDupl = Gene at position of iterator is duplicated, NOTE iterator stays on that gene     A  C* C  D  E
																				  	  							  																						^
			    TandemDupl = A sequence of genes is duplicated, NOTE iterator stays on that gene     		  A  C  C* D* E* C  D  E
																				  	  							           																						 ^
			    etc....
																																																				*********/
bool Genome::MutateGenome(float gene_mob, float loss, float dupl, float tdupl, float tdel, float inv, float gtn, float ntg, float gendisc)
{
	iter i = StringOfPearls->begin();
	bool mutated = FALSE;
	bool V = FALSE;
	while(i!=StringOfPearls->end())
		{
			/**** OTHER MUTATIONS ****/
				float mobmut_chance = genrand_real1();
				if(V) cout<< "mobmut_chance = " << mobmut_chance << endl;

				if(mobmut_chance < gene_mob) 									// LS: move to par file - DONE
				{
					MutateMobility(i), mutated=TRUE;
				}
				
				float mut_chance = genrand_real1 ();
				if(V) cout <<"mut_chance=" << mut_chance << endl;
				if(V)cout  <<"chance of A mutation happening = "<<loss+dupl+tdupl+tdel+inv+gtn+ntg << endl;
	
				if (mut_chance > loss+dupl+tdupl+tdel+inv+gtn+ntg)
				{
					i++;
				}
				else if(mut_chance < loss)
				{
					i = GeneLoss(i), mutated=TRUE;
				}
				else if(mut_chance < loss+dupl) 						// ELSE if, so these mutations are mutually exclusive on the level of a pearl (LS: note, why not happen at the same time?)
				{
					i = GeneDupl(i), mutated=TRUE;
				}
				else if(mut_chance < loss+dupl+tdupl)
				{
					i = TandemDupl(i), mutated=TRUE;
				}
				else if(mut_chance < loss+dupl+tdupl+tdel)
				{
					i = TandemDel(i), mutated=TRUE;
				}
				else if(mut_chance < loss+dupl+tdupl+tdel+inv)
				{
					i = Invert(i), mutated=TRUE;
				}
				else if(mut_chance < loss+dupl+tdupl+tdel+inv+gtn)
				{
					if(!IsNoncoding(*i)) i = Gene_To_NonCoding(i), mutated=TRUE;
				}
				else if(mut_chance < loss+dupl+tdupl+tdel+inv+gtn+ntg)
				{
					if(IsNoncoding(*i))  i = NonCoding_To_Gene(i), mutated=TRUE;
				}

		}
		

		
		Create_Gene_Lists();
		// DUPLICATING VIA CHUNKS (analogous to via extracellular DNA)
		




	bool discovered_gene = GeneDiscovery(gendisc);
	if(discovered_gene) mutated=TRUE;

	//if(V) cout << mutated << endl;
	//if(V) if(mutated) cout << this->ListContent(NULL) << endl << endl; // This is more or less arbitrary.. sometimes it prints the genome of a child

	return mutated;
}


bool Genome::TransposonDynamics(list<Pearl*> *StringOfPearls_scanned,float rate_mult, float break_chance)
{	
	bool V = false;
	if(V) cout << "-= START TRANSPOSONDYNAMICS =- " << endl;
	if(!Viable())
	{
		if(V) cout << "This is not a viable host cell, so transposons are unable to replicate too" << endl;
		return false;
	}
	bool mutated = false;
	int num_jumps = 0;
	float rate = rate_mult;
	if(rate < 0.0) rate = 0.0;
	int max_num_jumps = 10;
	bool hgt = false;
	if(StringOfPearls_scanned == NULL)
	{
		 StringOfPearls_scanned = this->StringOfPearls;
		 if(V) cout << "Using this stringOfPearls, so this is a intragenomic transposition event" << endl;
	}
	else
	{
		if(V) cout << "Using StringOfPearls passed in ARG 1, so this is a eDNA Fragment" << endl;
		hgt = true;
	}

	
	std::list<Pearl*> StringOfPearlsTemp;
	iter i = StringOfPearls_scanned->begin();
	CopyPartOfGenomeToTempList(i,StringOfPearls_scanned->size(),StringOfPearlsTemp);


	/// TRANSPOSON EITHER JUMPING FROM eDNA to GENOME, OR DUPLICATING WITHIN GENOME 
	if(V) cout << "Scanned DNA fragment: " << endl << ListContent(&StringOfPearlsTemp) << endl;
	if(V) cout << "For possible insertion into: " << endl << ListContent(StringOfPearls) << endl;
	i = StringOfPearlsTemp.begin();
	Pearl *pearl;
	int count = 0;
	while(i!=StringOfPearlsTemp.end())
	{		
		// if(V)cout << count << endl;
		count++;
		if(IsTransposase(*i))
		{							
			double mobility = GetPearlMobility(*i);
			//if(genrand_real1() < pow(mobility,2)*rate)
			if(genrand_real1() < mobility*rate)
			{
				num_jumps++;
				pearl=(*i);
				
				mutated = true;
				iter ii = StringOfPearls->begin();
				
				int randpos = (int)(genrand_real1()*(StringOfPearls->size()-1));
				if(V)cout << randpos << "," << count << endl;
				
				advance(ii,randpos);
				if(*i == *ii) continue; // Transposon inserting into itself, so just skip the entire thing since the genome remains the same. 
				if(V)cout << "Inserting into position " << randpos << endl;
				
				// Make a copy of the transposon
				Transposase *tra;	
				tra = new Transposase(0);
				if(hgt) tra->num_horizontal_transfers_++;
				else tra->num_jumps_++;
				
				tra->mobility = mobility;
				
				
				// Where it inserted, the existing gene might be broken (chance = break_chance. action = delete, replace with non-coding DNA)
				bool insert_in_gene = false;
				if(genrand_real1() < break_chance) insert_in_gene = true;
				if(insert_in_gene)
				{						
					delete *ii;
					ii = StringOfPearls->erase(ii);											
					Noncoding *nc;									
					nc = new Noncoding(0);							// Let 
					nc->mobility = 0.0;								// Newly introduced 
					ii = StringOfPearls->insert(ii,nc);
								
				}
				ii = StringOfPearls->insert(ii,tra);
				// Insert a copy of the transposon at the determined random posotion.									
			
				if(V) cout << "Insertion happened. New genome: " << endl << ListContent(StringOfPearls) << endl;				
				if(V) cout << "Insertion happened. New scange: " << endl << ListContent(StringOfPearls_scanned) << endl;				
				if(!Viable())
				{
					if(V)cout << "Breaking out of transposon dynamics since genome is already inviable due to HK genes being broken." << endl;
					compstrength = 0.0;
					return true;
				}
			}
			//  if(num_jumps > max_num_jumps) break;
		}

		delete *i;
		i = StringOfPearlsTemp.erase(i);			
		
		// i++;
	}
	StringOfPearlsTemp.clear();
	if(V) cout << "Scanned DNA fragment:  " << endl << ListContent(&StringOfPearlsTemp) << endl;
	if(V) cout << "Genome aFter scanning: " << endl << ListContent(StringOfPearls) << endl;

	if(mutated) 
	{
		Create_Gene_Lists();
		CalculateCompStrength();
		return true;
	}			
	return false;
}



void Genome::MutateMobility(iter ii)
{
	if(IsTransposase(*ii))
		dynamic_cast<Transposase *>(*ii)->mobility = min(1.0,max(0.0,( (dynamic_cast<Transposase *>(*ii)->mobility) + (2.0*genrand_real1()-1.0) * 0.1)));
}

Genome::iter Genome::GeneLoss(iter ii)
{
	delete *ii;
	ii = StringOfPearls ->erase(ii);
	genomesize_--;
	return ii;
}

Genome::iter Genome::GeneDupl(iter ii)
{
	Pearl *Pearl;
	Pearl=(*ii)->clone();
	StringOfPearls->insert(ii, Pearl);
	genomesize_++;
	return ii;
}

Genome::iter Genome::TandemDel(iter ii)
{
	bool V = FALSE;

	int size = genrand_int(2,StringOfPearls->size()/4);
	int deleted = 0;

	while(ii!=StringOfPearls->end())
	{
		delete *ii;
		ii = StringOfPearls ->erase(ii);
		genomesize_--;
		deleted++;
		if (deleted == size) return ii;
	}

	return ii;
}

Genome::iter Genome::TandemDupl(iter ii)
{
	bool V = FALSE;
	std::list<Pearl*> StringOfPearlsTemp;
	int size = genrand_int(2,StringOfPearls->size()/4);
	Pearl *pearl;

	CopyPartOfGenomeToTempList(ii,size,StringOfPearlsTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it

	if(V)
	{
		cout << "Before" << endl;
		cout << ListContent(NULL) << endl;
		cout << "Made temp stukje dna" << endl;
		cout << ListContent(&StringOfPearlsTemp) << endl;
	}
	iter t= StringOfPearlsTemp.begin();
	while(t!=StringOfPearlsTemp.end())
	{
		if(V)cout << ListContent(&StringOfPearlsTemp) << endl;
		pearl=(*t);
		if(V)cout << "Inserting bead.." << endl;
		StringOfPearls->insert(ii, pearl);
		genomesize_++;
		if(V)cout << "..inserted bead" << endl;
		t++;
	}
	StringOfPearlsTemp.clear();
	if(V)cout << "After" << endl;
	if(V)cout << ListContent(NULL) << endl;
	return ii;
}

Genome::iter Genome::Invert(iter ii)
{
	bool V = FALSE;
	std::list<Pearl*> StringOfPearlsTemp;
	int size = genrand_int(2,StringOfPearls->size()/4);
	int deleted=0;
	Pearl *pearl;

	CopyPartOfGenomeToTempList(ii,size,StringOfPearlsTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it

	if(V)
	{
		cout << "Before" << endl;
		cout << ListContent(NULL) << endl;
		cout << "Made temp stukje dna" << endl;
		cout << ListContent(&StringOfPearlsTemp) << endl;
	}

	while(ii!=StringOfPearls->end())
	{
		delete *ii;
		ii = StringOfPearls ->erase(ii);
		genomesize_--;
		deleted++;
		if (deleted == size) break;	// Done deleting, next we insert it in the reverse order
	}


	rev_iter t;
	t = StringOfPearlsTemp.begin();				// Reverse iterator starts
	while(t!=StringOfPearlsTemp.end())
	{
		if(V)cout << ListContent(&StringOfPearlsTemp) << endl;
		pearl=(*t);
		if(V)cout << "Inserting bead.." << endl;
		StringOfPearls->insert(ii, pearl);
		ii--;
		genomesize_++;
		if(V)cout << "..inserted bead" << endl;
		t++;
	}
	StringOfPearlsTemp.clear();
	if(V)cout << "After" << endl;
	if(V)cout << ListContent(NULL) << endl;
	if(V)cout << endl << endl;
	return ii;
}

Genome::iter Genome::Gene_To_NonCoding(iter ii)
{
	bool V=FALSE;
	Noncoding * nc;
	
	HK * hk;
	int type = -1;		// Nonecoding type heeft altijd -1. Deze moet nooit gebruikt worden voor "echte" genen, maar het is dus good practice om deze op -1 te zetten zodat het absoluut fout gaat als het perongeluk in een normaal gen overerft
	double mob = GetPearlMobility(*ii);
	//if(V) cout << "Inheriting mobility"
	if(V)cout << "Befor: " << ListContent(NULL) << endl;
	nc = new Noncoding(type);	 // Make new noncoding pearl with a random type (random with respect to our toxin/antitoxins)
	nc->mobility = 0.0; 	 	// Inherit mobility from original gene
	delete *ii;						 // Delete the original gene
	ii = StringOfPearls ->erase(ii); //Erase from list
	ii = StringOfPearls->insert(ii,nc);
	if(V)cout << "After: " << ListContent(NULL) << endl;
	//genomesize_--;
	return ii;
}

Genome::iter Genome::NonCoding_To_Gene(iter ii)
{
	bool V=FALSE;
	Noncoding * nc;
	
	Transposase * tra;
	HK * hk;
	// RESISTANCE GENE DISCOVERY
	if(genrand_real1() < 0.5)
	{
		int type = genrand_int(0,HKgenes_-1);		// Counts from 0 because computers :) 
		double mob = GetPearlMobility(*ii);
		if(V)cout << "Befor: " << ListContent(NULL) << endl;
		double randnum = genrand_real1();	
		hk = new HK(type);	 				// Make new arg pearl with type of the non-coding gene (in place for future modularity)	
		hk->mobility = 0.0; 	 				// For now, only NC "flanks" can infer gene mobility
		delete *ii;						 		// Delete the original non-coding gene
		ii = StringOfPearls ->erase(ii); 		// Erase from list
		ii = StringOfPearls->insert(ii,hk);	// Insert new arg	
	}
	else
	{
		int type = 0;							// Counts from 0 because computers :) 
		double mob = GetPearlMobility(*ii);
		if(V)cout << "Befor: " << ListContent(NULL) << endl;
		double randnum = genrand_real1();	
		tra = new Transposase(type);	 				// Make new arg pearl with type of the non-coding gene.. do we want this? Dunno
		tra->mobility = genrand_real1(); 	 				// Random mobility
		delete *ii;						 		// Delete the original non-coding gene
		ii = StringOfPearls ->erase(ii); 		// Erase from list
		ii = StringOfPearls->insert(ii,tra);	// Insert new arg	
	}
	if(V)cout << "After: " << ListContent(NULL) << endl;
	return ii;
}

bool Genome::GeneDiscovery(double gendisc)
{

	bool V = FALSE;
	int type, randpos;
	iter ii;
	bool mutated = FALSE;

	if(V) cout << "GeneDiscovery called" << endl;
	if(genrand_real1()<gendisc)
	{
		if(V) cout << ListContent(NULL) << endl;
		if(V) cout << "Discovering nc with chance " << gendisc << endl;
		Noncoding *nc;
		ii = StringOfPearls->begin();
		randpos = (int)(genrand_real1()*genomesize_);
		if(V) cout << "Randpos: " << randpos << " van de " << genomesize_ << endl;
		advance(ii,randpos);
		type = genrand_int(0,HKgenes_-1);		// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY
		if(V) cout << "Type: " << type << " van de " << HKgenes_ << endl;
		nc = new Noncoding(type);						// Let toxin point to a new Toxin
        nc->mobility = genrand_real1();	// Newly introduced 
		ii = StringOfPearls->insert(ii,nc);  		     	// Add pointer to toxin to pearl-list
		if(V) cout << ListContent(NULL) << endl;
		genomesize_++;
		mutated=TRUE;
	}	

	//cout << "GeneDiscovery called" << endl;
	int nr_discovered = 0;
	while(genrand_real1()<gendisc) nr_discovered++;
	for(int i = 0; i < nr_discovered; i++)
	{
		if(V) cout << ListContent(NULL) << endl;
		if(V) cout << "Discovering transposon with chance " << gendisc << endl;
		Transposase *tra;
		ii = StringOfPearls->begin();
		randpos = (int)(genrand_real1()*genomesize_);
		if(V) cout << "Randpos: " << randpos << " van de " << genomesize_ << endl;
		advance(ii,randpos);
		type = genrand_int(0,HKgenes_-1);		// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY
						
		// Transposon discovered
		if(V) cout << "Type: " << type << " van de " << HKgenes_ << endl;
		tra = new Transposase(type);						// Let toxin point to a new Toxin
        tra->mobility = 0.75+0.25*genrand_real1();	// Newly introduced gets random mobility	
		ii = StringOfPearls->insert(ii,tra);  		     	// Add pointer to toxin to pearl-list
		
		genomesize_++;
		mutated=TRUE;
	}
	
	return mutated;
}


void Genome::Create_Gene_Lists()
{
	
	/*
	This function generates three things:
	1) A vector ('array') of size TOT_NR_TOXINS that holds whether this individual is resistance to toxins
	2) A list of toxins this individual produces.
	3) A list of household genes this individual produces
	*/

	bool V = false;

	list<Pearl*>* StringOfPearls = this->StringOfPearls;
	iter i = StringOfPearls->begin();
	HK* hk;
	Transposase* tra;
	NonEss* noness;

	if(V) cout << "Making lookup table" << endl;

	HKgenes.clear();
	transposases.clear();
	nonessential.clear();

	if(V) cout << "Before:" << ListContent(NULL) << endl;
	if(V) cout << endl;

	while(i!=StringOfPearls->end())
	{		
		if(IsHK(*i))
		{
			hk=dynamic_cast<HK *>(*i);
			HKgenes.push_back(hk->type);
		}
		if(IsTransposase(*i))
		{
			tra=dynamic_cast<Transposase *>(*i);
			transposases.push_back(tra->type);
		}
		if(IsNonEssential(*i))
		{
			noness=dynamic_cast<NonEss *>(*i);
			nonessential.push_back(noness->type);
		}
		//IsNonEssential
		i++;
	}

	if(V) cout << "After:" << ListContent(NULL) << endl;
	if(V) cout << "Tra:" << transposases.size() << endl;

}



string Genome::ListContent(list<Pearl*> *StringOfPearls, bool ignore_hk, bool ignore_noncoding, bool includemobility) // If *Pearl is not given, it prints the entire genome
{
	bool V = FALSE;
	string GenomeContent;
	if(StringOfPearls == NULL) StringOfPearls = this->StringOfPearls;
	iter i = StringOfPearls->begin();
	Transposase *tra;
	HK *hk;
	Noncoding * nc;
	NonEss * noness;

	if(V)
		cout << "Attempting to list genome with length " << StringOfPearls->size() << endl;
	// GenomeContent+="\033[00m---";																					//LS: klopt Dit?
	while(i!=StringOfPearls->end())
	{
		if(ignore_noncoding)
			if(IsNoncoding(*i)) {i++;continue;}
		if(ignore_hk)
			if(IsHK(*i)) {i++;continue;}
		stringstream stringtemp;
		stringtemp.precision(2);
		stringtemp << fixed;
		// if(i!=StringOfPearls->begin()) GenomeContent +="-";
		if(IsTransposase(*i))
		{
			tra=dynamic_cast<Transposase*> (*i);
			stringtemp << "\033[1;41mT";
			if(includemobility) stringtemp <<  ":" <<tra->mobility;
			stringtemp << "\033[00m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}		
		else if(IsHK(*i))
		{
			hk=dynamic_cast<HK *>(*i);
			// stringtemp << "\033[1;42mH" << hk->type;
			stringtemp << "\033[1;44mH";
			if(includemobility) stringtemp << ":" <<hk->mobility;
			stringtemp << "\033[00m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsNoncoding(*i))
		{			
			nc=dynamic_cast<Noncoding *>(*i);
			stringtemp << "\033[1;90mn";
			if(includemobility) stringtemp << ":" << nc->mobility;
			stringtemp << "\033[00m";					
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsNonEssential(*i))
		{			
			noness=dynamic_cast<NonEss *>(*i);
			stringtemp << "\033[1;42mA";
			if(includemobility) stringtemp << ":" << noness->mobility;
			stringtemp << "\033[00m";					
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else
		{
			GenomeContent+="He?";
		}
		i++;
	}
	// GenomeContent+="\033[49m---\n";
	return GenomeContent;
}


string Genome::ListContentShort(list<Pearl*> *StringOfPearls, bool ignore_hk, bool ignore_noncoding, bool includemobility) // If *Pearl is not given, it prints the entire genome
{
	bool V = FALSE;
	string GenomeContent;
	if(StringOfPearls == NULL) StringOfPearls = this->StringOfPearls;
	iter i = StringOfPearls->begin();
	Transposase *tra;
	HK *hk;
	Noncoding * nc;
	NonEss * noness;


	if(V)
		cout << "Attempting to list genome with length " << StringOfPearls->size() << endl;
	while(i!=StringOfPearls->end())
	{
		if(ignore_noncoding)
			if(IsNoncoding(*i)) {i++;continue;}
		if(ignore_hk)
			if(IsHK(*i)) {i++;continue;}
		stringstream stringtemp;
		stringtemp.precision(2);
		stringtemp << fixed;
		if(IsTransposase(*i))
		{
			tra=dynamic_cast<Transposase*> (*i);
			stringtemp << "T";
			if(includemobility) stringtemp <<  ":" <<tra->mobility;
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsHK(*i))
		{
			hk=dynamic_cast<HK *>(*i);
			stringtemp << "H";
			if(includemobility) stringtemp << ":" <<hk->mobility;
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsNoncoding(*i))
		{
			nc=dynamic_cast<Noncoding *>(*i);
			stringtemp << "n";			
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsNonEssential(*i))
		{			
			noness=dynamic_cast<NonEss *>(*i);
			stringtemp << "A";			
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else{
			GenomeContent+="He?";
		}
		i++;
	}
	// GenomeContent+="\033[49m---\n";
	return GenomeContent;
}


/*********
				Get genome stat functions
********/

// Below are all the functions that check whether a Pearl is of a certain type

bool Genome::IsTransposase(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(Transposase));	// typeid determines the class of an object at runtime
}

bool Genome::IsHK(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(HK));	// typeid determines the class of an object at runtime
}

bool Genome::IsNonEssential(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(NonEss));	// typeid determines the class of an object at runtime
}

bool Genome::IsNoncoding(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(Noncoding));	// typeid determines the class of an object at runtime
}

int Genome::GetPearlType(Pearl *Pearl) const
{
	if(IsHK(Pearl))
		return dynamic_cast<HK *>(Pearl)->type;
	else if(IsNoncoding(Pearl))
		return dynamic_cast<Noncoding *>(Pearl)->type;
	else if(IsTransposase(Pearl))
		return dynamic_cast<Transposase *>(Pearl)->type;
}

double Genome::GetPearlMobility(Pearl *Pearl) const
{
	if(IsHK(Pearl))
		return dynamic_cast<HK *>(Pearl)->mobility;
	else if(IsNoncoding(Pearl))
		return dynamic_cast<Noncoding *>(Pearl)->mobility;
	else if(IsTransposase(Pearl))
		return dynamic_cast<Transposase *>(Pearl)->mobility;
	else if(IsNonEssential(Pearl))
		return dynamic_cast<NonEss *>(Pearl)->mobility;	
	else
	{
		cout << "What?" << endl;
		exit(0);
	}
}
