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
	toxins.clear();
	resistance.clear();
	resistance_lookup.clear();
}

/**********
		GENOME GENERATION / COPYING
			GenerateGenome = random genome generation (i.o.w. only at start of simulation)
			CloneGenome = clone genome from parent (or part of genome)
			CopyPartOfGenome = After making a new (empty genome), this copies part of a GIVEN genome to this genome
																				*********/
void Genome::GenerateGenome(int init_nr_args, int total_nr_args, int init_nr_coregenes, int init_nr_noncoding, float gene_cost, float genome_size_cost, float init_mob)
{
	gene_cost_ = gene_cost;
	genome_size_cost_ = genome_size_cost;
	tot_nr_ARGs = total_nr_args;

	bool V = FALSE;
	int i,j;

	Toxin* toxin; // int A;
	ARG* arg;
	Transposase* tra;
	Antitoxin* antitoxin;
	Core* core;
	Noncoding* nc;
	transformant=0;

	//printf("%p\n", genomesize_);	// With this you can check if pointer segfaults, which means that it doesnt exist
	// Set all copy-number and genome size thingies to 0
	genomesize_= 0;
	coregenes_=0;

	if(V) cout << "Making new Pearl list" << endl;
	StringOfPearls=new list<Pearl*>(); 		// Make a new (empty) list of Pearl pointers

	for(i=0;i<init_nr_args;i++)					// Make this par-file dependent later (LS: already done?)
	{
		if(genrand_real1()<0.5)																//LS: so one out of every 50 times someone gets a toxin? so every ind has on average init_nr_toxins/50 toxins? (so many none)
		{
			int type = genrand_int(0,init_nr_args-1);		// Counts from 0 because computers :)
			arg = new ARG(type);						// Let toxin point to a new Toxin
			arg->mobility = 0.0;															// LS: move initial mobility to parfile
			StringOfPearls->push_back(arg);  		     	// Add pointer to toxin to pearl-list LIJST.PUSH_BACK()
			genomesize_++;
		}
	}

	for(i=0;i<init_nr_coregenes;i++)
	{																																			//LS: everyone gets (some?) coregenes?
		core = new Core(i);					// Let core point to this Core gene
		core->mobility = 0.0;
		StringOfPearls->push_back(core);  	    	// Add pointer to core to list
		genomesize_++;
		coregenes_++;																													//LS: when making the coregenes also the minimum necessity is set? (somewhere else you say coregenes_ is the minimum of coregenes necessary)
	}
	for(i=0;i<init_nr_noncoding;i++)
	{
		nc = new Noncoding(-1);					// Default the "type" of a non-coding part corresponds to a random type of toxin/antitoxin
		nc->mobility = 0.80*genrand_real1();
		StringOfPearls->push_back(nc);  	    	// Add pointer to core to list
		genomesize_++;
	}
	for(i=0;i<1;i++)
	{
		if(genrand_real1() < 0.00)
		{
			tra = new Transposase(-1);					// Default the "type" of a non-coding part corresponds to a random type of toxin/antitoxin
			tra->mobility = 0.0;
			StringOfPearls->push_back(tra);  	    	// Add pointer to core to list
			genomesize_++;
		}
	}

	if(V) cout << "VOOR:" << endl;
	if(V) cout << ListContent(NULL) << endl;

	// Shuffling for lists isn't very efficient, but the 3 lines below do the trick
	   vector< Pearl* > Vec( StringOfPearls->begin(), StringOfPearls->end() );
	   random_shuffle( Vec.begin(), Vec.end() );
	   StringOfPearls->assign( Vec.begin(), Vec.end() );


	if(V) cout << "NA:" << endl;
    if(V) cout << ListContent(NULL) << endl;
	if(V)cout << endl << endl;
	viable = TRUE;
	Create_Gene_Lists(total_nr_args);			// This creates a look-up structure to quickly figure our resistance and stuff. :)
	CalculateCompStrength();
//	sleep(1);
	generation_ = 1;
}

void Genome::CloneGenome(const Genome *parent)
{
	//double mob_start = GetPearlMobility(parent->StringOfPearls->front());
	//double mob_end = GetPearlMobility(parent->StringOfPearls->back());

	bool V = FALSE;
	StringOfPearls=new list<Pearl*>();
	if(V) cout << "Trying to clone genome of size " << parent->StringOfPearls->size() << endl;
	if(V) cout << ListContent(parent->StringOfPearls) << endl;
	CopyPartOfGenome(parent->StringOfPearls->begin(),parent->StringOfPearls->end());		// Copies full genome to this genome

	genomesize_=parent->genomesize_;
	generation_=parent->generation_;
	tot_nr_ARGs=parent->tot_nr_ARGs;
	coregenes_=parent->coregenes_;
	genome_size_cost_=parent->genome_size_cost_;
	gene_cost_=parent->gene_cost_;
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
		StringOfPearls->push_back(pearl);
		i++;
	}
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

bool Genome::FetchTransposon(list<Pearl*> &Fragment, list<Pearl*> &PearlListTransp, int begin_scan, int end_scan, bool cut)			// For inversions and duplications of multiple genes
{
	bool V = false;
	// First, let's find a putative transposon
	iter it = Fragment.begin();	
	iter it_end = Fragment.begin();
	advance(it, begin_scan);
	advance(it_end, end_scan);



	// int LHS_pos = begin_scan;
	// int RHS_pos = begin_scan;
	// int T_posit = begin_scan;
	
	
	bool found_transposase = FALSE;
	int Tpos = 0;
	int pos = 0;
	vector<int> LHS_sites;
	vector<int> RHS_sites;
	
	while(it!=it_end)
	{
		if(IsTransposase(*it))
		{
			found_transposase = TRUE;	
		
		}
		else
		{
			double thismob = GetPearlMobility(*it);
			bool inverted_repeat = thismob > 0.75;
			if(inverted_repeat)
			{
				if(found_transposase)
				{
					RHS_sites.push_back(pos+begin_scan);
				}
				else
				{
					LHS_sites.push_back(pos+begin_scan);					
				}
			}
			if(!found_transposase) Tpos++;
		}

		pos++;
		it++;
	}
		
	if(!found_transposase) return false;
	if(LHS_sites.size() < 1) return false;
	if(RHS_sites.size() < 1) return false;
	
	if(V) cout << endl << "Finding transposable element on"  << ListContent(&Fragment, FALSE, FALSE, TRUE);
	if(V) cout << "Between " << begin_scan << " and " << end_scan << endl;
	if(V) cout << "Position of tranposase: " << Tpos+begin_scan << endl;	
	if(V) cout << "These are the inverted repeat sites:" << endl;
	if(V) for(int i = 0; i<LHS_sites.size();i++) cout << LHS_sites[i] << ",";
	if(V) cout << endl;
	if(V) for(int i = 0; i<RHS_sites.size();i++) cout << RHS_sites[i] << ",";
	if(V) cout << endl;
	if(V) cout << "Picking: " << LHS_sites[LHS_sites.size()-1] << " and " << RHS_sites[0] << " as sites." << endl;
	
	int start_pos = LHS_sites[LHS_sites.size()-1];
	int end_pos = RHS_sites[0];
	
	//V=false;
	if(V)cout << "Putative transposon area found at " << start_pos << " ending at pos " << end_pos << endl;			
	if(V) cout << "Start" << endl;
	// If it iss... we go on

	
	iter tr_start = Fragment.begin();
	advance(tr_start,start_pos);
	
	iter tr_end = Fragment.begin();
	advance(tr_end,end_pos);

	// if(!IsNoncoding(*tr_start)) cout << "!?!" << endl;
	// else cout << "Fine!" << endl;
	// if(!IsNoncoding(*tr_end)) cout << "!?!" << endl;
	// else cout << "Fine!" << endl;
	//V = true;
	if(V) cout << "Getting mobility" << endl;
	double mob_start = GetPearlMobility(*tr_start);	
	double mob_end = GetPearlMobility(*tr_end);
	// cout << "Mob start: " << mob_start << endl;
	// cout << "Mob end: " << mob_end << endl;
	
	mob_start = (mob_start-0.75)/0.25;
	mob_end = (mob_end-0.75)/0.25;
	// cout << "Mob start2: " << mob_start << endl;
	// cout << "Mob end2: " << mob_end << endl;
	if(V) cout << "Rolling the dice" << endl;
	double chance = mob_start*mob_end;
	double dice = genrand_real1();

	if(V) cout << "Dice: " << dice << endl;

	if(dice < chance)
	{
		
		if(V)cout << "Assembling transposon..." << endl;
		if(V)cout << "From this frag: " << ListContent(&Fragment,FALSE, FALSE, TRUE) << endl;
		iter it = tr_start;
		tr_end++;
		Pearl *pearl;
		while(it!=tr_end)
		{
			pearl=(*it)->clone();			
			PearlListTransp.push_back(pearl);
			if(cut) 
			{	
				delete(*it);	
				it = Fragment.erase(it);
			}
			else
				it++;
		}
		if(genrand_real1() < 0.0001) cout << "Fetched the following transposon: " << ListContent(&PearlListTransp,FALSE, FALSE, TRUE) << endl;
		return true;
	}
	else 
	{
		return false;
	}
}

void Genome::CalculateCompStrength()
{
			bool V = FALSE;
 				list<int> uniquecores;
				if(coregenes_ > 0)
				{
					uniquecores = coregenes;
					uniquecores.sort();
					uniquecores.unique();
					if(V)
					{
						cout << "Found " << uniquecores.size() << " core genes. Need " << coregenes_ << endl;		// For clarity: coregenes_ is just a number to check the number of genes the species should minimally have, while without the hyphen (coregenes) is actually the list of all coregenes in that genome
					}
				}

				if(StringOfPearls->size() == 0)
				{
					compstrength = 1.0;
				}
				else if(coregenes_ == uniquecores.size())			// A) Check if all core genes exist in genome LS:!! nu mag je niet meer unieke coregenes hebben dan het minimum noodzakelijk?!?
				{
					if(V) cout << "Determine fitness because all core genes are present" << endl;					
					//compstrength = 1.0-(resistance.size()+toxins.size()+coregenes.size())*gene_cost_;
					//compstrength = 1.0-(coregenes.size()+resistance.size())*gene_cost_;
					compstrength = 1.0-(coregenes.size()+transposases.size()+resistance.size())*gene_cost_;
					if(V) cout << "Fitnessl lowered with " << (coregenes.size()+resistance.size())*gene_cost_  << endl;
					compstrength -= (StringOfPearls->size()*genome_size_cost_);
					if(V) cout << "Fitness l furthermore " << (StringOfPearls->size()*genome_size_cost_)  << endl;
					compstrength = max(compstrength,0.0);
				}
				else compstrength = 0.0;
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
					// if (V) cout <<"Mutates mobility with chance"<< gene_mob << endl;
					MutateMobility(i), mutated=TRUE;
					// if (V) cout << " Done mutating mobility" <<endl;
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
		

		
		Create_Gene_Lists(tot_nr_ARGs);
		// DUPLICATING VIA CHUNKS (analogous to via extracellular DNA)
		




	bool discovered_gene = GeneDiscovery(gendisc);
	if(discovered_gene) mutated=TRUE;

	//if(V) cout << mutated << endl;
	//if(V) if(mutated) cout << this->ListContent(NULL) << endl << endl; // This is more or less arbitrary.. sometimes it prints the genome of a child

	return mutated;
}


bool Genome::TransposonDynamics()
{
	/// TRANSPOSON EITHER JUMPING OR DUPLICATING THEMSELVES
	
	bool V = false;
	
	if(V) cout << "\n\nTRANSPOSON DYNAMIC START\n\nScanning for transposable element from " <<  endl;
	bool mutated = false;
	Create_Gene_Lists(tot_nr_ARGs);

		int start = 0;
		int end = 5; //genrand_int(5,5); 			// NOTE HARD CODED, ALSO CHANGE IN DEGRADE FUNCTIno!		
		int max_num_jumps = 5;
		int num_jumps = 0;
		if(end > StringOfPearls->size()) end = StringOfPearls->size();		
		while(end < StringOfPearls->size())
		{
			bool cut = genrand_real1() < 0.5;
			bool paste = genrand_real1() < 0.5;
			if(V)cout << "-----------------------------------------------------" << endl;
			if(V)cout << "Before int" << endl;
			if(V)cout << "Cut?" << (cut) << endl;
			if(V)cout << "Paste?" << (paste) << endl;
			if(V)cout << ListContent(NULL) << endl;
			
			bool integrated = IntegrateDNA(*StringOfPearls,tot_nr_ARGs,start,end,cut,paste); //FetchTransposon(PearlList_ToIntegrate, Transposon, scan_start, scan_end, true);
			if(integrated) mutated = true;
			genomesize_ = StringOfPearls->size();
			
			if(V)cout << "After int" << endl;
			if(V)cout << ListContent(NULL) << endl;
			start+=(end-start);
			end=start+5; // genrand_int(1,5);			
			if(end > StringOfPearls->size()) end = StringOfPearls->size();
			num_jumps++;
			if(num_jumps > max_num_jumps) 
			{
				compstrength = 0.0;
				break;				
			}
		}


	if(mutated) 
	{
		Create_Gene_Lists(tot_nr_ARGs);
		return true;
	}			

	return false;
	// DUPLICATING VIA CHUNKS (analogous to via extracellular DNA)

}



void Genome::MutateMobility(iter ii)
{
	if(IsNoncoding(*ii))
		dynamic_cast<Noncoding *>(*ii)->mobility = min(1.0,max(0.0,( (dynamic_cast<Noncoding *>(*ii)->mobility) + (2.0*genrand_real1()-1.0) * 0.2)));
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
	Toxin * tox;
	Antitoxin * antitox;
	Core * core;
	int type = -1;		// Nonecoding type heeft altijd -1. Deze moet nooit gebruikt worden voor "echte" genen, maar het is dus good practice om deze op -1 te zetten zodat het absoluut fout gaat als het perongeluk in een normaal gen overerft
	double mob = GetPearlMobility(*ii);
	//if(V) cout << "Inheriting mobility"
	if(V)cout << "Befor: " << ListContent(NULL) << endl;
	nc = new Noncoding(type);	 // Make new noncoding pearl with a random type (random with respect to our toxin/antitoxins)
	nc->mobility = genrand_real1(); 	 	// Inherit mobility from original gene
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
	Toxin * tox;
	ARG * arg;
	Antitoxin * antitox;
	Transposase * tra;
	Core * core;
	// RESISTANCE GENE DISCOVERY
	if(genrand_real1() < 0.5)
	{
		int type = genrand_int(0,tot_nr_ARGs-1);		// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY, this respresents our scope of "possible" types
		double mob = GetPearlMobility(*ii);
		if(V)cout << "Befor: " << ListContent(NULL) << endl;
		double randnum = genrand_real1();	
		arg = new ARG(type);	 				// Make new arg pearl with type of the non-coding gene.. do we want this? Dunno
	//	arg->mobility = mob; 	 				// Inherit mobility from original non-coding sequence. This we do want.. i guess.
		arg->mobility = 0.0; 	 				// For now, only NC "flanks" can infer gene mobility
		delete *ii;						 		// Delete the original non-coding gene
		ii = StringOfPearls ->erase(ii); 		// Erase from list
		ii = StringOfPearls->insert(ii,arg);	// Insert new arg	
	}
	else
	{
		int type = genrand_int(0,tot_nr_ARGs-1);		// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY, this respresents our scope of "possible" types
		double mob = GetPearlMobility(*ii);
		if(V)cout << "Befor: " << ListContent(NULL) << endl;
		double randnum = genrand_real1();	
		tra = new Transposase(type);	 				// Make new arg pearl with type of the non-coding gene.. do we want this? Dunno
	//	arg->mobility = mob; 	 				// Inherit mobility from original non-coding sequence. This we do want.. i guess.
		tra->mobility = 0.0; 	 				// For now, only NC "flanks" can infer gene mobility
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

	//if(V) cout << "GeneDiscovery called" << endl;
	if(genrand_real1()<gendisc)
	{
		if(V) cout << ListContent(NULL) << endl;
		if(V) cout << "Discovering resistance with chance " << gendisc << endl;
		ARG *arg;
		ii = StringOfPearls->begin();
		randpos = (int)(genrand_real1()*genomesize_);
		if(V) cout << "Randpos: " << randpos << " van de " << genomesize_ << endl;
		advance(ii,randpos);
		type = genrand_int(0,tot_nr_ARGs-1);		// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY
		if(V) cout << "Type: " << type << " van de " << tot_nr_ARGs << endl;
		arg = new ARG(type);						// Let toxin point to a new Toxin
        arg->mobility = 0.0;	// Newly introduced 
		ii = StringOfPearls->insert(ii,arg);  		     	// Add pointer to toxin to pearl-list
		if(V) cout << ListContent(NULL) << endl;
		genomesize_++;
		mutated=TRUE;
	}	
		//if(V) cout << "GeneDiscovery called" << endl;
	if(genrand_real1()<gendisc)
	{
		if(V) cout << ListContent(NULL) << endl;
		if(V) cout << "Discovering nc with chance " << gendisc << endl;
		Noncoding *nc;
		ii = StringOfPearls->begin();
		randpos = (int)(genrand_real1()*genomesize_);
		if(V) cout << "Randpos: " << randpos << " van de " << genomesize_ << endl;
		advance(ii,randpos);
		type = genrand_int(0,tot_nr_ARGs-1);		// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY
		if(V) cout << "Type: " << type << " van de " << tot_nr_ARGs << endl;
		nc = new Noncoding(type);						// Let toxin point to a new Toxin
        nc->mobility = genrand_real1();	// Newly introduced 
		ii = StringOfPearls->insert(ii,nc);  		     	// Add pointer to toxin to pearl-list
		if(V) cout << ListContent(NULL) << endl;
		genomesize_++;
		mutated=TRUE;
	}	

	//if(V) cout << "GeneDiscovery called" << endl;
	if(genrand_real1()<gendisc)
	{
		if(V) cout << ListContent(NULL) << endl;
		if(V) cout << "Discovering tra with chance " << gendisc << endl;
		Transposase *tra;
		ii = StringOfPearls->begin();
		randpos = (int)(genrand_real1()*genomesize_);
		if(V) cout << "Randpos: " << randpos << " van de " << genomesize_ << endl;
		advance(ii,randpos);
		type = genrand_int(0,tot_nr_ARGs-1);		// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY
		
		// UNCOMMENT BELOW TO MAKE TRANSPASE-DISCOVERY INCLUDING IRS
		// Noncoding *nc;
		// nc = new Noncoding(type);						// Let toxin point to a new Toxin
        // nc->mobility = 0.85+genrand_real1()*0.15;	// Newly introduced 
		// ii = StringOfPearls->insert(ii,nc);  		     	// Add pointer to toxin to pearl-list
		// if(genrand_real1() < 0.5)
		// {
		// 	nc = new Noncoding(type);						// Let toxin point to a new Toxin
        // 	nc->mobility = 0.0;	// Newly introduced 
		// 	ii = StringOfPearls->insert(ii,nc);  		     	// Add pointer to toxin to pearl-list
		// }
		
		if(V) cout << "Type: " << type << " van de " << tot_nr_ARGs << endl;
		tra = new Transposase(type);						// Let toxin point to a new Toxin
        tra->mobility = 0.0;	// Newly introduced 
		ii = StringOfPearls->insert(ii,tra);  		     	// Add pointer to toxin to pearl-list

		// UNCOMMENT BELOW TO MAKE TRANSPASE-DISCOVERY INCLUDING IRS
		// if(genrand_real1() < 0.5)
		// {
		// 	nc = new Noncoding(type);						// Let toxin point to a new Toxin
        // 	nc->mobility = 0.0;	// Newly introduced 
		// 	ii = StringOfPearls->insert(ii,nc);  		     	// Add pointer to toxin to pearl-list
		// }

		// type = genrand_int(0,tot_nr_ARGs-1);				// Counts from 0 because computers :) USING RESISTANCE SIZE = OKAY
		// nc = new Noncoding(type);							// Let toxin point to a new Toxin
        // nc->mobility = 0.85+genrand_real1()*0.15;			// Newly introduced 
		// ii = StringOfPearls->insert(ii,nc);  		     	// Add pointer to toxin to pearl-list

		//cout << ListContent(NULL) << endl;
		genomesize_++;
		mutated=TRUE;
	}
	
	return mutated;
}

bool Genome::IntegrateDNA(list<Pearl*> &PearlList_ToIntegrate, int total_nr_args, int scan_start, int scan_end, bool cut, bool paste)
{
	
			
		

	bool V = false;
  	if(V)
	{
		cout << "Startintegrating:" << endl << ListContent(NULL, FALSE, FALSE, FALSE) << endl;
	}
	bool transposon_found = false;
	std::list<Pearl*> Transposon;
	
	int tries = 1;
	while(!transposon_found)
	{	
		if(V) cout << "... scanning for transposable element from " << scan_start << " to " << scan_end << endl;  	
		transposon_found = FetchTransposon(PearlList_ToIntegrate, Transposon, scan_start, scan_end, cut); // cut = true to remove transposon from original fragment
		if(V)cout << "Found? --> " << transposon_found << endl;
		tries--;
		if(tries==0) break;
	}

	if(V) {
		 if(transposon_found) cout << "\n\n\n\n\n\n TRANSPOSON FOUND !!!! \n\n\n\n\n\n" << endl;
		 else cout << "No transposon found...." << endl;
	}
	if(transposon_found && paste)
	{
			
		if(V)
		{
			cout << "This is the transposon that was found to jump: " << endl;
			cout << ListContent(&Transposon, FALSE, FALSE, FALSE) << endl;
			//cout << ListContent(&PearlList_ToIntegrate) << endl;
		}

		Pearl *pearl;

		// let's determine where in the GENOME the transposon will insert
		iter ii = StringOfPearls->begin();
		int randpos = (int)(genrand_real1()*(StringOfPearls->size()-1));
		advance(ii,randpos);


		// Where it inserts, the existing gene is broken (delete, replace with two flanks that are non-coding)
		delete *ii;
		ii = StringOfPearls->erase(ii);

		// Insert NC (left side of the broken gene)
		Noncoding *nc;		
		nc = new Noncoding(0);						// Let 
        nc->mobility = genrand_real1();					// Newly introduced 
		ii = StringOfPearls->insert(ii,nc);

		if(V)
		{
			cout << "Integrating piece of DNA of " << Transposon.size() << " long " << endl;
			cout << "Into position: " << randpos << " of " << endl << ListContent(NULL, FALSE, FALSE, FALSE) << endl;
		}

		// Insert the transposon
		iter it = Transposon.begin();

		while(it!=Transposon.end())
		{
			pearl=(*it);
			if(V) cout << "Inserting bead to DNA.." << endl;
			StringOfPearls->insert(ii, pearl);
			
			if(V) cout << "Deleting bead from transposon" << endl;
			if(V) cout << "Done." << endl;
			it = Transposon.erase(it);
		}

		// Insert NC (right side of the broken gene)
		
		nc = new Noncoding(0);						// Let 
        nc->mobility = genrand_real1();					// Newly introduced 
		ii = StringOfPearls->insert(ii,nc);

		if(V) cout << "Done with splice." << endl;
		
		genomesize_ = StringOfPearls->size();
		Create_Gene_Lists(total_nr_args);			// This creates a look-up structure to quickly figure our resistance and stuff. :)
		CalculateCompStrength();
		if(V)cout << ListContent(NULL, FALSE, FALSE, FALSE) << endl;
		return TRUE;
	}
	else
	{
		if(V) if(!transposon_found) cout << "Not integrated because no transposon found " << endl;
		if(V) if(!paste) cout << "Not integrated because it didnt paste" << endl;
		return FALSE;
	}

}

void Genome::Create_Gene_Lists(int total_nr_args)
{
	
	/*
	This function generates three things:
	1) A vector ('array') of size TOT_NR_TOXINS that holds whether this individual is resistance to toxins
	2) A list of toxins this individual produces.
	3) A list of household genes this individual produces
	*/

	bool V = FALSE;

	list<Pearl*>* StringOfPearls = this->StringOfPearls;
	iter i = StringOfPearls->begin();
	Toxin* tox;
	Core* core;
	Antitoxin* antitox;
	ARG* arg;
	Transposase* tra;

	if(V) cout << "Making lookup table" << endl;
	resistance_lookup.assign(total_nr_args, 0);	// Initialise as 0's

	if(V) cout << endl;

	while(i!=StringOfPearls->end())
	{
		if(IsToxin(*i))
		{
			tox=dynamic_cast<Toxin *>(*i);
			toxins.push_back(tox->type);
		}
		if(IsARG(*i))
		{
			arg=dynamic_cast<ARG *>(*i);
			resistance.push_back(arg->type);
			resistance_lookup[arg->type] += 1;
		}
		if(IsCore(*i))
		{
			core=dynamic_cast<Core *>(*i);
			coregenes.push_back(core->type);
		}
		if(IsTransposase(*i))
		{
			tra=dynamic_cast<Transposase *>(*i);
			transposases.push_back(tra->type);
		}
		i++;
	}
	if(V)
	{
		for (vector<int>::const_iterator i = resistance_lookup.begin(); i != resistance_lookup.end(); ++i)
			cout << *i << ' ';
		cout << endl;
		for (list<int>::const_iterator i = toxins.begin(); i != toxins.end(); ++i)
			cout << *i << ' ';
		for (list<int>::const_iterator i = resistance.begin(); i != resistance.end(); ++i)
			cout << *i << ' ';
		cout << endl;
	}

}



string Genome::ListContent(list<Pearl*> *StringOfPearls, bool ignore_core, bool ignore_noncoding, bool includemobility) // If *Pearl is not given, it prints the entire genome
{
	bool V = FALSE;
	string GenomeContent;
	if(StringOfPearls == NULL) StringOfPearls = this->StringOfPearls;
	iter i = StringOfPearls->begin();
	Toxin *tox;
	Transposase *tra;
	Antitoxin *antitox;
	ARG *arg;
	Core *core;
	Noncoding * nc;

	if(V)
		cout << "Attempting to list genome with length " << StringOfPearls->size() << endl;
	GenomeContent+="\033[00m---";																					//LS: klopt Dit?
	while(i!=StringOfPearls->end())
	{
		if(ignore_noncoding)
			if(IsNoncoding(*i)) {i++;continue;}
		if(ignore_core)
			if(IsCore(*i)) {i++;continue;}
		stringstream stringtemp;
		stringtemp.precision(2);
		stringtemp << fixed;
		if(i!=StringOfPearls->begin()) GenomeContent +="-";
		if(IsToxin(*i))
		{
			tox=dynamic_cast<Toxin *>(*i);
			stringtemp << "\033[1;41m" << tox->type;
			if(includemobility) stringtemp << ":" << tox->mobility;
			stringtemp << "\033[00m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsARG(*i))
		{
			arg=dynamic_cast<ARG*> (*i);
			stringtemp << "\033[1;44mA" << arg->type;
			if(includemobility) stringtemp <<  ":" <<arg->mobility;
			stringtemp << "\033[00m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsTransposase(*i))
		{
			tra=dynamic_cast<Transposase*> (*i);
			stringtemp << "\033[1;41mT";
			if(includemobility) stringtemp <<  ":" <<tra->mobility;
			stringtemp << "\033[00m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsAntitoxin(*i))
		{
			antitox=dynamic_cast<Antitoxin*> (*i);
			stringtemp << "\033[1;44m" << antitox->type;
			if(includemobility) stringtemp <<  ":" <<antitox->mobility;
			stringtemp << "\033[00m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsCore(*i))
		{
			core=dynamic_cast<Core *>(*i);
			stringtemp << "\033[1;42mH" << core->type;
			if(includemobility) stringtemp << ":" <<core->mobility;
			stringtemp << "\033[00m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsNoncoding(*i))
		{
			nc=dynamic_cast<Noncoding *>(*i);
			if(nc->mobility < 0.75) 
			{
				stringtemp << "\033[1;90m";
				if(includemobility) stringtemp << "" << nc->mobility;
				else stringtemp << "C";
				stringtemp << "\033[00m";
			}			
			else if(nc->mobility < 0.875)
			{
				stringtemp << "\033[1;43m";
				if(includemobility) stringtemp << "" << nc->mobility;
				else stringtemp << "r";
				stringtemp << "\033[00m";
			}
			else
			{
				stringtemp << "\033[1;43m";
				if(includemobility) stringtemp << "" << nc->mobility;
				else stringtemp << "R";
				stringtemp << "\033[00m";
			}
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else{
			GenomeContent+="He?";
		}
		i++;
	}
	GenomeContent+="\033[49m---\n";
	return GenomeContent;
}

string Genome::ListMobility(list<Pearl*> *StringOfPearls) // If *Pearl is not given, it prints the entire genome
{
	bool V = FALSE;
	string GenomeContent;
	if(StringOfPearls == NULL) StringOfPearls = this->StringOfPearls;
	iter i = StringOfPearls->begin();

	Toxin *tox;
	Antitoxin *antitox;
	ARG *arg;
	Core *core;
	Noncoding *nc;
	if(V)
		cout << "Attempting to list genome with length " << StringOfPearls->size() << endl;
	GenomeContent+="\033[0m---";
	while(i!=StringOfPearls->end())
	{
		stringstream stringtemp;
		if(i!=StringOfPearls->begin()) GenomeContent +="-";
		if(IsToxin(*i))
		{
			tox=dynamic_cast<Toxin *>(*i);
			stringtemp << "\033[1;31m" << tox->mobility << "\033[0m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsARG(*i))
		{
			arg=dynamic_cast<ARG *>(*i);
			stringtemp << "\033[1;34m" << arg->mobility << "\033[0m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsAntitoxin(*i))
		{
			antitox=dynamic_cast<Antitoxin *>(*i);
			stringtemp << "\033[1;34m" << antitox->mobility << "\033[0m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsCore(*i))
		{
			core=dynamic_cast<Core *>(*i);
			stringtemp << "\033[1;33m" << core->mobility << "\033[0m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsNoncoding(*i))
		{
			nc=dynamic_cast<Noncoding *>(*i);
			stringtemp << "\033[1;90m" << nc->mobility << "\033[0m";
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		i++;
	}
	GenomeContent+="\033[0m---";
	//cout << GenomeContent << endl;
	return GenomeContent;
}

/*********
							Get genome stat functions
																	********/
int Genome::GetTotalNrToxins(int type)
{
	int total_number = 0;
	for (list<int>::const_iterator i = toxins.begin(); i != toxins.end(); ++i)
			if(*i == type)
				total_number++;
	return total_number;
}

int Genome::GetTotalNrDefense(int type)
{
	return resistance_lookup[type];
}

double Genome::GetFractionTEs()
{	
	if(StringOfPearls == NULL) StringOfPearls = this->StringOfPearls;
	iter i = StringOfPearls->begin();
	int num_mobile_TEs = 0;
	int num_genes = 0;
	while(i!=StringOfPearls->end())
	{
		num_genes++;
		if(IsNoncoding(*i))
		{
			if(GetPearlMobility(*i) > 0.75) num_mobile_TEs++;
		}
		i++;
	}
	return num_mobile_TEs;
}

// Below are all the functions that check whether a Pearl is of a certain type
bool Genome::IsToxin(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(Toxin));	// typeid determines the class of an object at runtime
}

bool Genome::IsTransposase(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(Transposase));	// typeid determines the class of an object at runtime
}

bool Genome::IsAntitoxin(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(Antitoxin));	// typeid determines the class of an object at runtime
}

bool Genome::IsARG(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(ARG));	// typeid determines the class of an object at runtime
}


bool Genome::IsCore(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(Core));	// typeid determines the class of an object at runtime
}

bool Genome::IsNoncoding(Pearl *Pearl) const
{
	return (bool)(typeid(*Pearl) == typeid(Noncoding));	// typeid determines the class of an object at runtime
}

int Genome::GetPearlType(Pearl *Pearl) const
{
	if(IsToxin(Pearl))
		return dynamic_cast<Toxin *>(Pearl)->type;
	else if(IsAntitoxin(Pearl))
		return dynamic_cast<Antitoxin *>(Pearl)->type;
	else if(IsARG(Pearl))
		return dynamic_cast<ARG *>(Pearl)->type;
	else if(IsCore(Pearl))
		return dynamic_cast<Core *>(Pearl)->type;
	else if(IsNoncoding(Pearl))
		return dynamic_cast<Noncoding *>(Pearl)->type;
}

double Genome::GetPearlMobility(Pearl *Pearl) const
{
	//cout << IsToxin(Pearl) << " " << IsAntitoxin(Pearl) << " " << IsCore(Pearl) << " "<<  IsNoncoding(Pearl) << endl;
	if(IsToxin(Pearl))
		return dynamic_cast<Toxin *>(Pearl)->mobility;
	else if(IsAntitoxin(Pearl))
		return dynamic_cast<Antitoxin *>(Pearl)->mobility;
	else if(IsTransposase(Pearl))
		return dynamic_cast<Transposase *>(Pearl)->mobility;
	else if(IsARG(Pearl))
		return dynamic_cast<ARG *>(Pearl)->mobility;
	else if(IsCore(Pearl))
		return dynamic_cast<Core *>(Pearl)->mobility;
	else if(IsNoncoding(Pearl))
		return dynamic_cast<Noncoding *>(Pearl)->mobility;
	else
	{
		cout << "What?" << endl;
		exit(0);
	}
}
/*
void Genome::Increment_HGT_num(Pearl *Pearl) const
{
	//cout << IsToxin(Pearl) << " " << IsAntitoxin(Pearl) << " " << IsCore(Pearl) << " "<<  IsNoncoding(Pearl) << endl;
	if(IsToxin(Pearl))
		dynamic_cast<Toxin *>(Pearl)->num_horizontal_transfers_++;
	else if(IsAntitoxin(Pearl))
		dynamic_cast<Antitoxin *>(Pearl)->num_horizontal_transfers_++;
	else if(IsCore(Pearl))
		dynamic_cast<Core *>(Pearl)->num_horizontal_transfers_++;
	else if(IsNoncoding(Pearl))
		dynamic_cast<Noncoding *>(Pearl)->num_horizontal_transfers_++;
}
*/
