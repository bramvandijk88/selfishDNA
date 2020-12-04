#ifndef FUNCTIONS_DATA_HH
#define FUNCTIONS_DATA_HH

#include "Brem-cash/cash2-s.hpp"
#include "functions_data.hpp"



extern int nrow;
extern int ncol;
extern string folder;

extern float birth, birthNON;
extern float gene_loss, gene_dupl, tandem_dupl, tandem_del, inversions, gene_to_noncoding, noncoding_to_gene, mut_rate_scaling, ab_penalty, gene_mob, gene_discovery, init_mob;

extern float gene_cost, genome_size_cost;

extern int total_nr_args;
extern vector<double> ABs_concentration;

extern float HGTcap = 1.0;






/* -------------------------------------------------------------------------------------*/
/* ---------------------------- BELOW = FUNCTION PROTOCOLS -----------------------------*/
/* -------------------------------------------------------------------------------------*/


void InitialiseGenomes(TYPE2** vibrios, TYPE2** e, int init_nr_args, int total_nr_args, int init_nr_coregenes, int init_nr_noncoding)
{
	bool V=FALSE;
	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			if(vibrios[row][col].val > 1)
			{
				if(V) cout << "Initialising genome for individual at " << row << " " << col << endl;					//LS: using " " just to space the numbers?
				vibrios[row][col].G = new Genome();
				vibrios[row][col].G->GenerateGenome(init_nr_args, total_nr_args, init_nr_coregenes, init_nr_noncoding, gene_cost, genome_size_cost, init_mob);
				vibrios[row][col].G->generation_ = 1;
				vibrios[row][col].val = genrand_int(52,255);
				if(V) cout << vibrios[row][col].G->ListContent(NULL) << endl;
			}
			else
			{
				if(V) cout << "\n\nVal is " << vibrios[row][col].val << " so not making genome\n\n";
			}
		}
	}
}

void InitialiseDNAPlane(TYPE2 **DNA)
{
	for (int row=1; row <= nrow; row++)
		for (int col=1; col <= ncol; col++)
			DNA[row][col].DNA = new eDNA();		// Make objects
}

void RefreshEnvironment(TYPE2** environment)
{
	bool V=TRUE;
	// PLACE
	int offsetx = genrand_int(-nrow/2,nrow/2);
	int offsety = genrand_int(-ncol/2,ncol/2);
	int width = genrand_int(5000,5000);
	int resource_density = genrand_int(300,300);
	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			//De if hieronder tekent een cirkel
			environment[row][col].val2 = 300;
			/*if( (pow((col-ncol/2),2) + pow((row-nrow/2),2) ) < width)
			{
				int rowdraw, coldraw;
				rowdraw = row+offsetx+nrow;
				coldraw = col+offsety+ncol;
				//if(environment[rowdraw % nrow][coldraw % ncol].val2 == 0)
				environment[rowdraw % nrow+1][coldraw % ncol+1].val2 = max(resource_density,environment[rowdraw % nrow+1][coldraw % ncol+1].val2);  // Val2 is the ACTUAL concentration of resource
			}*/
		}
	}
}

int GetSpeciesColour(TYPE2** vibrios, int row, int col){
	// Dit is eens HASH functie, die de tekst-representatie van species omzet naar een (vaak) unieke int van 2 - 241, wat dan in CASH als kleur wordt gebruikt
	string species = vibrios[row][col].G->ListContent(NULL, FALSE, TRUE, FALSE);
	unsigned int hash = 0;
	for(int i = 0; i < species.length(); i++){
		hash += (int)species[i];
		hash += (hash << 10);
		hash ^= (hash >> 6);
	}
	hash += (hash << 3);
	hash ^= (hash >> 11);
	hash += (hash << 15);
	hash = ((int)(hash%160))+90;
	species.clear();
	return hash;
}

void CalculateABEffect(TYPE2 **vibrios)
{
	/*
	This function calculates fitness and death rate
	*/
	bool V=FALSE;

	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			if(vibrios[row][col].val > 0)                    	//LS: also assigning fitness to dying individuals
			{
				double sum_killing = 0.0;
				for(int i=0;i<total_nr_args;i++)
					if(vibrios[row][col].G->resistance_lookup[i] == 0 && ABs_concentration[i] > 0) sum_killing += ABs_concentration[i]*0.02;
				vibrios[row][col].ab_effect = sum_killing;
			}
			else
			{
				vibrios[row][col].fitness = 0.0;											//LS: why necessary to assign an empty grid with a zero fitness?
			}
		}
	}

}

// Spill DNA chops the genome up in pieces and drops it onto the DNA plane
void SpillDNA(TYPE2 **vibrios, TYPE2** DNA, int row, int col)
{
	if(vibrios[row][col].G->StringOfPearls->size() > 0)
		DNA[row][col].DNA->FragmentiseGenome(vibrios[row][col].G,1.0);
}

// Every fragment has a chance to disappear --- LS: shouldn't this be dependent on their length?
void DegradateDNA(TYPE2** DNA, int row, int col, float degr)
{
	DNA[row][col].DNA->Degradate(degr);
}

// Applies HGT for an invididual
void DoHGT(TYPE2 **vibrios, TYPE2** DNA, int total_nr_toxins)
{
	bool V=FALSE;
	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			if(vibrios[row][col].val>1)
			{
				//Add/change: Loop over fragment list to give everyone a chance to integrate
				eDNA::fr_iter frits;
				//cout << DNA[row][col].DNA->Fragments->size() << endl;
				frits = DNA[row][col].DNA->Fragments->begin();
				while(frits != DNA[row][col].DNA->Fragments->end())
 				{

						//if(V)cout << "Local poolsize [" << row << "][" << col << "] " << DNA[row][col].DNA->Fragments->size() << endl;
						if (genrand_real1() < 0.05) 																						//LS: fixed chance of uptake
						{
							bool mutated = false;
							stringstream old;
							stringstream frag_old;
							stringstream frag_new;

							//V = FALSE;
							bool transposase = false;
								
							frag_old << vibrios[row][col].G->ListContent(&(*frits),FALSE,FALSE,FALSE) << endl;
							old << vibrios[row][col].G->ListContent(NULL,FALSE,FALSE,FALSE) << endl; 
							//if(vibrios[row][col].G->transposases.size() > 0) transposase = true;		
							Genome::iter i = frits->begin();					
							while(i!=frits->end()) { if(vibrios[row][col].G->IsTransposase(*i)) transposase = true; i++; }
							
							if(V)
							{
								cout << "Trying to HGT into genome: " << endl;
								cout << vibrios[row][col].G->ListContent(NULL) << endl << endl;								
								cout << "Has transposase? --> " << transposase << endl;
								cout << DNA[row][col].DNA->Fragments->size() << endl;
								cout << "Fragment inserted:" << endl;
								cout << vibrios[row][col].G->ListContent(&(*frits)) << endl; 
							}								
							if(transposase) mutated = vibrios[row][col].G->IntegrateDNA(*frits, total_nr_args, 0, frits->size()-1,true,true); 	
							frag_new << vibrios[row][col].G->ListContent(&(*frits),FALSE,FALSE,TRUE) << endl;
							//else cout << "Not integrating because no transposon" << endl;
							if(V)cout << "Did it integrate? --> " << mutated << endl;
							
							//if(mutated) V = TRUE;																
								
							//DNA[row][col].DNA->Fragments->erase(frits);
							frits = DNA[row][col].DNA->DeleteFragment(frits);
							if (DNA[row][col].DNA->Fragments->size() == 0) break;
							if(mutated) 
							{			
									V = false;																
									if(V)cout << "Done with HGT: " << endl ;
									if(V)cout << "OLD: " << old.str() << endl;
									if(V)cout << "frag_old: " << frag_old.str() << endl;
									if(V)cout << "frag_new: " << frag_new.str() << endl;
									if(V)cout << "NEW: " << vibrios[row][col].G->ListContent(NULL,FALSE,FALSE,FALSE) << endl << endl<< endl;	
									vibrios[row][col].G->Create_Gene_Lists(total_nr_args);
									vibrios[row][col].G->transformant = 3;
									V = false;
							}
						}
						else
						{
							frits++;
						}

					}
				}
			}
		}
}


void CleanUpTheDead(TYPE2 **vibrios, TYPE2 **edna)
{
	/*
	This function calculates fitness and death rate --LS: no this cleans up the dead
	*/
	bool V=FALSE;

	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			if(vibrios[row][col].val == 1)				// Val 1 are the cells that died last update
			{
				if(V) cout << "Cleaning up " << row << " " << col <<  endl;
				vibrios[row][col].fitness = 0.0;
				delete vibrios[row][col].G;
				vibrios[row][col].val = 0;				// Val 1 = prepared to die, val 0 = actually empty spot
			}
			//else
				//if(V) cout << "Not cleaning up " << vibrios[row][col].val << endl;
		}
	}
}

void CompeteAndReproduce(TYPE2 **vibrios, TYPE2 **e, int row, int col)
{

	bool V = FALSE;

	TYPE2* competitor;												// Placeholder for invading cell
	double rand = genrand_real1();									// Roulette wheel number
	double pinvade = 0.0;											// Probability for current competitor to invade (cumulative)
	double fitsum = birthNON;										// Sum of all fitnessvalues for competitors (LS: initialised with birthNON, adds up in loop)
	for(int i=1;i<9;i++)											// 1 - 8 are all surrounding cells in the Moore neighborhood (0 is self)
	{
		competitor = GetNeighborP(vibrios,row,col,i);
		if(competitor->val > 1) fitsum += competitor->G->compstrength;
	}
	if(V) cout << endl << "Starting competition selection" << endl;
	if(V) cout << "rand is " << rand << endl;
	if(V) cout << "fitsum is " << fitsum << endl;
	for(int i=1;i<9;i++)						// 1 - 9 are all surrounding cells in the Moore neighborhood (0 is self)
	{
		competitor = GetNeighborP(vibrios,row,col,i);			// Store neightbour i in temporary placeholder
		if(competitor->val > 1) 																			// LS: arrow because only the pointer to the individual is stored?
		{
			pinvade += (competitor->G->compstrength)/(fitsum+birthNON);
			if(V) cout << "pinvade for competitor " << i << " is "  << pinvade << endl; 				// LS: the output is not the actual chance of invading for that individual, but the cumulative chance
			if(rand <= pinvade)																																	// LS: sure that this cyan thing doesnt chance the properties of the parameter?
			{
					if(V) cout << "!! invaded !! " << competitor->val <<  endl;
					
					if(V) cout << "generation " << competitor->G->generation_ << endl;
					vibrios[row][col].val = competitor->val;		// Inherit val (in this case, this is always 1?)
					//environment[row][col].val2-=3;			// Extra bit of resource spent for this reproduction event
					//environment[row][col].val2 = max(environment[row][col].val2,0);
					// Copy genome to new individual
					vibrios[row][col].G = new Genome();
					vibrios[row][col].G->CloneGenome(competitor->G);
					bool mutated = vibrios[row][col].G->MutateGenome(gene_mob, gene_loss, gene_dupl, tandem_dupl, tandem_del, inversions, gene_to_noncoding, noncoding_to_gene, gene_discovery); // LS: why ' mutated' in Genome.cc but 'mutate' here
					if(V)cout << "1" << endl;
					bool integrated = vibrios[row][col].G->TransposonDynamics();
					if(integrated) 
					{
						vibrios[row][col].G->transformant = 2;	
					}
					//bool integrated = false;
					
					///bool integrated= false;
					if(mutated || integrated) vibrios[row][col].G->Create_Gene_Lists(total_nr_args);

					if(mutated || integrated) vibrios[row][col].G->CalculateCompStrength();
					else
						vibrios[row][col].G->compstrength = competitor->G->compstrength;							//fitness = parent fitness

					if(V) cout << "Genome of kid is: " << endl;
					if(V) cout << vibrios[row][col].G->ListContent(NULL) << endl;
					if(V) cout << "Generation of kid is: " << vibrios[row][col].G->generation_ << endl;
					if(V) cout << "Genomesize of kid is: " << vibrios[row][col].G->genomesize_ << endl;

					break;
			}
		}
	}
}

void AgeGenes(TYPE2 **vibrios)
{
	/*
	This function ages genes and resets  [HAS TO BE ADDED!]after transitions have been calcutated
	*/
	bool V=FALSE;

	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			if(vibrios[row][col].val > 0)                    	//LS: do we need dying individuals to have their age increased?
			{
				list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
				Genome::iter i = StringOfPearls->begin();
				while(i!=StringOfPearls->end())
				{
					if(vibrios[row][col].G->IsToxin(*i))
					{
						if (V) cout << "found toxin with age "<< dynamic_cast<Toxin*>(*i)->gene_age_ << endl; // IF I dynamic cast here, can i still change the value?
						dynamic_cast<Toxin*>(*i)->gene_age_ += 1;
						if (V) cout << "new age is "<< dynamic_cast<Toxin*>(*i)->gene_age_<<endl;
					}
					else if(vibrios[row][col].G->IsAntitoxin(*i))
	 				{
						if (V) cout << "found antitoxin with age "<< dynamic_cast<Antitoxin*>(*i)->gene_age_ << endl;
						dynamic_cast<Antitoxin*>(*i)->gene_age_ += 1;
						if (V) cout << "new age is "<<  dynamic_cast<Antitoxin*>(*i)->gene_age_<<endl;
					}
	  				else if(vibrios[row][col].G->IsCore(*i))
					{
						if (V) cout << "found coregene with age "<< dynamic_cast<Core*>(*i)->gene_age_ << endl;
						dynamic_cast<Core*>(*i)->gene_age_ += 1;
						if (V) cout << "new age is "<< dynamic_cast<Core*>(*i)->gene_age_<<endl;
					}
				i++;
				}
			}
		}
	}
}



void ColourClones(TYPE2 **vibrios)
{
	bool V=FALSE;

	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			if(vibrios[row][col].val > 1)				// Val 1 are the cells that died last update
			{
				vibrios[row][col].val = GetSpeciesColour(vibrios,row,col);
			}
		}
	}
}

#endif
