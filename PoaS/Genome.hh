#ifndef GenomeHeader
#define GenomeHeader

#include <typeinfo>
#include <math.h>
#include "Header.hh"
#include "Pearl.hh"
#include "Toxin.hh"
#include "ARG.hh"
#include "Antitoxin.hh"
#include "Core.hh"
#include "Transposase.hh"
#include "Noncoding.hh"

class Genome
{
	public:
		list<Pearl*>* StringOfPearls;
		typedef std::list<Pearl*>::iterator iter;
		typedef std::list<Pearl*>::iterator rev_iter;

		int generation_;
		int genomesize_;
		

		int coregenes_;					// Number of unique coregenes needed
		int tot_nr_ARGs;			// Number of toxins in the WOOOORLLLD
		double compstrength;   // competitive strength for local reproduction
		bool viable;
		float gene_cost_;						// Every genome knows how much genes cost, avoids all the extra arguments
		float genome_size_cost_;
		int transformant;

		vector<int> resistance_lookup;			// A vector for resistance-lookup
		list<int> resistance;					// A list of resistance genes
		list<int> toxins;					    // A list for produced toxins
		list<int> coregenes;			        // A list of all core genes
		list<int> transposases;			        // A list of all core genes
		

		Genome();
		~Genome();
		void GenerateGenome(int init_nr_args, int total_nr_args, int init_nr_coregenes, int init_nr_noncoding, float gene_cost, float genome_size_cost, float init_mob);

		// Reproduce / copy parts of genomes
		void CloneGenome(const Genome *parent);
		void CopyPartOfGenome(iter begin,iter end);
		void CopyPartOfGenomeToTempList(iter begin, int size, list<Pearl*> &PearlListTemp);
		bool FetchTransposon(list<Pearl*> &Fragment, list<Pearl*> &PearlListTemp, int begin_scan, int end_scan, bool cut);

		// mutations
		bool MutateGenome(float gene_mob, float loss, float dupl, float tdupl, float tdel, float inv, float gtn, float ntg, float gendisc);
		bool TransposonDynamics();
		void MutateMobility(iter ii);
		Genome::iter GeneLoss(iter ii);
		Genome::iter GeneDupl(iter ii);
		Genome::iter TandemDupl(iter ii);
		Genome::iter TandemDel(iter ii);
		Genome::iter Invert(iter ii);
		Genome::iter Gene_To_NonCoding(iter ii);
		Genome::iter NonCoding_To_Gene(iter ii);

		bool IntegrateDNA(list<Pearl*> &PearlList_ToIntegrateFrom, int,int,int,bool cut, bool paste);

		bool GeneDiscovery(double chance);

		void CalculateCompStrength();

		// Toxin / anti-toxin shizzle
		void Create_Gene_Lists(int nr_toxins); 			// Creates a vector of resistance genes, and a list of toxins, giving quick look-up of resistance
		int GetTotalNrToxins(int type);						// Get total number of toxins of type X in this genome
		int GetTotalNrDefense(int type);
		double GetFractionTEs();

		// Fucntions below check type of Pearl in this genome
		bool IsToxin(Pearl *Pearl) const;
		bool IsAntitoxin(Pearl *Pearl) const;
		bool IsARG(Pearl *Pearl) const;
		bool IsCore(Pearl *Pearl) const;
		bool IsNoncoding(Pearl *Pearl) const;
		bool IsTransposase(Pearl *Pearl) const;

		int GetPearlType(Pearl *Pearl) const;
		double GetPearlMobility(Pearl *Pearl) const;

		string ListContent(list<Pearl*> *pearllist, bool ignorecore=FALSE, bool ignore_nc=FALSE, bool includemobility=FALSE);	// When *pearllist is NULL, it will list the entire genome
		string ListMobility(list<Pearl*> *pearllist);	// When *pearllist is NULL, it will list the entire genome
};

#endif
