#ifndef GenomeHeader
#define GenomeHeader

#include <typeinfo>
#include <math.h>
#include "Header.hh"
#include "Pearl.hh"
#include "HK.hh"
#include "NonEss.hh"

#include "Transposase.hh"
#include "Noncoding.hh"

class Genome
{
	public:
		list<Pearl*>* StringOfPearls;
		Genome *parent;
		typedef std::list<Pearl*>::iterator iter;
		typedef std::list<Pearl*>::iterator rev_iter;

		int generation_;
		int time_infected_;
		int genomesize_;
		string rowcoltime;
		string GenomeAtBirth;		

		int HKgenes_;					// Number of unique HKgenes needed
		double compstrength;   // competitive strength for local reproduction
		double extra_death;
		bool viable;
		float gene_cost_;						// Every genome knows how much genes cost, avoids all the extra arguments
		float transp_cost_;
		float genome_size_cost_;
		int transformant;
		float fitness_effect_noness;

		list<int> HKgenes;			        // A list of all HK genes
		list<int> transposases;			        // A list of all HK genes
		list<int> nonessential;			        // A list of all HK genes
		

		Genome();
		~Genome();
		void GenerateGenome(int init_nr_HKgenes, int nr_noness, int init_nr_noncoding, int init_nr_tra, float gene_cost, float transp_cost, float genome_size_cost, float init_mob, float fitness_eff);

		// Reproduce / copy parts of genomes
		void CloneGenome(Genome *parent);
		void RecombineGenomes(Genome *parent,Genome *parent2);
		void CopyPartOfGenome(iter begin,iter end);
		void Recombine(Genome *parent,Genome *parent2);
		void CopyPartOfGenomeToTempList(iter begin, int size, list<Pearl*> &PearlListTemp);
		bool FetchTransposon(list<Pearl*> &Fragment, list<Pearl*> &PearlListTemp, int begin_scan, int end_scan, bool cut);

		// mutations
		bool MutateGenome(float gene_mob, float loss, float dupl, float tdupl, float tdel, float inv, float gtn, float ntg, float gendisc);
		bool TransposonDynamics(list<Pearl*> *pearllist,float rate, float break_chance);
		void MutateMobility(iter ii);
		Genome::iter GeneLoss(iter ii);
		Genome::iter GeneDupl(iter ii);
		Genome::iter TandemDupl(iter ii);
		Genome::iter TandemDel(iter ii);
		Genome::iter Invert(iter ii);
		Genome::iter Gene_To_NonCoding(iter ii);
		Genome::iter NonCoding_To_Gene(iter ii);


		bool IntegrateDNA(list<Pearl*> &PearlList_ToIntegrateFrom, int,int,int,bool cut, bool paste, bool hgt, float break_chance);

		bool GeneDiscovery(double chance);

		bool Viable();
		void CalculateCompStrength();

		void Create_Gene_Lists(); 			// Creates a vector of resistance genes, and a list of toxins, giving quick look-up of resistance

		// Fucntions below check type of Pearl in this genome
		
		bool IsHK(Pearl *Pearl) const;
		bool IsNoncoding(Pearl *Pearl) const;
		bool IsTransposase(Pearl *Pearl) const;
		bool IsNonEssential(Pearl *Pearl) const;

		int GetPearlType(Pearl *Pearl) const;
		double GetPearlMobility(Pearl *Pearl) const;

		string ListContent(list<Pearl*> *pearllist, bool ignoreHK=FALSE, bool ignore_nc=FALSE, bool includemobility=FALSE);	// When *pearllist is NULL, it will list the entire genome
		string ListContentShort(list<Pearl*> *pearllist, bool ignoreHK=FALSE, bool ignore_nc=FALSE, bool includemobility=FALSE);	// When *pearllist is NULL, it will list the entire genome
		string ListMobility(list<Pearl*> *pearllist);	// When *pearllist is NULL, it will list the entire genome
};

#endif
