// Import C libraries (Cash, math, iomanip) and custom functions (compete cells, degrade DNA, etc.)
#include "clibs.hpp"
#include "functions.hpp"
#include <curses.h>
#include <iomanip>

static TYPE2** Vibrios;	// The TYPE2** (cash grid) on which bacterial cells live. Named after vibrionacaea because they're cool. 
static TYPE2** DNA;		// The TYPE2** (cash grid) on which eDNA resides
static TYPE2** Genomesize;	// Cash grid (visualisation only, -x)
static TYPE2** Transformant; // Cash grid (visualisation only, -x)

// Simulation stuff
int myseed;
int myscale =2;

// BIRTH/DEATH AND FITNESS STUFF (see par.txt)
float death, birth, ab_penalty, birthNON, gene_cost, genome_size_cost;
						
// Mutations (see par.txt)
float gene_loss, gene_dupl, tandem_dupl, tandem_del, inversions, gene_to_noncoding, noncoding_to_gene, mut_rate_scaling, gene_mob, gene_discovery, fitness_effect_noness;

// Toxins and antitoxins stuff
int init_nr_noncoding, init_nr_args, init_nr_HKgenes, init_nr_noness;

float init_mob;

//eDNA
int diff;
float degr;

double genediffusion;


int dodisplay, opendisplay, summaryinterval, displayinterval, fieldsize, livegraphs;

// Saving/plotting vars below
bool makemovie = 1; int moviecounter = 0; bool mixing = 0; bool diffusing = 0; bool mutation = 1; int genedisplay = 1;
bool Graphs = 0; bool Output = 1; bool mobilitygathering = 1; bool getgenefrequencies = 1; bool extinct = false;

string ancfolder; 

string map_extension;
string dirbuf;
string folder;
string toxin_mobility_folder;
string antitoxin_mobility_folder;
string coregene_mobility_folder;
string frequencies_folder;

// Command line controllable options
bool mix = false; 
bool mixgrid = false;
bool dohgt = true;
bool dojump = true;
bool sexual = false;
float jumprate = 0.1;			//j
float uptakerate = 0.1;			//u
float transp_cost = 0.0;		//c
float break_chance = 0.0;		//b

int start_sge_influx = 0;		
double sge_influx = 0.0;
int stop_sge_influx = 0;

set<Genome*> Ancestors;

int dna_diff = 2;				

void Initial(void)
{

	/*
	Below is the parsing of options from par.txt, and also some alternative methods to define these options from the command line.
	*/
	string readOut;
	bool cmd_eDNAdiffusion = false;	// Keep track if user defined eDNA-diffusion on command-line
	bool outputgiven = FALSE;
	bool cmd_fieldsize = false;		// Keeps track if user defined field-size on command-line
	for(int i = 0; i < (int)argc_g; i++)
	{
		readOut = (char*)argv_g[i];
		if(readOut == "-Size") {fieldsize = atoi(argv_g[i+1]); cmd_fieldsize = true;}			//W,H
		if(readOut == "-Scale") {myscale = atoi(argv_g[i+1]);}						
		if(readOut == "-Diff") {genediffusion = atof(argv_g[i+1]); cmd_eDNAdiffusion = true;}	//D
		if(readOut == "-Jumprate") {jumprate = atof(argv_g[i+1]);}								//j	
		if(readOut == "-Sexual") {sexual = true;}		
		if(readOut == "-Cost") {transp_cost = atof(argv_g[i+1]);}								//c
		if(readOut == "-BreakChance") {break_chance = atof(argv_g[i+1]);}						//b
		if(readOut == "-Uptakerate") {uptakerate = atof(argv_g[i+1]);}							//u
		if(readOut == "-startInTra") {start_sge_influx = atoi(argv_g[i+1]);}				
		if(readOut == "-stopInTra") {stop_sge_influx = atoi(argv_g[i+1]);}
		if(readOut == "-rateInTra") {sge_influx = atof(argv_g[i+1]);}
		if(readOut == "-MixDNA") {mix = true; }	
		if(readOut == "-MixPop") {mixgrid = true; }
		if(readOut == "-noHGT") {dohgt = false; }
		if(readOut == "-noJump") {dojump = false; }

		readOut = (char*)argv_g[i];
		string map_extension ="" ;
		if(readOut == "-o")
		{
			outputgiven = TRUE;
			map_extension = argv_g[i+1];
			//folder = \"/mnt/d/Selfish_output/\" + currentDateTime() + \"_\" + map_extension;	// Paths to local storages
			folder = "/mnt/d/Selfish_output/" + map_extension;	// Paths to local storages
			ancfolder = folder+"/Ancestor_traces";
			string command = "rm -rf " + folder;
			system(command.c_str());
			command = "mkdir -p " + folder;
			system(command.c_str());
			command = "mkdir -p " + ancfolder;
			system(command.c_str());
		}

  }
  
	for(int i = 0; i < (int)argc_g; i++)
	{
		if(readOut == "-t")
		{
			readOut = (char*)argv_g[i];
			string out = folder+"/log.out";
			string outerr = folder+"/log.error";
			cout << "Redirecting log to:";
			cout << out << endl;
			freopen(out.c_str(),"a",stdout);	// Als niet werkt weer 'w' doen
			//int serr = dup(fileno(stderr));
			freopen(outerr.c_str(),"a",stderr);
		}
	}
	if(outputgiven==FALSE)
	{
		cout << "No output given, saved under NoName" << endl;
		folder = "/home/brem/ARG_output/" + currentDateTime() + "_" + "NoName";
	}
	


	string model_folder = folder+"/current_model_code";					// LS: Current model and parameters are saved
	string command_model = "mkdir -p " + model_folder;
	string command_copy = "cp -r ./* " + model_folder;
	system(command_model.c_str());
	system(command_copy.c_str());

	// HOW TO OPEN AND WRITE IN A FILE
  // ofstream modelfile;
  // modelfile.open((model_folder+"/Model.dat").c_str());
  // cout << " writing into model file" << endl;
  // modelfile << " Here comes the current version of the model of this run" << endl;

	
	ReadOptions("par.txt");								// Read file
 	InDat("%f", "death", (int*)&death);						// Get death rate from file
	InDat("%f", "birth", (int*)&birth);
	InDat("%f", "birthNON", (int*)&birthNON);
	InDat("%d", "seed", (int*)&myseed);
	InDat("%d", "init_nr_noncoding", (int*)&init_nr_noncoding);
	InDat("%d", "init_nr_HKgenes", (int*)&init_nr_HKgenes);
	InDat("%d", "init_nr_noness", (int*)&init_nr_noness);
	InDat("%f", "fitness_effect_noness", (int*)&fitness_effect_noness);
	
	InDat("%f", "init_mob", (int*)&init_mob);
	if(!cmd_fieldsize) InDat("%d", "fieldsize", (int*)&fieldsize);			// Get field size from file IF not already defined on command-line

	InDat("%f", "gene_cost", (int*)&gene_cost);
	InDat("%f", "genome_size_cost", (int*)&genome_size_cost);

	InDat("%f", "mut_rate_scaling", (int*)&mut_rate_scaling);
	InDat("%f", "gene_loss", (int*)&gene_loss);
	InDat("%f", "gene_mob", (int*)&gene_mob);
	InDat("%f", "gene_dupl", (int*)&gene_dupl);
	InDat("%f", "tandem_dupl", (int*)&tandem_dupl);
	InDat("%f", "tandem_del", (int*)&tandem_del);
	InDat("%f", "inversions", (int*)&inversions);
	InDat("%f", "gene_to_noncoding", (int*)&gene_to_noncoding);
	InDat("%f", "noncoding_to_gene", (int*)&noncoding_to_gene);
	InDat("%f", "gene_discovery", (int*)&gene_discovery);

	InDat("%f", "degr", (int*)&degr);
  	InDat("%d", "diff", (int*)&diff);

	InDat("%f", "dodisplay", (int*)&dodisplay);
	InDat("%f", "opendisplay", (int*)&opendisplay);
	InDat("%d", "displayinterval", (int*)&displayinterval);
	InDat("%d", "summaryinterval", (int*)&summaryinterval);

	
	ofstream file;
  	stringstream filepath;
  	filepath  << folder << "/command.txt";
  	string filename = filepath.str();
  	file.open (filename.c_str(), ios::app);
	MaxTime = 200000;
	for(int i = 0; i < (int)argc_g; i++)
	{
		readOut = (char*)argv_g[i];
		file << readOut << " ";
		if(readOut == "-Seed") {myseed = atoi(argv_g[i+1]);}										//s
		if(readOut == "-Disco") {gene_discovery = atof(argv_g[i+1]);}								
		if(readOut == "-Phi") {init_mob = atof(argv_g[i+1]);}										//init value of phi
		if(readOut == "-Noncoding") {init_nr_noncoding = atoi(argv_g[i+1]);}						
		if(readOut == "-HK") {init_nr_HKgenes = atoi(argv_g[i+1]);}	
		if(readOut == "-MaxTime") {MaxTime = atoi(argv_g[i+1]);}	
		if(readOut == "-Non") {birthNON = atof(argv_g[i+1]);}			
		if(readOut == "-Del") { gene_loss = atof(argv_g[i+1]); }									//mutations
                if(readOut == "-Dup") { gene_dupl = atof(argv_g[i+1]); }
                if(readOut == "-TDel") {tandem_del = atof(argv_g[i+1]); }
                if(readOut == "-TDup") {tandem_dupl = atof(argv_g[i+1]); }
                if(readOut == "-Noness") {init_nr_noness = atoi(argv_g[i+1]); }
                if(readOut == "-Noness_eff") {fitness_effect_noness = atof(argv_g[i+1]); }

	}
	file.close();
	// Code below takes seed from the parameter file, or chooses time as a seed if this parameter file contains a 0 as seed (this seed is written to a file named seed.txt)
	if(myseed==0)
    	ulseedG = time(0);
	else ulseedG = myseed;

	
	nrow = fieldsize;
	ncol = fieldsize;
	nplane = 4;
    nplanedisp = 4;
	scale = myscale;
	boundary = WRAP;
	//boundaryvalue2 = (TYPE2){0,0,0,0,0,0.,0.,0.,0.,0.,{'S','S','S','S','S','S','S','S'}};
}

void InitialPlane(void)
{
  float initdens;
  ReadOptions("par.txt");
  InDat("%f", "initial_density", (int*)&initdens);

  MakePlane(&Vibrios,&DNA,&Genomesize,&Transformant);

  InitialSet(Vibrios,1,0,2,initdens);								// Initial density of bacteria

  InitialiseGenomes(Vibrios, init_nr_HKgenes, init_nr_noness,init_nr_noncoding);							// Initialises all pearl-on-a-string genomes for LIVING cells
  InitialiseDNAPlane(DNA);																																									// Every pixel gets an empty eDNA list
  //cout << endl << endl << " -= Start of simulation =- " << endl << endl;


  ColorRGB(0,0,0,0);                    // Color 0 reserved for BLACK
  ColorRGB(1,50,50,50);                 // Color 1 reserved for ALMOST BLACK == dieing cell that must still be cleaned up
  ColorRGB(2,0,0,255);                  // Color 2 reserved for BLUE
  ColorRGB(3,255,0,0);                  // Color 3 reserved for RED
  ColorRGB(4,0,255,0);									// Color 4 reserved for GREEN
  // Color 5-10 still reserved for some other stuff :)
  ColorRGB(255,255,255,255);            // Color 255 reserved for WHITE
  for(int q = 11; q<= 50; q++) ColorRGB(q,((float)(q-11)/39.0)*(254),((float)(q-11)/39.0)*(125),0);	// Oranjeish gradient for resource
  for(int q = 51; q<= 90; q++) ColorRGB(q,0,((float)(q-51)/39.0)*(254),0);			// Black to Green colour gradient
  for(int q = 91; q<= 254; q++) ColorRGB(q,(int)ceil(genrand_real1()*200+55),(int)ceil(genrand_real1()*200+55),(int)ceil(genrand_real1()*200+55));    // Color 91 - 254 are "random" colours
  
}

void NextState(int row,int col)
{
	if(Vibrios[row][col].val>1)
	{	
		// Random cell death with rate <death>								
		// LS: Cells with too many transposases die, but note that this is not the mechanisms responsible for streamlining. 
		//if(genrand_real1()< death || Vibrios[row][col].G->transposases.size() > 1000)  		 
        //if(genrand_real1()< death || !Vibrios[row][col].G->Viable() ||  Vibrios[row][col].G->transposases.size() > 1000) 
		if(Vibrios[row][col].G->transposases.size()>0) Vibrios[row][col].G->time_infected_++;
		if(genrand_real1() < death || !Vibrios[row][col].G->Viable() )  		 		
		{
			Vibrios[row][col].val = 1;											// 1 = Temporary dead cell, which is cleaned up and set to 0 later
			if(dohgt) SpillDNA(Vibrios,DNA,row,col);			
        }	
		else
		{
			bool integrated = false;
			if(dojump) integrated = Vibrios[row][col].G->TransposonDynamics(NULL,jumprate,break_chance);			
			// bool integrated=false;
			if(integrated) 
			{
				Vibrios[row][col].G->transformant = 2;	
				Vibrios[row][col].G->Create_Gene_Lists();
				Vibrios[row][col].G->CalculateCompStrength();
			}
		}
	}
	else if(Vibrios[row][col].val==0){ // && genrand_real1() < birth )			// LS: with birth at 1 this always happens
        	CompeteAndReproduce(Vibrios,row,col);
	}	
}

void DiffuseDNA(int row, int col)
{
	std::vector< std::list<Pearl*> >::iterator fr_iter;
	int pos=0;
	bool V = FALSE;
		for (fr_iter = DNA[row][col].DNA->Fragments->begin(); fr_iter != DNA[row][col].DNA->Fragments->end(); fr_iter++)
		{			
			if(genrand_real1() < genediffusion) 
			{
				pos ++;
				if(V) cout << "Swapping eDNAfragment at pos " << pos << endl;
				int randneigh = genrand_int(1,8);
				TYPE2* neigh = GetNeighborP(DNA,row,col,randneigh);
				neigh->DNA->Fragments->push_back(*fr_iter);
				fr_iter = DNA[row][col].DNA->Fragments->erase(fr_iter);
				fr_iter--;
			}
		}
}


void Update(void)
{

	Synchronous(1,Vibrios);																//LS: note to self: holds next-state function/loop thingy
	if(Time%1000==0) AncestorTrace(Vibrios);
	
	CleanUpTheDead(Vibrios, DNA);
	
	
	Asynch_func(*DiffuseDNA);    		// Special diffusion stuff by Brem to diffuse (potentially) unique objects. W.i.p. Also see DiffuseParticles() below
	
	if(dohgt) 
	{
		for(int row=1; row<=nrow; row++)for(int col=1; col<=ncol; col++) DegradateDNA(DNA,row,col,degr); 
		if(Time > start_sge_influx && Time < stop_sge_influx) for(int row=1; row<=nrow; row++)for(int col=1; col<=ncol; col++) InfluxDNA(DNA,row,col,sge_influx); 
		DoHGT(Vibrios,DNA);
	}
		
		
	AgeGenes(Vibrios);
	if(mixgrid) PerfectMix(Vibrios);
	if(mix) PerfectMix(DNA);
	// if(Time>100000) dohgt=false;
	if(Time%displayinterval==0)
	{
		for(int row=1; row<=nrow; row++)
		for(int col=1; col<=ncol; col++)
		{
			if(Vibrios[row][col].val>0)
			{
				Genomesize[row][col].val = 11+min((float)30.0*(Vibrios[row][col].G->transposases.size()>0),(float)39.0);	
				Transformant[row][col].val = 51+min((float)0.5*(Vibrios[row][col].G->StringOfPearls->size()),(float)39.0);		
				//Transformant[row][col].val = Vibrios[row][col].G->transformant;		
			}
			else
			{
				Genomesize[row][col].val = 0;
				Transformant[row][col].val = 0;		
			}
			
			float size_sum = DNA[row][col].DNA->Fragments->size();
			
			DNA[row][col].val = 51+min((int)((float)size_sum*2),39);
		}
		//cout << Get_Sim_Stats(Vibrios,false);


		
		if(dodisplay)
		{
			ColourClones(Vibrios);
			Display(nplanedisp,Vibrios,DNA,Genomesize,Transformant);
		}
	}
	int datainterval= 50;
	if(Time%datainterval==0)
	{
		Get_Sim_Stats(Vibrios,true);		
		// PrintMostAbundant(Vibrios);
		// if(Time%5000<1000 && Time%5000>0) Dump_Grids(Vibrios,DNA);
	}	
	
	
	datainterval= 25;
	if(Time%datainterval==0)
	{
		//PrintMostAbundant(Vibrios);		
		Dump_Genomes(Vibrios);
		Dump_Grids(Vibrios,DNA);		
	}
	
	//if(Time == MaxTime || extinct)
	if(Time == MaxTime)
	{
		exit(0);
	} 
	#include "gotmouse.cpp" // For pausing the display / asking stats in interactive mode (-x)
}
