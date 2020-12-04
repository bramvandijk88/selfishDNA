#include "clibs.hpp"
#include "functions.hpp"
#include <curses.h>
#include <iomanip>

static TYPE2** Vibrios;
static TYPE2** Environment;
static TYPE2** DNA;
static TYPE2** Genomesize;
static TYPE2** Transformant;

// Simulation stuff
int myseed;
int myscale =2;

// BIRTH/DEATH AND FITNESS STUFF (see par.txt)
float death, birth, ab_penalty, birthNON, gene_cost, genome_size_cost;
int max_internal_resource;																		// LS: What's this, for cases of growth on resource uptake ability?

// Mutations (see par.txt)
float gene_loss, gene_dupl, tandem_dupl, tandem_del, inversions, gene_to_noncoding, noncoding_to_gene, mut_rate_scaling, gene_mob, gene_discovery;

// Toxins and antitoxins stuff
int total_nr_args, init_nr_noncoding, init_nr_args, init_nr_coregenes;
float init_mob;

//eDNA
int diff;
float degr;

double extradeath, genediffusion;

std::vector<double> ABs_concentration;
std::vector<double> ABs_intake;
std::vector<double> ABs_timer;




int dodisplay, opendisplay, summaryinterval, displayinterval, fieldsize, livegraphs;

// Saving/plotting vars below
bool makemovie = 1; int moviecounter = 0; bool mixing = 0; bool diffusing = 0; bool mutation = 1; int genedisplay = 1;
bool Graphs = 1; bool Output = 1; bool mobilitygathering = 1; bool getgenefrequencies = 1;

string dumpfolder = "/dumps";
string genedata = "/genedata";
string specfolder = "/speciesfreq";
string heatfolder = "/heatmaps";
string timefolder = "/timeframes";
string ancefolder = "/ancestor";
string map_extension;
string dirbuf;
string folder;
string toxin_mobility_folder;
string antitoxin_mobility_folder;
string coregene_mobility_folder;
string frequencies_folder;

bool mix = false; 
bool mixgrid = false;
bool dohgt = true;
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
		if(readOut == "-Size") {fieldsize = atoi(argv_g[i+1]); cmd_fieldsize = true;}
		if(readOut == "-Scale") {myscale = atoi(argv_g[i+1]);}

		if(readOut == "-Diff") {genediffusion = atoi(argv_g[i+1]); cmd_eDNAdiffusion = true;}
		if(readOut == "-Mix") {mix = true; }
		if(readOut == "-MixGrid") {mixgrid = true; }
		if(readOut == "-noHGT") {dohgt = false; }

		readOut = (char*)argv_g[i];
		string map_extension ="" ;
		if(readOut == "-o")
		{
			outputgiven = TRUE;
			map_extension = argv_g[i+1];
			folder = "/home/brem/ARG_output/" + currentDateTime() + "_" + map_extension;	// Paths to local storages
			string command = "mkdir -p " + folder;
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

	string datafolder = folder+"/data";														//LS: makes data folder
	string command_data = "mkdir -p " + datafolder;
	system(command_data.c_str());

	string model_folder = datafolder+"/Current_Model";					// LS: Current model and parameters are saved
	string command_model = "mkdir -p " + model_folder;
	string command_copy = "cp -r ./* " + model_folder;
	system(command_model.c_str());
	system(command_copy.c_str());

	// HOW TO OPEN AND WRITE IN A FILE
  // ofstream modelfile;
  // modelfile.open((model_folder+"/Model.dat").c_str());
  // cout << " writing into model file" << endl;
  // modelfile << " Here comes the current version of the model of this run" << endl;

	if (mobilitygathering){
	 	toxin_mobility_folder = datafolder+"/toxin_mobilities";					// LS: In here go files with mobility of every toxin (made the first datainterval timestep and appended the others)
		string command_toxmob = "mkdir -p " + toxin_mobility_folder;
		system(command_toxmob.c_str());

		antitoxin_mobility_folder = datafolder+"/antitoxin_mobilities"; // LS: In here go files with mobility of every antitoxin
		string command_antitoxmob = "mkdir -p " + antitoxin_mobility_folder;
		system(command_antitoxmob.c_str());

		coregene_mobility_folder = datafolder+"/coregene_mobilities"; // LS: In here go files with mobility of every coregene
		string command_coremob = "mkdir -p " + coregene_mobility_folder;
		system(command_coremob.c_str());
  }
  	if (getgenefrequencies){
	 	frequencies_folder = datafolder+"/gene_frequencies";					// LS: In here go files with frequency of every gene (made the first datainterval timestep and appended the others)
		string command_toxmob = "mkdir -p " + frequencies_folder;
		system(command_toxmob.c_str());
  	}

	ReadOptions("par.txt");								// Read file
 	InDat("%f", "death", (int*)&death);						// Get death rate from file
	InDat("%f", "birth", (int*)&birth);
	InDat("%f", "birthNON", (int*)&birthNON);
	InDat("%f", "ab_penalty", (int*)&ab_penalty);
	InDat("%d", "seed", (int*)&myseed);
	InDat("%d", "max_internal_resource", (int*)&max_internal_resource);
	InDat("%d", "total_nr_args", (int*)&total_nr_args);
	InDat("%d", "init_nr_args", (int*)&init_nr_args);
	InDat("%d", "init_nr_noncoding", (int*)&init_nr_noncoding);
	InDat("%d", "init_nr_coregenes", (int*)&init_nr_coregenes);
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

	// Code below takes seed from the parameter file, or chooses time as a seed if this parameter file contains a 0 as seed (this seed is written to a file named seed.txt)
	if(myseed==0)
    	ulseedG = time(0);
	else ulseedG = myseed;

	MaxTime = 1000000;
	nrow = fieldsize;
	ncol = fieldsize;
	nplane = 5;
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

  MakePlane(&Vibrios,&DNA,&Environment,&Genomesize,&Transformant);

  InitialSet(Vibrios,1,0,2,initdens);								// Initial density of bacteria
  for(int i=0;i<3;i++)
  	RefreshEnvironment(Environment);							// Puts 3 circles of resource on the Resource plane
  InitialiseGenomes(Vibrios, Environment, init_nr_args, total_nr_args, init_nr_coregenes, init_nr_noncoding);							// Initialises all pearl-on-a-string genomes for LIVING cells
  InitialiseDNAPlane(DNA);																																									// Every pixel gets an empty eDNA list
  cout << endl << endl << " -= Start of simulation =- " << endl << endl;


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
  
  ABs_concentration.assign(total_nr_args, 0.0);
  ABs_intake.assign(total_nr_args, 0);
   ABs_timer.assign(total_nr_args, 0);
}

void NextState(int row,int col)
{
	if(Vibrios[row][col].val>1){																// Strictly speaking > 91.
        if(genrand_real1()< death + Vibrios[row][col].ab_effect)  		 // LS: so note, deathrate value is only toxin induced deathrate
		{
			Vibrios[row][col].val = 1;											// 1 = Temporary dead cell, which is cleaned up and set to 0 later
			if(dohgt) SpillDNA(Vibrios,DNA,row,col);
        }	
	}
	else if(Vibrios[row][col].val==0){ // && genrand_real1() < birth )			// LS: with birth at 1 this always happens
        	CompeteAndReproduce(Vibrios,Environment,row,col);
	}
	

}

void DiffuseDNA(int row, int col)
{
	std::vector< std::list<Pearl*> >::iterator fr_iter;
	int pos=0;
	bool V = FALSE;
		for (fr_iter = DNA[row][col].DNA->Fragments->begin(); fr_iter != DNA[row][col].DNA->Fragments->end(); fr_iter++)
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


void Update(void)
{
  //if(genrand_real1() < 0.05)													// LS: move chance to par.txt?
	RefreshEnvironment(Environment);									// Add new resource particle with a certain probability

	//for(int i = 0; i<3;i++) ObstacleMargolus(Vibrios,Environment,1);	// Vibrios do not diffuse if on a resource particle
	

	CalculateABEffect(Vibrios);								// extra death due to toxin is calculated. Competitive strenght is calcutated if Genome changes. No explicit fitness
	Synchronous(1,Vibrios);																//LS: note to self: holds next-state function/loop thingy
	CleanUpTheDead(Vibrios, DNA);
	
	for(int i =0;i<genediffusion;i++)
		Asynch_func(*DiffuseDNA);    // Special diffusion stuff by Brem to diffuse (potentially) unique objects. W.i.p. Also see DiffuseParticles() below
	if(dohgt) 
	{
		for(int row=1; row<=nrow; row++)for(int col=1; col<=ncol; col++) DegradateDNA(DNA,row,col,degr); 
		DoHGT(Vibrios,DNA,total_nr_args);
	}
		
		
	AgeGenes(Vibrios);
	//if (Time > 600000) PerfectMix(Vibrios);
	if(mixgrid) PerfectMix(Vibrios);
	if(mix)
	{
	        PerfectMix(DNA);
	}


	
	
	
	if(Time%displayinterval==0)
	{
		for(int row=1; row<=nrow; row++)
		for(int col=1; col<=ncol; col++)
		{
			if(Vibrios[row][col].val>0)
			{
				Genomesize[row][col].val = 11+min((float)30.0*(Vibrios[row][col].G->transposases.size()>0),(float)39.0);		
				Transformant[row][col].val = Vibrios[row][col].G->transformant;		
			}
			else
			{
				Genomesize[row][col].val = 0;
				Transformant[row][col].val = 0;		
			}
			
			float size_sum = DNA[row][col].DNA->Fragments->size();
			
			DNA[row][col].val = 51+min((int)((float)size_sum*2),39);
		}
		cout << GetToxin_Antitoxin_Stats(Vibrios);
		if(dodisplay)
		{
			ColourClones(Vibrios);
			Display(nplanedisp,Vibrios,DNA,Genomesize,Transformant);
		}
	}
	if(Time%summaryinterval==0)
	{
		cout << "\n\n\nSummary at simulation time: " << Time << endl;
		cout << Get_Genome_Stats(Vibrios) << endl;
		
		PrintMostAbundant(Vibrios);
		cout << " " << endl;
		//cout << PrintAverageToxinMobilities(Vibrios);
		cout << PrintAverageARGMobilities(Vibrios);
		
		cout << endl << endl << endl;
		cout << PrintAverageCoreMobilities(Vibrios);
		cout << PrintAverageNCMobilities(Vibrios) << endl;
		cout << endl << endl << endl;

		dirbuf = folder + timefolder;
		DrawSlide2((char*)dirbuf.c_str(),Vibrios, DNA, Genomesize,Transformant);
		if(Time%100==0)
		{
			cout << ">\t";
			for (int i=0;i<total_nr_args;i++)
				cout << "f(ARG" << i << ")\t" << "[ARG" << i << "]\t";
			cout << "AVGfit\tPopsize\tFract_IRs\tFract_Transposase"<< endl;
		}
	}

//DATA STUFF

	int datainterval= 1000;
	

	if(Time%50==0)
	{
		for(int i=0;i<total_nr_args;i++) 
			{			
				if(ABs_intake[i] == 0 && genrand_real1() < 0.02 && Time > 100000)
				{
					ABs_intake[i] = 1;
					ABs_timer[i] = 1500+genrand_real1()*1000;
				}
				
			}
	}
	
	for(int i=0;i<total_nr_args;i++) 
	{
		
		if(ABs_concentration[i]>0)
		{
			ABs_concentration[i] *= 0.975;
		}
		if(ABs_intake[i]==1){ 
			ABs_concentration[i] += 0.025;
			ABs_timer[i]--;
		}
		if(ABs_timer[i] <= 0) {
				ABs_intake[i] = 0;	
				ABs_timer[i] = 0;
			}	
	}
		
		#include "gotmouse.cpp"


}
