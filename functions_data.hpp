extern int nrow;
extern int ncol;
extern int Time;
extern bool Graphs;

extern int total_nr_args;
extern int init_nr_coregenes;
extern string toxin_mobility_folder;
extern string folder;
extern string antitoxin_mobility_folder;
extern string coregene_mobility_folder;
extern string frequencies_folder;

extern vector<double> ABs_concentration;
/********* Declarations  *********/

string Get_Sim_Stats(TYPE2**, bool save);   // LS: why is the word string not purple?

string filename; ofstream myfile;

string Get_Sim_Stats(TYPE2 **vibrios,bool save)
{
    //cout << "Getting toxstats for " << total_nr_toxins << " toxins" << endl;
    vector<int> toxcount(total_nr_args,0);
    vector<int> argcount(total_nr_args,0);
    int pop_size = 0;
    double num_TEs = 0;
    double sumfit = 0.0;
    double sumtra = 0.0;
    double sumtrans = 0.0;
    double sumtrans2 = 0.0;
    double sum_genomesize = 0;
    double sum_hgt_hk = 0;
    double sum_hgt_arg = 0;
    double sum_hgt_tra = 0; 
    double sum_jump_hk = 0;
    double sum_jump_arg = 0;
    double sum_jump_tra = 0; 
    int num_nc_pearls=0;
    int num_IR_pearls=0;
    int num_hk_pearls=0;
    int num_arg_pearls=0;
    Noncoding* nc;

    double sum_nc_mobility = 0.0;   
    double sum_IR_mobility = 0.0;
    
    for(int row = 1; row<=nrow;row++){
        for(int col=1;col<=ncol;col++){
            if(vibrios[row][col].val >= 1)
            {
                pop_size++;
                num_TEs += vibrios[row][col].G->GetFractionTEs();
                sumtra += vibrios[row][col].G->transposases.size();
                sumtrans += vibrios[row][col].G->transformant==1;
                sumtrans2 += vibrios[row][col].G->transformant==2;
		            sumfit += vibrios[row][col].G->compstrength;
                sum_genomesize += vibrios[row][col].G->StringOfPearls->size();
                for(int i=0;i<total_nr_args;i++)
                {
                    toxcount[i] += vibrios[row][col].G->GetTotalNrToxins(i);                    
                    // argcount[i] += vibrios[row][col].G->resistance_lookup[i];
                    if(vibrios[row][col].G->resistance_lookup[i]>0) argcount[i]++;
                }

                Genome::iter i = vibrios[row][col].G->StringOfPearls->begin();
                while(i!=vibrios[row][col].G->StringOfPearls->end())
                {
                  if(vibrios[row][col].G->IsNoncoding(*i))
                  {                  
                    nc=dynamic_cast<Noncoding*> (*i);
                    sum_nc_mobility += nc->mobility;
                    num_nc_pearls ++;
                    if(nc->mobility > 0.75)
                    {
                      num_IR_pearls++;
                      sum_IR_mobility+=nc->mobility;
                    }
                  }
                  else if(vibrios[row][col].G->IsARG(*i))
                  {
                    sum_hgt_arg+= (*i)->num_horizontal_transfers_;
                    sum_jump_arg+= (*i)->num_jumps_;
                    num_arg_pearls++;
                  }
                  else if(vibrios[row][col].G->IsCore(*i))
                  {
                    sum_hgt_hk+= (*i)->num_horizontal_transfers_;
                    sum_jump_hk+= (*i)->num_jumps_;
                    num_hk_pearls++;
                  }
                  else if(vibrios[row][col].G->IsTransposase(*i))
                  {
                    sum_hgt_tra+= (*i)->num_horizontal_transfers_;
                    sum_jump_tra+= (*i)->num_jumps_;
                  }
                  if(save)
                  {
                    (*i)->num_horizontal_transfers_=0; // Reset counter for num transfers
                    (*i)->num_jumps_ = 0;
                  }
                  i++;
                }
            }
        }
    }

    stringstream genestats;
    if(pop_size == 0) {cout << "Shutting down simulation as no more cells are alive" << endl; exit(0); }
    
    float fract_TEs = (float)num_TEs/pop_size;
    float fract_arg = (float)num_arg_pearls/pop_size;
    float fract_tra = (float)sumtra/pop_size;
    float fract_nc = (float)num_nc_pearls/pop_size;
    float fract_trans = (float)sumtrans/pop_size;
    float fract_trans2 = (float)sumtrans2/pop_size;

    float fract_hgt_hk = num_hk_pearls > 0 ? (float)sum_hgt_hk/num_hk_pearls: 0.0;
    float fract_hgt_tra = sumtra > 0 ? (float)sum_hgt_tra/sumtra : 0.0;
    float fract_hgt_arg = num_arg_pearls > 0 ? (float)sum_hgt_arg/num_arg_pearls : 0.0;

    float fract_jump_hk = num_hk_pearls > 0 ? (float)sum_jump_hk/num_hk_pearls: 0.0;
    float fract_jump_tra = sumtra > 0 ? (float)sum_jump_tra/sumtra : 0.0;
    float fract_jump_arg = num_arg_pearls > 0 ? (float)sum_jump_arg/num_arg_pearls : 0.0;
  
    if(Time==0)
    {
      //genestats << ">";
      genestats << "Time\tAVG_g\tAVG_f\tPopsize\tFracHK\tFracIRs\tFracTra\tFracARG\tFractNC\tHGT_hk\tJump_hk\tHGT_arg\tJump_arg\tHGT_tra\tJump_tra\tPhi_NC\tPhi_IR\t";
      for (int i=0;i<total_nr_args;i++)
        genestats << "f(ARG" << i << ")\t" << "[ARG" << i << "]\t";
      genestats << "\n";
		}
    
    genestats << fixed << setprecision(4) << Time << "\t" << (float)sum_genomesize/pop_size << "\t";
    
    genestats << sumfit/pop_size << "\t" << pop_size << "\t" << (float)num_hk_pearls/pop_size << "\t" << fract_TEs << "\t" << fract_tra << "\t" << fract_arg << "\t" << fract_nc << "\t";
    genestats << fract_hgt_hk << "\t" << fract_jump_hk << "\t" << fract_hgt_arg << "\t" << fract_jump_arg << "\t" << fract_hgt_tra << "\t" << fract_jump_tra << "\t" << (float)sum_nc_mobility/num_nc_pearls << "\t" << (float)sum_IR_mobility/num_IR_pearls << "\t";
    for (int i=0;i<total_nr_args;i++)
    {
		  genestats << fixed << setprecision(4) << (float)argcount[i]/pop_size << "\t";    				
      genestats << fixed << setprecision(4) << ABs_concentration[i] << "\t";
    }
    
    genestats << endl;
    string returnstring;
    
    returnstring+=genestats.str();
    if(Graphs)
    {
      double values_1[16] = {(float)pop_size,nrow*ncol+100,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double values_2[16] = {(float)argcount[0]/pop_size,(float)argcount[1]/pop_size,(float)argcount[2]/pop_size,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double values_3[16] = {ABs_concentration[0],ABs_concentration[1],ABs_concentration[2],1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double values_4[16] = {fract_trans,fract_tra,fract_trans2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      for(int i =0;i<16;i++)
        if(i>total_nr_args-1) values_2[i] = 0.0;
      PlotArrays(values_1,values_2,values_3,values_4,"");
    }

	genestats.clear();	
	toxcount.clear();
	argcount.clear();
  if(save)
  {
    ofstream file;
    stringstream summary;
    summary  << folder << "/Summary_stats.dat";
    string filename = summary.str();
    file.open (filename.c_str(), ios::app);
    file << returnstring;
    file.close();
  }
  return returnstring;
}

void Dump_Genomes(TYPE2 **vibrios)
{
  bool V = FALSE;
  ofstream file;
  stringstream all_genomes;
  all_genomes  << folder << "/Genomes.dat";
  string filename = all_genomes.str();
  file.open (filename.c_str(), ios::app);
  if(Time==0) file << "Time\trow\tcol\tgenome" << endl;
  for(int row = 1; row<=nrow;row++){
  for(int col=1;col<=ncol;col++){
        if(vibrios[row][col].val > 1)
        {                    
          file << Time << "\t" << row << "\t" << col << "\t" << vibrios[row][col].G->ListContent(NULL,false,false,true) << endl;            
        }
      }
  }
  file.close();
}


void Dump_Grids(TYPE2 **vibrios, TYPE2 **dna)
{
  bool V = FALSE;
  ofstream file;
  stringstream filepath;
  filepath  << folder << "/NumTransposonGrid.dat";
  string filename = filepath.str();
  file.open (filename.c_str(), ios::app);

  ofstream file2;
  stringstream filepath2;
  filepath2  << folder << "/NumTransformantGrid.dat";
  filename = filepath2.str();
  file2.open (filename.c_str(), ios::app);

  ofstream file3;
  stringstream filepath3;
  filepath3  << folder << "/ValGrid.dat";
  filename = filepath3.str();
  file3.open (filename.c_str(), ios::app);

  ofstream file4;
  stringstream filepath4;
  filepath4  << folder << "/MobGrid.dat";
  filename = filepath4.str();
  file4.open (filename.c_str(), ios::app);

  ofstream file5;
  stringstream filepath5;
  filepath5  << folder << "/DNAGrid.dat";
  filename = filepath5.str();
  file5.open (filename.c_str(), ios::app);

  if(Time==0) file << "Time\trow\tcol\tnum_transposons" << endl;
  if(Time==0) file2 << "Time\trow\tcol\ttransformant_type" << endl;
  if(Time==0) file3 << "Time\trow\tcol\tval" << endl;
  if(Time==0) file4 << "Time\trow\tcol\tmobility" << endl;
  if(Time==0) file5 << "Time\trow\tcol\tnum_frags" << endl;
  for(int row = 1; row<=nrow;row++)
  {
  for(int col=1;col<=ncol;col++)
      {
        if(vibrios[row][col].val > 1)
        {
          float sum_mob = 0.0;
          int num_mob = 0;
          Genome::iter i = vibrios[row][col].G->StringOfPearls->begin();
          while(i!=vibrios[row][col].G->StringOfPearls->end()) 
          {
            if(vibrios[row][col].G->IsNoncoding(*i)) 
            {
              double mob = vibrios[row][col].G->GetPearlMobility(*i);
              if(mob > 0.75)
              {
                num_mob++;
                sum_mob+=mob; 
              }
            }
            i++;
          }
          float avg_mob = num_mob > 0 ? sum_mob/num_mob : 0.0;

          file << Time << "\t" << row << "\t" << col << "\t" << vibrios[row][col].G->transposases.size() << endl;
          file2 << Time << "\t" << row << "\t" << col << "\t" << vibrios[row][col].G->transformant << endl;
          file3 << Time << "\t" << row << "\t" << col << "\t" << vibrios[row][col].val << endl;
          file4 << Time << "\t" << row << "\t" << col << "\t" << avg_mob << endl;
          file5 << Time << "\t" << row << "\t" << col << "\t" << dna[row][col].DNA->Fragments->size() << endl;
        }
      }
  }
  file.close();
  file2.close();
  file3.close();
  file4.close();
  file5.close();
}

string Get_Genome_Stats(TYPE2 **vibrios)
{
    //cout << "Getting toxstats for " << total_nr_toxins << " toxins" << endl;

    int pop_size = 0;
    int sum_genome_size = 0;
    int biggest_genome = 0;
    int smallest_genome = 99999999;
    float av_genome_size =  0;

    for(int row = 1; row<=nrow;row++){
        for(int col=1;col<=ncol;col++){
            if(vibrios[row][col].val > 1)
            {
                pop_size++;
                sum_genome_size += vibrios[row][col].G->genomesize_;
                if (vibrios[row][col].G->genomesize_ > biggest_genome) biggest_genome = vibrios[row][col].G->genomesize_;
                if (vibrios[row][col].G->genomesize_ < smallest_genome) smallest_genome = vibrios[row][col].G->genomesize_; //LS: else if?
            }
        }
    }

    av_genome_size = (float)sum_genome_size/pop_size;

    stringstream genomestats;

    genomestats << setprecision(3);
    genomestats << "AVG: " << av_genome_size << "\t";
    genomestats << "MAX: " << biggest_genome << "\t";
    genomestats << "MIN: " << smallest_genome << "\t";
    genomestats << endl;

    string returnstring;
    returnstring+=genomestats.str();
    genomestats.clear();

    return returnstring;
}

void PrintMostAbundant(TYPE2** vibrios)
{
	bool V = FALSE;
  map<string,int> m;

	vector < string > all_species;                                                //LS: what does this code do? Makes a vector of type string named all_species
	if(V) cout << "Storing all species-strings in vector" << endl;
	for (int row=1; row <= nrow; row++)
	{
		for (int col=1; col <= ncol; col++)
		{
			if(vibrios[row][col].val > 1)
            {
				          if(vibrios[row][col].G->transposases.size()>0) all_species.push_back(vibrios[row][col].G->ListContent(NULL, FALSE, FALSE, FALSE));    //LS: ??
            }
		}
	}
	//if(V)
	//	for (vector<string>::iterator it=all_colours.begin(); it!=all_colours.end(); ++it)        //LS: to determine which color is who?
    //		cout << ' ' << *it;
	//cout << endl;
  int max=0;
  string most_abundant = "";
  string ditgenoom = "";
  for(int i = 0; i < all_species.size(); i++)
  {
    ditgenoom = all_species[i];
    m[ditgenoom]++;
    if(m[ditgenoom] > max)
    {
      most_abundant = ditgenoom;
      max = m[ditgenoom];
    }
    break; // ! Added because for now just a random genome is shown. 
  }

	if(V) cout << "Getting abundant" << endl;
  if(all_species.size()==0) 
  {
    cout << "No cell found with a transposon" << endl;
    return;
  }
    cout << most_abundant;
    cout << setprecision(5) << " and it's frequency in the population is " << (float)max/all_species.size() << endl;     //LS: even gets printed when this bool is FALSE, because not if(V) in same line?

	ditgenoom.clear();
	most_abundant.clear();
	m.clear();
    all_species.clear();
}

string PrintAverageToxinMobilities(TYPE2** vibrios)

  {
      bool V = FALSE;

      //DECLARE PLACEHOLDER FOR TOX POINTER
         Toxin* tox;
         vector<int> num_of_toxins;
         vector<float> sum_tox_mobility;
         num_of_toxins.assign(total_nr_args, 0);
         sum_tox_mobility.assign(total_nr_args, 0.0);
         int i;

      //int biggest_mobility = 0;
      //int smallest_mobility = 99999999;


      for(int row = 1; row<=nrow;row++){
          for(int col=1;col<=ncol;col++){
              if(vibrios[row][col].val > 1)
              {
                if (V) cout <<" finds alive vibrio" << endl;
                    // ITERATING STUFF
                list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
                Genome::iter i = StringOfPearls->begin();
                    while(i!=StringOfPearls->end())
                    {
                    if(vibrios[row][col].G->IsToxin(*i))
                		{
                      if (V) cout << "found toxin" << endl;
                		tox=dynamic_cast<Toxin*> (*i);
                    sum_tox_mobility[tox->type] += tox->mobility;
                    num_of_toxins[tox->type] ++;
                    }
                    i++;
                    }
          }
      }
  }

  stringstream toxmobstats;

  toxmobstats << setprecision(3);
  toxmobstats << "Avgmobtox:\t";
  for (int i=0;i<total_nr_args;i++)
  {
  //float gemiddelde = 1/0;
  //if(num_of_toxins[i]>0)
    float gemiddelde_tox = (float)sum_tox_mobility[i]/num_of_toxins[i];
    toxmobstats << gemiddelde_tox << "\t";
  }
  toxmobstats << endl;

  string returnstringtoxmob;
  returnstringtoxmob+=toxmobstats.str();

  toxmobstats.clear();
  num_of_toxins.clear();
  sum_tox_mobility.clear();

  return returnstringtoxmob;
}


string PrintAverageARGMobilities(TYPE2** vibrios)

  {
      bool V = FALSE;

      //DECLARE PLACEHOLDER FOR TOX POINTER
         ARG* arg;
         vector<int> num_of_args;
         vector<float> sum_arg_mobility;
         num_of_args.assign(total_nr_args, 0);
         sum_arg_mobility.assign(total_nr_args, 0.0);
         int i;

      //int biggest_mobility = 0;
      //int smallest_mobility = 99999999;


      for(int row = 1; row<=nrow;row++){
          for(int col=1;col<=ncol;col++){
              if(vibrios[row][col].val > 1)
              {
                if (V) cout <<" finds alive vibrio" << endl;
                    // ITERATING STUFF
                list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
                Genome::iter i = StringOfPearls->begin();
                    while(i!=StringOfPearls->end())
                    {
                    if(vibrios[row][col].G->IsARG(*i))
                		{
                      if (V) cout << "found toxin" << endl;
                		arg=dynamic_cast<ARG*> (*i);
                    sum_arg_mobility[arg->type] += arg->mobility;
                    //cout << arg->mobility << endl;
                    num_of_args[arg->type] ++;
                    }
                    i++;
              }
          }
      }
  }
  

  stringstream argmobstats;

  argmobstats << setprecision(3);
  argmobstats << "AvgARGmob:\t";
  for (int i=0;i<total_nr_args;i++)
  {
  //float gemiddelde = 1/0;
  //if(num_of_toxins[i]>0)
    float gemiddelde_arg = (float)sum_arg_mobility[i]/num_of_args[i];
    argmobstats << gemiddelde_arg << "\t";
  }
  argmobstats << endl;
  string returnstringargmob;
  returnstringargmob+=argmobstats.str();

  argmobstats.clear();
  num_of_args.clear();
  sum_arg_mobility.clear();

  return returnstringargmob;
}




string PrintAverageNCMobilities(TYPE2** vibrios)

  {
      bool V = FALSE;

      //DECLARE PLACEHOLDER FOR TOX POINTER
         Noncoding* nc;
        double sum_nc_mob = 0.0;
        int sum_nc = 0;
         int i;

      //int biggest_mobility = 0;
      //int smallest_mobility = 99999999;


      for(int row = 1; row<=nrow;row++){
          for(int col=1;col<=ncol;col++){
              if(vibrios[row][col].val > 1)
              {
                if (V) cout <<" finds alive vibrio" << endl;
                    // ITERATING STUFF
                list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
                Genome::iter i = StringOfPearls->begin();
                    while(i!=StringOfPearls->end())
                    {
                    if(vibrios[row][col].G->IsNoncoding(*i))
                		{
                     
                		nc=dynamic_cast<Noncoding*> (*i);
                    sum_nc_mob += nc->mobility;
                    //cout << arg->mobility << endl;
                    sum_nc ++;
                    }
                    i++;
              }
          }
      }
  }
  

  stringstream mobstats;

  mobstats << setprecision(3);
  mobstats << "AvgNCmob:\t";
  
  mobstats << sum_nc_mob / sum_nc << endl;
  string returnstring;
  returnstring+=mobstats.str();

  mobstats.clear();

  return returnstring;
}




string PrintAverageCoreMobilities(TYPE2** vibrios)
  {
      bool V = FALSE;

      //DECLARE PLACEHOLDER FOR TOX POINTER
         Core* core;
         vector<int> num_coregenes;
         vector<float> sum_core_mobility;
         num_coregenes.assign(init_nr_coregenes, 0);
         sum_core_mobility.assign(init_nr_coregenes, 0.0);
         int i;

      //int biggest_mobility = 0;
      //int smallest_mobility = 99999999;


      for(int row = 1; row<=nrow;row++){
          for(int col=1;col<=ncol;col++){
              if(vibrios[row][col].val > 1)
              {
                
                list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
                Genome::iter i = StringOfPearls->begin();
                    while(i!=StringOfPearls->end())
                    {
                    if(vibrios[row][col].G->IsCore(*i))
                		{
                		core=dynamic_cast<Core*> (*i);
                    sum_core_mobility[core->type] += core->mobility;
                    //cout << arg->mobility << endl;
                    num_coregenes[core->type] ++;
                    }
                    i++;
              }
          }
      }
  }

  stringstream coremobstats;

  coremobstats << setprecision(3);
  coremobstats << "AvgCOREmob:\t";
  for (int i=0;i<init_nr_coregenes;i++)
  {
  //float gemiddelde = 1/0;
  //if(num_of_toxins[i]>0)
    float gemiddelde_core = (float)sum_core_mobility[i]/num_coregenes[i];
    coremobstats << gemiddelde_core << "\t";
  }
  coremobstats << endl;
  string returnstringcoremob;
  returnstringcoremob+=coremobstats.str();

  coremobstats.clear();
  num_coregenes.clear();
  sum_core_mobility.clear();

  return returnstringcoremob;

  }

string GetAverageAntiToxinMobility(TYPE2** vibrios)  // UNDER CONSTRUCTION! trying to figure out how to get the average toxin mobility, so need every toxin from every grid

  {
      bool V = FALSE;
      //cout << "Getting toxstats for " << total_nr_toxins << " toxins" << endl;



      //DECLARE PLACEHOLDER FOR TOX POINTER
         Antitoxin* antitox;
         vector<int> num_of_antitoxins;
         vector<float> sum_antitox_mobility;
         num_of_antitoxins.assign(total_nr_args, 0);
         sum_antitox_mobility.assign(total_nr_args, 0.0);
         int i;

      //int biggest_mobility = 0;
      //int smallest_mobility = 99999999;


      for(int row = 1; row<=nrow;row++){
          for(int col=1;col<=ncol;col++){
              if(vibrios[row][col].val > 1)
              {
                if (V) cout <<" finds alive vibrio" << endl;
                    // ITERATING STUFF
                list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
                Genome::iter i = StringOfPearls->begin();
                    while(i!=StringOfPearls->end())
                    {
                    if(vibrios[row][col].G->IsAntitoxin(*i))
                		{
                      if (V) cout << "found antitoxin" << endl;
                		antitox=dynamic_cast<Antitoxin*> (*i);
                    sum_antitox_mobility[antitox->type] += antitox->mobility;
                    num_of_antitoxins[antitox->type] ++;
                    }
                    i++;

                  //if (vibrios[row][col].G->genomesize_ > biggest_genome) biggest_genome = vibrios[row][col].G->genomesize_;
                  //if (vibrios[row][col].G->genomesize_ < smallest_genome) smallest_genome = vibrios[row][col].G->genomesize_; //LS: else if?
              }
          }
      }
  }

  stringstream antitoxmobstats;

  antitoxmobstats << setprecision(3);
  antitoxmobstats << "Avgmobant:\t";
  for (int i=0;i<total_nr_args;i++)
  {
  //float gemiddelde = 1/0;
  //if(num_of_toxins[i]>0)
    float gemiddelde_antitox = (float)sum_antitox_mobility[i]/num_of_antitoxins[i];
  antitoxmobstats << gemiddelde_antitox << "\t";
}
  antitoxmobstats << endl;
  string returnstring_antitoxmob;
  returnstring_antitoxmob+=antitoxmobstats.str();
  antitoxmobstats.clear();
  num_of_antitoxins.clear();
  sum_antitox_mobility.clear();

  return returnstring_antitoxmob;
}

void WriteTACMobilities(TYPE2** vibrios)
{
  bool V = FALSE;
  Antitoxin* antitox;
  Toxin* tox;
  Core* core;
    vector< vector < float > > antitoxmobility; // LS: This worked before, why not anymore?
    antitoxmobility.resize(total_nr_args);
    vector< vector < float > > toxmobility;
    toxmobility.resize(total_nr_args);
	vector< vector < float > > coremobility;
	coremobility.resize(init_nr_coregenes);

    for(int row = 1; row<=nrow;row++){
      for(int col=1;col<=ncol;col++){
        if(vibrios[row][col].val > 1)
        {
          if (V) cout <<" finds alive vibrio" << endl;
          // ITERATING STUFF
          list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
          Genome::iter i = StringOfPearls->begin();
          while(i!=StringOfPearls->end())
          {
            if(vibrios[row][col].G->IsToxin(*i))
            {
              if (V) cout << "found toxin" << endl;
              tox=dynamic_cast<Toxin*> (*i);
              toxmobility[tox->type].push_back(tox->mobility);
            }
            else if(vibrios[row][col].G->IsAntitoxin(*i))
            {
              if (V) cout << "found antitoxin" << endl;
              antitox=dynamic_cast<Antitoxin*> (*i);
              antitoxmobility[antitox->type].push_back(antitox->mobility);
            }
			else if(vibrios[row][col].G->IsCore(*i))
            {
              if (V) cout << "found coregene" << endl;
              core=dynamic_cast<Core*> (*i);
              coremobility[core->type].push_back(core->mobility);
            }
            i++;

            //if (vibrios[row][col].G->genomesize_ > biggest_genome) biggest_genome = vibrios[row][col].G->genomesize_;
            //if (vibrios[row][col].G->genomesize_ < smallest_genome) smallest_genome = vibrios[row][col].G->genomesize_; //LS: else if?
          }
        }
      }
    }


    int nr_mobilities_kept = 200;

    for(int q = 0; q<total_nr_args; q++) // Looped over alle TYPE toxins (0-total_nr_toxins-1)
    {

      ofstream toxmobdatafile;

      stringstream path_tox;
      path_tox  << toxin_mobility_folder << "/Toxin" << q << ".dat";
      string filename = path_tox.str();
      toxmobdatafile.open (filename.c_str(), ios::app);

      //Shuffeling vectors and gaining only first 200 entries to not blow up the save data (RANDOM ENOUGH???)
      random_shuffle(toxmobility[q].begin(), toxmobility[q].end()); // is this correctly notated?
      int range = min((int)toxmobility[q].size(),nr_mobilities_kept);
      toxmobdatafile << Time << ", " ;

      for(int r = 0; r < range; r++)  // Loop over alle values in dit type
      {
      //  cout << " savin toxin mobilities of type " << q << " max amount saved" << nr_mobilities_kept << endl;
        toxmobdatafile << toxmobility[q][r];
        if(r < range-1) toxmobdatafile << ", ";

      }
      toxmobdatafile << endl;

      toxmobdatafile.close();
    }

    for(int q = 0; q<total_nr_args; q++) // Looped over alle TYPE antitoxins (0-total_nr_toxins-1)
    {
      ofstream antitoxmobdatafile;

      stringstream path_anti;
      path_anti  << antitoxin_mobility_folder << "/Antitoxin" << q << ".dat";
      string filename = path_anti.str();
      antitoxmobdatafile.open (filename.c_str(), ios::app);

      random_shuffle(antitoxmobility[q].begin(), antitoxmobility[q].end()); // is this correctly notated?
      int range = min((int)antitoxmobility[q].size(),nr_mobilities_kept);
      antitoxmobdatafile << Time << ", " ;

      for(int r = 0; r < range; r++)  // Loop over alle values in dit type
      {
        //cout << " savin toxin mobilities of type " << q << " max amount saved" << nr_mobilities_kept << endl;
        antitoxmobdatafile << antitoxmobility[q][r];
        if(r < range-1) antitoxmobdatafile << ", ";

      }
      antitoxmobdatafile << endl;

      antitoxmobdatafile.close();
    }
	for(int q = 0; q<init_nr_coregenes; q++) // Looped over alle TYPE coregenes (0-nr coregenes-1)
    {

      ofstream coremobdatafile;

      stringstream path_core;
      path_core  << coregene_mobility_folder << "/Coregene" << q << ".dat";
      string filename = path_core.str();
      coremobdatafile.open (filename.c_str(), ios::app);

      //Shuffeling vectors and gaining only first 200 entries to not blow up the save data (RANDOM ENOUGH???)
      random_shuffle(coremobility[q].begin(), coremobility[q].end()); // is this correctly notated?
      int range = min((int)coremobility[q].size(),nr_mobilities_kept);
      coremobdatafile << Time << ", " ;

      for(int r = 0; r < range; r++)  // Loop over alle values in dit type
      {
        //cout << " savin core mobilities of type " << q << " max amount saved" << nr_mobilities_kept << endl;
        coremobdatafile << coremobility[q][r];
        if(r < range-1) coremobdatafile << ", ";

      }
      coremobdatafile << endl;

      coremobdatafile.close();
    }

	antitoxmobility.clear(); // LS: This worked before, why not anymore?
	toxmobility.clear();
	coremobility.clear();
}

void  SampleFrequencies(TYPE2** vibrios)
{
bool V = FALSE;

if (V) cout << "at least i make it in the function" <<endl;

int nr_frequencies_kept = 1000;

Antitoxin* antitox;
Toxin* tox;
Core* core;

vector< vector < float > > antitoxfrequency;
antitoxfrequency.resize(total_nr_args);
vector< vector < float > > toxfrequency;
toxfrequency.resize(total_nr_args);
vector< vector < float > > corefrequency;
corefrequency.resize(init_nr_coregenes);

	for(int row = 1; row<=nrow;row++){
		for(int col=1;col<=ncol;col++){
			if(vibrios[row][col].val > 1)
			{
				if (V){
					cout << " gets in alive vibrio" << endl;
					cout << total_nr_args << endl;
					cout << init_nr_coregenes << "\n";
				}
			int toxcounter_[total_nr_args]= {0};
			int antcounter_[total_nr_args]= {0};
			int corcounter_[0];
			if (init_nr_coregenes >0) int corcounter_[init_nr_coregenes]= {0};
				if (V)cout <<" i made the counters" << endl;
			list<Pearl*>* StringOfPearls = vibrios[row][col].G->StringOfPearls;
            Genome::iter i = StringOfPearls->begin();
	            while(i!=StringOfPearls->end())
	            {
				if (V)cout << " on a pearl" << endl;
				bool V1 = FALSE;
					if(vibrios[row][col].G->IsToxin(*i))
		            {
		                if (V1) cout << "found toxin" << endl;
		                tox=dynamic_cast<Toxin*> (*i);
		                toxcounter_[tox->type] += 1;
						if (V1) cout << "cell" <<row<<" "<< col <<"added to toxcounter of type"<< tox->type <<"and is now " << toxcounter_[tox->type] << endl;
		            }
		            else if(vibrios[row][col].G->IsAntitoxin(*i))
		            {
		                if (V1) cout << "found antitoxin" << endl;
		                antitox=dynamic_cast<Antitoxin*> (*i);
		                antcounter_[antitox->type]+=1;
						if (V1) cout << "cell" <<row<<" "<< col << "added to antcounter of type"<< antitox->type <<"and is now " << antcounter_[antitox->type] << endl;
		            }
		  			else if(vibrios[row][col].G->IsCore(*i))
		            {
		                if (V1) cout << "found coregene" << endl;
		                core=dynamic_cast<Core*> (*i);
						corcounter_[core->type]+=1;
						if (V1) cout << "cell" <<row<<" "<< col <<"added to corcounter of type"<< core->type <<"and is now " << corcounter_[core->type] << endl;
		            }
	            i++;
			  	}

			for(int p = 0; p<total_nr_args; p++){
			toxfrequency[p].push_back(toxcounter_[p]);
			}
			for(int q = 0; q<total_nr_args; q++){
			antitoxfrequency[q].push_back(antcounter_[q]);
			}
			for(int r = 0; r<init_nr_coregenes; r++){
			corefrequency[r].push_back(corcounter_[r]);
			}
			}
		}
	}
	for(int q = 0; q<total_nr_args; q++) // Looped over alle TYPE antitoxins (0-total_nr_toxins-1)
	{
	  ofstream toxfreqdatafile;

	  stringstream path_toxfreq;
	  path_toxfreq  << frequencies_folder << "/Toxin" << q << ".dat";
	  string filename = path_toxfreq.str();
	  toxfreqdatafile.open (filename.c_str(), ios::app);

	  random_shuffle(toxfrequency[q].begin(), toxfrequency[q].end()); // is this correctly notated?
	  int range = min((int)toxfrequency[q].size(),nr_frequencies_kept);
	  toxfreqdatafile << Time << ", " ;

	  for(int r = 0; r < range; r++)  // Loop over alle values in dit type
	  {
		//cout << " savin toxin mobilities of type " << q << " max amount saved" << nr_mobilities_kept << endl;
		toxfreqdatafile << toxfrequency[q][r];
		if(r < range-1) toxfreqdatafile << ", ";

	  }
	  toxfreqdatafile << endl;

	  toxfreqdatafile.close();
	}
	for(int q = 0; q<total_nr_args; q++) // Looped over alle TYPE antitoxins (0-total_nr_toxins-1)
	{
	  ofstream antitoxfreqdatafile;

	  stringstream path_antfreq;
	  path_antfreq  << frequencies_folder << "/Antitoxin" << q << ".dat";
	  string filename = path_antfreq.str();
	  antitoxfreqdatafile.open (filename.c_str(), ios::app);

	  random_shuffle(antitoxfrequency[q].begin(), antitoxfrequency[q].end()); // is this correctly notated?
	  int range = min((int)antitoxfrequency[q].size(),nr_frequencies_kept);
	  antitoxfreqdatafile << Time << ", " ;

	  for(int r = 0; r < range; r++)  // Loop over alle values in dit type
	  {
		//cout << " savin toxin mobilities of type " << q << " max amount saved" << nr_mobilities_kept << endl;
		antitoxfreqdatafile << antitoxfrequency[q][r];
		if(r < range-1) antitoxfreqdatafile << ", ";

	  }
	  antitoxfreqdatafile << endl;

	  antitoxfreqdatafile.close();
	}
	for(int q = 0; q<init_nr_coregenes; q++) // Looped over alle TYPE antitoxins (0-total_nr_toxins-1)
	{
	  ofstream corefreqdatafile;

	  stringstream path_corfreq;
	  path_corfreq  << frequencies_folder << "/Core" << q << ".dat";
	  string filename = path_corfreq.str();
	  corefreqdatafile.open (filename.c_str(), ios::app);

	  random_shuffle(corefrequency[q].begin(), corefrequency[q].end()); // is this correctly notated?
	  int range = min((int)corefrequency[q].size(),nr_frequencies_kept);
	  corefreqdatafile << Time << ", " ;

	  for(int r = 0; r < range; r++)  // Loop over alle values in dit type
	  {
		//cout << " savin toxin mobilities of type " << q << " max amount saved" << nr_mobilities_kept << endl;
		corefreqdatafile << corefrequency[q][r];
		if(r < range-1) corefreqdatafile << ", ";

	  }
	  corefreqdatafile << endl;

	  corefreqdatafile.close();
	}

	antitoxfrequency.clear();
	toxfrequency.clear();
	corefrequency.clear();

}

const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[120];
    tstruct = *localtime(&now);
    sprintf(buf, "%d-%02d-%02d_%02dh%02dm%02ds", 1900+tstruct.tm_year, tstruct.tm_mon+1, tstruct.tm_mday, tstruct.tm_hour, tstruct.tm_min, tstruct.tm_sec);
    return buf;
}
