extern int nrow;
extern int ncol;
extern int Time;
extern bool Graphs;
extern bool extinct;

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
    int num_ne_pearls=0;
    int sum_time_infected = 0;
    int num_infected = 0;
    Noncoding* nc;
    Transposase* tra;

    double sum_mobility = 0.0;   
    //double sum_IR_mobility = 0.0;
    
    for(int row = 1; row<=nrow;row++){
        for(int col=1;col<=ncol;col++){
            if(vibrios[row][col].val >= 1)
            {
                pop_size++;
                sumtra += vibrios[row][col].G->transposases.size();
                if(vibrios[row][col].G->time_infected_ > 0) 
                {
                  sum_time_infected += vibrios[row][col].G->time_infected_;
                  num_infected++;
                }
                sumtrans += vibrios[row][col].G->transformant==1;
                sumtrans2 += vibrios[row][col].G->transformant==2;
		            sumfit += vibrios[row][col].G->compstrength;
                sum_genomesize += vibrios[row][col].G->StringOfPearls->size();

                Genome::iter i = vibrios[row][col].G->StringOfPearls->begin();
                while(i!=vibrios[row][col].G->StringOfPearls->end())
                {
                  if(vibrios[row][col].G->IsNoncoding(*i))
                  {                  
                    nc=dynamic_cast<Noncoding*> (*i);         
                    num_nc_pearls ++;
                    
                  }                  
                  else if(vibrios[row][col].G->IsHK(*i))
                  {
                    sum_hgt_hk+= (*i)->num_horizontal_transfers_;
                    sum_jump_hk+= (*i)->num_jumps_;
                    num_hk_pearls++;
                  }
                  else if(vibrios[row][col].G->IsNonEssential(*i))
                  {
                    num_ne_pearls++;
                  }
                  else if(vibrios[row][col].G->IsTransposase(*i))
                  {
                    tra=dynamic_cast<Transposase*> (*i);
                    sum_mobility+=tra->mobility;

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

    
    float fract_TEs = (float)num_TEs/pop_size;
    float fract_tra = (float)sumtra/pop_size;
    float fract_ne = (float)num_ne_pearls/pop_size;
    float fract_nc = (float)num_nc_pearls/pop_size;
    float fract_trans = (float)sumtrans/pop_size;
    float fract_trans2 = (float)sumtrans2/pop_size;

    float fract_hgt_hk = num_hk_pearls > 0 ? (float)sum_hgt_hk/num_hk_pearls: 0.0;
    float fract_hgt_tra = sumtra > 0 ? (float)sum_hgt_tra/sumtra : 0.0;

    float fract_jump_hk = num_hk_pearls > 0 ? (float)sum_jump_hk/num_hk_pearls: 0.0;
    float fract_jump_tra = sumtra > 0 ? (float)sum_jump_tra/sumtra : 0.0;
  
    if(Time==0)
    {
      //genestats << ">";
      genestats << "Time\tAVG_g\tAVG_f\tPopsize\tFracHK\tFracNE\tFracTra\tFractNC\tHGT_hk\tJump_hk\tHGT_tra\tJump_tra\tPhi\tTime_infected\t";
      
      genestats << "\n";
		}
    
    genestats << fixed << setprecision(6) << Time << "\t" << (float)sum_genomesize/pop_size << "\t";
    
    genestats << sumfit/pop_size << "\t" << pop_size << "\t" << (float)num_hk_pearls/pop_size << "\t" << fract_ne << "\t" << fract_tra << "\t" << fract_nc << "\t";
    genestats << fract_hgt_hk << "\t" << fract_jump_hk << "\t" << fract_hgt_tra << "\t" << fract_jump_tra << "\t" << (float)sum_mobility/sumtra << "\t" << (float)sum_time_infected/num_infected << "\t";

    genestats << endl;
    string returnstring;
    
    returnstring+=genestats.str();
    if(Graphs)
    {
      double values_1[16] = {(float)pop_size,nrow*ncol+100,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double values_2[16] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double values_3[16] = {0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double values_4[16] = {fract_trans,fract_tra,fract_trans2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      
      PlotArrays(values_1,values_2,values_3,values_4,"");
    }

	genestats.clear();	
	
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
  if(pop_size == 0 || fract_tra == 0.0) 
  { 
      // cout << "Shutting down simulation as no more cells are alive or transposons died out" << endl; 
      extinct = true;
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
          file << Time << "\t" << row << "\t" << col << "\t" << vibrios[row][col].G->ListContentShort(NULL,false,false,false) << endl;            
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

  ofstream file6;
  stringstream filepath6;
  filepath6  << folder << "/NumTransposonExt.dat";
  filename = filepath6.str();
  file6.open (filename.c_str(), ios::app);

  ofstream file7;
  stringstream filepath7;
  filepath7  << folder << "/Gsize.dat";
  filename = filepath7.str();
  file7.open (filename.c_str(), ios::app);
  

  if(Time==0) file << "Time\trow\tcol\tnum_transposons" << endl;
  if(Time==0) file2 << "Time\trow\tcol\ttransformant_type" << endl;
  if(Time==0) file3 << "Time\trow\tcol\tval" << endl;
  if(Time==0) file4 << "Time\trow\tcol\tmobility" << endl;  
  if(Time==0) file5 << "Time\trow\tcol\tnum_frags" << endl;
  if(Time==0) file6 << "Time\trow\tcol\tnum_transposons\tnum_pearls" << endl;
  if(Time==0) file7 << "Time\trow\tcol\tgenomesize" << endl;
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
          file7 << Time << "\t" << row << "\t" << col << "\t" << vibrios[row][col].G->StringOfPearls->size() << endl;
        }
        int nr_tra_ext = 0;
        int nr_pearls = 0;
        for (std::vector< std::list<Pearl*> >::iterator frit = dna[row][col].DNA->Fragments->begin(); frit != dna[row][col].DNA->Fragments->end(); frit++)
        {
          
          // cout << vibrios[row][col].G->ListContent(&(*frit)) << endl;
           std::list<Pearl*>::iterator it = frit->begin();
           while(it!=frit->end())
           {
            //  cout << "." << endl;
             if(dna[row][col].DNA->IsTransposase(*it)) nr_tra_ext++;
             nr_pearls++;
             it++;             
          }
        }
        file6 << Time << "\t" << row << "\t" << col << "\t" << nr_tra_ext << "\t" << nr_pearls << endl;
        
      }
  }
  file.close();
  file2.close();
  file3.close();
  file4.close();
  file5.close();
  file6.close();
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
				          all_species.push_back(vibrios[row][col].G->ListContent(NULL, FALSE, FALSE, FALSE));    //LS: ??
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
  
    cout << most_abundant;
    cout << setprecision(5) << " and it's frequency in the population is " << (float)max/all_species.size() << endl;     //LS: even gets printed when this bool is FALSE, because not if(V) in same line?

	ditgenoom.clear();
	most_abundant.clear();
	m.clear();
    all_species.clear();
}



const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[120];
    tstruct = *localtime(&now);
    sprintf(buf, "%d-%02d-%02d_%02dh%02dm%02ds", 1900+tstruct.tm_year, tstruct.tm_mon+1, tstruct.tm_mday, tstruct.tm_hour, tstruct.tm_min, tstruct.tm_sec);
    return buf;
}
