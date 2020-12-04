/* ====================================================================================
		 Hieronder staat getLine stukje, klik op display triggers
  =====================================================================================*/
if(Mouse() == 1){
	gotMouseBrem();

		if(getSummary==1)
			cout << "Summary unavailable due to december changes" << endl; // Summary unavail
			getSummary=0;
		if(changegene==1){
			char answer[100];
			printf("Change gene to...? (0-39)\n");
			fgets(answer,100,stdin);
			genedisplay = atoi(answer);		
			changegene=0;
			rewind(stdin);
		}
		if(switchmix==1) {
			if(!mixing){mixing = 1;cout << "Mixing switched ON" << endl;}
			else {mixing = 0;cout << "Mixing switched OFF" << endl;}
			switchmix =0;
		}
		if(switchdif==1) {
			if(!diffusing){diffusing = 1; cout << "Diffusion switched ON" << endl;}
			else{ diffusing = 0; cout << "Diffusion switched OFF" << endl;}
			switchdif = 0;
		}
		if(dumpfile==1) {
			//DUMP
		}
		if(slideshow ==1){
			makemovie = 1; 
			cout << "Saving movie in dir /timeframes" << endl;
			slideshow = 0;
		}
		if(smutation == 1){
			if(!mutation){mutation = 1; cout << "Mutation switched ON" << endl;}
			else{ mutation = 0; cout << "Mutation switched OFF" << endl;}
			smutation = 0;
		}
		if(makemovie && moviecounter < 1000){
                	DrawSlide(Vibrios, "timeframes");
                	moviecounter++;
	        }
       		else{
                	moviecounter = 0;
        	}

	}

        if(diffusing) MDiffusion(Vibrios);

        if(mixing) PerfectMix(Vibrios);
	
	// EINDE GET-LINE
