/***********************************************************************************************
This code simulates the progression of HIV patients undergoing HAART.  

v1.0 -- 3/23/05. Calibrated off the CHORUS database.  Developed in C by Steven Shechter over the course of a year.  This initial model was originally created by Scott Braithwaite in DecisionMaker software.  We decided to move it to C for the flexibility and dramatic reduction in run time.
v.1.1 -- 3/28/05. Made some changes/fixes to how we process patients who start with initial mutations.  Added #if RES_BASED_CHOICE
that determines if we want to make an informed choice of the first regimen (presumably based in real life on results	of resistance tests) by setting that to 1, or we can set it to 0 and not "optimize" the starting regimen for patients starting with muts.  In either case, we run such patients through code to determine if they start off with resistance.
v.1.2 -- 3/28/05.  Changed the way we obtain compliance.  We used to always have it as: Max(0, comp + (1-comp)*(pat->diligence - .5)). 
This implicitly assumed that comp was >= .5 and we wanted a symmetric uniform distribution around comp that didn't get cut off by the upper bound of 1.  Now that we are doing more testing of lower adherence levels, we should have a corresponding construct if comp < .5: comp + (0-comp)*(pat->diligence - .5).
v.1.3 -- 4/1/05. Accidentally removed part of code that keeps patients on haart even if they start haart and rise above the start threshold level, so I reinserted that.  Also, cleaned up code to handle the situation where patients stay on haart even if they have hit the max number of regimens or experience resistance to all possible regimens.  This used to be indicated by "usefitness", but I gave it a more general name of STAY_ON_HAART.  We now allow for continued buildup of mutations in the face of resistance.  This is implemented through changing the part that assigns a 0 mutation rate if the totalresnc = 3.  Now we assign 0 only if the totalnc = 0.  e.g., if there are 3 resistant mutations, we still consider mutations to acquire. Similar consideration was given in a change to the function Get_Num_Each_Type.
v.1.4 -- 5/4/05.  Implemented/Activated functionality to incorporate utilities, costs, and discounting.
************************************************************************************************/

/*This is the version for India. Generated in December 2014. QZ*/

/*This would be the second version for India. We add simplified competing risks to identify the causes of non-HIV deaths. 2/25/2015 QZ*/

//kimnew
#define _CRT_SECURE_NO_DEPRECATE 1	//eliminates printf deprecation warnings
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include "hiv.h"
#include "constants.h"
#include "distributions.h" 
#include "perlin.h"
#include "constantsHeader.h"
#include "output.h"

void Initialize(patient*);
void Setup();
double Get_Rate(double table[TABLE_SPACE][2], int, double);
void Process_Death(int type, patient *pat);
double HIV_Lookup(patient *, bool haart);
void Process_Patient(patient*);
void Check_For_Mutations(patient*);

void Change_Regimen(patient *pat);
void Get_Compliance(patient *pat);
double Get_Median(double arr[], int);
void QuickSort(double arr[], int, int);
void Process_Output();
void Process_Initial_Mutations(patient *pat);
void Swap( double[], int, int );
void Initialize_Vars();
int Get_Category(patient *pat, int type);
double Utility(patient *pat);
double Cost(patient *pat, int whichCost);
int ReadSensitivityFile ();
int ReadSensitivityFile_RIC ();
int ChooseRegimen_Rich(patient *pat);
int ChooseRegimen_Poor(patient *pat);
void determineResistanceProfile_Poor(patient *pat);
void printDrugsInReg(patient *pat);
void printNumRemaining(patient *pat);
void getVLdec(patient *pat);
int RandomNumber(int upperRange);
void decrementForIntolerance(patient *pat);
void getInitialAIDSRates();
void Get_Num_Each_Type_In_Array(patient *pat, int numInArray, int drugArray[], int numEachType[]);
int numNotResOrIntol(patient *pat, int type);
int numNotIntol(patient *pat, int type);
int numNotRes (patient *pat, int type);
void getCD4regressionideal(patient *pat);
double addNoise(double input, patient *pat, int type);
void setNextMonitorCycle(patient *pat);
bool metTrigger(int trigger, patient *pat);
void Least_Squared(double mod[], double clin[], int size, double *d);
bool notEnoughDrugsForInitialReg(patient *pat);
int countResDrugs(patient *pat);
double linearInterpolate(double x0, double y0, double  x1, double  y1, double x);
void readCommandLineArgs(int argc, char* argv[]);
double getVLTerm(patient *pat);
void putPatOnReg(patient *pat, int type[]);
int giveMeBestDrugInClass(patient *pat, int type, int numAssignedSoFar);
void checkForCrossResAndCrossClassRes(patient *pat, int mutType);
int getNumRes(drug *reg[]);
int getNumIntol(drug *reg[]);
double Get_HIV_Death_Rate(patient *pat);
void getstate(patient *pat, int mystate[]);
void recordChange();
void showstates();
void writeTransitionsToFile();
int getResistanceCategory(patient *pat);
void incrementWHOArray(patient *pat, int numOfSomething[3]);
void WHO_print_year_by_year_stats();
void initializeWHOstats();
void RyanWhiteSetupPatient();
void set_regimen_and_triggers_WHO_discordance();
void set_regimen_and_triggers_WHO();
void record_discordance_data(patient *pat);
void process_WHO_discordance_output();
string IntToStr( int n );
void collectRyanWhiteStatsAtPointsInTime(patient *pat);
void start_haart(patient *pat);

#if RETENTION_IN_CARE
void Check_Engagement (patient*);
void check_for_AIDS_defining_events_RIC (patient *pat);
void process_HIV_death_RIC(patient *pat);
void process_age_death_RIC(patient *pat);
void addOutreachCostsToTotalCostAtDeath(patient *pat);
#endif

void graph_trajectories(patient *pat);

//for VRT
void fillUnifInit();
void fillUnifSim();
void printUnif();


int main(int argc, char *argv[])
{
	patient *pat;
	
	Setup(); //opens files, sets RNG seed, reads in male, female non-HIV-related mort tables and 1 HIV-related mort table

	if (BATCH) readCommandLineArgs(argc, argv);

#ifdef VARYING_AGE_VL_CD4
	start_age = 50;
	start_cd4 = 600;
	start_hiv = 6.0;
	cd4_treat = 350;

	avgAge_SD = 0;
	avgCD4_SD = 0;
	avgHIV_SD = 0;
#endif

#ifdef GETSTATES_LOOP
	int res_cat;
#if COMMAND_LINE_RESISTANCE_ADHERENCE
	if (argc > 2)
	{
		res_cat = atoi(argv[1]); // Command line argument 1 is the resistance level.
		comp = atof(argv[2]); // Command line argument 2 is the compliance (adherence level)
	} else {
		cout << "Arguments were not passed. Exiting." << endl;
		exit(0);
	}
#else
	res_cat = RCAT;
#endif

	for (patient_starts_on_haart = 0; patient_starts_on_haart < NUM_INITIAL_HAART_STATES; patient_starts_on_haart++)  //Run with patients starting on Care, and patients starting on Treatment
	{
		for (int start_cd4_cat = 0; start_cd4_cat<5; start_cd4_cat++)		//actual cd4s get set based on categories below
		{
			for (start_hiv = 2.0; start_hiv<6.99; start_hiv++)				//these are the actual VLs
			{
			//for (int res_cat = 0; res_cat<8; res_cat++)				//commented out because I ran 8 times on the cluster for time efficiency
			//{
			//	/*Resistant classes 
			//	0: No resistance
			//	1: NNRTI only
			//	2: PI only
			//	3: NRTI only
			//	4: NNRTI and NRTI
			//	5: PI and NRTI
			//	6: NNRTI and PI 
			//	7: All 3 classes
			//	*/

				switch (res_cat)
				{
				case 0: mut_nnrti_nevirapine_start = 0; mut_nrti_nontam_start = 0; mut_pi_boosted_start = 0; break;	 				
				case 1: mut_nnrti_nevirapine_start = 999; mut_nrti_nontam_start = 0; mut_pi_boosted_start = 0; break;	
				case 2: mut_nnrti_nevirapine_start = 0; mut_nrti_nontam_start = 0; mut_pi_boosted_start = 999; break;	
				case 3: mut_nnrti_nevirapine_start = 0; mut_nrti_nontam_start = 999; mut_pi_boosted_start = 0; break;	
				case 4: mut_nnrti_nevirapine_start = 999; mut_nrti_nontam_start = 999; mut_pi_boosted_start = 0; break;	
				case 5: mut_nnrti_nevirapine_start = 0; mut_nrti_nontam_start = 999; mut_pi_boosted_start = 999; break;	
				case 6: mut_nnrti_nevirapine_start = 999; mut_nrti_nontam_start = 0; mut_pi_boosted_start = 999; break;	
				case 7: mut_nnrti_nevirapine_start = 999; mut_nrti_nontam_start = 999; mut_pi_boosted_start = 999; break;	
				}

				switch (start_cd4_cat)
				{
				case 0: start_cd4 = 25; break;	
				case 1: start_cd4 = 125; break;	
				case 2: start_cd4 = 275; break;
				case 3: start_cd4 = 425; break;	
				case 4: start_cd4 = 705; break;	
				}

				avgCD4_SD = 0;
				avgHIV_SD = 0;
#endif

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS 

	//This call is made to ensure that the Perlin noise package gets initialized in a predictable position within the random
	// number chain, so that that parallel results are the same as the sequential.
	PerlinNoise1D (19.728000000000002, 1.5, 5, 2);

	while (ReadSensitivityFile ())
	{
		/*Reinitializing all random numbers between runs so that results run in parallel on the cluster 
		match results run in series*/
		seed = initial_seed; 
		init_genrand(seed); //sets the initial the seed
		srand( (unsigned)seed );

#endif		

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS_RIC

		//If using the cluster to distribute th sensitivity over many nodes
	#ifdef USING_COMMAND_LINE_ARGS
			if (argc > 1)
			{
				RIC_command_line_sensitivity_row_num = atoi(argv[1]);
				cout << "row_num: "<< RIC_command_line_sensitivity_row_num << endl;
			}

			else
			{
				cout<<"No command line argument was received.  Exiting."<<endl;
				exit(0);
			}
	#endif

		//This call is made to ensure that the Perlin noise package gets initialized in a predictable position within the random
		// number chain, so that that parallel results are the same as the sequential.
		PerlinNoise1D (19.728000000000002, 1.5, 5, 2);

		while (ReadSensitivityFile_RIC ())
		{
			/*Reinitializing all random numbers between runs so that results run in parallel on the cluster
			match results run in series*/
			seed = initial_seed;
			init_genrand(seed); //sets the initial the seed
			srand( (unsigned)seed );

#endif

#ifdef TRIGGER_COMPARISON
		for (trig=1; trig<=4; trig++) 
		{
			//TEST

			/*for (monitor_interval=3; monitor_interval<=12; monitor_interval+=3)
			{
				for (vl_thresh = 0; vl_thresh < 4; vl_thresh++)
					{
						if (vl_thresh==0) vl_regimen_failure_1st = vl_regimen_failure_after_1st = 2.7;
						else if (vl_thresh==1) vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
						else if (vl_thresh==2) vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.7;
						else if (vl_thresh==3) vl_regimen_failure_1st = vl_regimen_failure_after_1st = 4.0;

						//Will stop the loop after one pass if we don't need to loop
						if (!(trig==T_VL || trig==T_NESTED)) vl_thresh=3;	
						*/
#endif

#ifdef TEST_COST_COMPARISON
			trig=T_VL;
			testCostCD4 = 11.20;
			testCostVL = 70.00;
			//for (testCostCD4=1.0; testCostCD4<=12.0001; testCostCD4+=1.0)
			for (testCostVL=0.0; testCostVL<=70; testCostVL+=10.0)
			{
				
#endif
	
#ifdef MONITORING_PAPER_SENSITIVITY
				
				for (cd4_treat=200; cd4_treat> -1000; cd4_treat -= 1000)
				{
					for (numAvailableRegimens=1; numAvailableRegimens<=3; numAvailableRegimens++)
					{
						for (trig=1; trig<=4; trig++) 
						{
							if (numAvailableRegimens!=1 || (numAvailableRegimens==1 && trig==4))
							{
								for (monitor_interval=3; monitor_interval<=12; monitor_interval+=3)
								{
									if (!(trig==4 && monitor_interval >3))
									{
										for (vl_thresh = 0; vl_thresh < 2; vl_thresh++)
										{
											if (monitor_interval!=9)
											{
												if (vl_thresh==0) vl_regimen_failure_1st = vl_regimen_failure_after_1st = 2.7;
												else if (vl_thresh==1) vl_regimen_failure_1st = vl_regimen_failure_after_1st = 4.0;

												//Will stop the loop after one pass if we don't need to loop
												if (!(trig==T_VL || trig==T_NESTED)) vl_thresh=999;

												//Will stop the loop after one pass with patient not going on treatment
												if (cd4_treat<0) vl_thresh = trig = monitor_interval = trig = numAvailableRegimens = 999;

#endif
												
#if WHO_MONITORING_STRATEGY_ANALYSIS
												for (cd4_treat=200; cd4_treat<=500; cd4_treat+=150)
												{
													for (scenario=0; scenario<18; scenario++)
													{
														//scenario = SCENARIO;
														
														initializeWHOstats();
														
														set_regimen_and_triggers_WHO();
														
														
#endif
															
#ifdef WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
														cd4_treat = 350; // CHECK ON THIS!!!
														for (scenario=0; scenario<10; scenario++)
														{
															//scenario = SCENARIO;
															cout <<"scenario: "<< scenario << endl;
															set_regimen_and_triggers_WHO_discordance();
															
#endif

#ifdef RYAN_WHITE_ANALYSIS
														for (interv=0; interv<NUM_ANALYSES; interv++)
														{
															//For these runs, we want to target ALL  propTargeted = (double)numTargetedInterv[interv] / (double)TOTAL_PATS_RYAN_WHITE;
															propTargeted = 1.0;
															numTargeted = int(propTargeted * PATIENTS);
#endif

  															Initialize_Vars();

															//begin simulating patients
															for (patnum = 0; patnum < PATIENTS; patnum++) 
															{
																if (KIMTEST) fprintf(kim_test, "Beginning patient %d.\n", patnum);

																if (patnum % 500 == 0) printf("%d",patnum ); 
																else if (patnum % 100 == 0) printf(".");
																fflush (stdout);        //make sure the . goes to the consold window

																if (VRT) fillUnifInit();		

																pat = (patient *) calloc(1, sizeof(patient));
																Initialize(pat);
																
																if (KIMTEST) fprintf (kim_test, "Baseline Age: %.2f, CD4: %.2f, VL: %.2f, luck1: %.2f\n", pat->start_age, pat->CD4baseline, pat->HIVbaseline, pat->luck1);

																if (RYAN_WHITE) RyanWhiteSetupPatient();

																//run through Markov cycle until stopping criteria are met. In general, there will be many recursive calls within Process_Patient, ending when a patient either dies HIV, dies Age, or some other stopping criteria are met.
																Process_Patient(pat);
																											
																free(pat);
															} //ends for each pat

															Process_Output();

#ifdef TRIGGER_COMPARISON
														}//}}
#endif
#ifdef TEST_COST_COMPARISON
				}
#endif

#ifdef MONITORING_PAPER_SENSITIVITY
								}}}}}}}}
#endif

#if WHO_MONITORING_STRATEGY_ANALYSIS
					}}
#endif

#ifdef	WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
			}
#endif

#ifdef RYAN_WHITE_ANALYSIS
							}
#endif

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS 
	}	//ends while (ReadSensitivityFile)
#endif 

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS_RIC
            }	//ends while (ReadSensitivityFile_RIC)
#endif

#ifdef GETSTATES_LOOP
			}}}//}
#endif

	//Outputs data to file for transmission model lookup tables
#if (GETSTATES) 
	writeTransitionsToFile();
#endif

	fclose(output);
	fprintf(stateprob, "EOF\n");
	return 0;

} //end main



/*************************************************************************************************
We use the convention that a patient enters here each period, we decide if the patient is on HAART, whether they will be compliant or not to each drug of their regimen, etc.  Then we update VL and CD4 counts and use these values along with updated age and other vars to determine prob of dying from HIV or age-related death or prob of mutations and resistance developing against the different drug types patient was prescribed.
*************************************************************************************************/
void Process_Patient(patient *pat)
{
	double asr_male_rate, asr_female_rate, temp;
	int i;
	
	while (1)
	{
beginNewCycle:
		
		pat->cyclenum++;
		
		if (VRT)
		{
			fillUnifSim();          // fill array of random values using new seed
		}

		//Note that starting a patient already on haart must come before the GETSTATES code to capture that the patient is on treatment at the very beginning
		if (patient_starts_on_haart && pat->cyclenum == 0) start_haart(pat);

		if (GETSTATES)
		{
			//Set the initial state at cycle 0
			if (pat->cyclenum == 0) 
				getstate(pat, cur_state);

			//save the new state annually after the first time
			else
			{
				if (pat->cyclenum % 365 == 0) 
				{
					//The cost over the past year is today's total cost minus the total cost one year ago
					costThisYearThisPat = total_cost - totalCostAtStartOfYearThisPat;
					
					getstate(pat, new_state);
					recordChange();
				}
			}
		}
		

		/***************************************************************************************************
		At beginning of the month, we determine if the patient will be on haart for this next month.  When the CD4real falls below the treatment threshold, we put the patient on haart (as long has the patient hasn't exhausted all regimens, or if he has but we "stay on haart" which indicates we will continue to put the patient on haart because living with resistant strains is better than living with the wild type.  Has to do with viral fitness concept).  
		Note that once patient falls below the CD4_TREAT, we keep the patient on haart, even if the CD4real 	subsequently rises above the treatment threshold.  We have the variable pat->started_haart because we consider the VL/CD4 dynamics to still be affected by a haart regimen after the patient stops that regimen (through the lag effect).  
		*****************************************************************************************************/     


#if RETENTION_IN_CARE
		if (pat->engaged_in_care)   //calculate time-in-care
		{
			pat->time_in_care++;
		}

		if (!pat->started_haart && pat->CD4real <= cd4_treat && pat->engaged_in_care && (pat->willInitiateHAARTAtCD4TreatThresh || pat->CD4real <= 200))
		{ //if not the first time LTFU, once engage back, start HAART immediately
#else
		if (!pat->started_haart && pat->CD4real <= cd4_treat)
		{
#endif
			start_haart(pat);
		}


		//Increment time on current regimen
		//The ( pat->numreg < MAX_REGS_MODEL ) is to make sure we don't exceed the size of the storage array
		if ( ( pat->numreg <= MAX_REGS_MODEL ) && pat->haart )
		{
			pat->cyclesOnReg[pat->numreg-1]++;
			pat->counterreg++;
		}

		pat->age = pat->start_age + (double)pat->cyclenum/(double)CYCLES_PER_YEAR;

		if (pat->CD4Peak < pat->CD4real) pat->CD4Peak = pat->CD4real;

		if ( CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_THERAPY)
		{
			if ( pat->started_haart && 
				((pat->cyclenum - pat->haart_start_time) == 1*CYCLES_PER_YEAR ||
				(pat->cyclenum - pat->haart_start_time) == 3*CYCLES_PER_YEAR ||
				(pat->cyclenum - pat->haart_start_time) == 5*CYCLES_PER_YEAR)
				&& (pat->cyclenum/CYCLES_PER_YEAR < TIME)) 
			{	
				//set "years" which is used as the array index below
				if ((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR ==1) years = 0;
				else if ((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR ==3) years = 1;
				else if ((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR ==5) years = 2;

				numAlive_1_3_5_YearsIntoHaart[years]++;
				ChangeInCD4_1_3_5_YearsIntoTherapy[years][numAlive_1_3_5_YearsIntoHaart[years]-1] = pat->CD4real - pat->CD4atStartOfHaart;
				ChangeInVL_1_3_5_YearsIntoTherapy[years][numAlive_1_3_5_YearsIntoHaart[years]-1] = pat->VLreal - pat->VLatStartOfHaart;			
			}
		}

		if (RYAN_WHITE && patIsTargeted && pat->started_haart && ((pat->cyclenum - pat->haart_start_time) == CYCLES_PER_YEAR))
		{	if (pat->VLreal <= 2.7) numSuppressed++;
			if (pat->VLatStartOfHaart - pat->VLreal >= 1) numDroppedAtLeastOneLogVL++;
		}

		if (RYAN_WHITE) collectRyanWhiteStatsAtPointsInTime(pat); 
		

		if ( CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE)
		{
			if ( pat->exhausted_clean_regs && 
				((pat->cyclenum - pat->haart_exhausted_clean_regs_time) == 1*CYCLES_PER_YEAR ||
				(pat->cyclenum - pat->haart_exhausted_clean_regs_time) == 3*CYCLES_PER_YEAR ||
				(pat->cyclenum - pat->haart_exhausted_clean_regs_time) == 5*CYCLES_PER_YEAR)
				&& (pat->cyclenum/CYCLES_PER_YEAR < TIME)) 
			{
				if ((pat->cyclenum - pat->haart_exhausted_clean_regs_time)/CYCLES_PER_YEAR ==1) years = 0;
				else if ((pat->cyclenum - pat->haart_exhausted_clean_regs_time)/CYCLES_PER_YEAR ==3) years = 1;
				else if ((pat->cyclenum - pat->haart_exhausted_clean_regs_time)/CYCLES_PER_YEAR ==5) years = 2;

				numAlive_1_3_5_YearsIntoSalvage[years]++;
				CD4_1_3_5_YearsIntoSalvage[years][numAlive_1_3_5_YearsIntoSalvage[years]-1] = pat->CD4real;
				VL_1_3_5_YearsIntoSalvage[years][numAlive_1_3_5_YearsIntoSalvage[years]-1] = pat->VLreal;			
			}

			if (pat->totalres >=2 && !pat->reached2DrugsResOnSalvage)
			{
				pat->reached2DrugsResOnSalvage = true;
				numReached2DrugsResOnSalvage++;
			}
		}

		if ( YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY || DEBUG_NUM_MUTS || MONITORING_PAPER)
		{
#ifdef CD4CHANGE
			if (pat->started_haart && ((pat->cyclenum - pat->haart_start_time) % CYCLES_PER_MONTH == 0) && (pat->cyclenum / CYCLES_PER_MONTH < 12 * TIME))
#else
			if ( pat->started_haart && ((pat->cyclenum - pat->haart_start_time) % CYCLES_PER_YEAR == 0) && (pat->cyclenum/CYCLES_PER_YEAR < TIME))
#endif
			{
#ifdef CD4CHANGE
				year = (pat->cyclenum - pat->haart_start_time) / CYCLES_PER_MONTH;
#else				
				year = (pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR;
#endif
				YearlyCD4[year][numAliveAtYearIntoHaart[year]] = pat->CD4real;

				//Embedding this because it also uses numAliveAtYearIntoHaart and has to come before numAliveAtYearIntoHaart.  This way, I can't mess up the order!
				if (MONITORING_PAPER)
				{
					if (pat->cyclenum - pat->haart_start_time == CYCLES_PER_YEAR * 5)
					{
						VLatYear5IntoHaart[numAliveAtYearIntoHaart[5]] = pat->VLreal;
						totalNumRegs5yearsIntoHaart_ofPatsAlive5YrsIn += min(pat->numreg, numAvailableRegimens);
						totalNumMuts5yearsIntoHaart_ofPatsAlive5YrsIn +=pat->mutCount;
					}

					if (pat->cyclenum - pat->haart_start_time == CYCLES_PER_YEAR * 10)
					{
						VLatYear10IntoHaart[numAliveAtYearIntoHaart[10]] = pat->VLreal;
						totalNumRegs10yearsIntoHaart_ofPatsAlive10YrsIn += min(pat->numreg, numAvailableRegimens);
						totalNumMuts10yearsIntoHaart_ofPatsAlive10YrsIn += pat->mutCount;
					}
				}

				numAliveAtYearIntoHaart[year]++;
			}
		}

		
		if ( MARK_ROBERTS_AHRQ )
		{
			if (pat->cyclenum<=112*7 &&
				(pat->cyclenum == 1*7 || 
				pat->cyclenum == 4*7 || 
				pat->cyclenum == 8*7 || 
				pat->cyclenum == 12*7 || 
				pat->cyclenum == 16*7 || 
				pat->cyclenum == 20*7 || 
				pat->cyclenum == 24*7 || 
				(pat->cyclenum - 24*7) % (8*7) == 0))

			{
				if (pat->cyclenum > 32*7 && pat->VLreal > 1.7) VLFail=true;
				//else if (pat->cyclenum <=32*7 && pat->HIVbaseline - pat->VLreal < 1.0) VLFail=true;
			}

			if (pat->cyclenum == TIMEPOINT*CYCLES_PER_YEAR)
			{
				tot_alive++;				
				if (pat->numreg > 1) tot_treatFail++;
				tot_regs+= pat->numreg;
				numAIDS += pat->has_AIDS;
				tot_CD4elev += pat->CD4real - pat->CD4atStartOfHaart;
				if (VLFail) numVLFail++;

				if (pat->VLreal <= 1.7) numVLsup++;	//below 1.7, the VL is considered supressed
				if (VLFail && pat->mutCount>0) numVLFailWithMuts++;
			}
		}

		if (WHO_MONITORING_STRATEGY_ANALYSIS) 
		{
			if (pat->haart && patRegCombo==0) incrementWHOArray(pat, patCyclesOn1stLineTherapy);
			else if (pat->haart && patRegCombo==1) incrementWHOArray(pat, patCyclesOn2ndLineTherapy);
			incrementWHOArray(pat, lifeCyclesLived);
			if (pat->haart && !pat->has_AIDS) incrementWHOArray(pat, lifeCyclesOnSuccessfulART);	//life years before occurrence of AIDS defining event
			else if (pat->has_AIDS) incrementWHOArray(pat, lifeCyclesWithAIDS);
		}
		
#ifdef WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
			record_discordance_data(pat);
#endif

            /************** CHECK THE STATUS OF ENGAGEMENT ********************/
#if RETENTION_IN_CARE

            Check_Engagement(pat);
			
			if (pat->VLreal < 3) total_time_VL_less1000++;	
#endif
            /************** END OF CHECKING ENGAGEMENT ************************/
            
            //Quality adjust and discount
            qamedian[patnum] +=qa_time;
            
            //gives value between 0 and 1, according to patient's CD4,that represents 1 utility-weighted month
            qa_time = Utility(pat);
            total_qa_surv_time += qa_time;
            total_disc_surv_time += 1 / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
            total_qa_disc_surv_time += qa_time / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
            

            cost = Cost(pat, TOTAL);

			individual_cost[patnum]+=cost;							//Used for Retention in Care

            total_cost += cost;
            total_disc_cost += cost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
            
            //breakdown costs into Drug, Lab, Care, and Hospitalization
            drug_cost = Cost(pat, DRUG);
            total_drug_cost += drug_cost;
            total_drug_disc_cost += drug_cost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
            
            care_cost = Cost(pat, CARE);
            total_care_cost += care_cost;
            total_care_disc_cost += care_cost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
     
			hospital_cost = Cost(pat, HOSP);
            total_hospital_cost += hospital_cost;
            total_hospital_disc_cost += hospital_cost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));

#if DETECT_VL
			alc_cost = Cost(pat, ALC_INTERV);
			total_alc_cost += alc_cost;
			total_alc_disc_cost += alc_cost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
#endif

#if RETENTION_IN_CARE
			//Note:  This is just the daily "finding/tracking" cost of somebody on outreach.  
			total_outreach_finding_cost += Cost(pat, OUTREACH);
			total_outreach_cost += Cost(pat, OUTREACH);
			
			total_RR_intervention_cost += Cost(pat, RISK_REDUCTION);
			total_SP_intervention_cost += Cost(pat, SECONDARY_PREVENTION);
#endif

		if (RYAN_WHITE && patIsTargeted) costRyanWhiteInterv += (double)costPerPersonPerYear[interv]/(double)CYCLES_PER_YEAR;

		if (pat->haart)
		{
			//obtain the compliance pattern for this next month
			//For each drug in the regimen, will the patient be compliant this month?  Yes/No
			Get_Compliance(pat);

			//get number of drugs virus is *not* susceptible to. If virus is resistant to a drug or the patient is not compliant to that drug, then it is ineffective at fighting HIV, which we make use of below
			pat->totalresnc = 0;
			pat->totalres = 0;																																										

			for (i = 0; i < DRUGS_IN_REG; i++)
			{
				pat->totalresnc += (int) Max(pat->reg[i]->res, pat->nc[i]);
				pat->totalres += pat->reg[i]->res;
			}

			//the following (fract) is a measure of the effectiveness of this regimen. It gives the fraction of drugs the virus is susceptible to.
			pat->vl_fract = (double)(DRUGS_IN_REG-pat->totalresnc)/DRUGS_IN_REG;
		}

		/*******************************************************************************************
		UPDATING PATIENT VIRAL LOAD
		-------------------------------------------------------------------------------------------
		Each patient starts with a baseline viral load.  We have this notion of an "idea" and a "real" level of VL that has to do with a lag effect of the real level approaching the ideal level.  See description below in CD4 change section and constants.h  If the patient is not on haart, the VL ideal stays at the baseline level. When the patient starts haart, the ideal drops by a VLdec amount (defined in Initialize()) subject to an adjustment based on the effectiveness of the regimen (which is based on a combination of adherence and resistance
		*******************************************************************************************/
		pat->VLideal = pat->HIVbaseline - (pat->haart*pat->vl_fract*pat->VLdec);

		pat->VLreal_before_noise = (pat->VLreal_before_noise - pat->VLideal) * exp(-rVLreal*CYCTIME_IN_MONTHS) + pat->VLideal;

		/*The VL will be constant when the patient is not on haart or when the pat is resistant to all drugs in the reg.  
		We add noise to the VL trajectory to be more realistic.*/
		if ((!pat->haart || (pat->reg[0]->res && pat->reg[1]->res && pat->reg[2]->res)) && USE_PERLIN_NOISE) pat->VLreal = addNoise(pat->VLreal_before_noise, pat, VL);
		else pat->VLreal = pat->VLreal_before_noise;

		pat->VLcat = Get_Category(pat, VL);

		if (KIMTEST && !ONLY_CORE_EVENTS) fprintf (kim_test, "\tCycle %d:  VLideal = %.2f    VLreal_after_noise = %.2f    VLcat = %d    ", pat->cyclenum, pat->VLideal, pat->VLreal, pat->VLcat);

		/*******************************************************************************************
		UPDATING PATIENT CD4 COUNT
		-------------------------------------------------------------------------------------------
		We consider CD4 count to be an aggregate of three types of CD4 dynamics: 1) CD4 changes when off of haart, 2) CD4 changes between regimens, and 3) CD4 changes within regimens.  For 2 and 3, we consider an "ideal" level and a "real" level. The ideal represents the level that the real level is approaching asymptotically (or with a lag effect).  This is why we do the real = real + (ideal-real)/adjust.  The denominator (adjust) is based on the real being X% of the ideal after Y months if everything else stayed the same.  See "constants.h" for further explanation. Note in the following that the offhaart level decreases when a patient is off haart, and freezes while a patient is on haart; the aggregated "between" effects are active from the time a patient starts haart until they die, regardless if they stop haart sometime in between; the "within effect" are also always active, but when a patient is off of haart, it is only an effect to the extent that the within effect approaches 0 asymptotically.  Note that the first change to the between level really takes effect after a patient changes regimen in Change_Regimen
		*******************************************************************************************/
		/*To clarify the following if statements, I will reference the following states:
		1)  Haart Naive (or 0% compliant)
		2)  On Haart - minimal resistance
		3)  On Haart - significant resistance (Plato equation)
		*/

		//initialize CD4_plato to an impossibly high number (for use below where we take the min)
		CD4_plato = 9999;

		// (1: Haart Naive or 0% compliant)  Decrement based on the Mellors paper:
		temp = CYCTIME_IN_MONTHS*Max(0, 1.78 + 2.8*(pat->VLreal - 3)) + DCD4NOHAART; 

		pat->CD4Mellors = pat->CD4real_before_noise - temp;

		pat->CD4Mellors = Max(0,pat->CD4Mellors);

		// (2:  On Haart, minimal resistance)
		if (pat->haart)
		{
			//sets the variable pat->CD4regressionideal
			getCD4regressionideal(pat);

			if (pat->CD4real_before_noise < pat->CD4regressionideal) pat->CD4regressionreal = (pat->CD4real_before_noise - pat->CD4regressionideal) * exp(-rCD4within_delta_real*CYCTIME_IN_MONTHS)+ pat->CD4regressionideal;

			if (pat->CD4real_before_noise >= pat->CD4regressionideal) 
			{
				delta = (pat->CD4real_before_noise - pat->CD4regressionideal) / CYCLES_PER_YEAR;

				pat->CD4regressionreal = pat->CD4real_before_noise - delta;
			}
		}


		// (3:  Pat is resistant to >= NUM_DRUGS_RES_FOR_PLATO drugs and their VL is below 4.0)
		// VL condition added on 4/20/09.  If the VL is above 4.0, the decline in CD4 is governed by PLATO.  
		// But PLATO does not impose any decline if the VL is below 4.0.  
		if (pat->haart && (pat->totalres >= NUM_DRUGS_RES_FOR_PLATO && pat->VLreal >= 4.0))
		{
			//determine the number of cycles PLATO has been applied consecutively so far
			//if (pat->exhausted_clean_regs) cycles_since_PLATO_start = pat->cyclenum - pat->haart_exhausted_clean_regs_time;
			//else cycles_since_PLATO_start = pat->cyclenum - pat->cycleMetResistanceLimit;

			//decrement based on the PLATO study:
			//On the following line, that used to be an && instead of an ||,
			//but we don't want patients CD4s to increase if the VL is less than 4.0
			//so I'm changing it to an ||.
			//Up for debate whether or not to use the two-year limit (cycles_since_PLATO_start >=CYCLES_PER_YEAR*2 && pat->VLreal < 4.0 )

			if (pat->VLreal >= 4.0) temp = ((double)-50 * pat->VLreal + 200)/CYCLES_PER_YEAR; 
			else temp = 0;


			CD4_plato = pat->CD4real_before_noise + temp;  
			CD4_plato = Max(0, CD4_plato);
		}

		//if the patient isn't on haart, the Mellors equation applies
		if (!pat->haart) pat->CD4real_before_noise = pat->CD4Mellors;

		if (pat->haart)
		{
			//If the patient has less than 50% compliance, the patient's CD4 is calculated using a combination
			//of CD4 on haart (the regression equation) and CD4 not on haart (Mellors).  At 0%, the patient
			//should be entirely using Mellors, at 50%, they should be entirely using the Regression.
			if (pat->comp < .5) pat->CD4weightedForComp = linearInterpolate(0, pat->CD4Mellors, 0.5, pat->CD4regressionreal, pat->comp);
			else pat->CD4weightedForComp = pat->CD4regressionreal;
		}

		//If the patient is on haart, we want the CD4real to be the lesser value of the CD4 
		//calculated calculated from the PLATO equation and the CD4 from the regression equation 
		//(and then waited for compliance) which may produce a lower value especially when the VL is close to 4.0.
		if (pat->haart) pat->CD4real_before_noise = Min (pat->CD4weightedForComp, CD4_plato);

		//Add noise
		if (USE_PERLIN_NOISE) pat->CD4real = Max(0, addNoise(pat->CD4real_before_noise, pat, CD4));
		else pat->CD4real = Max(0, pat->CD4real_before_noise);

		pat->CD4cat = Get_Category(pat, CD4);

		//kn  Some output to aid in debugging.  Used for graphing CD4 parameters over time.
		if (VLCD4GRAPH) graph_trajectories(pat);
		
		if (KIMTEST && !ONLY_CORE_EVENTS) fprintf (kim_test, "CD4offreal: %.1f    CD4regressionideal: %.1f    CD4real_before_noise: %.1f    CD4real: %.1f\n", pat->CD4offreal, pat->CD4regressionideal, pat->CD4real_before_noise, pat->CD4real);
		if (KIMTEST && !ONLY_CORE_EVENTS && pat->haart) fprintf (kim_test, "\t Resistances:  %d  %d  %d\n", pat->reg[0]->res, pat->reg[1]->res, pat->reg[2]->res);


		/***********************************************************************************************
		Now that we have updated the VL and CD4 levels we check if the patient dies of hiv or non-hiv-related causes by the end of this month.  Then if we determine the patient is still alive, we continue on to check to see if resistant mutations developed.
		***********************************************************************************************/

		pat->deathrateHIV = Get_HIV_Death_Rate(pat);

		//This was just a test.   We decided not to go with this code.  Leaving in for now, commented out
		//getSmoothedHIVDeathRate(pat);

#if RETENTION_IN_CARE
        check_for_AIDS_defining_events_RIC (pat);
#else
		//if the patient does not yet have AIDS, check for AIDS
		if (!pat->has_AIDS || !pat->AIDSeventReg)
		{
			//The rate of an AIDS-defining event is AIDS_RATE_MULTIPLIER * pat->deathrateHIV
			pat->pAIDS =  1 - exp(-(AIDS_EVENT_MULTIPLIER * pat->deathrateHIV) * CYCTIME);

			if ((!VRT ? uniform(&seed) : unifSim[0]) <= pat->pAIDS)
			{
				if (!pat->has_AIDS) 
				{
					pat->has_AIDS = true;
					if (DEBUG_AIDS && (int)(pat->cyclenum/CYCLES_PER_YEAR) <= 5) 
						totalAIDSOrDeath++; 

					if (WHO_MONITORING_STRATEGY_ANALYSIS) incrementWHOArray(pat, numEverDevelopedAIDS);

				}

				//Patient has had an AIDS-defining event during the current regimen
				pat->AIDSeventReg = true;
				if (KIMTEST) fprintf (kim_test, "\nPatient had AIDS-defining event at cycle: %d\n\n", pat->cyclenum);
			}
		}
#endif

		//Adjust the deathrateHIV from the table for AIDS or non-AIDS
		if (pat->has_AIDS) pat->deathrateHIV = MORTMULTHIV * AIDS_ADJUST * COINFECTION_MULT *pat->deathrateHIV;

		else pat->deathrateHIV = MORTMULTHIV * NON_AIDS_ADJUST * COINFECTION_MULT *pat->deathrateHIV;

		pat->pDieHIV = 1 - exp(-pat->deathrateHIV*CYCTIME);

		//TEST NEW
				//TEST NEW
		/*if (pat->cyclenum == DIE_AT_DAY) pat->pDieHIV = 1.0;
		else pat->pDieHIV = 0;
		*/

		if ((!VRT ? uniform(&seed) : unifSim[1]) <= pat->pDieHIV)	
		{
			Process_Death(HIV, pat);		//patient dies of HIV
			goto end;
		}


		//determine if patient will die from non-hiv-related causes
#ifdef COMPETING_RISKS_1
		for(int ic=0; ic<NUM_MORT_CONDITIONS+1; ic++)
		{
			asr_male_rate = Get_Rate(male_age_table[ic], num_male_entries, pat->age);
			asr_female_rate = Get_Rate(female_age_table[ic], num_female_entries, pat->age);

			pat->deathrate_nonHIV = MALE*asr_male_rate + FEMALE*asr_female_rate;
			//modify deathrate for haart toxicity
			if (pat->haart) pat->deathrate_nonHIV *= haart_tox;

			//there may be increased toxicity over the first three months (.25 years)of being on haart
			//This is currently set to 10.0 for resource poor, per Constantin, and 1.0 fo resource rich (no effect)
			if (pat->haart && (pat->cyclenum - pat->haart_start_time) <= TIME_HIGH_MORT*(double)CYCLES_PER_YEAR) pat->deathrate_nonHIV *= TOX_MULT_1ST_3_MONTHS;

			pat->pDieAge = 1 - exp(-pat->deathrate_nonHIV * CYCTIME);

			if ((!VRT ? uniform(&seed) : unifSim[2]) <= pat->pDieAge)
			{
				num_death_condition[ic]++;
				Process_Death(AGE, pat);		//patient dies of age related causes
				goto end;
			}
		}
#else
		asr_male_rate = Get_Rate(male_age_table, num_male_entries, pat->age);
		asr_female_rate = Get_Rate(female_age_table, num_female_entries, pat->age);

		//We don't specify sex of patient, so adjust deathrate for proportion male/female in the cohort
		pat->deathrate_nonHIV = MALE*asr_male_rate + FEMALE*asr_female_rate;

#ifdef MONITORING_PAPER_SENSITIVITY
		if (floor((double)runnum/2) == 0) pat->deathrate_nonHIV *= sens[0][runnum % 2];
#endif
		
		//If calibrating resource_rich, we want to use the VACS hiv-negative mortality value of 0.01526 
		//(veterans have higher mortality than population average)
		if (RESOURCE_RICH && CALIBRATE)
		{
			//use the VACS hiv-neg mortality rate because it's higher
			pat->deathrate_nonHIV = max (VACS_HIV_NEG_MORT_RATE, pat->deathrate_nonHIV);
			pat->deathrate_nonHIV *= VACS_HIV_NEG_MORT_MULT;
		}

		//Turn up the HIV-neg mortality for resource poor
		if (!RESOURCE_RICH && CALIBRATE)
		{
			pat->deathrate_nonHIV *= RES_POOR_HIV_NEG_MORT_MULT;
		}

		//modify deathrate for haart toxicity
		if (pat->haart) pat->deathrate_nonHIV *= haart_tox;

		//there may be increased toxicity over the first three months (.25 years)of being on haart
		//This is currently set to 10.0 for resource poor, per Constantin, and 1.0 fo resource rich (no effect)
		if (pat->haart && (pat->cyclenum - pat->haart_start_time) <= TIME_HIGH_MORT*(double)CYCLES_PER_YEAR) pat->deathrate_nonHIV *= TOX_MULT_1ST_3_MONTHS;

		//TEST NEW
		//pat->pDieAge = 0;

		pat->pDieAge = 1 - exp(-pat->deathrate_nonHIV * CYCTIME);

		//used for LOOP analysis for Scott
		//Turns off non-hiv mortality but caps life expectancy at 100 years
		if (removeNonHIVMort) 
		{
			if (pat->cyclenum == 100 * CYCLES_PER_YEAR) pat->pDieAge = 1.0;	
			else pat->pDieAge = 0;
		}

		if ((!VRT ? uniform(&seed) : unifSim[2]) <= pat->pDieAge) 
		{
			Process_Death(AGE, pat);		//patient dies of age related causes
			goto end;
		}
#endif //COMPETING_RISKS_1

		//used for resource poor calibration
		if (CALIBRATE)	
		{
			//When calibrating resource poor, determine if patient was censored
			if (pat->CD4real < CENSOR_THRESHOLD) pat->censorRate = CENSOR_RATE;
			else pat->censorRate = 0;

			pat->pCensor = 1 - exp(-pat->censorRate * CYCTIME);

			if ((!VRT ? uniform(&seed) : unifSim[13]) <= pat->pCensor)
			{
				Process_Death(CENSORED, pat);		//patient is censored
				goto end;
			}
		}


#if DETECT_VL
		if ((!pat->interv_on) && pat->haart && pat->cyclenum == pat->nextMonCycle && (pat->cyclenum - pat->haart_start_time) >= 6 * CYCLES_PER_MONTH)
		{
			if (pat->VLreal >= VL_THRESHOLD)
			{
				pat->interv_on = true; //apply intervention
#if ALC_INTERVENTION
				if (uniform(&seed) <= PROP_CURED_ALC)
				{
					pat->alc = false;
					pat->comp = NEWCOMPALC;
				}
#else
#if BOTH_AD_ALC
				pat->comp = NEWCOMPAD;
				if (uniform(&seed) <= PROP_CURED_ALC)
				{
					pat->alc = false;
					pat->comp = NEWCOMPALC;
				}
#else
				pat->comp = NEWCOMPAD;
#endif
#endif
			}
		}
#endif

		/************************************************************************************************
		patient may stop the regimen for his own reasons...
		************************************************************************************************/

		if (pat->haart)
		{
			intol_rate = stop_reg;

			//if they've previously been intolerant, they are more likely to be intolerant again
			if (pat->historyOfIntolerance) intol_rate = intol_rate * intolMultiplier;

			pat->pstopreg = 1 - exp (-intol_rate*CYCTIME);

			//if patient stops regimen, then we send him to Change_Regimen to change the HIV treatment
			//New - we only change the regimen if pat has not exhausted all regimens.  
			if ((!VRT ? uniform(&seed) : unifSim[3]) <= pat->pstopreg)
			{	
				//If the patient has not already exhausted all regimens, change to a new reg
				if (!pat->exhausted_clean_regs)
				{
					if (REASONS_FOR_REG_CHANGES)
					{
						for (int reg=0; reg<3; reg++) if (pat->numreg-1==reg) regChangeIntol[reg]++;
					}

					if (VLCD4GRAPH)	regChangeIntolFlag = true;

					if (MARK_ROBERTS_AHRQ && pat->numreg==1 && pat->cyclenum <= TIMEPOINT * CYCLES_PER_YEAR) 
						tot_intolFirstReg++;

					if (!pat->historyOfIntolerance) pat->historyOfIntolerance = true;

					if (KIMTEST) fprintf (kim_test, "\nChanging Regimen due to intolerance at cycle %d.  \n", pat->cyclenum);  

					decrementForIntolerance(pat);

					//mark that patient changed reg due to intol (didn't fail)
					lastRegFailed = false;

					Change_Regimen(pat);
					goto beginNewCycle;
				}
			}
		}	

		/*************************************************************************************
		kn new v2.0
		If the patient's VL has risen above vl_regimen_failure or the CD4 meets WHO criteria for reg failure,
		the regimen has failed and a new regimen must be chosen.  
		Note:  The VL threshold for regimen failure is subject to change after the 1st time the patient has had regimen failure
		**************************************************************************************/

		if (pat->haart && pat->cyclenum == pat->nextMonCycle || (DO_CLINICAL_VISIT_EVERY_SIX_MONTHS_FOR_WHO && pat->cyclenum == sixMonthsSinceLastVisit)) 
		{
			//maxCD4onReg is used to trigger a regimen change if running resource-poor
			if (pat->haart && pat->cyclenum == pat->nextMonCycle && pat->CD4real > pat->maxCD4onReg) pat->maxCD4onReg = pat->CD4real;

			if (DO_CLINICAL_VISIT_EVERY_SIX_MONTHS_FOR_WHO) 
			{
				incrementWHOArray(pat, numPatVisits);
				if (KIMTEST) fprintf (kim_test, "Clinic visit at cycle: %d\n", pat->cyclenum); 
				sixMonthsSinceLastVisit += CYCLES_PER_YEAR/2;
			}

			//If patient has been on regimen for at least six months (per WHO Guidelines) and we meet the reg change trigger
			//Note:  Don't check the trigger unless it's been six months because the VL threshold changes after the first failure
			if (pat->counterreg >= .5 * CYCLES_PER_YEAR)
			{
				//Note:  We use (sixMonthsSinceLastVisit - CYCLES_PER_YEAR/2.0) because we already incremented sixMonthsSinceLastVisit above.  Now we need the original value again.
				if ((pat->cyclenum == pat->nextMonCycle && metTrigger(trig, pat)) ||
					(DO_CLINICAL_VISIT_EVERY_SIX_MONTHS_FOR_WHO && pat->cyclenum == (sixMonthsSinceLastVisit-CYCLES_PER_YEAR/2.0) && metTrigger(T_CLINICAL, pat)))
				{
					if (firstCheckUnderRegDone) //if this is not the first check, check whether CD4<100;
					{
						if (pat->CD4real < 100)
							pat->lastLess100 = true;
						else
							pat->lastLess100 = false;
					}
					if (!firstCheckUnderRegDone) firstCheckUnderRegDone = true; //if this is the first check, the value of CD4 does not matter, but remember we have done the first check
					//censor patients who would have failed - used for resource poor
					if (CALIBRATE) 
					{
						pat->pCensor = PROB_CENSOR_INSTEAD_OF_FAIL;

						if ((!VRT ? uniform(&seed) : unifSim[13]) <= pat->pCensor)
						{
							Process_Death(CENSORED, pat);		//patient is censored
							goto end;
						}
					}

					if (REASONS_FOR_REG_CHANGES)
					{
						for (int reg=0; reg<3; reg++) if (pat->numreg-1==reg) regChangeFail[reg]++;
					}

					if (VLCD4GRAPH) regChangeFailFlag = true;

					//If pat wasn't resistant, they must have been intolerant
					//If they aren't already intolerant to any drugs in the reg, decrement for intolerance
					pat->totalintol=0;
					for (int i=0; i<DRUGS_IN_REG; i++) if (pat->reg[i]->intol) pat->totalintol++;
					if (pat->totalres==0 && pat->totalintol==0) decrementForIntolerance(pat);

					//mark that patient failed reg (used in Change_Regimen)
					lastRegFailed = true;

					Change_Regimen(pat);
					goto beginNewCycle;
				}
			}

			setNextMonitorCycle(pat);

		}//end section on regimen failure


		//if patient is on haart, we'll have to check for mutations and then possibly resistance
		//if not on haart, we'll just send them right to the next cycle
		if (pat->haart)
		{
			Check_For_Mutations(pat);	
		}
		
	}//end while (1)

end:;
}//end Process_Patient


void Process_Death(int type, patient *pat)
{
	double dur, durHaart = 9999;
	int year;

	if (KIMTEST && type == HIV) fprintf (kim_test, "PATIENT DIED from HIV at cycle %d.\n\n", pat->cyclenum);
	if (KIMTEST && type == AGE) fprintf (kim_test, "PATIENT DIED from Age at cycle %d.\n\n", pat->cyclenum);
	if (KIMTEST && type == CENSORED) fprintf (kim_test, "PATIENT was censored at cycle %d.\n\n", pat->cyclenum);

	if (type == HIV)
	{
		total_hiv_deaths++;
		hivdeaths[int(pat->cyclenum/CYCLES_PER_YEAR)]++;
		
#if RETENTION_IN_CARE 
	process_HIV_death_RIC(pat);
#endif
	}
	
	else if (type == AGE)
	{
		total_age_deaths++;
		agedeaths[int(pat->cyclenum/CYCLES_PER_YEAR)]++;

#ifdef COMPETING_RISKS_2
		double pDieCon[NUM_MORT_CONDITIONS];
		double asr_male_rate, asr_female_rate;
		for (int ic = 0; ic < NUM_MORT_CONDITIONS; ic++)//get the proportion of death in each listed condition
		{
			asr_male_rate = Get_Rate(male_CR_table[ic], num_male_entries, pat->age);
			asr_female_rate = Get_Rate(female_CR_table[ic], num_female_entries, pat->age);

			pDieCon[ic] = MALE*asr_male_rate + FEMALE*asr_female_rate;
		}
		
		double chance = uniform(&seed2);  //draw from 0 to 1;
		double pbefore = 0;
		double pafter = 0;
		bool bDieCon;
		bDieCon = false;
		for (int ic = 0; ic < NUM_MORT_CONDITIONS; ic++) //check in which range the random number falls
		{
			pafter += pDieCon[ic];  //upper bound
			if (chance <= pafter && chance > pbefore)
			{
				num_death_condition[ic]++;
				bDieCon = true;
				break;
			}
			pbefore = pafter;
		}
		if(!bDieCon) //not die of any condition
			num_death_condition[NUM_MORT_CONDITIONS]++;  //put to other causes

#endif

#if RETENTION_IN_CARE
     process_age_death_RIC(pat);       
#endif

	}
	else if (type == CENSORED)
	{
		total_censored++;
	}
	else
	{
		printf("type in Process_Death not valid\n");
		exit(0);
	}

#if RETENTION_IN_CARE
	if (outreach_interv_enabled && pat->outreach)  //if this pat died during the outreach process
	{
		pat->outreach = false;
		pat->deathUnkown = true;
		total_num_stopoutreach++; //outreaching is terminated

		addOutreachCostsToTotalCostAtDeath(pat);
	}
#endif

	if (DEBUG_NUM_MUTS) totalResDrugsAtDeath += countResDrugs(pat);

	//if the patient has AIDS, they were already counted
	if (DEBUG_AIDS && !pat->has_AIDS && int(pat->cyclenum/CYCLES_PER_YEAR) <= 5) 
		totalAIDSOrDeath++;

	if ((KAPLAN_MEIER_SURVIVAL_CURVE) && pat->cyclenum < (TIME*CYCLES_PER_YEAR) )
	{
		if (type == HIV)				SurvivalPlot[0][pat->cyclenum]++;
		else if (type == AGE)			SurvivalPlot[1][pat->cyclenum]++;
		else if (type == CENSORED)		SurvivalPlot[2][pat->cyclenum]++;
	}

	dur = (double) pat->cyclenum;
	if ( pat->started_haart )
		durHaart = pat->cyclenum - pat->haart_start_time; 

	total_surv_time += dur;		//gives cum total years since patient entered model

	timeTillDeath[numDeaths] = (double) dur;
	numDeaths++;


	if (MORTALITY_1_3_5_YEARS_AFTER_THERAPY )
	{

		//	find out if person died within 1, 3, 5, 10, 15 ... years
		//  Kim changed the "<=" in each if to just "<" to match the Kaplan Meier Curve data.
		if ( dur < CYCLES_PER_YEAR )
			numDeadYear[0]++;
		if ( pat->started_haart && durHaart < CYCLES_PER_YEAR )
			numDeadYearHaart[0]++;    

		if ( dur < 3.0 * CYCLES_PER_YEAR )
			numDeadYear[1]++;
		if ( pat->started_haart && durHaart < 3.0 * CYCLES_PER_YEAR )
			numDeadYearHaart[1]++;

		for ( year = 1; year <=20; year++ )
		{
			if ( dur < (year*5*CYCLES_PER_YEAR) )   
				numDeadYear[year+1]++;
			if ( pat->started_haart && durHaart < ( year*5*CYCLES_PER_YEAR ))
				numDeadYearHaart[year+1]++;
		}

		if ( dur > maxMonths )
			maxMonths = dur;
	}

	if (TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3 || REG_INFO_FOR_SENSITIVITY )
	{
		for (int i=0; i<3; i++)
		{
			if ( pat->numreg > i+1 )			//Patient fully completed reg 1, etc.
			{
				cyclesToCompleteReg[i][numPeopleCompletedReg[i]++] = pat->cyclesOnReg[i];	//for median		
				total_completed_reg_cycles[i] += pat->cyclesOnReg[i];						//for mean
			}
		}
	}

	if (KAPLAN_MEIER_CURVE_OF_TIME_TO_TREATMENT_FAILURE_REG_1_2_3 || KAPLAN_MEIER_SURVIVAL_CURVE)
	{
		int regStart=0;

		for (int i=0; pat->cyclesOnReg[i]!=0 && i<3 && i<MAX_REGS; i++)
		{
			int j;

			//Mark all of the days during the regimen as "survived"
			for (j=0; j<=pat->cyclesOnReg[i] && j < (TIME*CYCLES_PER_YEAR) ; j++)
			{
				tttfDat[j][i].survived++;
			}

			//Mark this time j at the end of the reg as failed 
			//unless that time is NOW in which case the 
			//patient died during a regimen and we censor that patient
			//We need j-1 because it got incremented at the end of the for loop above
			//And we need numreg-1 because numreg starts counting from 1 but i counts from 0
#ifdef TTTF
#ifdef WHITE_CARD
			if (pat->haart && pat->counterreg == j - 1 && pat->numreg - 1 == i) tttfDat[j][i].censored++;	//patient died during a reg, censor them
#else
			if (pat->haart && pat->counterreg == j - 1 && pat->numreg - 1 == i) tttfDat[j][i].failed++;	//patient died during a reg, in calibration, count as failure.
#endif
#else
			if (pat->haart && pat->counterreg==j-1 && pat->numreg-1==i) tttfDat[j][i].censored++;	//patient died during a reg, censor them
#endif
			else tttfDat[j][i].failed++;							//patient died later, they must have failed this reg

			regStart += pat->cyclesOnReg[i];
		}
	}

	if (KIMTEST) 
	{
		fprintf (kim_test, "\n\nReg\tCycle Reg Started\tCycles on Reg\n");
		int regStart=pat->haart_start_time;
		for (int i=0; pat->cyclesOnReg[i]!=0 && i<MAX_REGS; i++)
		{
			fprintf (kim_test, "%d\t%d\t\t%d\n", i, regStart, pat->cyclesOnReg[i]);
			regStart += pat->cyclesOnReg[i];
		}
		fprintf (kim_test, "\n\n\n");
	}

	if (VLCD4GRAPH && (CYCLES_PER_YEAR == 12 || pat->cyclenum %1 ==0 )) 
	{
		fprintf (kim_test2, "\n\nPat num: %d\n\n", patnum);
	}

	if (ACTIVE_HAART_INFO || REG_INFO_FOR_SENSITIVITY) 
	{	
		if (pat->exhausted_clean_regs)
		{
			total_time_on_haart_before_salvage += (pat->haart_exhausted_clean_regs_time - pat->haart_start_time);
			total_res_drugs_before_salvage += pat->num_res_drugs_before_salvage;
			total_num_regs_before_salvage += pat->num_regs_before_salvage;

			num_exhausted_clean_regs++;	
		}
	}

	if (VRT)
	{
		while (pat_rns_used < MAX_PAT_RNS)
			fillUnifSim();
	}

	if (GETSTATES && type == HIV)
	{

		/*Record the cost. 

		On a normal cycle, we just add the daily cost to
		the running total for the year.  However, we are trying to find
		the annual cost.  Since the patient has died, we want to extrapolate
		out the costs so far this year to cover the whole year.

		Note:  If the patient is killed on cycle 0 or is LTFU and we haven't yet calculated a cost, do so now to avoid
		a divide-by-zero error below.  And since we can't do "==0" because of floating point error, just check whether the cost
		is a very small number.  That one-day cost then has to be extrapolated out for the year, so multiply by cycles per year.
		*/

		/*NEW Record the cost.

		Normally, on the last day of the year, we record the costs accumulated for that year.
		When the patient dies, we extrapolate their costs for the year so far out for the rest of the year
		to get a cost that is representative of an entire year.
		A problem occurs if the patient dies on day 365 because the denominator in that extrapolation
		will be 0 (365 % CYCLES_PER_YEAR will be 0).  In that case, we still want to make an entry 
		because the patient dies, but we will take the cost accumulated that day and extrapolate them out for the year.

		Also note that we use (pat->cyclenum % CYCLES_PER_YEAR) + 1.0 because if we're on cycle num one, we have
		actually accumulated costs for cycle 0 and cycle 1.  (Two days.)
		*/

		if (pat->cyclenum % CYCLES_PER_YEAR == 0) costThisYearThisPat = (total_cost - totalCostAtStartOfYearThisPat) * CYCLES_PER_YEAR;	
		else costThisYearThisPat = (total_cost - totalCostAtStartOfYearThisPat) * (double) CYCLES_PER_YEAR / (double(pat->cyclenum % CYCLES_PER_YEAR) + 1.0);
		
#if RETENTION_IN_CARE 
		if (outreach_interv_enabled && pat->deathUnkown && (pat->CD4real < outreach_cd4))
		{	
			/* We add cost for outreach if the patient's death was unknown.  The cost is
			equal to the finding cost for half the normal outreach period plus
			any outreach initiation cost.
			The number of cycles of outreach added to a given year for the transmissio model lookup tables
			can't exceed the number of days in the year. */
			cycles_of_outreach_year_of_death = min((outreach_end_trigger / 2), (CYCLES_PER_YEAR - (pat->cyclenum % CYCLES_PER_YEAR)));
			costThisYearThisPat += outreach_finding_cost_poor * cycles_of_outreach_year_of_death;
			costThisYearThisPat += outreach_initiation_cost_poor;			//I didn't add a finishing cost because these patients won't be relinked.
		} 
#endif

		//The patient has died.  Use all zeros for cd4, vl, treat and res, but use 1 in new_state[4] to indicate patient has died
		//We originally had a call here to getstate(pat, cur_state);
		
		new_state[0] = 0;
		new_state[1] = 0;
		new_state[2] = 0;
		new_state[3] = 0;
		new_state[4] = 1;
		recordChange();
	}	

	if (WHO_MONITORING_STRATEGY_ANALYSIS) 
	{
		incrementWHOArray(pat, totalDeaths);
		if (type ==HIV) incrementWHOArray(pat, totalHIVDeaths);
	}
}

/**************************************************************************************************
* Initializes various patient vars
**************************************************************************************************/
void Initialize(patient *pat)
{
	int i;

	pat->cyclenum = -1;									//Because the first thing we do each cycle is increment the cyclenum.  Starts us on cycle 0.
	pat->started_haart = false;
	pat->quit_haart = 0;					//huh? (int)(MAX_MONTHS/CYCTIME)+5;
	pat->haart = 0;
	pat->exhausted_clean_regs = false;
	pat->numreg = 1;
	pat->counterreg = 0;	//number of cycles pat is on the current regimen
	pat->vl_regimen_failure = vl_regimen_failure_1st;
	pat->has_AIDS = false;
	pat->AIDSeventReg = false;	//initialize for first reg
	pat->arv_types = ARV_TYPES;
	pat->num_regs_before_salvage = 0;
	pat->num_res_drugs_before_salvage = 0;	
	pat->maxCD4onReg=0;
	pat->mutCount=0;
	pat->reached2DrugsResOnSalvage = false;

	pat->historyOfIntolerance = false;
	pat->nextMonCycle = TIME*CYCLES_PER_YEAR;				//initialize to waaaay off in the future
	pat->lastLess100 = false;

	for (int i=0; i<ARV_TYPES; i++)
	{
		switch (i)
		{
		case 0: strcpy(pat->className[i], "PI_Singular"); break;
		case 1: strcpy(pat->className[i], "PI_Boosted"); break;
		case 2: strcpy(pat->className[i], "TAM"); break;
		case 3: strcpy(pat->className[i], "NONTAM"); break;
		case 4: strcpy(pat->className[i], "Efavirenz"); break;
		case 5: strcpy(pat->className[i], "Nevirapine"); break;
		}
	}

	//patient starts out not being tolerant or resistant to any drugs
	pat->num_res[PI_SINGULAR] = 0;
	pat->num_res[PI_BOOSTED] = 0;
	pat->num_res[NRTI_TAM] = 0;
	pat->num_res[NRTI_NONTAM] = 0;
	pat->num_res[NNRTI_EFAVIRENZ] = 0;
	pat->num_res[NNRTI_NEVIRAPINE] = 0;

	//Set the number of drugs in each class
	pat->num_in_class[PI_SINGULAR] = NUM_PI_SINGULAR;
	pat->num_in_class[PI_BOOSTED] = NUM_PI_BOOSTED;
	pat->num_in_class[NRTI_TAM] = NUM_NRTI_TAM;
	pat->num_in_class[NRTI_NONTAM] = NUM_NRTI_NONTAM;
	pat->num_in_class[NNRTI_EFAVIRENZ] = NUM_NNRTI_EFAVIRENZ;
	pat->num_in_class[NNRTI_NEVIRAPINE] = NUM_NNRTI_NEVIRAPINE;

	//initialize number of mutations against each drug type	
	pat->mut_arv[PI_SINGULAR] = mut_pi_singular_start;
	pat->mut_arv[PI_BOOSTED] = mut_pi_boosted_start;
	pat->mut_arv[NRTI_TAM] = mut_nrti_tam_start;
	pat->mut_arv[NRTI_NONTAM] = mut_nrti_nontam_start;
	pat->mut_arv[NNRTI_EFAVIRENZ] = mut_nnrti_efavirenz_start;
	pat->mut_arv[NNRTI_NEVIRAPINE] = mut_nnrti_nevirapine_start;

	//Min statements are unnec for our current constants
	//Since the above comment says they are unncessary, I'm removing them!  KAN v2.0
	//gives prob of a mutation being resistant to each ARV type
	pat->pmutres[PI_SINGULAR] = adjmutres*pmutres_pi_singular_est;  //Min(.99, adjmutres*pmutres_pi_singular_est);
	pat->pmutres[PI_BOOSTED] = pmutres_pi_boosted_est;   //Min(.99, adjmutres*pmutres_pi_boosted_est);
	pat->pmutres[NRTI_TAM] = pmutres_nrti_tam_est;   //Min(.99, adjmutres*pmutres_nrti_tam_est);
	pat->pmutres[NRTI_NONTAM] = pmutres_nrti_nontam_est;   //Min(.99, adjmutres*pmutres_nrti_nontam_est);
	pat->pmutres[NNRTI_EFAVIRENZ] = pmutres_nnrti_efavirenz_est;   //Min(.99, adjmutres*pmutres_nnrti_efavirenz_est);
	pat->pmutres[NNRTI_NEVIRAPINE] = pmutres_nnrti_nevirapine_est;   //Min(.99, adjmutres*pmutres_nnrti_nevirapine_est);

	//gives prob of a resistant strain being cross resistant to other drugs in the same class
	pat->pcrossres[PI_SINGULAR] = pcrossres_pi_singular_est;   //Min(.99, adjcrossres*pcrossres_pi_singular_est);
	pat->pcrossres[PI_BOOSTED] = pcrossres_pi_boosted_est;   //Min(.99, adjcrossres*pcrossres_pi_boosted_est);
	pat->pcrossres[NRTI_TAM] = pcrossres_nrti_tam_est;   //Min(.99, adjcrossres*pcrossres_nrti_tam_est);
	pat->pcrossres[NRTI_NONTAM] = pcrossres_nrti_nontam_est;   //Min(.99, adjcrossres*pcrossres_nrti_nontam_est);
	pat->pcrossres[NNRTI_EFAVIRENZ] = pcrossres_nnrti_efavirenz_est;   //Min(.99, adjcrossres*pcrossres_nnrti_efavirenz_est);
	pat->pcrossres[NNRTI_NEVIRAPINE] = pcrossres_nnrti_nevirapine_est;   //Min(.99, adjcrossres*pcrossres_nnrti_nevirapine_est);

	pat->crossClassType[PI_SINGULAR]=PI_BOOSTED;
	pat->crossClassType[PI_BOOSTED]=PI_SINGULAR;
	pat->crossClassType[NRTI_TAM]=NRTI_NONTAM; 
	pat->crossClassType[NRTI_NONTAM]=NRTI_TAM; 
	pat->crossClassType[NNRTI_EFAVIRENZ]=NNRTI_NEVIRAPINE; 
	pat->crossClassType[NNRTI_NEVIRAPINE]=NNRTI_EFAVIRENZ; 

	pat->pCrossClassRes[PI_SINGULAR]=pcrossres_pisingular_boosted; 
	pat->pCrossClassRes[PI_BOOSTED]=pcrossres_piboosted_singular;  
	pat->pCrossClassRes[NRTI_TAM]=pcrossres_tam_nontam;  
	pat->pCrossClassRes[NRTI_NONTAM]=pcrossres_nontam_tam;  
	pat->pCrossClassRes[NNRTI_EFAVIRENZ]=pcrossres_efav_nevir;  
	pat->pCrossClassRes[NNRTI_NEVIRAPINE]=pcrossres_nevir_efav;  

	//Given the overall mutation rate, this is the probability of a mutation occurring in this class
	pat->mutRateClass[PI_SINGULAR] = mut_rate_pi;
	pat->mutRateClass[PI_BOOSTED] = mut_rate_pi;
	pat->mutRateClass[NRTI_TAM] = mut_rate_nrti;
	pat->mutRateClass[NRTI_NONTAM] = mut_rate_nrti;
	pat->mutRateClass[NNRTI_EFAVIRENZ] = mut_rate_nnrti;
	pat->mutRateClass[NNRTI_NEVIRAPINE] = mut_rate_nnrti;

	//initialize the drugs array (note that this has to be done after the num_in_class array is set above!
	for (int c=0; c<pat->arv_types; c++)
	{
		for (int d=0; d<pat->num_in_class[c]; d++)
		{
			drugs[c][d].drugNum=d;
			drugs[c][d].classNum=c;
			drugs[c][d].res=false;
			drugs[c][d].intol=false;
		}
	}

	pat->init_reg[0] = INITIAL_REG_DRUG1;
	pat->init_reg[1] = INITIAL_REG_DRUG2;
	pat->init_reg[2] = INITIAL_REG_DRUG3;

	//If we're running the calibration for RESOURCE_RICH, we want to emulate
	//the drug regimens from the VACS population used in Scott's 2007 AIDS paper
	//We set the proportion of patients on each reg according to the paper
	//Since we don't have 3NRTI Regs or "Other" in our model, we recalate the 
	//proportions according to the regs we do have
	//Prop Efav = 1140/(6397-(517+500)) = .212
	//Prop Nevir = 512/(6397-(517+500)) = .0953
	//Prop Single PI = 3324 /(6397-(517+500)) = .618
	//Prop Boosted PI = 401 / (6397-(517+500)) = .0747
	if (RESOURCE_RICH && CALIBRATE)
	{
		if (patnum+1<PATIENTS*.212) pat->init_reg[0] = NNRTI_EFAVIRENZ;
		else if (patnum+1<PATIENTS*(.212+.0953)) pat->init_reg[0] = NNRTI_NEVIRAPINE;
		else if (patnum+1<PATIENTS*(.212+.0953+.618)) pat->init_reg[0] = PI_SINGULAR;
		else if (patnum+1>=PATIENTS*(.212+.0953+.618)) pat->init_reg[0] = PI_BOOSTED;
	}

#ifdef CD4CHANGE
	if (patnum + 1<PATIENTS*.494) pat->init_reg[0] = NNRTI_NEVIRAPINE;
	else pat->init_reg[0] = NNRTI_EFAVIRENZ;
#endif

	pat->start_age = Max(Min(MAX_AGE, start_age + avgAge_SD*(!VRT ? normal(0,1,&seed) : unifInit[0])), 0); 

#if CALIBRATE
#ifdef SURVIVAL
#ifdef WHITE_CARD
	pat->CD4baseline = Max(Min(MAX_CD4, start_cd4 + avgCD4_SD*(!VRT ? normal(0, 1, &seed) : unifInit[1])), 0);
#else
	pat->CD4baseline = Max(Min(MAX_CD4, uniform_a_b(CD4_LOW, CD4_HIGH, &seed)), 0);
#endif
#endif

#ifdef TTTF
#ifdef WHITE_CARD
	pat->CD4baseline = Max(Min(MAX_CD4, start_cd4 + avgCD4_SD*(!VRT ? normal(0, 1, &seed) : unifInit[1])), 0);
#else
	double ctemp = uniform(&seed);
	if (ctemp < CD4_PROP1) pat->CD4baseline = Max(Min(MAX_CD4, uniform_a_b(0, CD4_CUT1, &seed)), 0);
	else if (ctemp >= CD4_PROP1 && ctemp < (CD4_PROP2 + CD4_PROP1)) pat->CD4baseline = Max(Min(MAX_CD4, uniform_a_b(CD4_CUT1, CD4_CUT2, &seed)), 0);
	else if (ctemp >= (CD4_PROP2 + CD4_PROP1) && ctemp < (CD4_PROP1 + CD4_PROP2 + CD4_PROP3)) pat->CD4baseline = Max(Min(MAX_CD4, uniform_a_b(CD4_CUT2, CD4_CUT3, &seed)), 0);
	else pat->CD4baseline = Max(Min(MAX_CD4, uniform_a_b(CD4_CUT3, 500, &seed)), 0);
#endif
#endif

#ifdef CD4CHANGE
	pat->CD4baseline = Max(Min(MAX_CD4, start_cd4 + avgCD4_SD*(!VRT ? normal(0, 1, &seed) : unifInit[1])), 0);
#endif
#else
	pat->CD4baseline = Max(Min(MAX_CD4, start_cd4 + avgCD4_SD*(!VRT ? normal(0, 1, &seed) : unifInit[1])), 0);
#endif

	pat->HIVbaseline = Max(Min(MAX_HIV, start_hiv + avgHIV_SD*(!VRT ? normal(0,1,&seed) : unifInit[2])), 0);	
	pat->age = pat->start_age;

	//determine whether patient starts with AIDS
	pat->has_AIDS = false;
#if CALIBRATE
#ifdef SURVIVAL
	if (uniform(&seed) < START_WITH_AIDS) pat->has_AIDS = true;
#endif
#ifdef TTTF
#ifdef WHITE_CARD
	if (uniform(&seed) < START_WITH_AIDS) pat->has_AIDS = true;
#else
	if (pat->CD4baseline < CD4_CUT1) pat->has_AIDS = true;  //all CD4 0-50 are aids.
	else if (pat->CD4baseline < CD4_CUT2 && uniform(&seed) < START_WITH_AIDS) pat->has_AIDS = true;
#endif
#endif
#else
	if (pat->CD4baseline < 200 && (!VRT ? uniform(&seed) : unifInit[3]) < AIDSrate_less200) pat->has_AIDS = true;
	if (pat->CD4baseline >= 200 && (!VRT ? uniform(&seed) : unifInit[4]) < AIDSrate_greater200) pat->has_AIDS = true;
#endif

	if (KIMTEST && pat->has_AIDS) fprintf(kim_test, "Pat started with AIDS!\n");

	//initialize compliance to drug1, drug2, drug3
	
	for (i = 0; i < DRUGS_IN_REG; i++)
	{
		pat->nc[i] = 0;
	}

	pat->VLreal_before_noise = pat->VLideal = pat->VLreal = pat->HIVbaseline;
	pat->CD4offreal = pat->CD4Peak = pat->CD4real = pat->CD4real_before_noise = pat->CD4baseline;

	//the idea behind luck1 and luck2 is that luck1 represents an individual's tendency to have good/bad CD4 reactions
	//(this is random from patient to patient). luck2 represents regimen-specific randomness to good/bad CD4 changes
	//diligence allows a random variation around the average cohort compliance from patient to patient

	//VRT
	pat->luck1 = (!VRT ? normal(0,1,&seed) : unifInit[5]); //Min(2.0, Max(-2.0, normal(0, 1, &seed)));
	pat->luck2 = (!VRT ? normal(0,1,&seed) : unifInit[6]);
	pat->diligence = (!VRT ? uniform(&seed) : unifInit[7]);

	if (MARK_ROBERTS_AHRQ)
	{
		if ((double)patnum/(double)PATIENTS<=.38) comp = .62;
		else comp = 1.0;
		VLFail = false;
	}

	pat->comp = comp;
	
#if ALC_INTERVENTION
	pat->alc = true;
#if DETECT_VL		//alc with vl, patients starts at base comp
	pat->comp = comp;
#else
	if (uniform(&seed) <= PROP_CURED_ALC)
	{
		pat->alc = false;
		pat->comp = NEWCOMPALC;
	}
#endif
#endif

#if BOTH_AD_ALC
	pat->alc = true;
#endif

#if DETECT_VL
	pat->interv_on = false; //no intervention in the beginning
#endif


	for (i=0; i<MAX_REGS; i++)
	{
		pat->cyclesOnReg[i]=0;
	}

	//The noise pattern for each patient is looking the same.  
	//Add an offset to change the start of the sin wave so pat to pat is out of phase

	offset1 = RandomNumber(999);
	offset2 = RandomNumber(999);

	if (VLCD4GRAPH)
	{
		regChangeIntolFlag = false;
		regChangeFailFlag = false;
	}

	if (DEBUG_NUM_MUTS) pat->hasResistance=false;

	//if there are mutations at the start, we want to see if any of the accumulated mutations are resistant ones
	Process_Initial_Mutations(pat);

	//Used for generating the transmission model lookup tables with GETSTATES
	totalCostAtStartOfYearThisPat = total_cost;

	//Set to nonsense number that will never be correct
	if (DO_CLINICAL_VISIT_EVERY_SIX_MONTHS_FOR_WHO) sixMonthsSinceLastVisit = -999;

	//added for Retention In Care
#if RETENTION_IN_CARE
	individual_cost[patnum]=0;
	pat->engaged_in_care = true;
	pat->engaged_in_clinic = true;
	pat->LTFU = 0;
	pat->num_LTFU = 0;
	pat->num_cycles_ES1 = 0;
	pat->num_cycles_ES3 = 0;
	pat->num_cycles_ES4 = 0;
	pat->num_cycles_ART = 0;
    pat->time_in_care = 0;
	LTFUflag = false;						
	ES1flag = ES3flag = ES4flag = 0;			
	 
	//for OUTREACH_INTERVENTION
    pat->found = true;
	pat->outreach = false;
	pat->deathUnkown = false;
	LTFUfoundflag = false;					
	startoutreachflag = false;				
	stopoutreachflag = false;		

	//Determine if patient will initiate ART when their CD4 drops below the CD4_TREAT threshold
	if (uniform(&seed) <= prob_ART_init_at_CD4_treat) pat->willInitiateHAARTAtCD4TreatThresh = true;
	else pat->willInitiateHAARTAtCD4TreatThresh = false;

#endif
}

void Check_For_Mutations(patient *pat)
{
	int j, mutType;
	double mutMultES4 = 1;		//An increase in the mutation rate for patients in state ES4.  By default, this value will be 1 to have no effect

	//Check each drug in the patient's regimen for mutations.
	//For purposes of the model, patient can only develop mutations to drugs in their regimen.

	for (j = 0; j < DRUGS_IN_REG; j++)
	{
		//for simplicity below
		mutType = pat->reg[j]->classNum;

#if RETENTION_IN_CARE
        if (pat->num_cycles_ES4 > CYCLES_PER_MONTH)  // apply multiplier on mutation if in ES4>1 month and restart ART.
        {
            mutMultES4 *= mutation_mult_es4;
        }
#endif

		//VL is a surrogate for replication.  Thus, a higher viral load results in a higher prob of mutation
		pat->adjmutrate = pat->mutRateClass[mutType]*pow(factor, pat->VLideal - 2.31) * mutMultES4;

		//now we can get prob of mutation in a month/day (depending on the cycle time)
		pat->pMutate = 1 - exp(-pat->adjmutrate*CYCTIME);

		//determine if there is a virus mutation, if not, we just go back to cycle
		//Note:  
		//1) patient can only develop a mutation to a drug if they are compliant to that drug
		//2) if the switch is set to not allow the patient to develops muts to resistant drugs
		//	 we need to make sure the patient is not resistant to that drug.

		if (!pat->nc[j] &&
			(!VRT ? uniform(&seed) : unifSim[kludge_counter++]) <= pat->pMutate && //virus mutates	//VRT
			((!CAN_GET_MUTS_TO_RES_DRUGS && !pat->reg[j]->res) || CAN_GET_MUTS_TO_RES_DRUGS)) 
		{
			//Virus mutated!
			if (KIMTEST) 
			{
				fprintf (kim_test, "\nMutation at cycle %d against drug %s.%d.    ", pat->cyclenum, pat->className[mutType], pat->reg[j]->drugNum); 
			}

			if (DEBUG_NUM_MUTS && (pat->cyclenum - pat->haart_start_time) <= CYCLES_PER_YEAR) total_muts_all_pats++;

			pat->mutCount++;  // used in MARK_ROBERTS_AHRQ and MONITORING_PAPER

			if (ANNALS_REGS_MUTS)
			{
				if ((pat->cyclenum - pat->haart_start_time) <=CYCLES_PER_YEAR*5 ) totalNumMuts1st5years++;
				if ((pat->cyclenum - pat->haart_start_time) <=CYCLES_PER_YEAR*10) totalNumMuts1st10years++;
			}

			//increment number of mutations conferring possible resistance to this arv type
			pat->mut_arv[mutType]++;

#ifdef MONITORING_PAPER_SENSITIVITY
			if(floor((double)runnum/2) == 1) pat->pmutres[mutType] *= sens[1][runnum % 2];
#endif
			
			if ((!VRT ? uniform(&seed) : unifSim[kludge_counter++]) <= pat->pmutres[mutType])		//mut res to >= 1 drugs
			{
				//assign resistance to found drug in regimen
				pat->reg[j]->res = true;

				if (KIMTEST) 
				{
					fprintf (kim_test, "Mutation is resistant.\n"); 
					fprintf(kim_test, "Res profile: (%d %d %d)  \n", pat->reg[0]->res, pat->reg[1]->res, pat->reg[2]->res);
				}

				if (DEBUG_NUM_MUTS)
				{
					if ((pat->cyclenum - pat->haart_start_time) <= CYCLES_PER_YEAR) total_res_muts_all_pats++;
					if (!pat->hasResistance)
					{
						if ((pat->cyclenum - pat->haart_start_time) <=CYCLES_PER_YEAR*1 ) TotalPatsRes1st1years++;
						if ((pat->cyclenum - pat->haart_start_time) <=CYCLES_PER_YEAR*5 ) TotalPatsRes1st5years++;
						pat->hasResistance = true;
					}
				}

				//Now check for cross-resistance and cross-class resistance to other drugs
				checkForCrossResAndCrossClassRes(pat, mutType);

			}//end if mutation is resistant 

			else if (KIMTEST) fprintf (kim_test, "Mutation not resistant.\n"); 

		}//end of if mutation
	}//end of for loop
}



/*****************************************************************************************
This module takes in a patient who needs to change drug regimens.
********************************************************************************************/
void Change_Regimen(patient *pat)
{
	int reg=999;	//initialize to nonsense number

	reg = (RESOURCE_RICH ? ChooseRegimen_Rich(pat) : ChooseRegimen_Poor(pat));

	//set next cycle to monitor CD4 or VL
	setNextMonitorCycle(pat);

	//if the patient just ran out of regimens
	if (!pat->exhausted_clean_regs && reg == MAINTAIN_CURRENT_REG ) 
	{
		pat->exhausted_clean_regs = true;
		pat->haart_exhausted_clean_regs_time = pat->cyclenum;

		if (KIMTEST) fprintf(kim_test, "\n\n***Cycle: %d:   PATIENT EXHAUSTED ALL REGIMENS!***\n\n", pat->cyclenum);

		if (REG_INFO_FOR_SENSITIVITY || ACTIVE_HAART_INFO) 
		{
			pat->num_regs_before_salvage = pat->numreg;
			pat->num_res_drugs_before_salvage = countResDrugs(pat);
		}

		if (CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE)
		{
			CD4atStartOfSalvage[numStartedSalvage] = pat->CD4real;
			VLatStartOfSalvage[numStartedSalvage] = pat->VLreal;
			numStartedSalvage++;
		}

		//If the resource rich patient just ran out of regs, but they are resistant to 2 or more drugs, call this function again
		if (RESOURCE_RICH && pat->totalres>=2) reg = ChooseRegimen_Rich(pat);

		//If the resource poor patient just ran out of clean regs, but they are resistant to 1 or more drugs, call this function again
		if (!RESOURCE_RICH && pat->totalres>=1) reg = ChooseRegimen_Poor(pat);	
	}

	if (pat->exhausted_clean_regs && reg == MAINTAIN_CURRENT_REG) 
	{
		if (KIMTEST)
		{
			fprintf(kim_test, "*Cycle: %d:   Patient remaining on current regimen.  ", pat->cyclenum);
			printDrugsInReg(pat);
			fprintf(kim_test, "Res profile: (%d %d %d)\n", pat->reg[0]->res, pat->reg[1]->res, pat->reg[2]->res);
		}
	}

	if (reg == NEW_REG)
	{
		//We increase the regimen number
		//Note, in the case that we are changing the regimen due to initial mutatiions, we don't
		//want to increase the regimen number (since the patient never really went on that reg)
		if (!((pat->cyclenum == 0 && pat->counterreg == 0) || notEnoughDrugsForInitialReg(pat))) pat->numreg++;			//increase reg number

		pat->counterreg = 0;	//reset number of cycles on current reg
		pat->cyclesOnReg[pat->numreg - 1] = 0;	//-1 since we already incremented numreg

		//Get the new VLdec according to the drugs in the new regimen
		getVLdec(pat);

		//set start-of-reg values for CD4within
		pat->AgeRegBaseline = pat->age;

		//Only reset the CD4 baseline if the last regimen had failed (not for intol)
		if (lastRegFailed) pat->CD4RegBaseline = pat->CD4real;

		firstCheckUnderRegDone = false;

		//luck2 represents regimen-specific randomness.  Since we are changing the reg, change luck2.
		pat->luck2 = (!VRT ? normal(0, 1, &seed) : unifSim[4]);

		//with a new regimen, we consider the possibility that this patient's compliance pattern will differ for this next regimen
		pat->diligence = (!VRT ? uniform(&seed) : unifSim[5]);

		//initialize the max CD4 for this new regimen to the current CD4 count
		pat->maxCD4onReg = pat->CD4real;

		// now reset comp vars
		if (comp >= .5)
			pat->comp = comp + (1 - comp)*(pat->diligence - .5);
		else
			pat->comp = comp + (0 - comp)*(pat->diligence - .5);

#if DETECT_VL
#if ALC_INTERVENTION
		if (!pat->alc)//alc cured
		{
			if (NEWCOMPALC >= .5)
				pat->comp = NEWCOMPALC + (1 - NEWCOMPALC)*(pat->diligence - .5);
			else
				pat->comp = NEWCOMPALC + (0 - NEWCOMPALC)*(pat->diligence - .5);
		}
#else
#if BOTH_AD_ALC
		if (pat->interv_on && pat->alc)
		{
			if (NEWCOMPAD >= .5)
				pat->comp = NEWCOMPAD + (1 - NEWCOMPAD)*(pat->diligence - .5);
			else
				pat->comp = NEWCOMPAD + (0 - NEWCOMPAD)*(pat->diligence - .5);
		}
		if (pat->interv_on && !pat->alc)
		{
			if (NEWCOMPALC >= .5)
				pat->comp = NEWCOMPALC + (1 - NEWCOMPALC)*(pat->diligence - .5);
			else
				pat->comp = NEWCOMPALC + (0 - NEWCOMPALC)*(pat->diligence - .5);
		}

#else
		if (pat->interv_on)
		{
			if (NEWCOMPAD >= .5)
				pat->comp = NEWCOMPAD + (1 - NEWCOMPAD)*(pat->diligence - .5);
			else
				pat->comp = NEWCOMPAD + (0 - NEWCOMPAD)*(pat->diligence - .5);
		}
#endif
#endif
#else
#if ALC_INTERVENTION
		if (!pat->alc)//alc cured
		{
			if (NEWCOMPALC >= .5)
				pat->comp = NEWCOMPALC + (1 - NEWCOMPALC)*(pat->diligence - .5);
			else
				pat->comp = NEWCOMPALC + (0 - NEWCOMPALC)*(pat->diligence - .5);
		}
#endif
#endif

		//for the clinical trigger, reset whether the patient has had an AIDS-defining
		//even yet during this regimen
		pat->AIDSeventReg = false;

		if (KIMTEST) 
		{
			fprintf (kim_test, "Changing to Reg %d at Cycle %d.\n", pat->numreg, pat->cyclenum);
			fprintf (kim_test, "New Reg:  "); 
			printDrugsInReg(pat);
			fprintf(kim_test, "Res profile: (%d %d %d)  \n", pat->reg[0]->res, pat->reg[1]->res, pat->reg[2]->res);
			fprintf (kim_test, "Round diligence = %.2f:  Round compliance = %.2f\n", pat->diligence, pat->comp);
			fprintf (kim_test, "Assigned luck2: %.2f\n", pat->luck2);
			fprintf (kim_test, "VLdec is now %.1f\n", pat->VLdec);
		}

		if ( ANNALS_REGS_MUTS )
		{
			if (pat->cyclenum <= CYCLES_PER_YEAR * 5)	totalNumRegs1st5years++;
			if (pat->cyclenum <= CYCLES_PER_YEAR * 10)  totalNumRegs1st10years++;
		}
	}

}//end Change_Regimen


/******************************************************************************************************
opens files for writing, reads in mortality tables
*******************************************************************************************************/
void Setup
()
{
#ifdef MONITORING_PAPER_SENSITIVITY
	for (int A = 0; A < 5; A++) {
		sens[A][0] = .5;
		sens[A][1] = 1.5;
		if (A == 3) {
			sens[A][0] = -100;
			sens[A][1] = 100;
		}
	}	
#endif
	
	int i=0, ival;
	char sval1[30];
	FILE *hiv, *male_india, *female_india, *male_india_CR, *female_india_CR;
	string header;

	seed = initial_seed; 
	init_genrand(seed); //sets the initial the seed

#ifdef COMPETING_RISKS_2
	seed2 = initial_seed;
	init_genrand(seed2); //sets the initial the seed
#endif

	//create the file header containing settings from the constants.h file
	headerString = createConstantsHeader();

	//open files for writing to:
	if( (output = fopen( "output.txt", "w" )) == NULL )
	{
		printf("ERROR: The output file was not opened\n" );
		exit(0);
	}
	else fprintf (output, "%s", headerString.c_str());

	if( (file1 = fopen( "file1.txt", "w" )) == NULL )
	{
		printf("ERROR: The file was not opened\n" );
		exit(0);
	}	
	else fprintf (file1, "%s", headerString.c_str());

	if( (kim_test = fopen( "kim_test.txt", "w" )) == NULL )
	{
		printf("ERROR: The testing file was not opened\n" );
		exit(0);
	}
	else fprintf (kim_test, "%s", headerString.c_str());

	if( (kim_test2 = fopen( "kim_test2.txt", "w" )) == NULL )
	{
		printf("ERROR: The testing file was not opened\n" );
		exit(0);
	}
	else fprintf (kim_test2, "%s", headerString.c_str());

	//used for WHO_MONITORING_STRATEGY_ANALYSIS
	if( (who_monitoring = fopen( "who_monitoring.txt", "w" )) == NULL )
	{
		printf("ERROR: The who_monitoring file was not opened\n" );
		exit(0);
	}
	else fprintf (who_monitoring, "%s", headerString.c_str());


#ifdef ONE_WAY_SENSITIVITY_ANALYSIS
	senseFile.open ("SensitivityOutput.txt");
	{
		if (senseFile.fail()) 
		{
			cout << "Error opening SensitivityOutput.txt" << endl;
			exit (0);
		}
		else senseFile << headerString.c_str();
	}
#endif

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS_RIC
        senseFile.open ("SensitivityOutput_RIC.txt");
        {
            if (senseFile.fail())
            {
                cout << "Error opening SensitivityOutput_RIC.txt" << endl;
                exit (0);
            }
            else senseFile << headerString.c_str();
        }
#endif

	if((stateprob = fopen( "model_probabilities.txt", "w" )) == NULL )
	{
		printf("ERROR opening model_probabilities.txt\n" );
		exit(0);
	}

	//open files for reading from

	if( (hiv = fopen( "hiv_mort_table.txt", "r" )) == NULL )
	{
		printf("ERROR: The hiv mort table file was not opened\n" );
		exit(0);
	}

	if (!RESOURCE_RICH)
	{
#ifdef COMPETING_RISKS_1
		if( (male_india = fopen( "male_CR_table_india.txt", "r" )) == NULL )
		{
			printf("ERROR: The male age table file was not opened\n" );
			exit(0);
		}

		if( (female_india = fopen( "female_CR_table_india.txt", "r" )) == NULL )
		{
			printf("ERROR: The female age table file was not opened\n" );
			exit(0);
		}
#else
		if( (male_india = fopen( "male_age_table_india.txt", "r" )) == NULL )
		{
			printf("ERROR: The male age table file was not opened\n" );
			exit(0);
		}

		if( (female_india = fopen( "female_age_table_india.txt", "r" )) == NULL )
		{
			printf("ERROR: The female age table file was not opened\n" );
			exit(0);
		}
#ifdef COMPETING_RISKS_2
		if( (male_india_CR = fopen( "male_age_CR_prop_india.txt", "r" )) == NULL )
		{
			printf("ERROR: The male age table file was not opened\n");
			exit(0);
		}

		if ((female_india_CR = fopen("female_age_CR_prop_india.txt", "r")) == NULL)
		{
			printf("ERROR: The female age table file was not opened\n");
			exit(0);
		}
#endif
#endif
	}


	//read in HIV mortality table
	i=0;
	while (1)
	{
		if (i == TABLE_SPACE)
		{
			printf("WARNING: Out of table space for HIV mortality\n");
			exit(0);
		}
		fscanf(hiv, "%d %s\n", &ival, sval1);
		hiv_mort_table[i][0] = (double) ival;
		hiv_mort_table[i][1] = atof(sval1);
		i++;
		if (feof(hiv))
			break;
	}
	num_mort_entries = i;

#ifdef COMPETING_RISKS_1
	//read in male non-hiv age-related death table
	i = 0;
	while(1)
	{
		if (i == TABLE_SPACE)
		{
			printf("WARNING: Out of table space for male mortality\n");
			exit(0);
		}
		fscanf(male_india, "%d", &ival);
		for (int iy = 0; iy < NUM_MORT_CONDITIONS+1; iy++)
		{
			male_age_table[iy][i][0] = (double)ival;
			fscanf(male_india, "%s", sval1);
			male_age_table[iy][i][1] = atof(sval1);
		}
		fscanf(male_india, "\n");
		i++;
		if (feof(male_india))
			break;
	}
	num_male_entries = i;

	//read in female non-hiv age-related death table
	i = 0;
	while(1)
	{
		if (i == TABLE_SPACE)
		{
			printf("WARNING: Out of table space for female mortality\n");
			exit(0);
		}
		fscanf(female_india, "%d", &ival);
		for (int iy = 0; iy < NUM_MORT_CONDITIONS+1; iy++)
		{
			female_age_table[iy][i][0] = (double)ival;
			fscanf(female_india, "%s", sval1);
			female_age_table[iy][i][1] = atof(sval1);
		}
		fscanf(female_india, "\n");
		i++;
		if (feof(female_india))
			break;
	}
	num_female_entries = i;
#else
	//read in male non-hiv age-related death table
	i = 0;
	while(1)
	{
		if (i == TABLE_SPACE)
		{
			printf("WARNING: Out of table space for male mortality\n");
			exit(0);
		}
		fscanf(male_india, "%d %s\n", &ival, sval1);
		male_age_table[i][0] = (double) ival;
		male_age_table[i][1] = atof(sval1);
		i++;
		if (feof(male_india))
			break;
	}
	num_male_entries = i;

	//read in female non-hiv age-related death table
	i = 0;
	while(1)
	{
		if (i == TABLE_SPACE)
		{
			printf("WARNING: Out of table space for female mortality\n");
			exit(0);
		}
		fscanf(female_india, "%d %s\n", &ival, sval1);
		female_age_table[i][0] = (double) ival;
		female_age_table[i][1] = atof(sval1);
		i++;
		if (feof(female_india))
			break;
	}
	num_female_entries = i;

#ifdef COMPETING_RISKS_2
	//read in male non-hiv age-related condition based death proportion table
	i = 0;
	while (1)
	{
		if (i == TABLE_SPACE)
		{
			printf("WARNING: Out of table space for male mortality\n");
			exit(0);
		}
		fscanf(male_india_CR, "%d", &ival);
		for (int iy = 0; iy < NUM_MORT_CONDITIONS; iy++)
		{
			male_CR_table[iy][i][0] = (double)ival;
			fscanf(male_india_CR, "%s", sval1);
			male_CR_table[iy][i][1] = atof(sval1);
		}
		fscanf(male_india_CR, "\n");
		i++;
		if (feof(male_india_CR))
			break;
	}
//	num_male_entries = i;

	//read in female non-hiv age-related condition based death proportion table
	i = 0;
	while(1)
	{
		if (i == TABLE_SPACE)
		{
			printf("WARNING: Out of table space for female mortality\n");
			exit(0);
		}
		fscanf(female_india_CR, "%d", &ival);
		for (int iy = 0; iy < NUM_MORT_CONDITIONS; iy++)
		{
			female_CR_table[iy][i][0] = (double)ival;
			fscanf(female_india_CR, "%s", sval1);
			female_CR_table[iy][i][1] = atof(sval1);
	}
		fscanf(female_india_CR, "\n");
		i++;
		if (feof(female_india_CR))
			break;
	}
//	num_female_entries = i;

#endif
#endif

	//initialize for creating lookup tables for transmission models
#if GETSTATES
	//Loop through all elements of the statecount array to initialize it
	for (int cd4_start = 0; cd4_start < STCD4; cd4_start++)
	{
		for (int vl_start = 0; vl_start < STVL; vl_start++)
		{
			for (int treat_start = 0; treat_start < STTREAT; treat_start++)
			{
				for (int res_start = 0; res_start < STRESIST; res_start++)
				{
					for (int cd4_dest = 0; cd4_dest < STCD4; cd4_dest++)
					{
						for (int vl_dest = 0; vl_dest < STVL; vl_dest++)
						{
							for (int treat_dest = 0; treat_dest < STTREAT; treat_dest++)
							{
								for (int res_dest = 0; res_dest < STRESIST; res_dest++)
								{
									for (int dead_hiv = 0; dead_hiv < STDEAD; dead_hiv++)
									{
										statecount[cd4_start][vl_start][treat_start][res_start][cd4_dest][vl_dest][treat_dest][res_dest][dead_hiv] = 0;
										totalAnnualCost[cd4_start][vl_start][treat_start][res_start][cd4_dest][vl_dest][treat_dest][res_dest][dead_hiv] = 0;
									}
								}
							}
						}
					}
				}
			}
		}
	}
#endif

}

/********************************************************************************************************
This uses a binary search to find the index such that table[index][0] = lookup and table[index][1] = the rate
we want.  lookup for the hiv mort table is really an int and for the male and female age tables it's a double,
but the code covers both types.  In the case that this is an age mort table lookup, the code ends up doing a
linear interpolation between the two table values.
***********************************************************************************************************/
double Get_Rate(double table[TABLE_SPACE][2], int entries, double lookup)
{
	int first = 0, last = entries-1, middle, found = 0, index;
	double temp, slope;

	while (first < last)
	{
		middle = (int) ceil((first + last)/2.0);
		if (fabs(table[first][0] - lookup) < EPSILON)
		{
			index = first;
			found = 1;
			break;
		}
		if (fabs(table[last][0] - lookup) < EPSILON)
		{
			index = last;
			found = 1;
			break;
		}
		if (fabs(table[middle][0] - lookup) < EPSILON)
		{
			index = middle;
			found = 1;
			break;
		}

		//if we are in here for the male or female age lookup, then we
		//need to do linear interpolation
		if (first == last - 1)  //only can get here in the case of age lookup...do linear interp
		{
			temp = table[last][0] - table[first][0];
			slope = (table[last][1] - table[first][1]) / temp;
			temp = table[first][1];
			return (temp + slope*(lookup - table[first][0]));
		}

		if (lookup < table[middle][0])
			last = middle;
		else
			first = middle;
	}

	if (found == 0)
	{
		printf("ERROR: entry not found in table\n");
		exit(0);
	}
	else
	{
		return table[index][1];
	}
}

/*******************************************************************************************************
This procedure steps through each drug the patient is taking and determines if the patient is compliant
to that drug.  Once noncompliance is observed for any drug, there is an increased chance the patient is
noncompliant to the other drugs (this accounts for temporal clustering of compliance).  We do not currently
account for compliance as a function of side effects, pill burden, etc.
*******************************************************************************************************/
void Get_Compliance(patient *pat)
{
	int i, found_nc;
	double rand[DRUGS_IN_REG];

	pat->pNocomp = 1 - pat->comp;

	found_nc = 0;  //indicates if non-compliance is observed which triggers the use of the assoc comp

	for (i = 0; i < DRUGS_IN_REG; i++)
	{
		pat->nc[i] = 0;			//initializing all to false for next loop
		rand[i] = (!VRT ? uniform(&seed) : unifSim[10+i]);

		if (rand[i] <= pat->pNocomp)	
		{
			pat->nc[i] = 1;
			found_nc= 1; 
		}
	}

	if (KIMTEST && found_nc && !ONLY_CORE_EVENTS) 
	{
		fprintf (kim_test, "Patient not compliant to: "); 
		for (i = 0; i < DRUGS_IN_REG; i++)
		{
			if (pat->nc[i]) fprintf(kim_test, "%s ", pat->className[pat->reg[i]->classNum]);
		}
		fprintf (kim_test, "\n");
	}
}

/*******************************************************************************************************
Prints output.  The infrastructure for a variety of output is coded below, though most of it may be
commented out.  Uncomment according to what you want to see printed.  At the same time, you may need
to go into the rest of the code and uncomment parts related to the calculation of the output you want.
*******************************************************************************************************/
void Process_Output()
{
	int i, j, goodcurve = TRUE, num_alive;
	double median_surv_time, median_qalys, mean_qalys;
	double medianPriorHaart = 0, medianOnHaart = 0, medianPostHaart = 0;
	double meanPriorHaart = 0, meanOnHaart = 0, meanPostHaart = 0;
	int numPatientsOnReg1 = 0, numPatientsOnReg2 = 0, numPatientsOnReg3 = 0;
	double proportion_alive;
	static bool firstPass = true;
	double mean_time_before_salvage, mean_num_regs_before_salvage, mean_res_drugs_before_salvage;

	/*  Important!  Adding Retention in Care code here because it adds costs to the 
	total_cost which gets used to calculate the mean_cost below.  */

	//Retention in Care debug output
#if RETENTION_IN_CARE
	mean_ES1_time = (total_ES1_time/(double)CYCLES_PER_YEAR)/PATIENTS;
	mean_ES3_time = (total_ES3_time/(double)CYCLES_PER_YEAR)/PATIENTS;
	mean_ES4_time = (total_ES4_time/(double)CYCLES_PER_YEAR)/PATIENTS;
	mean_time_onART = (total_time_onART/(double)CYCLES_PER_YEAR)/PATIENTS;
	mean_LTFU = (double)total_LTFU/PATIENTS;
	mean_time_VL_less1000 = (total_time_VL_less1000/(double)CYCLES_PER_YEAR)/PATIENTS;
	cout<<"mean_ES1_time="<<mean_ES1_time<<endl;
	cout<<"mean_ES3_time="<<mean_ES3_time<<endl;
	cout<<"mean_ES4_time="<<mean_ES4_time<<endl;
	cout<<"mean_time_onART="<<mean_time_onART<<endl;
	cout<<"mean_LTFU="<<mean_LTFU<<endl;

	if (outreach_interv_enabled)
	{
		cout<<"total_num_outreach="<<total_num_outreach<<endl;
		mean_outreach_cost = total_outreach_cost/PATIENTS;
		mean_num_outreach = (double)total_num_outreach/PATIENTS;
		mean_num_stopoutreach=(double)total_num_stopoutreach/PATIENTS;
		cout<<"mean_num_outreach="<<mean_num_outreach<<endl;
		cout<<"mean_num_stopoutreach="<<mean_num_stopoutreach<<endl;
		if (mean_num_outreach!=mean_num_stopoutreach)
		{
			cout<<"CAUTION: mean_num_outreach!=mean_num_stopoutreach ERROR !!!!!!"<<endl;
			system("pause");
		}
	}

	if (risk_reduction_interv_enabled)
	{
		mean_RR_intervention_cost = total_RR_intervention_cost/PATIENTS;
	}

	if (secondary_prevention_interv_enabled)
	{
		mean_SP_intervention_cost = total_SP_intervention_cost/PATIENTS;
	}
#endif

	mean_surv_time = (total_surv_time/(double)CYCLES_PER_YEAR)/PATIENTS;
	mean_disc_surv_time = (total_disc_surv_time/(double)CYCLES_PER_YEAR)/PATIENTS;
	mean_qalys = (total_qa_surv_time/(double)CYCLES_PER_YEAR)/PATIENTS;
	mean_qa_disc_surv_time = (total_qa_disc_surv_time/(double)CYCLES_PER_YEAR)/PATIENTS;

	median_surv_time = Get_Median(timeTillDeath, numDeaths )/(double)CYCLES_PER_YEAR; 
	median_qalys = Get_Median(qamedian, numDeaths)/(double)CYCLES_PER_YEAR; 

	mean_cost = total_cost/PATIENTS;
	mean_disc_cost = total_disc_cost/PATIENTS;

	mean_lab_cost = total_lab_cost/PATIENTS;
	mean_disc_lab_cost = total_lab_disc_cost/PATIENTS;

	mean_drug_cost = total_drug_cost/PATIENTS;
	mean_disc_drug_cost = total_drug_disc_cost/PATIENTS;

	mean_care_cost = total_care_cost/PATIENTS;
	mean_disc_care_cost = total_care_disc_cost/PATIENTS;

	mean_hospital_cost = total_hospital_cost/PATIENTS;
	mean_disc_hospital_cost = total_hospital_disc_cost/PATIENTS;

	mean_alc_cost = total_alc_cost / PATIENTS;
	mean_disc_alc_cost = total_alc_disc_cost / PATIENTS;

	//Dos window output
	printf("\n  median_surv_time %.2f\tmean %.2f\tmedian_qalys %.2f\tmean %.2f", median_surv_time, mean_surv_time, median_qalys, mean_qalys);


	if (ACTIVE_HAART_INFO || REG_INFO_FOR_SENSITIVITY)
	{
		//values output below
		mean_time_before_salvage = ((double)total_time_on_haart_before_salvage/(double)num_exhausted_clean_regs)/(double)CYCLES_PER_YEAR;
		mean_num_regs_before_salvage = (double)total_num_regs_before_salvage / (double)num_exhausted_clean_regs;
		mean_res_drugs_before_salvage = (double)total_res_drugs_before_salvage / (double)num_exhausted_clean_regs;
	}

	if (TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3 || REG_INFO_FOR_SENSITIVITY )
	{
		for (i=0; i<3; i++)
		{
			//initialize
			mean_reg_time[i]=0;
			median_reg_time[i]=0;

			if (numPeopleCompletedReg[i])
			{
				mean_reg_time[i] = (total_completed_reg_cycles[i]/(double)CYCLES_PER_YEAR)/numPeopleCompletedReg[i];
				median_reg_time[i] = Get_Median(cyclesToCompleteReg[i], numPeopleCompletedReg[i])/(double)CYCLES_PER_YEAR;
			}
		} 	//output below
	}

	//In the following, prob1 is Prob(survive to time ti | survived up to time t(i-1)).  
	//Prob2 is the Kaplan-Meier estimator of the survival probability at time t(i).
	if (KAPLAN_MEIER_CURVE_OF_TIME_TO_TREATMENT_FAILURE_REG_1_2_3)
	{
		double tttfDat_prob1;		
		double tttfDat_prob2;
		double tttfDat_prevProb2[3];
		int tttfDat_total;

		int dontStop = 1;

		//initialize to some nonsense value greater than 0;
		tttfDat_prevProb2[0] = tttfDat_prevProb2[1] = tttfDat_prevProb2[2] = 999999;	

		fprintf (output, "\n\nData for Kaplan Meier Curve of Time To Treatment Failure\n\n\n cycle_num\tyears\tprop_on_reg1\tprop_on_reg2\tprop_on_reg3\n");

		for (int i=0; i<TIME*CYCLES_PER_YEAR && dontStop; i++)
		{	
			fprintf(output, "%d\t %f", i, i/double(CYCLES_PER_YEAR));

			dontStop = 0;

			for (int j=0; j<3; j++)
			{
				//Stop calculating when all have reached the event
				if (tttfDat_prevProb2[j] >0 )			
				{
					tttfDat_total = tttfDat[i][j].survived + tttfDat[i][j].censored + tttfDat[i][j].failed;		

					//We want to do this until we run out of patients alive.  But we run into a problem
					//if the last living patient was censored.  If the last patient was censored, the probability
					//of being alive doesn't change beyond that time.
					if (tttfDat_total !=0) 
					{
						tttfDat_prob1 = 1 - (double)tttfDat[i][j].failed/(double)tttfDat_total;
						dontStop = 1;
					}
					else tttfDat_prob1 = 1;

					//The first doesn't get multiplied
					if (i == 0) tttfDat_prob2 = tttfDat_prob1;					
					else tttfDat_prob2 = tttfDat_prob1 * tttfDat_prevProb2[j];

					tttfDat_prevProb2[j] = tttfDat_prob2;

					fprintf(output, "\t%f", tttfDat_prob2);

				}

				else fprintf(output, "\t0", 0);
			}

			fprintf(output, "\n");
		}
	}

	printf("\n  Start age: %.2f\n  CD4 Treatment: %.2f\n  Start HIV: %.2f\n  Toxicity: %.1f\n", start_age, cd4_treat, start_hiv, haart_tox);

	if (HIV_VS_AGE_DEATHS )
	{
		fprintf (output, "HIV_VS_AGE_DEATHS\n\n");

		fprintf(output, "Year\t HIV deaths\t Age deaths\t Num Censored\n");

		for ( i = 0 ; i < TIME ; i++ )
		{		
			fprintf(output, "%d\t %d\t %d\t %d\n", i, hivdeaths[i], agedeaths[i], agedeaths[i]);
		}
	}

	if (YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY )
	{
		fprintf (output, "\n\nYEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY\n\n\nYear\tmedian CD4\tmeanCD4\n\n");

#ifdef CD4CHANGE
		for (i = 0; i < 12*TIME && numAliveAtYearIntoHaart[i] != 0; i++)
#else
		for (i = 0; i < TIME && numAliveAtYearIntoHaart[i] != 0; i++)
#endif
		{
			double medianCD4, meanCD4, totalCD4;

			medianCD4 = Get_Median(YearlyCD4[i], numAliveAtYearIntoHaart[i]);

			//Calculate mean 
			totalCD4=0;
			for (j=0; j<numAliveAtYearIntoHaart[i]; j++) totalCD4 += YearlyCD4[i][j];
			meanCD4 = totalCD4/numAliveAtYearIntoHaart[i];

			fprintf (output, "%d\t%.2f\t%.2f\n", i, medianCD4, meanCD4);
		}
	}

	if (KAPLAN_MEIER_SURVIVAL_CURVE )
	{
		fprintf (output, "\n\nData for Kaplan Meier Curve (Proportion Alive)\n\n\nCycle\tYears\thiv_deaths\tage_deaths\tnum_censored\tnum_alive\n\n");

		num_alive = PATIENTS;
		for ( j = 0 ; (j < TIME*CYCLES_PER_YEAR) && (num_alive != 0); j++)
		{
			//SurvivalPlot[0][j] are the HIV deaths, [1][j] are the Age deaths, [2][j] are Censored
			num_alive = num_alive - (SurvivalPlot[0][j] + SurvivalPlot[1][j] + SurvivalPlot[2][j]);
			proportion_alive = num_alive/(double)PATIENTS;
			fprintf (output, "%d\t %.3f\t %d\t %d\t %d\t %d\n ", j, j/double(CYCLES_PER_YEAR), SurvivalPlot[0][j], SurvivalPlot[1][j], SurvivalPlot[2][j], num_alive);
		}
	}




#ifdef ONE_WAY_SENSITIVITY_ANALYSIS 
	if (REG_INFO_FOR_SENSITIVITY) 
	{
		senseFile << label << "\t" << start_age << "\t" << start_cd4 << "\t" << cd4_treat << "\t" << start_hiv << "\t" << haart_tox << "\t" << DELTA_UTILITY_WITH_HAART << "\t" << PATIENTS << "\t" << total_hiv_deaths << "\t" << total_age_deaths << "\t" <<  median_surv_time << "\t" << mean_surv_time << "\t" << median_qalys << "\t" << mean_qalys << "\t" << mean_time_before_salvage << "\t" << mean_num_regs_before_salvage << "\t" << mean_res_drugs_before_salvage << "\t" << mean_reg_time[0] << "\t" << mean_reg_time[1] << "\t" << mean_reg_time[2] << "\t" << num_exhausted_clean_regs;

		//Add columns at the end for number of initial mutations of each drugs class
		for (int drugclass=0; drugclass<6; drugclass++)
		{
			senseFile<<"\t"<<(double)total_res_after_initial_muts[drugclass]/PATIENTS;
		}
		senseFile<<endl;
	}
	else senseFile << label << "\t" << start_age << "\t" << start_cd4 << "\t" << cd4_treat << "\t" << start_hiv << "\t" << haart_tox << "\t" << DELTA_UTILITY_WITH_HAART << "\t" << PATIENTS << "\t" << total_hiv_deaths << "\t" << total_age_deaths << "\t" <<  median_surv_time << "\t" << mean_surv_time << "\t" << median_qalys << "\t" << mean_qalys << "\t" << endl;
#endif

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS_RIC
		senseFile << label <<"\t"<< outreach_interv_enabled<< "\t"<< risk_reduction_interv_enabled<< "\t"<<secondary_prevention_interv_enabled << "\t"<<start_cd4 <<"\t"<< start_age <<"\t"<< comp <<"\t"<< cd4_treat <<"\t"<< cost_reg1_poor <<"\t"<< cost_reg2_poor <<"\t"<< cost_of_care_poor <<"\t"<< hospital_cost_poor <<"\t"<< es1_to_es2_multiplier_for_sensitivity << "\t" << pre_ART_mult_for_prob_es1_to_es2 << "\t" << prob_es2_to_es4 << "\t" << prob_es5_to_es6 << "\t"<< mutation_mult_es4 << "\t"<<outreach_prob_identify<<"\t"<<outreach_prob_find<<"\t"<<outreach_prob_relink<<"\t"<<outreach_start_trigger<<"\t"<<outreach_end_trigger<<"\t"<<outreach_initiation_cost_poor<<"\t"<<outreach_finding_cost_poor<<"\t"<<outreach_finishing_cost_poor<<"\t"<<outreach_cd4<<"\t"<<risk_reduction_rr_prob_es1_to_es2<<"\t"<<risk_reduction_intervention_cost<<"\t"<<risk_reduction_cd4<<"\t"<<secondary_prevention_rr_prob_es1_to_es2<<"\t"<<secondary_prevention_intervention_cost<<"\t"<<secondary_prevention_cd4<<"\t"<< prob_ART_init_at_CD4_treat <<"\t"<< PATIENTS <<"\t"<<total_ES1_HIV_death_known<<"\t"<< total_ES1_HIV_death_unknown << "\t" << total_ES1_AGE_death_known<<"\t"<<total_ES1_AGE_death_unknown<<"\t"<<total_ES3_HIV_death<<"\t"<<total_ES3_AGE_death<<"\t"<<total_ES4_HIV_death<<"\t"<<total_ES4_AGE_death <<"\t"<<mean_ES1_time<<"\t"<<mean_ES3_time<<"\t"<<mean_ES4_time<<"\t"<<mean_time_onART<< "\t" << mean_surv_time << "\t"<<mean_disc_surv_time<<"\t"<<median_surv_time<<"\t"<< mean_qalys<<"\t"<<mean_qa_disc_surv_time <<"\t"<<median_qalys <<"\t"<<mean_time_VL_less1000<< "\t" <<mean_cost<<"\t"<<mean_disc_cost<<"\t"<<mean_lab_cost<<"\t"<<mean_drug_cost<<"\t"<<mean_care_cost<<"\t"<<mean_hospital_cost<<"\t"<<mean_outreach_cost<<"\t"<<mean_RR_intervention_cost<<"\t"<<mean_SP_intervention_cost<<"\t"<<total_num_outreach<<"\t"<<mean_num_outreach<<"\t"<<total_LTFU<<"\t"<<mean_LTFU<<"\t"<< endl;
#endif

	//General output

	if (firstPass) 
	{
		if (SURV_RATES) fprintf(file1, "Start age\t Age SD\t Start CD4\t CD4 SD\t CD4_treat\t Start_hiv\t HIV SD\t Tox\t Util_dec\t Num_patients\t HIV Deaths\t Age Deaths\t Median Survival Time\t Mean Survival Time\t Median QALYS\t Mean QALYS\t");

		if (CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_THERAPY) 
		{
			fprintf(file1, "Median change in CD4 count <num> years into therapy\t1\t3\t5\t");
			fprintf(file1, "Median change in log viral load <num> years into therapy\t1\t3\t5\t");
		}
		
		if (ACTIVE_HAART_INFO) fprintf (file1, "prop exhausted all regs\tmean years before salvage\t mean num regs before salvage\t mean res drugs before salvage\t");

		if (CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE) 
		{
			fprintf (file1, "Proportion who developed at least two resistant drugs on salvage\t");
			fprintf(file1, "Median CD4 at start of salvage\t");
			fprintf(file1, "Median CD4 count <num> years into salvage\t1\t3\t5\t");
			fprintf(file1, "Median VL at start of salvage\t");
			fprintf(file1, "Median VL <num> years into salvage\t1\t3\t5\t");
		}

		if (ANNALS_REGS_MUTS) fprintf(file1, "Mean number of muts at 5 yrs on Haart (any mut, denominator is total patients)\t Mean number of muts at 10 yrs on Haart (any mut, denominator is total patients)\t Mean number of regs at 5 yrs on Haart (any mut, denominator is total patients)\t Mean number of regs at 10 yrs on Haart (any mut, denominator is total patients)\t");

		if (REASONS_FOR_REG_CHANGES) fprintf(file1, "Reg1 Intol\t Reg1 Fail\t Reg2 Intol\t Reg2 Fail\t Reg3 Intol\t Reg3 Fail\t");

		if (DEBUG_NUM_MUTS ) fprintf (file1, "Mean num muts in 1st year on haart\t Mean num resistant muts in 1st year on haart\t Prop with resistance 1st year on Haart\t Prop with resistance 5 years on Haart\tMean num res drugs at death\t");


		if (TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3) fprintf(file1, "Mean tttf reg1\t Mean tttf reg2\t Mean tttf reg3\t");
		if (TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3) fprintf(file1, "Median tttf reg1\t Median tttf reg2\t Median tttf reg3\t");
		if (MORTALITY_1_3_5_YEARS_AFTER_THERAPY )
		{
			fprintf(file1, "Proportion who die within:\t 1yr\t 3yrs \t5yrs\t");
			fprintf(file1, "Proportion who die within <time> after starting HAART:\t 1yr\t 3yrs \t5yrs\t");
		}

		if ( MARK_ROBERTS_AHRQ ) fprintf (file1, "numVLsup\t tot_CD4elev\t numAIDS\t num with VL failure and muts\t tot_regs\t tot_treatFail\t tot_alive\t drug class\t time point\t VL Fail Threshold\t Num intol to reg1\t Num Virologic Failure\t");

		if (DEBUG_AIDS) fprintf (file1, "Prop AIDS or dead by year 5\t");

		if (MONITORING_PAPER) 
		{
			fprintf (file1, "Num regs available\tCan return to previous regs\tTrigger\tMonitor Interval (months)\t");
			fprintf (file1, "VL threshold (if applicable)\tCost of 3rd reg (2008 USD)\t");
			fprintf (file1, "Mean life years\tMean life years (Discounted)\tMean QALYs\tMean QALYs (Discounted)\t");
			fprintf (file1, "Mean cumulative mutations at Year 5\tMean cumulative mutations at Year 10\t");
			fprintf (file1, "Median CD4 at Year 5\tMedian CD4 at Year 10\tMedian viral load at Year 5\tMedian viral load at Year 10\t");
			fprintf (file1, "Mean number of regimens used at Year 5\tMean number of regimens used at Year 10\t");
			fprintf (file1, "Mean cost (2008 USD) \tMean disc cost (2008 USD)\t");
			fprintf (file1, "Mean lab cost (2008 USD) \tMean disc lab cost (2008 USD)\t");
			fprintf (file1, "Mean drug cost (2008 USD) \tMean disc drug cost (2008 USD)\t");
			fprintf (file1, "Mean care cost (2008 USD) \tMean disc care cost (2008 USD)\t");
			fprintf (file1, "Mean hospital cost (2008 USD) \tMean disc hospital cost (2008 USD)\t");
#if DETECT_VL
			fprintf(file1, "Mean alc intervention cost (2008 USD) \tMean disc alc intervention cost (2008 USD)\t");
#endif
		}

#ifdef COMPETING_RISKS_1
		fprintf(file1, "Num death of neonatal conditions\tNum death of Lower respiratory infections\tNum death of Diarrhoeal diseases\tNum death of Unintentional injuries\t");
		fprintf(file1, "Num death of Childhood-cluster diseases\tNum death of Meningitis\tNum death of Nutritional deficiencies\tNum death of Respiratory diseases\t");
		fprintf(file1, "Num death of Neurological conditions\tNum death of Digestive diseases\tNum death of Malignant neoplasms\tNum death of Acute hepatitis B\t");
		fprintf(file1, "Num death of Intentional injuries\tNum death of Cardiovascular diseases\tNum death of Maternal conditions\tNum death of Tuberculosis\t");
		fprintf(file1, "Num death of Kidney diseases\tNum death of Diabetes mellitus\tNum other age death\t");
#endif

#ifdef COMPETING_RISKS_2
		fprintf(file1, "Num death of neonatal conditions\tNum death of Lower respiratory infections\tNum death of Diarrhoeal diseases\tNum death of Unintentional injuries\t");
		fprintf(file1, "Num death of Childhood-cluster diseases\tNum death of Meningitis\tNum death of Nutritional deficiencies\tNum death of Respiratory diseases\t");
		fprintf(file1, "Num death of Neurological conditions\tNum death of Digestive diseases\tNum death of Malignant neoplasms\tNum death of Acute hepatitis B\t");
		fprintf(file1, "Num death of Intentional injuries\tNum death of Cardiovascular diseases\tNum death of Maternal conditions\tNum death of Tuberculosis\t");
		fprintf(file1, "Num death of Kidney diseases\tNum death of Diabetes mellitus\tNum other age death\t");
#endif

		if (WHO_MONITORING_STRATEGY_ANALYSIS)
		{
			fprintf (file1, "Num patient visits at Year 5\tNum patient visits at Year 10\tNum patient visits at Year 20\t");
			fprintf (file1, "Num CD4 tests at Year 5\tNum CD4 tests at Year 10\tNum CD4 tests at Year 20\t");
			fprintf (file1, "Num VL tests at Year 5\tNum VL tests at Year 10\tNum VL tests at Year 20\t");
			fprintf (file1, "Patient-years on 1st line therapy at Year 5\tPatient-years on 1st line therapy at Year 10\tPatient-years on 1st line therapy at Year 20\t");
			fprintf (file1, "Patient-years on 2nd line therapy at Year 5\tPatient-years on 2nd line therapy at Year 10\tPatient-years on 2nd line therapy at Year 20\t");
			fprintf (file1, "Total deaths at Year 5\tTotal deaths at Year 10\tTotal deaths at Year 20\t");
			fprintf (file1, "Total HIV deaths at Year 5\tTotal HIV deaths at Year 10\tTotal HIV deaths at Year 20\t");
			fprintf (file1, "Life-years lived at Year 5\tLife-years lived at Year 10\tLife-years lived at Year 20\t");
			fprintf (file1, "Life-years lived on successful ART at Year 5\tLife-years lived on successful ART at Year 10\tLife-years lived on successful ART at Year 20\t");
			fprintf (file1, "Life-years lived with AIDS at Year 5\tLife-years lived with AIDS at Year 10\tLife-years lived with AIDS at Year 20\t");
			fprintf (file1, "Num developed AIDS at Year 5\tNum developed AIDS at Year 10\tNum developed AIDS at Year 20\t");
		}

		if (RYAN_WHITE)
		{
			fprintf (file1, "Interv Num\tNum targeted Ryan White\tAdherence RR\tLinkage RR\tCost per person per year\tNum targeted in model run\tNum baseline comp\tNum targeted who are linked\tNum not targeted who are linked\tNum cd4 treat 50\tTotal cost of intervention\tNum targeted who are suppressed (VL < 500) at 1 yr\tNum targeted with >= 1 log drop in VL at 1 yr\t");		
			fprintf (file1, "Mean interv cost at year 2\tMean interv cost at year 5\tMean interv cost at year 10\tMean interv cost at year 20\t");
			fprintf (file1, "Mean life years at year 2\tMean life years at year 5\tMean life years at year 10\tMean life years at year 20\t");	
		}

		//add an endline
		fprintf(file1, "\n\n");
		firstPass = false;
	}

	if (SURV_RATES) fprintf(file1, "%.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.3f\t %d\t %d\t %d\t %.3f\t %.3f\t %.3f\t %.3f\t", start_age, avgAge_SD, start_cd4, avgCD4_SD, cd4_treat, start_hiv, avgHIV_SD, haart_tox, DELTA_UTILITY_WITH_HAART, PATIENTS, total_hiv_deaths, total_age_deaths, median_surv_time, mean_surv_time, median_qalys, mean_qalys);

	if (CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_THERAPY) 
	{
		fprintf(file1, "\t");
		for (years=0; years<3; years++) fprintf(file1, "%.2f\t ", Get_Median(ChangeInCD4_1_3_5_YearsIntoTherapy[years], numAlive_1_3_5_YearsIntoHaart[years]));
		fprintf(file1, "\t");
		for (years=0; years<3; years++) fprintf(file1, "%.2f\t ", Get_Median(ChangeInVL_1_3_5_YearsIntoTherapy[years], numAlive_1_3_5_YearsIntoHaart[years]));
	}
	
	if (ACTIVE_HAART_INFO) fprintf (file1, "%.3f\t%.2f\t %.2f\t %.2f\t", (double)num_exhausted_clean_regs/PATIENTS, mean_time_before_salvage, mean_num_regs_before_salvage, mean_res_drugs_before_salvage);

	if (CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE) 
	{
		fprintf (file1, "%.2f\t", (double)numReached2DrugsResOnSalvage/(double)PATIENTS);

		fprintf (file1, "%.2f\t", Get_Median(CD4atStartOfSalvage, numStartedSalvage));

		fprintf(file1, "\t");
		for (years=0; years<3; years++) fprintf(file1, "%.2f\t ", Get_Median(CD4_1_3_5_YearsIntoSalvage[years], numAlive_1_3_5_YearsIntoSalvage[years]));

		fprintf (file1, "%.2f\t", Get_Median(VLatStartOfSalvage, numStartedSalvage));

		fprintf(file1, "\t");
		for (years=0; years<3; years++) fprintf(file1, "%.2f\t ", Get_Median(VL_1_3_5_YearsIntoSalvage[years], numAlive_1_3_5_YearsIntoSalvage[years]));
	}

	if (ANNALS_REGS_MUTS) fprintf(file1, "%.2f\t %.2f\t %.2f\t %.2f\t", (double)totalNumMuts1st5years/(double)PATIENTS, (double)totalNumMuts1st10years/(double)PATIENTS,(double)totalNumRegs1st5years/(double)PATIENTS, (double)totalNumRegs1st10years/(double)PATIENTS);

	if (REASONS_FOR_REG_CHANGES) 
	{
		for (int reg=0; reg<3; reg++) 
		{
			if (regChangeIntol[reg]==0 && regChangeFail[reg]==0) fprintf(file1, "0\t0\t");
			else fprintf(file1, "%.2f\t%.2f\t", (double)regChangeIntol[reg]/(double)(regChangeIntol[reg]+regChangeFail[reg]), (double)regChangeFail[reg]/(double)(regChangeIntol[reg]+regChangeFail[reg]));
		}
	}

	if (DEBUG_NUM_MUTS ) fprintf (file1, "%.3f\t %.3f\t %.3f\t %.3f\t %.3f\t", (double)total_muts_all_pats/PATIENTS, (double)total_res_muts_all_pats/PATIENTS, (double)TotalPatsRes1st1years/(double)numAliveAtYearIntoHaart[1], (double)TotalPatsRes1st5years/(double)numAliveAtYearIntoHaart[5], (double)totalResDrugsAtDeath/(double)PATIENTS);

	if (TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3) fprintf(file1, "%.2f\t %.2f\t %.2f\t", mean_reg_time[0], mean_reg_time[1], mean_reg_time[2]);
	if (TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3) fprintf(file1, "%.2f\t %.2f\t %.2f\t", median_reg_time[0], median_reg_time[1], median_reg_time[2]);

	if (MORTALITY_1_3_5_YEARS_AFTER_THERAPY)
	{
		fprintf(file1, "\t%.4f\t %.4f\t %.4f\t", (double) numDeadYear[0]/PATIENTS, (double) numDeadYear[1]/PATIENTS, (double) numDeadYear[2]/PATIENTS);
		fprintf(file1, "\t %.4f\t %.4f\t %.4f\t", (double) numDeadYearHaart[0]/PATIENTS, (double) numDeadYearHaart[1]/PATIENTS, (double) numDeadYearHaart[2]/PATIENTS);
	}

	if ( MARK_ROBERTS_AHRQ ) fprintf (file1, "%d\t %f\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %.1f\t %d\t %d\t", numVLsup, tot_CD4elev, numAIDS, numVLFailWithMuts, tot_regs, tot_treatFail, tot_alive, INITIAL_REG_DRUG1, TIMEPOINT, (float)vl_regimen_failure_1st, tot_intolFirstReg, numVLFail);

	if (DEBUG_AIDS) fprintf (file1, "%.3f\t", (double)totalAIDSOrDeath/(double)PATIENTS);

	if (MONITORING_PAPER) 
	{
		fprintf (file1, "%d\t%d\t%d\t%d\t", numAvailableRegimens, poorCanGoBackOnPreviousRegs, trig, monitor_interval);
		fprintf (file1, "%.1f\t%.2f\t",  vl_regimen_failure_1st, cost_reg3_poor); 
		fprintf (file1, "%.3f\t%.3f\t%.3f\t%.3f\t", mean_surv_time, mean_disc_surv_time, mean_qalys, mean_qa_disc_surv_time);
		fprintf (file1, "%.2f\t%.2f\t", (double)totalNumMuts5yearsIntoHaart_ofPatsAlive5YrsIn/(double)numAliveAtYearIntoHaart[5],  (double)totalNumMuts10yearsIntoHaart_ofPatsAlive10YrsIn/(double)numAliveAtYearIntoHaart[10]);
		fprintf (file1, "%.2f\t%.2f\t%.2f\t%.2f\t", Get_Median(YearlyCD4[5], numAliveAtYearIntoHaart[5]), Get_Median(YearlyCD4[10], numAliveAtYearIntoHaart[10]), Get_Median(VLatYear5IntoHaart,numAliveAtYearIntoHaart[5]),  Get_Median(VLatYear10IntoHaart,numAliveAtYearIntoHaart[10]));
		fprintf (file1, "%.2f\t%.2f\t", (double)totalNumRegs5yearsIntoHaart_ofPatsAlive5YrsIn/(double)numAliveAtYearIntoHaart[5], (double)totalNumRegs10yearsIntoHaart_ofPatsAlive10YrsIn/(double)numAliveAtYearIntoHaart[10]);
		fprintf (file1, "%.2f\t%.2f\t", mean_cost, mean_disc_cost);
		fprintf (file1, "%.2f\t%.2f\t", mean_lab_cost, mean_disc_lab_cost);
		fprintf (file1, "%.2f\t%.2f\t", mean_drug_cost, mean_disc_drug_cost);
		fprintf (file1, "%.2f\t%.2f\t", mean_care_cost, mean_disc_care_cost);
		fprintf (file1, "%.2f\t%.2f\t", mean_hospital_cost, mean_disc_hospital_cost);
#if DETECT_VL
		fprintf(file1, "%.2f\t%.2f\t", mean_alc_cost, mean_disc_alc_cost);
#endif
	}

#ifdef COMPETING_RISKS_1
	for (int ix = 0; ix < NUM_MORT_CONDITIONS + 1; ix++)
		fprintf(file1, "%d\t", num_death_condition[ix]);
#endif

#ifdef COMPETING_RISKS_2
	for (int ix = 0; ix < NUM_MORT_CONDITIONS + 1; ix++)
		fprintf(file1, "%d\t", num_death_condition[ix]);
#endif

	if (WHO_MONITORING_STRATEGY_ANALYSIS)
	{	
		fprintf (file1, "%d\t%d\t%d\t", numPatVisits[WHO_5_YR_INDEX], numPatVisits[WHO_10_YR_INDEX], numPatVisits[WHO_20_YR_INDEX]);
		fprintf (file1, "%d\t%d\t%d\t", numCD4Tests[WHO_5_YR_INDEX], numCD4Tests[WHO_10_YR_INDEX], numCD4Tests[WHO_20_YR_INDEX]);
		fprintf (file1, "%d\t%d\t%d\t", numVLTests[WHO_5_YR_INDEX], numVLTests[WHO_10_YR_INDEX], numVLTests[WHO_20_YR_INDEX]);
		fprintf (file1, "%.3f\t%.3f\t%.3f\t", (double)patCyclesOn1stLineTherapy[WHO_5_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)patCyclesOn1stLineTherapy[WHO_10_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)patCyclesOn1stLineTherapy[WHO_20_YR_INDEX]/(double)CYCLES_PER_YEAR);
		fprintf (file1, "%.3f\t%.3f\t%.3f\t", (double)patCyclesOn2ndLineTherapy[WHO_5_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)patCyclesOn2ndLineTherapy[WHO_10_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)patCyclesOn2ndLineTherapy[WHO_20_YR_INDEX]/(double)CYCLES_PER_YEAR);
		fprintf (file1, "%d\t%d\t%d\t", totalDeaths[WHO_5_YR_INDEX], totalDeaths[WHO_10_YR_INDEX], totalDeaths[WHO_20_YR_INDEX]);
		fprintf (file1, "%d\t%d\t%d\t", totalHIVDeaths[WHO_5_YR_INDEX], totalHIVDeaths[WHO_10_YR_INDEX], totalHIVDeaths[WHO_20_YR_INDEX]);
		fprintf (file1, "%.3f\t%.3f\t%.3f\t", (double)lifeCyclesLived[WHO_5_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)lifeCyclesLived[WHO_10_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)lifeCyclesLived[WHO_20_YR_INDEX]/(double)CYCLES_PER_YEAR);
		fprintf (file1, "%.3f\t%.3f\t%.3f\t", (double)lifeCyclesOnSuccessfulART[WHO_5_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)lifeCyclesOnSuccessfulART[WHO_10_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)lifeCyclesOnSuccessfulART[WHO_20_YR_INDEX]/(double)CYCLES_PER_YEAR);
		fprintf (file1, "%.3f\t%.3f\t%.3f\t", (double)lifeCyclesWithAIDS[WHO_5_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)lifeCyclesWithAIDS[WHO_10_YR_INDEX]/(double)CYCLES_PER_YEAR, (double)lifeCyclesWithAIDS[WHO_20_YR_INDEX]/(double)CYCLES_PER_YEAR);
		fprintf (file1, "%d\t%d\t%d\t", numEverDevelopedAIDS[WHO_5_YR_INDEX], numEverDevelopedAIDS[WHO_10_YR_INDEX], numEverDevelopedAIDS[WHO_20_YR_INDEX]);
		
		WHO_print_year_by_year_stats();
	}

	if (RYAN_WHITE)
	{
		fprintf (file1, "%d\t%d\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%li\t%d\t%d\t", interv, numTargetedInterv[interv], adh_RR[interv], linkage_RR[interv], costPerPersonPerYear[interv], numTargeted, numBaselineComp, numTargetedLinked, numNotTargetedLinked, totalNumCD4_50, (long int)costRyanWhiteInterv, numSuppressed, numDroppedAtLeastOneLogVL);
		fprintf (file1, "%.2f\t%.2f\t%.2f\t%.2f\t", costAtYear[0]/(double)PATIENTS,costAtYear[1]/(double)PATIENTS,costAtYear[2]/(double)PATIENTS,costAtYear[3]/(double)PATIENTS);
		fprintf (file1, "%.2f\t%.2f\t%.2f\t%.2f\t", cyclesLivedAtYear[0]/((double)CYCLES_PER_YEAR*PATIENTS),cyclesLivedAtYear[1]/((double)CYCLES_PER_YEAR*PATIENTS),cyclesLivedAtYear[2]/((double)CYCLES_PER_YEAR*PATIENTS),cyclesLivedAtYear[3]/((double)CYCLES_PER_YEAR*PATIENTS));
	}

#ifdef WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
	process_WHO_discordance_output();
#endif

	//add an endline
	fprintf(file1, "\n");

	//flush the buffers so all data is written to output file or window
	fflush (NULL);
	fflush (stdout);


}

void Process_Initial_Mutations(patient *pat)
{
	int mutType, numMuts, drugNum;
	double pmutres;

	for (mutType=0; mutType<pat->arv_types; mutType++)
	{
		numMuts = pat->mut_arv[mutType];

		//Find the probability that at least one drug in the class is resistant
		pmutres = 1 - pow((1-pat->pmutres[mutType]), numMuts);

		if ((!VRT ? uniform(&seed) : unifSim[kludge_counter++]) <= pmutres)		
		{
			if (KIMTEST) fprintf (kim_test, "Initial mutation to class %s is resistant.  ", pat->className[mutType]);
			//Make sure there are drugs in the class, if so, pick one at random to mark as resistant
			if (pat->num_in_class[mutType] >0) 
			{
				drugNum = RandomNumber(pat->num_in_class[mutType]);
				drugs[mutType][drugNum].res = true;

				if (KIMTEST) fprintf (kim_test, "Choosing drug in class at random to mark as resistant:  %s.%d\n", pat->className[mutType], drugNum ); 
			}

			//Check for cross-resistance and cross-class resistance.
			checkForCrossResAndCrossClassRes(pat, mutType);
		}
	}

	if (REG_INFO_FOR_SENSITIVITY)
	{
		for (int drugclass=0; drugclass<ARV_TYPES; drugclass++)
		{
			total_res_after_initial_muts[drugclass] += pat->num_res[drugclass];
		}
	}
}

/*******************************************************************************************************
Initializes many counter variables to be 0.
*******************************************************************************************************/
void Initialize_Vars()
{ 
	int i, j;

	firstCheckUnderRegDone = false;

	numDeaths = 0;
	total_hiv_deaths = 0;
	total_age_deaths = 0;
	total_censored = 0;
	total_surv_time = 0;

#ifdef COMPETING_RISKS_1
	for (i = 0; i < NUM_MORT_CONDITIONS + 1; i++)
		num_death_condition[i] = 0;
#endif

#ifdef COMPETING_RISKS_2
	for (i = 0; i < NUM_MORT_CONDITIONS + 1; i++)
		num_death_condition[i] = 0;
#endif

	if (ANNALS_REGS_MUTS || MONITORING_PAPER)
	{
		totalNumMuts1st5years=totalNumMuts1st10years=0;
		totalNumRegs1st5years=totalNumRegs1st10years=PATIENTS;
	}

	if (MONITORING_PAPER) 
	{
		totalNumRegs5yearsIntoHaart_ofPatsAlive5YrsIn=0;
		totalNumRegs10yearsIntoHaart_ofPatsAlive10YrsIn=0;
		totalNumMuts5yearsIntoHaart_ofPatsAlive5YrsIn=0;
		totalNumMuts10yearsIntoHaart_ofPatsAlive10YrsIn=0;
	}

	total_qa_surv_time = total_qa_disc_surv_time = total_cost = total_disc_cost = total_disc_surv_time = 0;
	total_lab_cost = mean_lab_cost = total_lab_disc_cost = mean_disc_lab_cost = 0;
	care_cost = total_care_cost = mean_care_cost = total_care_disc_cost = mean_disc_care_cost = 0;
	drug_cost = total_drug_cost = mean_drug_cost = total_drug_disc_cost = mean_disc_drug_cost = 0;
	hospital_cost = total_hospital_cost = mean_hospital_cost = total_hospital_disc_cost = mean_disc_hospital_cost = 0;
	alc_cost = total_alc_cost = mean_alc_cost = total_alc_disc_cost = mean_disc_alc_cost = 0;

	num_exhausted_clean_regs = 0;

	if (ACTIVE_HAART_INFO || REG_INFO_FOR_SENSITIVITY) 
	{
		total_time_on_haart_before_salvage=0, total_res_drugs_before_salvage=0, total_num_regs_before_salvage=0;

		for (i=0; i<ARV_TYPES; i++)
		{
			total_res_after_initial_muts[i]=0;
		}
	}

	if (DEBUG_NUM_MUTS) 
	{
		total_muts_all_pats=0, total_res_muts_all_pats=0;
		TotalPatsRes1st1years = 0;
		TotalPatsRes1st5years = 0;
		totalResDrugsAtDeath = 0;
	}

	for(i=0; i<PATIENTS; i++)
	{
		qamedian[i] = 0;
	}

	for (j = 0; j < 22; j++ )
	{
		numDeadYear[j] = 0;
		numDeadYearHaart[j] = 0;
	}

	if (KAPLAN_MEIER_SURVIVAL_CURVE)
	{
		for (j=0; j<TIME; j++)
		{
			SurvivalPlot[0][j]=0;
			SurvivalPlot[1][j]=0;
			SurvivalPlot[2][j]=0;
		}
	}

	if (KAPLAN_MEIER_CURVE_OF_TIME_TO_TREATMENT_FAILURE_REG_1_2_3 || 
		TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3 ||
		REG_INFO_FOR_SENSITIVITY)
	{
		numPeopleCompletedReg[0] = numPeopleCompletedReg[1] = numPeopleCompletedReg[2] = 0;  //Kim added
		total_completed_reg_cycles[0] = total_completed_reg_cycles[1] = total_completed_reg_cycles[2] = 0;	//Kim
	}

	if (CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_THERAPY)
	{
		for (years=0; years<3; years++) numAlive_1_3_5_YearsIntoHaart[years] = 0;
	}

	if (CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE)
	{
		for (years=0; years<3; years++) numAlive_1_3_5_YearsIntoSalvage[years] = 0;
	}

	maxMonths = 0;

	if (HIV_VS_AGE_DEATHS )
	{
		for ( i = 0 ; i < TIME ; i++ )
		{
			hivdeaths[i] = 0;
			agedeaths[i] = 0;
		}
	}

	//determine the probability of developing a mutation to each drug type
	//note:  mut_rate_nrti + mut_rate_nnrti + mut_rate_pi = mut_rate;
	//note:  mut_rate_nrti + RATIO_NNRTI_TO_NRTI*mut_rate_nrti + RATIO_PI_TO_NRTI*pmut_nrti = mut_rate;
	mut_rate_nrti = mut_rate/(1.0 + RATIO_NNRTI_TO_NRTI + RATIO_PI_TO_NRTI);
	mut_rate_nnrti = RATIO_NNRTI_TO_NRTI * mut_rate_nrti;
	mut_rate_pi = RATIO_PI_TO_NRTI * mut_rate_nrti;

	if (MARK_ROBERTS_AHRQ)
	{
		//Good, but testing to get proportion supressed down
		mut_rate_nrti = .099*1.7;   //.091;
		mut_rate_nnrti = .69*1.7;
		mut_rate_pi = .099*1.7;   //.091;

	}


	getInitialAIDSRates();

	//calculate the exponential rates r used in the VLreal and CD4real calculations
	rVLreal = -log(1-1/VL_ADJUST);
	rCD4within_delta_real = -log(1-1/CD4_ADJUST_W);

	if (VRT)
	{
		if (CYCLE == MONTH) MAX_PAT_RNS = MAX_SEED_INIT + MAX_SEED_CYCLE*MAX_MONTHS;
		else if (CYCLE == DAY) MAX_PAT_RNS = MAX_SEED_INIT + MAX_SEED_CYCLE*((MAX_MONTHS-2)/12)*365;
	}

	if (KAPLAN_MEIER_CURVE_OF_TIME_TO_TREATMENT_FAILURE_REG_1_2_3 || KAPLAN_MEIER_SURVIVAL_CURVE)
	{
		for (int i = 0; i< TIME*CYCLES_PER_YEAR; i++)
		{
			for (int j = 0; j< 3; j++)
			{
				tttfDat[i][j].censored = 0;
				tttfDat[i][j].failed = 0;
				tttfDat[i][j].survived = 0;
			}
		}	
	}

	if ( YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY || DEBUG_NUM_MUTS || MONITORING_PAPER)
	{
		for (int i = 0; i < TIME; i++)
		{
			numAliveAtYearIntoHaart[i] = 0;
		}
	}

	if (REASONS_FOR_REG_CHANGES)
	{
		for (int reg=0; reg<3; reg++) 
		{
			regChangeIntol[reg]=0;
			regChangeFail[reg]=0;
		}
	}

	if ( MARK_ROBERTS_AHRQ )
	{
		numVLsup=0;
		tot_CD4elev=0;
		numAIDS=0;
		numVLFailWithMuts=0;
		tot_regs=0;
		tot_treatFail=0;
		tot_alive=0;
		tot_intolFirstReg=0;
		numVLFail=0;
	}

	if (DEBUG_AIDS) totalAIDSOrDeath = 0;

	numStartedSalvage = 0;
	numReached2DrugsResOnSalvage = 0;

	//used for Ryan White Analysis
	if (RYAN_WHITE) 
	{
		costRyanWhiteInterv = 0;

		//For verification of my algorithm
		numBaselineComp = 0;
		totalNumCD4_50 = 0;
		numTargetedLinked = 0;
		numNotTargetedLinked = 0;
		numSuppressed = 0;
		numDroppedAtLeastOneLogVL = 0;

		//initialize for the 4 time periods we are counting the following variables
		for (int i=0; i<NUM_TIME_PERIODS_RYAN_WHITE; i++)
		{
			costAtYear[i] = 0;		
			cyclesLivedAtYear[i] = 0;
		}
	}

#ifdef WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
	// Set counts for all categories to 0:
	for (int aids = 0; aids < 2; aids++) {
		for (int yrs = 0; yrs < 5; yrs++) {
			cumulative_years[aids][yrs].initialize();
		}
	}
#endif
        
#if RETENTION_IN_CARE
	total_LTFU=0;
	mean_LTFU =0;
	total_ES1_HIV_death_known = 0;
	total_ES1_HIV_death_unknown = 0;
	total_ES1_AGE_death_known = 0;
	total_ES1_AGE_death_unknown = 0;
	total_ES3_HIV_death = 0;
	total_ES3_AGE_death = 0;
	total_ES4_HIV_death = 0;
	total_ES4_AGE_death = 0;
	total_ES1_time = 0;
	total_ES3_time = 0;
	total_ES4_time = 0;
	total_time_onART = 0;
	total_time_VL_less1000 = 0;
	mean_ES1_time = 0;
	mean_ES3_time = 0;
	mean_ES4_time = 0;
	mean_time_onART = 0;
	mean_time_VL_less1000 = 0;

	//for OUTREACH_INTERVENTION
	total_outreach_cost = total_outreach_finding_cost = 0;
	total_num_outreach = 0;
	total_num_stopoutreach = 0;

	//for RISK_REDUCTION_INTERVENTION
	total_RR_intervention_cost=mean_RR_intervention_cost=0;

	//for SECONDARY_PREVENTION_INTERVENTION
	total_SP_intervention_cost=mean_SP_intervention_cost=0;

#endif
	
}
/*******************************************************************************************************
Obtains either the viral load category or the CD4 category  according to the "type" input parameter
*******************************************************************************************************/
int Get_Category(patient *pat, int type)
{
	if (type == VL)
	{
		if (pat->VLreal < 3.5)
			return 1;
		else if (pat->VLreal < 4.5)
			return 2;
		else if (pat->VLreal < 5.5)
			return 3;
		else
			return 4;
	}
	else //type == CD4
	{
		if (pat->CD4real < 50)
			return 1;
		else if (pat->CD4real < 200)
			return 2;
		else if (pat->CD4real < 350)
			return 3;
		else if (pat->CD4real < 500)
			return 4;
		else
			return 5;
	}
}

/*******************************************************************************************************
Swaps the (x+1)th and the (y+1)th element of the array "a".  Used in Quicksort algorithm.
*******************************************************************************************************/
void Swap( double a[], int x, int y )
{	
	double t = a[x]; 
	a[x] = a[y]; 
	a[y] = t; 
}

/*******************************************************************************************************
Sorts array a from highest to lowest value
l and r are the left and right array indices
*******************************************************************************************************/
void QuickSort(double a[], int l, int r )
{
	int i, j, /*v, *t,*/ind;

	double v, t;	//Kim changed v and t to doubles since our array is all doubles.

	if ( r > l )
	{
		ind = ((l+r)/2);

		v = a[ind]; i = l-1; j = r;

		//Swap( a, ind, r ); 
		t=a[ind];
		a[ind]=a[r];
		a[r]=t;

		for ( ;; )
		{
			while ( i < r && a[++i] < v ) {} 
			while ( j > l && a[--j] > v ) {}

			if ( i >= j ) break;
			//Swap( a, i, j );
			t=a[i];
			a[i]=a[j];
			a[j]=t;
		}

		//Swap( a, i, r );
		t=a[i];
		a[i]=a[r];
		a[r]=t;
		QuickSort( a, l, i-1 );
		QuickSort(a, i+1, r );
	}

	return;
}

/*******************************************************************************************************
Gets the median value of an array of values, a[], of size "size".  Calls Quicksort to sort the array
*****************************************************************************************************/
double Get_Median(double a[], int size)
{
	double median;
	int index1, index2;

	QuickSort(a, 0, size-1);

	//if size is even, average middle two values, if odd, take middle value
	if (size%2 == 0)		//size is even
	{
		index1 = size/2 - 1;
		index2 = size/2 ;
		//take average
		median = (a[index1] + a[index2])/2.0;
	}
	else
	{
		index1 = ( (size + 1)/2 ) - 1;
		median = a[index1];
	}

	return median;
}

/*******************************************************************************************************
Returns a value in [0,1] inidicating a utility weight that is based on various patient parameters.
Based on values from a Freedberg paper.
********************************************************************************************************/
double Utility(patient *pat)
{
	double retval = 0;

	if (pat->CD4real < 100)
		retval += UTILITY_CD4_LESS_THAN_100;
	else if (pat->CD4real < 200)
		retval += UTILITY_CD4_LESS_THAN_200;
	else 
		retval += UTILITY_CD4_ABOVE_200;

	if (pat->haart) retval += DELTA_UTILITY_WITH_HAART;

	return retval;
}

/*******************************************************************************************************
Returns a cost for the month.  Based on Sanders et al (NEJM 2005).  There they estimated the following annual
costs (which I convert to monthly values in the Cost() function):
HIV infection, CD4 > 500: $2,978
HIV infection, CD4 between 200 and 500: $5,096
HIV infection, CD4 < 200: $7,596
Cost of AIDS (irrespective of CD4 count): $10,998 (we are not implementing this yet)
3-drug therapy: $13,752
********************************************************************************************************/

double Cost(patient *pat, int whichCost)
{
	double retval = 0;

	if (RESOURCE_RICH)
	{

	if (pat->CD4real < 200)
		retval += 7596.0/(double)CYCLES_PER_YEAR;
	else if (pat->CD4real < 500)
		retval += 5096.0/(double)CYCLES_PER_YEAR;
	else
		retval += 2978.0/(double)CYCLES_PER_YEAR;

	if (pat->haart)
		retval += 13752.0/(double)CYCLES_PER_YEAR; 
	}

	else	//RESOURCE_POOR
	{
		

#if RETENTION_IN_CARE
			if(pat->engaged_in_care)
			{
#endif
				if (whichCost == TOTAL || whichCost == CARE)
				{
					//All patients incur a cost of care, regardless of whether they are on haart
					retval += cost_of_care_poor/(double)CYCLES_PER_YEAR;
				}

				if (whichCost == TOTAL || whichCost == HOSP)
				{
					//If patient has had an AIDS defining event and their CD4 drops below 200,
					//apply hospital costs.  
			//		if (pat->has_AIDS && pat->CD4real < 200) 
			//			retval += hospital_cost_poor/(double)CYCLES_PER_YEAR;
					//NEW WAY of calculating hospital cost. CD4 stratified. 04/30/2015 QZ.
					if(pat->CD4real < 200)
						retval += HOSPITAL_COST_POOR_LOW_CD4/(double)CYCLES_PER_YEAR;
					if(pat->CD4real >= 200 && pat->CD4real < 350)
						retval += HOSPITAL_COST_POOR_MEDIAN_CD4/(double)CYCLES_PER_YEAR;
					if(pat->CD4real >= 350)
						retval += HOSPITAL_COST_POOR_HIGH_CD4/(double)CYCLES_PER_YEAR;
				}
#if DETECT_VL
				if (whichCost == TOTAL || whichCost == ALC_INTERV)
				{
					//If patient is on alc intervention.  
					if (pat->interv_on)
						retval += COST_ALC_INTERVENTION / (double)CYCLES_PER_YEAR;
				}
#endif

#if RETENTION_IN_CARE
			}
#endif		
		

		if (whichCost == TOTAL || whichCost == DRUG)
		{
			//Add in drug costs
			//patRegCombo indicates which line regimen the patient is on
			//Note these were given in $/month, so we multiply by 12 then divide by CYCLES_PER_YEAR
			if (pat->haart)
			{
				if (patRegCombo==0) retval += cost_reg1_poor * 12.0 / (double)CYCLES_PER_YEAR;
				if (patRegCombo==1) retval += cost_reg2_poor * 12.0 / (double)CYCLES_PER_YEAR;
				if (patRegCombo==2) retval += cost_reg3_poor * 12.0 / (double)CYCLES_PER_YEAR;	
			}
		}

#if RETENTION_IN_CARE	
		if ((whichCost == TOTAL || whichCost == OUTREACH) && outreach_interv_enabled) //if patient dies, costs are calculated in function process_death
		{ 
			if (pat->outreach && (pat->CD4real <= outreach_cd4))
			{
				// if outreach condition is satisfied, accumulate finding cost
				retval += outreach_finding_cost_poor;
			}
		}

		if ((whichCost == TOTAL || whichCost == RISK_REDUCTION) && risk_reduction_interv_enabled)
		{
			//RR cost applys only if pat is in ES1
			if (pat->engaged_in_care && pat->engaged_in_clinic && (pat->CD4real <= risk_reduction_cd4)) 
			{   // if outreach condition is satisfied, accumulate intervention cost
				retval += risk_reduction_intervention_cost;
			}
		}

		if ((whichCost == TOTAL || whichCost == SECONDARY_PREVENTION) && secondary_prevention_interv_enabled)
		{
			//SP cost applys only if pat is in ES1
			if (pat->engaged_in_care && pat->engaged_in_clinic && (pat->num_cycles_ES4 > 0) && (pat->CD4real <= secondary_prevention_cd4)) 
			{   // if outreach condition is satisfied, accumulate intervention cost
				retval += secondary_prevention_intervention_cost;
			}
		}
#endif
	}
	
	return retval;
}

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS
int ReadSensitivityFile ()
{
	static bool moreLinesToRead = true;
	static bool firstRead = true;
	string headings;
	static ifstream infile;

	if (firstRead) 
	{
		infile.open ("OneWaySensitivity.txt");	
		if (infile.fail()) 
		{
			cout << "Error opening OneWaySensitivity.txt" << endl;
			exit(1);
		}
		getline(infile, headings);
		senseFile << "\tStart age\tStart CD4\tCD4_treat\tStart_hiv\tHaart_tox\tUtil_dec\tNum_patients\tHIV Deaths\tAge Deaths\tMedian Survival Time\tMean Survival Time\tMedian QALYS\tMean QALYS\tmean time before salvage\tmean num regs before salvage\tmean res drugs before salvage\ttttf1\ttttf2\ttttf3\tnum exhausted all regs\tnum res from initial muts PI_SING\tnum res from initial muts PI_BOOSTED\tnum res from initial muts NRTI_TAM\tnum res from initial muts NRTI_NONTAM\tnum res from initial muts NNRTI_EFAV\tnum res from initial muts NNRTI_NEVIR" << endl << endl;
		firstRead = false;
	}

	//Start reading first line
	infile >> cd4_treat;

	//If we're at the end of the file, return false;
	if (infile.eof()) return false;

	//If we're here, the must be a line, so continue reading it!
	infile >> start_cd4;
	infile >> avgCD4_SD;
	infile >> start_hiv;
	infile >> avgHIV_SD;
	infile >> start_age;
	infile >> avgAge_SD;
	infile >> mut_pi_singular_start;
	infile >> mut_pi_boosted_start;
	infile >> mut_nrti_tam_start;
	infile >> mut_nrti_nontam_start;
	infile >> mut_nnrti_efavirenz_start;
	infile >> mut_nnrti_nevirapine_start;

	infile >> pmutres_pi_singular_est;
	infile >> pmutres_pi_boosted_est;
	infile >> pmutres_nrti_tam_est;
	infile >> pmutres_nrti_nontam_est;
	infile >> pmutres_nnrti_efavirenz_est;
	infile >> pmutres_nnrti_nevirapine_est;

	infile >> pcrossres_pi_singular_est;
	infile >> pcrossres_pi_boosted_est;
	infile >> pcrossres_nrti_tam_est;	
	infile >> pcrossres_nrti_nontam_est;
	infile >> pcrossres_nnrti_efavirenz_est;
	infile >> pcrossres_nnrti_nevirapine_est;	

	infile >> mut_rate;
	infile >> comp;
	infile >> assoc_comp;
	infile >> start_with_aids;
	infile >> pcrossres_piboosted_singular;
	infile >> pcrossres_pisingular_boosted;
	infile >> pcrossres_tam_nontam;
	infile >> pcrossres_nontam_tam;
	infile >> pcrossres_efav_nevir;
	infile >> pcrossres_nevir_efav;


	getline (infile, label);

	if (!label.empty()) label = label.substr(1,label.size());	//Remove the leading tabs
	return true;
}
#endif


/*Possible regimens

NNRTI_EFAVIRENZ			NNRTI_NEVIRAPINE
TAM + NONTAM			TAM + NONTAM
2 NONTAM				2 NONTAM

PI_BOOSTED				PI_SINGULAR
TAM + NONTAM			TAM + NONTAM	
2 NONTAM				2 NONTAM
*/

/*
A patient starts out on their initial regimen.  When they become intolerant or meet the trigger to change the reg, we 
choose a new regimen.  Initially, we choose a new regimen with no resistant or intolerant drugs.  If a new "clean" reg
(one that contains no resistant or intolerant drugs) cannot be built, the patient stays on their last regimen.  
They can change off of that regimen only if they have developed resistance to two or more drugs in the regimen.  
If this happens, we put them on a reg with the least number of resistant drugs regardless of intolerance.
*/

int ChooseRegimen_Rich(patient *pat) 
{
	int tReg[3]={0,0,0};					//tReg stands for test regimen!
	int remReg[NUM_CLASS_COMBINATIONS][3];	//remReg stands for "remember reg"!
	int regimenIndex;
	int numGoodRegs = 0;
	int pickRegOutOfHat;
	int numEachTypeinReg[ARV_TYPES];
	int numRes = 3;
	int numIntol = 3;
	int leastNumRes;
	int leastNumIntol = 9999;			//initialize to very high num
	int greatestVLDec = -1;				//initialize to very low num
	int typeofleastnumres;
	bool retVal;	
	int n1, n2, n3;
	bool addRegToList;

	//initialize
	leastNumRes = pat->totalres;

	//Note that the regimen list is in order of most effective VLdec to least effective.
	for (regimenIndex = 0; regimenIndex < NUM_CLASS_COMBINATIONS; regimenIndex++)
	{	
		switch (regimenIndex) 
		{
		case (0): tReg[0]= NNRTI_EFAVIRENZ; tReg[1]= NRTI_TAM; tReg[2]= NRTI_NONTAM;  break;
		case (1): tReg[0]= NNRTI_EFAVIRENZ; tReg[1]= NRTI_NONTAM; tReg[2]= NRTI_NONTAM;  break;
		case (2): tReg[0]= PI_BOOSTED; tReg[1]= NRTI_TAM; tReg[2]= NRTI_NONTAM;  break;
		case (3): tReg[0]= PI_BOOSTED; tReg[1]= NRTI_NONTAM; tReg[2]= NRTI_NONTAM;  break;
		case (4): tReg[0]= NNRTI_NEVIRAPINE; tReg[1]= NRTI_TAM; tReg[2]= NRTI_NONTAM;  break;
		case (5): tReg[0]= NNRTI_NEVIRAPINE; tReg[1]= NRTI_NONTAM; tReg[2]= NRTI_NONTAM;  break;
#if CALIBRATE
		case (6): tReg[0]= PI_SINGULAR; tReg[1]= NRTI_TAM; tReg[2]= NRTI_NONTAM;  break;
		case (7): tReg[0]= PI_SINGULAR; tReg[1]= NRTI_NONTAM; tReg[2]= NRTI_NONTAM;  break;	
#endif
		}

		//Get the number of drugs of each class in the regimen being considered (This is set in array numEachTypeinReg[])
		Get_Num_Each_Type_In_Array(pat, DRUGS_IN_REG, tReg, numEachTypeinReg);

		//initialize 
		addRegToList = false;

		//first scenario:  pat is still on active haart
		if (!pat->exhausted_clean_regs)
		{
			//We make sure there are enough drugs of each type that the pat is not res or intolerant to
			n1 = numNotResOrIntol(pat, tReg[0]);
			n2 = numNotResOrIntol(pat, tReg[1]);
			n3 = numNotResOrIntol(pat, tReg[2]);

			if((n1 >= numEachTypeinReg[tReg[0]] &&
				n2 >= numEachTypeinReg[tReg[1]] &&
				n3 >= numEachTypeinReg[tReg[2]] )
				&& 	

				//Save the first eligible reg (because it will have the highest VLdec)
				//And also save any reg with the same first drug

				//**Just added the "OR CALIBRATE" because during calibration we will have a different first regimen
				//and subsequent regimens will be randomly chosen from all those not resistant or intolerant

				(numGoodRegs==0 || 
				(numGoodRegs > 0 && !CALIBRATE && tReg[0] == remReg[0][0])
				|| CALIBRATE))
			{
				numGoodRegs++;
				addRegToList = true;
			}
		}

		//pat has exhausted all regs
		else  
		{

			//We make sure there are enough drugs of each type that the pat is not intolerant to
			//We ignore intolerance if the patient is resistant to two or more drugs
			if(	(pat->totalres<2 &&
				numNotIntol(pat, tReg[0]) >= numEachTypeinReg[tReg[0]] &&
				numNotIntol(pat, tReg[1]) >= numEachTypeinReg[tReg[1]] &&
				numNotIntol(pat, tReg[2]) >= numEachTypeinReg[tReg[2]] )
				||
				pat->totalres>=2)
			{
				numRes = 
					max(0, (numEachTypeinReg[tReg[0]] - numNotRes(pat, tReg[0]))) + 
					max(0, (numEachTypeinReg[tReg[1]] - numNotRes(pat, tReg[1]))*(tReg[0]!=tReg[1])) +
					max(0, (numEachTypeinReg[tReg[2]] - numNotRes(pat, tReg[2]))*(tReg[2]!=tReg[1] && tReg[2]!=tReg[0]));

				if (numRes < leastNumRes)
				{
					numGoodRegs=1;				//this will be the first good reg - it gets array index 0
					leastNumRes = numRes;
					typeofleastnumres = tReg[0];
					addRegToList = true;
				}

				else if (numRes == leastNumRes && numRes < pat->totalres && 
					(CALIBRATE || (numRes == leastNumRes && tReg[0] == typeofleastnumres)))
				{
					numGoodRegs++;
					addRegToList = true;
				}
			}
		}

		if (addRegToList)
		{
			remReg[numGoodRegs-1][0] = tReg[0];
			remReg[numGoodRegs-1][1] = tReg[1];
			remReg[numGoodRegs-1][2] = tReg[2];
		}
	}

	//This chunk randomly chooses a regimen out of all the regimens saved.
	if (numGoodRegs > 0)	
	{
		retVal = NEW_REG;				//a new regimen was available		
		pickRegOutOfHat = (!VRT ? RandomNumber(numGoodRegs) : (int)unifSim[6]);

		//put patient on that specific regimen
		putPatOnReg(pat, remReg[pickRegOutOfHat]);
	}

	else retVal = MAINTAIN_CURRENT_REG;		//no reg available

	return retVal;	
}



int ChooseRegimen_Poor(patient *pat) 
{	
	drug *tReg[3][DRUGS_IN_REG];		//tReg stands for test regimen, NUM_CLASS_COMBINATIONS=3
	int numResCurrentReg=0, numIntolCurrentReg=0;
	int numRes, numIntol, bestReg=999, bestNumRes, bestNumIntol;
	int retVal=999;

	//needed a constant expression for the size of tReg above.  Previously used a #define, but 
	//switched to a variable, "numAvailableRegs" during the MonitorPaper Analysis.
	if (numAvailableRegimens>3)
	{
		cout<<"Need to change size of tReg array in ChooseRegimen_Poor\n"<<endl;
		exit(0);
	}

	//The regimens are in the order that they would be given in AMPATH
	//The initial regimen is NEVIRAPINE, NRTI_NONTAM, NRTI_TAM
	if (numAvailableRegimens >0){ tReg[0][0]= &drugs[INITIAL_REG_DRUG1][0]; tReg[0][1]= &drugs[INITIAL_REG_DRUG2][0]; tReg[0][2]= &drugs[INITIAL_REG_DRUG3][0];}
	if (numAvailableRegimens >1){ tReg[1][0]= &drugs[PI_BOOSTED][0]; tReg[1][1]= &drugs[NRTI_TAM][1]; tReg[1][2]= &drugs[NRTI_NONTAM][1];}
	if (numAvailableRegimens >2){ tReg[2][0]= &drugs[PI_BOOSTED][1]; tReg[2][1]= &drugs[NRTI_NONTAM][2]; tReg[2][2]= &drugs[NRTI_NONTAM][3];}

	//If the patient hasn't exhausted all clean regs, just put them on the next reg
	//Note that since we start the pat on pat->numreg==1, then the index of tReg for their NEXT regimen is 1 (or numreg!)
	if (!pat->exhausted_clean_regs)
	{
		if (pat->numreg<numAvailableRegimens)
		{	
			pat->reg[0] = tReg[pat->numreg][0];			
			pat->reg[1] = tReg[pat->numreg][1];
			pat->reg[2] = tReg[pat->numreg][2];

			retVal = NEW_REG;
			patRegCombo = pat->numreg;						//used in MONITORING_PAPER
		} 
		else 
			retVal = MAINTAIN_CURRENT_REG;
	}
	
	else if (poorCanGoBackOnPreviousRegs) //pat has exhausted all clean regs
	{
		//First, count the number res and intol in current regimen
		for (int i=0; i<DRUGS_IN_REG; i++)
		{
			numIntolCurrentReg += pat->reg[i]->intol;
			numResCurrentReg += pat->reg[i]->res;
		}

		//initialize the "best" reg to be the current one
		bestNumRes = numResCurrentReg;
		bestNumIntol = numIntolCurrentReg;

		//Then, count the number res and intol in each possible regimen
		//If the test regimen is better than the current regimen, save the test regimen as "best"
		for (int t=0; t<numAvailableRegimens; t++)
		{
			numIntol = 0;
			numRes = 0;

			//First, count the number res and intol in each possible regimen
			for (int d=0; d<DRUGS_IN_REG; d++)
			{
				numIntol += tReg[t][d]->intol;
				numRes += tReg[t][d]->res;
			}

			if (numRes < bestNumRes ||
				(numRes == bestNumRes && numIntol < numIntolCurrentReg)) 
			{
				bestReg=t;
				bestNumRes=numRes;
				bestNumIntol=numIntol;
			}
		}

		if (bestReg!=999) 
		{
			retVal = NEW_REG;

			//put patient on the best reg
			pat->reg[0] = tReg[bestReg][0];
			pat->reg[1] = tReg[bestReg][1];
			pat->reg[2] = tReg[bestReg][2];

			patRegCombo = bestReg;			//So we know which reg the pat went back to, if they returned to a previous reg

			if (KIMTEST) fprintf (kim_test, "The most favorable regimen is reg %d.\n", bestReg);
		}

		else retVal = MAINTAIN_CURRENT_REG;
	}

	//Patient had exhausted all regs, and wasn't allowed to go back on a previous reg
	else retVal = MAINTAIN_CURRENT_REG;			
		
	if (retVal==999) 
	{
		cout<<"Error in function: ChooseRegimen_Poor"<<endl;
		exit(0);
	}

	return retVal;	
}

void printDrugsInReg(patient *pat)
{
	for (int i=0; i<DRUGS_IN_REG; i++)
	{
		if (KIMTEST) fprintf (kim_test, "%s.%d  ", pat->className[pat->reg[i]->classNum], pat->reg[i]->drugNum);
	}
	if (KIMTEST) fprintf (kim_test, "\n");
}

void printNumRemaining(patient *pat)
{
	fprintf (kim_test,"\nRes\tIntol\tClass Name\n");

	for (int c=0; c<pat->arv_types; c++)
	{
		for (int d=0; d<pat->num_in_class[c]; d++)
		{
			fprintf (kim_test, "%d\t%d\t%s.%d\n", drugs[c][d].res, drugs[c][d].intol, pat->className[c], d );
		}
		if (pat->num_in_class[c] !=0) fprintf (kim_test, "\n");		//if there were no drugs in the class, prevent two endlines
	}
}

//The values for the VL decrements come from previous research.  Those values are true for an hiv_baseline of 4.5.  Thus, we scale the VL dec according to the patient's HIVbaseline.
void getVLdec(patient *pat)
{
	switch (pat->reg[0]->classNum)
		{
		case (PI_SINGULAR): 
			pat->VLdec = VL_DEC_PI_SINGULAR * pat->HIVbaseline / 4.5 + DVLHAART;  
			break;	
		case (PI_BOOSTED): 
			pat->VLdec = VL_DEC_PI_BOOSTED * pat->HIVbaseline / 4.5 + DVLHAART;  
			break;
		case (NNRTI_EFAVIRENZ): 
			pat->VLdec = VL_DEC_EFAVIRENZ * pat->HIVbaseline / 4.5 + DVLHAART;  
			break;
		case (NNRTI_NEVIRAPINE): 
			pat->VLdec = VL_DEC_NEVIRAPINE * pat->HIVbaseline / 4.5 + DVLHAART;  
			break;
		}
	pat->VLdec *= VL_DEC_MULT; //LL: VL_DEC_MULT is lately added. See Jason's email (April 07, 2015 10:24am)

	if (MARK_ROBERTS_AHRQ)
	{
		pat->VLdec *= MARK_VL_DEC_MULTIPLIER;
	}
#ifdef MONITORING_PAPER_SENSITIVITY
	if(floor((double)runnum/2) == 2) pat->VLdec *= sens[2][runnum % 2];
#endif
}

//Example:  RandomNumber(8) will return a random number between 0 and 7.
int RandomNumber(int upperRange)
{
	//this line has been replaced by the following:  return (rand() % upperRange);
	return (int)floor(uniform_a_b(0, upperRange-.0000000001, &seed));
}

void decrementForIntolerance(patient *pat)
{
	int j, k;

	//Before we put the patient on a new regimen, remove the drugs the patient is intolerant to from the drugs remaining in each class.  Part of the time, we will decrement for all, otherwise, we choose a drug at random to decrement for.

	//First case - we decrement for all drugs in the regimen
	if ((!VRT ? uniform(&seed) : unifSim[7]) <= INTOLERANCE_DING_FOR_ALL)
	{
		for (j = 0; j < DRUGS_IN_REG; j++)
		{
			pat->reg[j]->intol=true;
		}

		if (KIMTEST)
		{
			fprintf (kim_test,"Patient was intolerant to ALL drugs.  Decrementing for all.\n "); 									
			printDrugsInReg(pat);
			printNumRemaining(pat);
		}
	}

	else	//only decrementing for one drug in the regimen.  choose it at random
	{
		k = (!VRT ? RandomNumber(3) : (int)unifSim[8]);

		if (!pat->reg[k]->intol)
		{
			pat->reg[k]->intol= true;

			if (KIMTEST)
			{
				fprintf (kim_test,"Patient was intolerant.  Choosing one drug at random to decrement: %s.%d\n", pat->className[pat->reg[k]->classNum], pat->reg[k]->drugNum); 										
				printNumRemaining(pat);
			}
		}
		else if (KIMTEST) fprintf (kim_test,"Chose one drug at random to decrement but pat was already intol: %s.%d\n", pat->className[pat->reg[k]->classNum], pat->reg[k]->drugNum);
	}
}


/*The function getInitialAIDSRates() is used when START_WITH_AIDS (the percentage of patients starting with AIDS) 
is set to a non-zero value.  If some portion of patients are starting the model with AIDS, we want to distribute 
those patients across CD4 strata according to the hiv mortality table.  We are using "middle of the road" data for 
the lookup table:  Patient not on haart, 40<age<50, 2nd year, 4.5<=VL<5.5.  Using these values, gives us two strata for CD4.

2213 (CD4 < 50)		0.2241
2223 (50<=CD4<=200) 0.2241
2233 (200<=CD4<=350)0.027
2243 (350<=CD4<=500)0.027
2253 (CD4>=500)		0.027

Thus, we have two strata:  CD4<=200 and CD4>200.  This explains the significance of the 200 in the function below.  
If this table changes, the function should be updated accordingly.

NumPats CD4<=200 * AIDSrate<=200 + NumPats CD4>200 * AIDSrate>200 = START_WITH_AIDS * Total NumPats

This is best illustated with an example.

Ex:  

AVG_CD4 = 250
AVG_CD4_SD = 50
START_WITH_AIDS = .2

From normal distribution,	P(<200) = .1587
P(>200) = .8413

The ratio of AIDS rates (<200 : >200) = 0.2241:0.027

The percentage of patients with CD4 <= 200 starting with AIDS must be 
.2241*AIDS_EVENT_MULTIPLIER*x
The percentage of patients with CD4 >200 starting with AIDS must be 
.027*AIDS_EVENT_MULTIPLIER*x

(Note that the AIDS_EVENT_MULTIPLIER isn't really necessary here since the ratios would be the same without it.)

.1587N((.2241*AIDS_EVENT_MULTIPLIER)*x) + .8413N((.027*AIDS_EVENT_MULTIPLIER)*x) = .2N

The Ns cancel leaving,

.1587((.2241*AIDS_EVENT_MULTIPLIER)*x) + .8413((.027*AIDS_EVENT_MULTIPLIER)*x) = .2
*/
void getInitialAIDSRates()
{

	double x;
	double prob_less_200;			//Given uniform distribution
	double prob_greater_200;		//Given uniform distribution

	double p_aids_less_200 = .2241 * AIDS_EVENT_MULTIPLIER;
	double p_aids_greater_200 = .027 * AIDS_EVENT_MULTIPLIER;

	if (avgCD4_SD != 0)
	{
		if (start_cd4 > 200) prob_less_200 = 1-cdf(-(200-start_cd4)/avgCD4_SD);
		else prob_less_200 = cdf((200-start_cd4)/avgCD4_SD);

		prob_greater_200 = 1-prob_less_200;

		x = start_with_aids / (prob_less_200*p_aids_less_200 + 
			prob_greater_200 * p_aids_greater_200);

		AIDSrate_less200 =  p_aids_less_200 * x;
		AIDSrate_greater200 = p_aids_greater_200 * x;

		if (AIDSrate_less200 > 1) 
		{
			AIDSrate_less200 =  1;
			AIDSrate_greater200 = (start_with_aids - prob_less_200) / prob_greater_200;
		}
	}

	else AIDSrate_less200 = AIDSrate_greater200 = start_with_aids;
}

void Get_Num_Each_Type_In_Array(patient *pat, int numInArray, int drugArray[], int numEachType[])
{
	int i, j;

	//initialize
	for (i = 0; i < ARV_TYPES; i++)
		numEachType[i] = 0;

	for (i = 0; i < numInArray; i++)
	{
		//for each drug, we find the arv type of this drug and increment the
		//appropriate counter
		for (j = 0; j < ARV_TYPES; j++)
		{
			if (drugArray[i] == j) 
				numEachType[j]++;
		}
	}
}

int numNotResOrIntol(patient *pat, int type)
{
	int numNotResOrIntol=0;

	for (int i=0; i<pat->num_in_class[type]; i++)
	{
		if (!drugs[type][i].intol && !drugs[type][i].res) numNotResOrIntol++;
	}
	return numNotResOrIntol;
}

int numNotIntol (patient *pat, int type)
{
	int numNotIntol=0;

	for (int i=0; i<pat->num_in_class[type]; i++)
	{
		if (!drugs[type][i].intol) numNotIntol++;
	}
	return numNotIntol;
}

int numNotRes (patient *pat, int type)
{
	int numNotRes=0;

	for (int i=0; i<pat->num_in_class[type]; i++)
	{
		if (!drugs[type][i].res) numNotRes++;
	}
	return numNotRes;
}

void getCD4regressionideal(patient *pat)
{	
	static double total_ended_round1;
	static double total_begin_round2;
	static int num_ended_round1;
	static int num_begin_round2;

	//ONE REG
	pat->CD4regressionideal = 

		pat->CD4atStartOfHaart +

		CD4_MULT*(LUCK_MULTIPLIER * (WEIGHT_PAT_LUCK*pat->luck1 + WEIGHT_ROUND_LUCK*pat->luck2) * AVG_CD4_SD
		+
		(pat->numreg>=1)*(
		+
		(37-35)*(((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR)>=0)		
		+
		(33-10)*(((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR)>=1)		
		+
		(20-10)*(((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR)>=2)
		+
		14*(((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR)>=3)
		+		
		10*(((pat->cyclenum - pat->haart_start_time)/CYCLES_PER_YEAR)>=4)
		+
		0*(pat->AgeRegBaseline<30) 
		+
		26*(pat->AgeRegBaseline>=30 && pat->AgeRegBaseline<40)
		+
		18*(pat->AgeRegBaseline>=40 && pat->AgeRegBaseline<50)
		+
		16*(pat->AgeRegBaseline>=50 && pat->AgeRegBaseline<60)
		+
		6*(pat->AgeRegBaseline>=60 && pat->AgeRegBaseline<70)
		+
		-15*(pat->AgeRegBaseline>=70)
		+ 
		1*(pat->CD4RegBaseline<=200)
		+
		0*(pat->CD4RegBaseline>200 && pat->CD4RegBaseline<=350)
		+ 
		-6*(pat->CD4RegBaseline>350 && pat->CD4RegBaseline<=500)
		+
		-48*(pat->CD4RegBaseline>500)
		+
		3.3		//average of med class coefficients
		+
		0*(pat->comp>0.8)  
		+
		0*(pat->comp>0.6 && pat->comp<=0.8)
		+
		-6*(pat->comp>0.4 && pat->comp<=0.6)
		+
		-31*(pat->comp>0.2 && pat->comp<=0.4)
		+
		-51*(pat->comp<=0.2)
		+	
		//We decided not to use this line because it benefited patients who started with higher Viral Loads.  (Patients who started with high VLs experience more of a drop in VL, so got more of a benefit in cd4 bump)
        //-43*(pat->VLreal - pat->VLatStartOfHaart)
		//The getVLTerm function scales the size of the VL term for the possible size of VL increase
		43*getVLTerm(pat)
		));
#ifdef MONITORING_PAPER_SENSITIVITY
	if(floor((double)runnum/2) == 3) pat->CD4regressionideal += sens[3][runnum % 2];
#endif	
	pat->CD4regressionideal = Max (0, pat->CD4regressionideal);		//Don't let it fall below zero
}


double addNoise (double input, patient *pat, int type)
{	

	//These values are much the result of trial and error and make the CD4 look realistic.
	if (type==CD4)
	{
		//if (CYCLE==DAY) return (input + 100 * PerlinNoise1D(pat->cyclenum*.002*12,1.5,5,2));	//12 is the value that looks good and looks like the VL
		if (CYCLE==DAY) return (input + 100 * PerlinNoise1D((pat->cyclenum +offset1)*.002*5,1.5,5,2));		//this line give smoother noise than the line above
		else if (CYCLE==MONTH) return (input + 200 * PerlinNoise1D(pat->cyclenum*.002,1.5,5,2));
	}

	else if (type==VL)
	{
		if (CYCLE==DAY) return (input + 1 * PerlinNoise1D((pat->cyclenum +offset2)*.002*12,1.5,5,2));  //the 333 is an offset so the noise for the CD4 doesn't look the same as the VL
		else if (CYCLE==MONTH) return (input + .25* PerlinNoise1D(pat->cyclenum*.002,1.5,5,2));
	}

	else return (0);

}

void fillUnifInit()
{
	int index = 0;
	int i;

	for ( i = 0; i < 3; i++ )
	{
		unifInit[index++] = normal(0,1,&seed);		//1-3
	}

	unifInit[index++] = uniform(&seed);				//4
	unifInit[index++] = uniform(&seed);				//5

	unifInit[index++] = normal(0, 1, &seed);		//6
	unifInit[index++] = normal(0, 1, &seed);		//7
	unifInit[index++] = uniform(&seed);				//8

	pat_rns_used = MAX_SEED_INIT;

	printf ("First 3 nums for new pat:  %f %f %f \n", unifInit[0], unifInit[1], unifInit[2]);
}

void fillUnifSim()
{
	int index = 0;
	kludge_counter = 0;

	unifSim[0] = uniform(&seed);		
	unifSim[1] = uniform(&seed);		
	unifSim[2] = uniform(&seed);	
	unifSim[3] = uniform(&seed);	
	unifSim[4] = normal(0, 1, &seed);
	unifSim[5] = uniform(&seed);
	unifSim[6] = 0;
	unifSim[7] = uniform(&seed);	
	unifSim[8] = RandomNumber(3);
	unifSim[9] = int(monitor_interval/12.0*CYCLES_PER_YEAR);
	unifSim[10] = uniform(&seed);		
	unifSim[11] = uniform(&seed);		
	unifSim[12] = uniform(&seed);	


	kludge_counter = 13;

	for (index = 0; index <501; index++)
	{
		unifSim[13+index] = uniform(&seed);	
	}


	pat_rns_used += MAX_SEED_CYCLE - 2;	//VRT  Kim moved this out of the above loop!

	//if (pat_rns_used < 40) printUnif();
}


void printUnif()
{
	int i;
	for ( i = 0; i < (MAX_SEED_INIT); i++ )
	{
		printf("Value %d in Unif Init array: %lf\n",i,unifInit[i]);
	}

	for ( i = 0; i < (MAX_SEED_CYCLE); i++ )
	{
		printf("Value %d in Unif Sim array: %lf\n",i,unifSim[i]);
	}
}


void setNextMonitorCycle(patient *pat)
{
	int cyclesInInterval;
	int rand;

	//interval is in months and CYCTIME is cycles/month, result of division is cycles
	cyclesInInterval = int(monitor_interval/12.0*CYCLES_PER_YEAR);
	rand = (!VRT ? RandomNumber(cyclesInInterval): (int)unifSim[9]);

	pat->nextMonCycle = int(cyclesInInterval/2) + rand + pat->cyclenum;

#if (WHO_MONITORING_STRATEGY_ANALYSIS || defined(WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS))
		//For scenario 16 of this analysis, we want a one-time VL check, so set next monitor time to an impossibly distant date.
		if (WHO_MONITORING_STRATEGY_ANALYSIS && scenario==16 && pat->haart_start_time != pat->cyclenum) pat->nextMonCycle = 99999999;
		else pat->nextMonCycle = pat->cyclenum + cyclesInInterval;
#endif
}


//The costs in the function metTrigger are From Kara Wools-Kaloustian, per her e-mail dated 11/18/2009.  
//At AMPATH the cost of a CD4 test is $11.20, and VL is $70.
bool metTrigger(int trigger, patient *pat)
{
	bool retVal = false;				//indicates whether trigger to change regimen has been met
	const double CD4differential = 50.0;
	double testcost = 0;					//iniitalize the cost of any tests

	switch (trigger)
	{
	
	case T_VL:

		//Set the cost of the test
		testcost = testCostVL;

		if (pat->VLreal>=pat->vl_regimen_failure)
		{
			if (KIMTEST) fprintf(kim_test, "\nRegimen failed at cycle %d.  VL above threshold of %.1f.  VL=%.1f\n", pat->cyclenum, pat->vl_regimen_failure, pat->VLreal);

			//if this is the first time (vl_regimen_failure is still set to 1st), 
			//reset the vl_regimen_failure to after_1st
			if (pat->vl_regimen_failure == vl_regimen_failure_1st) 
				pat->vl_regimen_failure = vl_regimen_failure_after_1st;

			retVal = true;
		}

		if (WHO_MONITORING_STRATEGY_ANALYSIS)
		{
			incrementWHOArray(pat, numVLTests);

			if (trig != T_NESTED)
			{
				//Note - Scott says to implement "confirmatory retesting" where any positive VL test is followed up by another VL test
				//We ONLY do this if the regimen trigger is VL, NOT if it's NESTED.
				if (retVal == true) incrementWHOArray(pat, numVLTests);
			}
		}

		break;


	case T_CD4:

		//Set the cost of the test
		testcost = testCostCD4;

		//		if ((pat->CD4real < .5 * pat->maxCD4onReg ||
		//			pat->CD4real <= (pat->CD4RegBaseline - CD4differential) ||
		//			pat->CD4real < 100) && 
		//			!(CALIBRATE && pat->numreg == MAX_REGS))
		/*in March 2015, we tried to calibrate for RIC, Jason suggested that the new guideline for regimen change is different, so we change our criteria accordingly*/
		if ((pat->CD4real <= (pat->CD4RegBaseline - CD4differential) ||
			pat->CD4real < 100 && pat->lastLess100) &&
			!(CALIBRATE && pat->numreg == MAX_REGS))
		{
			if (KIMTEST) fprintf(kim_test, "\nRegimen failed at cycle %d.  CD4real: %.1f, peak CD4 on round: %.1f, round baseline cd4: %.1f\n",  pat->cyclenum, pat->CD4real, pat->maxCD4onReg, pat->CD4RegBaseline);

			retVal = true;
		}

		if (WHO_MONITORING_STRATEGY_ANALYSIS) incrementWHOArray(pat, numCD4Tests);

		break;

		
	case T_NESTED:

		//If CD4 failure criteria is met, then check VL
		if (trigger==T_NESTED && metTrigger(T_CD4, pat))
		{
			if (metTrigger(T_VL, pat))
			{
				if (KIMTEST) fprintf(kim_test, "\nRegimen failed at cycle %d with nested trigger\n",  pat->cyclenum);

				//Note, I don't add in the costs here because the costs are added when I check each individual trigger.

				retVal = true;
			}
		}

		break;


	case T_CLINICAL: 

		if (pat->AIDSeventReg)
		{
			if (KIMTEST) fprintf(kim_test, "\nRegimen failed at cycle %d with clinical trigger\n",  pat->cyclenum);

			retVal = true;
		}

		break;
	}

	//add the costs of any tests to the total test costs
	total_cost += testcost;
	total_disc_cost += testcost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));

	total_lab_cost += testcost;
	total_lab_disc_cost += testcost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));

	return retVal;
}

//giving this a similar call to kolmogorov
//d is the meanSquareError;

void Least_Squared(double mod[], double clin[], int size, double *d)
{
	double sum=0;
	double dist;

	for (int i=0; i<size; i++)
	{
		dist=mod[i]-clin[i];
		sum+=dist*dist;
	}

	//*d is what we're returning

	*d = sqrt(sum/double(size));
}


//Sometimes with testing I want to set the number of drugs in a certain drug class to 0.
//This function checks to see if there are enough drugs to make the initial regimen.
bool notEnoughDrugsForInitialReg(patient *pat)
{
	bool changeFirstReg=false;

	Get_Num_Each_Type_In_Array(pat, DRUGS_IN_REG, pat->init_reg, numEachTypeinReg);

 	for (int i = 0; i < DRUGS_IN_REG; i++)
	{
		if (numEachTypeinReg[pat->init_reg[i]] > pat->num_in_class[pat->init_reg[i]]) changeFirstReg = true;
	}

	return changeFirstReg;
}

int countResDrugs(patient *pat)
{
	int numResDrugs = 0;

	for (int c = 0; c < ARV_TYPES; c++)
	{
		for (int d=0; d<pat->num_in_class[c]; d++)
		{
			numResDrugs += drugs[c][d].res;
		}
	}

	return numResDrugs;
}

//Not currently used but leaving it in just in case!
/*void getSmoothedHIVDeathRate(patient *pat)
{
	/*  From GetCategory()

	if (pat->VLreal < 3.5)
	return 1;
	else if (pat->VLreal < 4.5)
	return 2;
	else if (pat->VLreal < 5.5)
	return 3;
	else
	return 4;
	*/

/*	long int hiv_index, hiv_index1, hiv_index2;
	double deathrateHIV1, deathrateHIV2;

	//if patient is in category	1 (VL < or 4, use the value in the table (the regular procedure)
	if (pat->VLcat == 1 || pat->VLcat == 4)
	{
		hiv_index = (long)HIV_Lookup(pat);
		pat->deathrateHIV = Get_Rate(hiv_mort_table, num_mort_entries, hiv_index);
	}

	//if patient is in category 2 or 3, use linear interpolation
	else
	{
		if (pat->VLcat == 2) 
		{
			hiv_index1 = (long)HIV_Lookup(pat);
			hiv_index2 = hiv_index1 + 1;
		}

		else if (pat->VLcat ==3)
		{
			hiv_index2 = (long)HIV_Lookup(pat);
			hiv_index1 = hiv_index2 - 1;
		}

		else 
		{
			cout<<"Error:  pat->VLcat was not a value between 1-4"<<endl;
			exit(0);
		}


		deathrateHIV1 = Get_Rate(hiv_mort_table, num_mort_entries, hiv_index1);
		deathrateHIV2 = Get_Rate(hiv_mort_table, num_mort_entries, hiv_index2);

		pat->deathrateHIV = linearInterpolate(3.5, deathrateHIV1, 5.5, deathrateHIV2, pat->VLreal);
	}
}*/

//from (y-y0)/(x-x0) = (y1-y0)/(x1-x0) we get
//y = y0 + (x-x0)(y1-y0)/(x1-x0)
double linearInterpolate(double x0, double y0, double x1, double y1, double x)
{
	double y;		//the return value

	y = y0 + (x-x0)*(y1-y0)/(x1-x0);

	return y;
}

/*******************************************************************************************************
Used if we want to do batch runs
*******************************************************************************************************/
void readCommandLineArgs(int argc, char* argv[])
{
#if BATCH
	//This doesn't seem to work right.  Needs testing before using
	start_age = atof(argv[1]);
	start_cd4 = atof(argv[2]);
	start_hiv = atof(argv[3]);
	cd4_treat = atof(argv[4]);
	initial_seed = atoi(argv[5]);

	avgAge_SD = 0;
	avgHIV_SD = 0;
	avgCD4_SD = 0;
	
#endif
}

/*This function is used to calculate the viral load term for the CD4 regression 
equation.  The initial VL term from Joyce's regression was 
-43*(pat->VLreal - pat->VLatStartOfHaart).  This term resulted in the 
individual patient trajectories showing a good response to VL.  However, patients 
with higher baseline viral loads were resulting in a higher CD4 and thus having longer
survivals.  So instead of using the original term, we're taking the VL decrease seen 
in VACS (at that VLfract, drugclass, etc.), comparing it to the maximum possible
VL decrease (perfect adherence, no resistance, on Efavirenz) for patients with that 
baseline VL, and then scaling that term for the current VL decrease as a proportion of 
max possible decrease.

Here is the logic I e-mailed to Scott:

First, some values: 

VACS mean VL baseline:  4.07
Efavirenz VL decrement (unadjusted):  3.09
PI_Boosted VL decrement (unadjusted): 2.68 

Equation for decrease in VL load = VLfract * VLdec (unadjusted) * VLbaseline / 4.5

Let’s say that in the VACS data, the average decrease in VL load on therapy was 1.58 (we can discuss what that true value should be.)

The maximum decrease in VL would be if the patient were on Efavirenz, and their VL_fract was 1.0. 

So the maximum decrease was
= 1.0 (3.09) (4.07) / 4.5
= 2.79

1.58 as a proportion of the maximum is 1.58/2.79 = 0.566

Now let’s supposed we are running a patient with a VL baseline of 6.0, on PI_Boosted, with a VL_fract of 1/3.

Their decrease in VL would be 

= 1/3 * (2.68)(6.0)/4.5
= 1.19

Their maximum possible decrease in VL (if VLfract = 1.0 and pat were on Efavirenz) would be 

= 1.0 * (3.09) (6.0)/4.5
= 4.12

This patient’s VL decrease as a proportion of the maximum is 1.19 / 4.12 = .289

What value do we use for the decrease in VL load in the CD4 equation?

= 1.58 (.289) / .566
= 0.807

Note:  The above value of 1.58 was a guesstimate.  The actual mean decrease 
in VL for VACS is calculated using a weighted average of the data
in the 2007 AIDS paper by Scott, calculated to be 1.348.

1.348 = .182(2.07) + .082(1.44) + .512(1.12) + .065(1.58) + .075(1.22) + .084(1.02)

Thus, 1.348 should replace the 1.58 in the examples above.
*/



double getVLTerm(patient *pat)
{
	const double meanVLdecreaseVACS = 1.348;		//from 2007 AIDS paper by R. Braithwaite
	const double meanBaselineVLVACS = 4.07;			//mean VL baseline from Joyce's VACS data
	double VLdecreasePat, maxVLdecreaseVACS, maxVLdecreasePat, propOfMaxVACS, propOfMaxPat;
	double dec;	

	//The maximum possible decrease in Viral Load would occur at a vl_fract of 1.0
	//(indicating perfect adherence, no resistance) and on Efavirenz (which has the 
	//highest VLdec.  The 4.5 comes from the equation used in the getVLdec function.
	maxVLdecreaseVACS = 1.0 * VL_DEC_EFAVIRENZ * meanBaselineVLVACS / 4.5;
	maxVLdecreasePat  = 1.0 * VL_DEC_EFAVIRENZ * pat->HIVbaseline   / 4.5;	

	//Calculate the drop in VL from start of haart (VLatStartOfHaart should be close
	//to pat->HIVbaseline but could vary somewhat because of imposed noise
	VLdecreasePat = pat->VLatStartOfHaart - pat->VLreal;

	//Calculate the decrease in VL as a proportion of the maximum possible decrease
	propOfMaxVACS = meanVLdecreaseVACS / maxVLdecreaseVACS;
	propOfMaxPat  = VLdecreasePat      / maxVLdecreasePat;

	//We take the average decrement seen by the VACS patients and scale up or down
	//according to the ratio of the proportion of max decrease seen by the current patient
	//to the average proportion of max decrease seen by the current patient.
	//dec is the value that should replace (pat->VLatStartOfHaart - pat->VLreal)
	//in the CD4 regression equation
	dec = meanVLdecreaseVACS / propOfMaxVACS * propOfMaxPat;

	return dec;
}

//put patient on that specific regimen
void putPatOnReg(patient *pat, int *type)
{
	int drugNum;

	for (int i=0; i<DRUGS_IN_REG; i++)
	{
		drugNum = giveMeBestDrugInClass(pat, *(type+i),  i);
		pat->reg[i]=&drugs[type[i]][drugNum];
	}
}

/*For a given drug class, this function returns the next available "best" drug in the class.  
The order of preference is:
1)  Neither res nor intol
2)  Intol but not res
3)  Res but not intol
4)  Both res and intol
*/
int giveMeBestDrugInClass(patient *pat, int type,  int numAssignedSoFar)
{
	int drugFound=NONE_AVAILABLE;

	for (int i=0; i<pat->num_in_class[type] && drugFound==NONE_AVAILABLE; i++)
	{
		if (!drugs[type][i].res && !drugs[type][i].intol) 
		{
			drugFound=i;

			for (int j=0; j<numAssignedSoFar; j++) 
			{
				if (type==pat->reg[j]->classNum && i==pat->reg[j]->drugNum) drugFound=NONE_AVAILABLE;
			}
		}

	}

	if (drugFound==NONE_AVAILABLE)
	{
		for (int i=0; i<pat->num_in_class[type] && drugFound==NONE_AVAILABLE; i++)
		{
			if (!drugs[type][i].res) 
			{
				drugFound=i;

				for (int j=0; j<numAssignedSoFar; j++) 
				{
					if (type==pat->reg[j]->classNum && i==pat->reg[j]->drugNum) drugFound=NONE_AVAILABLE;
				}
			}
		}
	}

	if (drugFound==NONE_AVAILABLE)
	{
		for (int i=0; i<pat->num_in_class[type] && drugFound==NONE_AVAILABLE; i++)
		{
			if (!drugs[type][i].intol)
			{
				drugFound=i;

				for (int j=0; j<numAssignedSoFar; j++) 
				{
					if (type==pat->reg[j]->classNum && i==pat->reg[j]->drugNum) drugFound=NONE_AVAILABLE;
				}
			}
		}
	}

	if (drugFound==NONE_AVAILABLE) drugFound=0;

	return drugFound;
}

void checkForCrossResAndCrossClassRes(patient *pat, int mutType)
{
	for (int c=0; c<pat->arv_types; c++)
	{
		for (int d=0; d<pat->num_in_class[c]; d++)
		{
			//Check for cross-resistance
#ifdef MONITORING_PAPER_SENSITIVITY
			if(floor((double)runnum/2) == 4) pat->pcrossres[mutType] *= sens[4][runnum % 2];		
#endif			
			if (c==mutType && 
				!drugs[c][d].res &&
				(!VRT ? uniform(&seed) : unifSim[kludge_counter++]) <= pat->pcrossres[mutType])
			{
				//assign resistance to found drug
				drugs[c][d].res = true;

				if (KIMTEST) 
				{
					fprintf (kim_test, "Developed cross resistance to %s.%d\n", pat->className[mutType], d); 
					if (pat->haart) fprintf(kim_test, "Res profile: (%d %d %d)  \n", pat->reg[0]->res, pat->reg[1]->res, pat->reg[2]->res);
				}
			}

			//Check for cross class resistance 
			if (c == pat->crossClassType[mutType] && 
				!drugs[c][d].res && 
				(!VRT ? uniform(&seed) : unifSim[kludge_counter++]) <= pat->pCrossClassRes[mutType] )
			{

				//assign resistance to found drug
				drugs[c][d].res = true;

				if (KIMTEST) 
				{	
					fprintf (kim_test, "Developed (cross-class) resistance to %s.%d\n", pat->className[pat->crossClassType[mutType]], d);  					
					if (pat->haart) fprintf(kim_test, "Res profile: (%d %d %d)  \n", pat->reg[0]->res, pat->reg[1]->res, pat->reg[2]->res);
				}
			}
		}
	}

	if (KIMTEST) printNumRemaining(pat);
}

//A single call to HIV_Lookup and then Get_Rate() has been replaced by the following function.
//For patients with low compliance (< 0.5), we combined the off-haart mortality with on-haart 
//mortality.  At a compliance of 0, the patient's mortality should be the off-haart mortality and
//at 50% comp, the patient should have the full on-haart mortality.  
double Get_HIV_Death_Rate(patient *pat)
{
	double deathRateHaart, deathRateOffHaart, deathrateHIV, hiv_index;

	if (!pat->haart)
	{
		hiv_index = HIV_Lookup(pat, false);
		deathrateHIV = Get_Rate(hiv_mort_table, num_mort_entries, hiv_index);
	}

	else
	{
		hiv_index = HIV_Lookup(pat, true);
		deathRateHaart = Get_Rate(hiv_mort_table, num_mort_entries, hiv_index);
		if (pat->comp >= .5) deathrateHIV = deathRateHaart;
		else
		{
			hiv_index = HIV_Lookup(pat, false);
			deathRateOffHaart = Get_Rate(hiv_mort_table, num_mort_entries, hiv_index);
			deathrateHIV = linearInterpolate(0, deathRateOffHaart, 0.5, deathRateHaart, pat->comp);
		}
	}

	return deathrateHIV;
}


/*************************************************************************************
* This obtains the index to use to look up the HIV death rate in the hiv mortality table.
* The factors that affect this index are whether the patient is on HAART or not, her age
* category, her length of time since starting haart, her CD4 cat, and her VL cat.

2/24/2010 This function has been modified to received two variables, pat and haart.  When patients have a 
low compliance (< 0.5), we want to find their mortality rate both on and off therapy and then interpolate
so that the patient gets the full benefit of haart starting at the 
*************************************************************************************/
double HIV_Lookup(patient *pat, bool haart)
{
	/*********************************************************************************************************
	* We are changing to the following form which gives everyone the mortality benefit of being on haart
	* even before starting haart.  As Scott mentioned in an e-mail, "this is a conservative assumption that
	* biases the analyses in favor of starting later."  This has to do with "how we really don't know the
	* non-haart mortality rates for those with favorable CD4 and viral loads is."
	* CHANGE: We are no longer doing this.  That was just for a particular set of analyses.  Now, if a pat
	* is not on haart, they get no mortality benefit (first condition).  If they are on haart, then it depends
	* on whether or not there is resistance to all drugs and whether or not we are condisering viral fitness.
	* If there is not complete resistance to all regimens, then the benefit remains.  If there is res to all, and
	* we are considering reduced virulence of resistant strains, then the benefit remains. Otherwise, if there
	* is res to all and we do not consider a weaker strain, we remove the benefit of being on haart at that point.

	1st digit:							4th digit:
	1 Patient is on Haart				1 CD4 < 50
	0 Patient is not on Haart			2 50 <= CD4 < 200
	3 200 <= CD4 < 350
	2nd digit:							4 350 <= CD4 < 500
	1 age < 40							5 CD4 >=500
	2 40 < age < 50					
	3 age >= 50						5th digit:
	1 VL < 3.5
	3rd digit:							2 3.5 <= VL < 4.5
	1 1st year							3 4.5 <= VL < 5.5
	2 2nd year							4 VL > 5.5
	3 3rd year
	**************************************************************************/
	int total = 0;

	if (haart) total = 10000;	//In it's original form, this line was as follows, but was changed since we no longer use stay_on_haart:  total = 10000 - 10000*(pat->exhausted_clean_regs==1)*(pat->stay_on_haart==0);
	else total=0;

	if ( pat->age < 39.99 )
		total += 1000;
	else if ( pat->age < 49.99)
		total += 2000;
	else
		total += 3000;

	/*knv1  
	if (pat->cyclenum < CYCLES_PER_YEAR)
	total += 100;
	else if (pat->cyclenum < 2 * CYCLES_PER_YEAR)
	total += 200;
	else
	*/
	total += 300;

	total += 10*pat->CD4cat;
	total += pat->VLcat;
	return (double)total;
}



void writeTransitionsToFile()
{
	int cd4_start, vl_start, onTreat_start, res_start;	
	int cd4_dest, vl_dest, onTreat_dest, res_dest, deadOfHIV_dest;
	bool newline = 0;
	const bool verboseprint = 0;
	string cd4_cat[5] = {"CD4 < 50","50 <= CD4 < 200","200 <= CD4 < 350","350 <= CD4 < 500","500 <= CD4" };
	string VL_cat[5] = {"VL < 2.5","2.5 <= VL < 3.5","3.5 <= VL < 4.5","4.5 <= VL < 5.5","5.5 <= VL" }, OnTreat[2] = {"","Yes"};
	string living[2] = {"","Dead"}, resistStat[8] = {"No resistance","NNRTI only","PI only","NRTI only","NNRTI and NRTI","PI and NRTI","NNRTI and PI","All 3 classes"};
	const bool debugTransitions = false;
	double averageAnnualCost;

/*	The resistance strata defined:
	0: No resistance
	1: NNRTI only
	2: PI only
	3: NRTI only
	4: NNRTI and NRTI
	5: PI and NRTI
	6: NNRTI and PI
	7: All 3 classes
*/	
	if (debugTransitions) cout << "\tCD4\tVL\tTreat\tResis\t\tCD4\tVL\tTreat\tResis\tAlive\tNumber"<< endl;

	if (verboseprint) {
		fprintf(stateprob, "CD4\tViral Load\tTreatment\tResistance\t==Transition to ==>\tCD4\tViral Load\tTreatment\tResistance\tLiving?\tCount\n\n");
	}

	for (cd4_start = 0; cd4_start < STCD4; cd4_start++) {
		for (vl_start = 0; vl_start < STVL; vl_start++) {
			for (onTreat_start = 0; onTreat_start < STTREAT; onTreat_start++) {
				for (res_start = 0; res_start < STRESIST; res_start++) {
					for (cd4_dest = 0; cd4_dest < STCD4; cd4_dest++) {
						for (vl_dest = 0; vl_dest < STVL; vl_dest++) {
							for (onTreat_dest = 0; onTreat_dest < STTREAT; onTreat_dest++) {
								for (res_dest = 0; res_dest < STRESIST; res_dest++) {
									for (deadOfHIV_dest = 0; deadOfHIV_dest < STDEAD;  deadOfHIV_dest++) {
										if (statecount[cd4_start][vl_start][onTreat_start][res_start][cd4_dest][vl_dest][onTreat_dest][res_dest][deadOfHIV_dest] !=0) {

											if (debugTransitions) 
											{
												cout << "row\t "<<cd4_start << "\t" << vl_start << "\t" << onTreat_start << "\t" << res_start;
												cout <<"\tColumn:\t"<< cd4_dest<< "\t"<< vl_dest << "\t" << onTreat_dest << "\t" << res_dest << "\t" << deadOfHIV_dest << "\t";
												cout <<statecount[cd4_start][vl_start][onTreat_start][res_start][cd4_dest][vl_dest][onTreat_dest][res_dest][deadOfHIV_dest];
												cout << endl;
											}

											averageAnnualCost = totalAnnualCost[cd4_start][vl_start][onTreat_start][res_start][cd4_dest][vl_dest][onTreat_dest][res_dest][deadOfHIV_dest] / statecount[cd4_start][vl_start][onTreat_start][res_start][cd4_dest][vl_dest][onTreat_dest][res_dest][deadOfHIV_dest];
											fprintf (stateprob, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t\n", cd4_start, vl_start, onTreat_start, res_start, cd4_dest, vl_dest, onTreat_dest, res_dest, deadOfHIV_dest, statecount[cd4_start][vl_start][onTreat_start][res_start][cd4_dest][vl_dest][onTreat_dest][res_dest][deadOfHIV_dest], averageAnnualCost);

											newline = 1;
										}
									}
								}
							}
						}
					}
					if (newline) 
					{
						newline = 0;
						fprintf(stateprob, "\n");
					}
				}
			}
		}
	}
}




void getstate(patient *pat, int mystate[])
{
	/*
	mystate[0] = CD4 level, 
	mystate[1] = Viral load level, 
	mystate[2] = On treatment, 
	mystate[3] = Resistance category
	*/ 

	if (pat->CD4real < 50) {
		mystate[0] = 0;
	} else if (pat->CD4real < 200) {
		mystate[0] = 1;
	} else if (pat->CD4real < 350) {
		mystate[0] = 2;
	} else if (pat->CD4real < 500) {
		mystate[0] = 3;
	} else {
		mystate[0] = 4;
	}
	
	if (pat->VLreal < 2.5) {
		mystate[1] = 0;
	} else if (pat->VLreal < 3.5) {
		mystate[1] = 1;
	} else if (pat->VLreal < 4.5) {
		mystate[1] = 2;
	} else if (pat->VLreal < 5.5) {
		mystate[1] = 3;
	} else {
		mystate[1] = 4;
	}
	
	mystate[2] = pat->haart; // On treatment or not

	mystate[3] = getResistanceCategory(pat); // resistance category

	mystate[4] = 0;	//if we're calling this function, patient must be alive!  

/*
	if (!headerdone) {
		cout << "\nCD4 categories\n 0 if <50, 1 if < 200, 2 if < 350, 3 if < 500, 4 if >=500\n";
		cout << "VL categories\n 0 if < 2.5, 1 if < 3.5, 2 if < 4.5, 3 if < 5.5\n";
		cout << "Treatment: 0 = not on treatment, 1 = on treatment\n";
		cout << "Resistance: Number of drugs resistant\n";
		cout << "Alive: 0 = living, 1 = dead from other causes, 2 = dead of AIDS\n";
		cout << "\ncyclenum\tCD4\tVL\tTreatment\tResistance\t"<< "cd4state\tVLstate\tTreatmentstate\tResistancestate\n";
		headerdone = 1;
	}
	cout <<pat->cyclenum/365 <<"\t"<< pat->CD4real << "\t"<< pat->VLreal <<"\t"<< pat->haart << "\t"<< pat->totalres;
	cout << "\t"<<mystate[0]<< "\t"<< mystate[1]<< "\t"<< mystate[2]<< "\t"<< mystate[3]<< "\n";
*/	
}

void recordChange()
{
	int i;

	//Record another change from source bucket to destination bucket!
	statecount[cur_state[0]][cur_state[1]][cur_state[2]][cur_state[3]][new_state[0]][new_state[1]][new_state[2]][new_state[3]][new_state[4]]++;

	//Record annual cost over the past year
	//NEWtotalAnnualCost[cur_state[0]][cur_state[1]][cur_state[2]][cur_state[3]][new_state[0]][new_state[1]][new_state[2]][new_state[3]][new_state[4]] += costSoFarThisYearThisPat;
	totalAnnualCost[cur_state[0]][cur_state[1]][cur_state[2]][cur_state[3]][new_state[0]][new_state[1]][new_state[2]][new_state[3]][new_state[4]] += costThisYearThisPat;

	//Reset parameters for next year
	for (i = 0; i < 5; i++) 
	{
		cur_state[i] = new_state[i];
	}

	totalCostAtStartOfYearThisPat = total_cost;
}



//determines the resistance category for the lookup table for the transmission model
int getResistanceCategory(patient *pat)
{
	bool resToClass[3];			//Indices 0=NNRTI, 1=PI, 2=NRTI
	int numResClasses = 0;
	int retval;
	
	//initialize resistance to false
	for (int c=0; c<3; c++) resToClass[c]=false;
	
	for (int c=0; c<pat->arv_types; c++)
	{
		for (int d=0; d<pat->num_in_class[c]; d++)
		{
			//if patient is resistant, mark that we have resistance to that class
			if (drugs[c][d].res) 
			{
				if (!resToClass[0] && (c==NNRTI_NEVIRAPINE || c==NNRTI_EFAVIRENZ)) 
				{
					resToClass[0]=true;
					numResClasses++;
				}
				if (!resToClass[1] && (c==PI_BOOSTED || c==PI_SINGULAR))
				{
					resToClass[1]=true;
					numResClasses++;
				}
				if (!resToClass[2] && (c==NRTI_TAM || c==NRTI_NONTAM))
				{
					resToClass[2]=true;
					numResClasses++;
				}
			}
		}
	}
	
	/*Resistant classes 
	 0: No resistance
	 1: NNRTI only
	 2: PI only
	 3: NRTI only
	 4: NNRTI and NRTI
	 5: PI and NRTI
	 6: NNRTI and PI
	 7: All 3 classes
	 */
	
	switch (numResClasses)
	{
		case 0: 
			retval = 0; break;
		case 1: 
			if (resToClass[0]) {retval = 1; break;} 
			if (resToClass[1]) {retval = 2; break;} 
			if (resToClass[2]) {retval = 3; break;} 
		case 2: 
			if (resToClass[0] && resToClass[2]) {retval = 4; break;} 
			if (resToClass[1] && resToClass[2]) {retval = 5; break;} 
			if (resToClass[0] && resToClass[1]) {retval = 6; break;} 
		case 3: 
			retval=7; break;
	}	
	
	return retval;
}

//used for WHO_MONITORING_STRATEGY_ANALYSIS
void incrementWHOArray(patient *pat, int numOfSomething[WHO_NUM_TIME_PERIODS])
{
	int timeIndex=0;	//0:  Within 5 years, 1:  Within 10 years, 2:  Within 20 years
	int yearIndex;
	
	//for the 5,10,20 year summaries
	for (int t=5; t<=20; t=t*2)
	{
		if ((double)pat->cyclenum/(double)CYCLES_PER_YEAR <= t) 
		{
			numOfSomething[WHO_YEARLY_OUTPUT_NUM_YEARS + timeIndex]++;
		}

		timeIndex++;
	}

	//for the yearly data
	yearIndex = int((double)pat->cyclenum/(double)CYCLES_PER_YEAR);

	if (yearIndex < WHO_YEARLY_OUTPUT_NUM_YEARS) numOfSomething[yearIndex]++;

}

//used for WHO_MONITORING_STRATEGY_ANALYSIS
void WHO_print_year_by_year_stats()
{
	
	fprintf (who_monitoring, "\n\nScenario: %d\n", scenario);
	fprintf (who_monitoring, "Year\t");
	fprintf (who_monitoring, "Num patient visits\t");
	fprintf (who_monitoring, "Num CD4 tests\t");
	fprintf (who_monitoring, "Num VL tests\t");
	fprintf (who_monitoring, "Patient-years on 1st line therapy\t");
	fprintf (who_monitoring, "Patient-years on 2nd line therapy\t");
	fprintf (who_monitoring, "Total deaths\t");
	fprintf (who_monitoring, "Total HIV deaths\t");
	fprintf (who_monitoring, "Life-years lived\t");
	fprintf (who_monitoring, "Life-years lived on successful ART\t");
	fprintf (who_monitoring, "Life-years lived with AIDS\t");
	fprintf (who_monitoring, "Num developed AIDS\t");
	fprintf (who_monitoring, "\n");


	for (int t=0; t<WHO_YEARLY_OUTPUT_NUM_YEARS; t++)
	{
		fprintf (who_monitoring, "%d\t", t+1);
		fprintf (who_monitoring, "%d\t", numPatVisits[t]);                                                                                 
		fprintf (who_monitoring, "%d\t", numCD4Tests[t]); 
		fprintf (who_monitoring, "%d\t", numVLTests[t]);        
		fprintf (who_monitoring, "%.3f\t", (double)patCyclesOn1stLineTherapy[t]/(double)CYCLES_PER_YEAR);                                                                                       
		fprintf (who_monitoring, "%.3f\t", (double)patCyclesOn2ndLineTherapy[t]/(double)CYCLES_PER_YEAR);     
		fprintf (who_monitoring, "%d\t", totalDeaths[t]);
		fprintf (who_monitoring, "%d\t", totalHIVDeaths[t]);
		fprintf (who_monitoring, "%.3f\t", (double)lifeCyclesLived[t]/(double)CYCLES_PER_YEAR);
		fprintf (who_monitoring, "%.3f\t", (double)lifeCyclesOnSuccessfulART[t]/(double)CYCLES_PER_YEAR);
		fprintf (who_monitoring, "%.3f\t", (double)lifeCyclesWithAIDS[t]/(double)CYCLES_PER_YEAR);
		fprintf (who_monitoring, "%d\t", numEverDevelopedAIDS[t]);
		fprintf (who_monitoring, "\n");
	}
	
}

//used for WHO_MONITORING_STRATEGY_ANALYSIS
void initializeWHOstats()
{
	for (int i=0; i<WHO_NUM_TIME_PERIODS; i++)					//We have arrays of 3 for 5, 10 and 20 year data
	{
		numPatVisits[i]=0;
		numCD4Tests[i]=0;
		numVLTests[i]=0;
		patCyclesOn1stLineTherapy[i]=0;
		patCyclesOn2ndLineTherapy[i]=0;
		totalDeaths[i]=0;
		totalHIVDeaths[i]=0;
		lifeCyclesLived[i]=0;
		lifeCyclesOnSuccessfulART[i]=0;	//life years before occurrence of AIDS defining event
		lifeCyclesWithAIDS[i]=0;
		numEverDevelopedAIDS[i]=0;
	}
}

void RyanWhiteSetupPatient()
{
	//The following patients are targeted by the intervention
	//They have improved compliance and a greater portion are linked to care
	if (patnum < numTargeted) 
	{
		patIsTargeted=true;
		comp = min(baselineComp * adh_RR[interv], 1.0);
		numLinked = int((min(baselineLinkage * linkage_RR[interv], 1.0) * numTargeted));

		//if patient is NOT linked to care, give them a cd4_treat of 50 (reflects showing up in treatment once showing symptoms)
		if (patnum >= numLinked) cd4_treat=50;
		else cd4_treat = CD4_TREAT;
	}
	else
	{
		patIsTargeted=false;
		comp = baselineComp;
		numLinked = int(baselineLinkage * (PATIENTS - numTargeted));

		//if patient is NOT linked to care, give them a cd4_treat of 50 (reflects showing up in treatment once showing symptoms)
		if (patnum >= numTargeted+numLinked) cd4_treat=50;
		else cd4_treat = CD4_TREAT;
	}

	//Checking that I did it right
	if (comp==baselineComp) 
		numBaselineComp++;
	if (cd4_treat < 51 && cd4_treat >49) 
		totalNumCD4_50++;
	if (patIsTargeted && cd4_treat==CD4_TREAT)
		numTargetedLinked++;
	if (!patIsTargeted && cd4_treat==CD4_TREAT)
		numNotTargetedLinked++;
}


void set_regimen_and_triggers_WHO_discordance()
{
	switch (scenario)
	{
		case 0:			//Scenario B0
			numAvailableRegimens = 1;
			trig = T_CLINICAL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 0.0;
			monitor_interval = 0;
			break;
		case 1:			//Scenario B1
			numAvailableRegimens = 2;
			trig = T_CLINICAL;
			monitor_interval = 6;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 0.0;
			break;
		case 2:			//Scenario B2
			numAvailableRegimens = 2;
			trig = T_CD4;
			monitor_interval = 6;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 0.0;
			break;
		case 3:			//Scenario B4
			numAvailableRegimens = 2;
			trig = T_NESTED;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 6;
			break;
		case 4:			//Scenario B6
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 6;
			break;
		case 5:			//Scenario B7
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 12;
			break;
		case 6:			//Scenario B8
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 36;
			break;
		case 7:			//Scenario A1
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 2.7;
			monitor_interval = 6;
			break;
		case 8:			//Scenario A3
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.7;
			monitor_interval = 6;
			break;
		case 9:			//Scenario A4
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 4.0;
			monitor_interval = 6;
			break;
	}
	
}


void set_regimen_and_triggers_WHO()
{
	switch (scenario)
	{
		case 0:			//Scenario B0
			numAvailableRegimens = 1;
			trig = T_CLINICAL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 0.0;
			monitor_interval = 0;
			break;
		case 1:			//Scenario B1
			numAvailableRegimens = 2;
			trig = T_CLINICAL;
			monitor_interval = 6;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 0.0;
			break;
		case 2:			//Scenario B2
			numAvailableRegimens = 2;
			trig = T_CD4;
			monitor_interval = 6;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 0.0;
			break;
		case 3:			//Scenario B4
			numAvailableRegimens = 2;
			trig = T_NESTED;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 2.7;
			monitor_interval = 6;
			break;
		case 4:			//Scenario B4
			numAvailableRegimens = 2;
			trig = T_NESTED;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 6;
			break;
		case 5:			//Scenario B4
			numAvailableRegimens = 2;
			trig = T_NESTED;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.7;
			monitor_interval = 6;
			break;
		case 6:			//Scenario B6
			numAvailableRegimens = 2;
			trig = T_NESTED;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 4.0;
			monitor_interval = 6;
			break;
		case 7:			//Scenario B6
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 2.7;
			monitor_interval = 6;
			break;
		case 8:			//Scenario B6
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 6;
			break;
		case 9:			//Scenario B6
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.7;
			monitor_interval = 6;
			break;
		case 10:			//Scenario B7
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 4.0;
			monitor_interval = 6;
			break;
		case 11:			//Scenario B7
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 2.7;
			monitor_interval = 12;
			break;
		case 12:			//Scenario B7
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 12;
			break;
		case 13:			//Scenario B7
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.7;
			monitor_interval = 12;
			break;
		case 14:			//The penciled in scenario B8
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 4.0;
			monitor_interval = 12;
			break;
			
		case 15:			//New
			numAvailableRegimens = 2;
			trig = T_NESTED;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 6;
			break;
			
		case 16:			//New
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 6;
			break;
			
		case 17:			//New
			numAvailableRegimens = 2;
			trig = T_VL;
			vl_regimen_failure_1st = vl_regimen_failure_after_1st = 3.0;
			monitor_interval = 36;
			break;
	}
	
}

#ifdef WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
void record_discordance_data(patient *pat)
{
	int time_category;
	const int timepoint[4] = {CYCLES_PER_YEAR*0, CYCLES_PER_YEAR*5,CYCLES_PER_YEAR*10,CYCLES_PER_YEAR*15};

	for (time_category = 0; time_category < 4; time_category++)
	{
		if (pat->cyclenum < timepoint[time_category]) 	// cumulative_years will increment the current AIDS status/time category.
			cumulative_years[pat->has_AIDS][time_category].increment_category(pat);
		
		if (pat->cyclenum >= timepoint[time_category] && pat->cyclenum < (timepoint[time_category] + CYCLES_PER_YEAR)) { // In case they mean "during that year"
			alt_cumulative_years[pat->has_AIDS][time_category].increment_category(pat);
		}
	}
	
}

void process_WHO_discordance_output()
{
	FILE *who_discordance; //used for WHO_DISCORDANCE_ANALYSIS
	string headstring;
	headstring = "Year\tCD_100A: cumulative person-years lived with CD4 <100  and VL ³5000 , no stage 3 or 4 events\tCD_100B: cumulative person-years lived with CD4 <100 and VL < 5000, ³ 1000, no stage 3 or 4 events\tCD_100C: cumulative person-years lived with CD4 <100 and VL <1000, ³ 500, no stage 3 or 4 events\tCD_100D: cumulative person-years lived with CD4 <100 and VL <500, no stage 3 or 4 events\t";
	headstring += "CD_BLPA: cumulative person-years lived with CD4 <50% of peak or below baseline and VL ³5000 , no stage 3 or 4 events\tCD_BLPB: cumulative person-years lived with CD4 <50% of peak or below baseline and VL < 5000, ³ 1000, no stage 3 or 4 events\tCD_BLPC: cumulative person-years lived with CD4 <50% of peak or below baseline and VL <1000, ³ 500, no stage 3 or 4 events\tCD_BLPD: cumulative person-years lived with CD4 <50% of peak or below baseline and VL <500, no stage 3 or 4 events\t";
	headstring += "CD_100A_S4: cumulative person-years lived with CD4 <100  and VL ³5000 , with stage 4 event\tCD_100B_S4: cumulative person-years lived with CD4 <100 and VL < 5000, ³ 1000 , with stage 4 event\tCD_100C_S4: cumulative person-years lived with CD4 <100 and VL <1000, ³ 500 , with stage 4 event\tCD_100D_S4: cumulative person-years lived with CD4 <100 and VL < 500 , with stage 4 event\t";
	headstring += "CD_BLPA_S4: cumulative person-years lived with CD4 <50% of peak or below baseline and VL ³5000 , with stage 4 event\tCD_BLPB_S4: cumulative person-years lived with CD4 <50% of peak or below baseline and VL < 5000, ³ 1000 ,  with stage 4 event\tCD_BLPC_S4: cumulative person-years lived with CD4 <50% of peak or below baseline and VL <1000, ³ 500 , with stage 4 event\tCD_BLPD_S4: cumulative person-years lived with CD4 <50% of peak or below baseline and VL < 500 , with stage 4 event\t";
	
	
	string file_name = "WHO_Discordance"+IntToStr(scenario)+".txt"; // With 10 scenarios there will be 10 WHO_Discorcdance files.
	if( (who_discordance = fopen( file_name.c_str(), "w" )) == NULL )
	{
		printf("ERROR: The file %s was not opened\n", file_name.c_str());
		exit(0);
	}
	
	const string scenario_name[10] = {"B0", "B1", "B2", "B4", "B6", "B7", "B8", "A1", "A3", "A4"};
	fprintf(who_discordance, "Scenario %s\n",scenario_name[scenario].c_str());
	fprintf (who_discordance, "%s", headstring.c_str());
	fprintf(who_discordance, "\t\t");
	fprintf (who_discordance, "%s", headstring.c_str());
	fprintf(who_discordance, "\n");
	
	int years_in_category[4] = {0,5,10,15};
	// For each number of years output a row:
	for (int yr = 0; yr < 4; yr++)
	{
		fprintf(who_discordance, "%d\t",years_in_category[yr]);
		// Output non-stage 4, then stage 4 results:
		for (int aids_event = 0; aids_event < 2; aids_event++)
		{
			fprintf(who_discordance, "%f\t%f\t%f\t%f\t", (double)cumulative_years[aids_event][yr].cd4_L100_vl_G5000/CYCLES_PER_YEAR , (double)cumulative_years[aids_event][yr].cd4_L100_vl_L5000G1000/CYCLES_PER_YEAR, (double)cumulative_years[aids_event][yr].cd4_L100_vl_L1000G500/CYCLES_PER_YEAR, (double)cumulative_years[aids_event][yr].cd4_L100_vl_L500/CYCLES_PER_YEAR);
			fprintf(who_discordance, "%f\t%f\t%f\t%f\t",(double)cumulative_years[aids_event][yr].cd4_L50p_base_vl_G5000/CYCLES_PER_YEAR, (double)cumulative_years[aids_event][yr].cd4_L50p_base_vl_L5000G1000/CYCLES_PER_YEAR, (double)cumulative_years[aids_event][yr].cd4_L50p_vl_L1000G500/CYCLES_PER_YEAR, (double)cumulative_years[aids_event][yr].cd4_L50p_vl_L500/CYCLES_PER_YEAR);
		}
		fprintf(who_discordance, "\t\t%d\t",years_in_category[yr]);
		for (int aids_event = 0; aids_event < 2; aids_event++) {
			fprintf(who_discordance, "%f\t%f\t%f\t%f\t", (double)alt_cumulative_years[aids_event][yr].cd4_L100_vl_G5000/CYCLES_PER_YEAR , (double)alt_cumulative_years[aids_event][yr].cd4_L100_vl_L5000G1000/CYCLES_PER_YEAR, (double)alt_cumulative_years[aids_event][yr].cd4_L100_vl_L1000G500/CYCLES_PER_YEAR, (double)alt_cumulative_years[aids_event][yr].cd4_L100_vl_L500/CYCLES_PER_YEAR);
			fprintf(who_discordance, "%f\t%f\t%f\t%f\t",(double)alt_cumulative_years[aids_event][yr].cd4_L50p_base_vl_G5000/CYCLES_PER_YEAR, (double)alt_cumulative_years[aids_event][yr].cd4_L50p_base_vl_L5000G1000/CYCLES_PER_YEAR, (double)alt_cumulative_years[aids_event][yr].cd4_L50p_vl_L1000G500/CYCLES_PER_YEAR, (double)alt_cumulative_years[aids_event][yr].cd4_L50p_vl_L500/CYCLES_PER_YEAR);

		}
		fprintf(who_discordance, "\n");
	}

	fclose(who_discordance);
}

string IntToStr(int n)
{
	std::ostringstream result;
	result << n;
	return result.str();
}

#endif

void collectRyanWhiteStatsAtPointsInTime(patient *pat)
{
	int num_years;

	for (int time_idx = 0; time_idx < NUM_TIME_PERIODS_RYAN_WHITE; time_idx++)
	{
		switch (time_idx)
		{
		case 0: num_years = 2; break;
		case 1: num_years = 5; break;
		case 2: num_years = 10; break;
		case 3: num_years = 20; break;
		}

		if (pat->cyclenum < num_years*CYCLES_PER_YEAR)
		{
			costAtYear[time_idx] += costRyanWhiteInterv;		
			cyclesLivedAtYear[time_idx]++;
		}
	}
}

/*Puts a patient on haart and initializes all relevent parameters.  
This code used to be in Process_Patient but was put into a function so that we could easily call it if 
we want to start patients in the model already on treatment.  We need to do this now when we generate
rate files for the transmission model.
*/
void start_haart(patient *pat)
{
	pat->started_haart = true;
	pat->haart = true;
	
	//Assign patient to initial regimen.  If the user has specified a situation where there aren't
	//enough drugs available to make up the initial regimen, print a message and exit.
	if (notEnoughDrugsForInitialReg(pat))
	{
		cout<<"User specified impossible initial regimen.  Exiting."<<endl;
		fprintf (kim_test, "User specified impossible initial regimen.  Exiting.\n");
		exit(0);
	}
	else
	{
		pat->reg[0]=&drugs[pat->init_reg[0]][0];
		pat->reg[1]=&drugs[pat->init_reg[1]][0];
		pat->reg[2]=&drugs[pat->init_reg[2]][0];
		
		patRegCombo=0;			//see explanation of bestReg in hiv.h
	}
	
	//set start-of-reg values for CD4within
	pat->AgeRegBaseline = pat->age;
	pat->CD4RegBaseline = pat->CD4real;
	pat->VLatStartOfHaart = pat->VLreal;
	pat->CD4atStartOfHaart = pat->CD4real;
	
	if (KIMTEST) fprintf (kim_test, "Start the patient on haart!  Cycle: %d, luck2: %.2f, Regimen:  ", pat->cyclenum, pat->luck2);
	printDrugsInReg(pat);
	
	// set VL decrement based on initial regimen
	getVLdec(pat);
	if (KIMTEST) fprintf (kim_test, "VLdec is now %.1f\n", pat->VLdec);
	
	if (KIMTEST) fprintf(kim_test, "Res profile: (%d %d %d)  \n", pat->reg[0]->res, pat->reg[1]->res, pat->reg[2]->res);
	
	pat->haart_start_time = pat->cyclenum;
	
	setNextMonitorCycle(pat);					//figure out when to monitor CD4 or VL
	
	if (DO_CLINICAL_VISIT_EVERY_SIX_MONTHS_FOR_WHO)
	{
		incrementWHOArray(pat, numPatVisits);
		if (KIMTEST) fprintf (kim_test, "Clinic visit at cycle: %d\n", pat->cyclenum);
		sixMonthsSinceLastVisit = pat->haart_start_time + CYCLES_PER_YEAR/2;
	}
	
	//When starting with initial muts, it's possible that the patient was initially resistant to drugs in the initial reg.
	//If this is the case, immediately change the regimen.  Only do this for research rich, because
	//with resource poor, we don't have the resistance testing to know.
	if (RESOURCE_RICH && (pat->reg[0]->res || pat->reg[1]->res|| pat->reg[2]->res))
	{
		if (KIMTEST) fprintf (kim_test, "Patient resistant to drug in initial reg.  Changing reg.\n");
		
		Change_Regimen(pat);
	}
}
    
void graph_trajectories(patient *pat)
{
	if (pat->counterreg==1 && !pat->exhausted_clean_regs) eventIndicator=1000;
	else if (pat->counterreg==1 && pat->exhausted_clean_regs) eventIndicator=200;
	else if (pat->counterreg==2 && pat->numreg>1)
	{
		if (regChangeFailFlag) eventIndicator=-100;
		else if (regChangeIntolFlag) eventIndicator=-50;
		//reset reason-for-reg-change parameters
		regChangeFailFlag = false;
		regChangeIntolFlag = false;
	}
	else if (!pat->haart && pat->cyclenum == (pat->quit_haart + 1)) eventIndicator=100;
	else eventIndicator = 0;

#if RETENTION_IN_CARE

	if (outreach_interv_enabled)
	{
		if (ES1flag)
		{
			eventIndicator_RIC = 100;
			ES1flag = false;
		}
		else if (ES3flag)
		{
			eventIndicator_RIC = 300;
			ES3flag = false;
		}
		else if (ES4flag)
		{
			eventIndicator_RIC = 400;
			ES4flag = false;
		}
		else if (LTFUflag)
		{
			eventIndicator_RIC = -50;
			LTFUflag = false;
		}
		else if (startoutreachflag)
		{
			eventIndicator_RIC = -100;
			startoutreachflag = false;
		}
		else if (stopoutreachflag)
		{
			eventIndicator_RIC = -150;
			stopoutreachflag = false;
		}
		else if (LTFUfoundflag)
		{
			eventIndicator_RIC = -200;
			LTFUfoundflag = false;
		}
		else eventIndicator_RIC = 0;

		fprintf (kim_test2, "%f\t %f\t %f\t %d\t %d\n", pat->CD4RegBaseline, pat->CD4real, pat->VLreal, eventIndicator, eventIndicator_RIC);
	}
#else 
	fprintf (kim_test2, "%f\t %f\t %f\t %f\t %f\t %f\t %d \n", pat->CD4RegBaseline, pat->CD4regressionideal, pat->CD4real_before_noise, pat->CD4real, pat->HIVbaseline-((pat->haart)?pat->VLdec:0), pat->VLreal, eventIndicator);
#endif
}
        
        
        
        
#if RETENTION_IN_CARE
/* pat->LTFU: as long as pat is not in ES1, it's LTFU
pat->found: (1) pat goes back to ES1 (2)found by outreach intervention
pat->outreach: outreach process?
*/
void Check_Engagement (patient *pat)
{
	double random, random2;

	/************************** Assign transition probabilities ************************/
	/* prob_es1_to_es2 changes based on time in care*/
	if (pat->time_in_care <= 6*CYCLES_PER_MONTH)
	{
		prob_es1_to_es2 = PROB_ES1_TO_ES2_6MON;
	}
	else if (pat->time_in_care <= 12*CYCLES_PER_MONTH)
	{
		prob_es1_to_es2 = PROB_ES1_TO_ES2_12MON;
	}
	else
	{
		prob_es1_to_es2 = PROB_ES1_TO_ES2_12MON_MORE;
	}
	/* if not yet on ART, apply this multiplier */
	if(!pat->haart)
	{
		prob_es1_to_es2 *= pre_ART_mult_for_prob_es1_to_es2;
	}

	//If running the sensitivity analysis, we want to vary the prob_es1_to_es2
	prob_es1_to_es2 *= es1_to_es2_multiplier_for_sensitivity;

	if (risk_reduction_interv_enabled && (pat->CD4real <= risk_reduction_cd4) ) prob_es1_to_es2 *= risk_reduction_rr_prob_es1_to_es2;

	if (secondary_prevention_interv_enabled && (pat->num_cycles_ES4 >0) && (pat->CD4real <= secondary_prevention_cd4)) 
		prob_es1_to_es2 *= secondary_prevention_rr_prob_es1_to_es2;

	prob_es3_to_es4 = prob_es1_to_es2 * prob_es2_to_es4;
	/************************** Finish Assignning transition probabilities ************************/

	random = (uniform(&seed));

	if (pat->haart == true)       //accumulate num of cycles on ART
	{
		pat->num_cycles_ART++;
		total_time_onART++;
	}
	if (pat->engaged_in_care && pat->engaged_in_clinic)  // engagement status ES1: pat is in care and in clinic
	{
		pat->num_cycles_ES1++;
		total_ES1_time++;
		// engagement status ES2: pat is not in clinic (transition state)
		if (random < prob_es1_to_es2)
		{
			pat->LTFU = true;                    // patient is lost-follow-up
			LTFUflag = true;					 //marks the day patient became LTFU
			pat->num_LTFU++;					 //increment the number of times this patient is LTFU
			total_LTFU++;						 //increment the total number of times all patients are LTFU
			pat->cyclenum_LTFU = pat->cyclenum;  //ccycle num pat last showed up in clinic
			pat->engaged_in_clinic = false;      // lost clinic engagement

			if (outreach_interv_enabled) pat->found = false;    //at the time point of LTFU, impossible found immediately

			random2 =(uniform(&seed));
			if (random2 <= prob_es2_to_es4)
			{
				ES4flag = true;
				pat->engaged_in_care = false;    //go to ES4, else will go to ES3
				pat->haart = false;
				pat->time_in_care = 0;                //disengaged from care, reset time-in-care to 0
			}
			else
			{
				ES3flag = true;
			}

			if (KIMTEST)
			{
				fprintf (kim_test, "Patient is LTFU at cycle %d.\n", pat->cyclenum);
				if (ES3flag) fprintf (kim_test, "Patient is ES3 (in care in another clinic).\n");
				if (ES4flag) fprintf (kim_test, "Patient is ES4 (care disengaged).\n");
			}
        }
	}

	if (outreach_interv_enabled)
	{
		if (!(pat->engaged_in_care && pat->engaged_in_clinic))                // if patient is not on state ES1, start to outreach patient if triggers are met.
		{
			if ((pat->LTFU)&&(!pat->outreach)&&(!pat->found)&&((pat->cyclenum - pat->cyclenum_LTFU) == outreach_start_trigger) && (pat->CD4real <= outreach_cd4))
			{
				pat->outreach = true;   // start outreaching
				startoutreachflag = true;
				total_num_outreach++;   // increase total amount of outreaching over all population by 1
				if (KIMTEST) fprintf (kim_test, "Initiating Outreach at cycle %d.  Incrementing num outreach to %d.\n", pat->cyclenum, total_num_outreach);


				total_outreach_cost += outreach_initiation_cost_poor;
					
				//add the outreach cost to the total cost
				total_cost += outreach_initiation_cost_poor;
				total_disc_cost += outreach_initiation_cost_poor / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
			}
			if (pat->outreach)
			{
				if ((pat->cyclenum-pat->cyclenum_LTFU) >= (outreach_start_trigger+outreach_end_trigger)) //if we still cannot find this pat when time limit is reached, stop outreach.
				{
					pat->outreach = false;
					stopoutreachflag = true;
					total_num_stopoutreach++;
					if (KIMTEST) fprintf (kim_test, "Stopping outreach at cycle %d.\n", pat->cyclenum);

					total_outreach_cost += outreach_finishing_cost_poor;

					//add the outreach cost to the total cost
					total_cost += outreach_finishing_cost_poor;
					total_disc_cost += outreach_finishing_cost_poor / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
				}
				else if (pat->found == true)
				{
					if ((pat->engaged_in_care && !pat->engaged_in_clinic)||(uniform(&seed) < outreach_prob_relink)) //if patient is found on ES3, no relinkage, back to ES1 directly. Otherwise, go through relinkage process.
					{
						ES1flag = true;                        //back to ES1
						pat->engaged_in_care = true;           
						pat->engaged_in_clinic = true;
						pat->LTFU = false;                     //not LTFU anymore
						pat->outreach = false;                 //outreach process finished
						stopoutreachflag = true;
						total_num_stopoutreach++;
						if (!pat->haart && pat->started_haart)
						{
							pat->haart = true;
							setNextMonitorCycle(pat);
						}
					}
				}
				else if (pat->found == false)
				{
					if (uniform(&seed) < outreach_prob_find) // Probability of finding. If cannot find, keep finding until meet end trigger
					{
						pat->found = true;     // found this pat (1, relink-not LTFU anymore; 2, not relink, found but still LTFU)
						LTFUfoundflag = true;   // indicates found of this pat
						if (KIMTEST) fprintf (kim_test, "Found patient at cycle %d.", pat->cyclenum);

					}
				}
			}
		}
	}	//ends if (outreach_interv_enabled)

	if (pat->engaged_in_care && !pat->engaged_in_clinic) // ES 3: pat is engaged in care but not engaged in clinic
	{
		pat->num_cycles_ES3++;
		total_ES3_time++;

		if (uniform(&seed) < prob_es3_to_es4) // pat will flow from ES3 to ES4, and be disengaged from care (ES4)
		{
			ES4flag = true;
			pat->engaged_in_care = false;
			pat->haart = false;
			pat->time_in_care = 0; // reset time-in-care = 0 
		}
	}
	if ((!pat->engaged_in_care) && (!pat->engaged_in_clinic)) // ES 4: pat is not in care (neither in care nor in clinic)
	{
		pat->num_cycles_ES4++;
		total_ES4_time++;
	}
}

void check_for_AIDS_defining_events_RIC (patient *pat)
{
	//check for AIDS every cycle even if patient already has AIDS
	//The rate of an AIDS-defining event is AIDS_RATE_MULTIPLIER * pat->deathrateHIV
	pat->pAIDS =  1 - exp(-(AIDS_EVENT_MULTIPLIER * pat->deathrateHIV) * CYCTIME);

	if (uniform(&seed) <= pat->pAIDS)
	{
		if (!pat->has_AIDS)
		{
			pat->has_AIDS = true;
			if (DEBUG_AIDS && (int)(pat->cyclenum/CYCLES_PER_YEAR) <= 5)
				totalAIDSOrDeath++;
		}

		//Patient has had an AIDS-defining event during the current regimen
		pat->AIDSeventReg = true;
		if (KIMTEST) fprintf (kim_test, "\nPatient had AIDS-defining event at cycle: %d\n\n", pat->cyclenum);

		//Patient has a new AIDS-defining event and disengaged from care goes back to care.
		if ((!pat->engaged_in_care) && (!pat->engaged_in_clinic)) // ES 4: pat is not in care (neither in care nor in clinic)
		{
			if (uniform(&seed) < prob_es4_to_es1_AIDSpat)
			{
				ES1flag =true;
				pat->engaged_in_care = true;
				pat->engaged_in_clinic = true;
				pat->LTFU = false;

				if (outreach_interv_enabled)
				{
					if (!pat->found)
					{
						pat->found = true;
						LTFUfoundflag = true; // indicates this pat is found
					}
					if (pat->outreach) // pat is back by himself, stop outreaching if already in process
					{
						pat->outreach = false;
						stopoutreachflag = true;
						total_num_stopoutreach++;
					}
				}

				if (pat->started_haart)
				{
					pat->haart = true;
					setNextMonitorCycle(pat);
				}
			}
			else
			{
				ES3flag = true;
				pat->engaged_in_care =true;
				if (pat->started_haart)
				{
					pat->haart = true;
					setNextMonitorCycle(pat);
				}
			}
		}
	}
}

void process_HIV_death_RIC(patient *pat)
{
	double random;

	//Even if the patient in engaged in care and clinic, there is a possibility
	//that their death will be unknown
	if(pat->engaged_in_care&&pat->engaged_in_clinic) //ES1
	{
		random = (uniform(&seed));
		if (random < prob_es5_to_es6 )
		{
			total_ES1_HIV_death_known++;
		}
		else
		{
			total_ES1_HIV_death_unknown++;
			total_LTFU++;
			if (KIMTEST) fprintf (kim_test, "HIV death unknown.  Incrementing total_LTFU to %d.\n", total_LTFU);

			if (outreach_interv_enabled && (pat->CD4real <= outreach_cd4))  //death unknown on state 1 consumes entire outreach effort from initiation to finishing
			{
				if (pat->outreach == true)
				{
					cout<<"ERROR: it's impossible for ES1 patients to be in outreach process"<<endl;
				}
				pat->outreach=true;  //if there's an unknown death, there is an outreach effort. Here turn this on in order to calculate total_outreach_finding_cost.
				total_num_outreach++;

				if (KIMTEST) fprintf (kim_test, "Initiating Outreach at HIV death at cycle %d.  Incrementing num outreach to %d.\n", pat->cyclenum, total_num_outreach);
			}

		}
	}
	else if (pat->engaged_in_care&&!pat->engaged_in_clinic) //ES3
	{
		total_ES3_HIV_death++;
	}
	else if (!pat->engaged_in_care&&!pat->engaged_in_clinic) //ES4
	{
		total_ES4_HIV_death++;
	}
}

void process_age_death_RIC(patient *pat)
{
	double random;

	//Even if the patient in engaged in care and clinic, there is a possibility
	//that their death will be unknown
	if(pat->engaged_in_care&&pat->engaged_in_clinic) //ES1
	{
		random = (uniform(&seed));
		if (random < prob_es5_to_es6 )
		{
			total_ES1_AGE_death_known++;
		}
		else
		{
			total_ES1_AGE_death_unknown++;
			total_LTFU++;
			if (KIMTEST) fprintf (kim_test, "Age death unknown.  Incrementing total_LTFU to %d.\n", total_LTFU);

			if (outreach_interv_enabled && (pat->CD4real < outreach_cd4)) //death unknown on state 1 consumes entire outreach effort from initiation to finishing
			{
				if (pat->outreach == true)
				{
					cout<<"ERROR: it's impossible for ES1 patients in outreach process"<<endl;
				}
				pat->outreach=true;  //if there's an unknown death, there is an outreach effort. Here turn this on in order to calculate total_outreach_finding_cost.
				total_num_outreach++;
				
				if (KIMTEST) fprintf (kim_test, "Initiating Outreach at Age Death at cycle %d.  Incrementing num outreach to %d.  ", pat->cyclenum, total_num_outreach);
			}
		}
	}
	else if (pat->engaged_in_care&&!pat->engaged_in_clinic) //ES3
	{
		total_ES3_AGE_death++;
	}
	else if (!pat->engaged_in_care&&!pat->engaged_in_clinic) //ES4
	{
		total_ES4_AGE_death++;
	}
}

//If the patient dies and their death is unknown, we add half the cost of normal outreach.
void addOutreachCostsToTotalCostAtDeath(patient *pat)
{
	double outreachCost;
	
	outreachCost = (outreach_finding_cost_poor * outreach_end_trigger/2) 
		+ outreach_initiation_cost_poor; // assume we spend half of maximum outreach duration (outreach_end_trigger) for tracking unknown death.  No finishing cost because patient isn't relinked.
	
	//Add this cost to the total outreach cost
	total_outreach_cost += outreachCost;

	//add the costs of outreach to the total cost
	total_cost += outreachCost;
	total_disc_cost += outreachCost / pow( (1+DISC_RATE), ((double)pat->cyclenum/(double)CYCLES_PER_YEAR));
}

#endif


#ifdef ONE_WAY_SENSITIVITY_ANALYSIS_RIC
int ReadSensitivityFile_RIC ()
{
	static bool moreLinesToRead = true;
	static bool firstRead = true;
	string headings;
	static ifstream infile;

#ifdef USING_COMMAND_LINE_ARGS
	static int lineNum = 0;					//used for running Sensitivity on cluster
#endif

	//For now, to hold the two cost parameters

	if (firstRead)
	{
		infile.open ("OneWaySensitivity_RIC.txt");

		if (infile.fail())
		{
			cout << "Error opening OneWaySensitivity_RIC.txt" << endl;
			exit(1);
		}
		getline(infile, headings);

		senseFile << "\tOutreach Intervention On\tRisk Reduction Intervention On\tSecondary Prevention Intervention On\tInitial CD4 (mean)\tInitial age (mean) \t Comp \t CD4 treat\t Cost of 1st line ART\t Cost of 2nd line ART\tCost of outpatient care annually\tHospitalization Cost\tUnderlying risk of non-retention multiplier (rate ES1->ES2)\tPreART non-retention multiplier\tProbability of disengagement from care given disengagement from clinic (ES4 given ES2)\tProbability of ascertainment of death (ES6 given ES5)\tMultiplier on mutation rate if in ES4 > 1 month and relink to care\tProbability of Identification of disengaged from clinic (pI)\tProbability of finding following identification (pF)\tProbability of re-linkage following successful tracking (pL)\tStart “trigger” for outreach\tStop “trigger” for outreach\tOutreach Identification Cost\tOutreach Tracking Cost\tOutreach Relinking Cost\tOutreach CD4 Threshold\tRisk Reduction rr on ES1->ES2 \tRisk Reduction Interv Cost\t Risk Reduction CD4 Threshold \t Secondary Prevention rr on ES1->ES2 \t Secondary Prevention Interv cost \t sp CD4 threshold \t prob ART init at CD4 treat \t Patients \t Total ES1 HIV death known \t Total ES1 HIV death unknown \t total ES1 AGE death known \t total ES1 AGE death unknown \t total ES3 HIV death \t total ES3 AGE death \t total ES4 HIV death \t total ES4 AGE death \t mean ES1 time (years) \t mean ES3 time (years) \t mean ES4 time (years) \t mean time on ART (years) \t mean surv time (years) \t mean disc surv time (years) \t median surv time (years) \t mean QALYs\tmean discounted QALYs\tmedian qalys \tmean time VL<1000 (years) \t mean total cost \t mean disc total cost \t mean lab cost \t mean drug cost \t mean care cost \t mean hospital cost \t mean outreach cost (per patient)\t mean RR interv cost \t mean SP interv cost \ttotal num times outreach initiated \t mean num times outreach initiated per patient \t total LTFU events \t mean num LTFU events per patient" <<endl<<endl;

		firstRead = false;
	}

#ifdef USING_COMMAND_LINE_ARGS
	else if (lineNum == RIC_command_line_sensitivity_row_num) return false;

	while (lineNum<RIC_command_line_sensitivity_row_num)
	{
#endif
		//Start reading first line
		infile >> outreach_interv_enabled;


		//If we're at the end of the file, return false;
		if (infile.eof()) return false;

		//If we're here, the must be a line, so continue reading it!
		infile >> risk_reduction_interv_enabled;
		infile >> secondary_prevention_interv_enabled;
		infile >> start_cd4;
		infile >> start_age;
		infile >> comp;
		infile >> cd4_treat;
		infile >> cost_reg1_poor;
		infile >> cost_reg2_poor;
		infile >> cost_of_care_poor;
		infile >> hospital_cost_poor;
		infile >> es1_to_es2_multiplier_for_sensitivity;
		infile >> pre_ART_mult_for_prob_es1_to_es2;
		infile >> prob_es2_to_es4;
		infile >> prob_es5_to_es6;
		infile >> mutation_mult_es4;

		//for OUTREACH_INTERVENTION
		infile >> outreach_prob_identify;
		infile >> outreach_prob_find;
		infile >> outreach_prob_relink;
		infile >> outreach_start_trigger;
		infile >> outreach_end_trigger;
		infile >> outreach_initiation_cost_poor;
		infile >> outreach_finding_cost_poor;
		infile >> outreach_finishing_cost_poor;  
		infile >> outreach_cd4;

		//for RISK_REDUCTION_INTERVENTION
		infile >> risk_reduction_rr_prob_es1_to_es2;
		infile >> risk_reduction_intervention_cost;
		infile >> risk_reduction_cd4;

		//for SECONDARY_PREVENTION_INTERVENTION
		infile >> secondary_prevention_rr_prob_es1_to_es2;
		infile >> secondary_prevention_intervention_cost;
		infile >> secondary_prevention_cd4;

		//the probability that each patient will start on treatment at the specified threshold instead of 200.
		infile >> prob_ART_init_at_CD4_treat;

		getline (infile, label);

		if (!label.empty()) label = label.substr(1,label.size());	//Remove the leading tabs

#ifdef USING_COMMAND_LINE_ARGS
		lineNum++;
	}
#endif

	return true;
}
#endif  //end ONE_WAY_SENSITIVITY_ANALYSIS_RIC




