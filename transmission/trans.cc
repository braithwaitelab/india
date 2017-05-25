#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;		//required now for cout to work!
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <algorithm>			//required now for min/max to work
#include "constants.h"
#include "trans.h"
#include "time.h"
#include "intervention.h"
#include "distributions.h"
#pragma warning(disable:4127)

void initialize();
void setup();
double convertToTimeStep(double rate);
void Runge_Kutta();
void Diff(double Pop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES], double &newInfections, double &HIV_DeathRate, double &newInfectionsAlc, double &numCircumsizedThisCycle, double &newInfectionsIDU);
void outputData();
void calculateBetaandLambda(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES]);
void balancePartnerships();
void initializeArrays(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES]);
void calculateBalanceMatrix(double c_matrix[NUM_SEX][NUM_ACT][NUM_AGE][NUM_ACT][NUM_AGE], bool errorFlag);
inline bool isEqual(double x, double y);
void determineProbabilityMatrix();
void calculateAgingBirthsDeaths(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES], double &numInfectedBirths, double &numCircumsizedThisCycle);
void printMatrix(double matrix[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE]);
void readHIVprogrates(string &off_care_path, string &on_care_root);
void printNonZeroCompartments(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES]);
void calculateCD4andVLratesNewInfecteds();
double normal_cdf(double x);
void calculateCompromiseConcurrencies(double Pop[NUM_STATUS][NUM_SEX][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES]);
double compromise (double a, double b);
void set_age_change_rates(int age_of_sexual_debut);
void scaleUpFertilityRates();
void createBasketOutput();
void setBaselineProbs();
double convertProbabilityToRate(double prob);
void initialize_interventions();
double getBaselineRRproduct(int p, int status, int sex, int orient, int act, int alc, int idu);
double getBasketCost(int utilization_or_prepurchased, bool discounted_or_normal);
double getTotalDeadAIDS();
void calcPropRes(double Pop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES]);
void setBaseLineProbsTestedLinkedandAdherence();
double linearInterpolate(double x0, double y0, double x1, double y1, double x);
double getAdhStratifiedHIVProgRate(int status, int curr_cd4, int curr_vl, int curr_res, int dest_cd4, int dest_vl,  int destOnTreatment, int destRes, int dead, int sex, int orient, int act, int alc, int idu);
double getAdhStratifiedCostCareTreat(int status, int curr_cd4, int curr_vl, int curr_res, int sex, int orient, int act, int alc, int idu);
string IntToStr(int n);
string DoubleToStr(double n);
void show_rates();
void setTotalPeopleinModel(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES]);
void calculateCost();
void read_population();
void save_population();
void rampUpValuesDuringCalibration();
void printNumInEachStatus();
void processInterventionsThatMovePeople();
double getAlphaMult(int sex, int orient, int act, int alc, int idu, int vl, int p_status, int p_sex, int p_orient, int p_act, int p_alc, int p_idu, int p_vl);
void setPropOrientActAlcIDU();
void calculate_propIDUinAct(); 

void change_prevalence(bool high);
void movehalf(int sx);
void vary_interv_effectsizes(int run);
void vary_baseline_probabilities(int run);
void vary_interv_additive(int num);
void vary_adherence(int num);
void vary_epsilon_gamma(int num);
void backup_arrays();
void restore_arrays();
void vary_yintercept(int num);
void vary_proportion_alcoholic(bool high);
void test_alc_baselines(int num);
void vary_activity_classesFREQ(bool high);
void vary_durations(bool high);
void vary_concurrencies(bool high);
void chris_sensitivity(int onewayrun);
void postProcessBaselineProbsForAlcAnalysis();
void postProcessBaselineProbsForOptimizingAlcoholIntervention();
void setAlcoholicProportionsAtStartOfRun();
void choose_rate_files();
double getProportionOfWomenWhoGiveBirthInAYear();
void setPropWomenOnPMTCT();
void LifeYears_and_QALYS();
double Utility(int status, int cd4);
double get_discounted_value(double val_to_discount);
void initializeVarsForAlcoholSurfaceAnalysis();
void preCalculateAdhStratHIVProgRateAndCostArrays();	//EFFICIENCY!
void precalc_getAlphaMult();							//EFFICIENCY!
void check_HIVprogrates();								//Error checking of rate files
bool rate_is_zero(int curr_cd4, int curr_vl, int curr_onTreat, int curr_resistance, int adh);  //Error checking of rate files
void findNumMenCircumsizedByIntervAtStartOfRun();
void postProcessBaselineProbsForOptimizingPreventionPorfolio();
void alcoholIntervention();
void initialize_willPartnerWithMatrix();
void print_time();
void preventBalanceError();
double lognormal(double mu, double lower, double upper, long int seed, int max_value);
void zero_all_baskets();
void do_loop(int loops_remaining);

int main(int argc, char * const argv[])
{
	print_time();			//Print the time the run starts
	
	setup();
	atLeastOneInterventionTurnedOn = false;

#if RUN_ALL_USING_SCRIPT
	// Running each intervention separately using baseline effect sizes and costs.

	if (RUN_NUM != 999)
	{
		basket[RUN_NUM] = 1;
	}
#endif

#if OPTIMIZING_PREVENTION_PORTFOLIO
    
	int run_num = RUN_NUM, i;
    
	// Turn all interventions off.
	for (i = 0; i < NUM_INTERVENTIONS-9; i++)  basket[i+9] = 0;
	
	// Read in the run number from the command line (if using command line arguments)
	// Note, we usually AREN'T.  We're usually using RUN_NUM from constants.h
	if (argc > 1)
	{
		run_num = atoi(argv[1]);
		cout << "run_num: "<< run_num << endl;
	}
    
    
	//Chris' bit shifting magic.  Use the bits in how the value of the run number is stored to set
	//whether each intervention should be turned on or off.
	for (i = 0; i < NUM_INTERVENTIONS-9; i++)
	{
		if (run_num % 2 == 1)
		{
			basket[i+9] = true;
			if (run_in_optimization[i + 9] == 0) //check if the intervention is included in the analysis
			{
				cout << "Intervention number " << i+9 << " not included in analysis.  Exiting." << endl;
				exit(0);
			}
		}
		run_num = run_num>>1;
	}
#endif
    
#if COST_SENSITIVITY
	// Running each intervention using baseline effect sizes, but varying costs between high and low values.

	int intervention_number = RUN_NUM/2;// Each intervention number has a high and low, so there are two RUN_NUM values for each
	basket[intervention_number] = 1;
	if(RUN_NUM % 2 > 0)
	{
		cout << "Testing intevention "<< intervention_number << " HIGH\n";
		high_low = "High "+IntToStr(intervention_number);
		interv_costs_per_person[intervention_number] = interv_costs_per_person_high[intervention_number];
	} else
	{
		cout << "Testing intevention "<< intervention_number << " LOW\n";
		high_low = "Low "+IntToStr(intervention_number);
		interv_costs_per_person[intervention_number] = interv_costs_per_person_low[intervention_number];
	}
#endif

#if INTERVENTION_EFFECT_SIZE_SENSITIVITY
	// Running each intervention separately using high and then low values for effect sizes.

	int interv_sense = RUN_NUM;

	if (argc > 1)
	{
		interv_sense = atoi(argv[1]);
		cout << "interv_sense: "<< interv_sense << endl;
	}

	int intervention_number = interv_sense/2;// Each intervention number has a high and low, so there are two RUN_NUM values for each
	basket[intervention_number] = 1;
	if(interv_sense % 2 > 0)
	{
		cout << "Testing intevention "<< intervention_number << " HIGH\n";
		high_low = "High "+IntToStr(intervention_number)+" path(s): ";
		for (int path_num = 0; path_num < NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE; path_num++)
		{
			interv_effect_size[intervention_number][path_num] = interv_effect_size_high[intervention_number][path_num];			
			if (interv_effect_size[intervention_number][path_num] != 1) {
				high_low += "P " + IntToStr(path_num) + ": " + DoubleToStr(interv_effect_size[intervention_number][path_num])+" ";
			}
		}
	} else
	{
		cout << "Testing intevention "<< intervention_number << " LOW\n";
		high_low = "Low "+IntToStr(intervention_number)+" path(s): ";
		for (int path_num = 0; path_num < NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE; path_num++)
		{
			interv_effect_size[intervention_number][path_num] = interv_effect_size_low[intervention_number][path_num];
			if (interv_effect_size[intervention_number][path_num] != 1) {
				high_low += "P " + IntToStr(path_num) + ": " + DoubleToStr(interv_effect_size[intervention_number][path_num])+" ";
			}
		}
		cout << endl;;
	}

#endif

#if OPTIMIZING_ALCOHOL_INTERVENTION	

	int condomRR_index, adhereRR_index, stiRR_index, prop_alc_index, effect_size_index;

	double RR_condoms[2] = {baselineRR_raw[0][5], 10};
	double RR_nonadherence[2] = {baselineRR_raw[3][5], 10};
	double RR_STIs[2] = {baselineRR_raw[5][5], 10};

	double prop_alcoholic_M[3] = {PROP_ALCOHOL_M, (PROP_ALCOHOL_M + 0.35) / 2.0, 0.35};	//Baseline, midpoint of baseline and .35, .35
	double prop_alcoholic_F[3] = {PROP_ALCOHOL_F, (PROP_ALCOHOL_F + 0.35) / 2.0, 0.35};	//Baseline, midpoint of baseline and .35, .35

	double alc_effect_size[3] = {1, 0.55, 0.1};

	for (int i = 0; i < NUM_INTERVENTIONS; i++)
	{
		basket[i] = 0; // Turn all interventions off.
	}

	basket[3] = 1; // Turn on intervention 3: alcohol issues.	

	//These loops are commented out and constants used above so that we can break this up to run on 24 nodes on the cluster. See script: optalcohol.sh

	/*
	for (condomRR_index = 0; condomRR_index < 2; condomRR_index++)
	{
	for (adhereRR_index = 0; adhereRR_index < 2; adhereRR_index++)
	{
	for (stiRR_index = 0; stiRR_index < 2; stiRR_index++)
	{
	for (prop_alc_index = 0; prop_alc_index < 3; prop_alc_index++)
	{	
	*/
	for (effect_size_index = 0; effect_size_index < 3; effect_size_index++)
	{
		if (alc_effect_size[effect_size_index] == 1)
		{
			is_baseline_run = true;
		} 
		else
		{
			is_baseline_run = false;
		}

		//Commented out for use with script
		/*
		prop_alcohol[M] = prop_alcoholic_M[prop_alc_index];
		prop_alcohol[F] = prop_alcoholic_F[prop_alc_index];

		condomRR = RR_condoms[condomRR_index];
		adhereRR = RR_nonadherence[adhereRR_index];
		stiRR = RR_STIs[stiRR_index];

		*/

		prop_alcohol[M] = prop_alcoholic_M[ALCOHOL_VALUE];
		prop_alcohol[F] = prop_alcoholic_F[ALCOHOL_VALUE];

		condomRR = RR_condoms[CONDOM_VALUE];
		adhereRR = RR_nonadherence[ADHERENCE_VALUE];
		stiRR = RR_STIs[STI_VALUE];

		interv_effect_size[3][10] = alc_effect_size[effect_size_index];
#endif 
		
			if (RUNNING_MULTIPLE_BASKETS) basket_num = -1;  //initialize to -1 to calculate with no interventions
			else basket_num = 0;


#ifdef RUN_ALL_NO_SCRIPT

			for (int interv_num=0; interv_num < NUM_INTERVENTIONS; interv_num++)
			{
				for (int i=0; i<NUM_INTERVENTIONS; i++) basket[i]=0;

				//Turn on the next intervention
				if (basket_num != -1) basket[interv_num]=1;
#endif

#if ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION

				baselineRR_raw_backup_condoms = baselineRR_raw[0][5];			//the alcohol RR on the condom pathway
				baselineRR_raw_backup_nonadherence = baselineRR_raw[3][5];		//the alcohol RR on the nonadherence pathway
				baselineRR_raw_backup_STIs = baselineRR_raw[5][5];				//the alcohol RR on the STI pathway

				run = RUN_NUM;   // Before scripted for cluster was this:  for (run=0; run<10; run++) {

				//initialize all relevent relative risks to 1 for alcohol compartments
				baselineRR_raw[0][5] = 1;		//condoms
				baselineRR_raw[3][5] = 1;		//adherence
				baselineRR_raw[5][5] = 1;		//STI's

				// Initialize all interventions to off.
				for (int bsk = 0; bsk < NUM_INTERVENTIONS; bsk++) basket[bsk] = 0;

				//Turn on the alcohol intervention for the even numbered runs and off for the odd numbered runs
				if (run%2==0) basket[3]=1;
				else basket[3]=0;
#endif

#ifdef TESTING_LINKAGE_ADHERENCE

				double baselineValForLater[3] = {baselineRR_raw[1][0], baselineRR_raw[2][0], baselineRR_raw[3][0]}; 

				for (int scenario = 0; scenario <=3; scenario++)
				{
					for (int baseline_prob_magnitude = -1; baseline_prob_magnitude <= 1; baseline_prob_magnitude ++)
					{
						num_infections_no_interv = 0;			//Since we aren't running interventions, this variable gets reused to hold the number of new infections.
						
						switch (baseline_prob_magnitude) 
						{
						case -1:	//testing, linkage and adherence baseline probs set to baselines
							baselineRR_raw[1][0] = baselineValForLater[0];	//testing	
							baselineRR_raw[2][0] = baselineValForLater[1];	//linkage
							baselineRR_raw[3][0] = baselineValForLater[2];	//adherence
							break;
						case 0:		//testing, linkage and adherence baseline probs set to 0s
							baselineRR_raw[1][0]= 0;	//testing	
							baselineRR_raw[2][0]= 0;	//linkage
							baselineRR_raw[3][0]= 0;	//adherence
							break; 
						case 1:		//testing, linkage and adherence baseline probs set to 1s
							baselineRR_raw[1][0]= 1;	//testing
							baselineRR_raw[2][0]= 1;	//linkage
							baselineRR_raw[3][0]= 1;	//adherence
							break; 
						}

						switch (scenario) 
						{
						case 0:	
							turn_Off_HIV_Mortality = true;
							turn_Off_Transmission_On_Treatment = false;
							cd4_treat_threshold = 10000;
							break;
						case 1:	
							turn_Off_HIV_Mortality = false;
							turn_Off_Transmission_On_Treatment = false;
							cd4_treat_threshold = 10000;
							break;
						case 2:	
							turn_Off_HIV_Mortality = false;
							turn_Off_Transmission_On_Treatment = true;
							cd4_treat_threshold = 10000;
							break;
						case 3:	
							turn_Off_HIV_Mortality = false;
							turn_Off_Transmission_On_Treatment = false;
							cd4_treat_threshold = 200;
							break;
						}
#endif

#if ALCOHOL_SURFACE
					//Get the indices into the arrays for the cost and effect sizes of the alcohol intervention
					//as well as for the RR modifier on the pathways affected by alcohol use
					
					if (argc == 4) // Using command line arguments
					{
						cost_index = atoi(argv[1]);
						effect_index = atoi(argv[2]);
						rr_index = atoi(argv[3]);
					} 
					else if(argc == 1) // Using constants
					{
						cost_index = SURF_COST;
						effect_index = SURF_EFFECT;
						rr_index = SURF_RR;
					} 
					else // Something's wrong!
					{ 
						cout << "argc is: "<< argc << " and we can only use 1 or 4. Exiting.\n";
						exit(0);							
					}

					initializeVarsForAlcoholSurfaceAnalysis();
#endif

#if CARE_PLUS
	// Turn all interventions off.
	for (int i = 0; i < NUM_INTERVENTIONS; i++)	basket[i] = 0;

	// Turn on intervention 7: Care+.
	basket[7] = 1;
#endif

#if COST_EFFECTIVENESS_SENSITIVITY || PROBABILISTIC_RUN
  

	if (argc == 2)  //the name of the model plus the run_num value
	{
		// Read in the run number from the command line
		run_num = atoi(argv[1]);

		//Set the initial seed so that each run will have a different
		//value of the seed
		seed = (long int)run_num;

		cout << "seed = " << seed << endl;
	}

	else // Something's wrong!
	{
		cout << "argc is: " << argc << " and it should be 2. Exiting.\n";
		exit(0);
	}
#endif
	
#if COST_EFFECTIVENESS_SENSITIVITY || PROBABILISTIC_RUN
	init_genrand(seed); //sets the initial seed

	/*Find the effectiveness and cost for each of the interventions we are testing using a distribution
	Use lognormal for both cost and relative risk based on Briggs, Claxton, Sculpher 
	Decision Modelling for Health Economic Evalulation page 108*/
	
	interv_effect_size[10][10] = lognormal(.36, .15, .82, seed, 1); //0.36 range: 0.15-0.82
	interv_effect_size[19][0] = lognormal(.92, .85, .98, seed, 1); //0.92 range: 0.85-0.98
	interv_effect_size[19][5] = lognormal(.64, .44, .89, seed, 1); //0.64 range: 0.44-0.89
	interv_effect_size[20][0] = lognormal(.97, .94, .99, seed, 1); //0.97 range: 0.94-0.99
	interv_effect_size[20][5] = lognormal(.81, .68, .95, seed, 1); //0.81 range: 0.68-0.95
	interv_effect_size[21][0] = lognormal(.94, .9, .99, seed, 1); //0.94 range: 0.90-0.99
	interv_effect_size[21][5] = lognormal(.71, .54, .92, seed, 1); //0.71 range: 0.54-0.92
	interv_effect_size[22][0] = lognormal(.97, .94, .99, seed, 1); //0.97 range 0.94-0.99
	interv_effect_size[22][5] = lognormal(.78, .59, 1.0, seed, 1); //0.78 range: 0.59-1
	interv_effect_size[23][3] = lognormal(.67, .53, .84, seed, 1); //0.67 range: 0.53-0.84
	interv_effect_size[24][3] = lognormal(.75, .58, .96, seed, 1); //0.75 range: 0.58-0.96
    


	//cost range between 0.5 - 1.5 of the mean
	interv_costs_per_person[10] = uniform_a_b(3.28, 9.84, &seed); //6.56 range: 3.28-9.84
	interv_costs_per_person[19] = uniform_a_b(7.38, 22.14, &seed); //14.76 range: 7.38-22.14
	interv_costs_per_person[20] = uniform_a_b(0.82, 2.46, &seed); //1.64 range: 0.82-2.46
	interv_costs_per_person[21] = uniform_a_b(2.952, 8.856, &seed); //5.904 range: 2.952-8.856
	interv_costs_per_person[22] = uniform_a_b(3.335, 10.005, &seed); //6.67 range: 3.335-10.005
	interv_costs_per_person[23] = uniform_a_b(3.28, 9.84, &seed); //6.56 range: 3.28-9.84
	interv_costs_per_person[24] = uniform_a_b(1.23, 3.69, &seed); //2.46 range: 1.23-3.69

	

	/**************find the basket file we are using based on an input **************/ 
	zero_all_baskets();
	//basket_count=0;

		//Make sure the user set the number of baskets to 2^NUM_INTERV_TO_COMBINE
		if (TOTAL_NUM_BASKETS != pow(2.0, NUM_INTERV_TO_COMBINE))
		{
			cout << "Adjust TOTAL_NUM_BASKETS!!!  Should be " << pow(2.0, NUM_INTERV_TO_COMBINE) << "." << endl;
			exit(0);
		}

		//basket_count and loops_remaining must be initialized as follows before calling do_loop
		//Because do_loop is recursive, these can't be initialized inside the function
		basket_count = 0; 
		basket_use = BASKET_USE;
		loops_remaining = NUM_INTERV_TO_COMBINE;

		//The function do_loop magically fills in the array_of_baskets with the 1s and 0s for each basket
		do_loop(loops_remaining - 1);

	//just for debugging, print out all baskets
	/*
	for (basket_count = 0; basket_count<TOTAL_NUM_BASKETS; basket_count++)
	{
		cout << "basket " << basket_count << endl;

			for (int i = 0; i<NUM_INTERVENTIONS; i++)
			{
				cout << array_of_baskets[basket_count][i] << " ";
			}
			cout << endl;
		}
	*/
		//Based on input basket_use (range from 1-128) determine which combination will be run. 
		for (int i = 0; i < NUM_INTERVENTIONS; i++)
		{
			basket[i] = array_of_baskets[basket_use][i];
		}
#endif 

#if PROBABILISTIC_RUN
                        
                        prob_pull[84] = interv_effect_size[10][10] ;
                        prob_pull[85] = interv_effect_size[19][0] ;
                        prob_pull[86] = interv_effect_size[19][5] ;
                        prob_pull[87] = interv_effect_size[20][0] ;
                        prob_pull[88] = interv_effect_size[20][5] ;
                        prob_pull[89] = interv_effect_size[21][0] ;
                        prob_pull[90] = interv_effect_size[21][5] ;
                        prob_pull[91] = interv_effect_size[22][0] ;
                        prob_pull[92] = interv_effect_size[22][5] ;
                        prob_pull[93] = interv_effect_size[23][3] ;
                        prob_pull[94] = interv_effect_size[24][3] ;
                        
                        
                        
                        //cost range between 0.5 - 1.5 of the mean
                        prob_pull[95] = interv_costs_per_person[10] ;
                        prob_pull[96] = interv_costs_per_person[19] ;
                        prob_pull[97] = interv_costs_per_person[20] ;
                        prob_pull[98] = interv_costs_per_person[21] ;
                        prob_pull[99] = interv_costs_per_person[22] ;
                        prob_pull[100] = interv_costs_per_person[23] ;
                        prob_pull[101] = interv_costs_per_person[24] ;
                        
		prob_pull[1] = prop_alcohol[0] = uniform_a_b(prob_values[1][0], prob_values[1][1], &seed); //prob_num 1 
		prob_pull[2] = prop_alcohol[1] = uniform_a_b(prob_values[2][0], prob_values[2][1], &seed); //prob_num 2
		prob_pull[3] = uniform_a_b(prob_values[3][0], prob_values[3][1], &seed); 

		 //PROP_GAY_M PROP_BI_M
		prob_pull[4] = propOrientInPop[0][1] = uniform_a_b(prob_values[4][0], prob_values[4][1], &seed); //prob_num 4
		prob_pull[5] = propOrientInPop[0][2] = uniform_a_b(prob_values[5][0], prob_values[5][1], &seed); //prob_num 5
		propOrientInPop[0][0] = 1 - propOrientInPop[0][2] - propOrientInPop[0][1];

		prob_pull[6] = uniform_a_b(prob_values[6][0], prob_values[6][1], &seed); //prob_num 6;

		prob_pull[7] = pmtct_meds = uniform_a_b(prob_values[7][0], prob_values[7][1], &seed); //prob_num 7

		double dur_multiplier = uniform_a_b(0.5, 1.5, &seed);  //MULTIPLIER TO KEEP RANK ORDER between 0.5 and 1.5
		prob_pull[24] = dur[1] *= dur_multiplier; 
		prob_pull[25] = dur[2] *= dur_multiplier;
		prob_pull[26] = dur[3] *= dur_multiplier;

		prob_pull[27] = conc_given[2] = uniform_a_b(prob_values[27][0], prob_values[27][1], &seed); //prob_num 27
		prob_pull[28] = conc_given[3] = uniform_a_b(prob_values[28][0], prob_values[28][1], &seed); //prob_num 28
		
		prob_pull[29] = uniform_a_b(prob_values[29][0], prob_values[29][1], &seed);
		prob_pull[30] = uniform_a_b(prob_values[30][0], prob_values[30][1], &seed);
		prob_pull[31] = uniform_a_b(prob_values[31][0], prob_values[31][1], &seed);
		prob_pull[32] = normal(CD4SD_F, prob_values[32][3], &seed);
		prob_pull[33] = uniform_a_b(prob_values[33][0], prob_values[33][1], &seed);
		
		double utility_mult = uniform_a_b(-0.05, 0.05, &seed); //MULTIPLIER TO KEEP RANK ORDER between -0.05 and 0.05
		prob_pull[37] = UTILITY_CD4_LESS_THAN_50 + utility_mult;
		prob_pull[38] = UTILITY_CD4_50_TO_200 + utility_mult;
		prob_pull[39] = UTILITY_CD4_ABOVE_200 + utility_mult;
		
		prob_pull[40] = uniform_a_b(prob_values[40][0], prob_values[40][1], &seed);
		prob_pull[41] = alpha_raw[1][0] = normal(alpha_raw[1][0], prob_values[41][2], &seed); //prob_num 41
		prob_pull[42] = alpha_raw[0][1] = normal(alpha_raw[0][1], prob_values[42][2], &seed); //prob_num 42
		prob_pull[44] = alpha_raw[0][0] = normal(alpha_raw[0][0], prob_values[44][2], &seed); //prob_num 44
		prob_pull[43] = alpha_IDU_raw = normal(alpha_IDU_raw, prob_values[43][2], &seed); //prob_num 43
		prob_pull[45] = uniform_a_b(prob_values[45][0], prob_values[45][1], &seed); 
		prob_pull[46] = uniform_a_b(prob_values[46][0], prob_values[46][1], &seed);
		prob_pull[47] = uniform_a_b(prob_values[47][0], prob_values[47][1], &seed);
		
		//condom nonuse pathway in baselineRR_raw matrix 
		prob_pull[48] = baselineRR_raw[0][0] = uniform_a_b(prob_values[48][0], prob_values[48][1], &seed); //prob_num 48
		prob_pull[49] = baselineRR_raw[0][1] = lognormal(baselineRR_raw[0][1], prob_values[49][0], prob_values[49][1], seed, 999); //prob_num 49 ***CHECK***
		prob_pull[50] = baselineRR_raw[0][2] = lognormal(baselineRR_raw[0][2], prob_values[50][0], prob_values[50][1], seed, 999); //prob_num 50
		prob_pull[51] = baselineRR_raw[0][5] = lognormal(baselineRR_raw[0][5], prob_values[51][0], prob_values[51][1], seed, 999); //prob_num 51
		prob_pull[52] = baselineRR_raw[0][6] = lognormal(baselineRR_raw[0][6], prob_values[52][0], prob_values[52][1], seed, 999); //prob_num 52 ***CHECK***
		prob_pull[53] = baselineRR_raw[0][7] = uniform_a_b(prob_values[53][0], prob_values[53][1], &seed); //prob_num 53
		prob_pull[55] = baselineRR_raw[0][9] = lognormal(baselineRR_raw[0][9], prob_values[55][0], prob_values[55][1], seed, 999); //prob_num 55
		prob_pull[56] = baselineRR_raw[0][10] = lognormal(baselineRR_raw[0][10], prob_values[56][0], prob_values[56][1], seed, 999); //prob_num 56
		prob_pull[57] = baselineRR_raw[0][12] = lognormal(baselineRR_raw[0][12], prob_values[57][0], prob_values[57][1], seed, 999); //prob_num 57

		//HIV testing pathway in baseline raw matrix 
		prob_pull[58] = baselineRR_raw[1][0] = uniform_a_b(prob_values[58][0], prob_values[58][1], &seed); //prob_num 58
		prob_pull[59] = baselineRR_raw[1][2] = lognormal(baselineRR_raw[1][2], prob_values[59][0], prob_values[59][1], seed, 999); //prob_num 59
		prob_pull[60] = baselineRR_raw[1][3] = lognormal(baselineRR_raw[1][3], prob_values[60][0], prob_values[60][1], seed, 999); //prob_num 60
		prob_pull[61] = baselineRR_raw[1][5] = lognormal(baselineRR_raw[1][5], prob_values[61][0], prob_values[61][1], seed, 999); //prob_num 61
		prob_pull[62] = baselineRR_raw[1][6] = lognormal(baselineRR_raw[1][6], prob_values[62][0], prob_values[62][1], seed, 999); //prob_num 62
		prob_pull[63] = baselineRR_raw[1][8] = lognormal(baselineRR_raw[1][8], prob_values[63][0], prob_values[63][1], seed, 999); //prob_num 63

		//Probability of ART nonadherence in baseline raw matrix 
		prob_pull[64] = baselineRR_raw[3][0] = uniform_a_b(prob_values[64][0], prob_values[64][1], &seed); //prob_num 64
		prob_pull[65] = baselineRR_raw[3][7] = lognormal(baselineRR_raw[3][7], prob_values[65][0], prob_values[65][1], seed, 999); //prob_num 65 ***CHECK***
		prob_pull[66] = baselineRR_raw[3][8] = lognormal(baselineRR_raw[3][8], prob_values[66][0], prob_values[66][1], seed, 999); //prob_num 66 ***CHECK***

		//Untreated STI pathway in baseline raw matrix 
		prob_pull[67] = baselineRR_raw[5][0] = uniform_a_b(prob_values[67][0], prob_values[67][1], &seed); //prob_num 67
		prob_pull[68] = baselineRR_raw[5][1] = uniform_a_b(prob_values[68][0], prob_values[68][1], &seed); //prob_num 68
		prob_pull[69] = baselineRR_raw[5][6] = lognormal(baselineRR_raw[5][6], prob_values[69][0], prob_values[69][1], seed, 999); //prob_num 69 ***CHECK***
		prob_pull[70] = baselineRR_raw[5][7] = uniform_a_b(prob_values[70][0], prob_values[70][1], &seed); //prob_num 70
		prob_pull[71] = baselineRR_raw[5][8] = lognormal(baselineRR_raw[5][8], prob_values[71][0], prob_values[71][1], seed, 999); //prob_num 71 ***CHECK***

		//Uncircumcised pathway in baseline raw matrix 
		prob_pull[72] = baselineRR_raw[7][0] = uniform_a_b(prob_values[72][0], prob_values[72][1], &seed); //prob_num 72
		prob_pull[73] = baselineRR_raw[7][9] = lognormal(baselineRR_raw[7][9], prob_values[73][0], prob_values[73][1], seed, 999); //prob_num 73 ***CHECK***
		prob_pull[74] = baselineRR_raw[7][10] = lognormal(baselineRR_raw[7][10], prob_values[74][0], prob_values[74][1], seed, 999); //prob_num 74
		prob_pull[75] = baselineRR_raw[7][11] = lognormal(baselineRR_raw[7][11], prob_values[75][0], prob_values[75][1], seed, 999); //prob_num 75
		prob_pull[76] = baselineRR_raw[7][12] = lognormal(baselineRR_raw[7][12], prob_values[76][0], prob_values[76][1], seed, 999); //prob_num 76 ***CHECK***

		 //min PROP_GAY_F PROP_BI_F 
		prob_pull[77] = propOrientInPop[1][1] = uniform_a_b(prob_values[77][0], prob_values[77][1], &seed);//prob_num 77
		prob_pull[78] = propOrientInPop[1][2] = uniform_a_b(prob_values[78][0], prob_values[78][1], &seed); //prob_num 78
		propOrientInPop[1][0] = 1 - propOrientInPop[1][2] - propOrientInPop[1][1];


		//probability of mother to child transmission
		double motherToChild_mult = uniform_a_b(0.5, 1.5, &seed); //MULTIPLIER TO KEEP RANK ORDER between 0.5 and 1.5
		prob_pull[79] = motherToChild[0] *= motherToChild_mult; 
		prob_pull[80] = motherToChild[1] *= motherToChild_mult;
		prob_pull[81] = motherToChild[2] *= motherToChild_mult;
		prob_pull[82] = motherToChild[3] *= motherToChild_mult;
		prob_pull[83] = motherToChild[4] *= motherToChild_mult;
                    
                        
#endif 
                        
#if ONE_WAY_SENSITIVITY_ANALYSIS
	if (RUN_NUM==3 || RUN_NUM==4) prop_alcohol[0]=sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 5 || RUN_NUM==6) prop_alcohol[1]=sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 9 || RUN_NUM==10) 
	{ //max PROP_GAY_M 
		propOrientInPop[0][0] = 1 - PROP_BI_M - sensitivity_values[RUN_NUM];
		propOrientInPop[0][1]= sensitivity_values[RUN_NUM]; 
	}
	else if (RUN_NUM == 11 || RUN_NUM==12)
	{ //min PROP_BI_M
		propOrientInPop[0][0] = 1 - sensitivity_values[RUN_NUM] - PROP_GAY_M; 
		propOrientInPop[0][2]=  sensitivity_values[RUN_NUM];
	}
	else if (RUN_NUM==15 || RUN_NUM==16) pmtct_meds = sensitivity_values[RUN_NUM];

	else if (RUN_NUM == 155 || RUN_NUM == 156)
	{ //min PROP_GAY_F
		propOrientInPop[1][0] = 1 - PROP_BI_F - sensitivity_values[RUN_NUM];
		propOrientInPop[1][1] = sensitivity_values[RUN_NUM];
	}
	else if (RUN_NUM == 157 || RUN_NUM == 158)
	{ //min PROP_BI_F
		propOrientInPop[1][0] = 1 - sensitivity_values[RUN_NUM] - PROP_GAY_F;
		propOrientInPop[1][2] = sensitivity_values[RUN_NUM];
	}
	else if (RUN_NUM == 49 || RUN_NUM==50) dur[1] = sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 51 || RUN_NUM==52) dur[2] = sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 53 || RUN_NUM==54) dur[3] = sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 95 || RUN_NUM == 96) conc_given[1]= sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 55 || RUN_NUM == 56) conc_given[2]= sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 57 || RUN_NUM == 58) conc_given[3]= sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 81 || RUN_NUM == 82) alpha_raw[1][0]=sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 83 || RUN_NUM == 84) alpha_raw[0][1]=sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 87 || RUN_NUM == 88) alpha_raw[0][0]=sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 85 || RUN_NUM == 86) alpha_IDU_raw = sensitivity_values[RUN_NUM]; 
	//condom nonuse pathway in baselineRR_raw matrix 
	else if (RUN_NUM == 97 || RUN_NUM == 98) baselineRR_raw[0][0] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 99 || RUN_NUM == 100) baselineRR_raw[0][1] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 101 || RUN_NUM == 102) baselineRR_raw[0][2] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 103 || RUN_NUM == 104) baselineRR_raw[0][5] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 105 || RUN_NUM == 106) baselineRR_raw[0][6] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 107 || RUN_NUM == 108) baselineRR_raw[0][7] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 109 || RUN_NUM == 110) baselineRR_raw[0][8] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 111 || RUN_NUM == 112) baselineRR_raw[0][9] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 113 || RUN_NUM == 114) baselineRR_raw[0][10] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 115 || RUN_NUM == 116) baselineRR_raw[0][12] = sensitivity_values[RUN_NUM];
	//HIV testing pathway in baseline raw matrix 
	else if (RUN_NUM == 117 || RUN_NUM == 118) baselineRR_raw[1][0] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 119 || RUN_NUM == 120) baselineRR_raw[1][2] = sensitivity_values[RUN_NUM]; 
	else if (RUN_NUM == 121 || RUN_NUM == 122) baselineRR_raw[1][3] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 123 || RUN_NUM == 124) baselineRR_raw[1][5] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 125 || RUN_NUM == 126) baselineRR_raw[1][6] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 127 || RUN_NUM == 128) baselineRR_raw[1][8] = sensitivity_values[RUN_NUM];
	//Probability of ART nonadherence in baseline raw matrix 
	else if (RUN_NUM == 129 || RUN_NUM == 130) baselineRR_raw[3][0] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 131 || RUN_NUM == 132) baselineRR_raw[3][7] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 133 || RUN_NUM == 134) baselineRR_raw[3][8] = sensitivity_values[RUN_NUM];
	//Untreated STI pathway in baseline raw matrix 
	else if (RUN_NUM == 135 || RUN_NUM == 136) baselineRR_raw[5][0] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 137 || RUN_NUM == 138) baselineRR_raw[5][1] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 139 || RUN_NUM == 140) baselineRR_raw[5][6] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 141 || RUN_NUM == 142) baselineRR_raw[5][7] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 143 || RUN_NUM == 144) baselineRR_raw[5][8] = sensitivity_values[RUN_NUM];
	//Uncircumcised pathway in baseline raw matrix 
	else if (RUN_NUM == 145 || RUN_NUM == 146) baselineRR_raw[7][0] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 147 || RUN_NUM == 148) baselineRR_raw[7][9] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 149 || RUN_NUM == 150) baselineRR_raw[7][10] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 151 || RUN_NUM == 152) baselineRR_raw[7][11] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 153 || RUN_NUM == 154) baselineRR_raw[7][12] = sensitivity_values[RUN_NUM];
	//probability of mother to child transmission
	else if (RUN_NUM == 159 || RUN_NUM == 160) motherToChild[0] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 161 || RUN_NUM == 162) motherToChild[1] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 163 || RUN_NUM == 164) motherToChild[2] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 165 || RUN_NUM == 166) motherToChild[3] = sensitivity_values[RUN_NUM];
	else if (RUN_NUM == 167 || RUN_NUM == 168) motherToChild[4] = sensitivity_values[RUN_NUM];
#endif 



				outfile<<"Int num: ";
				cout<<"Int num: ";

				if (basket_num==-1) 
				{
					outfile<<"No intervention";
					cout<<"No intervention";
				}

				//print the list of interventions in the current basket
				for (int i=0; i<NUM_INTERVENTIONS; i++) 
				{ 
					if (basket[i]) outfile<<i<<" "; 
					if (basket[i]) cout<<i<<" ";
					if (atLeastOneInterventionTurnedOn == false && basket[i]==true) atLeastOneInterventionTurnedOn = true;
				}	

				outfile<<endl;
				cout<<endl;

				num_infections[basket_num] = 0;

				if (!atLeastOneInterventionTurnedOn)
				{
					outfile<<"Interventions turned off"<<endl;
					cout<<endl<<"Interventions turned off"<<endl<<endl;
				}

				initialize();

				cout<<endl<<"**************AFTER MOVING PEOPLE (IF AT ALL)******************"<<endl<<endl;
				startAndEnd<<endl<<"**************AFTER MOVING PEOPLE (IF AT ALL)******************"<<endl<<endl;
				printNumInEachStatus();

				for (t = 0; t < STOPTIME; t+=DT)
				{						
					cout << "Current Time is: " << t << endl;
	
					if (CALIBRATE && t >= 14.83 && cd4_treat_threshold != CD4_TREAT_THRESHOLD_NOV_2011) // in Nov, 2011 the treatment threshold changes to 350
					{
						cd4_treat_threshold = CD4_TREAT_THRESHOLD_NOV_2011;
						choose_rate_files();
					}

					initializeArrays(Comp);					//used for compensateForDifferentialMortality and calculateBetaandLambda
					balancePartnerships();
					setTotalPeopleinModel(Comp);			//Checks for any problems with the population such as a negative value in a compartment
					LifeYears_and_QALYS();
					alcoholIntervention();					//Apply the SBIRT intervention to any alcs every cycle
					calculateBetaandLambda(Comp);			//TEST is it right to call this with Pop?

					//During the calibration phase, the rate of testing and the rate of being linked to care is changing over time
					//So we are recalculating the baseline probabilities during the calibration phase from 1997 - 2014
					//If we are not running in calibration mode, we should only calculate these values once at the beginning of the model run
					//This happens in the initialize() function.
					if (CALIBRATE && t <= NUM_YEARS_CALIBRATION + 0.0001)
					{
						rampUpValuesDuringCalibration();
						setBaselineProbs();						
						setBaseLineProbsTestedLinkedandAdherence();		//Get the rates from INF to DET and from DET to CARE, also the adherence
						if (isEqual(t, NUM_YEARS_CALIBRATION) && SAVE_POPULATION) save_population();
					}
					calcPropRes(Comp);						//determine proportion in each resistance category (use propRes to store total number in each category for now)

					if (PRINT_CYCLE_BY_CYCLE_OUTPUT) outputData();

					Runge_Kutta();	

                    if (!atLeastOneInterventionTurnedOn) num_infections_no_interv += incidence * DT;	//Keep a running count of infections averted.  'incidence' is the annual incidence calculated in Runge Kutta.  
					else num_infections[basket_num] += incidence * DT;

					num_new_infections_alcohol += alcohol_incidence * DT;
					num_new_infections_IDU += IDU_incidence*DT; 

					//Calculate the cost of the interventions and Care and Treatment
					calculateCost();

					sumOfNumInfectionsPerInfected += incidence * DT / totalinfectedsinmodel;
					sumOfNumDeathsPerInfected += HIVdeathrate * DT / totalinfectedsinmodel;
					sumOfCostOfCareAndTreatPerInfected += costOfCareAndTreatThisCycle / totalinfectedsinmodel; //Note that costOfCareAndTreatThisCycle already has a *DT term in it

#if PRINT_NUM_IN_EACH_STRATA_EVERY_YEAR
					if (isEqual(floor(t+.00000001),t) && t>0 && t < STOPTIME) printNumInEachStatus();
#endif
				}
				
				printNumInEachStatus();

				if (PRINT_CYCLE_BY_CYCLE_OUTPUT) outputData();

				if (basket_num != -1) createBasketOutput();
				else totalDeadAIDS_no_interv = getTotalDeadAIDS();

#ifdef RUN_ALL_NO_SCRIPT
				if (basket_num == -1) interv_num--;
				basket_num++;					//increment the basket count
			}
#endif

#if ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION
			//Before scripted was this:	}
#endif

#if OPTIMIZING_ALCOHOL_INTERVENTION

		///*
		//}}}}
		//*/
	}
#endif

#ifdef TESTING_LINKAGE_ADHERENCE
}}
#endif

#if DEBUG_SENSITIVITY && ONE_WAY_SENSITIVITY_ANALYSIS
cout << "run num: " << RUN_NUM << "\tsensitivity name " << sensitivity_name[RUN_NUM] << "\tsensitivity value " << sensitivity_values[RUN_NUM] << endl;

if (RUN_NUM == 3 || RUN_NUM == 4) cout << "prop_alcohol[0] " << prop_alcohol[0] << endl;
else if (RUN_NUM == 5 || RUN_NUM == 6) cout << "prop_alcohol[1] " << prop_alcohol[1] << endl;
else if (RUN_NUM == 9 || RUN_NUM == 10) cout << "propOrientInPop[0][0] " << propOrientInPop[0][0] << "\tpropOrientInPop[0][1]" << propOrientInPop[0][1] << endl;
else if (RUN_NUM == 11 || RUN_NUM == 12) cout << "propOrientInPop[0][0] " << propOrientInPop[0][0] << "\tpropOrientInPop[0][2]" << propOrientInPop[0][2] << endl;
else if (RUN_NUM == 49 || RUN_NUM == 50) cout << "dur[1]  " << dur[1] << endl;
else if (RUN_NUM == 51 || RUN_NUM == 52) cout << "dur[2]  " << dur[2] << endl;
else if (RUN_NUM == 53 || RUN_NUM == 54) cout << "dur[3]  " << dur[3] << endl;
else if (RUN_NUM == 95 || RUN_NUM == 96) cout << "conc_given[1] " << conc_given[1] << endl;
else if (RUN_NUM == 55 || RUN_NUM == 56) cout << "conc_given[2] " << conc_given[2] << endl;
else if (RUN_NUM == 57 || RUN_NUM == 58) cout << "conc_given[3] " << conc_given[3] << endl;
else if (RUN_NUM == 81 || RUN_NUM == 82) cout << "alpha_raw[1][0]  " << alpha_raw[1][0] << endl;
else if (RUN_NUM == 83 || RUN_NUM == 84) cout << "alpha_raw[0][1]  " << alpha_raw[0][1] << endl;
else if (RUN_NUM == 87 || RUN_NUM == 88) cout << "alpha_raw[1][1]  " << alpha_raw[1][1] << endl;
else if (RUN_NUM == 85 || RUN_NUM == 86) cout << "alpha_IDU_raw  " << alpha_IDU_raw << endl;
//condom nonuse pathway in baselineRR_raw matrix 
else if (RUN_NUM == 97 || RUN_NUM == 98) cout << "baselineRR_raw[0][0]  " << baselineRR_raw[0][0] << endl;
else if (RUN_NUM == 99 || RUN_NUM == 100) cout << "baselineRR_raw[0][1]  " << baselineRR_raw[0][1] << endl;
else if (RUN_NUM == 101 || RUN_NUM == 102) cout << "baselineRR_raw[0][2]  " << baselineRR_raw[0][2] << endl;
else if (RUN_NUM == 103 || RUN_NUM == 104) cout << "baselineRR_raw[0][5]  " << baselineRR_raw[0][5] << endl;
else if (RUN_NUM == 105 || RUN_NUM == 106) cout << "baselineRR_raw[0][6]  " << baselineRR_raw[0][6] << endl;
else if (RUN_NUM == 107 || RUN_NUM == 108) cout << "baselineRR_raw[0][7]  " << baselineRR_raw[0][7] << endl;
else if (RUN_NUM == 109 || RUN_NUM == 110) cout << "baselineRR_raw[0][8]  " << baselineRR_raw[0][8] << endl;
else if (RUN_NUM == 111 || RUN_NUM == 112) cout << "baselineRR_raw[0][9]  " << baselineRR_raw[0][9] << endl;
else if (RUN_NUM == 113 || RUN_NUM == 114) cout << "baselineRR_raw[0][10]  " << baselineRR_raw[0][10] << endl;
else if (RUN_NUM == 115 || RUN_NUM == 116) cout << "baselineRR_raw[0][12]  " << baselineRR_raw[0][12] << endl;
//HIV testing pathway in baseline raw matrix 
else if (RUN_NUM == 117 || RUN_NUM == 118) cout << "baselineRR_raw[1][0]  " << baselineRR_raw[1][0] << endl;
else if (RUN_NUM == 119 || RUN_NUM == 120) cout << "baselineRR_raw[1][2]  " << baselineRR_raw[1][2] << endl;
else if (RUN_NUM == 121 || RUN_NUM == 122) cout << "baselineRR_raw[1][3]  " << baselineRR_raw[1][3] << endl;
else if (RUN_NUM == 123 || RUN_NUM == 124) cout << "baselineRR_raw[1][5]  " << baselineRR_raw[1][5] << endl;
else if (RUN_NUM == 125 || RUN_NUM == 126) cout << "baselineRR_raw[1][6]  " << baselineRR_raw[1][6] << endl;
else if (RUN_NUM == 127 || RUN_NUM == 128) cout << "baselineRR_raw[1][8]  " << baselineRR_raw[1][8] << endl;
//Probability of ART nonadherence in baseline raw matrix 
else if (RUN_NUM == 129 || RUN_NUM == 130) cout << "baselineRR_raw[3][0]  " << baselineRR_raw[3][0] << endl;
else if (RUN_NUM == 131 || RUN_NUM == 132) cout << "baselineRR_raw[3][7]  " << baselineRR_raw[3][7] << endl;
else if (RUN_NUM == 133 || RUN_NUM == 134) cout << "baselineRR_raw[3][8]  " << baselineRR_raw[3][8] << endl;
//Untreated STI pathway in baseline raw matrix 
else if (RUN_NUM == 135 || RUN_NUM == 136) cout << "baselineRR_raw[5][0]  " << baselineRR_raw[5][0] << endl;
else if (RUN_NUM == 137 || RUN_NUM == 138) cout << "baselineRR_raw[5][1]  " << baselineRR_raw[5][1] << endl;
else if (RUN_NUM == 139 || RUN_NUM == 140) cout << "baselineRR_raw[5][6]  " << baselineRR_raw[5][6] << endl;
else if (RUN_NUM == 141 || RUN_NUM == 142) cout << "baselineRR_raw[5][7]  " << baselineRR_raw[5][7] << endl;
else if (RUN_NUM == 143 || RUN_NUM == 144) cout << "baselineRR_raw[5][8]  " << baselineRR_raw[5][8] << endl;
//Uncircumcised pathway in baseline raw matrix 
else if (RUN_NUM == 145 || RUN_NUM == 146) cout << "baselineRR_raw[7][0]  " << baselineRR_raw[7][0] << endl;
else if (RUN_NUM == 147 || RUN_NUM == 148) cout << "baselineRR_raw[7][9]  " << baselineRR_raw[7][9] << endl;
else if (RUN_NUM == 149 || RUN_NUM == 150) cout << "baselineRR_raw[7][10]  " << baselineRR_raw[7][10] << endl;
else if (RUN_NUM == 151 || RUN_NUM == 152) cout << "baselineRR_raw[7][11]  " << baselineRR_raw[7][11] << endl;
else if (RUN_NUM == 153 || RUN_NUM == 154) cout << "baselineRR_raw[7][12]  " << baselineRR_raw[7][12] << endl;
//sex behavior proportions
else if (RUN_NUM == 17 || RUN_NUM == 18) cout << "propSexOrientInAct[M][STRAIGHT][0] " << propSexOrientInAct[M][STRAIGHT][0] << "\tpropSexOrientInAct[M][STRAIGHT][1] " << propSexOrientInAct[M][STRAIGHT][1] << endl; 
else if (RUN_NUM == 19 || RUN_NUM == 20) cout << "propSexOrientInAct[M][STRAIGHT][1] " << propSexOrientInAct[M][STRAIGHT][1] << "\tpropSexOrientInAct[M][STRAIGHT][2] " << propSexOrientInAct[M][STRAIGHT][2] << endl;
else if (RUN_NUM == 21 || RUN_NUM == 22) cout << "propSexOrientInAct[M][STRAIGHT][1] " << propSexOrientInAct[M][STRAIGHT][1] << "\tpropSexOrientInAct[M][STRAIGHT][3] " << propSexOrientInAct[M][STRAIGHT][3] << endl;
else if (RUN_NUM == 23 || RUN_NUM == 24) cout << "propSexOrientInAct[F][STRAIGHT][0] " << propSexOrientInAct[F][STRAIGHT][0] << "\tpropSexOrientInAct[F][STRAIGHT][1] " << propSexOrientInAct[F][STRAIGHT][1] << endl;
else if (RUN_NUM == 25 || RUN_NUM == 26) cout << "propSexOrientInAct[F][STRAIGHT][1] " << propSexOrientInAct[F][STRAIGHT][1] << "\tpropSexOrientInAct[F][STRAIGHT][2] " << propSexOrientInAct[F][STRAIGHT][2] << endl;
else if (RUN_NUM == 27 || RUN_NUM == 28) cout << "propSexOrientInAct[F][STRAIGHT][1] " << propSexOrientInAct[F][STRAIGHT][1] << "\tpropSexOrientInAct[F][STRAIGHT][3] " << propSexOrientInAct[F][STRAIGHT][3] << endl;
else if (RUN_NUM == 29 || RUN_NUM == 30) cout << "propSexOrientInAct[M][GAY][1] " << propSexOrientInAct[M][GAY][1] << "\tpropSexOrientInAct[M][GAY][2] " << propSexOrientInAct[M][GAY][2] << endl;
else if (RUN_NUM == 31 || RUN_NUM == 32) cout << "propSexOrientInAct[M][GAY][2] " << propSexOrientInAct[M][GAY][2] << "\tpropSexOrientInAct[M][GAY][3] " << propSexOrientInAct[M][GAY][3] << endl;
else if (RUN_NUM == 33 || RUN_NUM == 34) cout << "propSexOrientInAct[F][GAY][0] " << propSexOrientInAct[F][GAY][0] << "\tpropSexOrientInAct[F][GAY][1] " << propSexOrientInAct[F][GAY][1] << endl;
else if (RUN_NUM == 35 || RUN_NUM == 36) cout << "propSexOrientInAct[F][GAY][1] " << propSexOrientInAct[F][GAY][1] << "\tpropSexOrientInAct[F][GAY][2] " << propSexOrientInAct[F][GAY][2] << endl;
else if (RUN_NUM == 37 || RUN_NUM == 38) cout << "propSexOrientInAct[F][GAY][1] " << propSexOrientInAct[F][GAY][1] << "\tpropSexOrientInAct[F][GAY][3] " << propSexOrientInAct[F][GAY][3] << endl;
else if (RUN_NUM == 39 || RUN_NUM == 40) cout << "propSexOrientInAct[M][BI][1] " << propSexOrientInAct[M][BI][1] << "\tpropSexOrientInAct[M][BI][2] " << propSexOrientInAct[M][BI][2] << endl;
else if (RUN_NUM == 41 || RUN_NUM == 42) cout << "propSexOrientInAct[M][BI][2] " << propSexOrientInAct[M][BI][2] << "\tpropSexOrientInAct[M][BI][3] " << propSexOrientInAct[M][BI][3] << endl;
else if (RUN_NUM == 43 || RUN_NUM == 44) cout << "propSexOrientInAct[F][BI][0] " << propSexOrientInAct[F][BI][0] << "\tpropSexOrientInAct[F][BI][1] " << propSexOrientInAct[F][BI][1] << endl;
else if (RUN_NUM == 45 || RUN_NUM == 46) cout << "propSexOrientInAct[F][BI][1] " << propSexOrientInAct[F][BI][1] << "\tpropSexOrientInAct[F][BI][2] " << propSexOrientInAct[F][BI][2] << endl;
else if (RUN_NUM == 47 || RUN_NUM == 48) cout << "propSexOrientInAct[F][BI][1] " << propSexOrientInAct[F][BI][1] << "\tpropSexOrientInAct[F][BI][3] " << propSexOrientInAct[F][BI][3] << endl;
else if (RUN_NUM == 7 || RUN_NUM == 8) cout << "propIDUinACT[0][IDU] " << propIDUinAct[0][IDU] << "\tpropIDUinACT[1][IDU] " << propIDUinAct[1][IDU] << "\tpropIDUinACT[2][IDU] " << propIDUinAct[2][IDU] << "\tpropIDUinACT[3][IDU] " << propIDUinAct[3][IDU] << endl;

else if (RUN_NUM == 155 || RUN_NUM == 156) cout << "propOrientInPop[1][0] " << propOrientInPop[1][0] << "\tpropOrientInPop[1][1]" << propOrientInPop[1][1] << endl;
else if (RUN_NUM == 157 || RUN_NUM == 158) cout << "propOrientInPop[1][0] " << propOrientInPop[1][0] << "\tpropOrientInPop[1][2]" << propOrientInPop[1][2] << endl;

else if (RUN_NUM == 159 || RUN_NUM == 160) cout << "motherToChild[0] " << motherToChild[0] << endl; 
else if (RUN_NUM == 161 || RUN_NUM == 162) cout << "motherToChild[1] " << motherToChild[1] << endl;
else if (RUN_NUM == 163 || RUN_NUM == 164) cout << "motherToChild[2] " << motherToChild[2] << endl;
else if (RUN_NUM == 165 || RUN_NUM == 166) cout << "motherToChild[3] " << motherToChild[3] << endl;
else if (RUN_NUM == 167 || RUN_NUM == 168) cout << "motherToChild[4] " << motherToChild[4] << endl;

#endif 

	print_time();		//print the time the run completes

	return (0);
}




void initialize()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res, p_sex;
	incidence=0;
	totalDeadAIDS = 0;		//used in createBasketOutput
	sumOfNumInfectionsPerInfected = 0;			//used to generate number of infections per infected person
	sumOfNumDeathsPerInfected = 0;
	sumOfCostOfCareAndTreatPerInfected = 0;		//used to generate cost of Care and Treatment per infected person
	num_new_infections_alcohol = 0;		//used to track total number of alcohol
	num_new_infections_IDU = 0; 

	costOfCareAndTreat = 0;					//the cost of care and treatment that comes from the progression model lookup tables
	discounted_costOfCareAndTreat = 0;
	
	//initialize survival time variables
	//We will maintain LYs, QALYs, and disc LYs and disc QALYs for 0)  All pats, 1)  All HIV+  2)  All On Care + Treat
	for (int group=0; group<3; group++)
	{
		total_life_years[group] = total_qa_life_years[group] = total_disc_life_years[group] = total_disc_qa_life_years[group] = 0;
	}

	//Used for balancing partnerships.  Only needs to be done once.
	initialize_willPartnerWithMatrix();

	//set up the array for proportion of people with/without alcohol abuse
	//Prop_alc used to be prop_alc[NUM_ALC][NUM_SEX] but changed to be over all the elements of the "who"
	for (status=0; status<NUM_ALIVE_STATE; status++)
	{
		for (sex=0; sex<NUM_SEX; sex++)
		{
			for (act=0; act<NUM_ACT; act++)
			{
				prop_alc[status][sex][act][ALC] = prop_alcohol[sex];
				prop_alc[status][sex][act][NONALC] = 1 - prop_alcohol[sex];
			}
		}
	}

	//initialize the prop_orient_act array
	propSexOrientInAct[M][STRAIGHT][0] = PROP_ACT_0_STRAIGHT_M;
	propSexOrientInAct[M][STRAIGHT][1] = PROP_ACT_1_STRAIGHT_M;
	propSexOrientInAct[M][STRAIGHT][2] = PROP_ACT_2_STRAIGHT_M;
	propSexOrientInAct[M][STRAIGHT][3] = PROP_ACT_3_STRAIGHT_M;

	propSexOrientInAct[F][STRAIGHT][0] = PROP_ACT_0_STRAIGHT_F;
	propSexOrientInAct[F][STRAIGHT][1] = PROP_ACT_1_STRAIGHT_F;
	propSexOrientInAct[F][STRAIGHT][2] = PROP_ACT_2_STRAIGHT_F;
	propSexOrientInAct[F][STRAIGHT][3] = PROP_ACT_3_STRAIGHT_F;

	propSexOrientInAct[M][GAY][0] = PROP_ACT_0_GAY_M;
	propSexOrientInAct[M][GAY][1] = PROP_ACT_1_GAY_M;
	propSexOrientInAct[M][GAY][2] = PROP_ACT_2_GAY_M;
	propSexOrientInAct[M][GAY][3] = PROP_ACT_3_GAY_M;

	propSexOrientInAct[F][GAY][0] = PROP_ACT_0_GAY_F;
	propSexOrientInAct[F][GAY][1] = PROP_ACT_1_GAY_F;
	propSexOrientInAct[F][GAY][2] = PROP_ACT_2_GAY_F;
	propSexOrientInAct[F][GAY][3] = PROP_ACT_3_GAY_F;

	propSexOrientInAct[M][BI][0] = PROP_ACT_0_BI_M;
	propSexOrientInAct[M][BI][1] = PROP_ACT_1_BI_M;
	propSexOrientInAct[M][BI][2] = PROP_ACT_2_BI_M;
	propSexOrientInAct[M][BI][3] = PROP_ACT_3_BI_M;

	propSexOrientInAct[F][BI][0] = PROP_ACT_0_BI_F;
	propSexOrientInAct[F][BI][1] = PROP_ACT_1_BI_F;
	propSexOrientInAct[F][BI][2] = PROP_ACT_2_BI_F;
	propSexOrientInAct[F][BI][3] = PROP_ACT_3_BI_F;

#if ONE_WAY_SENSITIVITY_ANALYSIS
		//Straight male proportions
		if (RUN_NUM == 17 || RUN_NUM==18 )
		{
			propSexOrientInAct[M][STRAIGHT][0] = sensitivity_values[RUN_NUM]; 						 
			propSexOrientInAct[M][STRAIGHT][1] = 1- PROP_ACT_2_STRAIGHT_M - PROP_ACT_3_STRAIGHT_M - sensitivity_values[RUN_NUM]; 
		}
		else if (RUN_NUM == 19 || RUN_NUM==20 )
		{	 
			propSexOrientInAct[M][STRAIGHT][1] = 1- PROP_ACT_0_STRAIGHT_M - PROP_ACT_3_STRAIGHT_M - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[M][STRAIGHT][2] = sensitivity_values[RUN_NUM]; 
		}
		
		else if (RUN_NUM == 21 || RUN_NUM==22 )
		{				 
			propSexOrientInAct[M][STRAIGHT][1] = 1- PROP_ACT_0_STRAIGHT_M - PROP_ACT_2_STRAIGHT_M - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[M][STRAIGHT][3] = sensitivity_values[RUN_NUM]; 
		}

		//Straight female proportions
		else if (RUN_NUM == 23 || RUN_NUM==24 )
		{
			propSexOrientInAct[F][STRAIGHT][0] = sensitivity_values[RUN_NUM]; 						 
			propSexOrientInAct[F][STRAIGHT][1] = 1- PROP_ACT_2_STRAIGHT_F - PROP_ACT_3_STRAIGHT_F - sensitivity_values[RUN_NUM]; 
		}
		else if (RUN_NUM == 25 || RUN_NUM==26 )
		{	 
			propSexOrientInAct[F][STRAIGHT][1] = 1- PROP_ACT_0_STRAIGHT_F - PROP_ACT_3_STRAIGHT_F - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[F][STRAIGHT][2] = sensitivity_values[RUN_NUM]; 
		}
		
		else if (RUN_NUM == 27 || RUN_NUM==28 )
		{				 
			propSexOrientInAct[F][STRAIGHT][1] = 1- PROP_ACT_0_STRAIGHT_F - PROP_ACT_2_STRAIGHT_F - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[F][STRAIGHT][3] = sensitivity_values[RUN_NUM]; 
		}
		 
		//Gay male proportions 
		else if (RUN_NUM == 29 || RUN_NUM==30 )
		{	 
			propSexOrientInAct[M][GAY][1] = sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[M][GAY][2] = 1- PROP_ACT_0_GAY_M - PROP_ACT_3_GAY_M - sensitivity_values[RUN_NUM]; 
		}
		
		else if (RUN_NUM == 31 || RUN_NUM == 32 )
		{				 
			propSexOrientInAct[M][GAY][2] = 1- PROP_ACT_0_GAY_M - PROP_ACT_1_GAY_M - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[M][GAY][3] = sensitivity_values[RUN_NUM]; 
		}

		//Gay female proportions
		else if (RUN_NUM == 33 || RUN_NUM==34 )
		{
			propSexOrientInAct[F][GAY][0] = sensitivity_values[RUN_NUM]; 						 
			propSexOrientInAct[F][GAY][1] = 1- PROP_ACT_2_GAY_F - PROP_ACT_3_GAY_F - sensitivity_values[RUN_NUM]; 
		}
		else if (RUN_NUM == 35 || RUN_NUM==36 )
		{	 
			propSexOrientInAct[F][GAY][1] = 1- PROP_ACT_0_GAY_F - PROP_ACT_3_GAY_F - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[F][GAY][2] = sensitivity_values[RUN_NUM]; 
		}
		
		else if (RUN_NUM == 37 || RUN_NUM==38 )
		{				 
			propSexOrientInAct[F][GAY][1] = 1- PROP_ACT_0_GAY_F - PROP_ACT_2_GAY_F - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[F][GAY][3] = sensitivity_values[RUN_NUM]; 
		}

		//Bi male proportions
		else if (RUN_NUM == 39 || RUN_NUM==40 )
		{	 
			propSexOrientInAct[M][BI][1] = sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[M][BI][2] = 1- PROP_ACT_0_BI_M - PROP_ACT_3_BI_M - sensitivity_values[RUN_NUM]; 
		}
		
		else if (RUN_NUM == 41 || RUN_NUM==42 )
		{				 
			propSexOrientInAct[M][BI][2] = 1- PROP_ACT_0_BI_M - PROP_ACT_1_BI_M - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[M][BI][3] = sensitivity_values[RUN_NUM]; 
		}

		//Bi female proportions
		else if (RUN_NUM == 43 || RUN_NUM==44 )
		{
			propSexOrientInAct[F][BI][0] = sensitivity_values[RUN_NUM]; 						 
			propSexOrientInAct[F][BI][1] = 1- PROP_ACT_2_BI_F - PROP_ACT_3_BI_F - sensitivity_values[RUN_NUM]; 
		}
		else if (RUN_NUM == 45 || RUN_NUM==46 )
		{	 
			propSexOrientInAct[F][BI][1] = 1- PROP_ACT_0_BI_F - PROP_ACT_3_BI_F - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[F][BI][2] = sensitivity_values[RUN_NUM]; 
		}
		
		else if (RUN_NUM == 47 || RUN_NUM==48 )
		{				 
			propSexOrientInAct[F][BI][1] = 1- PROP_ACT_0_BI_F - PROP_ACT_2_BI_F - sensitivity_values[RUN_NUM]; 
			propSexOrientInAct[F][BI][3] = sensitivity_values[RUN_NUM]; 
		}
#endif
		
#if PROBABILISTIC_RUN
		//Straight male proportions
		prob_pull[8] = propSexOrientInAct[M][STRAIGHT][0] = normal(PROP_ACT_0_STRAIGHT_M, prob_values[8][2], &seed); //prob_num 8
		prob_pull[9] = propSexOrientInAct[M][STRAIGHT][2] = uniform_a_b(prob_values[9][0], prob_values[9][1], &seed); //prob_num 9
		prob_pull[10] = propSexOrientInAct[M][STRAIGHT][3] = uniform_a_b(prob_values[10][0], prob_values[10][1], &seed); //prob_num 10
		propSexOrientInAct[M][STRAIGHT][1] = 1 - propSexOrientInAct[M][STRAIGHT][0] - propSexOrientInAct[M][STRAIGHT][2] - propSexOrientInAct[M][STRAIGHT][3];
	
		//Straight female proportions
		prob_pull[11] = propSexOrientInAct[F][STRAIGHT][0] = uniform_a_b(prob_values[11][0], prob_values[11][1], &seed); //prob_num 11
		prob_pull[12] = propSexOrientInAct[F][STRAIGHT][2] = uniform_a_b(prob_values[12][0], prob_values[12][1], &seed); //prob_num 12
		prob_pull[13] = propSexOrientInAct[F][STRAIGHT][3] = uniform_a_b(prob_values[13][0], prob_values[13][1], &seed); //prob_num 13
		propSexOrientInAct[F][STRAIGHT][1] = 1 - propSexOrientInAct[F][STRAIGHT][0] - propSexOrientInAct[F][STRAIGHT][2] - propSexOrientInAct[F][STRAIGHT][3];

		//Gay male proportions 
		prob_pull[14] = propSexOrientInAct[M][GAY][1] = normal(PROP_ACT_1_GAY_M, prob_values[14][2], &seed); //prob_num 14
		prob_pull[15] = propSexOrientInAct[M][GAY][3] = uniform_a_b(prob_values[15][0], prob_values[15][1], &seed); //prob_num 15
		propSexOrientInAct[M][GAY][2] = 1 - propSexOrientInAct[M][GAY][0] - propSexOrientInAct[M][GAY][1] - propSexOrientInAct[M][GAY][3];

		//Gay female proportions
		prob_pull[16] = propSexOrientInAct[F][GAY][0] = uniform_a_b(prob_values[16][0], prob_values[16][1], &seed); //prob_num 16
		prob_pull[17] = propSexOrientInAct[F][GAY][2] = uniform_a_b(prob_values[17][0], prob_values[17][1], &seed); //prob_num 17
		prob_pull[18] = propSexOrientInAct[F][GAY][3] = uniform_a_b(prob_values[18][0], prob_values[18][1], &seed); //prob_num 18
		propSexOrientInAct[F][GAY][1] = 1 - propSexOrientInAct[F][GAY][0] - propSexOrientInAct[F][GAY][2] - propSexOrientInAct[F][GAY][3];

		//Bi male proportions
		prob_pull[19] = propSexOrientInAct[M][BI][1] = normal(PROP_ACT_1_BI_M, prob_values[19][2], &seed); //prob_num 19
		prob_pull[20] = propSexOrientInAct[M][BI][3] = uniform_a_b(prob_values[20][0], prob_values[20][1], &seed); //prob_num 20
		propSexOrientInAct[M][BI][2] = 1 - propSexOrientInAct[M][BI][0] - propSexOrientInAct[M][BI][1] - propSexOrientInAct[M][BI][3];

		//Bi female proportions
		prob_pull[21] = propSexOrientInAct[F][BI][0] = uniform_a_b(prob_values[21][0], prob_values[21][1], &seed); //prob_num 21
		prob_pull[22] = propSexOrientInAct[F][BI][2] = uniform_a_b(prob_values[22][0], prob_values[22][1], &seed); //prob_num 22
		prob_pull[23] = propSexOrientInAct[F][BI][3] = uniform_a_b(prob_values[23][0], prob_values[23][1], &seed); //prob_num 23
		propSexOrientInAct[F][BI][1] = 1 - propSexOrientInAct[F][BI][0] - propSexOrientInAct[F][BI][2] - propSexOrientInAct[F][BI][3];
    
#if DEBUG_SENSITIVITY
    
    for (int i = 0; i < NUM_PROB_TESTS + NUM_PROB_COST_EFFECT; i++)
    {
        cout << probabilistic_name[i] << "\t" << prob_pull[i] << endl;
    }
    
#endif
    
#endif

	calculate_propIDUinAct();

	//p_sex is the sex doing the infecting
	for (sex=0; sex<NUM_SEX; sex++)
	{
		for (p_sex=0; p_sex<NUM_SEX; p_sex++)
		{
			for (vl=0; vl<NUM_VL; vl++)
			{
				alpha[sex][p_sex][vl] = alpha_raw[sex][p_sex] * alpha_vl_mults[vl];
				alpha_IDU[vl] = alpha_IDU_raw * alpha_vl_mults[vl];
			}
		}
	}


	//Initialize all states to zero, some will be overwritten with new values below.
	LOOP_THROUGH_ALL_COMPARTMENTS

		Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res]=0;

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	setPropOrientActAlcIDU();

	//If we aren't running in calibrate mode, the first thing we want to do is load the saved population
	if (!CALIBRATE)  
	{	
		read_population();					
		printNumInEachStatus();
	}
	
	//Set up the initial population for 1997
	if (CALIBRATE)
	{
		LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS

		if (status == SUSC && cd4 == 0 && vl == 0 && res == 0) Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] = INITIAL_POP * propOrientActAlcIDU[status][sex][orient][act][alc][idu] * prop_sex_age[sex][age] * (1 - prop_inf_sex_age[sex][age]);
		else if (status > SUSC && act == 3)
		{
			if (status == INF && res == 0) Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] = INITIAL_POP * propOrientActAlcIDU[status][sex][orient][act][alc][idu] * prop_sex_age[sex][age] * prop_inf_sex_age[sex][age] * (1 - PROP_DETECTED) * init_INF_cd4_vl_distr[vl][cd4] / (propSexOrientInAct[sex][orient][act]); 
			else if (status == HVL && res == 0) Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] = INITIAL_POP * propOrientActAlcIDU[status][sex][orient][act][alc][idu]* prop_sex_age[sex][age] * prop_inf_sex_age[sex][age] * init_HVL_cd4_vl_distr[vl][cd4]/(propSexOrientInAct[sex][orient][act]);
			else if (status == DET && res == 0) Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] = INITIAL_POP * propOrientActAlcIDU[status][sex][orient][act][alc][idu] * prop_sex_age[sex][age] * prop_inf_sex_age[sex][age] * PROP_DETECTED * init_INF_cd4_vl_distr[vl][cd4]/(propSexOrientInAct[sex][orient][act]);
		}

			//TEST (this was not part of the original population setup, but the code crashes with no people in the detected state
			//So, FIX THIS.  WE need to determine the proportion of people in detected vs. INF, vs. care

		END_LOOP_THROUGH_ALL_COMPARTMENTS

			//Now that we've allocated, add up the total people in all compartments
			setTotalPeopleinModel(Comp);

		//Check that the number of people in all the compartments equals the initial population
		//Did we do the allocation correctly?
		if (!isEqual(totalPeopleInModel, INITIAL_POP))
		{
			cout<<"The number of people allocated over all compartments does not equal INITIAL_POP.  Total people is "<<totalPeopleInModel<<" INITIAL_POP is "<<INITIAL_POP<<" Exiting."<<endl;
			exit(0);
		}
	}

	//TEST - do we still need this here???
	setTotalPeopleinModel(Comp);
	
#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM==1 || RUN_NUM==2)) 
		set_age_change_rates(int (sensitivity_values[RUN_NUM])); 
#else 
	set_age_change_rates(int (AGE_SEXUAL_DEBUT));
#endif 


	//Now initialize functions, etc. that relate to running interventions and the interventions that move people between compartments (alcohol and fewer partners)
	initialize_interventions();					

	//Note, we need the baseline probabilities for the rate from inf to care, so we need to call this function even if we aren't running interventions.
	if (CALIBRATE) rampUpValuesDuringCalibration();

	setBaselineProbs();	

	//For circumcision intervention (interv #4).  If interv is on, we immediately do a virtual mass-circumcision of adult men according to the interv effect size.
	if (basket[4] == true) findNumMenCircumsizedByIntervAtStartOfRun();

#if ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION
	//Note: This must come after setBaselineProbs() and before setBaseLineProbsTestedLinkedandAdherence() 
	//and before processInterventionsThatMovePeople().
	postProcessBaselineProbsForAlcAnalysis();
#endif

#if OPTIMIZING_PREVENTION_PORTFOLIO && (SET_RR_ALC > 0) // SET_RR_ALC will be 1 for high and 2 for low
	const double newbaselineRRsettings[2][3] = { // High and low settings for Condoms, Adherence, STI
		{2, 3.5, 2.5},
		{1.14, 1.67, 1.36}
	};
	condomRR = newbaselineRRsettings[SET_RR_ALC-1][0]; // Condoms
	adhereRR = newbaselineRRsettings[SET_RR_ALC-1][1]; // Adherence
	stiRR = newbaselineRRsettings[SET_RR_ALC-1][2]; // STI
	//Note: This must come after setBaselineProbs() and before setBaseLineProbsTestedLinkedandAdherence()
	//and before processInterventionsThatMovePeople().
	postProcessBaselineProbsForOptimizingAlcoholIntervention();
#endif

#if OPTIMIZING_ALCOHOL_INTERVENTION

	setAlcoholicProportionsAtStartOfRun();

	cout<<endl<<"**************AFTER reallocating alcoholics ******************"<<endl<<endl;
	startAndEnd<<endl<<"**************AFTER reallocating alcoholics ******************"<<endl<<endl;
	printNumInEachStatus();

	//Note: This must come after setBaselineProbs() and before setBaseLineProbsTestedLinkedandAdherence() 
	//and before processInterventionsThatMovePeople().
	postProcessBaselineProbsForOptimizingAlcoholIntervention();
#endif

#if ALCOHOL_SURFACE
	//Note: This must come after setBaselineProbs() and before setBaseLineProbsTestedLinkedandAdherence() 
	//and before processInterventionsThatMovePeople().
	postProcessBaselineProbsForOptimizingAlcoholIntervention();
#endif

	setPropWomenOnPMTCT();
	choose_rate_files();

	//For setBaseLineProbsTestedLinkedandAdherence(), make sure that initialize_interventions is called before this function (sets the effect sizes)
	//also make sure that setBaselineProbs is called before this function
	setBaseLineProbsTestedLinkedandAdherence();				

#if OPTIMIZING_PREVENTION_PORTFOLIO && SET_ALC_PROPORTIONS
	setAlcoholicProportionsAtStartOfRun();
#endif

	processInterventionsThatMovePeople();					//make sure that initialize_interventions is called before this function (sets the effect sizes)

	//EFFICIENCY!  Calling the function to pre calculate values of HIVPrograte so we don't repeat the same interpolation (for the right level of adherence) again and again.
	if (EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS) preCalculateAdhStratHIVProgRateAndCostArrays();

	//EFFICIENCY!  Calling the function to pre calculate alphaMult values so that we don't calculate the same values repeatedly.
	if (EFFICIENCY_PRECALCULATE_ALPHA_MULT_VALUES) precalc_getAlphaMult();
}

/*

Within the Diff function, here are the possible transfers people can make from each compartment:

if (status==SUSC) 
{
//1)  People will transfer out with new infections
}
if (status==HVL)
{
//1) People will transfer in with new infections
//2) People will transfer out at a steady rate
}

if (status==INF)
{
//1) Transfer in from HVL
//2)  Transfer out by going on care
//3)  In by progress (inner loop)
//4)  Out by progress (inner loop)
//5)  Out by dying of HIV
}

if (status==CARE) 
{
//1) Transfer in from INF
//2)  Transfer out by going on treatment (inner loop)
//3)  In by progress (inner loop)
//4)  Out by progress (inner loop)
//5)  Out by dying of HIV
}

if (status==TREAT) 
{
//1)  Transfer in from CARE 
//2)  Transfer out by loss to follow up
//3)  In by progress (inner loop)
//4)  Out by progress (inner loop)
//5)  Out by dying of HIV
}
*/

void Diff(double Pop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES], double &newInfections, double &HIV_DeathRate, double &newInfectionsAlc, double &numCircumsizedThisCycle, double &newInfectionsIDU)
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	int vl_d, cd4_d, treat_d, res_d, dead_d;
	double trans_amount;
	double insum=0, outsum=0;
	double numInfectedBirths = 0;		//Really, the point in time infected births per year
	double newInfectedAdults = 0;		//Really, the point in time number of new infected adults per year

	double trans_in[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];
	double trans_out[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];


	//initialize arrays to zero
	memcpy(trans_in, zerod_array, sizeof(trans_in));
	memcpy(trans_out, zerod_array, sizeof(trans_out));

	calculateAgingBirthsDeaths(Pop, numInfectedBirths, numCircumsizedThisCycle);		//This uses the proportion resistant, so make sure to call after calcPropRes();

	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS  

		//TEST (Keep this in the loop, dummy!)
		//if (age>=3) lambda[sex][orient][act][age]=0.1;

		//People will transfer from SUSC to HVL with new infections (only adults are sexually active) 
		if (status==SUSC && vl==0 && cd4==0 && res==0 && age >= ADULTHOOD) 
		{
			trans_amount = lambda[sex][orient][act][age][alc][idu] * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;
			if (DEBUG_TRANSITIONS &&  trans_amount>0) printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_out", trans_amount);

			//Add new infections to calculate incidence (note that this will only calculate ADULT new infections)
			newInfectedAdults += trans_amount;	

			//Add new infections to calculate alcohol related incidence (note that this will only calculate ADULT new infections - ok, since kids don't drink!)
			if (alc == ALC || alc == ALC_TREATED_NOT_CURED) newInfectionsAlc += trans_amount;
			if (idu == IDU) newInfectionsIDU += trans_amount; 

			//Apportion the newly infecteds into VL/CD4/Res groups
			for (cd4_d=0; cd4_d<NUM_CD4; cd4_d++)
			{
				for (vl_d=0; vl_d<NUM_VL; vl_d++)
				{	
					for (res_d=0; res_d<NUM_RES; res_d++)
					{
						if (new_Infections_Go_Immediately_On_Treat)
							trans_in[TREAT][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d] += trans_amount * probcd4vl[sex][cd4_d][vl_d] * propRes[res_d];
						else //normal behavior
							trans_in[HVL][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d] += trans_amount * probcd4vl[sex][cd4_d][vl_d] * propRes[res_d];

						if (DEBUG_TRANSITIONS && trans_in[HVL][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d]>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", HVL, sex, orient, act, age, alc, idu, vl_d, cd4_d, res_d, "trans_in", trans_in[HVL][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d]);
					}
				}
			}
		}

		//All susceptibles should have a vl, cd4 and res of 0.  Print error if otherwise.
		else if (status==SUSC && (vl!=0 || cd4!=0 || res!=0) && Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res]>0)
		{
			cout<<"There were susceptibles who did not have VL, CD4 and Res values of 0.  Exiting.\n";
			exit(0);
		}

		//Transfer to INF from HVL
		if (status==HVL)
		{
			trans_amount = HVL_TO_INF * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;            
			trans_in[INF][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;

			if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_out", trans_amount);
			if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", INF, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_in", trans_amount);
		}

		//Transfer to Detected from INF
		if (status==INF)
		{	
			trans_amount = inf_to_det_rate[sex][orient][act][alc][idu]  * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;
			trans_in[DET][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;

			if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_out", trans_amount);
			if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", DET, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_in", trans_amount);
		}

		//Transfer to Care from Detected
		if (status==DET)
		{
			trans_amount = det_to_care_rate[sex][orient][act][alc][idu] * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;            
			trans_in[CARE][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;

			if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_out", trans_amount);
			if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", CARE, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_in", trans_amount);
		}

		if (!turn_Off_HIV_Mortality)	//Sometimes we want to turn off HIV mortality while testing the model
		{
			//All infected states (not including initial infection) can transfer out by dying of HIV
			if (status>=INF)		//(You have to be infected to die of HIV)
			{
				//If not calibrating the precalculated values to avoid calling getHIVProgRate() repeatedly for the same settings
				//the indices are ... preCalcAdhStratifiedHIVProgRate[status][cd4][vl][res][cd4_d][vl_d][treat_d][res_d][dead_d][sex][act][alc] = getAdhStratifiedHIVProgRate(status, cd4, vl, res, cd4_d, vl_d, treat_d, res_d, dead_d, sex, act, alc)
	
				//TEST
				//if (status>=CARE)
				//{
				//cout<<preCalcAdhStratifiedHIVProgRate[status][cd4][vl][res][0][0][0][0][1][sex][act][alc]<<endl;
				//cout<<getAdhStratifiedHIVProgRate(status, cd4, vl, res, 0, 0, 0, 0, 1, sex, act, alc)<<endl;
				//}

				if (EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS)	
					trans_amount = preCalcAdhStratifiedHIVProgRate[status][cd4][vl][res][0][0][0][0][1][sex][orient][act][alc][idu] * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
				else 
					trans_amount = getAdhStratifiedHIVProgRate(status, cd4, vl, res, 0, 0, 0, 0, 1, sex, orient, act, alc, idu) * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			
				//Mortality adjustment for <5 yr olds who are HIV infected and not in treatment
				if (age==0 && status<=CARE)
				{
					#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM==13 || RUN_NUM==14 )) || PROBABILISTIC_RUN
						#if PROBABILISTIC_RUN
							trans_amount *= prob_pull[6]; //prob_num 6;
						#else
							trans_amount *=sensitivity_values[RUN_NUM];
						#endif 
					#else 
						trans_amount *= MORT_MULT_CHILDREN_NOT_ON_TREAT;
					#endif
				}

				trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;            
				trans_in[DEAD_HIV][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;

				HIV_DeathRate += trans_amount;

				if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_out", trans_amount);
				if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", DEAD_HIV, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_in", trans_amount);
			}
		}

		//Now process people with HIV moving to other alive states
		if (status >= INF)
		{
			for (vl_d=0; vl_d<NUM_VL; vl_d++)
			{
				for (cd4_d=0; cd4_d<NUM_CD4; cd4_d++)
				{
					for (res_d=0; res_d<NUM_RES; res_d++)
					{
						for (treat_d=0; treat_d<2; treat_d++)
						{
							dead_d=0;		//for all of these destinations, people don't arrive dead (we did that above)


							//If not calibrating the precalculated values to avoid calling getHIVProgRate() repeatedly for the same settings
							if (EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS)
								trans_amount = preCalcAdhStratifiedHIVProgRate[status][cd4][vl][res][cd4_d][vl_d][treat_d][res_d][dead_d][sex][orient][act][alc][idu] * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
							else 
								trans_amount = getAdhStratifiedHIVProgRate(status, cd4, vl, res, cd4_d, vl_d, treat_d, res_d, dead_d, sex, orient, act, alc, idu) * Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];

							//Here are people transferring out
							trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] += trans_amount;
							if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, "trans_out", trans_amount);

							//These people stay within INF or DET, but can change their vl and cd4
							if (status==INF || status==DET)
							{
								if (treat_d==0)
								{
									trans_in[status][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d] += trans_amount;
									if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl_d, cd4_d, res_d, "trans_in", trans_amount);
								}
								else if (trans_amount > 0)			
								{
									cout << "People shouldn't be able to transfer from INF or DET to TREAT.  Exiting." <<endl;
									exit (0);
								}
							}

							//People in Care can either go on Treatment, or stay in Care
							if (status==CARE) 
							{
								//Transfer out by going on treatment
								if (treat_d==1)
								{
									trans_in[TREAT][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d] += trans_amount;
									if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", TREAT, sex, orient, act, age, alc, idu, vl_d, cd4_d, res_d, "trans_in", trans_amount);
								}

								//Stay in Care, but change CD4/VL/res
								else  //this means treat_d == 0
								{
									trans_in[status][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d] += trans_amount;
									if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl_d, cd4_d, res_d, "trans_in", trans_amount);
								}
							}

							if (status==TREAT) 
							{
								//Stay on TREAT, but change CD4/VL/res
								if (treat_d==1)
								{
									trans_in[status][sex][orient][act][age][alc][idu][vl_d][cd4_d][res_d] += trans_amount;
									if (DEBUG_TRANSITIONS && trans_amount>0) printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\n", status, sex, orient, act, age, alc, idu, vl_d, cd4_d, res_d, "trans_in", trans_amount);
								}
								else if (trans_amount > 0)
								{
									cout << "People shouldn't be able to transfer anywhere else from here (since we've already processed people moving to DEAD)." <<endl;
									exit (0);
								}
							}

						}}}}}

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	if (DEBUG_TRANSITIONS) cout<<"Transition Totals: "<<endl;

	LOOP_THROUGH_ALL_COMPARTMENTS

		dPop[status][sex][orient][act][age][alc][idu][vl][cd4][res] = trans_in[status][sex][orient][act][age][alc][idu][vl][cd4][res] - trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] + age_births_deaths[status][sex][orient][act][age][alc][idu][vl][cd4][res];

		if (DEBUG_TRANSITIONS)
		{
			if (trans_in[status][sex][orient][act][age][alc][idu][vl][cd4][res] !=0) printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\ttrans_in\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, trans_in[status][sex][orient][act][age][alc][idu][vl][cd4][res]);
			if (trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res] !=0) printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\ttrans_out\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res]);
			//if (dPop[status][sex][orient][act][age][alc][idu][vl][cd4][res]!=0) printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\tdiff\t%f\n", status, sex, orient, act, age, alc, idu, vl, cd4, res, dPop[status][sex][orient][act][age][alc][idu][vl][cd4][res]);

			insum += trans_in[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			outsum += trans_out[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		}

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	if (DEBUG_TRANSITIONS) cout<<"IN "<<insum<<"  OUT "<<outsum<<endl<<endl;

	//Really, the point in time number of new infections per year
	newInfections = numInfectedBirths + newInfectedAdults;		
}

void balancePartnerships()
{
	int sex, orient, act, age;
	int p_sex, p_orient, p_age, p_act;
	double numPartnerships1stPartner, numPartnerships2ndPartner; 
	bool mutualPartnershipsOffered;

	determineProbabilityMatrix();

	for (sex=0; sex<NUM_SEX; sex++)
	{
		for (orient=0; orient<NUM_ORIENT; orient++)
		{
			for (act=FIRST_ACTIVE_CLASS; act<NUM_ACT; act++)
			{
				for (age=ADULTHOOD; age<NUM_AGE; age++)
				{
					for (p_sex=0; p_sex<NUM_SEX; p_sex++)
					{
						for (p_orient=0; p_orient<NUM_ORIENT; p_orient++)
						{
							for (p_act=FIRST_ACTIVE_CLASS; p_act<NUM_ACT; p_act++)
							{
								for (p_age=ADULTHOOD; p_age<NUM_AGE; p_age++)
								{	
									//compromised duration
									dur_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = 
										dur_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = 
										compromise(dur[act], dur[p_act]);

									//num partnerships offered to each group of opposite sex = N * conc * p
									numPartnerships1stPartner = nSexOrientActivityAge[sex][orient][act][age] * conc_given[act] * prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
									numPartnerships2ndPartner = nSexOrientActivityAge[p_sex][p_orient][p_act][p_age] * conc_given[p_act] * prob_matrix[p_sex][p_orient][p_act][p_age][sex][orient][act][age];

									//If either side is offering zero partnerships (maybe because the groups can't partner with each other),
									//then don't try to balance.  It will result in ugly divide-by-zero errors.
									if (numPartnerships1stPartner > 0 && numPartnerships2ndPartner > 0) mutualPartnershipsOffered=true;
									else mutualPartnershipsOffered = false;

									if (mutualPartnershipsOffered)
									{
										//update the number of partnerships based on the compromise
										num_partnerships_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = 
											num_partnerships_matrix[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = 
											compromise(numPartnerships1stPartner, numPartnerships2ndPartner);

										//update the concurancy matrix based on the compromised number of parternships 
										conc_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = num_partnerships_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] / (nSexOrientActivityAge[sex][orient][act][age] * prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age]);
										conc_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = num_partnerships_matrix[p_sex][p_orient][p_act][p_age][sex][orient][act][age] / (nSexOrientActivityAge[p_sex][p_orient][p_act][p_age] * prob_matrix[p_sex][p_orient][p_act][p_age][sex][orient][act][age]);

										if (FREQ_VARIES_WITH_CONC && act!=3)//if (FREQ_VARIES_WITH_CONC)
										{
											//As the number of partners increases, the frequency of sexual acts will also increase but not multiplicatively. 
											//For example, if someone has 1 partner with a frequency of 3, a person with 10 partners will not likely 
											//have a frequency of 30 (10x3).  Therefore, a linear equation was used to determine how frequency varies with 
											//partner number, where partner number is along the x axis and frequency along the y. 
											desired_overall_freq[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = SLOPE * (conc_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age] <= 4 ? conc_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age] : 4) + yintercept;
											desired_overall_freq[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = SLOPE * (conc_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age] <= 4 ? conc_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age] : 4) + yintercept;
										}
										else
										{
											desired_overall_freq[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = FREQ;
											desired_overall_freq[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = FREQ;
										}

										//The frequency of sexual acts with each partner (An individual with r concurrent partners will have a frequency of Freq/r for each partner)
										freq_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = desired_overall_freq[sex][orient][act][age][p_sex][p_orient][p_act][p_age] / conc_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
										freq_per_partnership[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = desired_overall_freq[p_sex][p_orient][p_act][p_age][sex][orient][act][age] / conc_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age];

										//two individuals in a partnership must experience the same number of acts within the partnership (freq of acts) X (duration of relationship)
										acts_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = freq_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age] * dur_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
										acts_per_partnership[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = freq_per_partnership[p_sex][p_orient][p_act][p_age][sex][orient][act][age] * dur_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age];

										//now compromise on the total number of acts per partnership
										acts_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age] =
											acts_per_partnership[p_sex][p_orient][p_act][p_age][sex][orient][act][age] =
											compromise(acts_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age], acts_per_partnership[p_sex][p_orient][p_act][p_age][sex][orient][act][age]);

										//update the freq_per_partnership matrix based on the compromised number of acts with each partner
										freq_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = acts_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age] / dur_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
										freq_per_partnership[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = acts_per_partnership[p_sex][p_orient][p_act][p_age][sex][orient][act][age] / dur_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age];

										//calculate the partner change rate (number of concurrent partners per year) / (duration of each partnership)
										change_rate_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = conc_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age] / dur_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
										change_rate_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age] = conc_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age] / dur_comp[p_sex][p_orient][p_act][p_age][sex][orient][act][age];
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if (debug)
	{
		cout<<endl<<"Duration of partnerships after compromise (first men, then women)"<<endl;
		printMatrix(dur_comp);

		cout<<endl<<"Num partnerships offered after compromise (first men, then women)"<<endl;
		printMatrix(num_partnerships_matrix);

		cout<<endl<<"Concurrency after compromise (first men, then women)"<<endl;
		printMatrix(conc_comp);

		cout<<endl<<"Desired overall frequency (acts per year) - cumulative over all partners (freq_per_partnership varies with number of concurrent partners) (first men, then women)"<<endl;
		printMatrix(desired_overall_freq);

		cout<<endl<<"Frequency (acts per partnership per year) after compromise (first men, then women)"<<endl;
		printMatrix(freq_per_partnership);

		cout<<endl<<"Acts per partnership after compromise (first men, then women)"<<endl;
		printMatrix(acts_per_partnership);

		cout<<endl<<"Change rates after compromise (first men, then women)"<<endl;
		printMatrix(change_rate_comp);
	}
}


void outputData()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	double plwha = 0, deadhiv_this_cycle = 0, adult_prevalence;
	static bool firstPass=1;
	static double dead_HIV_cumulative = 0;

	double numSuscWomen, numInfWomen, numDiedAIDSWomen, numDiedAgeWomen;
	double numSuscMen, numInfMen, numDiedAIDSMen, numDiedAgeMen;
	double numSuscChild, numInfChild, numDiedAIDSChild, numDiedAgeChild, numOnTreat;
	double numSexbyAct[NUM_SEX][NUM_ACT]={0};
	double totalIDU, numInfIDU; 

	numSuscWomen=0; numInfWomen=0; numDiedAIDSWomen=0; numDiedAgeWomen=0;
	numSuscMen=0; numInfMen=0; numDiedAIDSMen=0; numDiedAgeMen=0;
	numSuscChild=0; numInfChild=0; numDiedAIDSChild=0; numDiedAgeChild=0;
	numOnTreat = 0; totalIDU = 0, numInfIDU=0;


	//print column headers
	if (firstPass) 
	{
		outfile<<"Years\t";
		outfile<<"Total Susceptible\tIncidence\tAdult Prevalence\tPLWHA\tTotal Dead HIV\tSusc Men\tSusc Women\tSusc Children\tInf Men\tInf Women\tInf Children\tDied AIDS Men\tDied AIDS Women\tDied AIDS Children\tDied Age Men\tDied Age Women\tDied Age Children\tIDU_incidence\tNum New IDU Infections\tOn Treat\tTotal IDU\tTotal Inf IDU\t";

		if (OUTPUT_BY_ACT_GROUP)
		{
			for (sex=0; sex<NUM_SEX; sex++)
			{
				for (act=0; act<NUM_ACT; act++)
				{
					outfile<<"\t"<<(sex==0?"Men":"Women")<< " Act "<<act;
				}
			}
		}

		outfile<<endl;
		firstPass=0;
	}


	LOOP_THROUGH_ALL_COMPARTMENTS	

		//Find the total number of sexually active women

		if (status==SUSC)
		{
			if (age>=ADULTHOOD) 
			{
				if (sex==0) numSuscMen+= Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
				else numSuscWomen += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];

				if (OUTPUT_BY_ACT_GROUP) numSexbyAct[sex][act] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}
			else numSuscChild += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		}

		else if (status>=HVL && status<=TREAT)
		{
			if (age>=ADULTHOOD) 
			{
				if (sex==0) numInfMen+= Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
				else numInfWomen += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}
			else numInfChild += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			if (idu == 1)
			{
				numInfIDU += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}
			if (status == TREAT)
			{
				numOnTreat += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}
		}

		else if (status==DEAD_HIV)
		{
			if (age>=ADULTHOOD) 
			{
				if (sex == 0) numDiedAIDSMen += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
				else numDiedAIDSWomen += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}
			else numDiedAIDSChild += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		}

		else if (status==DEAD_AGE)
		{
			if (age>=ADULTHOOD) 
			{
				if (sex == 0) numDiedAgeMen += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
				else numDiedAgeWomen += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}
			else numDiedAgeChild += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		}

		if (idu==1) totalIDU+=Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	//Create data for calibration
	adult_prevalence = (numInfMen + numInfWomen)/(numInfMen + numInfWomen + numSuscMen + numSuscWomen);
	plwha	= numInfMen + numInfWomen + numInfChild;
	deadhiv_this_cycle = (numDiedAIDSMen + numDiedAIDSWomen + numDiedAIDSChild - dead_HIV_cumulative);
	dead_HIV_cumulative = numDiedAIDSMen + numDiedAIDSWomen + numDiedAIDSChild;

	/*********************************START OUTPUT HERE*********************************************/

	outfile<<t<<"\t"<<(numSuscMen+numSuscWomen+numSuscChild)<<"\t"<<incidence/(numSuscMen+numSuscWomen+numSuscChild)<<"\t"<< adult_prevalence << "\t" << plwha << "\t" << deadhiv_this_cycle << "\t";
	outfile<<numSuscMen<<"\t"<<numSuscWomen<<"\t"<<numSuscChild<<"\t";
	outfile<<numInfMen<<"\t"<<numInfWomen<<"\t"<<numInfChild<<"\t";
	outfile<<numDiedAIDSMen<<"\t"<<numDiedAIDSWomen<<"\t"<<numDiedAIDSChild<<"\t";
	outfile << numDiedAgeMen << "\t" << numDiedAgeWomen << "\t" << numDiedAgeChild << "\t";
	outfile << IDU_incidence << "\t" << num_new_infections_IDU << "\t" << numOnTreat << "\t" << totalIDU << "\t" << numInfIDU << "\t";

	if (OUTPUT_BY_ACT_GROUP)
	{
		for (sex=0; sex<NUM_SEX; sex++)
		{
			for (act=0; act<NUM_ACT; act++)
			{
				outfile<<"\t"<<numSexbyAct[sex][act];
			}
		}
	}

	outfile<<endl;
}

void setup()
{
	int sex, orient, act, alc, idu, vl, status;

	//initialize zerod_array which is a very large array of all 0's used to initialize other arrays in the model
	//Note, this needs to come first because it is used in function calls below.  
	//Note also that this array has the same indices as preCalculatedAlphMultArray which is the biggest array in the model

	for (sex=0; sex<NUM_SEX; sex++){
		for (orient=0; orient<NUM_ORIENT; orient++){
			for (act=0; act<NUM_ACT; act++){
				for (alc=0; alc<NUM_ALC; alc++){
					for (idu=0; idu<NUM_IDU; idu++){
						for (vl=0; vl<NUM_VL; vl++){
							for (status = 0; status < CARESTATUS; status++){
								for (sex=0; sex<NUM_SEX; sex++){
									for (orient=0; orient<NUM_ORIENT; orient++){
										for (act=0; act<NUM_ACT; act++){
											for (alc=0; alc<NUM_ALC; alc++){
												for (idu=0; idu<NUM_IDU; idu++){
													for (vl=0; vl<NUM_VL; vl++){
														zerod_array[sex][orient][act][alc][idu][vl][status][sex][orient][act][alc][idu][vl] = 0;
													}}}}}}}}}}}}};

	//sets the outfile stream to output 4 digits after the decimal point
	outfile.precision(7);
	outfile << fixed;

	cout.precision(7);
	cout << fixed;

	basketfile.precision(7);
	basketfile << fixed;

	startAndEnd.precision(7);
	startAndEnd << fixed;

	//open files for writing to
	outfile.open ("output.txt");
	if (outfile.fail())
	{
		cout<< "Could not open output.txt" << endl;
		exit(0);
	}

	//open files for writing to
	basketfile.open ("basket.txt");
	if (basketfile.fail())
	{
		cout<< "Could not open basket.txt" << endl;
		exit(0);
	}

	//open files for writing to
	startAndEnd.open ("startAndEndStats.txt");
	if (startAndEnd.fail())
	{
		cout<< "Could not open startAndEndStats.txt" << endl;
		exit(0);
	}

	my_pop_file.open("Population.txt");
	if (my_pop_file.fail())
	{
		cout << "Could not open Population.txt" << endl;
		exit(0);
	}
	my_pop_file.precision(25); //Set high precision because this data will be used to populate the array for future runs.


	//Calculate the rates at which new infecteds will go into the different VL/CD4 buckets based
	//on a CD4 and VL distribution
	calculateCD4andVLratesNewInfecteds();

	//The following are just to make test scenarios easier
	turn_Off_HIV_Mortality = TURN_OFF_HIV_MORTALITY;
	turn_Off_Transmission_On_Treatment = TURN_OFF_TRANSMISSION_ON_TREATMENT;
	new_Infections_Go_Immediately_On_Treat = NEW_INFECTIONS_GO_IMMEDIATELY_ON_TREAT;
	cd4_treat_threshold = CD4_TREAT_THRESHOLD;
}

void Runge_Kutta()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	//Runge Kutta requires 4 passes
	double dPop1[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];
	double dPop2[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];
	double dPop3[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];
	double dPop4[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];
	double tmpPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];

	double numAnnualNewInfections[4] = {0, 0, 0, 0}; 
	double numAnnualNewHIVDeaths[4] = {0, 0, 0 , 0};
	double numAnnualNewInfections_alc[4] = {0, 0, 0, 0}; // Temporary storage for the numbers of alcohol and NONALC related infections.
	double numCircumsizedAnnually[4] = {0, 0, 0, 0};
	double numAnnualNewInfections_IDU[4] = { 0, 0, 0, 0 }; 
	int round;	//which call to diff (out of 4)


	if (debug) cout<<"Beginning Diff 1 at time "<<t<<endl<<endl;

	round=0;
	Diff(Comp, numAnnualNewInfections[round], numAnnualNewHIVDeaths[round], numAnnualNewInfections_alc[round], numCircumsizedAnnually[round], numAnnualNewInfections_IDU[round]);
	LOOP_THROUGH_ALL_COMPARTMENTS
		dPop1[status][sex][orient][act][age][alc][idu][vl][cd4][res] = dPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		tmpPop[status][sex][orient][act][age][alc][idu][vl][cd4][res] = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] + DT * dPop1[status][sex][orient][act][age][alc][idu][vl][cd4][res] / 2;
	END_LOOP_THROUGH_ALL_COMPARTMENTS

		if (debug) cout<<"Beginning Diff 2 at time "<<t<<endl<<endl;

	round=1;

	Diff(tmpPop, numAnnualNewInfections[round], numAnnualNewHIVDeaths[round], numAnnualNewInfections_alc[round], numCircumsizedAnnually[round], numAnnualNewInfections_IDU[round]);
	LOOP_THROUGH_ALL_COMPARTMENTS
		dPop2[status][sex][orient][act][age][alc][idu][vl][cd4][res] = dPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		tmpPop[status][sex][orient][act][age][alc][idu][vl][cd4][res] = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] + DT * dPop2[status][sex][orient][act][age][alc][idu][vl][cd4][res] / 2;
	END_LOOP_THROUGH_ALL_COMPARTMENTS

		if (debug) cout<<"Beginning Diff 3 at time "<<t<<endl<<endl;

	round=2;

	Diff(tmpPop, numAnnualNewInfections[round], numAnnualNewHIVDeaths[round], numAnnualNewInfections_alc[round], numCircumsizedAnnually[round], numAnnualNewInfections_IDU[round]);
	LOOP_THROUGH_ALL_COMPARTMENTS
		dPop3[status][sex][orient][act][age][alc][idu][vl][cd4][res] = dPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		tmpPop[status][sex][orient][act][age][alc][idu][vl][cd4][res] = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] + DT * dPop3[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	END_LOOP_THROUGH_ALL_COMPARTMENTS

		if (debug) cout<<"Beginning Diff 4 at time "<<t<<endl<<endl;

	round=3;

	Diff(tmpPop, numAnnualNewInfections[round], numAnnualNewHIVDeaths[round], numAnnualNewInfections_alc[round], numCircumsizedAnnually[round], numAnnualNewInfections_IDU[round]);
	LOOP_THROUGH_ALL_COMPARTMENTS
		dPop4[status][sex][orient][act][age][alc][idu][vl][cd4][res] = dPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] + (dPop1[status][sex][orient][act][age][alc][idu][vl][cd4][res] / 6 + dPop2[status][sex][orient][act][age][alc][idu][vl][cd4][res] / 3 + dPop3[status][sex][orient][act][age][alc][idu][vl][cd4][res] / 3 + dPop4[status][sex][orient][act][age][alc][idu][vl][cd4][res] / 6) * DT;
	END_LOOP_THROUGH_ALL_COMPARTMENTS

	incidence = numAnnualNewInfections[0]/6 + numAnnualNewInfections[1]/3 + numAnnualNewInfections[2]/3 + numAnnualNewInfections[3]/6;
	alcohol_incidence = (numAnnualNewInfections_alc[0]/6 + numAnnualNewInfections_alc[1]/3 + numAnnualNewInfections_alc[2]/3 + numAnnualNewInfections_alc[3]/6);
	IDU_incidence = (numAnnualNewInfections_IDU[0] / 6 + numAnnualNewInfections_IDU[1] / 3 + numAnnualNewInfections_IDU[2] / 3 + numAnnualNewInfections_IDU[3] / 6); 
	HIVdeathrate = numAnnualNewHIVDeaths[0]/6 + numAnnualNewHIVDeaths[1]/3 + numAnnualNewHIVDeaths[2]/3 + numAnnualNewHIVDeaths[3]/6;
	numCircumsizedByIntervThisCycle = (numCircumsizedAnnually[0]/6 + numCircumsizedAnnually[1]/3 + numCircumsizedAnnually[2]/3 + numCircumsizedAnnually[3]/6) * DT;
		return;
}


void initializeArrays(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES])
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	//First initialize all to zero
	num_risky_idu_users = 0;
	double total_adults = 0; 
	memcpy(nSexOrientActivityAge, zerod_array, sizeof(nSexOrientActivityAge));
	memcpy(num_partnerships_offered_by_group, zerod_array, sizeof(num_partnerships_offered_by_group));

	memcpy(numSexuallyActiveWomen,zerod_array, sizeof(numSexuallyActiveWomen));
	memcpy(numAbstinentWomen, zerod_array, sizeof(numAbstinentWomen));


	for (status=0; status<=INDEX_MAX_ACTIVE_STATE; status++)
	{
		for (sex=0; sex<NUM_SEX; sex++)
		{
			for (orient=0; orient<NUM_ORIENT; orient++)
			{
				for (act=0; act<NUM_ACT; act++)
				{
					for (age=ADULTHOOD; age<NUM_AGE; age++)
					{
						for (alc=0; alc<NUM_ALC; alc++)
						{
							for (idu=0; idu<NUM_IDU; idu++)
							{
								for (vl=0; vl<NUM_VL; vl++)
								{
									for (cd4=0; cd4<NUM_CD4; cd4++)
									{
										for (res=0; res<NUM_RES; res++)
										{
											if (act >= FIRST_ACTIVE_CLASS)
											{
												//Finds N*c offered for each row in the mixing matrix, 
												//Note that this does not take any assortativeness into account
												num_partnerships_offered_by_group[sex][orient][act][age] += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res] * conc_given[act];

												if (sex==F) numSexuallyActiveWomen[age] += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];	

												nSexOrientActivityAge[sex][orient][act][age] += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
											}
											else	//Abstinent
											{
												if (sex==F) numAbstinentWomen[age] += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
											}

											if (idu == 1) num_risky_idu_users += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
											total_adults += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
#if KELLYTEST
	double prop_idu_check = num_risky_idu_users / total_adults;
	cout << "num_risky_idu_users: " << num_risky_idu_users << endl;
	cout << "total adults: " << total_adults << endl;
	cout << "proportion idu: " << prop_idu_check << endl;
#endif 
}



void determineProbabilityMatrix()
{
	int sex, orient, act, age;
	int p_sex, p_orient, p_age, p_act;
	double proportional_matrix[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
	double assortative_matrix[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
	double row_sum;
	double total_partnerships_avail_to_group_fully_proportional[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
	double total_partnerships_avail_to_group_fully_assort[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];

	preventBalanceError();

	for (sex=0; sex<NUM_SEX; sex++)
	{
		for (orient=0; orient<NUM_ORIENT; orient++)
		{
			for (act=FIRST_ACTIVE_CLASS; act<NUM_ACT; act++)
			{
				for (age=ADULTHOOD; age<NUM_AGE; age++)
				{
					//initialize
					total_partnerships_avail_to_group_fully_proportional[sex][orient][act][age] = 0;
					total_partnerships_avail_to_group_fully_assort[sex][orient][act][age] = 0;	

					for (p_sex=0; p_sex<NUM_SEX; p_sex++)
					{
						for (p_orient=0; p_orient<NUM_ORIENT; p_orient++)
						{
							for (p_act=FIRST_ACTIVE_CLASS; p_act<NUM_ACT; p_act++)
							{
								for (p_age=ADULTHOOD; p_age<NUM_AGE; p_age++)
								{
									//Add up the N*c for the proportional case and the assortative case
									total_partnerships_avail_to_group_fully_proportional[sex][orient][act][age] += num_partnerships_offered_by_group[p_sex][p_orient][p_act][p_age] * willPartnerWith[sex][orient][p_sex][p_orient];
									total_partnerships_avail_to_group_fully_assort[sex][orient][act][age] += num_partnerships_offered_by_group[p_sex][p_orient][p_act][p_age] * willPartnerWith[sex][orient][p_sex][p_orient] * (age==p_age) * (act==p_act);
								}	
							}
						}
					}

					//Initialize, for making sure rows of probability matrix sum to 1
					row_sum = 0;

					for (p_sex=0; p_sex<NUM_SEX; p_sex++)
					{
						for (p_act=FIRST_ACTIVE_CLASS; p_act<NUM_ACT; p_act++)
						{
							for (p_orient=0; p_orient<NUM_ORIENT; p_orient++)
							{
								for (p_age=ADULTHOOD; p_age<NUM_AGE; p_age++)
								{
									proportional_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] =
										num_partnerships_offered_by_group[p_sex][p_orient][p_act][p_age]
										* willPartnerWith[sex][orient][p_sex][p_orient] 
										/ total_partnerships_avail_to_group_fully_proportional[sex][orient][act][age];

									assortative_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = 
										num_partnerships_offered_by_group[p_sex][p_orient][p_act][p_age] 
									* (age==p_age) * (act==p_act)
									* willPartnerWith[sex][orient][p_sex][p_orient] 
									/ total_partnerships_avail_to_group_fully_assort[sex][orient][act][age];

									#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM == 61 || RUN_NUM == 62)) || PROBABILISTIC_RUN
										#if PROBABILISTIC_RUN //prob_run 30
											prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = 
												prob_pull[30] * proportional_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] +
												(1 - prob_pull[30]) * assortative_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
										#else
											prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = 
												sensitivity_values[RUN_NUM] * proportional_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] +
												(1-sensitivity_values[RUN_NUM]) * assortative_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
										#endif
									#else 
										prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] = 
											EPS * proportional_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] +
											(1-EPS) * assortative_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
									#endif

									//the sum of probabilities in every row should equal 1
									row_sum += prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age];
								}	
							}
						}
					}

					if (!isEqual(row_sum, 1.0))
					{
						cout<<endl<<endl<<"***Error in probability matrix.  Row not summing to one!***  Exiting."<<endl<<endl;
						exit(0);
					}
				}
			}
		}
	}

	//cout<<"Probability Matrix:"<<endl;
	//printMatrix(prob_matrix);
}





//Stole this snippet from http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
inline bool isEqual(double x, double y)
{
	const double maxRelativeError = .000001;
	float relativeError;

	if (x == y) return true;


	if (fabs(y) > fabs(x)) relativeError = (float)fabs((x - y) / y);

	else relativeError = (float)fabs((x - y) / x);

	if (relativeError <= maxRelativeError) return true;

	return false;
} 

void calculateBetaandLambda(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES])
{
	int sex, orient, act, age, alc, vl, idu;
	int p_status, p_sex, p_orient, p_act, p_age, p_vl, p_cd4, p_res, p_alc, p_idu;
	double thisBeta, thisAlpha = 0;
	vl=0; //because all susceptibles have a vl=0

	//initialize lambda array to all 0's
	memcpy(lambda, zerod_array, sizeof(lambda));

	if (debug) cout<<"Lambda: First men, then women"<<endl;

	//Initialize all lambda to 0. This loop was moved to initializeArrays() for efficiency.

	//Calculate lambda for every combination of sex/activity/age/alcohol
	for (sex=0; sex<NUM_SEX; sex++)
	{
		for (orient=0; orient<NUM_ORIENT; orient++)
		{
			if (debug) cout<<endl<<"orient: "<<orient<<endl;

			for (act=0; act<NUM_ACT; act++)
			{
				if (debug) cout<<"Act: "<<act<<"\t";

				for (age=ADULTHOOD; age<NUM_AGE; age++)
				{
					for (alc=0; alc<NUM_ALC; alc++)
					{
						for (idu = 0; idu<NUM_IDU; idu++)
						{
							for (p_sex=0; p_sex<NUM_SEX; p_sex++)
							{
								for (p_orient=0; p_orient<NUM_ORIENT; p_orient++)
								{
									for (p_act = 0; p_act < NUM_ACT; p_act++)
									{
										for (p_age = ADULTHOOD; p_age < NUM_AGE; p_age++)
										{
											//Prevent divide-by-zero error in following calculation
											if (!isEqual(nSexOrientActivityAge[p_sex][p_orient][p_act][p_age], 0))
											{
												//loop through the sexually active infective states (State 0 is susceptible, so start with 1)
												for (p_status = HVL; p_status <= INDEX_MAX_ACTIVE_STATE; p_status++)
												{
													for (p_alc = 0; p_alc < NUM_ALC; p_alc++)
													{
														for (p_idu = 0; p_idu < NUM_IDU; p_idu++)
														{
															for (p_vl = 0; p_vl < NUM_VL; p_vl++)
															{
																for (p_cd4 = 0; p_cd4 < NUM_CD4; p_cd4++)
																{
																	for (p_res = 0; p_res < NUM_RES; p_res++)
																	{
																		/*****************Sexual Transmission*******************/

																		//KIMTEST
																		//There is no reason to calculate sexual transmission if the probability of these two 
																		//groups partnering is zero.  Also, we can assume the transmissbility of women to women is zero.  
																		if (!(prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] == 0) &&
																			!(sex==F && p_sex==F))
																		{
																			//For sexual transmission, we only want to consider the sexually active classes
																			if (act >= FIRST_ACTIVE_CLASS && p_act >= FIRST_ACTIVE_CLASS)
																			{
																				if (p_status == HVL)
																				{
																					switch (p_vl)
																					{
																						// No need to do the math each time!  case (0) alpha = (1-1.5/2.5) * alpha[0] + (1.5/2.5*2.0/3.0) * alpha[1] + (1.5/2.5*1.0/3.0)* alpha[2]; break;
																					case (0) : thisAlpha = .4 * alpha[sex][p_sex][0] + .4 * alpha[sex][p_sex][1] + .2 * alpha[sex][p_sex][2]; break;
																					case (1) : thisAlpha = .5 * alpha[sex][p_sex][2] + .5 * alpha[sex][p_sex][3]; break;
																					case (2) : thisAlpha = .5 * alpha[sex][p_sex][3] + .5 * alpha[sex][p_sex][4]; break;
																					case (3) : thisAlpha = alpha[sex][p_sex][4]; break;
																					case (4) : thisAlpha = alpha[sex][p_sex][4]; break;
																					}
																				}
																				else
																				{
																					thisAlpha = alpha[sex][p_sex][p_vl];
																				}

																				//For testing purposes, we sometimes want to eliminate all chances of transmission on treatment
																				if (turn_Off_Transmission_On_Treatment && p_status == TREAT) thisAlpha = 0;
																				if (TURN_OFF_SEXUAL_TRANSMISSION == 1 || TURN_OFF_ALL_TRANSMISSION) thisAlpha = 0;

																				thisAlpha *= SEX_ALPHA_SCALE;

																				//getAlphaMult gets the multiplier that is a results of any/all internventions that are turned on
																				//EFFICIENCY! 
																				//Note that we are pre-calculated the alpha multipliers because there are a fixed number of combinations that exist and thus a fixed number of return values.

																				//VL is now in here for targeting the alcohol intervention by VL
																				if (EFFICIENCY_PRECALCULATE_ALPHA_MULT_VALUES)
																					thisAlpha *= preCalculatedAlphMultArray[sex][orient][act][alc][idu][vl][p_status][p_sex][p_orient][p_act][p_alc][p_idu][p_vl];

																				//VL is now in here for targeting the alcohol intervention by VL
																				else thisAlpha *= getAlphaMult(sex, orient, act, alc, idu, vl, p_status, p_sex, p_orient, p_act, p_alc, p_idu, p_vl);
																				thisBeta = 1 - (pow((1 - thisAlpha), acts_per_partnership[sex][orient][act][age][p_sex][p_orient][p_act][p_age]));

																				lambda[sex][orient][act][age][alc][idu] += change_rate_comp[sex][orient][act][age][p_sex][p_orient][p_act][p_age] *
																					prob_matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age] * thisBeta * myPop[p_status][p_sex][p_orient][p_act][p_age][p_alc][p_idu][p_vl][p_cd4][p_res]
																				/ nSexOrientActivityAge[p_sex][p_orient][p_act][p_age];
																			}
																		}

																		/*****************IDU Transmission*******************/
																		//lambda due to IDU
																		//check that there are some risky idu users to prevent divide by zero error
																		if (idu == 1 && p_idu == 1 && num_risky_idu_users>0)
																		{
																			if (p_status == HVL)
																			{
																				switch (p_vl)
																				{
																				case (0) : thisAlpha = .4 * alpha_IDU[0] + .4 * alpha_IDU[1] + .2 * alpha_IDU[2]; break;
																				case (1) : thisAlpha = .5 * alpha_IDU[2] + .5 * alpha_IDU[3]; break;
																				case (2) : thisAlpha = .5 * alpha_IDU[3] + .5 * alpha_IDU[4]; break;
																				case (3) : thisAlpha = alpha_IDU[4]; break;
																				case (4) : thisAlpha = alpha_IDU[4]; break;
																				}
																			}
																			else
																			{
																				thisAlpha = alpha_IDU[p_vl];
																			}
																			
																			#if (ONE_WAY_SENSITIVITY_ANALYSIS  && (RUN_NUM==93 || RUN_NUM==94)) || PROBABILISTIC_RUN
																				#if PROBABILISTIC_RUN
																					thisAlpha *= prob_pull[47];
																				#else
																					thisAlpha *= sensitivity_values[RUN_NUM];
																				#endif
																			#else 
																				thisAlpha *= DRUG_ALPHA_SCALE;
																			#endif 

																			#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM >=89 && RUN_NUM <= 92)) || PROBABILISTIC_RUN
																				#if PROBABILISTIC_RUN //prob_run 45 and 46
																					
																				thisBeta = 1 - (pow((1 - thisAlpha), prob_pull[45] / prob_pull[46]));
																					lambda[sex][orient][act][age][alc][idu] += thisBeta *
																						myPop[p_status][p_sex][p_orient][p_act][p_age][p_alc][p_idu][p_vl][p_cd4][p_res]
																						* prob_pull[46] / num_risky_idu_users;
																				#else
																					if  (RUN_NUM == 91 || RUN_NUM == 92) 
																					{
																								thisBeta = 1 - (pow((1 - thisAlpha), (double)SHARED_INJECTIONS_PER_YEAR / sensitivity_values[RUN_NUM]));
																								lambda[sex][orient][act][age][alc][idu] += thisBeta *
																									myPop[p_status][p_sex][p_orient][p_act][p_age][p_alc][p_idu][p_vl][p_cd4][p_res]
																									* sensitivity_values[RUN_NUM] / num_risky_idu_users;
																					}
																					else if (RUN_NUM == 89 || RUN_NUM == 90)
																					{
																							thisBeta = 1 - (pow((1 - thisAlpha), (double)sensitivity_values[RUN_NUM]/IDU_PARTNER_CHANGE_RATE));
																							lambda[sex][orient][act][age][alc][idu] += thisBeta *
																								myPop[p_status][p_sex][p_orient][p_act][p_age][p_alc][p_idu][p_vl][p_cd4][p_res]
																								* IDU_PARTNER_CHANGE_RATE/ num_risky_idu_users;
																					}
																				#endif 
																			#else 
																				thisBeta = 1 - (pow((1 - thisAlpha), (double)SHARED_INJECTIONS_PER_YEAR / IDU_PARTNER_CHANGE_RATE));
																				lambda[sex][orient][act][age][alc][idu] += thisBeta *
																					myPop[p_status][p_sex][p_orient][p_act][p_age][p_alc][p_idu][p_vl][p_cd4][p_res]
																				* IDU_PARTNER_CHANGE_RATE / num_risky_idu_users;
																			#endif 
																		}
																		
																		/*********************end IDU Transmission***********************/
																	}
																}
															}
														}
													}
												}
											}
											
										}
									}
								}
							}
						}
					}
					if (debug) cout<<lambda[sex][orient][act][age][alc][idu]<<" ";
				}
			}
			if (debug) cout<<endl;
		}
		if (debug) cout<<endl;
	}
	if (debug) cout<<endl<<endl;
}



void calculateAgingBirthsDeaths(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES], double &numInfectedBirths, double &numCircumsizedThisCycle)
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	double healthyBirths = 0;
	double infectedBirths[NUM_RES];		//the number of infected births stratified by resistance
	double birthAllocation[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];

	double births;
	double age_in;
	double age_out;
	double age_mort;

	double allocation;
	double total_infected_births = 0;
	double total_uninfected_allocation = 0;
	double total_infected_allocation_INF = 0;
	double total_infected_allocation_TREAT = 0;

	numCircumsizedThisCycle = 0;

	double moms_adj_HVL[NUM_VL];
	double moms_in_vl[NUM_VL];  //mother's VL after taking into account any mtct prevention meds

	//initialize the number of deaths to zero
	memcpy(age_deaths, zerod_array, sizeof(age_deaths));

	//initialize birth allocation to 0 (because may not set all values of array if FIRST_ACTIVE_CLASS is not zero
	memcpy(birthAllocation, zerod_array, sizeof(birthAllocation));

	//initialize the number of infected births in each resistance category to 0
	memcpy(infectedBirths, zerod_array, sizeof(infectedBirths));

	//initialize  to 0
	memcpy(moms_adj_HVL, zerod_array, sizeof(moms_adj_HVL));


	scaleUpFertilityRates();

	if (!TURN_OFF_BIRTHS)	//Sometimes for testing, we want no births
	{
		//Calculate the total number of births across all women in all age and activity classes

		for (status = 0; status < NUM_ALIVE_STATE; status++)
		{
			for (orient = 0; orient < NUM_ORIENT; orient++)
			{
				for (act=FIRST_ACTIVE_CLASS; act<NUM_ACT; act++)
				{
					for (age=ADULTHOOD; age<NUM_AGE; age++)
					{
						for (alc=0; alc<NUM_ALC; alc++)
						{
							for (idu = 0; idu < NUM_IDU; idu++)
							{
								for (cd4 = 0; cd4 < NUM_CD4; cd4++)
								{
									for (res = 0; res < NUM_RES; res++)
									{
										if (status >= HVL && status <= CARE)
										{
											for (int viral_load = 0; viral_load < NUM_VL; viral_load++)
											{
												if (status == HVL)
												{
													if (viral_load == 0) moms_adj_HVL[viral_load] = 0;
													if (viral_load == 1) moms_adj_HVL[viral_load] = 2.0 / 3.0 * myPop[status][F][orient][act][age][alc][idu][viral_load - 1][cd4][res];
													if (viral_load == 2 || viral_load == 3) moms_adj_HVL[viral_load] = 1.0 / 3.0 * myPop[status][F][orient][act][age][alc][idu][viral_load - 2][cd4][res] + 2.0 / 3.0 * myPop[status][F][orient][act][age][alc][idu][viral_load - 1][cd4][res];
													if (viral_load == 4) moms_adj_HVL[viral_load] = myPop[status][F][orient][act][age][alc][idu][viral_load][cd4][res] + myPop[status][F][orient][act][age][alc][idu][viral_load - 1][cd4][res] + 1.0 / 3.0 * myPop[status][F][orient][act][age][alc][idu][viral_load - 2][cd4][res];
												}
												else moms_adj_HVL[viral_load] = myPop[status][F][orient][act][age][alc][idu][viral_load][cd4][res];
											}
										}

										for (vl = 0; vl < NUM_VL; vl++)
										{
											//the "1" is for women
											//if the mother is healthy (not infected)
											if (status == SUSC || TURN_OFF_ALL_TRANSMISSION) healthyBirths += myPop[status][F][orient][act][age][alc][idu][vl][cd4][res] * fertRateScaledUp[age];

											else if (status == TREAT)
											{
												healthyBirths += (1 - motherToChild[vl]) * myPop[status][F][orient][act][age][alc][idu][vl][cd4][res] * fertRateScaledUp[age];
												infectedBirths[res] += motherToChild[vl] * myPop[status][F][orient][act][age][alc][idu][vl][cd4][res] * fertRateScaledUp[age];
											}

											else
											{
												/*From the progression model
												pat->VLdec = VL_DEC_NEVIRAPINE * pat->HIVbaseline / 4.5
												where VL_DEC_NEVIRAPINE = 2.22
												the average VL setpoint for women = 4.27(from (Saathoff, 2010))
												Thus, the VL decrement for an average women on pregnancy meds = 2.22*4.27/4.5 = 2.1065
												This is about equivalent to moving two VL groups down
												*/

												if (vl == 0) moms_in_vl[vl] = moms_adj_HVL[vl] + pmtct_meds * moms_adj_HVL[vl + 1] + pmtct_meds * moms_adj_HVL[vl + 2];
												if (vl == 1 || vl == 2) moms_in_vl[vl] = (1 - pmtct_meds) * moms_adj_HVL[vl] + pmtct_meds * moms_adj_HVL[vl + 2];
												if (vl == 3 || vl == 4) moms_in_vl[vl] = (1 - pmtct_meds) * moms_adj_HVL[vl];

												healthyBirths += (1 - motherToChild[vl]) * moms_in_vl[vl] * fertRateScaledUp[age];
												infectedBirths[res] += motherToChild[vl] * moms_in_vl[vl] * fertRateScaledUp[age];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}



		//Now allocate those births among the sex, age, and activity groups
		//A proportion of babies, defined by "motherToChild" are both infected and enter the first infected state
		//We divide the total healthy and infected births by two because half will be boys, half girls

		for (sex=0; sex<NUM_SEX; sex++)
		{
			for (orient=0; orient<NUM_ORIENT; orient++)
			{
				for (act=0; act<NUM_ACT; act++)
				{
					for (alc=0; alc<NUM_ALC; alc++)
					{
						for (idu = 0; idu<NUM_IDU; idu++)
						{
							//Healthy births
							//propInAct //proportion of men/women of each status/sex/alc combination in each activity group (may have been modified by an intervention)
							allocation = healthyBirths / 2.0 * propOrientActAlcIDU[SUSC][sex][orient][act][alc][idu];
							birthAllocation[SUSC][sex][orient][act][alc][idu][0][0][0] = allocation;

							//Add up the total SUSC allocation (for QA checking)
							total_uninfected_allocation += allocation;

							for (vl = 0; vl < NUM_VL; vl++)
							{
								for (cd4 = 0; cd4 < NUM_CD4; cd4++)
								{
									for (res = 0; res < NUM_RES; res++)
									{
										//Infected births go into "Infected" by default, but some proportion go immediately on Treat

										//Infected births NOT immediately on treatment
										allocation = (1 - PROP_INFECTED_INFANTS_START_IN_TREAT) * infectedBirths[res] / 2.0 * propOrientActAlcIDU[INF][sex][orient][act][alc][idu] * probcd4vl[sex][cd4][vl];
										birthAllocation[INF][sex][orient][act][alc][idu][vl][cd4][res] = allocation;

										//Add up the total INF allocation (for QA checking)
										total_infected_allocation_INF += allocation;

										//Infected births IMMEDIATELY on treatment
										allocation = PROP_INFECTED_INFANTS_START_IN_TREAT * infectedBirths[res] / 2.0 * propOrientActAlcIDU[TREAT][sex][orient][act][alc][idu]* probcd4vl[sex][cd4][vl];
										birthAllocation[TREAT][sex][orient][act][alc][idu][vl][cd4][res] = allocation;

										//Add up the total TREAT allocation (for QA checking)
										total_infected_allocation_TREAT += birthAllocation[TREAT][sex][orient][act][alc][idu][vl][cd4][res];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	

	/**************** TEST THAT BIRTH ALLOCATIONS ADD UP!!! *******************/

	//Add up the total number of infected births
	for (res=0; res<NUM_RES; res++) total_infected_births += infectedBirths[res];
	
	//Check that the number of infected births allocated equals the number of infected babies born
	if (!isEqual(total_infected_births, (total_infected_allocation_INF + total_infected_allocation_TREAT)))
	{
		cout<<"The number of infected births allocated does not equal the number of infected babies born.  Exiting."<<endl;
		exit(0);
	}

	//Check that the number of healthy births allocated equals the number of healthy babies born
	if (!isEqual(healthyBirths, total_uninfected_allocation))
	{
		cout<<"The number of healthy births allocated does not equal the number of healthy babies born.  Exiting."<<endl;
		exit(0);
	}
	/****************************************************************************/

	//Now determine the change in population in each compartment including
	//births, deaths, aging in, and aging out
	LOOP_THROUGH_ALL_COMPARTMENTS

		if (status < NUM_ALIVE_STATE)
		{
			if (age == 0)
			{	
				births = birthAllocation[status][sex][orient][act][alc][idu][vl][cd4][res];
				age_in = 0;
				age_out = age_out_rate[age] * myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}

			else
			{	
				births = 0;
				age_in = age_out_rate[age - 1] * myPop[status][sex][orient][act][age - 1][alc][idu][vl][cd4][res];
				age_out = age_out_rate[age] * myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];

				//The circumcision intervention is intervention #4.  
				//The circumcision path is pathway #7.
				//The number of men circumsized by the intervention is the number suitable for circumcision
				//multiplied by one minus the intervention effect size.
				
				//ANIK
				if (basket[4] && age==ADULTHOOD) numCircumsizedThisCycle += (age_in * interv[4].who[status][sex][orient][act][alc][idu][vl])
					* baselineProb[7][status][sex][orient][act][alc][idu]
					* (1 - interv_effect_size[4][7]);
			}

			age_mort = ageMortRate[sex][age] * myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];

			//TEST  No births, no aging out of population
			//if (age==NUM_AGE-1) age_out=0;
			//births = 0; 

			age_births_deaths[status][sex][orient][act][age][alc][idu][vl][cd4][res] = births + age_in - age_out - age_mort;

			//For the age death terminal state
			age_deaths[sex][orient][act][age][alc][idu][vl][cd4][res] += age_mort;

			if (TURN_ON_DEMOGRAPHICS==0) 
			{
				age_births_deaths[status][sex][orient][act][age][alc][idu][vl][cd4][res] = 0;
				age_deaths[sex][orient][act][age][alc][idu][vl][cd4][res] = 0;
			}
		}

		if (status == DEAD_AGE) age_births_deaths[status][sex][orient][act][age][alc][idu][vl][cd4][res] = age_deaths[sex][orient][act][age][alc][idu][vl][cd4][res];

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	//Add up the total number of infected births (used to calculate overall incidence)
	for (res=0; res<NUM_RES; res++) numInfectedBirths += infectedBirths[res];
}

void printMatrix(double matrix[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE])
{
	int sex, orient, act, age, p_sex, p_orient, p_act, p_age;

	cout<<endl;

	for (sex=0; sex<NUM_SEX; sex++)
	{
		cout<<"Sex: "<<sex<<endl;
		for (orient=0; orient<NUM_ORIENT; orient++)
		{
			for (act=FIRST_ACTIVE_CLASS; act<NUM_ACT; act++)
			{
				for (age=ADULTHOOD; age<NUM_AGE; age++)
				{
					for (p_sex=0; p_sex<NUM_SEX; p_sex++)
					{
						for (p_orient=0; p_orient<NUM_ORIENT; p_orient++)
						{
							for (p_act=FIRST_ACTIVE_CLASS; p_act<NUM_ACT; p_act++)
							{
								for (p_age=ADULTHOOD; p_age<NUM_AGE; p_age++)
								{
									cout<<matrix[sex][orient][act][age][p_sex][p_orient][p_act][p_age]<<" ";
								}
							}
						}
					}
					cout<<endl;
				}
			}
		}
	}
}


//In the prograte tables, the order of columns is:
//curr_cd4  curr_vl  curr_onTreat  curr_res dest_cd4  dest_vl  dest_onTreat  dest_res dead  count  rate
//Also, reminder:  double HIVprograte[NUM_ADHERENCE][CARESTATUS][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_ALIVEORDEAD];

void readHIVprogrates(string &off_care_path, string &on_care_root)
{
	int curr_cd4, curr_vl, curr_onTreat, curr_resistance;
	int dest_cd4, dest_vl, dest_onTreat, dest_resistance, dead, thecount;
	float cost;
	string filepath, nameend=".txt";
	double therate;
	char header_line[1000];


	memcpy(HIVprograte, zerod_array, sizeof(HIVprograte));
	memcpy(costCareTreat, zerod_array, sizeof(costCareTreat));

	//Open the off care rate file
	if((stateprob = fopen( off_care_path.c_str(), "r" )) == NULL )
	{
		printf("ERROR opening %s\n", off_care_path.c_str() );
		exit(0);
	}
	cout << "Reading off care rate file: "<< off_care_path << endl;

	// We need to ignore the first line of the file. That line is a human readable header.
	fgets(header_line, 1000, stateprob); 

	//Read in the offCare model rates
	do 
	{
		fscanf(stateprob,"%d%d%d%d%d%d%d%d%d%d%lf%f", &curr_cd4, &curr_vl, &curr_onTreat, &curr_resistance, &dest_cd4, &dest_vl, &dest_onTreat, &dest_resistance, &dead, &thecount, &therate, &cost);

		//read in both the rates and costs 
		HIVprograte[OFFCARE_INDEX][0/* off care */][curr_cd4][curr_vl][curr_onTreat][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][dead] = (double)therate;
		costCareTreat[OFFCARE_INDEX][0/* off care */][curr_cd4][curr_vl][curr_onTreat][curr_resistance] = (double)cost;

	} while (!feof(stateprob));

	fclose(stateprob);


	//Loop over the adherence-stratified OnCare lookup tables
	for (int filenum = 0; filenum < NUM_ADHERENCE-1; filenum++) 
	{
		filepath = on_care_root + IntToStr(filenum) + nameend;

		//Open the file
		if((stateprob = fopen( filepath.c_str(), "r" )) == NULL )
		{
			printf("ERROR opening %s\n", filepath.c_str());
			exit(0);
		}
		
		// We need to ignore the first line of the file. That line is a human readable header.
		fgets(header_line, 1000, stateprob); 
		
		//Read in the onCare model rates from that file
		do {
			fscanf(stateprob,"%d%d%d%d%d%d%d%d%d%d%lf%f", &curr_cd4, &curr_vl, &curr_onTreat, &curr_resistance, &dest_cd4, &dest_vl, &dest_onTreat, &dest_resistance, &dead, &thecount, &therate, &cost);

			HIVprograte[filenum/*adherence level*/][1/* on care */][curr_cd4][curr_vl][curr_onTreat][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][dead] = (double)therate;	
			costCareTreat[filenum/*adherence level*/][1/* on care */][curr_cd4][curr_vl][curr_onTreat][curr_resistance] = (double)cost;

		} while (!feof(stateprob));

		
		fclose(stateprob);

		cout << "Reading on care rate file: "<< filepath.c_str() << endl;

	}

	// We need to be sure that every compartment with living people has at least one rate from that starting compartment. 
	//This function will report an error and quit if this is not the case.

	check_HIVprogrates();

	//For testing, prints the rates to a file
	//show_rates();

	//Apply the HIV mortality rate multiplier
	for (int onCare = 0; onCare < CARESTATUS; onCare++) {
		for (curr_cd4 = 0; curr_cd4 < NUM_CD4; curr_cd4++) {
			for (curr_vl = 0; curr_vl < NUM_VL; curr_vl++) {
				for (curr_resistance = 0 ; curr_resistance < NUM_RES; curr_resistance++) {
					for (dest_cd4 = 0; dest_cd4 < NUM_CD4; dest_cd4++) {
						for (dest_vl = 0; dest_vl < NUM_VL; dest_vl++) {
							for (dest_onTreat = 0; dest_onTreat < NUM_TREATORNOT; dest_onTreat++) {
								for (dest_resistance = 0 ; dest_resistance < NUM_RES; dest_resistance++) {
									for(int adh = 0; adh < NUM_ADHERENCE; adh++){
										//We found during calibration the the off treatment deaths peaked too early and therefore created two different 
										//HIV mortality multipliers for pre and post ART.  
										HIVprograte[adh][onCare][curr_cd4][curr_vl][0/*Off Treat*/][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][1/*dead*/] *= HIV_MORT_MULT_OFF_TREAT;
										HIVprograte[adh][onCare][curr_cd4][curr_vl][1/*On Treat*/][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][1/*dead*/] *= HIV_MORT_MULT_ON_TREAT;

									}
								}
							}
						}
					}
				}
			}
		}
	}

	//For testing, prints the rates to a file
	//show_rates();
}




string IntToStr(int n)
{
	std::ostringstream result;
	result << n;
	return result.str();
}

string DoubleToStr(double n)
{
	std::ostringstream result;
	result << n;
	return result.str();
}

void show_rates()
{
	int Adherence, curr_cd4 = 0, curr_vl = 0, curr_onTreat = 0, curr_resistance = 0;	
	int dest_cd4 = 0, dest_vl = 1, dest_onTreat = 1, dest_resistance = 0, dead = 0;
	int oncare;
	string filename;
	ofstream testoutfile;
	testoutfile.precision(10);

	for (Adherence = 0; Adherence < NUM_ADHERENCE; Adherence++) {
		filename = "TEST_Adherence_" + IntToStr(Adherence) + "_output.txt";
		testoutfile.open(filename.c_str());
		for (oncare = 0; oncare < 2; oncare++) {
			for (curr_cd4 = 0; curr_cd4 < NUM_CD4; curr_cd4++) {
				for (curr_vl = 0; curr_vl < NUM_VL; curr_vl++) {
					for (curr_onTreat = 0; curr_onTreat < 2; curr_onTreat++) {
						for (curr_resistance = 0; curr_resistance < NUM_RES; curr_resistance++) {
							for (dest_cd4 = 0; dest_cd4 < NUM_CD4; dest_cd4++) {
								for (dest_vl = 0; dest_vl < NUM_VL; dest_vl++) {
									for (dest_onTreat = 0; dest_onTreat < 2; dest_onTreat++) {
										for (dest_resistance = 0; dest_resistance < NUM_RES; dest_resistance++) {
											for (dead = 0; dead < 2; dead++) {
												if (HIVprograte[Adherence][oncare][curr_cd4][curr_vl][curr_onTreat][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][dead] != 0) {
													testoutfile << curr_cd4 << "\t" << curr_vl << "\t" << curr_onTreat << "\t" << curr_resistance << "\t" << dest_cd4 << "\t" << dest_vl << "\t" << dest_onTreat << "\t" << dest_resistance << "\t" << dead << "\t\t" << HIVprograte[Adherence][oncare][curr_cd4][curr_vl][curr_onTreat][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][dead] << endl;

												}
											}
										}
									}
								}
							}
						}
					}
				}				
			}
		}
		testoutfile.close();
	}
	//exit(11);
}

void printNonZeroCompartments(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES])
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	cout<<"t="<<t<<endl;
	cout<<"status\tsex\tact\tage\tvl\tcd4\tres\n\n";

	LOOP_THROUGH_ALL_COMPARTMENTS
		if (!isEqual(myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res],0))
		{
			printf ("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", status, sex, orient, act, age, idu, vl, cd4, res);
			cout << myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res]<<endl;
		}
	END_LOOP_THROUGH_ALL_COMPARTMENTS
	cout<<endl;
}

void calculateCD4andVLratesNewInfecteds()
{
	int cd4, vl, sex;
	double mean_cd4, mean_cd4sd, mean_vl, mean_vlsd;
	double cd4probs[NUM_CD4], vlprobs[NUM_VL];	//the probability of going into each cd4 or vl bucket

	for (sex=0; sex<NUM_SEX; sex++)
	{
		if (sex==M)
		{
			mean_cd4 = CD4MEAN_M;
			mean_cd4sd = CD4SD_M;
			mean_vl = VLMEAN_M;
			mean_vlsd = VLSD_M;
		}
		else //(SEX==F)
		{
			mean_cd4 = CD4MEAN_F;
			mean_cd4sd = CD4SD_F;
			mean_vl = VLMEAN_F;
			mean_vlsd = VLSD_F;
		}

		if (mean_cd4sd == 0) mean_cd4sd = 0.0001;			//Need a non-zero CD4SD
		if (mean_vlsd == 0) mean_vlsd = 0.000001;			//Need a non-zero VLSD

		 //change values for sensitivity analysis
#if ONE_WAY_SENSITIVITY_ANALYSIS
		if (sex == F  && (RUN_NUM == 63 || RUN_NUM == 64)) 
		{
				mean_cd4 = sensitivity_values[RUN_NUM]; 
				mean_cd4sd = sensitivity_values[RUN_NUM+106]; //run number 169/170 were used for cd4_std_F
		}
		else if (sex == F && (RUN_NUM == 65 || RUN_NUM == 66)) mean_vl = sensitivity_values[RUN_NUM]; 
		else if (sex == M  && (RUN_NUM == 67 || RUN_NUM == 68)) 
		{
				mean_cd4 = sensitivity_values[RUN_NUM]; 
				mean_cd4sd = sensitivity_values[RUN_NUM+104]; //run number 171/172  were used for cd4_std_M
		}
		else if (sex == M && (RUN_NUM == 69 || RUN_NUM == 70)) mean_vl = sensitivity_values[RUN_NUM];
#endif 

#if PROBABILISTIC_RUN

			mean_cd4 = prob_pull[31]; 
			mean_vl = prob_pull[33];
			mean_cd4sd = prob_pull[32];

#endif 
		//Calculate the probability of entering each CD4 bucket.  Note that the first and last buckets are special cases.
		cd4probs[0] = normal_cdf((cd4pts[0] - mean_cd4)/mean_cd4sd);

		for (cd4 = 1; cd4 < NUM_CD4-1; cd4++) {
			cd4probs[cd4] = normal_cdf((cd4pts[cd4] - mean_cd4)/mean_cd4sd) - normal_cdf((cd4pts[cd4-1]-mean_cd4)/mean_cd4sd) ;
		}

		cd4probs[NUM_CD4 - 1] = 1 -  normal_cdf((cd4pts[cd4-1]-mean_cd4)/mean_cd4sd);

		//Calculate the probability of entering each VL bucket.  Note that the first and last buckets are special cases.

		vlprobs[0] = normal_cdf((vlpts[0] - mean_vl)/mean_vlsd);

		for (vl = 1; vl < NUM_VL-1; vl++) 
		{
			vlprobs[vl] = normal_cdf((vlpts[vl] - mean_vl)/mean_vlsd) - normal_cdf((vlpts[vl-1] - mean_vl)/mean_vlsd);
		}
		vlprobs[vl] = 1-normal_cdf((vlpts[vl-1] - mean_vl)/mean_vlsd);

		//Now multiply the cd4 and vl probabilities to get the probability of entering each VL/CD4 combination category

		if (debug) cout<<"Rates of entering each VL/CD4 state upon initial infection: cd4 is down, vl is across"<<endl;

		for (cd4 = 0; cd4 < NUM_CD4; cd4++) 
		{
			for (vl = 0; vl < NUM_VL; vl++) 
			{
				probcd4vl[sex][cd4][vl] = cd4probs[cd4] * vlprobs[vl];
				if (debug) cout << probcd4vl[sex][cd4][vl] <<"\t";
			}
			if (debug) cout << endl;
		}
		if (debug) cout << endl<<endl;
	}
}

//Taken from:  http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
double normal_cdf( double x)
{
	const double b1 =  0.319381530;
	const double b2 = -0.356563782;
	const double b3 =  1.781477937;
	const double b4 = -1.821255978;
	const double b5 =  1.330274429;
	const double p  =  0.2316419;
	const double c  =  0.39894228;

	if(x >= 0.0) {
		double t = 1.0 / ( 1.0 + p * x );
		return (1.0 - c * exp( -x * x / 2.0 ) * t *
			( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
	}
	else {
		double t = 1.0 / ( 1.0 - p * x );
		return ( c * exp( -x * x / 2.0 ) * t *
			( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
	}
}

double compromise (double a, double b)
{
	if (COMPROMISE == MINIMUM) return min(a,b);

	else if (COMPROMISE == MAXIMUM) return max(a,b);

	else if (COMPROMISE == MEAN) return ((a+b)/2.0);

	cout << "COMPROMISE does not equal MINIMUM, MAXIMUM, or MEAN!\n";
	exit(1);
}

void set_age_change_rates(int age_of_sexual_debut)
{
	for (int age = 0; age < NUM_AGE; age++) {
		if (AGE_CHANGE_RATE > 0)			
		{
			switch (age) 
			{
			case 2:
				age_out_rate[age] = 1.0 / (age_of_sexual_debut - 10);
				break;
			case 3:
				age_out_rate[age] = 1.0 / (20 - age_of_sexual_debut);
				break;
			default:
				age_out_rate[age] = AGE_CHANGE_RATE; 
				break;
			}
		}
		else age_out_rate[age] = 0;
	}
}

//If there is an abstinent class, we want to scale up the fertility rates 
//so that the proper fertility is applied to the population
void scaleUpFertilityRates()
{
	double fertilityRate;

	for (int age = ADULTHOOD; age < NUM_AGE; age++) 
	{
		//We put in an array of fertility rates from 1997 to year 2014.  When we are calibrating, we should
		//use the rate from the appropriate year.  When we are running for real, we should use the 2010 value.
		if (CALIBRATE && int(t) < NUM_YEARS_CALIBRATION) fertilityRate = fertRate[int(t)][age];
		else fertilityRate = fertRate[NUM_YEARS_CALIBRATION-1][age];

		//We only do this if there is an abstinent class
		if (FIRST_ACTIVE_CLASS!=0)
		{
			fertRateScaledUp[age] = fertilityRate * 
				(numAbstinentWomen[age] + numSexuallyActiveWomen[age]) / numSexuallyActiveWomen[age];
		}

		else	//no scaling up necessary, there is no abstinent class
		{
			fertRateScaledUp[age] = fertilityRate;
		}
	}
}


//UPDATE THIS!!!
/*The indices for the baselineRR array.  
0:  Overall
1:  Female versus male, RR
2:  Same versus opposite, RR
3:  Both versus opposite, RR
4:  Risk "0" versus Risk "1", RR
5:  Risk "2" versus Risk "1", RR
6:  IDU "1" versus IDU "0", RR
7:  HIV+ infected, known versus HIV-, RR
8:  HIV+ infected, unknown versus HIV-, RR
9:	HIV+ infected, known, treated versus HIV-, RR
*/

/* The pathways
0: More condom use						(mult on sex alpha)
1: Fewer partners						(changes distribution of initial population in sexual risk groups)
2: Less IDU								(changes proportion of population in IDU group)
- prop in idu group					Note:  proportion in idu category 1 should be Path2*(1-Path3)
3: Less risk per IDU act
- prop who bleach needles			(changes proportion of population in IDU group)
4: Earlier ARV							(changes rate to care)
- rate of being tested for HIV		Note:  rate to care = Path4*(1-Path5)  
5: Earlier ARV							(changes rate to care)
- probability of failing to be linked to care once diagnosed	
6: More effective ARV 					(changes lookup table)
- prob of nonadherence
7: Prophylactic ARV						(mult on sex alpha)
8: STD tx								(mult on sex alpha)
9: microbicides							(mult on sex alpha)
10:circumcision							(mult on sex alpha)
*/

void initialize_interventions() 
{
	int n, p;
	int status, sex, orient, act, alc, idu, vl;

	double proportionOfWomenWhoGiveBirthInAYear = getProportionOfWomenWhoGiveBirthInAYear();


	//Initialize the interventions, set the cost, the "who" and the effect sizes
	for (n = 0; n < NUM_INTERVENTIONS; n++) 
	{	
		//Initialize each intervention and assign the intervention costs
		interv[n].initialize();	
		interv[n].overall_effect_size = overall_effect_size[n];
		interv[n].cost_per_person = interv_costs_per_person[n];

		//Fill in the effect sizes for the intervention/pathway combinations
		for (p=0; p<NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE; p++)
		{
			interv[n].effect_size[p]=interv_effect_size[n][p];

			//Sanity check to make sure all effect sizes are between 0 and 1.
			if (interv[n].effect_size[p] > 1 || interv[n].effect_size[p] < 0)
			{
				cout<<endl<<"The effect size of intervention "<<n<<", pathway "<<p<<" is NOT between 0 and 1.  The value is "<<interv[n].effect_size[p]<<".  Exiting."<<endl;
				exit(0);
			}
		}

		//Fill in the "Who" array for who the intervention applies to
		//Note that at this time (9/29/2014) VL targeting only works for the alcohol intervention (further code changes are needed to apply it to other pathways).
		LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS

			switch (n)
			{
				//interventions are in the order of spreadsheet
				//Currently all interventions apply only to adults
				
				case 0: if (sex == F && act == 3 ) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1;	break;					

				// Applies to adults who are not on Care or Treat (Note, the CD4 threshold portion 
				//	of this intervention applies to ALL, but the rest applies only to those not on 
				//	Care or Treat and that is how we apply cost to this intervention.)
				case 1: if(status < CARE ) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				// Applies to all Pregnant women
				case 2: if (sex==F) interv[n].who[status][sex][orient][act][alc][idu][vl] = proportionOfWomenWhoGiveBirthInAYear;  break;
				
				// SBIRT Applies to all alc adults
				case 3:  if (alc==ALC) 
					if (status >= SBIRT_TARGET_START && status <= SBIRT_TARGET_STOP && vl >= SBIRT_VL_TARGET_THRESHOLD) 
						interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; 
					break;  
				
				// Applies to all HIV- adult men
				case 4: if(status == SUSC && sex == M ) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1;  break;

				// Applies to all adults
				case 5: interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				// Applies to all adults
				case 6: interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				// PrEP All HIV negative adults in sex risk class 2+3 (classes 0+1 excluded)
				case 7: if (status == SUSC && act > 1 ) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				// Care+ Applies to all known HIV positive
				case 8: if (status >= DET && status <= TREAT) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Alcohol Individual Short: applies only to Male ALCS with HIV on TREAT  
				case 9: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc == ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break; 
				
				//Alcohol Individual Long: applies only to Male ALCS with HIV on TREAT  
				case 10: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc == ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Alcohol Group Short: applies only to Male ALCS with HIV on TREAT   
				case 11: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc == ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Alcohol Group Long: applies only to Male ALCS with HIV on TREAT   
				case 12: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc == ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Alcohol Community: cost applies to everyone, intervention applies only to ALCS (fixed cost for entire population)
				case 13: interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Depression Individual Short: applies only to Male ALCS with HIV on TREAT  
				case 14: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = PROP_ALC_DEPRESSED; break;

				//Depression Individual Long: applies only to Male ALCS with HIV on TREAT 
				case 15: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = PROP_ALC_DEPRESSED; break;

				//Depression Group Short: applies only to Male ALCS with HIV on TREAT  
				case 16: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = PROP_ALC_DEPRESSED; break;

				//Depression Group Long: applies only to Male ALCS with HIV on TREAT  
				case 17: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = PROP_ALC_DEPRESSED; break;

				//Sex Individual Short: applies only to Male ALCS with HIV on TREAT 
				case 18: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Sex Individual Long: applies only to Male ALCS with HIV on TREAT 
				case 19: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Sex Individual Short: applies only to Male ALCS with HIV on TREAT 
				case 20: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Sex Individual Long: applies only to Male ALCS with HIV on TREAT 
				case 21: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Sex Community: applies to all adults (fixed cost for entire population)
				case 22:  interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Weekly SMS adherence   
				case 23: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Brief adherence counseling
				case 24: if (status >= INDIA_ALC_TARGET_START && status <= INDIA_ALC_TARGET_STOP && alc >= ALC && sex == M) interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;

				//Mass media HIV education campaigns 
				case 25:  interv[n].who[status][sex][orient][act][alc][idu][vl] = 1; break;
			}
	
		END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS
	}

	//Move the high concurrency people and IDU users (just once!) if are applying those interventions.  
	//This function will need to be rewritten once we know how conc reduction will work
	//Apply_Conc_and_IDU_Interventions();
}



void createBasketOutput()
{
	static bool firstpass=true;
	string the_interventions;
	
	if (firstpass)
	{
#if ALCOHOL_SURFACE
		basketfile<<"Surface run\tCost per person\tEffect size\tRR Mult\tCondom RR\tAdherence RR\tSTI RR\t";
#endif

#if	OPTIMIZING_ALCOHOL_INTERVENTION
		basketfile<<"baseline run?\tIntervention effect size\tproportion Males w/alcohol issues\tproportion Females w/alcohol issues\t";
#endif

#if ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION || OPTIMIZING_ALCOHOL_INTERVENTION
		basketfile<<"alcohol RR of condom nonuse\talcohol RR of nonadherence\talcohol RR of STIs\tAlcohol Intervention turned on\t";
#endif

		if (RUNNING_MULTIPLE_BASKETS) 
		{
			basketfile<<"Num new infections with no interventions: "<<num_infections_no_interv<<endl;
			basketfile<<"Num dead of AIDS with no interventions: "<<totalDeadAIDS_no_interv<<endl;
		}


#ifdef RUN_OPTIMISTIC_COMBOS_FOR_SCOTT

		basketfile<<"Effect sizes equal 0\t";
#endif

#ifdef RUN_JASON_B
		basketfile<<"Hazard Strata\t";
		basketfile<<"Baseline testing frequency\t";
#endif

#if CARE_PLUS
	basketfile <<"CD4 treatment threshold\tCondom Path Effect Size\tLinkage Path Effect Size\tAdherence Path Effect Size\tIntervention Cost\t";
#endif

#if ONE_WAY_SENSITIVITY_ANALYSIS
		basketfile<<"Sensitivity name\tSensitivity Value\tRun Number\t";
#endif		
		
		basketfile<<"Intervention(s)\t";
		
		basketfile<<"Total basket cost Prepurchased\tTotal discounted basket cost Prepurchased\t";
		basketfile<<"Total basket cost Utilization\tTotal discounted basket cost Utilization\t";
		basketfile<<"Total basket cost mixing cost types\tTotal discounted basket cost mixing cost types\t";
		basketfile<<"Cost of Care and Treatment\tDiscounted Cost of Care and Treatment\t";

		basketfile<<"Total Life Years (All)\tTotal Discounted Life Years (All)\tTotal QALYs (All)\tTotal Discounted QALYs (All)\t";
		basketfile<<"Total Life Years (HIV+)\tTotal Discounted Life Years (HIV+)\tTotal QALYs (HIV+)\tTotal Discounted QALYs (HIV+)\t";
		basketfile<<"Total Life Years (Care+Treat)\tTotal Discounted Life Years (Care+Treat)\tTotal QALYs (Care+Treat)\tTotal Discounted QALYs (Care+Treat)\t";
		

		basketfile<<"Num new infections over "<<t<<" years\t";
		basketfile<<"Num infections averted over "<<t<<" years"<<"\t";

		basketfile<<"Num died of AIDS over "<<t<<" years"<<"\t";

		basketfile<<"Mean num new infections per infected per year\t";
		basketfile<<"Mean num HIV deaths per infected per year\t";
		basketfile<<"Cost of care and treatment per infected per year\t";

#if	OPTIMIZING_ALCOHOL_INTERVENTION
		basketfile<<"Proportion new infections alcohol related\t";
#endif

		basketfile<<"Proportion Detected on treatment\tProportion HIV+ on treatment\tProportion HIV+ detected\t";

		for (int i=0; i<NUM_INTERVENTIONS; i++) basketfile<<"Int "<<i<<" Pre\tInt "<<i<<" Ut\t"; 

#if PROBABILISTIC_RUN /*output all of the probabilistic run info*/
		for (int i = 0; i < NUM_PROB_TESTS + NUM_PROB_COST_EFFECT; i++)
		{
			basketfile << probabilistic_name[i] << "\t";
		}

#endif 

		basketfile<<endl;
		
		firstpass = false;
	}

#if ALCOHOL_SURFACE
	basketfile<<cost_index<<"_"<<effect_index<<"_"<<rr_index <<"\t"<< interv_costs_per_person[3] <<"\t" << interv_effect_size[3][10] << "\t" << alcohol_rr_modifier[rr_index] <<"\t"<< condomRR <<"\t" << adhereRR<<"\t"<< stiRR << "\t";// "Cost per person for alcohol intervention.\tAlcohol effect size on pathway 10, the pathway that reduces substance abuse.\tRelative Risk 0 - condoms\t3 - adherence\t5 - Untreated STI\t";
#endif
	
#if	OPTIMIZING_ALCOHOL_INTERVENTION
	basketfile<< (is_baseline_run ? "Yes":"No") << "\t" << interv_effect_size[3][10] << "\t" << prop_alcohol[M] <<"\t"<<prop_alcohol[F] <<"\t";
#endif

#if ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION
	basketfile<<baselineRR_raw[0][5]<<"\t"<<baselineRR_raw[3][5]<<"\t"<<baselineRR_raw[5][5]<<"\t"<<(basket[3]? "YES":"no")<<"\t";
#endif

#if OPTIMIZING_ALCOHOL_INTERVENTION
	basketfile<<condomRR<<"\t"<<adhereRR<<"\t"<<stiRR<<"\t"<<(basket[3]? "YES":"no")<<"\t";
#endif

#ifdef RUN_JASON_B
	basketfile<<hazard<<"\t";
	basketfile<<baselineRR_raw[4][0]<<"\t";
#endif

#if CARE_PLUS
	basketfile<<cd4_treat_threshold<<"\t";
	basketfile<<interv_effect_size[7][0]<<"\t"<<interv_effect_size[7][2]<<"\t"<<interv_effect_size[7][3]<<"\t";
	basketfile<<interv_costs_per_person[7]<<"\t";
#endif


#if ONE_WAY_SENSITIVITY_ANALYSIS
	basketfile << sensitivity_name[RUN_NUM] << "\t" << sensitivity_values[RUN_NUM] << "\t" << RUN_NUM << "\t";
#endif


	string used_interventions;
	used_interventions = "";
	
	int numintvs_used = 0, last_used = 0;
	for (int intv = 0; intv < NUM_INTERVENTIONS; intv++)
	{
		if (basket[intv])
		{
			last_used = intv;
			used_interventions += IntToStr(intv) + ": ";
			numintvs_used++;
		}
	}
	
	if (numintvs_used == 1)
	{
		used_interventions += intervention_name[last_used];
	} else if(numintvs_used == 0)
	{
		used_interventions += "none";
	}
	
	basketfile << used_interventions << "\t";

	basketfile<<getBasketCost(PP, NORM)<<"\t"<<getBasketCost(PP, DISC)<<"\t";
	basketfile<<getBasketCost(UT, NORM)<<"\t"<<getBasketCost(UT, DISC)<<"\t";
	basketfile<<getBasketCost(COMBO, NORM)<<"\t"<<getBasketCost(COMBO, DISC)<<"\t";
	basketfile<<costOfCareAndTreat<<"\t"<<discounted_costOfCareAndTreat<<"\t";
	
	for (int group=0; group<3; group++)
	{
		basketfile<<total_life_years[group]<<"\t"<<total_disc_life_years[group]<<"\t"<<total_qa_life_years[group]<<"\t"<<total_disc_qa_life_years[group]<<"\t";
	}

#if OPTIMIZING_ALCOHOL_INTERVENTION
	if (is_baseline_run) num_infections_no_interv = num_infections[basket_num];
#endif

	if (atLeastOneInterventionTurnedOn)
	{
		basketfile<<num_infections[basket_num]<<"\t";
		basketfile<<num_infections_no_interv - num_infections[basket_num]<<"\t";
	}

	if (!atLeastOneInterventionTurnedOn)
	{
		basketfile<<num_infections_no_interv<<"\t";
		basketfile<<"N/A\t";
	}

	basketfile<<getTotalDeadAIDS()<<"\t";
	basketfile <<sumOfNumInfectionsPerInfected / STOPTIME<<"\t";
	basketfile <<sumOfNumDeathsPerInfected / STOPTIME<<"\t";
	basketfile <<sumOfCostOfCareAndTreatPerInfected / STOPTIME<<"\t";

#if OPTIMIZING_ALCOHOL_INTERVENTION
	basketfile<<num_new_infections_alcohol/num_infections[basket_num]<<"\t";
#endif

	basketfile<<prop_detected_on_treatment<<"\t"<<prop_infected_on_treatment<<"\t"<<prop_infected_detected<<"\t";

	for (int i=0; i<NUM_INTERVENTIONS; i++) 
	{
		if (basket[i]) basketfile<<interv[i].total_prepurchased_cost<< "\t" << interv[i].total_utilization_cost << "\t";
		else basketfile<<"-\t-\t";
	}

#if SHOW_EFFECT_SIZES
	the_interventions.clear();
	for (int bsk = 0; bsk < NUM_INTERVENTIONS; bsk++)
	{
		basketfile<<intervention_name[bsk] <<"\t";
		for (int thepath = 0; thepath < NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE; thepath++)
		{
			if (basket[bsk])
			{
				basketfile << interv_effect_size[bsk][thepath]<<"\t";
			} else {
				basketfile<<"-\t";
			}
		}
	}
#endif

#if PROBABILISTIC_RUN /*output all of the probabilistic run info*/
	for (int i = 0; i < NUM_PROB_TESTS + NUM_PROB_COST_EFFECT; i++)
	{
		basketfile << prob_pull[i] << "\t"; 
	}

#endif 


	basketfile<<endl;
}


void setBaselineProbs()
{
	int p;
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	double popInCombo, totalPop;
	double rr[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];
	double sumOfNumxRR, x;
	bool setThisOneManually;

	for (p = 0; p < NUM_PATH_AND_CHAR; p++)
	{
		sumOfNumxRR=0;
		totalPop=0;
		x=0;			//While x doesn't need intiialization (it's not cumulative), we get a "potential uninitialized variable error" if it's not, because there are whole statuses of people during calibration that contain no people.

		for (status=0; status<NUM_ALIVE_STATE; status++)
		{
			for (sex=0; sex<NUM_SEX; sex++)
			{
				for (orient=0; orient<NUM_ORIENT; orient++)
				{
					for (act=0; act<NUM_ACT; act++)
					{
						for (alc=0; alc<NUM_ALC; alc++)
						{
							for (idu=0; idu<NUM_IDU; idu++)
							{
								popInCombo=0;

								//See comment below
								setThisOneManually=false;		

								//These are special cases where the reference compartment is actually the one that is N/A
								//So for instance, for the linkage and adherence pathways (2 and 3), they shouldn't apply to Susceptibles
								//And the microbicides pathway shouldn't apply to men
								//For these comparments, we set them to N/A "manually", hence the variable above
								if (!((p==2 && status==SUSC) || (p==3 && status==SUSC) || (p==6 && sex==M)))
								{
									for (age=ADULTHOOD; age<NUM_AGE; age++)
									{
										for (vl = 0; vl < NUM_VL; vl++)
										{
											for (cd4 = 0; cd4 < NUM_CD4; cd4++)
											{
												for (res = 0; res < NUM_RES; res++)
												{
													popInCombo += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
												}
											}
										}
									}
								}


								//Reference characteristic needs to be N/A
								else setThisOneManually=true;

								if (setThisOneManually) rr[status][sex][orient][act][alc][idu] = NOT_APPLICABLE;

								else rr[status][sex][orient][act][alc][idu] = getBaselineRRproduct(p, status, sex, orient, act, alc, idu);

								if (int(rr[status][sex][orient][act][alc][idu]) != NOT_APPLICABLE)
								{
									sumOfNumxRR += popInCombo * rr[status][sex][orient][act][alc][idu];

									totalPop += popInCombo;
								}
							}
						}
					}
				}
			}
		}


		if (sumOfNumxRR > 0) x = baselineRR_raw[p][0] * totalPop / sumOfNumxRR;
		//We thought this line would prevent an error, but it prevent calibration from working where there are no people to begin with in some of the compartments.  else {cout<<"x was not greater than sumOfNumxRR.  Exiting.";  exit(0);}

		for (status = 0; status<NUM_ALIVE_STATE; status++){
			for (sex = 0; sex<NUM_SEX; sex++){
				for (orient = 0; orient<NUM_ORIENT; orient++){
					for (act = 0; act<NUM_ACT; act++){
						for (alc = 0; alc<NUM_ALC; alc++){
							for (idu = 0; idu<NUM_IDU; idu++){

								//Just added the "|| sumOfNumxRR==0".  We just ran into a situation where the only state applicable in the pathway
								//had nobody in it.  Thus, calculating x gave a divide by zero error.  If there isn't anybody available
								//with which to calculate the baseline probabililties for that pathway (for example, nobody is detected at the start of the model
								//so the pathway for Probability of failing to be linked to care once diagnosed has no people to apply it to)
								//then set the baseline RR to be NOT_APPLICABLE
								if (rr[status][sex][orient][act][alc][idu] == NOT_APPLICABLE  || sumOfNumxRR==0)
								{
									baselineProb[p][status][sex][orient][act][alc][idu] = NOT_APPLICABLE;
								}

								else 
								{
									baselineProb[p][status][sex][orient][act][alc][idu] = rr[status][sex][orient][act][alc][idu] * x;

									//Make sure that the baseline RR isn't greater than 1.0.
									if (baselineProb[p][status][sex][orient][act][alc][idu] > 1) baselineProb[p][status][sex][orient][act][alc][idu] = 1.0;
								}

								//TEST
								//if (p==3 /*&& status==SUSC && sex==1 && orient==0 && risk==2 && idu==0 && age==1*/)
								//if (p==0 || p==3 || p==5 && baselineProb[p][status][sex][act][alc] != NOT_APPLICABLE) 
								//cout<<"path "<<p<<" status "<<status<<" sex "<<sex<<" act "<<act<<" alc "<<alc<<"   "<<baselineProb[p][status][sex][act][alc]<<endl;

							}
						}
					}
				}
			}
		}
	}


	//Check
	/*for (p = 0; p < NUM_PATH_AND_CHAR; p++)
	{
	double checkTotal=0;
	double totalPeople=0;

	for (status=0; status<NUM_ALIVE_STATE; status++)
	{
	for (sex=0; sex<NUM_SEX; sex++)
	{
	for (act=0; act<NUM_ACT; act++)
	{
	for (alc=0; alc<NUM_ALC; alc++)
	{
	popInCombo=0;

	for (age=ADULTHOOD; age<NUM_AGE; age++)
	{
	for (vl=0; vl<NUM_VL; vl++)
	{
	for (cd4=0; cd4<NUM_CD4; cd4++)
	{
	for (res=0; res<NUM_RES; res++)
	{
	popInCombo += Comp[status][sex][orient][act][age][alc][vl][cd4][res];
	}
	}
	}
	}

	//TEST if (getBaselineRRproduct(p, status, sex, act, alc) != NOT_APPLICABLE)
	if (baselineProb[p][status][sex][act][alc] != NOT_APPLICABLE)
	{
	checkTotal += popInCombo * baselineProb[p][status][sex][act][alc];
	totalPeople += popInCombo;
	}			
	}
	}
	}
	}

	cout<<"Recalculated baseline probability for pathway"<<p<<": "<<checkTotal/totalPeople<<endl;
	}*/

	/*p=2;
	for (status=DET; status<NUM_ALIVE_STATE; status++)
	{
	for (sex=0; sex<NUM_SEX; sex++)
	{
	for (orient=0; orient<NUM_ORIENT; orient++)
	{
	for (act=0; act<NUM_ACT; act++)
	{
	for (alc=0; alc<NUM_ALC; alc++)
	{
	for (vl=0; vl<NUM_VL; vl++)
	{
	for (idu=0; idu<NUM_IDU; idu++)
	{



	if (baselineProb[p][status][sex][orient][act][alc][idu] != NOT_APPLICABLE) 
	{
	//if (status==DET && sex==F && orient==STRAIGHT && act==2 && alc == 0 && vl==0 && idu==0)
	cout<<"path "<<p<<" status "<<status<<" sex "<<sex<<" orient "<<orient<<" act "<<act<<" alc "<<alc<<" idu "<<idu<<"   "<<baselineProb[p][status][sex][orient][act][alc][idu]<<endl;
	}
	END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS*/


}

double convertProbabilityToRate(double prob)
{
	double rate;

	if (prob < 0 || prob > 1) 
	{
		cout<<"Error in convertProbabilityToRate().  Probability was not a value between 0 and 1.  Exiting."<<endl;
		exit(0);
	}

	else if (prob == 1) rate = RATE_WHEN_PROB_EQUALS_ONE;
	else if (prob == 0) rate = 0;		//Without this line, the rate calculated below resulted in a rate of -0, which seemed to behave ok, but looked weird.  So we added in this special case just to be safe.

	else rate = -log(1-prob)/1;			//the "1" is for one year

	return rate;
}


double getBaselineRRproduct(int p, int status, int sex, int orient, int act, int alc, int idu)
{
	double product = 1;
	for (int rr=1; rr<NUM_BASELINE_RR && product != NOT_APPLICABLE; rr++)
	{
		switch (rr)
		{
		case 1:  if (sex == F) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;			
		case 2:  if (orient == GAY)
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;			
		case 3:  if (orient == BI)
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;			
		case 4:  if (act == 0) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;			
		case 5:  if (act == 2) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;			
		case 6:  if (act == 3) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;			
		case 7:  if (alc == ALC || alc == ALC_TREATED_NOT_CURED) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;		
		case 8:  if (idu == IDU) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;		
		case 9:  if (status == DET) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;	
		case 10:  if (status == CARE) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;	
		case 11:  if (status==HVL || status==INF) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;
		case 12:  if (status==TREAT) 
				 {if (baselineRR_raw[p][rr] != NOT_APPLICABLE) product *= baselineRR_raw[p][rr]; else product = NOT_APPLICABLE;} break;		
		}
	}

	//TEST
	//if (p==4 && product !=1 && product !=-999) 
	//cout<<"stop"<<endl;
	return product;
}

//getBasketCost takes two arguments.  
//The first, utilization_or_prepurchased can be either PP (prepurchased),  UT (utilization), or COMBO which will use the
//costType array in intervention.h to determine which costType (PP or UT) to use for that intervention.
//The second, discounted_or_normal can be either DISC (discounted) or NORM (normal)
//The function returns the appropriate total cost for the basket of interventions

double getBasketCost(int utilization_or_prepurchased, bool discounted_or_normal)
{
	double total_basket_cost = 0;
	double cost;

	//Loop through all interventions
	for (int i=0; i<NUM_INTERVENTIONS; i++)
	{
		//if the intervention is turned on
		if (basket[i]==true) 
		{
			cost = 0;

			//if we want the prepurchased cost
			if (utilization_or_prepurchased == PP || (utilization_or_prepurchased == COMBO && costType[i]==PP))
			{
				switch ((int)discounted_or_normal) // Xcode displays a warning when using a boolean in a switch statement, so we caste to an int
				{
				case (NORM): cost = interv[i].total_prepurchased_cost;  break;
				case (DISC): cost = interv[i].discounted_total_prepurchased_cost;  break;
				}
			}

			//if we want the utilization cost
			else if (utilization_or_prepurchased == UT || (utilization_or_prepurchased == COMBO && costType[i]==UT))
			{
				switch ((int)discounted_or_normal) // Xcode displays a warning when using a boolean in a switch statement, so we caste to an int
				{
				case (NORM): cost = interv[i].total_utilization_cost;  break;
				case (DISC): cost = interv[i].discounted_total_utilization_cost;  break;
				}
			}

			//illegal value was passed
			else if (utilization_or_prepurchased > COMBO) 
			{
				cout<<"Illegal costType value.  Exiting.";  
				exit(0);
			}

			total_basket_cost += cost;
		}		
	}

	return total_basket_cost;
}

double getTotalDeadAIDS()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	//Add up the total number of AIDS deaths so far
	status=DEAD_HIV;
	for (sex=0; sex<NUM_SEX; sex++)	{
		for (orient=0; orient<NUM_ORIENT; orient++)	{
			for (act=0; act<NUM_ACT; act++)	{
				for (age=0; age<NUM_AGE; age++)	{
					for (alc = 0; alc<NUM_ALC; alc++){
						for (idu = 0; idu<NUM_IDU; idu++){
							for (vl=0; vl<NUM_VL; vl++){
								for (cd4=0; cd4<NUM_CD4; cd4++){
									for (res = 0; res < NUM_RES; res++)
									{
										totalDeadAIDS += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
									}}}}}}}}}

	return totalDeadAIDS;
}

void calcPropRes(double Pop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES])
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	double totalInf = 0;							//the total number of infected persons in model

	//initialize arrays to zero
	memcpy(propRes, zerod_array, sizeof(propRes));

	//determine proportion in each resistance category (use propRes to store total number in each category for now)
	LOOP_THROUGH_ALL_COMPARTMENTS

		//This used to say status>=INF, thought that was wrong, so changed it
		if (status>=HVL && status<=TREAT && age>=ADULTHOOD) 
		{
			propRes[res] += Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			totalInf += Pop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		}

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	//divide by total number infected and sexually active to get proportion in each resistance category
	//Only do this if the number of infectives is not zero to prevent a divide by zero error!
	if (totalInf) {for (res=0; res<NUM_RES; res++) propRes[res] /= totalInf;}
	//If there isn't anybody infected in the population, set all new infecteds to get no resistance (all go into the first resistance category which represents no resistance)
	else propRes[0]=1;			
}

//This function finds the rates from INF to Det and from Det to Care
void setBaseLineProbsTestedLinkedandAdherence()
{
	int status, sex, orient, act, alc, idu, vl;
	double final_effect_size, effect_size;
	double adjustedBaselineProb;
	int i, p;

	//Repeat for paths 1,2 and 3 (rates of being tested and of being linked to care and adherence)
	for (p=1; p<=3; p++)
	{	
		//Loops through the relative risk categories on the baseline probabilities
		LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS 

			//intialize final effect size to no effect
			final_effect_size = 1;

			//Find the cumulative effect over all interventions that activate that pathway
			for (i=0; i < NUM_INTERVENTIONS; i++)
			{
				//If an intervention is turned on that activates this pathway and the effect sizes is not 1 (implying there actually is an effect)
				if (basket[i] == true &&
					!(baselineProb[p][status][sex][orient][act][alc][idu] == NOT_APPLICABLE) &&
					!(isEqual(interv[i].effect_size[p], 1)))
				{

					//To find the actual effect size for an intervention, a who of '0' should have no effect (resulting in an effect size of 1)
					//A who of '1' applies the effect to everyone in the compartment, hence the full effect size applies
					//A who of some value in between gets partial effect, so interpolate!
					//linearInterpolate(x0, y0, x1, y1, x)  applied as....
					//linearInterpolate(who 0, effect size 1, who 1, full intervention effect size, who for this intervention

					effect_size = linearInterpolate(0, 1, 1, interv[i].effect_size[p], interv[i].who[status][sex][orient][act][alc][idu][vl]);

					if (MODE == OPTIMISTIC) final_effect_size *= effect_size;			//Optimistic

					else final_effect_size = min(final_effect_size, effect_size);		//Pessimistic

				}
			}

			if (baselineProb[p][status][sex][orient][act][alc][idu] != NOT_APPLICABLE)
			{
				adjustedBaselineProb = baselineProb[p][status][sex][orient][act][alc][idu] * final_effect_size;

				switch (p)
				{	
					//The rate of moving from INF to Detected (pathway 1 is prob of NOT being tested).  Convert to rate.
				case (1): 
					if (status==INF) 
					{
						if (adjustedBaselineProb < 0 || adjustedBaselineProb > 1)
						{
							cout<<"adjustedBaselineProb (case 1) is not a value between 0 and 1.  Exiting"<<endl;
							exit(0);
						}

						inf_to_det_rate[sex][orient][act][alc][idu] = convertProbabilityToRate(1 - adjustedBaselineProb); 
					}
					break;

					//The rate of moving from Detected to CARE (pathway 2 is prob of NOT being linked).  Convert to rate.
				case (2): 
					if (status==DET) 
					{
						if (adjustedBaselineProb < 0 || adjustedBaselineProb > 1)
						{
							cout<<"adjustedBaselineProb (case 2) is not a value between 0 and 1.  Exiting"<<endl;
							exit(0);
						}

						det_to_care_rate[sex][orient][act][alc][idu] = convertProbabilityToRate(1 - adjustedBaselineProb);
					}
					break;

					//The probability of being adherent (pathway 3). 
				case (3):  
					if (status==TREAT || status==CARE)			//(Remember, we need to use the adherence-stratified lookup tables if the status is TREAT or CARE
					{
						if (adjustedBaselineProb < 0 || adjustedBaselineProb > 1)
						{
							cout<<"adjustedBaselineProb (case 3) is not a value between 0 and 1.  Exiting"<<endl;
							exit(0);
						}

						adherence[sex][orient][act][alc][idu] = 1 - adjustedBaselineProb;
					}
					break;
				}
			}

		END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS
	}
}

double linearInterpolate(double x0, double y0, double x1, double y1, double x)
{
	double y;				//the return value

	if (x==x0) y = y0;	//prevent divide-by-zero error
	else y = y0 + (x-x0)*(y1-y0)/(x1-x0);

	return y;
}

//For reference:  HIVprograte[NUM_ADHERENCE][CARESTATUS][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_ALIVEORDEAD];		
double getAdhStratifiedHIVProgRate(int status, int curr_cd4, int curr_vl, int curr_res, int dest_cd4, int dest_vl, int destOnTreatment, int destRes, int destDead, int sex, int orient, int act, int alc, int idu)
{
	double retval;
	int progFloor=0;
	int progCeiling=0;
	double adhere;		//adherence, but "adherence[status][sex][act][alc]" is a global.  Couldn't reuse the whole word.

	/* Very important to know when accessing the progression rates.
	For people transitioning to dead of hiv, we set all of their destination variables (other than dead) to 0.

	If this modification has not been made at the time the function is called, we make sure the appropriate tweaks
	are made here.
	*/

	if (destDead==1) 
	{
		dest_cd4 = 0;
		dest_vl = 0;
		destOnTreatment = 0;
		destRes = 0;
	}


	/***********************end of special cases ************************/

	//If on care, make sure to use the appropriate adherence level in the array
	if (ONCARE(status))
	{
		adhere = adherence[sex][orient][act][alc][idu];

		if (adhere <0 || adhere > 1) {cout<<"Error!  Adherence is not a value between 0 and 1.  Exiting."<<endl;  exit(0);}

		//Note that people who aren't yet on treatment, don't have a level of adherence (it's NOT_APPLICABLE)
		//So we need to pick which on-care table we are going to use.  Theoretically, they should all have the
		//same level of downward progression for people who are not yet on treatment, so let's use the 
		//average adherence in the progression model
		//if (status==CARE || isEqual(adherence, NOT_APPLICABLE)) adherence = 0.63;
		//else adherence = 0.63;

		// Using a cast to integer instead of floor/ceiling to save time.  (Used to say progFloor = (int) floor( 10 * adhere); progCeiling = (int) ceil( 10 * adhere);)
		progFloor = (int)( 10 * adhere);
		progCeiling = min (progFloor+1, 10);	//In general, progCeiling will be one index higher, except in the case of perfect adherence - there are only 11 adherence files, 10 is the highest numbered file!

		retval = linearInterpolate(progFloor/10.0,
			HIVprograte[progFloor][1/*On Care*/][curr_cd4][curr_vl][ONTREAT(status)][curr_res][dest_cd4][dest_vl][destOnTreatment][destRes][destDead],
			progCeiling/10.0,
			HIVprograte[progCeiling][1/*On Care*/][curr_cd4][curr_vl][ONTREAT(status)][curr_res][dest_cd4][dest_vl][destOnTreatment][destRes][destDead], 
			adhere);

		//TEST
		//cout<<HIVprograte[progFloor][1/*On Care*/][curr_cd4][curr_vl][ONTREAT(status)][curr_res][dest_cd4][dest_vl][destOnTreatment][destRes][destDead];
	}

	//If not on care, make sure to use OFFCARE_INDEX
	else
		retval = HIVprograte[OFFCARE_INDEX][0/*Not on Care*/][curr_cd4][curr_vl][ONTREAT(status)][curr_res][dest_cd4][dest_vl][destOnTreatment][destRes][destDead];
	
	if (retval<0) {cout<<"**************ERROR Retval < 0!!!!***********************"<<endl; exit (0); }

	return retval;
}

//For reference:  costCareTreat[NUM_ADHERENCE][CARESTATUS][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES];		
double getAdhStratifiedCostCareTreat(int status, int curr_cd4, int curr_vl, int curr_res, int sex, int orient, int act, int alc, int idu)
{
	double retval;
	int progFloor=0;
	int progCeiling=0;
	double adhere;		//adherence, but "adherence[status][sex][act][alc]" is a global.  Couldn't reuse the whole word.

	//If on care, make sure to use the appropriate adherence level in the array
	if (ONCARE(status))
	{
		adhere = adherence[sex][orient][act][alc][idu];
		if (adhere <0 || adhere > 1) {cout<<"Error!  Adherence is not a value between 0 and 1.  Exiting."<<endl;  exit(0);}

		//Note that people who aren't yet on treatment, don't have a level of adherence (it's NOT_APPLICABLE)
		//So we need to pick which on-care table we are going to use.  Theoretically, they should all have the
		//same level of downward progression for people who are not yet on treatment, so let's use the 
		//average adherence in the progression model
		//if (status==CARE || isEqual(adherence, NOT_APPLICABLE)) adherence = 0.63;
		//else adherence = 0.63;

		// Using a cast to integer instead of floor/ceiling to save time.  (Used to say progFloor = (int) floor( 10 * adhere); progCeiling = (int) ceil( 10 * adhere);)
		progFloor = (int)( 10 * adhere);
		progCeiling = min (progFloor+1, 10);	//In general, progCeiling will be one index higher, except in the case of perfect adherence - there are only 11 adherence files, 10 is the highest numbered file!

		retval = linearInterpolate(progFloor/10.0,
			costCareTreat[progFloor][1/*On Care*/][curr_cd4][curr_vl][ONTREAT(status)][curr_res],
			progCeiling/10.0,
			costCareTreat[progCeiling][1/*On Care*/][curr_cd4][curr_vl][ONTREAT(status)][curr_res], 
			adhere);
	}

	//If not on care, make sure to use OFFCARE_INDEX
	else
		retval = costCareTreat[OFFCARE_INDEX][0/*Not on Care*/][curr_cd4][curr_vl][ONTREAT(status)][curr_res];
	
		if (retval<0) {cout<<"**************ERROR Retval < 0!!!!***********************"<<endl; exit (0); }

	return retval;
}


void setTotalPeopleinModel(double myPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES])
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	totalPeopleInModel=0, totalinfectedsinmodel=0;

	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS 

		totalPeopleInModel += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];

		if (status>SUSC && status<=TREAT)  totalinfectedsinmodel += myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		
		if (myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res] < 0) 
		{
			cout<<"Negative Pop!"<<myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res]<<endl;
			printf ("%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n", t, status, sex, orient, act, age, alc, idu, vl, cd4, res, myPop[status][sex][orient][act][age][alc][idu][vl][cd4][res]);
			exit(0);
		}

	END_LOOP_THROUGH_ALL_COMPARTMENTS 

}


//  Old 1 patient test function
/*void calculateCost()
{
	int sex, act, age, alc, vl, cd4, res;
	double costThisCycle;
	costOfCareAndTreatThisCycle = 0;	//initialize for use below
	double effect_size;					//the intervention effect size after taking into account the "who"
	double totalNumIncludingMoved;
	cout << "status\tcd4\tvl\tres\tsex\tact\talc\tCost\n";

	for (int res=0; res<NUM_RES; res++)
	{
		for (int status = INF; status<=TREAT; status++)
		{
			for (int cd4=0; cd4<NUM_CD4; cd4++)
			{
				for (int vl=0; vl < NUM_VL; vl++)
				{
					// ****************TEST*********************
					costOfCareAndTreat = 0;
					costOfCareAndTreatThisCycle = 0;

					sex=1;
					act=3;
					age=4;
					alc=0;
					Comp[status][sex][orient][act][age][alc][vl][cd4][res]=1;


					//Now calculate the cost of care and treatment using the data from the progression model
					//Note that include the cost of care and treatment for adults as well as children

					if (status>=INF && Comp[status][sex][orient][act][age][alc][vl][cd4][res] > 0) // change CARE to INF
					{

						//TEST
						// Need to loop over only source, not destination
						// Formula should be costOfCareAndTreatThisCycle += Comp[status][sex][orient][act][age][alc][vl][cd4][res] * preCalcAdhStratifiedCostCareTreat[status][cd4][vl][res][sex][act][alc] * DT

						if (EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS) {
							costOfCareAndTreatThisCycle += Comp[status][sex][orient][act][age][alc][vl][cd4][res] *
								preCalcAdhStratifiedCostCareTreat[status][cd4][vl][res][sex][act][alc] *
								DT;
						} 
						else
							costOfCareAndTreatThisCycle += Comp[status][sex][orient][act][age][alc][vl][cd4][res] *
							getAdhStratifiedCostCareTreat(status, cd4, vl, res, sex, act, alc) *
							DT;

						cout << status<< "\t" <<cd4<< "\t" <<vl<< "\t" <<res<< "\t" <<sex<< "\t" << act << "\t" << alc << "\t" << costOfCareAndTreatThisCycle << endl;


					}

					costOfCareAndTreat += costOfCareAndTreatThisCycle;
					discounted_costOfCareAndTreat += get_discounted_value(costOfCareAndTreatThisCycle);

				}}}}
}
*/

void calculateCost()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
    double costThisCycle;
	costOfCareAndTreatThisCycle = 0;	//initialize for use below
	double effect_size;					//the intervention effect size after taking into account the "who"
	double totalNumIncludingMoved;
	bool firstPass=true;				//Used for circumcision intervention below.


	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS

		for (int i=0; i<NUM_INTERVENTIONS; i++)
		{
			if (basket[i]==true)
			{
				//Only apply the cost of the interventions to adults
				if (age>=ADULTHOOD)
				{
					//Initialize here to prevent compile warnings (or inadvertently using the costThisCycle from the last cycle!)
					costThisCycle = 0;

					//The SBIRT intervention doesn't belong inside the if(who...) because
					//the people who have been treated have been moved out of the targeted
					//compartments and into ALC_TREATED_CURED and ALC_TREATED_NOT_CURED
					if (i==3 || (i>8 && i<13)) 
					{
						if (alc == ALC_TREATED_CURED || alc == ALC_TREATED_NOT_CURED)
						{
							costThisCycle = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res]
							* interv[i].cost_per_person
								* DT;
						}
					}

					//We do "else if" so that we don't do the SBIRT interv (interv 3 again here.)
					else if (interv[i].who[status][sex][orient][act][alc][idu][vl]>0)
					{
						/*Special case for intervention 0:  Condom distribution, risk reduction counseling, and STI treatment for female CSW
						We want to apply to cost to the number of women targeted in this group BEFORE some were moved to lower activity groups
						We don't have to modify the other interventions that move people because it just so happens that they all apply to all adults
						so the costing comes out the same.*/
						if (i==0)
						{
							//To find the actual effect size for an intervention, a who of '0' should have no effect (resulting in an effect size of 1)
							//A who of '1' applies the effect to everyone in the compartment, hence the full effect size applies
							//A who of some value in between gets partial effect, so interpolate!
							//linearInterpolate(x0, y0, x1, y1, x)  applied as....
							//linearInterpolate(who 0, effect size 1, who 1, full intervention effect size, who for this intervention

							effect_size = linearInterpolate(0, 1, 1, interv[i].effect_size[13], interv[i].who[status][sex][orient][act][alc][idu][vl]);

							//Find the total number of people in the compartment including those who were moved out due to interventions
							totalNumIncludingMoved = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] / effect_size;

							costThisCycle = totalNumIncludingMoved
								* interv[i].who[status][sex][orient][act][alc][idu][vl] 
							* interv[i].cost_per_person
								* DT;
						}

						//Special case for circumcision intervention (intervention #4)
						//We don't need to do this for every compartment, so it comes before the 
						//loops below. By this point, we know how many people have been circumsized by the intervention
						//so we only want to calculate the cost once.
						else if (i==4)
						{
							if (firstPass)		
							{
								//At time 0, we add the cost of circumcising all adult HIV- men, according to effect size
								if (t==0) costThisCycle = numCircumsizedByIntervAtStart * interv[i].cost_per_person;
								costThisCycle += numCircumsizedByIntervThisCycle * interv[i].cost_per_person;
								firstPass = false;
							}
						}

						//Normal default case
						else
						{
							costThisCycle = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res]
							* interv[i].who[status][sex][orient][act][alc][idu][vl] 
							* interv[i].cost_per_person
								* DT;
						}
					}//end of if "who"


					//Add the cost this cycle to the running total
					interv[i].total_prepurchased_cost += costThisCycle;
					interv[i].discounted_total_prepurchased_cost += get_discounted_value(costThisCycle);
			
				}//end of if "adulthood"
			}//end of if "basket[i]==true"
			

			/******************** Now calculate the utilization cost ******************************/
			interv[i].total_utilization_cost = interv[i].total_prepurchased_cost * (1 - interv[i].overall_effect_size);
			interv[i].discounted_total_utilization_cost = get_discounted_value(interv[i].total_utilization_cost);
		}

		//Now calculate the cost of care and treatment using the data from the progression model
		//Note that include the cost of care and treatment for adults as well as children

		//This used to say "if (status >=CARE )" but we changed CARE to INF to include
		//hospital costs for people not yet in care
		if (status>=INF && Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] > 0)
		{
			if (EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS) {
				costOfCareAndTreatThisCycle += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] *
					preCalcAdhStratifiedCostCareTreat[status][cd4][vl][res][sex][orient][act][alc][idu] *
					DT;
			} 
			else
				costOfCareAndTreatThisCycle += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] *
					getAdhStratifiedCostCareTreat(status, cd4, vl, res, sex, orient, act, alc, idu) *
					DT;
		}

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	costOfCareAndTreat += costOfCareAndTreatThisCycle;
	discounted_costOfCareAndTreat += get_discounted_value(costOfCareAndTreatThisCycle);
}

void save_population()
{ 
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	ofstream my_pop_file;

	cout<<"Writing population file"<<endl<<endl;

	my_pop_file.open ("Population.txt");
	my_pop_file.precision(25); //Set high precision because this data will be used to populate the array for future runs.

	my_pop_file <<"status\tsex\torient\tact\tage\talc\tidu\tvl\tcd4\tres\tpop"<< endl; // Header

	// Loop through entire population array and output its indices and values to file, one line per value.
	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS
		my_pop_file << status << "\t" << sex << "\t" << orient << "\t" << act << "\t" << age << "\t" << alc<< "\t" << idu << "\t" << vl << "\t" << cd4 << "\t" << res << "\t" << Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] << endl;
	END_LOOP_THROUGH_ALL_COMPARTMENTS

	my_pop_file.close(); 
}

void movehalf(int sex)
{
	int status, orient, act, age, alc, idu, vl, cd4, res;
	int oppSex = (sex + 1) % 2;
	double men = 0, women = 0;
	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS
		Comp[status][oppSex][orient][act][age][alc][idu][vl][cd4][res] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] * 0.5;
		Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] *= 0.5;
		if (sex==M) men += Comp[status][M][orient][act][age][alc][idu][vl][cd4][res];
		else if (sex==F) women += Comp[status][F][orient][act][age][alc][idu][vl][cd4][res];
	END_LOOP_THROUGH_ALL_COMPARTMENTS
	cout <<"Proportion male: " << men/(men + women) << endl;
}

void change_prevalence(bool high) 
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	double pop=0, infpop = 0, numbertomove;
	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS
		pop += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	if (status) infpop += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	END_LOOP_THROUGH_ALL_COMPARTMENTS
		cout << "Starting prevalence: "<< (infpop/pop) << endl;
	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS
		if (status != SUSC) 
		{
			numbertomove = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] * .9;
			if (high) 
			{ // Increase prevelance by 90%
				Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] += numbertomove;
				Comp[SUSC][sex][orient][act][age][alc][idu][0][0][0] -= numbertomove;
			} 
			else 
			{ // Reduce prevelance by 90%
				Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] -= numbertomove;
				Comp[SUSC][sex][orient][act][age][alc][idu][0][0][0] += numbertomove;
			}
		}
		END_LOOP_THROUGH_ALL_COMPARTMENTS
			pop = 0;
		infpop = 0;
		LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS
			pop += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		if (status != SUSC)	infpop += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		END_LOOP_THROUGH_ALL_COMPARTMENTS
			cout << "NEW PREVALENCE = " << (infpop/pop) << endl;
}

void read_population() 
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	double num;
	char buff[256];
	FILE *my_file;
	string pop_file = RESOURCE_FILE_PATH;
	pop_file += POPULATION_FILE;

	if((my_file = fopen(pop_file.c_str(), "r" )) == NULL )
	{
		printf("ERROR opening %s\n", pop_file.c_str() );
		exit(0);
	}

	fgets(buff, 100, my_file); // Skip the first line beacause it's just the header.
	
	while (fscanf(my_file,"%d%d%d%d%d%d%d%d%d%d%le", &status, &sex, &orient, &act, &age, &alc, &idu, &vl, &cd4, &res, &num) == 11) 
	{ 
		// Read indices and value
		Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] = num; // Assign value to population array.	
	}
	fclose(my_file);
}

//During the calibration period from 1997-2014, we increase the probability of being tested linearly over that time period.
//We also linearly increase the probability of being linked to care from 2003 to 2014.  (Before 2003, we assume no linkage to care.)
//We also linearly increase the proportion of women taking PMTCT meds from 1997 to 2011.
void rampUpValuesDuringCalibration()
{
	double x0, y0, x1, y1, x;

	/**************** First determine the probability of not being tested for HIV **************/

	x0 = 1997;													//the "time" at the start of calibration
	y0 = prob_not_being_tested_for_hiv_start_of_calibration;	//the prob of not being tested at the start of calibration
	x1 = 2014;													//the "time" at the end of caibration
	y1 = prob_not_being_tested_for_hiv_end_of_calibration;		//the prob of not being tested at the end of calibration
	x = max(t,0.0)+1997;										//how much time has passed since the start of the run (note that at first pass, t seems to be -0.02.  I don't remember why we did this.  But don't let t go less than 0.

	//baselineRR_raw[1][0] is the probability of not being tested for HIV
	baselineRR_raw[1][0] = linearInterpolate(x0, y0, x1, y1, x);

	/**************** Next determine the probability of failing to be linked to care **************/

	if (t < NUM_YEARS_NO_LINKAGE_TO_CARE_FROM_1997)	baselineRR_raw[2][0] = 1;		//from 1997 - 2003, we assume no linkage to care
	else		
	{
		x0 = 1997 + NUM_YEARS_NO_LINKAGE_TO_CARE_FROM_1997;	 //the "time" at the start of the ramp up period for people transitioning to care.  Previous to this time, nobody transferred to Care
		y0 = 1;									//the prob of failing to be diagnosed at the start of the ramp up period
		x1 = 2014;								//the "time" at the end of caibration ramp up period
		y1 = prob_failing_to_be_linked_to_care_once_diagnosed_end_of_calibration;		//the prob of failing to be diagnosed at the end of the ramp up period
		x = t+1997;								//how much time has passed since the start of the ramp up period

		//baselineRR_raw[2][0] is the probability of failing to be linked to care
		baselineRR_raw[2][0] = linearInterpolate(x0, y0, x1, y1, x);
	}

	/********** Ramp up the proportion of pregnant women on PMTCT meds **************/

	x0 = 1997;													//the "time" at the start of calibration
	y0 = PROP_PREGNANT_WOMAN_ON_PMTCT_1997;						//the prop using PMTCT at the start of calibration
	x1 = 2014;	//the "time" at the end of caibration
	
	//the prop using PMTCT at the end of calibration
	y1 = PROP_PREGNANT_WOMAN_ON_PMTCT_2014;			

	x = max(t,0.0)+1997;										//how much time has passed since the start of the run (note that at first pass, t seems to be -0.02.  I don't remember why we did this.  But don't let t go less than 0.

	//set the proportion of pregnant women on PMTCT meds
	pmtct_meds = linearInterpolate(x0, y0, x1, y1, x);
}



void printNumInEachStatus()
{
	double total_in_status[NUM_STATUS] = { 0 };
	double total_in_activity[NUM_ACT] = { 0 };
	double total_alc_adults[NUM_ALC] = { 0 };
	double total_alc_children[NUM_ALC] = { 0 };
	double total_idu[NUM_IDU] = { 0 };
	double numVLcat[NUM_VL] = { 0 };
	double numCD4cat[NUM_CD4] = { 0 };
	double numRescat[NUM_RES] = { 0 };
	double vl_on_treat[NUM_VL] = { 0 };
	double vl_not_on_treat[NUM_VL] = { 0 };
	double total_sex_alc[NUM_SEX][NUM_ALC] = { 0 };
	double total_infected_alc = { 0 };
	double total_in_sex[NUM_SEX] = { 0 };
	double total_inf_act[NUM_ACT]={0}; 
	double total_inf_sex_age[NUM_SEX][NUM_AGE]={0}; 
	double total_sex_age[NUM_SEX][NUM_AGE]={0}; 
	double total=0, infected=0, prop_inf=0;
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	static std::streambuf *sb = std::cout.rdbuf();

	//initialize totals to zero
	memcpy(total_in_status, zerod_array, sizeof(total_in_status));
	memcpy(total_in_activity, zerod_array, sizeof(total_in_activity));
	memcpy(total_alc_children, zerod_array, sizeof(total_alc_children));
	memcpy(total_alc_adults, zerod_array, sizeof(total_alc_adults));
	memcpy(total_idu, zerod_array, sizeof(total_idu));
	memcpy(numVLcat, zerod_array, sizeof(numVLcat));
	memcpy(numCD4cat, zerod_array, sizeof(numCD4cat));
	memcpy(numRescat, zerod_array, sizeof(numRescat));
	memcpy(vl_on_treat, zerod_array, sizeof(vl_on_treat));
	memcpy(vl_not_on_treat, zerod_array, sizeof(vl_not_on_treat));
	memcpy(total_sex_alc, zerod_array, sizeof(total_sex_alc));
	memcpy(total_in_sex, zerod_array, sizeof(total_in_sex));
	total_infected_alc = 0;

	LOOP_THROUGH_ALL_COMPARTMENTS

	total_in_status[status] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	total_in_activity[act] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	total_in_sex[sex] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	total_sex_alc[sex][alc] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	total_idu[idu] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	total_sex_age[sex][age]+= Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	if (age >= ADULTHOOD) total_alc_adults[alc] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
	else total_alc_children[alc] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];

		if (status>=INF && status<=TREAT)
		{
			numVLcat[vl] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			numCD4cat[cd4] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			numRescat[res] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			if (status == TREAT) vl_on_treat[vl] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			else vl_not_on_treat[vl] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			if (alc == ALC || alc == ALC_TREATED_NOT_CURED) total_infected_alc += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			
			if (age >= ADULTHOOD) 
			{
				total_inf_act[act]+= Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
				total_inf_sex_age[sex][age]+= Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
				infected+= Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			}
		}

		if (age >= ADULTHOOD) total += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];

	END_LOOP_THROUGH_ALL_COMPARTMENTS

	//Calculate some ratios
	if (total_in_status[DET] + total_in_status[CARE] + total_in_status[TREAT] == 0) prop_detected_on_treatment = -999; 
	else prop_detected_on_treatment = total_in_status[TREAT]/(total_in_status[DET]+total_in_status[CARE]+total_in_status[TREAT]);
	prop_infected_on_treatment = total_in_status[TREAT]/(total_in_status[HVL]+total_in_status[INF]+total_in_status[DET]+total_in_status[CARE]+total_in_status[TREAT]);
	prop_infected_detected = (total_in_status[DET]+total_in_status[CARE]+total_in_status[TREAT])/(total_in_status[HVL]+total_in_status[INF]+total_in_status[DET]+total_in_status[CARE]+total_in_status[TREAT]);
	
	//prints the output first to an output file, then to the console
	for (int output_style=0; output_style<2; output_style++)
	{
		if (output_style==0) std::cout.rdbuf(startAndEnd.rdbuf());
		else std::cout.rdbuf(sb);

		cout<<"\n\nPop at time:  "<<t<<endl;

		for (int status=0; status<NUM_STATUS; status++)
		{
			cout<<"status "<<status<<": \t"<<total_in_status[status]<<endl;
		}
		cout<<endl;

		for (int sex=0; sex<NUM_SEX; sex++)
		{
			cout<<"sex "<<sex<<": \t"<<total_in_sex[sex]<<endl;
		}
		cout<<endl;

		for (int act=0; act<NUM_ACT; act++)
		{
			cout<<"act "<<act<<": \t"<<total_in_activity[act]<<endl;
		}
		cout<<endl;

		for (int alc=0; alc<NUM_ALC; alc++)
		{
			cout<<"alc adults "<<alc<<": \t"<<total_alc_adults[alc]<<endl;
		}
		cout<<endl;

		for (int alc=0; alc<NUM_ALC; alc++)
		{
			cout<<"alc children "<<alc<<": \t"<<total_alc_children[alc]<<endl;
		}
		cout<<endl;

		for (int idu = 0; idu<NUM_IDU; idu++)
		{
			cout << "idu " << idu << ": \t" << total_idu[idu] << endl;
		}
		cout << endl;

		for (int act = 0; act<NUM_ACT; act++)
		{
			cout << "infected in act " << act << ": \t" << total_inf_act[act] << endl;
		}
		cout << endl;
		
		/*for (int sex = 0; sex<NUM_SEX; sex++)
		{
			for (int age = ADULTHOOD; age<NUM_AGE; age++)
			{
				prop_inf=total_inf_sex_age[sex][age]/total_sex_age[sex][age]; 
				cout << "prop infected sex " << sex << " age " << age << ": \t" << prop_inf << endl;
			}
		}*/

		cout << "prop infected total " << infected/total << endl; 
		cout << "Male % alc: \t" << (100*total_sex_alc[M][1]/(total_sex_alc[M][0] + total_sex_alc[M][1])) << endl;
		cout << "Female % alc: \t" << (100*total_sex_alc[F][1]/(total_sex_alc[F][0] + total_sex_alc[F][1])) << endl;
		
		cout<<endl;

		cout<<"alc adults infected:\t"<<total_infected_alc<<endl;
		cout<<"alc adults susceptible:\t"<<total_alc_adults[ALC] + total_alc_adults[ALC_TREATED_NOT_CURED] - total_infected_alc<<endl;

		cout<<endl;

		for (int vl=0; vl<NUM_VL; vl++)
		{
			cout<<"vl "<<vl<<": \t"<<numVLcat[vl]<<endl;
		}

		cout<<endl;

		for (int vl=0; vl<NUM_VL; vl++)
		{
			cout<<"vl (on treatment) "<<vl<<": \t"<<vl_on_treat[vl]<<endl;
		}

		cout<<endl;

		for (int vl=0; vl<NUM_VL; vl++)
		{
			cout<<"vl (not on treatment) "<<vl<<": \t"<<vl_not_on_treat[vl]<<endl;
		}

		cout<<endl;

		for (int cd4=0; cd4<NUM_CD4; cd4++)
		{
			cout<<"cd4 "<<cd4<<": \t"<<numCD4cat[cd4]<<endl;
		}
		cout<<endl;

		for (int res=0; res<NUM_RES; res++)
		{
			cout<<"res "<<res<<": \t"<<numRescat[res]<<endl;
		}
		cout<<endl;

		if (prop_detected_on_treatment == -999) cout<< "Prop detected on treatment: \t"<< "None in DET, CARE, or TREAT"<<endl;
		else cout<< "Prop detected on treatment: \t"<< prop_detected_on_treatment <<endl;
		cout<< "Prop infected on treatment: \t"<< prop_infected_on_treatment<<endl;
		cout<< "Prop infecteds detected: \t"<< prop_infected_detected<<endl;
		cout<<endl;

		cout<<"Total people in model: "<<total<<endl;
	}
}

//Note that the alcohol intervention now has it's own function (even though it moves people) as it no longer fits the pattern.
void processInterventionsThatMovePeople()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	int i, p;

	double final_effect_size[NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE][NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL], effect_size;

	//Start with the first pathway that moves people (which are located at the end of the list of pathways) and process until the end of the list of pathways
	for (p=NUM_PATHWAYS; p<NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE; p++)
	{	
		LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS 

			//intialize final effect size to no effect
			final_effect_size[p][status][sex][orient][act][alc][idu][vl] = 1;

			//Find the cumulative effect over all interventions that activate that pathway
			for (i=0; i < NUM_INTERVENTIONS; i++)
			{
				//If an intervention is turned on that activates this pathway and the effect sizes is not 1 (implying there actually is an effect)
				if (basket[i] == true && !(isEqual(interv[i].effect_size[p], 1)))
				{
					/*NOTE:  Ideally, we would be able to apply these interventions according to the "who", but it's a little 
					complicated (I think) because propInAct is only over sex and activity.  So more code changes would have to
					take place than just implementing the "who" as below.  Also, note that this line of code doesn't currently
					work because the variables status, sex, act and alc never get set in this function.	
					effect_size = linearInterpolate(0, 1, 1, interv[i].effect_size[p], interv[i].who[status][sex][act][alc]);
					To really implement this, you would end up with an effect size array over the same elements as the "who" array"*/

					//To find the actual effect size for an intervention, a who of '0' should have no effect (resulting in an effect size of 1)
					//A who of '1' applies the effect to everyone in the compartment, hence the full effect size applies
					//A who of some value in between gets partial effect, so interpolate!
					//linearInterpolate(x0, y0, x1, y1, x)  applied as....
					//linearInterpolate(who 0, effect size 1, who 1, full intervention effect size, who for this intervention

					//ANIK
					effect_size = linearInterpolate(0, 1, 1, interv[i].effect_size[p], interv[i].who[status][sex][orient][act][alc][idu][vl]);

					if (MODE == OPTIMISTIC) final_effect_size[p][status][sex][orient][act][alc][idu][vl] *= effect_size;			//Optimistic

					else final_effect_size[p][status][sex][orient][act][alc][idu][vl] = min(final_effect_size[p][status][sex][orient][act][alc][idu][vl], effect_size);		//Pessimistic
				}
			}

		END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS

	} //End of for loop over pathways to determine final_effect_size array



	//Loop through the pathways again and apply the final_effect_sizes found above
	for (p=NUM_PATHWAYS; p<NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE; p++)
	{	
		switch (p)
		{

		//Pathway 8: Fewer partners (Applies to group(s) with concurrent partners, but not CSWs)
		case (8): 

			//For both men and women, move people from the higher activity groups into the next lower activity group according to effect size

			//Note that at this time (9/29/2014) VL targeting only works for the alcohol intervention (it will NOT work for this pathway).
			//We made a change to allow intervention targeting to include VL status.  Hence, final_effect_size now has an index for VL.  
			//But propInAct does not vary by VL status.  
			//If in the future the effect sizes vary by VL we need to build this into the propInAct variable, which could be a mess

			//We have to set a value for VL, so we set it to 0.  Ok, because the effect sizes for VL should all be the same for now.
			vl = 0;

			for (status=0; status<NUM_ALIVE_STATE; status++)
			{
				for (sex=0; sex<NUM_SEX; sex++)
				{
					for (orient=0; orient<NUM_ORIENT; orient++)
					{
						for (alc=0; alc<NUM_ALC; alc++)
						{
							for (idu=0; idu<NUM_IDU; idu++)
							{
								propOrientActAlcIDU[status][sex][orient][1][alc][idu] = propOrientActAlcIDU[status][sex][orient][1][alc][idu] + propOrientActAlcIDU[status][sex][orient][2][alc][idu] * (1 - final_effect_size[p][status][sex][orient][2][alc][idu][vl]);
								propOrientActAlcIDU[status][sex][orient][2][alc][idu] = propOrientActAlcIDU[status][sex][orient][2][alc][idu] * final_effect_size[p][status][sex][orient][2][alc][idu][vl];
							}
						}
					}
				}
			}

			LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS

				if (act==1) 
				{
					Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] + Comp[status][sex][orient][2][age][alc][idu][vl][cd4][res] * (1 - final_effect_size[p][status][sex][orient][2][alc][idu][vl]);
				}
				else if (act==2) 
				{
					Comp[status][sex][orient][2][age][alc][idu][vl][cd4][res] = Comp[status][sex][orient][2][age][alc][idu][vl][cd4][res] * final_effect_size[p][status][sex][orient][2][alc][idu][vl];
				}

				END_LOOP_THROUGH_ALL_COMPARTMENTS

					break;

		//Pathway 9: Reduction in transactional sex or intergenerational partnering
		case (9): 

			break;

		//Pathway 13:  Fewer partners (just high activity females) 
		//For highest activity women, move people to the next lower activity group according to effect size
		case (13):

			//Note that at this time (9/29/2014) VL targeting only works for the alcohol intervention (it will NOT work for this pathway).
			//We made a change to allow intervention targeting to include VL status.  Hence, final_effect_size now has an index for VL.  
			//But propInAct does not vary by VL status.  
			//If in the future the effect sizes vary by VL we need to build this into the propInAct variable, which could be a mess

			//We have to set a value for VL, so we set it to 0.  Ok, because the effect sizes for VL should all be the same for now.
			vl = 0;

			for (status=0; status<NUM_ALIVE_STATE; status++)
			{
				for (orient=0; orient<NUM_ORIENT; orient++)
				{
					for (alc=0; alc<NUM_ALC; alc++)
					{
						for (idu=0; idu<NUM_IDU; idu++)
						{
							propOrientActAlcIDU[status][F][orient][2][alc][idu] = propOrientActAlcIDU[status][F][orient][2][alc][idu] + propOrientActAlcIDU[status][F][orient][3][alc][idu] * (1 - final_effect_size[p][status][F][orient][3][alc][idu][vl]);
							propOrientActAlcIDU[status][F][orient][3][alc][idu] = propOrientActAlcIDU[status][F][orient][3][alc][idu] * final_effect_size[p][status][F][orient][3][alc][idu][vl];
						}	
					}
				}
			}

			LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS

				if (sex==F)
				{
					if (act==2)
					{
						Comp[status][sex][orient][2][age][alc][idu][vl][cd4][res] += Comp[status][sex][orient][3][age][alc][idu][vl][cd4][res] * (1 - final_effect_size[p][status][sex][orient][3][alc][idu][vl]);
					}
					else if (act==3)
					{
						Comp[status][sex][orient][3][age][alc][idu][vl][cd4][res] = Comp[status][sex][orient][3][age][alc][idu][vl][cd4][res] * final_effect_size[p][status][sex][orient][3][alc][idu][vl];
					}
				}

				END_LOOP_THROUGH_ALL_COMPARTMENTS

					break;

		} //end of switch

	}//End of for loop
}


double getAlphaMult(int sex, int orient, int act, int alc, int idu, int vl, int p_status, int p_sex, int p_orient, int p_act, int p_alc, int p_idu, int p_vl)
{
	double  baseline_Prob_susc, baseline_Prob_infector, who_susc, who_infector, min_this_path, thisEffectSize, alphaMult;

	double final_alphaMult = 1;		//initialize to no effect
	double alphaMult_susc = 1;		//initialize to no effect
	double alphaMult_infector = 1;	//initialize to no effect

	for (int p=0; p<NUM_PATHWAYS; p++)
	{
		//Do only for alpha multiplier pathways
		if ((p==0 || p>=4))		
		{
			baseline_Prob_susc = baselineProb[p][SUSC][sex][orient][act][alc][idu];
			baseline_Prob_infector = baselineProb[p][p_status][p_sex][p_orient][p_act][p_alc][p_idu];
			
			/*********************************** Special Case*************************************/
			//Because FSW are more likey to use condoms with their clients compared to their regular partner, we had to add a multiplier 
			//can modify the baselineRR values in the getAlphaMult function.  The current value in the matrix is for a FSW partnering with her 
			//regular partner (ACT1), so use this when they are partnering with a client (ACT2 and ACT3) 
			
			if (p==0) //only for condom use pathway
			{
				#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM == 59 || RUN_NUM == 60)) || PROBABILISTIC_RUN
					#if PROBABILISTIC_RUN
						if (sex == F && act == 3 && p_act > 1 && p_sex == M) baseline_Prob_susc *= prob_pull[29]; 
						else if (sex == M && act > 1 && p_act == 3 && p_sex == F) baseline_Prob_infector *= prob_pull[29];
					#else
						if (sex == F && act == 3 && p_act > 1 && p_sex == M) baseline_Prob_susc *= sensitivity_values[RUN_NUM];
						else if (sex == M && act > 1 && p_act == 3 && p_sex == F) baseline_Prob_infector *= sensitivity_values[RUN_NUM]; 
					#endif 
				#else
					if (sex == F && act == 3 && p_act > 1 && p_sex == M) baseline_Prob_susc *= FSW_CONDOM_MULTIPLIER; 
					else if (sex == M && act > 1 && p_act == 3 && p_sex == F) baseline_Prob_infector *= FSW_CONDOM_MULTIPLIER; 
				#endif
			}
			/************************************ End Special Case *********************************/

			min_this_path = 1.0;			//initialize to no effect for this pathway

			for (int i=0; i<NUM_INTERVENTIONS; i++)
			{
				if (basket[i]==true) thisEffectSize = interv[i].effect_size[p]; 
				else thisEffectSize = 1.0;

				//ANIK
				who_susc = interv[i].who[SUSC][sex][orient][act][alc][idu][vl];
				who_infector = interv[i].who[p_status][p_sex][p_orient][p_act][p_alc][p_idu][p_vl]; 

				alphaMult_susc = (baseline_Prob_susc - (baseline_Prob_susc*(1-thisEffectSize)*who_susc))*1 
					+ (1-baseline_Prob_susc + baseline_Prob_susc*(1-thisEffectSize)*who_susc)*alpha_sex_mult[p];

				alphaMult_infector = (baseline_Prob_infector - (baseline_Prob_infector*(1-thisEffectSize)*who_infector))*1 
					+ (1-baseline_Prob_infector + baseline_Prob_infector*(1-thisEffectSize)*who_infector)*alpha_sex_mult[p];

				switch (p)
				{
					//Condoms:  Take mean of resultant alpha mults
				case (0):	
					alphaMult = (alphaMult_susc + alphaMult_infector) / 2.0;
					if (baseline_Prob_susc == NOT_APPLICABLE || baseline_Prob_infector == NOT_APPLICABLE) {cout<<"ERROR 0!!!"<<endl; exit(0);}
					break;	

					//PreP:  Use susceptible's alpha mult
				case (4): 
					alphaMult = alphaMult_susc;	
					if (baseline_Prob_susc == NOT_APPLICABLE) {cout<<"ERROR 4!!!"<<endl; exit(0);}
					break;

					//Untreated STI:  Use max of resultant alpha mults
				case (5): 
					alphaMult = max(alphaMult_susc, alphaMult_infector);	
					if (baseline_Prob_susc == NOT_APPLICABLE || baseline_Prob_infector == NOT_APPLICABLE) {cout<<"ERROR 5!!!"<<endl; exit(0);}
					break;	

					//Microbicides:  If susceptible is a woman, use her alpha mult, else use 1.0.
				case (6): 
					alphaMult = sex==F?alphaMult_susc:1.0;	
					if (sex==F && baseline_Prob_susc == NOT_APPLICABLE) {cout<<"ERROR 6!!!"<<endl; exit(0);}
					break;	

					//Circumcision:  If susceptible is a man, use his alpha mult, else use 1.0.
				case (7): 
					alphaMult = sex==M?alphaMult_susc:1.0;	
					if (sex==M && baseline_Prob_susc == NOT_APPLICABLE) {cout<<"ERROR 7!!!"<<endl; exit(0);}
					break;		

				default: 
					cout<<"Pathway error.  Exiting."<<endl;  
					exit(0);
					break;
				}

				if (alphaMult < 0)
				{
				cout<<"*****************  Uh oh!  alphaMult is less than zero! - it's "<<alphaMult<<endl;
				exit (0);
				}
			
				if (MODE == OPTIMISTIC) final_alphaMult *= alphaMult;

				else min_this_path = min(min_this_path, alphaMult);
			}

			if (MODE == PESSIMISTIC) final_alphaMult *= min_this_path;
		}
	}
	return final_alphaMult;
}

#if ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION

//We call this function after the baseline probabilities have already been calculated with all alcohol RR's set to 1.  In this function, we now multiple the baseline prob by the RR
void postProcessBaselineProbsForAlcAnalysis()
{
	int path, status, sex, act, alc;

	LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS 

		//if (status==TREAT) cout<<"before: run "<<run<<" path "<<path<<" status "<<status<<" sex "<<sex<<" act "<<act<<" alc "<<alc<<"   "<<baselineProb[path][status][sex][act][alc]<<endl;

		switch (run)
	{

		case(2):	
		case(3):
			path = 0;
			baselineRR_raw[path][5] = baselineRR_raw_backup_condoms; 

			if (baselineProb[path][status][sex][act][alc] != NOT_APPLICABLE && alc==1) baselineProb[path][status][sex][act][alc] *= baselineRR_raw[path][5];
			baselineProb[path][status][sex][act][alc] = min(baselineProb[path][status][sex][act][alc], 1.0);		//Can't have a probability greater than 1

			path = 3;
			baselineRR_raw[path][5] = baselineRR_raw_backup_nonadherence; 

			if (baselineProb[path][status][sex][act][alc] != NOT_APPLICABLE && alc==1) baselineProb[path][status][sex][act][alc] *= baselineRR_raw[path][5];
			baselineProb[path][status][sex][act][alc] = min(baselineProb[path][status][sex][act][alc], 1.0);		//Can't have a probability greater than 1

			path = 5;
			baselineRR_raw[path][5] = baselineRR_raw_backup_STIs; 

			if (baselineProb[path][status][sex][act][alc] != NOT_APPLICABLE && alc==1) baselineProb[path][status][sex][act][alc] *= baselineRR_raw[path][5];
			baselineProb[path][status][sex][act][alc] = min(baselineProb[path][status][sex][act][alc], 1.0);		//Can't have a probability greater than 1

			break;

		case(4):	
		case(5):
			path = 0;
			baselineRR_raw[path][5] = baselineRR_raw_backup_condoms; 
			break;
		case(6):	
		case(7):
			path = 3;
			baselineRR_raw[path][5] = baselineRR_raw_backup_nonadherence; 
			break;
		case(8):	
		case(9):  
			path = 5;
			baselineRR_raw[path][5] = baselineRR_raw_backup_STIs; 
			break;
	}


	if (run>=4 && baselineProb[path][status][sex][act][alc] != NOT_APPLICABLE && alc==1) baselineProb[path][status][sex][act][alc] *= baselineRR_raw[path][5];
	baselineProb[path][status][sex][act][alc] = min(baselineProb[path][status][sex][act][alc], 1.0);		//Can't have a probability greater than 1

	END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS 
}
#endif


void setPropOrientActAlcIDU()
{
	int status, sex, orient, act, idu, alc;

	//proportion of status, sex, orient, age, alc and idu in each activity group
	//propOrientActAlcIDU[status][sex][orient][alc][idu][act]
	//note that propOrientActAlcIDU was originally just propOrientActAlcIDU[sex][act] but we expanded it
	//because interventions could be applied differently according to status, orient, idu, etc. according to the "who"
	for (status=0; status<NUM_ALIVE_STATE; status++)
	{
		for (sex=0; sex<NUM_SEX; sex++)
		{
			for (orient=0; orient<NUM_ORIENT; orient++)
			{
				for (act=0; act<NUM_ACT; act++)
				{
					for (alc=0; alc<NUM_ALC; alc++)
					{
						for (idu=0; idu<NUM_IDU; idu++)
						{
							propOrientActAlcIDU[status][sex][orient][act][alc][idu] = propSexOrientInAct[sex][orient][act] * propOrientInPop[sex][orient] * propIDUinAct[act][idu] * prop_alc[status][sex][act][alc];
						}
					}
				}
			}
		}
	}
}



#if OPTIMIZING_ALCOHOL_INTERVENTION || ALCOHOL_SURFACE || (OPTIMIZING_PREVENTION_PORTFOLIO && (SET_RR_ALC > 0))
void postProcessBaselineProbsForOptimizingAlcoholIntervention()
{
	int path, status, sex, act, alc;
	double mult;

	int path_array[3] = {0,3,5};

	for (int path_index=0; path_index <3; path_index++)
	{

		LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS

			//set the path to condoms, stis, or adherence
			path = path_array[path_index];

		if (baselineProb[path][status][sex][act][alc] != NOT_APPLICABLE && alc==1) 
		{
			switch (path)
			{
			case 0:  mult = condomRR / baselineRR_raw[path][5];  break;
			case 3:  mult = adhereRR / baselineRR_raw[path][5];  break;
			case 5:  mult = stiRR / baselineRR_raw[path][5];  break;
			default:  {cout << "Error in postProcessBaselineProbsForOptimizingAlcoholIntervention().  Exiting\n";  exit(0);}
			}

			//Make sure new baseline prob isn't greater than 1!
			baselineProb[path][status][sex][act][alc] = min(baselineProb[path][status][sex][act][alc] * mult, 1.0);		//Can't have a probability greater than 1
		}

		END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS
	}
}
#endif

void setAlcoholicProportionsAtStartOfRun()
{
	int status, sex, orient, act, age, idu, vl, cd4, res;

	double total;

	for (status=0; status<NUM_ALIVE_STATE; status++)
	{
		for (sex=0; sex<NUM_SEX; sex++)
		{
			for (orient=0; orient<NUM_ORIENT; orient++)
			{
				for (act=0; act<NUM_ACT; act++)
				{
					for (age=0; age<NUM_AGE; age++)
					{
						for (idu = 0; idu<NUM_IDU; idu++)
						{
							for (vl=0; vl<NUM_VL; vl++)
							{
								for (cd4=0; cd4<NUM_CD4; cd4++)
								{
									for (res = 0; res < NUM_RES; res++)
									{
										//determine the proportion alcoholic for each compartment at the start of the model

										//prevent divide by zero errors in the denominator
										total = Comp[status][sex][orient][act][age][ALC][idu][vl][cd4][res] + Comp[status][sex][orient][act][age][NONALC][idu][vl][cd4][res];

										//reallocate to desired settings
										Comp[status][sex][orient][act][age][ALC][idu][vl][cd4][res] = prop_alcohol[sex] * total;
										Comp[status][sex][orient][act][age][NONALC][idu][vl][cd4][res] = (1 - prop_alcohol[sex]) * total;
									}}}}}}}}}
}

void choose_rate_files()
{
	// readHIVprogrates() adds file numbers and file extensions to read adherence files.
	string root_of_oncare = RESOURCE_FILE_PATH; 
	string path_to_off_care = RESOURCE_FILE_PATH;

	switch (cd4_treat_threshold) 
	{
	case 200:	root_of_oncare += CD4_200_ONCARE; break;
	case 350:	root_of_oncare += CD4_350_ONCARE; break;
	case 500:	root_of_oncare += CD4_500_ONCARE; break;
	case 10000:	root_of_oncare += CD4_10K_ONCARE; break;

	default:
		{
			cout << "Sorry, CD4 treatment threshold ("<< cd4_treat_threshold << ") is not a valid value." << endl;
			exit(0);
		}
		break;
	}

	path_to_off_care += OFF_CARE_RATES;

	for (int intervention = 0; intervention < NUM_INTERVENTIONS; intervention++)
	{
		if (basket[intervention] && interv[intervention].effect_size[CD4_THRESH_PATHWAY] == CD4_TREAT_THRESH_10000)
		{
			root_of_oncare = RESOURCE_FILE_PATH;
			root_of_oncare += CD4_10K_ONCARE;
			break; // If any intervention activates this pathway we can break out of the loop.
		} 
		else if (basket[intervention] && interv[intervention].effect_size[CD4_THRESH_PATHWAY] != CD4_TREAT_THRESH_DEFAULT)
		{
			cout << "The effect size for pathway "<< CD4_THRESH_PATHWAY<< " must be 0 or 1."<< endl;
			exit(0);
		}
	}

	//Now call readHIVprogrates to read in the proper rate files.
	readHIVprogrates(path_to_off_care, root_of_oncare);
}

double getProportionOfWomenWhoGiveBirthInAYear()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;

	double totalWomenInAgeGroup[NUM_AGE], totalWomen=0; 

	double fertilityRate;
	double retval = 0;

	//initialize the number of women in each age strata to 0
	memcpy(totalWomenInAgeGroup, zerod_array, sizeof(totalWomenInAgeGroup));

	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS 

		if (sex==F && age >= ADULTHOOD) 
		{
			totalWomenInAgeGroup[age] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
			totalWomen += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res];
		}

	END_LOOP_THROUGH_ALL_COMPARTMENTS 


	for (int age = ADULTHOOD; age < NUM_AGE; age++) 
	{
		//We put in an array of fertility rates from 1997 to year 2014.  When we are calibrating, we should
		//use the rate from the appropriate year.  When we are running for real, we should use the 2014 value.
		fertilityRate = fertRate[NUM_YEARS_CALIBRATION-1][age];

		//Add up the proportion of women in each group times their fertility rate to get the total proportion of women
		//who give birth in a year.
		retval += totalWomenInAgeGroup[age] * fertilityRate / totalWomen;

	}

	return retval;
}

void setPropWomenOnPMTCT()
{
	for (int i=0; i < NUM_INTERVENTIONS; i++)
	{
		//Pathway 11 is the PMTCT pathway
		if (basket[i] && interv[i].effect_size[11] > pmtct_meds)
		{
			pmtct_meds = interv[i].effect_size[11];
		}
	}
}

#if ALCOHOL_SURFACE
void initializeVarsForAlcoholSurfaceAnalysis()
{
	// Turn all interventions off.
	for (int i = 0; i < NUM_INTERVENTIONS; i++)
	{
		basket[i] = 0;
	}

	// Turn on intervention 3: alcohol issues.
	basket[3] = 1;

	//Set the cost of the alcohol intervention
	interv_costs_per_person[3] = alcohol_cost[cost_index];

	//Set the effect size of the alcohol intervention on the pathway that moves people out of alcohol compartments
	interv_effect_size[3][10] = alcohol_effect_size[effect_index];

	//We are varying the alcohol RRs for the three pathways that are affected by alcohol use.
	condomRR = alcohol_rr_modifier[0][rr_index] * baselineRR_raw[0][5];
	adhereRR = alcohol_rr_modifier[1][rr_index] * baselineRR_raw[3][5];
	stiRR = alcohol_rr_modifier[2][rr_index] * baselineRR_raw[5][5];
}
#endif

void LifeYears_and_QALYS()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	// At each time slice we add to a running sum the number of life years, discounted life years, QALYs and disc QALYs.
	//We will maintain LYs, QALYs, and disc LYs and disc QALYs for 0)  All pats, 1)  All HIV+  2)  All On Care + Treat

	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS

		for (int group=0; group<3; group++)
		{
			if (group == ALL_PATS ||
				(group == ALL_HIV_POS && status >= HVL) ||
				(group == ALL_CARE_AND_TREAT && status >= CARE))
			{
				total_life_years[group] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] * DT;
				total_disc_life_years[group] += get_discounted_value(Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] * DT);

				total_qa_life_years[group] += Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] * DT * Utility(status, cd4);
				total_disc_qa_life_years[group] += get_discounted_value(Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] * DT * Utility(status, cd4));
			}
		}

	END_LOOP_THROUGH_ALL_COMPARTMENTS
}

/*******************************************************************************************************
 Returns a value in [0,1] inidicating a utility weight that is based on various patient parameters.
 Based on values from a Freedberg paper.
 ********************************************************************************************************/
double Utility(int status, int cd4)
{
	double retval = 0;
	
	if (status == SUSC)
	{
		return 1;
	}
	
	if (cd4<1)
	{
		#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM==71 || RUN_NUM==72)) || PROBABILISTIC_RUN
			#if PROBABILISTIC_RUN
				retval += prob_pull[37]; 
			#else
				retval += sensitivity_values[RUN_NUM];
			#endif 
		#else 
			retval += UTILITY_CD4_LESS_THAN_50; 
		#endif 
	}
	else if (cd4 < 2)
	{
		#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM==73 || RUN_NUM==74)) || PROBABILISTIC_RUN
			#if PROBABILISTIC_RUN
				retval += prob_pull[38]; 
			#else
				retval += sensitivity_values[RUN_NUM];
			#endif 
		#else 
			retval += UTILITY_CD4_50_TO_200;
		#endif 
	}
	else 
	{
		#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM==75 || RUN_NUM==76)) || PROBABILISTIC_RUN
			#if PROBABILISTIC_RUN
				retval += prob_pull[39]; 
			#else 
				retval += sensitivity_values[RUN_NUM];
			#endif
		#else 
			retval += UTILITY_CD4_ABOVE_200;
		#endif 
	}

	if (status == TREAT)
	{
		#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM==77 || RUN_NUM==78)) || PROBABILISTIC_RUN
			#if PROBABILISTIC_RUN
				retval += prob_pull[40]; 
			#else
				retval += sensitivity_values[RUN_NUM];
			#endif 
		#else 
			retval += DELTA_UTILITY_WITH_HAART;
		#endif 
	}

	
	return retval;
}

double get_discounted_value(double val_to_discount)
{
	#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM == 79 || RUN_NUM == 80 ))  
		return val_to_discount/pow(1+sensitivity_values[RUN_NUM], t);
	#else 
		return val_to_discount/pow(1+DISC_RATE, t);
	#endif 
}


//*****************************EFFICIENCY*********************************
void precalc_getAlphaMult()
{
	// This function precalculates values of getAlphaMult() for all possible calls to that function. This keeps the function from being called more than once for each combination of inputs.
	// Values are stored in the array AlphaMult_pre which has indices matching the parameters for getAlphaMult().
	int sex, orient, act, alc, idu, vl, p_status, p_sex, p_orient, p_act, p_alc, p_idu, p_vl; 

	//initialize to zeros
	memcpy(preCalculatedAlphMultArray, zerod_array, sizeof(preCalculatedAlphMultArray));

	for (sex = 0; sex < NUM_SEX; sex++)
	{
		for (orient = 0; orient < NUM_ORIENT; orient++)
		{
			for (act = FIRST_ACTIVE_CLASS; act < NUM_ACT; act++)
			{
				for (alc = 0; alc < NUM_ALC; alc++)
				{
					for (idu = 0; idu < NUM_IDU; idu++)
					{
						for (vl = 0; vl < NUM_VL; vl++)
						{
							for (p_status = SUSC; p_status < NUM_ALIVE_STATE; p_status++)
							{
								for (p_sex = 0; p_sex < NUM_SEX; p_sex++)
								{
									for (p_orient = 0; p_orient < NUM_ORIENT; p_orient++)
									{
										for (p_act = FIRST_ACTIVE_CLASS; p_act < NUM_ACT; p_act++)
										{
											for (p_alc = 0; p_alc < NUM_ALC; p_alc++)
											{
												for (p_idu = 0; p_idu < NUM_IDU; p_idu++)
												{
													for (p_vl = 0; p_vl < NUM_VL; p_vl++)
													{
														preCalculatedAlphMultArray[sex][orient][act][alc][idu][vl][p_status][p_sex][p_orient][p_act][p_alc][p_idu][p_vl] 
														= getAlphaMult(sex, orient, act, alc, idu, vl, p_status, p_sex, p_orient, p_act, p_alc, p_idu, p_vl);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void preCalculateAdhStratHIVProgRateAndCostArrays()
{
	// This function precalculates values of getHIVProgRate() for all possible calls to that function. This keeps the function from being called more than once for each combination of inputs.
	// Values are stored in the array HIV_Prograt_pre which has indices matching the parameters for getHIVProgRate() with the exception that RATE and COST are represented by 0 and 1 as the first index of the array.
	int status, cd4, vl, res, cd4_d, vl_d, treat_d, res_d, dead_d, sex, orient, act, alc, idu;

	for (status = SUSC; status < NUM_ALIVE_STATE; status++)
	{
		for (sex = M; sex <= F; sex++)
		{
			for (orient = 0; orient < NUM_ORIENT; orient++)
			{
				for (act = 0; act < NUM_ACT; act++)
				{
					for (alc = 0; alc < NUM_ALC; alc++)
					{
						for (idu = 0; idu < NUM_IDU; idu++)
						{ 
							for (cd4 = 0; cd4 < NUM_CD4; cd4++)
							{
								for (vl = 0; vl < NUM_VL; vl++)
								{
									for (res = 0; res < NUM_RES; res++)
									{
										//The cost array is not dependent on destination compartments
										preCalcAdhStratifiedCostCareTreat[status][cd4][vl][res][sex][orient][act][alc][idu] = getAdhStratifiedCostCareTreat(status, cd4, vl, res, sex, orient, act, alc, idu);

										for (dead_d = 0; dead_d < 2; dead_d++)
										{
											for (cd4_d = 0; cd4_d < NUM_CD4; cd4_d++)
											{
												for (vl_d = 0; vl_d < NUM_VL; vl_d++)
												{
													for (treat_d = 0; treat_d < 2; treat_d++)
													{
														for (res_d = 0; res_d < NUM_RES; res_d++)
														{
															preCalcAdhStratifiedHIVProgRate[status][cd4][vl][res][cd4_d][vl_d][treat_d][res_d][dead_d][sex][orient][act][alc][idu] = getAdhStratifiedHIVProgRate(status, cd4, vl, res, cd4_d, vl_d, treat_d, res_d, dead_d, sex, orient, act, alc, idu);
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

// We need to be sure that every compartment with living people has at least one rate from that starting compartment. This function will report an error and quit if this is not the case.
void check_HIVprogrates()
{
	int curr_cd4, curr_vl, curr_onTreat, curr_resistance, adh;
	bool has_problems_so_quit = false;
	
	//All adherence levels includes the off care file, as well
	for (adh = 0; adh < NUM_ADHERENCE; adh++)
	{
		for (curr_cd4 = 0; curr_cd4 < NUM_CD4; curr_cd4++)
		{
			for (curr_vl = 0; curr_vl < NUM_VL; curr_vl++)
			{
				// When adh is 11 (OFFCARE_INDEX) we're using the off care rate file, and that has no one on treatment. 
				//Hence: (adh==11 ? 1 : NUM_TREATORNOT) will prevent us from looking for on treatment rates in the off care file
				for (curr_onTreat = 0 ; curr_onTreat < (adh==OFFCARE_INDEX ? 1 : NUM_TREATORNOT); curr_onTreat++)
				{
					for (curr_resistance = 0 ; curr_resistance < NUM_RES; curr_resistance++)
					{
						if (rate_is_zero(curr_cd4, curr_vl, curr_onTreat, curr_resistance, adh))
						{
							//We only want to print the header once.
							//has_problems_so_quit will still be false here after we've found the first problem
							//So use it as a trigger to print the header
							if (!has_problems_so_quit)
							{
								has_problems_so_quit = true;
								cout <<"\nThe following start settings have no rates\nCD4\tVL\tTreat\tResistance\tAdh\n";
							}
						
							cout << curr_cd4 << "\t" << curr_vl << "\t" << curr_onTreat << "\t" << curr_resistance << "\t" << adh << endl;
						}
					}
				}
			}
		}
	}
	if (has_problems_so_quit)
	{
		exit(0);
	}
}

bool rate_is_zero(int curr_cd4, int curr_vl, int curr_onTreat, int curr_resistance, int adh)
{
	int dest_cd4, dest_vl, dest_onTreat, dest_resistance, dead;

	for (dest_cd4 = 0; dest_cd4 < NUM_CD4; dest_cd4++)
	{
		for (dest_vl = 0; dest_vl < NUM_VL; dest_vl++)
		{
			for (dest_onTreat = 0; dest_onTreat < NUM_TREATORNOT; dest_onTreat++)
			{
				for (dest_resistance = 0 ; dest_resistance < NUM_RES; dest_resistance++)
				{
					for (dead = 0; dead< NUM_ALIVEORDEAD; dead++)
					{
						if ((adh < OFFCARE_INDEX && HIVprograte[adh][1/*On Care*/][curr_cd4][curr_vl][curr_onTreat][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][dead] > 0) 
							|| (adh == OFFCARE_INDEX && HIVprograte[adh][0 /*Off Care*/][curr_cd4][curr_vl][curr_onTreat][curr_resistance][dest_cd4][dest_vl][dest_onTreat][dest_resistance][dead] > 0))
						{
							return false;
						}
					}
				}
			}
		}
	}
	return true;
}

void findNumMenCircumsizedByIntervAtStartOfRun()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	double baselineNumCircumsized;
	double baselineProbCircumsized;
	double adultSuscMales=0;

	numCircumsizedByIntervAtStart = 0;

	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS

		if (age >= ADULTHOOD)
		{
			//We are only applying the program to Susceptible adult males
			//The circumcision intervention is the 4th interv.
			//ANIK
			adultSuscMales = Comp[status][sex][orient][act][age][alc][idu][vl][cd4][res] * interv[4].who[status][sex][orient][act][alc][idu][vl];

			//Circumcision is the 7th pathway
			baselineProbCircumsized = 1 - baselineProb[7][status][sex][orient][act][alc][idu];	

			//Of all of the males who we may wish to circumcize with the program (male, Susc, HIV negative)
			//the baselineNumCircumsized is the number of those who are already circumsized.
			baselineNumCircumsized = baselineProbCircumsized * adultSuscMales;

			//The circumcision intervention is intervention #4.  
			//The circumcision path is pathway #7.
			//The number of men circumsized by the intervention is the number suitable for circumcision
			//multiplied by one minus the intervention effect size. 
			numCircumsizedByIntervAtStart += ( adultSuscMales - baselineNumCircumsized) * (1 - interv_effect_size[4][7]);
		}

	END_LOOP_THROUGH_ALL_COMPARTMENTS 
}

void alcoholIntervention()
{
	int status, sex, orient, act, age, alc, idu, vl, cd4, res;
	int i, p;

	double num_alcs;		//temporary storage variable
	double final_effect_size[NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE][NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL], effect_size;

	p = 10;		//Pathway 10 is the unhealthy alcohol use pathway.	

	LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS 

		//intialize final effect size to no effect
		final_effect_size[p][status][sex][orient][act][alc][idu][vl] = 1;

		//Find the cumulative effect over all interventions that activate that pathway
		for (i=0; i < NUM_INTERVENTIONS; i++)
		{
			//If an intervention is turned on that activates this pathway and the effect sizes is not 1 (implying there actually is an effect)
			if (basket[i] == true && !(isEqual(interv[i].effect_size[p], 1)))
			{
				//To find the actual effect size for an intervention, a who of '0' should have no effect (resulting in an effect size of 1)
				//A who of '1' applies the effect to everyone in the compartment, hence the full effect size applies
				//A who of some value in between gets partial effect, so interpolate!
				//linearInterpolate(x0, y0, x1, y1, x)  applied as....
				//linearInterpolate(who 0, effect size 1, who 1, full intervention effect size, who for this intervention

				effect_size = linearInterpolate(0, 1, 1, interv[i].effect_size[p], interv[i].who[status][sex][orient][act][alc][idu][vl]);

				if (MODE == OPTIMISTIC) final_effect_size[p][status][sex][orient][act][alc][idu][vl] *= effect_size;			//Optimistic

				else final_effect_size[p][status][sex][orient][act][alc][idu][vl] = min(final_effect_size[p][status][sex][orient][act][alc][idu][vl], effect_size);		//Pessimistic
			}
		}

	END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS

	//Apply the final_effect_sizes found above

	LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS

	if (age >= ADULTHOOD)
	{
		//We use the alc==ALC because we want to move people between alc compartments once, not for
		//every value of "alc" in the loop.  
		//We test for a final_effect_size < 1 so that we only enter this code for the
		//targeted group.  (Note that if the group is targeted, the final_effect_size
		//would be less than 1.)
		if (alc==ALC && final_effect_size[p][status][sex][orient][act][ALC][idu][vl] < 1)
		{
			//Save the original number of alcs because we are going to change that 
			//value before we need to use it.  
			num_alcs = Comp[status][sex][orient][act][age][ALC][idu][vl][cd4][res];

			//Successfully treated alcs are moved to ALC_TREATED_CURED
			Comp[status][sex][orient][act][age][ALC_TREATED_CURED][idu][vl][cd4][res] += num_alcs * (1 - final_effect_size[p][status][sex][orient][act][ALC][idu][vl]);

			//Unsuccessfully treated alcs are moved to ALC_TREATED_NOT_CURED
			Comp[status][sex][orient][act][age][ALC_TREATED_NOT_CURED][idu][vl][cd4][res] += num_alcs * final_effect_size[p][status][sex][orient][act][ALC][idu][vl];

			//Since we attempt treatment on all alcs, none are left in the ALC compartment.
			Comp[status][sex][orient][act][age][ALC][idu][vl][cd4][res] = 0;
		}
	}

	END_LOOP_THROUGH_ALL_COMPARTMENTS
}


void initialize_willPartnerWithMatrix()
{
	int sex, orient, p_sex, p_orient;

	//initialize the array to zero
	memcpy(willPartnerWith, zerod_array, sizeof(willPartnerWith));

	if (debug) 
	{
		cout<<"willPartnerWith array: "<<endl<<endl;
		cout<<"MS\tMG\tMB\tWS\tWG\tWB\t"<<endl;
	}

	for (sex=0; sex<NUM_SEX; sex++)
	{
		for (orient=0; orient<NUM_ORIENT; orient++)
		{
			for (p_sex=0; p_sex<NUM_SEX; p_sex++)
			{
				for (p_orient=0; p_orient<NUM_ORIENT; p_orient++)
				{		
					//Gay with Gay or Gay with Bi
					if (sex==p_sex && orient==GAY && (p_orient==GAY || p_orient==BI)) 					
						willPartnerWith[sex][orient][p_sex][p_orient] = true;

					//Straight with Straight or Straight with Bi
					if (sex!=p_sex && orient==STRAIGHT && (p_orient==STRAIGHT || p_orient==BI)) 					
						willPartnerWith[sex][orient][p_sex][p_orient] = true;

					//Bi with Bi or Bi with Gay or Bi with Straight
					if (orient==BI && (p_orient==BI || (sex!=p_sex && p_orient==STRAIGHT) || (sex==p_sex && p_orient==GAY))) 					
						willPartnerWith[sex][orient][p_sex][p_orient] = true;

					if (debug) cout<< willPartnerWith[sex][orient][p_sex][p_orient]<<"\t";
				}
			}
			if (debug) cout<<endl;
		}
	}
	if (debug) cout<<endl<<endl;
}

void print_time()
{
	time_t now;
	struct tm *current;
	now = time(0);
	current = localtime(&now);
	cout << "Time   hour: " << current->tm_hour << "  mins: " << current->tm_min << "  sec: " << current->tm_sec << endl << endl;
}

void preventBalanceError()
{
	if (PROP_GAY_M==1 || PROP_GAY_F==1 || PROP_STRAIGHT_M==1 || PROP_STRAIGHT_F==1)
	{
		cout<<endl<<"The model will not run with any gender being all straight or all gay.  Exiting."<<endl;
		exit(0);
	}
}

void calculate_propIDUinAct()
{

	int sex, orient, act, age;
	double propSexInPop[NUM_SEX]; //starting proportion of people in each sex
	double propActInPop[NUM_ACT]; //starting proportions of people in each activity group 

	/*initialize the proportion of injection drug users within each activity group using the proportion of
	IDU in the entire population (PROP_IDU), the proportion of IDU in each activity group (PROP_ACT_0_IDU...)
	and the proportion of people starting in each activity group (using the PROP_GAY_M, PROP_BI_M etc...)
	for calculations see "Dropbox (NYU CEDS)/SOLVE HIV Modeling/India Models/Transmission code/propIDUinAct_calculations.xlsx"*/
	//1. Zero out matrices
	for (sex = 0; sex < NUM_SEX; sex++) propSexInPop[sex] = 0;
	for (act = 0; act < NUM_ACT; act++) propActInPop[act] = 0;

	//2. Find the proportion of males and females in population
	for (sex = 0; sex < NUM_SEX; sex++)
	{
		for (age = 0; age < NUM_AGE; age++)
		{
			propSexInPop[sex] += prop_sex_age[sex][age];
		}
	}

	//3. Find the proportion of people in the entire population in each activity group
	for (act = 0; act < NUM_ACT; act++)
	{
		for (sex = 0; sex < NUM_SEX; sex++)
		{
			for (orient = 0; orient < NUM_ORIENT; orient++)
			{
				propActInPop[act] += propSexOrientInAct[sex][orient][act] * propSexInPop[sex] * propOrientInPop[sex][orient];
			}
		}
	}

	//4. Find the proportion of IDU in each activity group
#if (ONE_WAY_SENSITIVITY_ANALYSIS && (RUN_NUM==7 || RUN_NUM==8)) || PROBABILISTIC_RUN
	#if PROBABILISTIC_RUN
			
			propIDUinAct[0][IDU] = (PROP_ACT_0_IDU  * prob_pull[3]) / propActInPop[0];
			propIDUinAct[1][IDU] = (PROP_ACT_1_IDU  * prob_pull[3]) / propActInPop[1];
			propIDUinAct[2][IDU] = (PROP_ACT_2_IDU  * prob_pull[3]) / propActInPop[2];
			propIDUinAct[3][IDU] = (PROP_ACT_3_IDU  * prob_pull[3]) / propActInPop[3];
	#else
			propIDUinAct[0][IDU] = (PROP_ACT_0_IDU  * sensitivity_values[RUN_NUM]) / propActInPop[0];
			propIDUinAct[1][IDU] = (PROP_ACT_1_IDU  * sensitivity_values[RUN_NUM]) / propActInPop[1];
			propIDUinAct[2][IDU] = (PROP_ACT_2_IDU  * sensitivity_values[RUN_NUM]) / propActInPop[2];
			propIDUinAct[3][IDU] = (PROP_ACT_3_IDU  * sensitivity_values[RUN_NUM]) / propActInPop[3];
	#endif
#else 
		//max proportion of IDU 
		propIDUinAct[0][IDU] = (PROP_ACT_0_IDU  * PROP_IDU ) / propActInPop[0];
		propIDUinAct[1][IDU] = (PROP_ACT_1_IDU  * PROP_IDU ) / propActInPop[1];
		propIDUinAct[2][IDU] = (PROP_ACT_2_IDU  * PROP_IDU ) / propActInPop[2];
		propIDUinAct[3][IDU] = (PROP_ACT_3_IDU  * PROP_IDU ) / propActInPop[3];
#endif 

	propIDUinAct[0][NON_IDU] = 1 - propIDUinAct[0][IDU];
	propIDUinAct[1][NON_IDU] = 1 - propIDUinAct[1][IDU];
	propIDUinAct[2][NON_IDU] = 1 - propIDUinAct[2][IDU];
	propIDUinAct[3][NON_IDU] = 1 - propIDUinAct[3][IDU];

	//if the proportion of IDU in the population making up that activity group is greater than
	//the total proportion of all people in that activity group, give an error and either reduce the 
	//PROP_ACT_IDU level for that activity group or decrease the PROP_IDU level
	for (act = 0; act < NUM_ACT; act++)
	{
		if (propIDUinAct[act][IDU] > 1)
		{
			cout << "The number of IDU in activity group " << act << " is higher than the number of people in that group. Exiting." << endl;
			exit(0);
		}
	}
}

double lognormal(double mu, double lower, double upper, long int seed_input, int max_val)
{
	double se = (upper - lower) / (1.96 * 2); //convert the CI to a SE 
	double v = pow(se, 2); 
	double log_mu1 = log(mu / sqrt(1 + v / pow(mu, 2))); 
	double log_mu = log(pow(mu, 2) / sqrt(v + pow(mu, 2)));
	//double log_std = sqrt(log(pow((se / mu), 2) + 1));
	double log_se = (log(upper) - log(lower)) / (1.96 * 2); 
	double temp = normal(log_mu, log_se, &seed_input); 
	double exp_result = exp(temp); 
    if (exp_result > max_val) exp_result = max_val; 
	return exp_result; 
}


#if COST_EFFECTIVENESS_SENSITIVITY || PROBABILISTIC_RUN
//The next two functions are for finding the optimum basket

//Turn off all interventions in each basket array
void zero_all_baskets()
{
	for (basket_count = 0; basket_count<TOTAL_NUM_BASKETS; basket_count++)
	{
		for (int i = 0; i<NUM_INTERVENTIONS; i++)
		{
			array_of_baskets[basket_count][i] = 0;
		}
	}
}

void do_loop(int loops_remaining)
{
	for (loop_var[loops_remaining] = 0; loop_var[loops_remaining]<2; loop_var[loops_remaining]++)
	{
		if (loops_remaining >0) do_loop(loops_remaining - 1);

		else
		{
			for (int l = 0; l<NUM_INTERV_TO_COMBINE; l++)
			{
				array_of_baskets[basket_count][intervs_to_combine[l]] = loop_var[l];
			}

			basket_count++;
		}
	}
}

#endif










