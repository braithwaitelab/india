//References to "A&G" refer to 
//Garnett GP, Anderson RM. Balancing sexual partnerships in an age and activity stratified model 
//of HIV transmission in heterosexual populations. IMA J Math Appl Med Biol. 1994;11(3):161-92.

double Comp[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];		//the number in each compartment
double dPop[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];		//used in differentiating

double lambda[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU];				//the force of infection for someone in a given
																//sex, age, activity and alc group

bool debug=DEBUGG;

double pInfectious;				//probability partner is infectious
double birth_rate;				//birth rate
double treat_rate;				//rate of movement from infected to treated state

double t = -DT;					//the current step, initialized for one time step early so
//when we print the initial conditions, the time will
//show up properly

//in order to balance the partnerships, the partner change rate is modified.
//conc_given_expanded puts the given partner change rate matrix in the same form as the 
//later form, just to make it easier to switch between the two.
double conc_given_expanded[NUM_SEX][NUM_ACT][NUM_AGE][NUM_ACT][NUM_AGE];

//the probability matrix
double prob_matrix[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
bool willPartnerWith[NUM_SEX][NUM_ORIENT][NUM_SEX][NUM_ORIENT];	//with one sex/orientation partner with another sex/orientation

double num_partnerships_offered_by_group[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];

double nSexOrientActivityAge[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];		

double balance_matrix[NUM_ACT][NUM_AGE][NUM_ACT][NUM_AGE];

ofstream outfile, basketfile, startAndEnd, my_pop_file;


double fertRateScaledUp[NUM_AGE];			//fertility rates scaled up so that abstinent group doesn't drag them down

double relChangeRateAct[NUM_ACT];			//relative partner acquisition rate across activity groups

double relChangeRateAge[NUM_AGE];			//relative partner acquisition rate across age groups
//Note, we only include adults here

double age_births_deaths[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];

double age_deaths[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];			//the number of age deaths from all compartments this time period
double hiv_deaths[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_ALC][NUM_IDU][NUM_VL][NUM_CD4][NUM_RES];			//the number of AIDS deaths from all compartments this time period

//Chris's stuff for reading in state transitions
int statecount[NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_ALIVEORDEAD];
double HIVprograte[NUM_ADHERENCE][CARESTATUS][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_ALIVEORDEAD];		
double costCareTreat[NUM_ADHERENCE][CARESTATUS][NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES];  //Costs are now dependent only on the source compartment, not the source-destination combination
int	staterowsum[NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES],statecolsum[NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES][NUM_ALIVEORDEAD];
bool headerdone = 0;
bool state_is_used[NUM_CD4][NUM_VL][NUM_TREATORNOT][NUM_RES];
FILE *stateprob, *stateprobhuman, *sumsofrows, *sumsofcols, *graphs, *theunion;


//variables for apportioning newly infecteds into vl/cd4 buckets using a distribution
double cd4pts[NUM_CD4-1] = {50, 200, 350, 500};	//these values are taken from our CD4 categories
double vlpts[NUM_VL-1] = {2.5, 3.5, 4.5, 5.5};	//these values are taken from our VL categories
double probcd4vl[NUM_SEX][NUM_CD4][NUM_VL];					//the probability of going into each cd4/vl combination bucket


double dur[NUM_ACT] = {NOT_APPLICABLE, DUR_ACT_1, DUR_ACT_2, DUR_ACT_3};	//the duration of partnerships before compromise (meaningless for abstinent group!)	
double conc_given[NUM_ACT] = {NOT_APPLICABLE, MEDIAN_CONC_PART_ACT_1, MEDIAN_CONC_PART_ACT_2, MEDIAN_CONC_PART_ACT_3};	//the concurrency of partnerships before compromise (meaningless for abstinent group!)				
	
double conc_comp[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];	//concurrency modified after balancing according to c_current
double dur_comp[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
double num_partnerships_matrix[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
double desired_overall_freq[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];		//desired frequency (can vary with concurrency)
double freq_per_partnership[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
double acts_per_partnership[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];
double change_rate_comp[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_AGE];

//For zeroing out arrays in InitializeArrays().  //Note also that this array has the same indices as preCalculatedAlphMultArray which is the biggest array in the model
double zerod_array[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL][NUM_ALIVE_STATE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL] = { 0 };		//This comment is outdated, but I'm leaving it and it's code in because it can't hurt.  It seems that without this "={0}", the model is seg-faulting on the cluster!

double age_out_rate[NUM_AGE];					//This is the rate of aging out of one age group and into the next
double numSexuallyActiveWomen[NUM_AGE];
double numAbstinentWomen[NUM_AGE];

//for calibration
//the proportion of people to start in each age/sex group
//From "TRANSMISSION MODEL INPUTS", tab "pop breakdown"
double prop_sex_age[NUM_SEX][NUM_AGE] = { { 0.0740, 0.0674, 0.0694, 0.0485, 0.0665, 0.0385, 0.0385, 0.0385, 0.0385, 0.0384 }, { 0.0702, 0.0641, 0.0655, 0.0455, 0.0615, 0.0350, 0.0350, 0.0350, 0.0350, 0.0350 } };

//the proportion in each age,sex that are infected
//From "TRANSMISSION MODEL INPUTS", tab "prop inf at start"
double prop_inf_sex_age[NUM_SEX][NUM_AGE] = { { 0.000107, 0.000107, 0.000107, 0.000107, 0.002036, 0.004670, 0.006857, 0.005679, 0.004393, 0.005143 }, { 0.000750, 0.000750, 0.000750, 0.000750, 0.001821, 0.003000, 0.004821, 0.002464, 0.002036, 0.001821}}; 

//The initial vl/cd4 distribution
//From "TRANSMISSION MODEL INPUTS", tab "init_INF_cd4_vl_distr"
double init_INF_cd4_vl_distr[NUM_VL][NUM_CD4] = { {0.000000, 0.000000, 0.003968, 0.003968, 0.031746}, {0.000000, 0.000000, 0.012698, 0.012698, 0.073016}, {0.005952, 0.017857, 0.038095, 0.038096, 0.176191}, {0.013095, 0.039286, 0.061905, 0.061904, 0.161905}, {0.016667, 0.050000, 0.059524, 0.034524, 0.036905} };
double init_HVL_cd4_vl_distr[NUM_VL][NUM_CD4] = { {0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 0.000000, 0.008333, 0.008333}, {0.000000, 0.000000, 0.000000, 0.016667, 0.016667}, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000} };

//end of calibration vars

//from (Garcia, 1999)
double motherToChild[NUM_VL] = {0.000, 0.148, 0.338, 0.604, 0.725};  

//the first sex is the gender of the person we're calculating Beta for
double alpha_raw[NUM_SEX][NUM_SEX] = { { (0.23 * (0.008312 + 0.0062) / 2), .00042 }, { 0.00081, 0 } }; //{{[M][M], [M][F]}, {[F][M], [F][F]}} Taken from NYC model
double alpha_vl_mults[NUM_VL] = {.03, 0.33, 1.16, 1.57, 1.60};
double alpha[NUM_SEX][NUM_SEX][NUM_VL]; //the second sex is the one doing the infecting


//Proportion of mothers on mother-to-child prevention meds
double pmtct_meds = PMTCT_DEFAULT;


//fertility rates per capita (per female), (early groups are too young to be sexually active)
// From TRANSMISSION MODEL INPUTS fertRate tab 
//These rates run from 1997-2014 (18 years)
double fertRate[18][NUM_AGE] ={ 
{ 0, 0, 0, 0.0880, 0.2300, 0.1740, 0.0960, 0.0470, 0.0190, 0.0050 },
{ 0, 0, 0, 0.0853, 0.2277, 0.1717, 0.0940, 0.0457, 0.0189, 0.0051 },
{ 0, 0, 0, 0.0815, 0.2253, 0.1698, 0.0917, 0.0441, 0.0181, 0.0049 },
{ 0, 0, 0, 0.0778, 0.2230, 0.1680, 0.0895, 0.0425, 0.0174, 0.0048 },
{ 0, 0, 0, 0.0740, 0.2206, 0.1661, 0.0872, 0.0409, 0.0199, 0.0046 },
{ 0, 0, 0, 0.0722, 0.2184, 0.1641, 0.0854, 0.0395, 0.0152, 0.0046 },
{ 0, 0, 0, 0.0666, 0.2160, 0.1623, 0.0826, 0.0377, 0.0151, 0.0044 },
{ 0, 0, 0, 0.0628, 0.2136, 0.1604, 0.0803, 0.0361, 0.0143, 0.0042 },
{ 0, 0, 0, 0.0591, 0.2113, 0.1586, 0.0781, 0.0345, 0.0136, 0.0041 },
{ 0, 0, 0, 0.0553, 0.2089, 0.1567, 0.0758, 0.0329, 0.0128, 0.0039 },
{ 0, 0, 0, 0.5060, 0.2066, 0.1550, 0.0734, 0.0312, 0.0112, 0.0038 },
{ 0, 0, 0, 0.0479, 0.2043, 0.1529, 0.0712, 0.0297, 0.0113, 0.0037 }, 
{ 0, 0, 0, 0.0441, 0.2019, 0.1510, 0.0689, 0.0281, 0.0105, 0.0035 }, 
{ 0, 0, 0, 0.0404, 0.1996, 0.1492, 0.0667, 0.0265, 0.0098, 0.0034 }, 
{ 0, 0, 0, 0.0366, 0.1972, 0.1473, 0.0644, 0.0249, 0.0090, 0.0032 }, 
{ 0, 0, 0, 0.0329, 0.1949, 0.1454, 0.0621, 0.0233, 0.0083, 0.0031 }, 
{ 0, 0, 0, 0.0292, 0.1926, 0.1435, 0.0598, 0.0217, 0.0075, 0.0030 }, 
{ 0, 0, 0, 0.0254, 0.1902, 0.1416, 0.0575, 0.0201, 0.0067, 0.0028 }
}; 

//age related mortality rates
//From Sewankambo, 2000, AIDS (Table 1, p 2394)
double ageMortRate[NUM_SEX][NUM_AGE]={					
{0.048, 0.001, 0.001, 0.001, 0.002, 0.002, 0.003, 0.004, 0.006, 0.008},
{0.05, 0.001, 0.001, 0.001, 0.002, 0.002, 0.002, 0.002, 0.003, 0.004}};

double totalDeadAIDS_no_interv, totalDeadAIDS;	

// the proportion of men/women in each status/sex/activity combo with alcohol abuse or IDU
//note that for the prop alc and prop non-alc for each status/sex/activity combo should sum to 1.
double prop_alc[NUM_STATUS][NUM_SEX][NUM_ACT][NUM_ALC];
double prop_idu[NUM_STATUS][NUM_SEX][NUM_ACT][NUM_IDU];

//the proportion of the sexually active infected population in each resistance category
double propRes[NUM_RES];

//the rates of flow from the infected to detected state and from the detected to on care state
double inf_to_det_rate[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];
double det_to_care_rate[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];
double adherence[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];

double totalkidsinmodel, totaladultsinmodel;
double totalPeopleInModel, totalinfectedsinmodel;

double incidence, HIVdeathrate, alcohol_incidence, num_new_infections_alcohol, IDU_incidence, num_new_infections_IDU;

double costOfCareAndTreatThisCycle;
double costOfCareAndTreat, discounted_costOfCareAndTreat;

//We will maintain LYs, QALYs, and disc LYs and disc QALYs for 0)  All pats, 1)  All HIV+  2)  All On Care + Treat
double total_life_years[3], total_qa_life_years[3], total_disc_life_years[3], total_disc_qa_life_years[3];


//The following are just to make test scenarios easier
	bool turn_Off_HIV_Mortality;
	bool turn_Off_Transmission_On_Treatment;
	bool new_Infections_Go_Immediately_On_Treat;
	int cd4_treat_threshold;


//Kendall Bryant version
#if KENDALL_BRYANT_ALCOHOL_ANALYSIS
double analysis[3][3]= {
	{1.286,1.19,2.30},
	{2.33,1.67,3.10},
	{1.72,1.17,5.96}
}; // analysis
double analysis_effects[7] = {.5,.6,.7,.8,.9,.95, 1}; // analysis
#endif
//Variables needed for one way sensitivity
double yintercept = Y_INTERCEPT;


double sumOfNumInfectionsPerInfected, sumOfCostOfCareAndTreatPerInfected, sumOfNumDeathsPerInfected;

double num_infections[MAX_NUM_BASKETS];		//the number of new infections in each basket being run
int basket_num = NOT_APPLICABLE;			//initialize to not applicable (N/A gets used as a trigger in running some loops of baskets)
double num_infections_no_interv = 0;		//The number of new infections if no interventions are run.  This doesn't work when we put it in the initialize() function.  Leave here.
double infected_counts[NUM_STATUS];

#if ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION
int run;
double baselineRR_raw_backup_condoms;
double baselineRR_raw_backup_nonadherence;
double baselineRR_raw_backup_STIs;
#endif

bool atLeastOneInterventionTurnedOn;

double prop_alcohol[NUM_SEX] = {PROP_ALCOHOL_M, PROP_ALCOHOL_F};		//The proportion of people with unhealthy alcohol use
double propOrientInPop[NUM_SEX][NUM_ORIENT] ={ { PROP_STRAIGHT_M, PROP_GAY_M, PROP_BI_M }, { PROP_STRAIGHT_F, PROP_GAY_F, PROP_BI_F } }; //proportion with each orientation (straight/gay/bi) in males and females

double propOrientActAlcIDU[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];	//proportion of men/women of each status in orient/act/alc/idu combination (may be modified by an intervention)
double propSexOrientInAct[NUM_SEX][NUM_ORIENT][NUM_ACT];  //proportion of men/women in each orientation in each activity group
double propIDUinAct[NUM_ACT][NUM_IDU];   //Within each activtiy group, the proportion who are IDU users. These numbers are made up!  make sure we fill in correct values


#if OPTIMIZING_ALCOHOL_INTERVENTION || (OPTIMIZING_PREVENTION_PORTFOLIO && (SET_RR_ALC > 0))
double condomRR, adhereRR, stiRR;
double is_baseline_run;
#endif

double prop_detected_on_treatment, prop_infected_on_treatment, prop_infected_detected;

#if ALCOHOL_SURFACE
double alcohol_cost[5] = {1,10,20,50,100};
double alcohol_effect_size[5] = {0.95, 0.75, 0.55, 0.35, 0.15};

//RR modifiers for ALCOHOL_SURFACE_ANALYSIS "Values come from spreadsheet New values for surface RR mults 2013.03.28.xlsx"
double alcohol_rr_modifier[3][12] = {
	{1.0, 0.777, 0.817, 0.858, 0.898, 0.938, 0.978, 1.018, 1.058, 1.098, 1.138, 1.178}, //Condom
	{1.0, 0.429, 0.672, 0.915, 1.157, 1.400, 1.643, 1.886, 2.129, 2.372, 2.615, 2.858}, //Adherence
	{1.0, 0.580, 1.145, 1.710, 2.274, 2.839, 3.404, 3.968, 4.533, 5.098, 5.662, 6.227}  //STI
};

double condomRR, adhereRR, stiRR;
int cost_index, effect_index, rr_index;
#endif

//EFFICIENCY pre-calculated arrays
double preCalculatedAlphMultArray[NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL][NUM_ALIVE_STATE][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL];
double preCalcAdhStratifiedHIVProgRate[NUM_ALIVE_STATE][NUM_CD4][NUM_VL][NUM_RES][NUM_CD4][NUM_VL][2/*On treat or not*/][NUM_RES][2/*Alive or Dead*/][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];
double preCalcAdhStratifiedCostCareTreat[NUM_ALIVE_STATE][NUM_CD4][NUM_VL][NUM_RES][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];

//For circumcision intervention
double numCircumsizedByIntervAtStart, numCircumsizedByIntervThisCycle;

//Injection Drug Use variables
double alpha_IDU_raw = .0036;								//From Cardo, 1997;
double alpha_IDU[NUM_VL];
double num_risky_idu_users;

#if COST_EFFECTIVENESS_SENSITIVITY  || ONE_WAY_SENSITIVITY_ANALYSIS || PROBABILISTIC_RUN
long int seed;
int run_num;
#endif

#if COST_EFFECTIVENESS_SENSITIVITY || PROBABILISTIC_RUN
//to figure out which of the baskets to run based on an input 
int intervs_to_combine[NUM_INTERV_TO_COMBINE] = {10, 19, 20, 21, 22, 23, 24 };  
int loop_var[NUM_INTERV_TO_COMBINE];
bool array_of_baskets[TOTAL_NUM_BASKETS][NUM_INTERVENTIONS];
int basket_count;
int basket_use; 
int loops_remaining;

#endif 

/***************************Settings for One Way Sensitivity**********************/
#if ONE_WAY_SENSITIVITY_ANALYSIS 
string sensitivity_name[NUM_SENSITIVITY_TESTS] = { "Baseline",
"Age sexual debut, min", "Age sexual debut, max", //(1-2)
"Proportion of ALCS male, min", "Proportion of ALCS male, max", //(3-4)
"Proportion of ALCS female, min", "Proportion of ALCS female, max", //(5-6)
"Proportion IDU, min", "Proportion IDU, max", //(7-8)
"Proportion gay men, min", "Proportion gay men, max", //(9-10)
"Proporition bi males, min", "Proportion bi males, max", //(11-12)
"Mortality multiplier <5 year olds, min", "Mortality multiplier <5 years old, max", //(13-14)
"Proportion pregnant women on PMTCT in 2014, min", "Proportion pregnant women on PMTCT in 2014, max", //(15-16) 
"Proportion straight males in ACT0, min", "Proportion straight males in ACT0, max", //(17-18)
"Proportion straight males in ACT2, min", "Proportion straight males in ACT2, max", //(19-20)
"Proportion straight males in ACT3, min", "Proportion straight males in ACT3, max", //(21-22)
"Proportion straight females in ACT0, min", "Proportion straight females in ACT0, max", //(23-24)
"Proportion straight females in ACT2, min", "Proportion straight females in ACT2, max", //(25-26)
"Proportion straight females in ACT3, min", "Proportion straight females in ACT3, max", //(27-28)
"Proportion gay males in ACT1, min", "Proportion gay males in ACT2, max", //(29-30)
"Proportion gay males in ACT3, min", "Proportion gay males in ACT3, max", //(31-32)
"Proportion gay females in ACT0, min", "Proportion gay females in ACT0, max", //(33-34)
"Proportion gay females in ACT2, min", "Proportion gay females in ACT2, max", //(35-36)
"Proportion gay females in ACT3, min", "Proportion gay females in ACT3, max", //(37-38)
"Proportion bi males in ACT1, min", "Proportion bi males in ACT1, max",  //(39-40)
"Proportion bi males in ACT3, min", "Proportion bi males in ACT3, max", //(41-42)
"Proportion bi females in ACT0, min", "Proportion bi females in ACT0, max", //(43-44)
"Proportion bi females in ACT2, min", "Proportion bi females in ACT2, max", //(45-46)
"Proportion bi females in ACT3, min", "Proportion bi females in ACT3, max", //(47-48)
"Duration of ACT1 partnership, min", "Duration of ACT1 partnership, max", //(49-50)
"DUration of ACT2 partnership, min", "Duration of ACT2 partnership, max", //(51-52)
"Duration of ACT3 partnership, min", "Duration of ACT3 partnership, max", //(53-54)
"Median number concurrent partnerships ACT2, min", "Median number concurrent partnerships ACT2, max", //55-56
"Median number concurrent partnerships ACT3, min", "Median number concurrent partnerships ACT3, max", //57-58
"FSW condom multiplier, min", "FSW condom multiplier, max", //59-60
"EPS, min", "EPS, max", //61-62
"CD4 mean and std females, min", "CD4 mean and std females, max", //63-64
"VL mean females, min", "VL mean females, max", //65-66
"CD4 mean and std males, min", "CD4 mean and std males, max", //67-68
"VL mean males, min", "VL mean males, max", //69-70
"Utility CD4 less than 50, min", "Utility CD4 less than 50, max", //71-72
"Utility CD4 50 to 200, min", "Utility CD4 50-200, max", //73-74
"Utility CD4 above 200, min", "Utility CD4 above 200, max", //75-76
"Delta utility with HAART, min", "Delta utility with HAART, max", //77-78
"Discount rate, min", "Discount rate, max", //79-80
"alpha_raw male infecting female, min", "alpha_raw male infected female, max", //81-82
"alpha_raw female infecting male, min", "alpha_raw female infecting male, max", //83-84
"alpha IDU raw, min", "alpha IDU raw, max", //85-86
"alpha raw male infecting male, min", "alpha_raw male infecting male, max", //87-88
"number of shared injections, min", "number of shared injections, max", //89-90
"idu partner change rate, min", "idu partner change rate, max", //91-92
"drug alpha scale, min", "drug alpha scale, max", //93-94
"Median number concurrent partnerships ACT1, min", "Median number concurrent partnerships ACT1, max", //95-96
"proportion of condom nonuse in the general population, min", "proportion of condom nonuse in the general population, max", //97-98
"RR of condon nonuse female vs. male, min", "RR of condon nonuse female vs. male, max", //99-100
"RR of condom nonuse gay vs. straight, min","RR of condom nonuse gay vs. straight, max", //101-102
"RR of condom nonuse risk 2 vs. risk 1, min", "RR of condom nonuse risk 2 vs. risk 1, max", //103-104
"RR of condom nonuse risk3 vs. risk 1, min", "RR of condom nonuse risk 3 vs. risk 1, max", //105-106
"RR of condom nonuse alcohol/drug use vs. not, min", "RR of condom nonuse  alcohol/drug use vs. not, max", //107-108
"RR of condom nonuse IDU 1 vs. IDU 0, min", "RR of condom nonuse IDU 1 vs. IDU 0, max", //109-110
"RR of condom nonuse HIV inf detected vs HIV-, min", "RR of condom nonuse HIV inf detected vs HIV-, max", //111-112
"RR of condom nonuse HIV inf on care vs HIV-, min", "RR of condom nonuse HIV inf on care vs HIV-, max", //113-114
"RR of condom nonuse HIV inf treated vs HIV-, min", "RR of condom nonuse HIV inf treated vs HIV-, max", //115-116
"proportion of not being tested for HIV in the general population, min", "proportion of not being tested for HIV in the general population, max", //117-118
"RR of not being tested for HIV gay vs. straight, min", "RR of not being tested for HIV gay vs. straight, max", //119-120
"RR of not being tested for HIV bi vs. straight, min", "RR of not being tested for HIV bi vs. straight, max", //121-122
"RR of not being tested for HIV risk 2 vs. risk 1, min", "RR of not being tested for HIV risk 2 vs. risk 1, max", //123-124
"RR of not being tested for HIV risk3 vs. risk 1, min", "RR of not being tested for HIV risk 3 vs. risk 1, max", //125-126
"RR of not being tested for HIV IDU 1 vs. IDU 0, min", "RR of not being tested for HIV IDU 1 vs. IDU 0, max", //127-128
"proportion of ART nonadherence in the general population, min", "proportion of ART nonadherence in the general population, max", //129-130
"RR of ART nonadherence alcohol/drug use vs. not, min", "RR of ART nonadherence alcohol/drug use vs. not, max", //131-132
"RR of ART nonadherence IDU 1 vs. IDU 0, min", "RR of ART nonadherence IDU 1 vs. IDU 0, max", //133-134
"proportion of untreated STI in the general population, min", "proportion of untreated STI in the general population, max", //135-136
"RR of untreated STI female vs. male, min", "RR of untreated STI female vs. male, max", //137-138
"RR of untreated STI risk3 vs. risk 1, min", "RR of untreated STI risk 3 vs. risk 1, max", //139-140
"RR of untreated STI alcohol/drug use vs. not, min", "RR of untreated STI alcohol/drug use vs. not, max", //141-142
"RR of untreated STI IDU 1 vs. IDU 0, min", "RR of untreated STI IDU 1 vs. IDU 0, max", //143-144
"proportion of not being circumcised in the general population, min", "proportion of not being circumcised in the general population, max", //145-146
"RR of not being circumcised HIV inf detected vs HIV-, min", "RR ofnot being circumcised HIV inf detected vs HIV-, max", //147-148
"RR of not being circumcised HIV inf on care vs HIV-, min", "RR of not being circumcised HIV inf on care vs HIV-, max", //149-150
"RR of not being circumcised HIV inf unknown vs HIV-, min", "RR of not being circumcised HIV inf unknown vs HIV-, max", //151-152
"RR of not being circumcised HIV inf treated vs HIV-, min", "RR of not being circumcised HIV inf treated vs HIV-, max", //153-154

"Proportion gay women, min", "Proportion gay women, max", //155-156
"Proporition bi women, min", "Proportion bi women, max", //157-158
"probability of mother to child transmission (VL0), min", "probability of mother to child transmission (VL0), max", //159-160 
"probability of mother to child transmission (VL1), min", "probability of mother to child transmission (VL1), max", //161-162
"probability of mother to child transmission (VL2), min", "probability of mother to child transmission (VL2), max", //163-164
"probability of mother to child transmission (VL3), min", "probability of mother to child transmission (VL3), max", //165-166
"probability of mother to child transmission (VL4), min", "probability of mother to child transmission (VL4), max", //167-168
"CD4 std females, min", "CD4 std females, max", //CD4SD_F 169-170
"CD4 std males, min", "CD4 std males, max"//CD4SD_M 172-173

};

double sensitivity_values[NUM_SENSITIVITY_TESTS] = { 0, //RUN 0
17, 21, //AGE_SEXUAL_DEBUT (RUNS 1-2)
0.0925, 0.2775, //PROP_ALCOHOL_M (RUNS 3-4)
0.0096, 0.0224, //PROP_ALCOHOL_F (RUNS 5-6)
0.00004715, 0.00007685, //PROP_IDU (RUNS 7-8)
0.00059322, 0.00096694, //PROP_GAY_M (RUNS 9-10)
0.000394884, 0.000487116,  //PROP_BI_M (RUNS 11-12)
3, 9, //MORT_MULT_CHILDREN_NOT_ON_TREAT (RUNS 13-14)
0.74, 0.94, //PROP_PREGNANT_WOMAN_ON_PMTCT_2014 (RUNS 15-16)

0.21, 0.41, //PROP_ACT_0_STRAIGHT_M (RUNS 17-18)
0.1298, 0.21164, //PROP_ACT_2_STRAIGHT_M (RUNS 19-20)
0.02, 0.0326, //PROP_ACT_3_STRAIGHT_M (RUNS 21-22)

0.1, 0.42, //PROP_ACT_0_STRAIGHT_F (RUNS 23-24)
0.047, 0.07661, //PROP_ACT_2_STRAIGHT_F (RUNS 25-26)
0.00290753, 0.00473927, //PROP_ACT_3_STRAIGHT_F (RUNS 27-28)

0.192, 0.2652, //PROP_ACT_1_GAY_M (RUNS 29-30)
0.0971, 0.1583, //PROP_ACT_3_GAY_M (RUNS 31-32)

0.1, 0.42, //PROP_ACT_0_GAY_F (RUNS 33-34)
0.047, 0.07661, //PROP_ACT_2_GAY_F (RUNS 35-36)
0.00290753, 0.00473927, //PROP_ACT_3_GAY_F (RUNS 37-38)

0.192, 0.2652, //PROP_ACT_1_BI_M (RUNS 39-40)
0.0971, 0.1583, //PROP_ACT_3_BI_M (RUNS 41-42)

0.1, 0.42, //PROP_ACT_0_BI_F (RUNS 43-44)
0.047, 0.07661, //PROP_ACT_2_BI_F (RUNS 45-46)
0.00290753, 0.00473927, //PROP_ACT_3_BI_F (RUNS 47-48)

15, 45, //DUR_ACT_1 (RUNS 49-50)
0.5, 1.5, //DUR_ACT_1 (RUNS 51-52)
0.25, 0.75,//DUR_ACT_1 (RUNS 53-54)
1.5, 4.5,//MEDIAN_CONC_PART_ACT_2 (RUNS 55-56)
5, 15,  //MEDIAN_CONC_PART_ACT_3 (RUNS 57-58)
0.042, 0.126,  //FSW_CONDOM_MULTIPLIER (RUNS 59-60)
0.05, 0.5, //EPS (RUNS 61-62)
294, 994, //CD4MEAN_F (RUNS 63-64)
4, 5, //VLMEAN_F (RUNS 65-66)
294, 994, //CD4MEAN_M (RUNS 67-68)
4, 5, //VLMEAN_M (RUNS 69-70)
0.74,0.84, //UTILITY_CD4_LESS_THAN_50 (RUNS 71-72) 
0.80,0.90, //UTILITY_CD4_LESS_THAN_50 (RUNS 73-74) 
0.89,0.99,//UTILITY_CD4_ABOVE_200 (RUNS 75-76) 
-0.0265,-0.0795, //DELTA_UTILITY_WITH_HAART (RUNS 77-78)
0.02,0.1, //DISC_RATE (RUNS 79-80)

0.0004, 0.0012, //alpha_raw[F][M] (runs 81-82)
0.00021, 0.00063, //alpha_raw[m][f] (runs 83-84)
0.0018, 0.0054, //alpha_IDU_raw (runs 85-86)
0.000834, 0.002503, //alpha_raw[M][M] (runs 87-88)
54, 150, //SHARED_INJECTION_PER_YEAR (runs 89-90)
2.5, 7.5, //IDU_PARTNER_CHANGE_RATE (runs 91-92)
0.125, 0.375, //DRUG_ALPHA_SCALE (runs 93-94)
MEDIAN_CONC_PART_ACT_1*0.5, MEDIAN_CONC_PART_ACT_1*1.5, //MEDIAN_CONC_PART_ACT_2 (RUNS 95-96)

0.73, 0.815, //proportion of condom nonuse  97-98
1.063, 1.113,  //RR of condon nonuse female vs. male 99-100
0.602, 0.658, //RR of condom nonuse gay vs. straight 101-102
0.7, 0.834, //RR of condom nonuse risk 2 vs. risk 1, 103-104
1.086, 1.188, //RR of condom nonuse risk3 vs. risk 1,  105-106
1, 1.58, //RR of condom nonuse alcohol/drug use vs. not 107-108
0.549, 0.697, //RR of condom nonuse IDU 1 vs. IDU 0 109-110
0.4, 0.54, //RR of condom nonuse HIV inf detected vs HIV- 111-112
0.4, 0.54, //RR of condom nonuse HIV inf on care vs HIV- 113-114
0.4, 0.54, //RR of condom nonuse HIV inf treated vs HIV- 115-116

0.96, 0.99, //proportion of not being tested for HIV in the general population 117-118
0.276, 0.349, //RR of not being tested for HIV gay vs. straight 119-120
0.276, 0.349, //RR of not being tested for HIV bi vs. straight 121-122
0.793, 0.819, //RR of not being tested for HIV risk 2 vs. risk 1 123-124
0.221, 0.251, //RR of not being tested for HIV risk3 vs. risk 1 125-126
0.539, 0.666, //RR of not being tested for HIV IDU 1 vs. IDU 0 127-128

0.26, 0.364, //proportion of ART nonadherence in the general population 129-130
1.166, 3.499, //RR of ART nonadherence alcohol/drug use vs. not 131-132
1, 3, //RR of ART nonadherence IDU 1 vs. IDU 0 133-134

0.06, 0.096, //proportion of untreated STI in the general population 135-136
0.5, 1.5, //RR of untreated STI female vs. male 137-138
7.377, 10.328, //RR of untreated STI risk3 vs. risk 1 139-140
1.40, 2.046, //RR of untreated STI alcohol/drug use vs. not 141-142
1.223, 1.630, //RR of untreated STI IDU 1 vs. IDU 0 143-144

0.683, 0.917, //proportion of not being circumcised in the general population 145-146
1.163, 6.667, //RR of not being circumcised HIV inf detected vs HIV- 147-148
1.163, 6.667, //RR of not being circumcised HIV inf on care vs HIV-, 149-150
1.163, 6.667,//RR of not being circumcised HIV inf unknown vs HIV- 151-152
1.163, 6.667, //RR of not being circumcised HIV inf treated vs HIV- 153-154

0.0009, 0.00147, //prop gay women 155-156
0.000188, 0.000232, //prop bi women 157-158

0, 0, //motherToChild[0] 159-160
0.074, 0.222, //motherToChild[1] 161-162
0.169, 0.507, //motherToChild[2] 163-164
0.302, 0.906, //motherToChild[3] 165-166
0.3625, 1.0875, //motherToChild[4] 167-168

65,585,//CD4SD_F 169-170
65,585//CD4SD_M 172-173

};
#endif 

/***************************Settings for Probabilistic Run**********************/
#if PROBABILISTIC_RUN
string probabilistic_name[NUM_PROB_TESTS + NUM_PROB_COST_EFFECT] = {
"Age sexual debut", //0
"Proportion of ALCS male", //1
"Proportion of ALCS female", //2
"Proportion IDU", //3
"Proportion gay men", //4
"Proporition bi males", //5
"Mortality multiplier <5 year olds", //6
"Proportion pregnant women on PMTCT in 2014", //7
"Proportion straight males in ACT0", //8
"Proportion straight males in ACT2", //9
"Proportion straight males in ACT3", //10
"Proportion straight females in ACT0", //11
"Proportion straight females in ACT2", //12
"Proportion straight females in ACT3", //13
"Proportion gay males in ACT1,", //14
"Proportion gay males in ACT3", //15
"Proportion gay females in ACT0", //16
"Proportion gay females in ACT2", //17
"Proportion gay females in ACT3", //18
"Proportion bi males in ACT1", //19
"Proportion bi males in ACT3", //20
"Proportion bi females in ACT0", //21
"Proportion bi females in ACT2", //22
"Proportion bi females in ACT3", //23
"Duration of ACT1 partnership", //24
"DUration of ACT2 partnership", //25
"Duration of ACT3 partnership", //26
"Median number concurrent partnerships ACT2", //27
"Median number concurrent partnerships ACT3", //28
"FSW condom multiplier", //29
"EPS", //30
"CD4 mean males/females", //31, 32
"CD4 std males/females",
"VL mean males/females", //33
"blank", //34, 35
"blank",
"blank", //36
"Utility CD4 less than 50", //37
"Utility CD4 50 to 200", //38
"Utility CD4 above 200", //39
"Delta utility with HAART", //40
"alpha_raw male infecting female", //41
"alpha_raw female infecting male", //42
"alpha IDU raw", //43
"alpha raw male infecting male", //44
"number of shared injections", //45
"idu partner change rate", //46
"drug alpha scale", //47
"proportion of condom nonuse in the general population", //48
"RR of condon nonuse female vs. male", //49
"RR of condom nonuse gay vs. straight", //50
"RR of condom nonuse risk 2 vs. risk 1", //51
"RR of condom nonuse risk3 vs. risk 1", //52
"RR of condom nonuse alcohol/drug use vs. not", //53
"RR of condom nonuse IDU 1 vs. IDU 0", //54
"RR of condom nonuse HIV inf detected vs HIV-", //55
"RR of condom nonuse HIV inf on care vs HIV-", //56
"RR of condom nonuse HIV inf treated vs HIV-", //57
"proportion of not being tested for HIV in the general population", //58
"RR of not being tested for HIV gay vs. straight", //59
"RR of not being tested for HIV bi vs. straight,", //60
"RR of not being tested for HIV risk 2 vs. risk 1", //61
"RR of not being tested for HIV risk3 vs. risk 1", //62
"RR of not being tested for HIV IDU 1 vs. IDU 0", //63
"proportion of ART nonadherence in the general population", //64
"RR of ART nonadherence alcohol/drug use vs. not", //65
"RR of ART nonadherence IDU 1 vs. IDU 0", //66
"proportion of untreated STI in the general population", //67
"RR of untreated STI female vs. male", //68
"RR of untreated STI risk3 vs. risk 1", //69
"RR of untreated STI alcohol/drug use vs. not", //70
"RR of untreated STI IDU 1 vs. IDU 0", //71
"proportion of not being circumcised in the general population", //72
"RR of not being circumcised HIV inf detected vs HIV-", //73
"RR of not being circumcised HIV inf on care vs HIV-", //74
"RR of not being circumcised HIV inf unknown vs HIV-", //75
"RR of not being circumcised HIV inf treated vs HIV-", //76
"Proportion gay women", //77
"Proporition bi women", //78
"probability of mother to child transmission (VL0)", //79
"probability of mother to child transmission (VL1)", //80
"probability of mother to child transmission (VL2)", //81
"probability of mother to child transmission (VL3)", //82
"probability of mother to child transmission (VL4)", //83
"intervention 10 effect size, path 10",//84
"intervention 19 effect size, path 0", //85
"intervention 19 effect size, path 5", //86
"intervention 20 effect size path 0", //87
"intervention 20 effect size, path 5", //88
"intervention 21 effect size path 0", //89
"intervention 21 effect size, path 5", //90
"intervention 22 effect size path 0", //91
"intervention 22 effect size, path 5", //92
"intervention 23 effect size, path 3", //93
"intervention 24 effect size, path 3", //94
"intervention 10 cost", //95
"intervention 19 cost", //96
"intervention 20 cost", //97
"intervention 21 cost", //98
"intervention 22 cost", //99
"intervention 23 cost", //100
"intervention 24 cost" //101

};

double prob_pull[NUM_PROB_TESTS + NUM_PROB_COST_EFFECT] ; //create an empty matrix for all of the data pulls


double prob_values[NUM_PROB_TESTS][3] = { //min, max, std
		{ 17, 21, 1.02 }, //AGE_SEXUAL_DEBUT 0
		{ 0.0925, 0.2775, 0.047}, //PROP_ALCOHOL_M 1
		{ 0.0096, 0.0224, 0.0033}, //PROP_ALCOHOL_F 2
		{ 0.00004715, 0.00007685, 0.00000758 }, //PROP_IDU 3
		{ 0.00059322, 0.00096694, 0.00009534 }, //PROP_GAY_M 4
		{ 0.00028602, 0.00035283, 0.00001704 },  //PROP_BI_M 5
		{ 3, 9, 1.531 }, //MORT_MULT_CHILDREN_NOT_ON_TREAT 6
		{ 0.74, 0.94, 0.051 }, //PROP_PREGNANT_WOMAN_ON_PMTCT_2014 7
		{ 0.21, 0.41, 0.051 }, //PROP_ACT_0_STRAIGHT_M 8
		{ 0.1298435, 0.21164491, 0.02086771 }, //PROP_ACT_2_STRAIGHT_M 9
		{ 0.02, 0.0326, 0.00321429 }, //PROP_ACT_3_STRAIGHT_M 10
		{ 0.1, 0.42, 0.08163 }, //PROP_ACT_0_STRAIGHT_F 11
		{ 0.047, 0.07661, 0.00755 }, //PROP_ACT_2_STRAIGHT_F 12
		{ 0.00290753, 0.00473927, 0.00046728 }, //PROP_ACT_3_STRAIGHT_F 13
		{ 0.192, 0.2652, 0.01867 }, //PROP_ACT_1_GAY_M 14
		{ 0.0971, 0.1583, 0.0156 }, //PROP_ACT_3_GAY_M 15
		//GAY females same as STRAIGHT females
		{ 0.1, 0.42, 0.08163 }, //PROP_ACT_0_GAY_F 16
		{ 0.047, 0.07661, 0.00755}, //PROP_ACT_2_GAY_F 17
		{ 0.00290753, 0.00473927, 0.00046728 }, //PROP_ACT_3_GAY_F 18
		//BI males same as GAY males
		{ 0.192, 0.2652, 0.01867 }, //PROP_ACT_1_BI_M 19
		{ 0.0971, 0.1583, 0.0156 }, //PROP_ACT_3_BI_M 20
		//BI females same as STRAIGHT females
		{ 0.1, 0.42, 0.08163 }, //PROP_ACT_0_BI_F 21
		{ 0.047, 0.07661, 0.00755 }, //PROP_ACT_2_BI_F 22
		{ 0.00290753, 0.00473927, 0.00046728 }, //PROP_ACT_3_BI_F 23

		{ 15, 45, 7.65 }, //DUR_ACT_1 24
		{ 0.5, 1.5, 0.255 }, //DUR_ACT_1 25
		{ 0.25, 0.75, 0.12755 },//DUR_ACT_1 26
		{ 1.5, 4.5, 0.765 },//MEDIAN_CONC_PART_ACT_2 27
		{ 5, 15, 2.551 },  //MEDIAN_CONC_PART_ACT_3 28
		{ 0.042, 0.126, 0.0214 },  //FSW_CONDOM_MULTIPLIER 29
		{ 0.05, 0.5, 0.275 }, //EPS 30

		{ 294, 994, 178.57 }, //CD4MEAN_F 31
		{ 65, 585, 132.65 },//CD4SD_F 32
		{ 4, 5, 0.255 }, //VLMEAN_F 33
		{ 294, 994, 178.57 }, //CD4MEAN_M 34
		{ 65, 585, 132.65 },//CD4SD_M 35
		{ 4, 5, 0.255 }, //VLMEAN_M 36

		{ 0.74, 0.84, 0.0255 }, //UTILITY_CD4_LESS_THAN_50 37
		{ 0.80, 0.90, 0.0255 }, //UTILITY_CD4_LESS_THAN_50 38
		{ 0.89, 0.99, 0.0255 },//UTILITY_CD4_ABOVE_200 39
		{ -0.0265, -0.0795, -0.0135 }, //DELTA_UTILITY_WITH_HAART 40

		{ 0.000405, 0.001215, 0.00020663 }, //alpha_raw[F][M] 41
		{ 0.00021, 0.00063, 0.00010714 }, //alpha_raw[m][f] 42
		{ 0.0018, 0.0054, 0.0009184 }, //alpha_IDU_raw 43
		{ 0.000834, 0.002503, 0.0004257 }, //alpha_raw[M][M] 44

		{ 54, 150, 24.5 }, //SHARED_INJECTION_PER_YEAR 45
		{ 2.5, 7.5, 1.276 }, //IDU_PARTNER_CHANGE_RATE 46
		{ 0.125, 0.375, 0.0638 }, //DRUG_ALPHA_SCALE 47

		{ 0.73, 0.815, 0 }, //proportion of condom nonuse  48
		{ 1.063, 1.113, 0.0117 },  //RR of condon nonuse female vs. male 49
		{ 0.602, 0.658, 0.02297 }, //RR of condom nonuse gay vs. straight 50
		{ 0.7, 0.834, 0.0448 }, //RR of condom nonuse risk 2 vs. risk 1, 51
		{ 1.086, 1.188, 0.0228 }, //RR of condom nonuse risk3 vs. risk 1,  52
		{ 1, 1.58, 0.1167 }, //RR of condom nonuse alcohol/drug use vs. not 53
		{ 0.549, 0.697, 0.061 }, //RR of condom nonuse IDU 1 vs. IDU 0 54
		{ 0.4, 0.54, 0.0766 }, //RR of condom nonuse HIV inf detected vs HIV- 55
		{ 0.4, 0.54, 0.0766 }, //RR of condom nonuse HIV inf on care vs HIV- 56
		{ 0.4, 0.54, 0.0766 }, //RR of condom nonuse HIV inf treated vs HIV- 57

		{ 0.96, 0.99, 0.0078 }, //proportion of not being tested for HIV in the general population 58
		{ 0.276, 0.349, 0.0598  }, //RR of not being tested for HIV gay vs. straight 59
		{ 0.276, 0.349, 0.0598 }, //RR of not being tested for HIV bi vs. straight 60
		{ 0.793, 0.819, 0.00814 }, //RR of not being tested for HIV risk 2 vs. risk 1 61
		{ 0.221, 0.251, 0.0318 }, //RR of not being tested for HIV risk3 vs. risk 1 62
		{ 0.539, 0.666, 0.05399 }, //RR of not being tested for HIV IDU 1 vs. IDU 0 63

		{ 0.26, 0.364, 0.0859 }, //proportion of ART nonadherence in the general population 64
		{ 1.166, 3.499, 0.2803 }, //RR of ART nonadherence alcohol/drug use vs. not 65
		{ 1, 3, 0.2803 }, //RR of ART nonadherence IDU 1 vs. IDU 0 66

		{ 0.06, 0.096, 0.115 }, //proportion of untreated STI in the general population 67
		{ 0.5, 1.5, 0.28 }, //RR of untreated STI female vs. male 68
		{ 7.377, 10.328, 0.0858 }, //RR of untreated STI risk3 vs. risk 1 69
		{ 1.40, 2.046, 0.0968 }, //RR of untreated STI alcohol/drug use vs. not 70
		{ 1.223, 1.630, 0.0732 }, //RR of untreated STI IDU 1 vs. IDU 0 71

		{ 0.683, 0.917, 0.0752 }, //proportion of not being circumcised in the general population 72
		{ 1.163, 6.667, 0.4455 }, //RR of not being circumcised HIV inf detected vs HIV- 73
		{ 1.163, 6.667, 0.4455 }, //RR of not being circumcised HIV inf on care vs HIV-, 74
		{ 1.163, 6.667, 0.4455 },//RR of not being circumcised HIV inf unknown vs HIV- 75
		{ 1.163, 6.667, 0.4455 }, //RR of not being circumcised HIV inf treated vs HIV- 76

		{ 0.00059322, 0.00096694, 0.00009534 }, //prop gay women 77
		{ 0.000188, 0.000232, 0.0000112 }, //prop bi women 78

		{ 0, 0, 0 }, //motherToChild[0] 79
		{ 0.074, 0.222, 0.0378 }, //motherToChild[1] 80
		{ 0.169, 0.507, 0.0862 }, //motherToChild[2] 81
		{ 0.302, 0.906, 0.154 }, //motherToChild[3] 82
		{ 0.3625, 1.0875, 0.1849 } //motherToChild[4] 83



};
#endif 





