#include "constants.h"
#include "output.h"


int patRegCombo;	//Patient may return to a previous reg, if so, we need to know 
					//if it was reg 0, reg 1, reg 2, etc. for costing purposes
					//This value refers to the reg number in the list of regs in ChooseRegimen_Poor()



double totalCD4=0;
double meanCD4;
bool lastRegFailed;

bool firstCheckUnderRegDone; //if the first check of CD4 after being under current reg is done.

//used for DEBUG_AIDS
int totalAIDSOrDeath;

//used for DEBUG_NUM_MUTS
int totalResDrugsAtDeath;

//variable for LOOP analysis
int removeNonHIVMort = 0;		//default to NOT removing it!

//variables for data for MARK_ROBERTS_AHRQ
int numVLsup;			//number of patients with viral load supressed
double tot_CD4elev;		//total of all patients CD4 elevation from baseline
int numAIDS;			//number of patients with AIDS
int tot_regs;		//total regs over all pats
int tot_treatFail;	//total patients with treatment failure of initial reg
int tot_alive;		//total alive at certain time
int numVLFail;		//num patients who experience VL failure
int numVLFailWithMuts;	//the number of patients with virological failure and one or more res muts
bool VLFail;
int tot_intolFirstReg;	//total patients who changed first reg due to intol

int TotalPatsRes1st1years, TotalPatsRes1st5years;

int trig = TRIG;											//The type of regimen change trigger to use (values in hiv.h)

//for VRT
double unifInit[MAX_SEED_INIT];
double unifSim[MAX_SEED_CYCLE];
long int MAX_PAT_RNS;
int pat_rns_used;
int kludge_counter;

//for VL/CD4 trajectories
bool regChangeFailFlag;
bool regChangeIntolFlag;

int cycles_since_PLATO_start;
double intol_rate;

int regChangeIntol[3];										//Counts the number of reg changes due to intolerance
int regChangeFail[3];											//Coun the number of reg changes due to reg failure

//bool regFailedLastCycle, patIntolLastCycle;					//Used to provide more info to the trajectory graphs

string headerString;

//for debugging the number of mutations with DEBUG_NUM_MUTS
long int total_muts_all_pats, total_res_muts_all_pats;

int numEachTypeinReg[ARV_TYPES];

double rVLreal;
double rCD4within_delta_real;
int SurvivalPlot[3][TIME * (int)CYCLES_PER_YEAR];				//3D array: [0] is for hiv deaths, [1] is age deaths, [2] is censored
int hivdeaths[TIME];
int agedeaths[TIME];
double ChangeInCD4_1_3_5_YearsIntoTherapy[3][PATIENTS];				
double ChangeInVL_1_3_5_YearsIntoTherapy[3][PATIENTS];				
double CD4_1_3_5_YearsIntoSalvage[3][PATIENTS];			
double VL_1_3_5_YearsIntoSalvage[3][PATIENTS];	
double CD4atStartOfSalvage[PATIENTS];
double VLatStartOfSalvage[PATIENTS];
double cyclesToCompleteReg[3][PATIENTS];					//Kim added.  Used in calculating time to treatment failure for regs 1-3.
double timeTillDeath[PATIENTS];
int totalNumRegs1st5years, totalNumRegs1st10years;
int totalNumMuts1st5years, totalNumMuts1st10years;
int numDeaths;
int numPeopleCompletedReg[3];								//Added by Kim
int numDeadYear[22];
int numDeadYearHaart[22];
//TEST
extern long int initial_seed = INITIAL_SEED;
long int seed;
#ifdef COMPETING_RISKS_2
long int seed2;
#endif

int num_exhausted_clean_regs;								//Used in TIME_RES_DRUGS_BEFORE_SALVAGE
#ifdef CD4CHANGE
double YearlyCD4[600][PATIENTS];							//Used in YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY
int numAliveAtYearIntoHaart[600];							//Used in YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY AND MONITORING_PAPER
#else
double YearlyCD4[TIME][PATIENTS];							//Used in YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY
int numAliveAtYearIntoHaart[TIME];							//Used in YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY AND MONITORING_PAPER
#endif
double VLatYear5IntoHaart[PATIENTS];						//Used in MONITORING_PAPER
double VLatYear10IntoHaart[PATIENTS];						//Used in MONITORING_PAPER
int year;													//Used in YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY
int years;													//Used in CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_THERAPY (and SALVAGE)
int numStartedSalvage;										//Used in CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE
int numReached2DrugsResOnSalvage;							//Used in CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE
int totalNumRegs5yearsIntoHaart_ofPatsAlive5YrsIn;							//Used in MONITORING_PAPER
int totalNumRegs10yearsIntoHaart_ofPatsAlive10YrsIn;						//Used in MONITORING_PAPER
int totalNumMuts5yearsIntoHaart_ofPatsAlive5YrsIn;							//Used in MONITORING_PAPER
int totalNumMuts10yearsIntoHaart_ofPatsAlive10YrsIn;							//Used in MONITORING_PAPER

int total_hiv_deaths, total_age_deaths, total_censored;
double total_time_on_haart_before_salvage;					//Making this variable a double because at 1M runs, the value it holds exceeds the space allocated for an integer
int total_num_regs_before_salvage, total_res_drugs_before_salvage;
int total_res_after_initial_muts[ARV_TYPES];
int num_mort_entries;
int num_male_entries;
int num_female_entries;
int patnum;

double maxMonths;
int vl_monitor_interval;
int cd4_monitor_interval;
double CD4_previous, delta;
int eventIndicator;

double CD4_plato;			//used in calculating CD4

int offset1, offset2;

double hiv_mort_table[TABLE_SPACE][2];
#ifdef COMPETING_RISKS_1
double male_age_table[NUM_MORT_CONDITIONS + 1][TABLE_SPACE][2];
double female_age_table[NUM_MORT_CONDITIONS + 1][TABLE_SPACE][2];
int num_death_condition[NUM_MORT_CONDITIONS + 1];//number of death of each condition
#else
double male_age_table[TABLE_SPACE][2];
double female_age_table[TABLE_SPACE][2];
#ifdef COMPETING_RISKS_2
double male_CR_table[NUM_MORT_CONDITIONS][TABLE_SPACE][2];
double female_CR_table[NUM_MORT_CONDITIONS][TABLE_SPACE][2];
int num_death_condition[NUM_MORT_CONDITIONS + 1];//number of death of each condition
#endif
#endif
double vacs_hiv_neg_mort[NUM_YEARS_VACS_HIV_NEG];									//This is the mortality calculated from VACS hiv-negatives
double total_surv_time;
double mean_surv_time;
double mean_reg_time[3], median_reg_time[3];				//used in TIME_TO_TREATMENT_FAILURE
double total_completed_reg_cycles[3];						//used in KAPLAN_MEIER_CURVE_OF_TIME_TO_TREATMENT_FAILURE_REG_1_2_3

double Max(double x, double y) { return x > y ? x: y; }
double Min(double x, double y) { return x < y ? x: y; }

double start_age = AVG_AGE; 
double start_cd4 = AVG_CD4;
double start_hiv = AVG_HIV;
double cd4_treat = CD4_TREAT;

double avgAge_SD = AVG_AGE_SD;
double avgHIV_SD = AVG_HIV_SD;
double avgCD4_SD = AVG_CD4_SD;

double mut_rate_pi;
double mut_rate_nrti;
double mut_rate_nnrti;

double stop_reg = STOP_REG;;
double mortMultHIV = MORTMULTHIV;
double intolMultiplier = INTOL_MULTIPLIER;

int mut_pi_singular_start = MUT_PI_SINGULAR_START;
int mut_pi_boosted_start = MUT_PI_BOOSTED_START;
int mut_nrti_tam_start = MUT_NRTI_TAM_START;
int mut_nrti_nontam_start = MUT_NRTI_NONTAM_START;
int mut_nnrti_efavirenz_start = MUT_NNRTI_EFAVIRENZ_START;
int mut_nnrti_nevirapine_start = MUT_NNRTI_NEVIRAPINE_START;

double pmutres_pi_singular_est = PMUTRES_PI_SINGULAR_EST;
double pmutres_pi_boosted_est = PMUTRES_PI_BOOSTED_EST;
double pmutres_nrti_tam_est = PMUTRES_NRTI_TAM_EST;
double pmutres_nrti_nontam_est = PMUTRES_NRTI_NONTAM_EST;
double pmutres_nnrti_efavirenz_est = PMUTRES_NNRTI_EFAVIRENZ_EST;
double pmutres_nnrti_nevirapine_est = PMUTRES_NNRTI_NEVIRAPINE_EST;

double pcrossres_pi_singular_est = PCROSSRES_PI_SINGULAR_EST;
double pcrossres_pi_boosted_est	= PCROSSRES_PI_BOOSTED_EST;
double pcrossres_nrti_tam_est = PCROSSRES_NRTI_TAM_EST;
double pcrossres_nrti_nontam_est = PCROSSRES_NRTI_NONTAM_EST;
double pcrossres_nnrti_efavirenz_est = PCROSSRES_NNRTI_EFAVIRENZ_EST;
double pcrossres_nnrti_nevirapine_est = PCROSSRES_NNRTI_NEVIRAPINE_EST;

double pcrossres_piboosted_singular = PCROSSRESPI_BOOSTED_SINGULAR;  
double pcrossres_pisingular_boosted = PCROSSRESPI_SINGULAR_BOOSTED;  
double pcrossres_tam_nontam = PCROSSRESNRTI_TAM_NONTAM;  
double pcrossres_nontam_tam = PCROSSRESNRTI_NONTAM_TAM;  
double pcrossres_efav_nevir = PCROSSRESNNRTI_EFAV_NEVIR;  
double pcrossres_nevir_efav = PCROSSRESNNRTI_NEVIR_EFAV; 

double comp = COMP;
double mut_rate = MUT_RATE;
double factor = FACTOR;
double adjmutres = ADJMUTRES;
double adjcrossres = ADJCROSSRES;

double start_with_aids = START_WITH_AIDS;
double AIDSrate_less200, AIDSrate_greater200;

double haart_tox = HAART_TOX;		//Kim added.

double testCostVL = TESTCOSTVL;
double testCostCD4 = TESTCOSTCD4;
double cost_of_care_poor = COST_OF_CARE_POOR;        //LL added
//double hospital_cost_poor = HOSPITAL_COST_POOR;      //LL added
double cost_reg1_poor = COST_REG1_POOR;                  //LL added
double cost_reg2_poor = COST_REG2_POOR;                  //LL added
double cost_reg3_poor = COST_REG3_POOR;

int monitor_interval = MONITOR_INTERVAL;
double vl_regimen_failure_1st = VL_REGIMEN_FAILURE_1ST;					//for MONITORING_PAPER
double vl_regimen_failure_after_1st = VL_REGIMEN_FAILURE_AFTER_1ST;		//for MONITORING_PAPER
int vl_thresh;															//for MONITORING_PAPER loop
int numAvailableRegimens = NUM_CLASS_COMBINATIONS;						//only used for resource poor for MONITORING_PAPER
int poorCanGoBackOnPreviousRegs = POOR_CAN_GO_BACK_ON_PREVIOUS_REGIMENS;  
int patient_starts_on_haart = PATIENT_STARTS_ON_HAART; // Indicates whether patient will be initialized to be on HAART when patient is created.

double qa_time, mean_qa_surv_time, total_qa_surv_time, mean_qa_disc_surv_time, total_qa_disc_surv_time, total_disc_surv_time, mean_disc_surv_time;
double cost, total_cost, mean_cost, total_disc_cost, mean_disc_cost,individual_cost[PATIENTS];
double drug_cost, total_drug_cost, mean_drug_cost, total_drug_disc_cost, mean_disc_drug_cost;
double total_lab_cost, mean_lab_cost, total_lab_disc_cost, mean_disc_lab_cost;
double hospital_cost, total_hospital_cost, mean_hospital_cost, total_hospital_disc_cost, mean_disc_hospital_cost;
double care_cost, total_care_cost, mean_care_cost, total_care_disc_cost, mean_disc_care_cost;
double alc_cost, total_alc_cost, mean_alc_cost, total_alc_disc_cost, mean_disc_alc_cost; //costs related to alc intervention


int numAlive_1_3_5_YearsIntoHaart[3];  //Used in CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_THERAPY 
int numAlive_1_3_5_YearsIntoSalvage[3];  //Used in CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE	




//double Util_dec;
double qamedian[PATIENTS];

FILE *output, *file1, *kim_test, *kim_test2, *best_fit, *stateprob;
FILE *who_monitoring; //used for WHO_MONITORING_STRATEGY_ANALYSIS
ofstream senseFile;
string label;

/*
class changes {
	public:
	int state[5][5][5][5]3][3][4][4][3][3];
	int rowcd4, rowvl, rowtrt, rowres, rowlfe;
	
	void movto(int cd4, int vl, int trt, int res, int lfe)
	{
		state[rowcd4][cd4][rowvl][vl][rowtrt][trt][rowres][res][rowlfe][lfe]++;
	}
	
	
};*/

//TEST
typedef struct 
{
	int drugNum;
	int classNum;
	bool res;
	bool intol;
}drug;

drug drugs[MAX_ARV_TYPES][MAX_DRUGS_IN_CLASS];		//the array of all drugs

typedef struct {
	double start_age;
	double CD4baseline;
	double HIVbaseline;

    int cyclesOnReg[MAX_REGS_MODEL];    // the number of cycles on a given therapy 0 = haart1, 1 = haart2, 2 = haart3, etc.

	int num_regs_before_salvage;		//used for REG_INFO_FOR_SENSITIVITY
	int num_res_drugs_before_salvage;	//used for REG_INFO_FOR_SENSITIVITY

	int arv_types;					//the number of different drug classes

	int mut_arv[MAX_ARV_TYPES];		//gives number of mutations against each arv type
	int num_in_class[MAX_ARV_TYPES];	//a counter for the number of drugs in each drug class
	int num_res[MAX_ARV_TYPES];		//a counter for the number of drugs in each class the patient is resistant to
	double pmutres[MAX_ARV_TYPES];
	double pcrossres[MAX_ARV_TYPES];
	int crossClassType[MAX_ARV_TYPES];
	double pCrossClassRes[MAX_ARV_TYPES];
	double mutRateClass[MAX_ARV_TYPES];	//Given the overall mutation rate, this is the probability of a mutation occurring in this class
	char className[MAX_ARV_TYPES][20];

	drug *reg[3];					//The patient's regimen.  The pointers point to the drugs in the drug array
	int init_reg[3];				//holds the drug types of the initial regimen

	int nc[DRUGS_IN_REG];		//indicates if pat is compliant to drug 1, 2, 3

	bool haart;					//indicates if pat is on haart or not
	bool started_haart;
	int quit_haart;
    int CD4cat;					
	int VLcat;

	int totalresnc;				//total number of drugs being taken that virus is resistant to or pat is not compliant to.  
	int totalres;				//total number of drugs in reg that virus is resistant to
	int totalintol;				//total number of drugs in reg that pat is intol to
	int numreg;
	bool exhausted_clean_regs;	//indicates if patient is resistant to all possible regimens

	int cyclenum;				//I use this to mean number of times through "Cycle" branch, however
								//DM increases this anytime it goes through a branch such as OnePI, DR, etc.
	int counterreg;				//gives number of cycles on this reg
	int haart_start_time; 
	int haart_exhausted_clean_regs_time;
	double vl_regimen_failure;	//the VL threshold for regimen failure
	bool has_AIDS;				//bool indicating whether an AIDS defining event has ever occurred
	bool AIDSeventReg;			//bool indicating whether an AIDS defining event has occured during the current regimen


	double CD4real;				//gives the real CD4 level that accounts for non-inst change to the cd4 ideal.
								//set up so patient's CD4 moves a fraction of the way to the CD4 ideal
	double CD4real_before_noise;//the smooth CD4real - we add noise to make this look more natural
	double CD4offreal;
	double CD4regressionideal;	//the CD4 as calculated by Joyce's regression equation
	double CD4regressionreal;

	double CD4weightedForComp;

	double CD4Mellors;			//the CD4 off haart as calculated by the Mellors equation
	double VLreal;
	double VLreal_before_noise;				
	double VLideal;
	double VLdec;

	bool historyOfIntolerance;		//1 if patient has ever been intolerant, else 0
	bool hasResistance;				//says if pat has developed resistance to any drugs whatsoever
	
	int nextMonCycle;				//Number of cycles from now to check VL or CD4 again

	//values at start of each regimen - used for CD4within
	double AgeRegBaseline;
	double CD4RegBaseline;
	double CD4atStartOfHaart;
	double VLatStartOfHaart;

	double censorRate;			//used only for calibrating resource poor
	double pCensor;				//used only for calibrating resource poor
	double deathrateHIV;
	double deathrate_nonHIV;
	double pDieHIV;
	double pDieAge;
	double pAIDS;
	double age;
	double pNocomp;
	double pMutate;
	double adjmutrate;
	double vl_fract;
	double luck1;
	double luck2;
	double diligence;
	double comp;
#if ALC_INTERVENTION
	bool alc;  //is alc or not?
#endif
#if BOTH_AD_ALC
	bool alc;  //is alc or not?
#endif
#if DETECT_VL
	bool interv_on; //on intervention?
#endif

	double pstopreg;
	double maxCD4onReg;						//used for cd4 reg change trigger
	int mutCount;				//The number of mutations a patient has, regardless of resistance
	bool reached2DrugsResOnSalvage;								//"  "  Says whether patient was ever stuck with two resistant drugs on salvage (no better regimen and on Plato for the long slide down)
	double CD4Peak;				//The maximum CD4 ever reached on ART

	bool lastLess100; //on last check, CD4 is less than 100

	//LL: following are added for Retention In Care 
#if RETENTION_IN_CARE
	bool engaged_in_care;
	bool engaged_in_clinic;
	bool LTFU;						//binary: is this pat Lost-Follow-Up? yes/no
	int num_LTFU;					//how many times LTFU
	double num_cycles_ES1;
	double num_cycles_ES3;          //how many cycles in state E3
	double num_cycles_ES4;          //how many cycles in state E4
	double num_cycles_ART;          //how many cycles on ART;
	double cyclenum_LTFU;				//time point patient shows up in clinic last time
	double time_in_care;      // time_in_care= start of model+cycle time; should resert to 0 if patient becomes disengaged from care (ES4)
	
	//for OUTREACH_INTERVENTION
	bool found;               // case 1: no outreach intervention: back by himeself; Case 2: with outreach intervention: Have we found the pat? Note: stopoutreach != found
	bool outreach;            // Is the pat in outreach process?
	bool deathUnkown;  // Did the patient die during outreach?

	bool willInitiateHAARTAtCD4TreatThresh;		//if TRUE, pat starts on HAAART at CD4_TREAT, else will start HAART at cd4=200
#endif


} patient;

typedef patient *pat_ptr;




struct tttfData
{
	int failed;
	int censored;
	int survived;
};

tttfData tttfDat[TIME*(int)CYCLES_PER_YEAR][3];		


// used for generating transition model lookup table
int statecount[STCD4][STVL][STTREAT][STRESIST][STCD4][STVL][STTREAT][STRESIST][STDEAD]; // cd4, viral load, treatment, resistant, live/deadAids  (we are only recording deaths due to HIV)
int cur_state[5], new_state[5], cost_state[5];		//Each "state" is made up of the five following categories {cd4, vl, treat yes/no, res, alive/dead of hiv}														// holds cd4, vl, treat, res and alive/dead	 
double totalAnnualCost[STCD4][STVL][STTREAT][STRESIST][STCD4][STVL][STTREAT][STRESIST][STDEAD];
double totalCostAtStartOfYearThisPat;
double costThisYearThisPat;

#ifdef MONITORING_PAPER_SENSITIVITY
double sens[5][2];
int runnum = RUN_NUM;
#endif

//used for WHO_MONITORING_STRATEGY_ANALYSIS

//We want to gather stats at 5, 10, 20 years
#define WHO_YEARLY_OUTPUT_NUM_YEARS 20
#define WHO_5_YR_INDEX 20
#define WHO_10_YR_INDEX 21
#define WHO_20_YR_INDEX 22
#define WHO_NUM_TIME_PERIODS 23

int scenario;
int numPatVisits[WHO_NUM_TIME_PERIODS];
int numCD4Tests[WHO_NUM_TIME_PERIODS];
int numVLTests[WHO_NUM_TIME_PERIODS];
int patCyclesOn1stLineTherapy[WHO_NUM_TIME_PERIODS];
int patCyclesOn2ndLineTherapy[WHO_NUM_TIME_PERIODS];
int totalDeaths[WHO_NUM_TIME_PERIODS];
int totalHIVDeaths[WHO_NUM_TIME_PERIODS];
int lifeCyclesLived[WHO_NUM_TIME_PERIODS];
int lifeCyclesOnSuccessfulART[WHO_NUM_TIME_PERIODS];	//life years before occurrence of AIDS defining event
int lifeCyclesWithAIDS[WHO_NUM_TIME_PERIODS];
int numEverDevelopedAIDS[WHO_NUM_TIME_PERIODS];
int sixMonthsSinceLastVisit;
//end of WHO_MONITORING_STRATEGY_ANALYSIS


//used for RYAN_WHITE_ANALYSIS
#ifdef RYAN_WHITE_ANALYSIS
#define RYAN_WHITE 1
#else
#define RYAN_WHITE 0
#endif
#define NUM_ANALYSES 8	//First analysis is a regular run with no changes from baseline
#define TOTAL_PATS_RYAN_WHITE 17535
int costPerPersonPerYear[NUM_ANALYSES] = {0, 2181, 1279, 6825, 1276, 4426, 2212, 3057};
int numTargetedInterv[NUM_ANALYSES] = {TOTAL_PATS_RYAN_WHITE, 2701, 3033, 1556, 3581, 5371, 2141, 986};
double adh_RR[NUM_ANALYSES] = {1, /*Old 1.19*/ 1.28, /*Old 1.12*/1.5, 1.093, 1, /*Old 1*/1.5, 1, 1};
double linkage_RR[NUM_ANALYSES] = {1, 2.0, /*1.28*/1.7, 1.124, 1.21, 1.23, 1.275, 1.25};
bool patIsTargeted;
int interv;
double baselineComp = .63;
double baselineLinkage = .749;
int numTargeted;
double costRyanWhiteInterv;
double propTargeted;
int numLinked;
int numBaselineComp, totalNumCD4_50, numTargetedLinked, numNotTargetedLinked;
int numSuppressed;
int numDroppedAtLeastOneLogVL;

//Will hold costs, life years, etc. at years 2,5,10,20 and lifetime
#define NUM_TIME_PERIODS_RYAN_WHITE 4
double costAtYear[NUM_TIME_PERIODS_RYAN_WHITE];		
double cyclesLivedAtYear[NUM_TIME_PERIODS_RYAN_WHITE];

//End variables used for RYAN_WHITE_ANALYSIS

#ifdef WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
struct discordance_stages {
	//Variables are named to correspond to the columns in the spreadsheet:
	//cd4 Less than 100 viral load Greater than 5000 is cd4_L100_vl_G5000
	//cd4 Less than 50% peak viral load Greater than 5000 is cd4_L50p_base_vl_G5000
	int cd4_L100_vl_G5000, cd4_L100_vl_L5000G1000, cd4_L100_vl_L1000G500, cd4_L100_vl_L500;
	int cd4_L50p_base_vl_G5000, cd4_L50p_base_vl_L5000G1000, cd4_L50p_vl_L1000G500, cd4_L50p_vl_L500;

	discordance_stages()
	{
		cd4_L100_vl_G5000 = cd4_L100_vl_L5000G1000 = cd4_L100_vl_L1000G500 =  cd4_L100_vl_L500 = 0;
		cd4_L50p_base_vl_G5000 = cd4_L50p_base_vl_L5000G1000 = cd4_L50p_vl_L1000G500 = cd4_L50p_vl_L500 = 0;
	}
	
	void initialize()
	{
		discordance_stages();
	}
	
	void increment_category(patient *pat)
	{
		if (pat->CD4real < 100)
		{
			if (pat->VLreal > 3.7) // Greater than 5000
			{
				cd4_L100_vl_G5000++;
			} else if(pat->VLreal > 3) // Greater than 1000
			{
				cd4_L100_vl_L5000G1000++;
			} else if(pat->VLreal > 2.7) // Greater than 500
			{
				cd4_L100_vl_L1000G500++;
			} else // Less than or equal to 500
			{
				cd4_L100_vl_L500++;
			}
		}
		
		if (pat->CD4real < pat->CD4baseline || pat->CD4real < pat->CD4Peak/2.0)
		{
			if (pat->VLreal > 3.7) // Greater than 5000
			{
				cd4_L50p_base_vl_G5000++;
			} else if(pat->VLreal > 3) // Greater than 1000
			{
				cd4_L50p_base_vl_L5000G1000++;
			} else if(pat->VLreal > 2.7) // Greater than 500
			{
				cd4_L50p_vl_L1000G500++;
			} else // Less than or equal to 500
			{
				cd4_L50p_vl_L500++;
			}
		}
	}
	
};

discordance_stages cumulative_years[2][4]; // has_AIDS, time
discordance_stages alt_cumulative_years[2][4];
#endif

#if RETENTION_IN_CARE
int eventIndicator_RIC;
bool LTFUflag;
bool ES1flag, ES3flag, ES4flag;		// flags are for graphing trajectories - indicates an event happened
int total_ES1_HIV_death_known;
int total_ES1_HIV_death_unknown;
int total_ES1_AGE_death_known;
int total_ES1_AGE_death_unknown;
int total_ES3_HIV_death;
int total_ES3_AGE_death;
int total_ES4_HIV_death;
int total_ES4_AGE_death;
double total_ES1_time;          // cumulative amount of time for entire population on E1
double total_ES3_time;
double total_ES4_time;
double total_time_onART;
double total_time_VL_less1000;      // total time with VL <1000 (Progression)
double mean_ES1_time;           // mean_E1_time = total_E1_time/#patients
double mean_ES3_time;
double mean_ES4_time;
double mean_time_onART;
double mean_time_VL_less1000;       // mean time with VL <1000 (Progression)
int total_LTFU;
double mean_LTFU;
double pre_ART_mult_for_prob_es1_to_es2 = PRE_ART_MULT_FOR_PROB_ES1_TO_ES2;
double prob_es1_to_es2;  //not assign value here cus it varies based on time-in-care
double prob_es2_to_es4 = PROB_ES2_TO_ES4;
double prob_es3_to_es4;  //prob_es3_to_es4 = prob_es1_to_es2 * prob_es2_to_es4. do in hiv_toxicity.cc
double prob_es5_to_es6 = PROB_ES5_TO_ES6;
double prob_es4_to_es1_AIDSpat = PROB_ES4_TO_ES1_AIDSPAT; // pat returns to care only if he develops AIDS defining event
double mutation_mult_es4 = MUTATION_MULT_ES4;
double es1_to_es2_multiplier_for_sensitivity = 1;				//Used for RIC sensitivity (1 is no effect)

int RIC_command_line_sensitivity_row_num;	//used for running the RIC sensitivity on the cluster

//for OUTREACH_INTERVENTION
bool outreach_interv_enabled = OUTREACH_INTERVENTION;			//on or off
bool LTFUfoundflag;
bool stopoutreachflag;
bool startoutreachflag;
//double base_multiplier_for_outreach = BASE_MULTIPLIER_FOR_OUTREACH;
double outreach_prob_identify = OUTREACH_PROB_IDENTIFY;       //Id patient
double outreach_prob_find = OUTREACH_PROB_FIND;               //Finding/reaching
double outreach_prob_relink = OUTREACH_PROB_RELINK;           //Relinking
int outreach_start_trigger = OUTREACH_START_TRIGGER;
int outreach_end_trigger = OUTREACH_END_TRIGGER;
double outreach_initiation_cost_poor = OUTREACH_INITIATION_COST_POOR;
double outreach_finding_cost_poor = OUTREACH_FINDING_COST_POOR;    //
double outreach_finishing_cost_poor = OUTREACH_FINISHING_COST_POOR;
double total_outreach_cost, total_disc_outreach_cost, /*total_outreach_finding_cost,*/ mean_outreach_cost;
int total_num_outreach, total_num_stopoutreach;
double mean_num_outreach, mean_num_stopoutreach; // this is to test whether we stop outreach effort correctly
double total_outreach_finding_cost;		//temporary!
double outreach_cd4 = OUTREACH_CD4;

//for RISK_REDUCTION_INTERVENTION
bool risk_reduction_interv_enabled = RISK_REDUCTION_INTERVENTION;					//on or off
double risk_reduction_rr_prob_es1_to_es2 = RISK_REDUCTION_RR_PROB_ES1_TO_ES2;
double risk_reduction_intervention_cost = RISK_REDUCTION_INTERVENTION_COST;
double total_RR_intervention_cost, mean_RR_intervention_cost;
double risk_reduction_cd4 = RISK_REDUCTION_CD4;

//for SECONDARY_PREVENTION_INTERVENTION
bool secondary_prevention_interv_enabled = SECONDARY_PREVENTION_INTERVENTION;
double secondary_prevention_rr_prob_es1_to_es2 = SECONDARY_PREVENTION_RR_PROB_ES1_TO_ES2;
double secondary_prevention_intervention_cost = SECONDARY_PREVENTION_INTERVENTION_COST;
double total_SP_intervention_cost, mean_SP_intervention_cost;
double secondary_prevention_cd4 = SECONDARY_PREVENTION_CD4;

double prob_ART_init_at_CD4_treat = PROB_ART_INITIATION_AT_CD4_TREAT;
#endif
