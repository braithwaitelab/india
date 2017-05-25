#define STOPTIME 20//for calibration
#define CALIBRATE 0 							//run calibration for 14 years to generate new population file

/*******************************Specific Analyses********************************/
//#define RUN_ALL_NO_SCRIPT						//runs each intervention, one at a time
//#define TESTING_LINKAGE_ADHERENCE

#define RUN_ALL_USING_SCRIPT 0					//runs each intervention, one at a time
#define OPTIMIZING_PREVENTION_PORTFOLIO 0		//runs all possible combinations of interventions
#define COST_SENSITIVITY 0
#define INTERVENTION_EFFECT_SIZE_SENSITIVITY 0	//varies effect sizes within plausible range
#define ONE_WAY_SENSITIVITY_ANALYSIS 0
#define ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION 0		//Isolates effect of condom use, stis, and adherence
#define OPTIMIZING_ALCOHOL_INTERVENTION 0				//Variables that import values from run script are toward the bottom of this page
#define ALCOHOL_SURFACE 0
#define CARE_PLUS 0
#define DELPHI 0 //use delphi results when equals 1 and use literature when equals 0 
#define COST_EFFECTIVENESS_SENSITIVITY 0 //probabilistic analysis of different basket combinations
#define PROBABILISTIC_RUN	1
/*******************************End Specific Analyses********************************/
#define SHOW_EFFECT_SIZES 0					//Creates extra output with effect sizes
#define PRINT_CYCLE_BY_CYCLE_OUTPUT 0
#define DEBUGG 0 
#define DEBUG_TRANSITIONS 0
#define PRINT_BALANCE_ERROR 1
#define OUTPUT_BY_ACT_GROUP 0
#define SAVE_POPULATION 0
#define PRINT_NUM_IN_EACH_STRATA_EVERY_YEAR 0
#define DEBUG_SENSITIVITY 0

//Note that these will get shut off automatically during calibration in code below
#define EFFICIENCY_PRECALCULATE_ALPHA_MULT_VALUES 1
#define EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS 1

//The following are just to make test scenarios easier0
#define TURN_ON_DEMOGRAPHICS 1
#define TURN_OFF_HIV_MORTALITY 0
#define TURN_OFF_TRANSMISSION_ON_TREATMENT 0
#define NEW_INFECTIONS_GO_IMMEDIATELY_ON_TREAT 0
#define TURN_OFF_SEXUAL_TRANSMISSION 0  
#define TURN_OFF_ALL_TRANSMISSION 0							//Turns off sexual and mother to child transmission
#define TURN_OFF_BIRTHS 0

#define DT .02



//The following provide the number of possible values of each category
#define NUM_STATUS 8					//the number of states 
#define NUM_ALIVE_STATE 6				//the number of states in which people are alive
#define NUM_ACT	4						//the number of activity groups (defined by durations, below)
#define NUM_AGE	10						//the number of 5-year age groups, 0 represents 0-4, 1 is 5-9 etc.
#define NUM_SEX	2						//the number of genders (M,F)
#define NUM_ORIENT 3					//the number of orientations (Straight, Gay, Bi)
#define NUM_VL 5						//the number of VL categories
#define NUM_CD4 5						//the number of CD4 categories
#define NUM_RES 2						//the number of resistance categories
#define NUM_ALC 4						//the number of alcohol categories, where "alcohol" is used to represent a person with alcohol abuse, mental health issues, or drug abuse.  '0' is not having it, '1' is having it
#define NUM_IDU 2						//the number of Injection Drug User categories.  0 is non-user, 1 is user  (Set NUM_IDU to 1 if model doesn't need IDU users.  Set to 2 if model does need to incorporate IDU users.)
#define INDEX_MAX_ACTIVE_STATE TREAT	//the array index of the last sexually active state
#define ADULTHOOD 3						//the index of the first age group containing sexually active adults
#define FIRST_ACTIVE_CLASS 1			//the index of the first non-abstinent activity group (so if you want an abstinent group, make it activity group 0!)
#define AGE_CHANGE_RATE (1.0/5.0)		//the rate at which people move in and out of the 5-year age groups

/*******Currently AGE_SEXUAL_DEBUT does not work as intended.  Needs to be fixed in subsequent version!! *******/
#define AGE_SEXUAL_DEBUT 19				//the age at which a person becomes sexually active

#define EPS	0.2 						//degree of assortative mixing (0=assortative, 1=proportionate)	

#define COMPROMISE MINIMUM //CT MEAN	//Determines how two groups will compromise during the balancing

#define PROP_GAY_M 0.0006
#define PROP_GAY_F 0.0006
#define PROP_BI_M 0.0004
#define PROP_BI_F 0.0004
#define PROP_STRAIGHT_M (1 - PROP_GAY_M - PROP_BI_M) 
#define PROP_STRAIGHT_F (1 - PROP_GAY_F - PROP_BI_F) 

#define PROP_ALCOHOL_M .185				//The initial proportion of men with alcohol abuse problems
#define PROP_ALCOHOL_F .02				//The initial proportion of women with alcohol abuse problems
#define PROP_ALC_DEPRESSED 0.57			//The proportion of ALCS who are depressed 

#if (NUM_IDU == 2)
//the initial proportion who are risky injection drug users (meaning they used shared or unclean needles)
#define PROP_IDU 0.0000499211/*proportion of IDU in population*/*0.26/*number of IDUs consistently sharing needles*/ 
#else 
#define PROP_IDU 0						//If the model doesn't need IDU (thus we've set NUM_IDU to 1), set the PROP_IDU to 0.
#endif

#define PROP_DETECTED 0				//the proportion of infecteds who are detected at the start of the calibration period in 1997

#define PROP_INFECTED_INFANTS_START_IN_TREAT	0.0		//Proportion of HIV infected infantas immediately started in "TREAT"
#define MORT_MULT_CHILDREN_NOT_ON_TREAT			6.0		//Mortality adjustment for <5 yr olds who are HIV infected and not in treatment
#define PROP_PREGNANT_WOMAN_ON_PMTCT_1997		0.0		//Proportion of pregnant woman on PMTCT in 1997
#define PROP_PREGNANT_WOMAN_ON_PMTCT_2014		.84		//Proportion of pregnant woman on PMTCT iat start of 2014 change to .8 when running intervention 2: Enhanced PMTCT/Immediate ART for all pregnant women
#define PROB_NOT_BEING_TESTED_FOR_HIV_START_OF_CALIBRATION 1  //Probability of not being tested for HIV in 1997
#define NUM_YEARS_NO_LINKAGE_TO_CARE_FROM_1997	6.0		//There was no linkage to care for this number of years after 1997.  Thus, this is the begining of the ramp-up period where the prob of transitioning to Care increases over time during calibration.

//Equation of line where x is concurrency and y is annual frequency
//This data yields a line where (conc=1, freq_per_partnership = 98.8) and (conc=4, freq_per_partnership = 145.6)
//Note that beyond four partnerships, we plateau and simply use the max value of 145.6
#define FREQ_VARIES_WITH_CONC 1					//Yes or No
#define SLOPE 15.6
#define Y_INTERCEPT 83.2

//Choose a version of the activity group allocation (proportions, duration, concurrency, etc.)


/***************Variables to set up activity groups*********************/

#define PROP_ACT_0_STRAIGHT_M	.2775		//Proportion of males who are abstinent
#define PROP_ACT_1_STRAIGHT_M	.5733		//Proportion of males who are in stable, monogamous relationships
#define PROP_ACT_2_STRAIGHT_M	.1292		//Proportion of males in activity group 2
#define PROP_ACT_3_STRAIGHT_M	.02			//Proportion of males in activity group 3

#define PROP_ACT_0_STRAIGHT_F	.26			//Proportion of females who are abstinent
#define PROP_ACT_1_STRAIGHT_F	.69			//Proportion of females who are in stable, monogamous relationships
#define PROP_ACT_2_STRAIGHT_F	.047		//Proportion of females in activity group 2
#define PROP_ACT_3_STRAIGHT_F	.003		//Proportion of females in activity group 3

#define PROP_ACT_0_GAY_M	0			//Proportion of gay males who are abstinent
#define PROP_ACT_1_GAY_M	.23			//Proportion of gay males who are in stable, monogamous relationships
#define PROP_ACT_2_GAY_M	.675		//Proportion of gay males in activity group 2
#define PROP_ACT_3_GAY_M	.095		//Proportion of gay males in activity group 3

#define PROP_ACT_0_BI_M		0			//Proportion of bi males who are abstinent
#define PROP_ACT_1_BI_M		.23			//Proportion of bi males who are in stable, monogamous relationships
#define PROP_ACT_2_BI_M		.675		//Proportion of bi males in activity group 2
#define PROP_ACT_3_BI_M		.095		//Proportion of bi males in activity group 3

#define PROP_ACT_0_GAY_F	.26			//Proportion of gay females who are abstinent
#define PROP_ACT_1_GAY_F	.69			//Proportion of gay females who are in stable, monogamous relationships
#define PROP_ACT_2_GAY_F	.047		//Proportion of gay females in activity group 2
#define PROP_ACT_3_GAY_F	.003		//Proportion of gay females in activity group 3

#define PROP_ACT_0_BI_F		.26			//Proportion of bi females who are abstinent
#define PROP_ACT_1_BI_F		.69			//Proportion of bi females who are in stable, monogamous relationships
#define PROP_ACT_2_BI_F		.047		//Proportion of bi females in activity group 2
#define PROP_ACT_3_BI_F		.003		//Proportion of bi females in activity group 3


#define PROP_ACT_0_IDU		.13		//Proportion of IDU who are abstinent 
#define PROP_ACT_1_IDU		.19		//Proportion of IDU who are in stable monogamous relationships
#define PROP_ACT_2_IDU		.34		//Proportion of IDU who are in activity group 2
#define PROP_ACT_3_IDU		.34		//Proportion of IDU who are in activity group 3


#define DUR_ACT_1		30.0	//Average duration of stable, monogamous partnerships
#define DUR_ACT_2		1.0		//Average duration of activity group 2
#define DUR_ACT_3		.5		//Average duration of activity group 3
#define MEDIAN_CONC_PART_ACT_1		1.0		//Median number of concurrent partnerships for activity group 1
#define MEDIAN_CONC_PART_ACT_2		3.0		//Median number of concurrent partnerships for activity group 2
#define MEDIAN_CONC_PART_ACT_3		10.0	//Median number of concurrent partnerships for activity group 3


#define FREQ	700		
/**********************************************************************/

//Calibration parameters
#define SEX_ALPHA_SCALE	3	
#define HIV_MORT_MULT_OFF_TREAT	.30
#define HIV_MORT_MULT_ON_TREAT .55
#define DRUG_ALPHA_SCALE 0.25

#define NUM_PATHWAYS 8
#define NUM_PATH_AND_CHAR 8
#define NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE 14
#define NUM_BASELINE_RR 13
#define NUM_INTERVENTIONS 26
#define CD4_THRESH_PATHWAY 12

#define MODE PESSIMISTIC//PESSIMISTIC//OPTIMISTIC

//for reading in the HIV progress rates
#define CARESTATUS 2					//two possibiltiies, on care or not
#define NUM_ALIVEORDEAD 2				//two possibilities, dead or alive
#define NUM_TREATORNOT 2				//two possibilities, on treatment or not
#define NUM_ADHERENCE 12				//the number of adherence strata
#define OFFCARE_INDEX 11				//This is the element of the adherence array used to hold the offcare rates.

#define NOT_APPLICABLE -999

#define SUSC 0							//Susceptible
#define HVL 1							//Infected, initial high VL state
#define INF 2							//Infected (not yet detected)
#define DET 3							//Infected and detected
#define CARE 4							//Infected, On Care
#define TREAT 5							//Infected, On Care, On Treatment
#define DEAD_HIV 6						//Dead of HIV
#define DEAD_AGE 7						//Dead of age

#define M 0								//Men
#define F 1								//Women

#define STRAIGHT 0						//Orientations
#define GAY 1
#define BI	2

#define INCIDENT 0						//for use in alpha array
#define PREVALENT 1
#define LATE_STAGE 2

#define PESSIMISTIC 0					//Settings for how to combine multiple interventions
#define OPTIMISTIC 1

//To make intervention pathway 12 (cd4 treatment threshold) clearer
#define CD4_TREAT_THRESH_DEFAULT 1
#define CD4_TREAT_THRESH_10000 0

//To make intervention pathway 11 (PMTCT) clearer
#define PMTCT_DEFAULT   PROP_PREGNANT_WOMAN_ON_PMTCT_2014

//To make Kim's life easier
#define LOOP_THROUGH_ALL_COMPARTMENTS  for (status=0; status<NUM_STATUS; status++){for (sex=0; sex<NUM_SEX; sex++){for (orient=0; orient<NUM_ORIENT; orient++){for (act=0; act<NUM_ACT; act++){for (age=0; age<NUM_AGE; age++){for (alc=0; alc<NUM_ALC; alc++){for (idu=0; idu<NUM_IDU; idu++){for (vl=0; vl<NUM_VL; vl++){for (cd4=0; cd4<NUM_CD4; cd4++){for (res=0; res<NUM_RES; res++){
#define LOOP_THROUGH_ALL_ALIVE_COMPARTMENTS  for (status=0; status<NUM_ALIVE_STATE; status++){for (sex=0; sex<NUM_SEX; sex++){for (orient=0; orient<NUM_ORIENT; orient++){for (act=0; act<NUM_ACT; act++){for (age=0; age<NUM_AGE; age++){for (alc=0; alc<NUM_ALC; alc++){for (idu=0; idu<NUM_IDU; idu++){for (vl=0; vl<NUM_VL; vl++){for (cd4=0; cd4<NUM_CD4; cd4++){for (res=0; res<NUM_RES; res++){
#define LOOP_THROUGH_ALL_ADULT_ALIVE_COMPARTMENTS for (status=0; status<NUM_ALIVE_STATE; status++){for (sex=0; sex<NUM_SEX; sex++){for (orient=0; orient<NUM_ORIENT; orient++){for (act=0; act<NUM_ACT; act++){for (age=ADULTHOOD; age<NUM_AGE; age++){for (alc=0; alc<NUM_ALC; alc++){for (idu=0; idu<NUM_IDU; idu++){for (vl=0; vl<NUM_VL; vl++){for (cd4=0; cd4<NUM_CD4; cd4++){for (res=0; res<NUM_RES; res++){

#define END_LOOP_THROUGH_ALL_COMPARTMENTS }}}}}}}}}}

#define LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS for (status=0; status<NUM_ALIVE_STATE; status++){for (sex=0; sex<NUM_SEX; sex++){for (orient=0; orient<NUM_ORIENT; orient++){for (act=0; act<NUM_ACT; act++){for (alc=0; alc<NUM_ALC; alc++){for (vl=0; vl<NUM_VL; vl++){for (idu=0; idu<NUM_IDU; idu++){
#define END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS }}}}}}}



//don't edit these Macros!
#define mOppSex(sexx) (sexx==M)?F:M		//returns the value of the partners sex

#define START 0							//used for indicating whether pop is being printed at start or end of given time
#define END 1

#define ONCARE(stat) ((stat>=CARE)?1:0)
#define ONTREAT(stat) ((stat>=TREAT)?1:0)

#define HVL_TO_INF  8.7		//Rate of going from HVL to INF.  Average period of initial infection is 6 weeks = 42days = .115 years.  Thus the rate is 1/.115 = 8.7 years.

//the distributions for VL and CD4 for new infecteds (different sets of values for male and female.  
#define CD4MEAN_F 644			
#define CD4SD_F	260					
#define VLMEAN_F 4.46
#define VLSD_F 0.99

#define CD4MEAN_M 644				
#define CD4SD_M 260			
#define VLMEAN_M 4.46
#define VLSD_M 0.99

//For running tests with just two groups and varying compromising mechanisms
//used in Compromise function
#define MINIMUM 0
#define MAXIMUM 1
#define MEAN 2

#define NONALC 0					//Never had unhealthy alcohol use
#define ALC 1						//Unhealthy alcohol use
#define ALC_TREATED_CURED 2			//Given an intervention and cured of unhealthy alcohol use
#define ALC_TREATED_NOT_CURED 3		//Given an intervention but not cured of unhealthy alcohol use

#define NON_IDU 0
#define IDU 1

#define SHARED_INJECTIONS_PER_YEAR 102			//Scott says vary shared injections per partnership from 10-100
#define IDU_PARTNER_CHANGE_RATE	5				//Range from 1-10

#define SBIRT_TARGET_START SUSC	 	//Tells alchol interv to target people at a status >= this.  (ex.  SUSC, HVL, INF, DET, CARE, TREAT, etc.)
#define SBIRT_TARGET_STOP TREAT		//Tells alcohol interv to target people at a status <= this. (ex.  SUSC, HVL, INF, DET, CARE, TREAT, etc.)
#define SBIRT_VL_TARGET_THRESHOLD 0  //VL cat at which alcohol interv targeting occurs (make sure to set this to 0 to target everyone, even if targeting >= SUSC!)
#define INDIA_ALC_TARGET_START TREAT //Tells the india individual and group interventions to target people at a status >=this (ex.  SUSC, HVL, INF, DET, CARE, TREAT, etc.)
#define INDIA_ALC_TARGET_STOP TREAT //Tells the india individual and group interventions to target people at a status <=this (ex.  SUSC, HVL, INF, DET, CARE, TREAT, etc.)


//For calculating LifeYears and QALYS
#define ALL_PATS 0
#define ALL_HIV_POS 1
#define ALL_CARE_AND_TREAT 2

//for calling getHIVProgRate which can now return a rate or a cost
//No longer needed, except for zero'd array - afraid to remove without testing first
#define RATE 0
#define COST 1

#define UT 0		//Utilization cost
#define PP 1		//Prepurchased cost
#define COMBO 2		//We want to use a mix of the two types of costs as defined in the costType[] array in interventions.h

#define NORM 0		//Normal not-discounted cost
#define DISC 1		//Discounted cost

#define INITIAL_POP 91496195		//Population of Maharashtra at start of calibration	

#define MAX_NUM_BASKETS 10000

#if defined(RUN_ALL_NO_SCRIPT) || defined(FIND_OPTIMUM_BASKET) /*|| defined(RUN_OPTIMISTIC_COMBOS_FOR_SCOTT)*/
#define RUNNING_MULTIPLE_BASKETS 1
#else 
#define RUNNING_MULTIPLE_BASKETS 0
#endif

#if RUN_TEST_INTERVENTIONS
#undef NUM_INTERVENTIONS
#define NUM_INTERVENTIONS 45
#endif

//Kendall Bryant Version
#define AVG_VL_DEC	2.41	//The average VL decrement from the adherence model =((2.68 + 3.09 + 2.22)/3 * 4.07 / 4.5)

#define RATE_WHEN_PROB_EQUALS_ONE 10

//lookup table files names, etc.
#define OFF_CARE_RATES "IN_100k_offcare_2015_05_06_0.txt"
#define CD4_200_ONCARE "IN_100k_oncare_200_TCD4_6mos_2regs_2015_05_06_"
#define CD4_350_ONCARE "IN_100k_oncare_350_TCD4_6mos_2regs_2015_05_06_"	
#define CD4_500_ONCARE "IN_100k_oncare_500_TCD4_6mos_2regs_2015_05_06_"	
#define CD4_10K_ONCARE "IN_100k_oncare_10000_TCD4_6mos_2regs_2015_05_15_"


#define RESOURCE_FILE_PATH "/ifs/data/hivmodels/India/pop_and_ratefiles/"
//#define RESOURCE_FILE_PATH "H:/Personal/Visual Studio 2013/Projects/India_pop_and_rate_files/"
//#define RESOURCE_FILE_PATH "C:/Users/kruggles7/Documents/Visual Studio 2012/Projects/India_pop_and_rate_files/"
//#define RESOURCE_FILE_PATH "/ifs/home/nucifk01/Rate_tables/EA/"

//Settings for Chris, Kelly and Kim's account on cluster
//#define RESOURCE_FILE_PATH "/ifs/home/toohec01/pop_and_rate_files/"
//#define RESOURCE_FILE_PATH "/ifs/home/toohec01/share_chris/pop_and_rate_files/"
//#define RESOURCE_FILE_PATH "/ifs/home/nucifk01/Rate_tables/EA/"
//#define RESOURCE_FILE_PATH "/ifs/home/nucifk01/Rate_tables/EA/revised_rates_2013_07_31/"
//#define RESOURCE_FILE_PATH "/home/kan25/pop_and_rate_files/"
//#define RESOURCE_FILE_PATH "/ifs/home/rugglk01/Braithwaite/Transmission_Model/India_pop_and_rate_files/"
//#define RESOURCE_FILE_PATH "/ifs/home/patela23/India_pop_and_rate_files/"

#define CD4_TREAT_THRESHOLD 200
#define CD4_TREAT_THRESHOLD_NOV_2011 350 
// Code will look for files named ON_CARE_ADHERENCE0.txt through ON_CARE_ADHERENCE10.txt


#define POPULATION_FILE "Population_2015_05_14.txt" //From 2015_04_15 calibration run


//These four are used in OPTIMIZING_ALCOHOL_INTERVENTION
#if OPTIMIZING_ALCOHOL_INTERVENTION
#define CONDOM_VALUE cdmv 
#define ADHERENCE_VALUE adhv
#define STI_VALUE stiv
#define ALCOHOL_VALUE alchv
#endif

//Used when running analyses with scripts.
//Each node on the cluster will get one numbered run.
#if ONE_WAY_SENSITIVITY_ANALYSIS || ISOLATING_EFFECTS_OF_ALCOHOL_INTERVENTION || RUN_ALL_USING_SCRIPT || INTERVENTION_EFFECT_SIZE_SENSITIVITY || COST_SENSITIVITY || OPTIMIZING_PREVENTION_PORTFOLIO || PROBABILISTIC_RUN
#define RUN_NUM 0
#endif



#define UTILITY_CD4_LESS_THAN_50 .79
#define UTILITY_CD4_50_TO_200 .85
#define UTILITY_CD4_ABOVE_200 .94
#define DELTA_UTILITY_WITH_HAART -.053
#define DISC_RATE .03

#if ALCOHOL_SURFACE
#define	SURF_COST 0
#define SURF_EFFECT 0
#define SURF_RR 1
#endif

//Turn off the HIV Prog Rate Efficiency array when calibrating (it will break the calibration!)
#if CALIBRATE==1
#undef EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS
#define EFFICIENCY_PRECALCULATE_HIV_PROG_RATE_AND_COST_ARRAYS 0
#undef EFFICIENCY_PRECALCULATE_ALPHA_MULT_VALUES
#define EFFICIENCY_PRECALCULATE_ALPHA_MULT_VALUES 0
#endif

#define SET_ALC_PROPORTIONS 0 	// Used with OPTIMIZING_PREVENTION_PORTFOLIO
#define SET_RR_ALC 0  			// Used with OPTIMIZING_PREVENTION_PORTFOLIO to change condoms baselineRR_raw[0][5] default(1.29), adherence baselineRR_raw[3][5] default(2.33), and STI baselineRR_raw[5][5] default(1.72). Use 1 for high and 2 for low settings.

//Because FSW are more likey to use condoms with their clients compared to their regular partner, we had to add a multiplier 
//can modify the baselineRR values in the getAlphaMult function.  The current value in the matrix is for a FSW partnering with her 
//regular partner (ACT1), so use this when they are partnering with a client (ACT2 and ACT3) 
#define FSW_CONDOM_MULTIPLIER 0.084 
#define NUM_YEARS_CALIBRATION 18

//To be deleted

//To be used for CarePlus intervention. Values should be changed through script if run on cluster, O/W, ERROR!
#define CONDOM_PATH_EFFECT 1
#define LINKAGE_PATH_EFFECT 1
#define ADHERENCE_PATH_EFFECT 1
#define CAREPLUS_INTERV_COST 10.90
#define PROP_ACT_MULT 0

#if COST_EFFECTIVENESS_SENSITIVITY || PROBABILISTIC_RUN
//When running the COST_EFFECTIVENESS_SENSITIVITY
#define NUM_INTERV_TO_COMBINE 7								//For determining the optimum basket, how many interventions are we going to potentially combine in one basket?
#define TOTAL_NUM_BASKETS 128								//Should equal NUM_INTERV_TO_COMBINE^2	
#define BASKET_USE bskt_count //input from the sh file CHANGE THIS BEFORE RUNNING! 
#define NUM_PROB_TESTS 84
#define NUM_PROB_COST_EFFECT 18
#endif 


#if ONE_WAY_SENSITIVITY_ANALYSIS
#define NUM_SENSITIVITY_TESTS  173
#endif
