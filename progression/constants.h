#include "output.h"

#define PATIENTS 10000//0//0

#define RESOURCE_RICH 0			//Set as 1 for resource rich, 0 for resource poor

#define CALIBRATE 0 			//There are some settings only used during calibration

//#define COMPETING_RISKS_1
//#define COMPETING_RISKS_2

#define CYCLE DAY				//length of cycle can be MONTH or DAY 

#define ALC_INTERVENTION	0  //intervention for alcohol
#define DETECT_VL	0       //start alcohol intervention for patients with detactable VL
#define BOTH_AD_ALC  0		//both ad and alc interventions on, only turn this on when DETECT_VL is on, if this is on, ALC_INTERVENTION needs to be off.
#define VL_THRESHOLD	2.7  //threshold to start intervention
#define PROP_CURED_ALC  .9//probability alc cured

#define VL_DEC_MULT  1.613        //LL: 100/62 see Jason's email (April 07, 2015 10:24am)
#define CD4_MULT     1.613            //LL: 100/62 see Jason's email (April 07, 2015 10:24am)

#if !RESOURCE_RICH //India is a resource poor country

#if CALIBRATE

//#define SURVIVAL
//#define TTTF
#define CD4CHANGE
#define WHITE_CARD

//these are universal params
#define CD4_TREAT				10000	//Use 10000 during calibration so patients start in HAART immediately 
#define STOP_REG				0.12//0.04	//rate is per year was at .019 for previous good cal. 
										//increased from 0.04 to 0.12 when we incorprated changes to reg change trig and CD4 multiplier
#define INTOL_MULTIPLIER		3.0		//if a patient has previously been intolerant, they're likely to be intol again
#define MORTMULTHIV				1.0  //1.0											
#define TRIG					T_CD4    
#define MONITOR_INTERVAL		6		//in units of months
#define CENSOR_THRESHOLD 300

#define CENSOR_RATE				0//0.5	//.5	//Used only for resource poor
#define PROB_CENSOR_INSTEAD_OF_FAIL    0//0.5	//.5  //Patients who were about to fail a regimen have this probability of being censored instead
#define COINFECTION_MULT		3.0//6.0  //5.0		
#define TOX_MULT_1ST_3_MONTHS	10//10		//Multiplies toxicity by this factor over the first 3 months of haart

#define TIME_HIGH_MORT  0.5//0.46//0.42   //5 months  0.25 //3 months

//now those specific to cohort
#ifdef SURVIVAL
//#define LOWCD4
#ifdef LOWCD4
#define CD4_HIGH	99     //for India calibration
#define CD4_LOW		0   
#define START_WITH_AIDS			1.0 
//#else
//#define CD4_HIGH	199   //for India calibration
//#define CD4_LOW		100 
//#define START_WITH_AIDS			.686 
#endif

#ifdef WHITE_CARD
#define AVG_CD4 214.18
#define AVG_CD4_SD 173.7
#define START_WITH_AIDS			0.11 
#endif

#define AVG_AGE 41.35	
//#define AVG_CD4 50   //actually not used in calibration
#define AVG_HIV 4.5				//Unknown.  Will try 4.5 to 5.0
#define AVG_AGE_SD 9.1 
//#define AVG_CD4_SD 10
#define AVG_HIV_SD 1.11			//Unkown.  Using same as VACS
#define MALE		.597
#define FEMALE		.403
#define COMP					.744		//use C_compliance when generating adherence stratified rate tables
#define NUM_CLASS_COMBINATIONS 2 //TEST //1//2//3				//The number of different possible class combinations
#endif

#ifdef TTTF
#define START_WITH_AIDS			0.678  //for CD4 50-99

#define AVG_AGE 36	
#define AVG_CD4 50 
#define AVG_HIV 5.18				//Unknown.  Will try 4.5 to 5.0
#define AVG_AGE_SD 10 
#define AVG_CD4_SD 10
#define AVG_HIV_SD 1.11			//Unkown.  Using same as VACS
#define MALE		.74
#define FEMALE		.26
#define COMP					.744		//use C_compliance when generating adherence stratified rate tables
#define NUM_CLASS_COMBINATIONS 3 //TEST //1//2//3				//The number of different possible class combinations

#ifndef WHITE_CARD
#define CD4_CUT1	50
#define CD4_CUT2	100
#define CD4_CUT3	200
#define CD4_PROP1	0.344
#define CD4_PROP2	0.183
#define CD4_PROP3	0.294
#endif

#endif

#ifdef CD4CHANGE
#ifdef WHITE_CARD
#define START_WITH_AIDS			0.11 
#define AVG_AGE 41.35	
#define AVG_CD4 214.48 
#define AVG_HIV 4.5				//Unknown.  Will try 4.5 to 5.0
#define AVG_AGE_SD 9.1 
#define AVG_CD4_SD 173.7
#define AVG_HIV_SD 1.11			//Unkown.  Using same as VACS
#define MALE		.597
#define FEMALE		.403
#define COMP					.74		//use C_compliance when generating adherence stratified rate tables
#else
#define START_WITH_AIDS			0.34  //0.175	//From Kara:  11% had a stage 4 diagnosis at the time of ART initiation.
#define AVG_AGE 42.4//33	
#define AVG_CD4 198.94//180 
#define AVG_HIV 5.18//4.5				//Unknown.  Will try 4.5 to 5.0
#define AVG_AGE_SD 8.06//6.7 
#define AVG_CD4_SD 144.04//211.2
#define AVG_HIV_SD 1.11			//Unkown.  Using same as VACS
#define MALE		.6188//.595
#define FEMALE		.3812//.405
#define COMP					1.0//.744		//use C_compliance when generating adherence stratified rate tables
#endif
#define NUM_CLASS_COMBINATIONS 2 //TEST //1//2//3				//The number of different possible class combinations
#endif

#else 
//not for calibration
#define CD4_TREAT				10000	//200 
#define STOP_REG				0.12//0.04	//rate is per year was at .019 for previous good cal. 0.12 is used for the version with new reg change trigger and CD4 mult
#define INTOL_MULTIPLIER		3.0		//if a patient has previously been intolerant, they're likely to be intol again
#define MORTMULTHIV				1.0  											
#define TRIG					T_CD4    
#define MONITOR_INTERVAL		6		//in units of months
#define CENSOR_THRESHOLD 300

#define CENSOR_RATE				0	//.5	//Used only for resource poor
#define PROB_CENSOR_INSTEAD_OF_FAIL    0	//.5  //Patients who were about to fail a regimen have this probability of being censored instead
#define COINFECTION_MULT		3.0		
#define TOX_MULT_1ST_3_MONTHS	10//50//15		//Multiplies toxicity by this factor over the first 3 months of haart

#define TIME_HIGH_MORT  0.5//0.42   //5 months  0.25 //3 months

#define START_WITH_AIDS			0.11  //From Kara:  11% had a stage 4 diagnosis at the time of ART initiation.

#define AVG_AGE 41.35	
#define AVG_CD4 214.48 
#define AVG_HIV 4.5				//Unknown.  Will try 4.5 to 5.0
#define AVG_AGE_SD 9.1 
#define AVG_CD4_SD 173.7
#define AVG_HIV_SD 1.11			//Unkown.  Using same as VACS
#define MALE		.597
#define FEMALE		.403
#define COMP		0.74//0.609//0.945//0.886//0.763//0.72//0.695//0.609		//use C_compliance when generating adherence stratified rate tables
#define NEWCOMPAD		0.945//0.945//0.886//0.763//0.72//0.695//0.609		//new level for ad intervention
#define NEWCOMPALC		0.969//0.969//0.934//0.855//0.824//0.806//0.74 //adherence level if cured alc

#define NUM_CLASS_COMBINATIONS 2 //TEST //1//2//3				//The number of different possible class combinations
#endif

//All patients will begin on this regimen
#define INITIAL_REG_DRUG1 NNRTI_NEVIRAPINE
#define INITIAL_REG_DRUG2 NRTI_NONTAM
#define INITIAL_REG_DRUG3 NRTI_TAM

//Number of drugs in each class
#define NUM_PI_SINGULAR 0 //Nelfinavir is the only unboosted PI
#define NUM_PI_BOOSTED 2 //0//1//2	
#define NUM_NRTI_TAM 2 //1//2
#define NUM_NRTI_NONTAM 4 //1//2//4 

#ifdef CD4CHANGE
#define NUM_NNRTI_EFAVIRENZ 1
#else 
#define NUM_NNRTI_EFAVIRENZ 0
#endif

#define NUM_NNRTI_NEVIRAPINE 1	

//the VL threshold for regimen failure before and after the first regimen failure 
//Used in MONITORING_PAPER
#define VL_REGIMEN_FAILURE_1ST 3.7
#define VL_REGIMEN_FAILURE_AFTER_1ST 3.7

#endif /******************************* end POOR *********************************/

#define PATIENT_STARTS_ON_HAART 0 		//If value is 1, patient is initialized to be on haart at the start of the model.
#define	VACS_HIV_NEG_MORT_RATE  0.01526	//The mortality rate of HIV negative patients in the VACS cohort
#define VACS_HIV_NEG_MORT_MULT	2.5		//We doubled the HIV negative mortality during calibration to
										//bring the survivals in line

#define RES_POOR_HIV_NEG_MORT_MULT 1.8  

#define POOR_CAN_GO_BACK_ON_PREVIOUS_REGIMENS	1 //We may want to allow a person in a resource poor setting to go back n previous regimens, only applies to resource poor

#define LUCK_MULTIPLIER	0.5						//to control the effect of patient luck on the CD4within


#define MUT_RATE  0.18 //0.18 per year equals 0.015 per month
#define FACTOR 3.3		//Used in determining the mutation rate based on viral load

#define DISC_RATE .03
#define HAART_TOX 1.0  //toxicity

#define DELTA_UTILITY_WITH_HAART -.053
#define UTILITY_CD4_LESS_THAN_100 .81
#define UTILITY_CD4_LESS_THAN_200 .87
#define UTILITY_CD4_ABOVE_200 .94

#define ARV_TYPES 6		//represents the number of classes: PI, NRTI, NNRTI (x2)
#define DRUGS_IN_REG 3  //3 indicates triple-drug therapy

//Probability that mutation results in resistance
#define PMUTRES_PI_SINGULAR_EST  .5		
#define PMUTRES_PI_BOOSTED_EST  .5	
#define PMUTRES_NRTI_TAM_EST  .5
#define PMUTRES_NRTI_NONTAM_EST  .5
#define PMUTRES_NNRTI_EFAVIRENZ_EST  .9
#define PMUTRES_NNRTI_NEVIRAPINE_EST  .9

//Probability that mutation results in cross-resistance (to other drugs within class)
#if RESOURCE_RICH							//All PIs are boosted except for Nelfinavir
#define PCROSSRES_PI_SINGULAR_EST 1.0 
#define PCROSSRES_PI_BOOSTED_EST  .24
#else
#define PCROSSRES_PI_SINGULAR_EST .22		//All PIs are unboosted.  This includes all PIs
#define PCROSSRES_PI_BOOSTED_EST .24//1.0		//except for Lopinavir and Saquinavir
#endif
#define PCROSSRES_NRTI_TAM_EST 1.0
#define PCROSSRES_NRTI_NONTAM_EST  .48
#define PCROSSRES_NNRTI_EFAVIRENZ_EST 1.0
#define PCROSSRES_NNRTI_NEVIRAPINE_EST 1.0

//Probability of cross-resistance between different classes
#define PCROSSRESPI_BOOSTED_SINGULAR .05
#define PCROSSRESPI_SINGULAR_BOOSTED .07
#define PCROSSRESNRTI_TAM_NONTAM .02
#define PCROSSRESNRTI_NONTAM_TAM .08
#define PCROSSRESNNRTI_EFAV_NEVIR 1.0
#define PCROSSRESNNRTI_NEVIR_EFAV .875

//ratio of probability of mutations between drug classes
#define RATIO_PI_TO_NRTI	1.0		//.2		
#define RATIO_NNRTI_TO_NRTI	1.0

#define MUT_PI_SINGULAR_START 0
#define MUT_PI_BOOSTED_START 0	
#define MUT_NRTI_TAM_START 0
#define MUT_NRTI_NONTAM_START 0
#define MUT_NNRTI_EFAVIRENZ_START 0
#define MUT_NNRTI_NEVIRAPINE_START 0

//the probability of an AIDS defining event occurring is this factor times the probability of
//death from the HIV lookup table.  (Source: ARTCC - Antiretroviral Therapy Cohort Collaboration)
#define AIDS_EVENT_MULTIPLIER	3.0

//The number of drugs the patient must be resistant to on all regimens before going on salvage
#define NUM_DRUGS_RES_FOR_PLATO 2

//The probability of death from HIV is adjusted for whether the pat has AIDS
//Source:  AIDS 2007, 21:1185-1197, 2.33 (2.01-2.69) 95% confidence interval
//The NON_AIDS_ADJUST was calculated based on a .310 initial AIDS rate in the CHORUS
//cohort.  Source:  Journal of Clinical Epidemiology 57 (2004) 89-97.
#define AIDS_ADJUST				2.33		
#define NON_AIDS_ADJUST			.401

//The proportion of time we decrement for ALL of the drugs in a patient's drug regimen vs. choosing one at random to decrement.
#define INTOLERANCE_DING_FOR_ALL .25	

//Switch determines whether it is possible to develop mutations against a drug in the regimen that the patient has already developed resistance to.	
#define CAN_GET_MUTS_TO_RES_DRUGS 0	

#define VL_DEC_PI_SINGULAR 1.84//1.78 //=(1.12 * 100/63)
#define VL_DEC_PI_BOOSTED 2.68//2.51	//=(1.58 * 100/63)
#define VL_DEC_EFAVIRENZ 3.09//3.29	//=(2.07 * 100/63)	
#define VL_DEC_NEVIRAPINE 2.22//2.29	//=(1.44 * 100/63)


#define WEIGHT_PAT_LUCK		.8
#define WEIGHT_ROUND_LUCK	.2

#define VL_ADJUST 1.5  //kn changed for v2.0.  Was 3.0 in 1.5c.
//the following is obtained so that if we approach something at the rate of (1/3.13) per month, then
//after 6 months we'll be 90% of the way there
#define CD4_ADJUST_W 3.13  //3.0 
//the following is obtained so that if we approach something at the rate of (1/13.51) per month, then
//after 30 months we'll be 90% of the way there

#define RES_BASED_CHOICE 0	//for patients with initial mutations prior to starting haart.
							//set to 1 if we want to pick a good initial regimen that makes use of the 
							//information on number of baseline mutations against each haart type

//These costs are from Anik's lit search.  

#define TESTCOSTVL 49.54
#define TESTCOSTCD4 6.32

#define COST_OF_CARE_POOR 132.18			//Annual cost of care for Resource Poor
#define HOSPITAL_COST_POOR_LOW_CD4 347.25			//Annual costs for hospitalization for CD4<200
#define HOSPITAL_COST_POOR_MEDIAN_CD4 40.85			//Annual costs for hospitalization for 200<CD4<350
#define HOSPITAL_COST_POOR_HIGH_CD4 5.72			//Annual costs for hospitalization for CD4>350

#define COST_REG1_POOR 11.86				//Monthly costs of reg1;
#define COST_REG2_POOR 49.27				//Monthly costs of reg2;
#define COST_REG3_POOR 119.55				//Monthly costs of reg3:  for sensitivity $113.43, $255.60, $1022.41 

#define COST_ALC_INTERVENTION	1 // a dollar for a year

#define MAX_MONTHS 2402
#define MAX_REGS	50		//No limit on number of regimens for resource rich

#define MAX_ARV_TYPES ARV_TYPES //50
#define MAX_DRUGS_IN_CLASS 8	

#define MAX_REGS_MODEL	MAX_REGS	//This is the number of possible regimens the model can store data for

/*From Rob's e-mail - Each one will produce a different stream of random numbers. 
We've been sticking to the 60076 for as long as I can remember.  
In case you're curious about where they came from, 60076 was Steven's zipcode 
when he was growing up and 15232 was my zipcode when i was growing up.  
60076 produces a much longer stream of randoms, so it's easier to just use that one. */

#define INITIAL_SEED 60076

#define MAX_AGE 80
#define MAX_CD4 1500
#define MAX_HIV	 8

#define TIME 50						//Number of years of model run

//for sensitivity analyses
#define DVLHAART 0 //baseline: 0 //-1, 1
#define DCD4NOHAART 0 //baseline: 0 //-200, 200
#define DCD4HAART	0 //baseline: 0 //-200, 200
#define ADJMUTRES 1 //baseline: .5 //0, 1 
#define ADJCROSSRES 1 //baseline: 1 //0, 10
#define MORTMULTNOHIV 1 //baseline: 1//.1, 10

#define PAT_CHANGEREG 1
#define RES_CHANGEREG 2
#define RES_ALLREG 3
#define DIE_AGE 4
#define DIE_HIV 5

#define VL 1
#define CD4 2

#define TRUE 1
#define FALSE 0

//used for breaking down costs
#define TOTAL 0
#define DRUG 1
#define LAB 2
#define CARE 3
#define HOSP 4
#define OUTREACH 5  //cost of OUTREACH_INTERVENTION
#define RISK_REDUCTION 6 //cost of RISK_REDUCTION_INTERVENTION
#define SECONDARY_PREVENTION 7 //cost of SECONDARY_PREVENTION
#define ALC_INTERV  8  //cost of alc intervention

#define TABLE_SPACE 400

#ifdef COMPETING_RISKS_1
#define NUM_MORT_CONDITIONS 18 //leading causes of death + other death
#endif

#ifdef COMPETING_RISKS_2
#define NUM_MORT_CONDITIONS 18 //leading causes of death + other death
#endif

#define CALIBRATION_YEARS_POOR 5		//We have calibration data for up to five years for resource poor (may be less)
#define NUM_YEARS_VACS_HIV_NEG 14		//We have 14 years worth of mortality data from vacs HIV negs
#define HIV 0	
#define AGE 1
#define CENSORED 2
#define EPSILON .00001
#define DAY 1
#define MONTH 2

//Regimen types used in function ChooseRegimen()
#define MAINTAIN_CURRENT_REG 0
#define NEW_REG 1

//Used in putPatOnReg
#define NONE_AVAILABLE -999

#define PI_SINGULAR 0
#define PI_BOOSTED 1			
#define NRTI_TAM 2
#define NRTI_NONTAM 3
#define NNRTI_EFAVIRENZ 4			
#define NNRTI_NEVIRAPINE 5

//Possible regimen change triggers (for comparison)
#define T_CD4		1
#define T_VL		2
#define T_NESTED	3			//If CD4 failure criteria is met, then check VL
#define T_CLINICAL	4

//For VRT
#define MAX_SEED_INIT 8 
#define MAX_SEED_CYCLE 510 


//All rate and time variables are in units of years unless otherwise specified.
//Agreed with Scott to round the number of days per year to 365 (ignoring leap years).
//365 and 365.0 should give the same results
//Note:  I've kept "CYCTIME_IN_MONTHS" because some legacy code was just too complicated 
//to change when switching over to a day cycle time.

#define DAYS_PER_YEAR			365		

#if CYCLE == DAY
	#define CYCTIME				1.0/DAYS_PER_YEAR 		
	#define CYCLES_PER_YEAR		DAYS_PER_YEAR	
	#define CYCTIME_IN_MONTHS	12.0/DAYS_PER_YEAR	//CYCTIME_IN_MONTHS is the number of months per cycle
													//For daily, .032877 (one day = .032877 months)  
#ifdef CD4CHANGE
#define CYCLES_PER_MONTH	30 
#else
#define CYCLES_PER_MONTH	DAYS_PER_YEAR / 12.0
#endif
#else							
	#define CYCTIME				1.0/12.0			//CYCTIME is the number of months per cycle
	#define CYCLES_PER_YEAR		12.0
	#define CYCTIME_IN_MONTHS	1.0					//	One cycle is equal to one month
#endif


#define MARK_ROBERTS_AHRQ 0

#if MARK_ROBERTS_AHRQ

#define AVG_AGE 38	
#define AVG_CD4 191
#define AVG_HIV 4.8
#define AVG_AGE_SD 0
#define AVG_CD4_SD 150
#define AVG_HIV_SD 1.0 
#define MALE		.80
#define FEMALE		.2
#define INITIAL_REG_DRUG1 NNRTI_EFAVIRENZ  //NNRTI_EFAVIRENZ or PI_BOOSTED
#define INITIAL_REG_DRUG2 NRTI_NONTAM
#define INITIAL_REG_DRUG3 NRTI_NONTAM
#define VL_DEC_PI_SINGULAR 2.59  //=((2.07 -.44)* 100/63)	
#define VL_DEC_PI_BOOSTED 2.86	//=((2.07 -.27)* 100/63)	
#define VL_DEC_EFAVIRENZ 3.29	//=((2.07 -0)* 100/63)		//Efav is the reference group
#define VL_DEC_NEVIRAPINE 2.90	//=((2.07 -.24)* 100/63)	

#define VL_REGIMEN_FAILURE_1ST 3.7
#define VL_REGIMEN_FAILURE_AFTER_1ST 3.7

#define STOP_REG				.099		//rate is per year	

#define MARK_VL_DEC_MULTIPLIER 1.575		//1.55 looks good

#define AIDS_EVENT_MULTIPLIER	7.5//10.0		//default is 3.0

#endif

//more for above analysis for Mark
#define MARK_VL_DEC_MULTIPLIER 1.575		
#define	TIMEPOINT 2				//num years at which to produce output (also part of MARK_ROBERTS_AHRQ)

//used for generating transmission model lookup table
#define STCD4 5
#define STVL 5
#define STTREAT 2
#define STRESIST 8
#define STDEAD 2



//For running the sensitivity analysis of the monitoring paper
//#define MONITORING_PAPER_SENSITIVITY

#define RUN_NUM 99			//can leave defined, won't have any effect

#ifdef MONITORING_PAPER_SENSITIVITY

#if RUN_NUM == 10
#undef COMP
#define COMP					.75
#endif

#if RUN_NUM == 11
#undef COMP
#define COMP					.95
#endif

#if RUN_NUM == 14
#undef AVG_AGE
#define AVG_AGE 20	
#define AVG_AGE_SD 0
#endif

#if RUN_NUM == 15
#undef AVG_AGE
#define AVG_AGE 50	
#define AVG_AGE_SD 0 
#endif

#if RUN_NUM == 16
#undef MUT_RATE
#define MUT_RATE  0.16
#endif

#if RUN_NUM == 17
#undef MUT_RATE
#define MUT_RATE  0.20
#endif

#if RUN_NUM == 20
#undef COST_OF_CARE_POOR
#define COST_OF_CARE_POOR .5*287.28
#endif

#if RUN_NUM == 18
#undef DELTA_UTILITY_WITH_HAART
#define DELTA_UTILITY_WITH_HAART  0
#endif

#if RUN_NUM == 19
#undef DELTA_UTILITY_WITH_HAART
#define DELTA_UTILITY_WITH_HAART  -.10
#endif

#if RUN_NUM == 21
#undef COST_OF_CARE_POOR
#define COST_OF_CARE_POOR 1.5*287.28
#endif

#if RUN_NUM == 22
#undef HOSPITAL_COST_POOR 
#define HOSPITAL_COST_POOR .5*390.27			//Annual costs for hospitalization for Resource Poor
#endif

#if RUN_NUM == 23
#undef HOSPITAL_COST_POOR 
#define HOSPITAL_COST_POOR 1.5*390.27			//Annual costs for hospitalization for Resource Poor
#endif

#if RUN_NUM == 24
#undef TESTCOSTVL
#define TESTCOSTVL .5*70.00
#endif

#if RUN_NUM == 25
#undef TESTCOSTVL
#define TESTCOSTVL 1.5*70.00
#endif

#if RUN_NUM == 26
#undef TESTCOSTCD4
#define TESTCOSTCD4 .5*11.20
#endif

#if RUN_NUM == 27
#undef TESTCOSTCD4
#define TESTCOSTCD4 1.5*11.20
#endif

#if RUN_NUM == 28
#undef COST_REG3_POOR
#define COST_REG3_POOR 113.43				//Monthly costs of reg3:  for sensitivity $113.43, $255.60, $1022.41 
#endif

#if RUN_NUM == 29
#undef COST_REG3_POOR
#define COST_REG3_POOR 1022.41				//Monthly costs of reg3:  for sensitivity $113.43, $255.60, $1022.41 
#endif

#if RUN_NUM == 30
#undef COST_REG2_POOR
#undef COST_REG3_POOR
#define COST_REG2_POOR 15.78				//Monthly costs of reg2;
#define COST_REG3_POOR 15.78				//Monthly costs of reg3:  for sensitivity $113.43, $255.60, $1022.41 
#endif

#if RUN_NUM == 31
#undef COST_REG2_POOR
#undef COST_REG3_POOR
#define COST_REG2_POOR 1.5 * 113.43				//Monthly costs of reg2;
#define COST_REG3_POOR 1022.41				//Monthly costs of reg3:  for sensitivity $113.43, $255.60, $1022.41 
#endif

#if RUN_NUM == 32
#undef MUT_NNRTI_NEVIRAPINE_START
#define MUT_NNRTI_NEVIRAPINE_START 1
#endif

#if	RUN_NUM == 33
#undef POOR_CAN_GO_BACK_ON_PREVIOUS_REGIMENS
#define POOR_CAN_GO_BACK_ON_PREVIOUS_REGIMENS 0 //We may want to allow a person in a resource poor setting to go back n previous regimens, only applies to resource poor
#endif

#endif // end of #ifdef MONITORING_PAPER_SENSITIVITY

/**************************** Retention in Care Parameters **********************************/

#define RETENTION_IN_CARE 0     //LL: add Retention In Care in progression model

#if RETENTION_IN_CARE
#define OUTREACH_INTERVENTION 1
#define RISK_REDUCTION_INTERVENTION 1
#define SECONDARY_PREVENTION_INTERVENTION 1

//Patients will start on ART at the CD4_TREAT threshold with a probability of PROB_ART_INITIATION_AT_CD4_TREAT
//If they don't start at CD4_TREAT, they will start when their CD4 drops below 200.
#define PROB_ART_INITIATION_AT_CD4_TREAT 1		

/*
engagement status (ES) ES1: pat is in care and in clinic
engagement status (ES) ES2: pat is not in clinic (transition state)
engagement status (ES) ES3: pat is not in clinic but is in care
engagement status (ES) ES4: pat is not in clinic and is not in care
engagement status (ES) ES5: pat is dead (transition state)
engagement status (ES) ES6: DEAD_KNOWN
engagement status (ES) Es7: DEAD_UNKNOWN
*/

//The following are daily probabilities
#define PRE_ART_MULT_FOR_PROB_ES1_TO_ES2 2      // if not yet eligible for ART, apply this multiplier to prob_es1_to_es2. based on Jason's 11/4/2013 email
#define PROB_ES1_TO_ES2_6MON       0.000828     // 0-6 months, per day, per Jason's 11/4/2013 email
#define PROB_ES1_TO_ES2_12MON      0.000399     //6-12 months, per day, per Jason's 11/4/2013 email
#define PROB_ES1_TO_ES2_12MON_MORE 0.000141     // >12 months, per day, per Jason's 11/4/2013 email
#define PROB_ES2_TO_ES4 0.27                    // doc name "Retention In Care Logic_11_2013" by Jason
#define PROB_ES5_TO_ES6 0.42                    //doc name "Retention In Care Logic_11_2013" by Jason
#define PROB_ES4_TO_ES1_AIDSPAT 0.5
#define MUTATION_MULT_ES4 3.37        //multiplier on mutation rate if in ES4>1 month and restarts ART. doc name "Retention In Care Logic_11_2013" by Jason

//for OUTREACH_INTERVENTION
//#define BASE_MULTIPLIER_FOR_OUTREACH 3	//by this point, multiplier_for_outreach = base_multiplier_for_outreach * P1 *P2 * P3
#define OUTREACH_PROB_IDENTIFY 1.0			//doc name "Retention In Care Logic_11_2013" by Jason
#define OUTREACH_PROB_FIND 0.0046			//doc name "Retention In Care Logic_11_2013" by Jason
#define OUTREACH_PROB_RELINK 0.0051			//doc name "Retention In Care Logic_11_2013" by Jason
#define OUTREACH_START_TRIGGER 90		    //90 days from disengagement of ES1. doc name "Retention In Care Logic_11_2013" by Jason
#define OUTREACH_END_TRIGGER 90             //90 days from the initiation of outreach (duration).  doc name "Retention In Care Logic_11_2013" by Jason
#define OUTREACH_INITIATION_COST_POOR 0     //one-time cost ??? no source yet
#define OUTREACH_FINDING_COST_POOR    0.15  //Daily cost;doc name "Retention In Care Logic_11_2013" by Jason
#define OUTREACH_FINISHING_COST_POOR  0     //one-time cost ??? no source yet
#define OUTREACH_CD4 10000					//Intervention will only kick in if patient CD4 is below this value			

//for RISK_REDUCTION_INTERVENTION
#define RISK_REDUCTION_RR_PROB_ES1_TO_ES2 0.6                 // doc name "Retention in care analysis_v1.0" by Jason
#define RISK_REDUCTION_INTERVENTION_COST 0.33  // $10/person/month; doc "Retention in care analysis_v1.0"
#define RISK_REDUCTION_CD4 10000			//Intervention will only kick in if patient CD4 is below this value			

//for SECONDARY_PREVENTION_INTERVENTION
#define SECONDARY_PREVENTION_RR_PROB_ES1_TO_ES2 0.6                 // Jason says use same as Risk Reduction for now
#define SECONDARY_PREVENTION_INTERVENTION_COST 0.33  // Jason says use same as Risk Reduction for now
#define SECONDARY_PREVENTION_CD4 10000		//Intervention will only kick in if patient CD4 is below this value			



#endif

/**************************** End of Retention in Care Parameters **********************************/


#if COMMAND_LINE_RESISTANCE_ADHERENCE
#undef COMP
#define COMP 0  // RCAT and COMP will not be used when using command line arguments, but the compiler will complain if they're not defined.
#define RCAT 0
#elif RUNNING_WITH_SCRIPT_ON_CLUSTER
#undef COMP
#define COMP					C_compliance//.85		//use C_compliance when generating adherence stratified rate tables
#define RCAT					R_resistance			// 0 - 7 Initial resistance category. for automation of runs, edit R_resistance and C_compliance.
#else
#define RCAT					0			// 0 - 7 Initial resistance category. for automation of runs, edit R_resistance and C_compliance.
#endif



