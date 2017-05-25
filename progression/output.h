//debug settings
#define KIMTEST 0
#define VLCD4GRAPH 0													//kim_test2
#define ONLY_CORE_EVENTS 1

//analyses
//#define ONE_WAY_SENSITIVITY_ANALYSIS
//#define ONE_WAY_SENSITIVITY_ANALYSIS_RIC
//#define USING_COMMAND_LINE_ARGS		//used for running RIC sensitivity on cluster
//define FIND_BEST_FIT 
//#define TRIGGER_COMPARISON
//#define TEST_COST_COMPARISON
//#define CALIBRATION
//#define VARYING_AGE_VL_CD4
#define WHO_MONITORING_STRATEGY_ANALYSIS 0
//#define WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS
//#define RYAN_WHITE_ANALYSIS
#define GENERATING_TABLES_FOR_WHO 0

#if WHO_MONITORING_STRATEGY_ANALYSIS
#define SCENARIO 0
#endif


//For generating the transition rate tables for the transmission model
#define NUM_INITIAL_HAART_STATES 2 // Patients can either be forced to start on HAART or not. 
#define GETSTATES 0		//Turns on the code to keep track of the transitions
//#define GETSTATES_LOOP //Loops through all possible combinations of VL/CD4/Res categories
#ifdef GETSTATES_LOOP
#define COMMAND_LINE_RESISTANCE_ADHERENCE 0 // For using the cloud commandline with parameters
#define RUNNING_WITH_SCRIPT_ON_CLUSTER 1
#endif


//Allowable behaviors (to help debug)
#define USE_PERLIN_NOISE 1
#define VRT	0	
#define BATCH 0	

//output settings - all to go file1
#define SURV_RATES 1												
#define TIME_TO_TREATMENT_FAILURE_OF_REGIMENS_1_2_3 1				
#define CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_THERAPY 1		
#define CHANGE_IN_VL_and_CD4_1_3_5_YEARS_INTO_SALVAGE 1					
#define ANNALS_REGS_MUTS 1											
#define ACTIVE_HAART_INFO 1											
#define REASONS_FOR_REG_CHANGES 1								
#define MORTALITY_1_3_5_YEARS_AFTER_THERAPY 1						
#define DEBUG_NUM_MUTS 1	
#define DEBUG_AIDS 1
#define MONITORING_PAPER 1

//more output settings
#define HIV_VS_AGE_DEATHS 0												

#ifdef CALIBRATION								//these go to the file "output.txt"						
#define KAPLAN_MEIER_SURVIVAL_CURVE 0								
#define KAPLAN_MEIER_CURVE_OF_TIME_TO_TREATMENT_FAILURE_REG_1_2_3 0
#define YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY 1
#else
#define KAPLAN_MEIER_SURVIVAL_CURVE 0									
#define KAPLAN_MEIER_CURVE_OF_TIME_TO_TREATMENT_FAILURE_REG_1_2_3 0	
#define YEARLY_CD4_AND_YEARLY_CHANGE_IN_CD4_AFTER_START_OF_THERAPY 0							
#endif

#ifdef ONE_WAY_SENSITIVITY_ANALYSIS
#define REG_INFO_FOR_SENSITIVITY	1									//output
#else 
#define REG_INFO_FOR_SENSITIVITY	0
#endif

#if WHO_MONITORING_STRATEGY_ANALYSIS || defined(WHO_MONITORING_STRATEGY_DISCORDANCE_ANALYSIS) || GENERATING_TABLES_FOR_WHO
#define DO_CLINICAL_VISIT_EVERY_SIX_MONTHS_FOR_WHO 1
#else
#define DO_CLINICAL_VISIT_EVERY_SIX_MONTHS_FOR_WHO 0
#endif
