#if ONE_WAY_SENSITIVITY_ANALYSIS
bool basket[NUM_INTERVENTIONS] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#else 
bool basket[NUM_INTERVENTIONS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, /*India specific starts here*/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
//bool basket[NUM_INTERVENTIONS] = {0, 1, 2, 3, 4, 5, 6 ,7, 8, /*India specific starts here*/9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25} //for sanity's sake
//bool basket[NUM_INTERVENTIONS] = {1, 1, 1, 1, 1, 1, 1, 1, 1 };
#endif



struct intervention 
{
	//ANIK
	double who[NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU][NUM_VL];	//gets set to "1" if the intervention is applied to everyone in that comparment, "0" if nobody in that compartment, "0.5" for half the people, etc.
	double cost_per_person;								//the cost of implementing the intervention per person
	double total_prepurchased_cost, total_utilization_cost;	//the total cost of this intervention so far during the model run
	double discounted_total_prepurchased_cost, discounted_total_utilization_cost;
	double effect_size[NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE];
	double overall_effect_size;							//used to help determine utilization cost
		
	int status, sex, orient, act, age, alc, idu, vl;	

	void initialize() 
	{
		cost_per_person=0;						//eventually gets set to interv_costs_per_person[NUM_INTERVENTIONS], shown below, in the function initialize_interventions()	
		total_prepurchased_cost=0;
		total_utilization_cost=0;
		overall_effect_size = 0;
		discounted_total_prepurchased_cost = 0;
		discounted_total_utilization_cost = 0;

		LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS
			who[status][sex][orient][act][alc][idu][vl] = 0;
		END_LOOP_THROUGH_ALL_INTERVENTION_COMPARTMENTS
	}
};

intervention interv[NUM_INTERVENTIONS];			//create an array of interventions

//the array to hold the baseline probability of each compartment over each pathway and characteristic
double baselineProb[NUM_PATH_AND_CHAR][NUM_STATUS][NUM_SEX][NUM_ORIENT][NUM_ACT][NUM_ALC][NUM_IDU];	


//The values of the baseline relative risks from the second tab of the "master" spreadsheet
//Note that there are some special cases hardcoded into the function setBaselineProbs()
double baselineRR_raw[NUM_PATH_AND_CHAR][NUM_BASELINE_RR] = {
/* 0 */{.73, 1.08, .63, 1, NOT_APPLICABLE, .77, 1.14, 1.29, .99, .47, .47, 1, .47},
/* 1 */{.98, 1, .32, .31, 1, .81, .23, 1, .60, NOT_APPLICABLE, NOT_APPLICABLE, 1, NOT_APPLICABLE},
/* 2 */{.70, 1, 1, 1, 1, 1, 1, 1, 1, 1, NOT_APPLICABLE, NOT_APPLICABLE, NOT_APPLICABLE},
/* 3 */{.26, 1, 1, 1, 1, 1, 1, 2.33, 2, NOT_APPLICABLE, 1, NOT_APPLICABLE, 1},
/* 4 */{1, 1, 1, 1, 1, 1, 1, 1, 1, NOT_APPLICABLE, NOT_APPLICABLE, 1, NOT_APPLICABLE},
/* 5 */{.06, 1, 1, 1, NOT_APPLICABLE, 1, 8.53, 1.72, 1.60, 1, 1, 1, 1},
/* 6 */{1, 1, 1, 1, NOT_APPLICABLE, 1, 1, 1, 1, 1, 1, 1, 1},
/* 7 */{.8, NOT_APPLICABLE, 1, 1, 1, 1, 1, 1, 1, 2.22, 2.22, 2.22, 2.22}
};

#if !DELPHI
string intervention_name[NUM_INTERVENTIONS] = {
	"Condom distribution, risk reduction counselling, and STI treatment for female CSW", "Expanded VCT and universal ART", "PMTCT", "Alcohol/SBIRT intervention", "Voluntary Male Circumcision", "Mass media HIV education campaigns", "Structural condom distribution", "PrEP", "Care Plus", 
	"Alcohol Individual Short","Alcohol Individual Long","Alcohol Group Short","Alcohol Group Long", "Alcohol Community", "Depression Individual Short", "Depression Individual Long", 
	"Depression Group Short ", "Depression Group Long", "Sex Individual Short", "Sex Individual Long", "Sex Group Short", "Sex  Group Long", "Sex Community", "Weekly SMS", "Brief Adherence Counseling", "Mass media HIV education campaigns"
};
#else 
string intervention_name[NUM_INTERVENTIONS] = {
	"Condom distribution, risk reduction counselling, and STI treatment for female CSW", "Expanded VCT and universal ART", "PMTCT", "Alcohol/SBIRT intervention", "Voluntary Male Circumcision", "Mass media HIV education campaigns", "Structural condom distribution", "PrEP", "Care Plus", 
	"Alcohol Individual", "Alcohol Individual Long", "Alcohol Group ", "Alcohol Group Long", "Alcohol Community", "Depression Individual ", "Depression Individual Long",
	"Depression Group ", "Depression Group Long", "Sex Individual", "Sex Individual Long", "Sex Group", "Sex  Group Long", "Sex Community", "Adherence Individual", "Adherence Group", "Mass media HIV education campaigns"
};
#endif 
string path_names[NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE] = {
	"0: Prop w/ condom nonuse", "1: Prob of not being tested", "2: Prob of not being linked", "3: Prob of nonadherence", "4: Prob of not using prophylactic ARV", "5: Prob of untreated STI", "6: Prob not using microbicides", "7: Prob of not being circumsized", "8: Pathway: Fewer partners", "9: Pathway: Reduction in transactional sex or intergenerational partnering", "10: Pathway: Less substance abuse", "11:  PMTCT", "12:  Pathway: Change in the ART start threshold", "13: Proportion with >1 concurrent partners"
};

//Normal settings.  The one-way-sensitivity settings are at the bottom of this file.
//the intervention costs per person
#if !DELPHI
double interv_costs_per_person[NUM_INTERVENTIONS] = { 290, 39.09, 64.86, 5.20, 58, 1, 5.20, 165, CAREPLUS_INTERV_COST, 1.64, 6.56, 1.64, 5.9, 0, 13.12, 36.08, 3.608, 7.216, 1.64, 14.76, 1.64, 5.904, 6.67, 6.56, 2.46, 1.33};
#else
double interv_costs_per_person[NUM_INTERVENTIONS] = { 290, 39.09, 64.86, 5.20, 58, 1, 5.20, 165, CAREPLUS_INTERV_COST, 2.73, 0, 0.98, 0, 0.005, 2.73, 0, 0.71, 0, 2.73, 0, 0.19, 0, 0, 2.19, 0.3, 0};
#endif 

#if	COST_SENSITIVITY
double interv_costs_per_person_high[NUM_INTERVENTIONS] = {435, 58.5, 97.4, 7.8, 69, 1.5, 10.4, 18.93};
double interv_costs_per_person_low[NUM_INTERVENTIONS] = {145, 19.5, 32.4, 2.6, 47, .5, 1.04, 2.86};
string high_low;
#endif

//Used for determining utilization cost
#if !DELPHI
double overall_effect_size[NUM_INTERVENTIONS] = { .51, .58, .25, .8, .41, .93, .81, .5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.93}; 
double run_in_optimization[NUM_INTERVENTIONS]={0, 0, 0, 0, 0, 0, 0, 0, 0,/*start india*/ 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0 };
#else
double overall_effect_size[NUM_INTERVENTIONS] = { .51, .58, .25, .8, .41, .93, .81, .5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; 
double run_in_optimization[NUM_INTERVENTIONS] = { 0, 0, 0, 0, 0, 0, 0, 0, 0,/*start india*/1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0};
#endif 

double interv_effect_size[NUM_INTERVENTIONS][NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE]= 
{
#if !DELPHI //literature results


	/*Pathways in comments
	0		1		2		3  4  5		6  7	8		9  10 11				12							13  */
	{ .51, 1, 1, 1, 1, .31, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, .89 },	//0     Condoms, risk reduction conseling, and STI treatment for female CSW
	{ 1, .59, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_10000, 1 },	//1     Expaned VCT and universal ART
	{ 1, .24, 1, 1, 1, 1, 1, 1, 1, 1, 1, .8, CD4_TREAT_THRESH_DEFAULT, 1 },	//2     PMTCT
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .55, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//3     Alcohol/SBIRT
	{ 1, 1, 1, 1, 1, 1, 1, .4, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//4     Circumcision
	{ .91, .93, 1, 1, 1, 1, 1, 1, .89, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//5     Mass media HIV education campaigns     
	{ .81, 1, 1, 1, 1, .69, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//6     Structural condom distribution
	{ 1, 1, 1, 1, .5, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //7 PrEP
	{ CONDOM_PATH_EFFECT, 1, LINKAGE_PATH_EFFECT, ADHERENCE_PATH_EFFECT, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//8     Care+ intervention
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .68, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //9 Alcohol Individual Short
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .34, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //10 Alcohol Individual Long
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .62, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //11 Alcohol Group Short
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .47, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //12 Alcohol Group Long
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //13 Alcohol Community
	{ 1, 1, 1, .84, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //14 Depression Individual Short
	{ 1, 1, 1, .62, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //15 Depression Individual Long
	{ 1, 1, 1, .81, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //16 Depression Group Short
	{ 1, 1, 1, .71, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //17 Depression Group Long
	{ .97, 1, 1, 1, 1, .84, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //18 Sex Individual Short
	{ .92, 1, 1, 1, 1, .64, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //19 Sex Individual Long
	{ .97, 1, 1, 1, 1, .81, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //20 Sex Group Short
	{ .94, 1, 1, 1, 1, .71, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //21 Sex Group Long
	{ .97, 1, 1, 1, 1, .84, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //22 Sex Community
	{ 1, 1, 1, .67, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //23 Adherence SMS
	{ 1, 1, 1, .75, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //24 Brief adherence counseling
	{ .91, .93, 1, 1, 1, 1, 1, 1, .89, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, .89 } //25 Mass media HIV education campaigns
	
#else
	/*Pathways in comments
	0		1		2		3  4  5		6  7	8		9  10 11				12							13  */
	{ .51, 1, 1, 1, 1, .31, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, .89 },	//0     Condoms, risk reduction conseling, and STI treatment for female CSW
	{ 1, .59, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_10000, 1 },	//1     Expaned VCT and universal ART
	{ 1, .24, 1, 1, 1, 1, 1, 1, 1, 1, 1, .8, CD4_TREAT_THRESH_DEFAULT, 1 },	//2     PMTCT
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .55, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//3     Alcohol/SBIRT
	{ 1, 1, 1, 1, 1, 1, 1, .4, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//4     Circumcision
	{ .91, .93, 1, 1, 1, 1, 1, 1, .89, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//5     Mass media HIV education campaigns     
	{ .81, 1, 1, 1, 1, .69, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//6     Structural condom distribution
	{ 1, 1, 1, 1, .5, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //7 PrEP
	{ CONDOM_PATH_EFFECT, 1, LINKAGE_PATH_EFFECT, ADHERENCE_PATH_EFFECT, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 },	//8     Care+ intervention
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .475, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //9 Alcohol Individual 
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //10 Alcohol Individual Long
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .8, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //11 Alcohol Group 
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //12 Alcohol Group Long
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .85, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //13 Alcohol Community
	{ 1, 1, 1, .65, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //14 Depression Individual 
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //15 Depression Individual Long
	{ 1, 1, 1, .87, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //16 Depression Group 
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //17 Depression Group Long
	{ .55, 1, 1, 1, 1, .55, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //18 Sex Individual 
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //19 Sex Individual Long
	{ .8, 1, 1, 1, 1, .8, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //20 Sex Group 
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //21 Sex Group Long
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //22 Sex Community
	{ 1, 1, 1, .475, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //23 Adherence Individual
	{ 1, 1, 1, .74, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 }, //24 Adherence Group 
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, PMTCT_DEFAULT, CD4_TREAT_THRESH_DEFAULT, 1 } //25 Mass media HIV education campaigns
#endif 
};


#if INTERVENTION_EFFECT_SIZE_SENSITIVITY
string high_low;
double interv_effect_size_low[NUM_INTERVENTIONS][NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE]= 
{
	/*Pathways in comments  
     0		1		2		3  4  5		6  7	8		9  10 11				12                  13  */
	{0.1,	1,		1,		1, 1, .01,	1, 1,	1,		1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//0
	{1,		.5,		1,		1, 1, 1,	1, 1,	1,		1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_10000,    1}, //1
	{1,		.1,		1,		1, 1, 1,	1, 1,	1,		1, 1, .9,               CD4_TREAT_THRESH_DEFAULT,  1},	//2
	{1,		1,		1,		1, 1, 1,	1, 1,	1,		1, .5,PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//3
	{1,		1,		1,		1, 1, 1,	1, .38,	1,		1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//4
	{0.67,	.8,		1,		1, 1, 1,	1, 1,	.8,		1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//5
	{.69,	1,      1,		1, 1, .52,	1, 1,	1,      1, 1, PMTCT_DEFAULT,	CD4_TREAT_THRESH_DEFAULT,  1   }	//6     Structural condom distribution

};

double interv_effect_size_high[NUM_INTERVENTIONS][NUM_PATHWAYS_INCLUDING_MOVING_PEOPLE]= 
{
	/*Pathways in comments  
     0		1		2		3  4  5		6  7	8		9  10 11				12                  13  */
	{0.9,	1,		1,		1, 1, .51,	1, 1,	1,		1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//0
	{1,		.68,	1,		1, 1, 1,	1, 1,	1,		1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_10000,    1}, //1
	{1,		.5,		1,		1, 1, 1,	1, 1,	1,		1, 1, .7,               CD4_TREAT_THRESH_DEFAULT,  1},	//2
	{1,		1,		1,		1, 1, 1,	1, 1,	1,		1, .9,PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//3
	{1,		1,		1,		1, 1, 1,	1, .45,	1,		1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//4
	{1.0,	1.0,	1,		1, 1, 1,	1, 1,	1.0,	1, 1, PMTCT_DEFAULT,    CD4_TREAT_THRESH_DEFAULT,  1},	//5
	{.94,	1,      1,		1, 1, .91,	1, 1,	1,      1, 1, PMTCT_DEFAULT,	CD4_TREAT_THRESH_DEFAULT,  1   }	//6     Structural condom distribution
};

#endif

//the alpha multiplier for each pathway (if 1, then no effect on alpha)
double alpha_sex_mult[NUM_PATHWAYS] = {.2, 1, 1, 1, .56, .6, .61, .41};


//These values will get overwritten during the calibration period.  We need to store them so that we can
//ramp up to this value over the length of the calibration period.
double prob_not_being_tested_for_hiv_end_of_calibration = baselineRR_raw[1][0];
double prob_failing_to_be_linked_to_care_once_diagnosed_end_of_calibration = baselineRR_raw[2][0];
double prob_not_being_tested_for_hiv_start_of_calibration = PROB_NOT_BEING_TESTED_FOR_HIV_START_OF_CALIBRATION;  //makes sense to put this in the same place


// Enables us to mix types of costs when adding up the total basket cost - prepurchased or utilization for each intervention (Optimising Prevention Portfolio of Prevention)
bool costType[NUM_INTERVENTIONS] = { PP, UT, PP, UT, UT, PP, UT, PP, PP, /*new India values start here*/ PP, PP, PP, PP, PP, PP, PP, PP, PP, PP, PP, PP, PP, PP, PP };
