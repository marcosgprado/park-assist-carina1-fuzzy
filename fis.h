// fis.h

#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstring>
#include <cstdlib>

using namespace std;

#ifndef ABS
# define ABS(x)   ( (x) > (0) ? (x): (-(x)) )
#endif
#ifndef MAX
# define MAX(x,y) ( (x) > (y) ? (x) : (y) )
#endif
#ifndef MIN
# define MIN(x,y) ( (x) < (y) ? (x) : (y) )
#endif
#define MF_PARA_N 4
#define STR_LEN 500
#define MF_POINT_N 101

#if (defined(MATLAB_MEX_FILE) && !defined(__SIMSTRUC__))
# define FREE mxFree
#else
# define FREE free
#endif

#define FREEMAT(mat,m) fisFreeMatrix(mat,m)
#define FREEARRAY(array) FREE(array)

void fisError(const char *msg);
void *fisCalloc(int num_of_x, int size_of_x);
void fisGetString2(char *target, double *array, int max_leng);
int compareString(const char *string1, char *string2);
FILE *fisOpenFile(char *file, const char *mode);

class Mf_Node
{
public:
	// atributos
	char label[STR_LEN];	/* MF name */
	char type[STR_LEN];		/* MF type */
	int nparams;			/* length of params field */
	double *params;			/* MF parameters */
	int userDefined;		/* 1 if the MF is user-defined */
	double (Mf_Node::*mfFcn)(double, double *); /* pointer to a mem. fcn */ 
	double value;		    /* for Sugeno only */
	double *value_array;	/* for Mamdani only, array of MF values */	
	
	// metodos
	double fisTriangleMf(double x, double *params);
	double fisTrapezoidMf(double x, double *params);
	double fisGaussianMf(double x, double *params);
	double fisGaussian2Mf(double x, double *params);
	double fisSigmoidMf(double x, double *params);
	double fisProductSigmoidMf(double x, double *params);
	double fisDifferenceSigmoidMf(double x, double *params);
	double fisGeneralizedBellMf(double x, double *params);
	double fisSMf(double x, double *params);
	double fisZMf(double x, double *params);
	double fisPiMf(double x, double *params);
	double fisAllMf(double x, double *params);

};

class Io_Node
{
	public:
	char name[STR_LEN];
	int mf_n;
	double bound[2];
	double value;
	Mf_Node **mf;

	
		
};


class Fis_Node{
// atributos
public:
	int handle;
	int load_param;
	char name[STR_LEN];
	char type[STR_LEN];
	char andMethod[STR_LEN];
	char orMethod[STR_LEN];
	char impMethod[STR_LEN];
	char aggMethod[STR_LEN];
	char defuzzMethod[STR_LEN];
	int userDefinedAnd;
	int userDefinedOr;
	int userDefinedImp;
	int userDefinedAgg;
	int userDefinedDefuzz;
	int in_n;
	int out_n;
	int rule_n;
	int **rule_list;
	double *rule_weight;
	int *and_or;	/* AND-OR indicator */
	double *firing_strength;
	double *rule_output;
	/* Sugeno: output for each rules */
	/* Mamdani: constrained output MF values of rules */
	
	double (Fis_Node::*andFcn)(double, double);
	double (Fis_Node::*orFcn)(double, double);
	double (Fis_Node::*impFcn)(double, double);
	double (Fis_Node::*aggFcn)(double, double);
	double (Fis_Node::*defuzzFcn)(void *fis, int m, double *mf, int numofpoints); // inserido os quatro parâmetros por diogo 07 junho 2010
	double *BigOutMfMatrix;	/* used for Mamdani system only */
	double *BigWeightMatrix;/* used for Mamdani system only */
	double *mfs_of_rule;	/* MF values in a rule */

	double *bias; /*bias, to be tuned when no rules are fired*/
	int isbias;

	Fis_Node *next;
	
	Io_Node **input;
	Io_Node **output;

// métodos
//	void fisError(const char *msg);
//	FILE *fisOpenFile(char *file, const char *mode);
//	void *fisCalloc(int num_of_x, int size_of_x);
	char **fisCreateMatrix(int row_n, int col_n, int element_size);
	void fisFreeMatrix(void **matrix, int row_n);
	double**fisCopyMatrix(double **source, int row_n, int col_n);
	void fisPrintMatrix(double **matrix, int row_n, int col_n);
	void fisPrintArray(double *array, int size);
	void fisPause();
	

	
	double fisMin(double x, double y);
	double fisMax(double x, double y);
	double fisProduct(double x, double y);
	double fisProbOr(double x, double y);
	double fisSum(double x, double y);
	double fisArrayOperation(double *array, int size, double (Fis_Node::*fcn)(double, double));
	double defuzzCentroid(void *fis1, int m, double *mf, int numofpoints);
	double defuzzBisector(void *fis1, int m, double *mf, int numofpoints);
	double defuzzMeanOfMax(void *fis1, int m, double *mf, int numofpoints);
	double defuzzSmallestOfMax(void *fis1, int m, double *mf, int numofpoints);
	double defuzzLargestOfMax(void *fis1, int m, double *mf, int numofpoints);
	
	void fisAssignMfPointer(Fis_Node *fis);
	void fisAssignFunctionPointer(Fis_Node *fis);
	
	void fisPrintData(Fis_Node *fis);

	void fisFreeIoList(Io_Node *io_list, int n);
	void fisFreeMfList(Mf_Node *mf_list, int n);
	void fisFreeFisNode(Fis_Node *fis);
	
	void fisComputeOutputMfValueArray(Fis_Node *fis, int numofpoints);
	
	void fisCheckDataStructure(Fis_Node *fis);
	Io_Node *fisBuildIoList(int node_n, int *mf_n);
	void fisBuildFisNode(Fis_Node *fis, double **fismatrix, int col_n, int numofpoints);
	
	void fisLoadParameter(Fis_Node *fis, double **fismatrix, int numofpoints);
	void fisLoadParameter1(Fis_Node *fis, double *para_array, int numofpoints);
	Fis_Node *fisMatchHandle(Fis_Node *head, int handle);
	Fis_Node *fisMatchName(Fis_Node *head, char *name);
	int fisFindMaxHandle(Fis_Node *head);
	int fisGetMfParaN(char *mfType);
	
	void fisComputeInputMfValue(Fis_Node *fis);
	void fisComputeTskRuleOutput(Fis_Node *fis);
	void fisComputeFiringStrength(Fis_Node *fis);
	double fisFinalOutputMf(Fis_Node *fis, int m, int n);
	void fisFinalOutputMf2(Fis_Node *fis, int m, double *aggMF, int numofpoints);
	void fisEvaluate(Fis_Node *fis, int numofpoints);
	void getFisOutput(double *input, Fis_Node *fis, double *output);
	
	double ** returnEmptyFismatrix(char *filename, int *row_n_p, int *col_n_p);
	double ** returnFismatrix(char *fis_file, int *row_n_p, int *col_n_p);
	double ** returnDataMatrix(char *filename, int *row_n_p, int *col_n_p);

	
		
};





