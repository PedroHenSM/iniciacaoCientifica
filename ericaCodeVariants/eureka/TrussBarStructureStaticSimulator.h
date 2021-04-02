/*
 * File:   TrussBarStructureStaticSimulator.h
 * Author: Heder S. Bernardino
 *
 * Created on September 2, 2011, 1:49 PM
 */

#ifndef TRUSSBARSTRUCTURESTATICSIMULATOR_H
#define	TRUSSBARSTRUCTURESTATICSIMULATOR_H

#include <fstream>
#include <tr1/memory>
#include <stdlib.h>

using namespace std;
using namespace std::tr1;

namespace truss_structure_simulator {

	class TrussBarStructureStaticSimulator {
    public:
		/* Constructor.
		 */
        TrussBarStructureStaticSimulator(double gamma, int numnp_cpp, int numeg_cpp, int nlcase_cpp, int modex_cpp, int* node_n_cpp, int** node_id_cpp, double** node_position_cpp, int* ll_cpp, int* loads_n_cpp, int** loads_node_cpp, int** loads_direction_cpp, double** loads_cpp, int elements_npar_1, int elements_npar_2, int elements_npar_3, int* element_id_cpp, double* element_cpp, double* element_extra_cpp, int* m_cpp, int* ii_cpp, int* jj_cpp, int* mtyp_cpp, int* kg_cpp);
		/* Constructor.
		 */
        TrussBarStructureStaticSimulator(const TrussBarStructureStaticSimulator& orig);
		/* Destructor.
		 */
        virtual ~TrussBarStructureStaticSimulator();
		/* Evaluate the structure.
		 * param - vector with the design variables.
		 * fobj - variable in which the weight of the structure will be written.
		 * v - vector in which the constraint violations will be written.
		 */
		double stap(double *param, double *fobj, double *v);
		/*
		 * Constant which adjust the amount of memory must be allocated.
		 * TODO -- make it adjust automatically
		 */
		static const int MEMORY_SPACE = 400000;
    private:
		/*
		 * Specific weight of the material.
		 */
		double gamma;
		//variables with data from input file
		int numnp_cpp;
		int numeg_cpp;
		int nlcase_cpp;
		int modex_cpp;
		int* node_n_cpp;
		int** node_id_cpp;
		double** node_position_cpp;
		int* ll_cpp;
		int* loads_n_cpp;
		int** loads_node_cpp;
		int** loads_direction_cpp;
		double** loads_cpp;
		int elements_npar_1;
		int elements_npar_2;
		int elements_npar_3;
		int* element_id_cpp;
		double* element_cpp;
		double* element_extra_cpp; //this information is not used by stap
		int* m_cpp;
		int* ii_cpp;
		int* jj_cpp;
		int* mtyp_cpp;
		int* kg_cpp;
		//replace read-write files un1 and un2 (ielmnt and iload)
		int* ielmnt_info1;
		int* ielmnt_info2;
		double** ielmnt;
		double** iload;
		int i_iload;
		//variables with data from input file -- end
		//Variables
		int* npar;
		int numnp, neq, nwk, numest, midest, maxest, mk;
		int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15;
		int n101, n102, n103, n104, n105, n106, n107;
		int ind, numeg;
		int mtot; // used to verify the limit of the 'memory' variable -- candidate to be removed
		int nfirst, nlast;
		int itwo; //responsible to the precision used (2) -- candidate to be removed
		int ng, modex;
		double* memory;
		//Methods
		int tensoes(double *tens1, double *tens2, double *tens3, int *l, double *tens);
		int desloc(double *d1, double *d2, double * d3, int *l, double *desl, int n_desloc_cpp);
		int endarea(double *area, int *nume, double * param);
		int colht(int *mht, int *nd, int *lm);
		int addban(double *a, int *maxa, double *s, int *lm, int *nd);
		int russ(int *id, double *x, double *y, double *z__, double *u, int *mht, double *e, double *area, int *lm, double *xyz, int *matp, double *tens, int *ncar);
		int error(int *n, int *i__);
		int truss(double *param, int *ncar);
		int elemnt(double *param, int *ncar);
		int stress(double *aa, double *param, int *ncar);
		int writed(double *disp, int *id, int *neq, int *numnp);
		int loadv(double *r__, int *neq);
		int colsol(double *a, double *v, int *maxa, int *nn, int *nwk, int *nnm, int *kkk);
		int assem(double *aa, double *param, int *ncar);
		int clear(double *a, int *n);
		int addres(int *maxa, int *mht);
		int elcal(double *param, int *ncar);
		int loads(double *r__, int *nod, int *idirn, double *fload, int *id, int *nload, int *neq, int loadCase);
		int input(int *id, double *x, double *y, double *z__, int *numnp, int *neq);
		int calcfobj(double *u, double *tens, double *xyz, double *area, int *nume, double *param, double *fobj, double *v, double *d1, double * d2, double *d3, double *tens1, double *tens2, double * tens3, int n_desloc_cpp, int nlcase);
    };

	// typedef shared_ptr<TrussBarStructureStaticSimulator> TrussBarStructureStaticSimulatorPtr; // TODO: MODIFIED. 'shared_ptr' is ambiguous
    typedef std::shared_ptr<TrussBarStructureStaticSimulator> TrussBarStructureStaticSimulatorPtr;


}

#endif	/* TRUSSBARSTRUCTURESTATICSIMULATOR_H */

