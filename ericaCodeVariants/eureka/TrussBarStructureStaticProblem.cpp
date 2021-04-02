/* 
 * File:   TrussBarStructureStaticProblem.cpp
 * Author: Heder S. Bernardino
 * 
 * Created on September 11, 2011, 10:11 AM
 */

#include "TrussBarStructureStaticProblem.h"
#include <iostream>
#include <string.h>
#include <cmath>

using namespace std;

namespace problem {

	TrussBarStructureStaticProblem::TrussBarStructureStaticProblem(int dimension, void* bounds, int maxNumberObjectiveFunctionEvaluations, int numberObjectives, int numberConstraints, int numberOfBars, double gamma, double stressConstraint, double displacementConstraint, /*int numberOfDisplacementConstraints,*/ string inputFileName, double lowerBound, double upperBound): Problem(dimension, bounds, maxNumberObjectiveFunctionEvaluations, numberObjectives, numberConstraints), numberOfBars(numberOfBars), gamma(gamma), stressConstraint(stressConstraint), displacementConstraint(displacementConstraint)/*, numberOfDisplacementConstraints(numberOfDisplacementConstraints)*/, areasAux(new double[ this->numberOfBars ]), numberOfDisplacements(-1) {
		double** b = new double*[this->dimension];
        for (int i = 0; i < dimension; i++) {
            b[i] = new double[2];
            b[i][0] = lowerBound;
            b[i][1] = upperBound;
        }
        this->bounds = b;
		
		this->readFile(inputFileName);
		this->numberOfDisplacements = this->getNumberOfDisplacements();
		this->stressDisplacementAux = new double[ (this->nlcase_cpp*(this->numberOfBars + this->numberOfDisplacements)) ];
		this->simulator = TrussBarStructureStaticSimulatorPtr(
			new TrussBarStructureStaticSimulator(
				this->gamma, this->numnp_cpp, this->numeg_cpp, this->nlcase_cpp, this->modex_cpp, this->node_n_cpp, this->node_id_cpp, node_position_cpp, this->ll_cpp, this->loads_n_cpp, this->loads_node_cpp, this->loads_direction_cpp, this->loads_cpp, this->elements_npar_1, this->elements_npar_2, this->elements_npar_3, this->element_id_cpp, this->element_cpp, this->element_extra_cpp, this->m_cpp, this->ii_cpp, this->jj_cpp, this->mtyp_cpp, this->kg_cpp
			)
		);
	}

	TrussBarStructureStaticProblem::TrussBarStructureStaticProblem(const TrussBarStructureStaticProblem& orig): Problem(orig), areasAux(new double[ this->numberOfBars ]), numberOfDisplacements(-1) {
		double** bOrig = (double**) orig.bounds;
        double** b = new double*[this->dimension];
        for (int i = 0; i < dimension; i++) {
            b[i] = new double[2];
            memcpy(b[i], bOrig[i], 2 * sizeof (double));
        }
        this->bounds = b;
	}

	TrussBarStructureStaticProblem::~TrussBarStructureStaticProblem() {
		//bounds
		double** b = (double**) this->bounds;
        for (int i = 0; i < this->dimension; i++) {
            delete[] b[i];
        }
        delete[] b;
		//aux 
		delete[] this->stressDisplacementAux;
		delete[] this->areasAux;
		//input data
		for(int i=0; i<this->numnp_cpp; i++) {
			delete[] this->node_id_cpp[i];
			delete[] this->node_position_cpp[i];
		}
		delete[] this->node_id_cpp;
		delete[] this->node_position_cpp;
		delete[] this->node_n_cpp;
		for(int i=0; i<this->nlcase_cpp; i++) {
			delete[] this->loads_node_cpp[i];
			delete[] this->loads_direction_cpp[i];
			delete[] this->loads_cpp[i];
		}
		delete[] this->ll_cpp;
		delete[] this->loads_n_cpp;
		delete[] this->loads_node_cpp;
		delete[] this->loads_direction_cpp;
		delete[] this->loads_cpp;
		delete[] this->element_id_cpp;
		delete[] this->element_cpp;
		delete[] this->element_extra_cpp;
		delete[] this->m_cpp;
		delete[] this->ii_cpp;
		delete[] this->jj_cpp;
		delete[] this->mtyp_cpp;
		delete[] this->kg_cpp;
	}
	
	void TrussBarStructureStaticProblem::evaluation(void* vector, void* values) {
        double* x = (double*) vector;
		//double x[] = {2.147422, 2.649941, 1.399858, 2.257593, 3.441807, 1.851235, 2.404124, 2.930071, 2.297117, 4.236111, 1.862823, 0.545913, 2.976073, 2.877205, 3.211396, 0.787109, 2.788300, 1.893071, 1.359410, 2.437747, 2.818608, 2.468314, 1.056329, 1.772348, 2.319242};
		this->fillAreasAux(x);
		
        double* val = (double*) values;
        double objectiveFunctionValue = 0;
        
        this->simulator->stap(this->areasAux, &objectiveFunctionValue, this->stressDisplacementAux);
        val[0] = objectiveFunctionValue;
		int i=0;
        for(; i< this->nlcase_cpp*this->numberOfBars; i++) {
            val[i+1] = fabs( this->stressDisplacementAux[i] )/this->stressConstraint - 1;
        }
		for(; i< this->numberConstraints; i++) {
            val[i+1] = fabs( this->stressDisplacementAux[i] )/this->displacementConstraint - 1;
        }
		
		/*int max = (this->nlcase_cpp*(this->numberOfBars + this->numberOfDisplacements));
		int base = this->nlcase_cpp*this->numberOfBars;
		int base2 = this->numberOfDisplacements;
		for(i = base; i< base+base2; i++) {
            cout << i << "\t" << (abs(this->stressDisplacementAux[i])/this->displacementConstraint - 1) << endl;
        }
		max = 0;*/
		
    }
	
	void TrussBarStructureStaticProblem::fillAreasAux(double* x) {
		
		int* grouping = this->getGrouping();
		
		if (grouping==NULL) {
			for(int i = 0; i < this->getDimension(); i++) {
				this->areasAux[i] = x[i];
			}
		} else {
			int j = 0;
			for(int i = 0; i < this->getDimension(); i++) {
				for(; j< grouping[i]; j++) {
					this->areasAux[j] = x[i];
				}
			}
		}
		
	}
	
	int TrussBarStructureStaticProblem::getNumberOfDisplacements() {
		if (this->numberOfDisplacements < 0) {
			this->numberOfDisplacements = 0;
			for(int i=0; i < this->numnp_cpp; i++) {
				for (int j = 0; j < 3; j++) {
					if (this->node_id_cpp[i][j] == 0) {
						numberOfDisplacements++;
					}
				}
			}
		}
		return numberOfDisplacements;
	}
	
	void TrussBarStructureStaticProblem::readFile(string fileName) {

		
		ifstream inputFile( fileName.c_str());
		string line;
		char* lineC;
		int i;
		if (inputFile.is_open()) {
			if (inputFile.good()) {
				getline(inputFile, line); //first line only presents information about the problem
				//first information
				getline(inputFile, line); 
				lineC = new char[line.length() + 1];
				strcpy(lineC, line.c_str());
				this->numnp_cpp = atoi(strtok(lineC, " "));
				this->numeg_cpp = atoi(strtok(NULL, " "));
				this->nlcase_cpp = atoi(strtok(NULL, " "));
				this->modex_cpp = atoi(strtok(NULL, " "));
				delete[] lineC;
				//cout << this->numnp_cpp << "\t" << this->numeg_cpp << "\t" << this->nlcase_cpp << "\t" << this->modex_cpp << endl;
				
				//nodes
				this->node_n_cpp = new int[this->numnp_cpp];
				this->node_id_cpp = new int*[this->numnp_cpp];
				this->node_position_cpp = new double*[this->numnp_cpp];
				i = 0;
				while (inputFile.good() && i < this->numnp_cpp) {
					this->node_id_cpp[i] = new int[3];
					this->node_position_cpp[i] = new double[3];
					getline(inputFile, line);
					lineC = new char[line.length() + 1];
					strcpy(lineC, line.c_str());
					this->node_n_cpp[i] = atoi(strtok(lineC, " "));
					this->node_id_cpp[i][0] = atoi(strtok(NULL, " "));
					this->node_id_cpp[i][1] = atoi(strtok(NULL, " "));
					this->node_id_cpp[i][2] = atoi(strtok(NULL, " "));
					this->node_position_cpp[i][0] = atof(strtok(NULL, " "));
					this->node_position_cpp[i][1] = atof(strtok(NULL, " "));
					this->node_position_cpp[i][2] = atof(strtok(NULL, " "));
					delete[] lineC;
					//cout << this->node_n_cpp[i] << "\t" << this->node_id_cpp[i][0] << "\t" << this->node_id_cpp[i][1] << "\t" << this->node_id_cpp[i][2] << "\t" << this->node_position_cpp[i][0] << "\t" << this->node_position_cpp[i][1] << "\t" << this->node_position_cpp[i][2] << endl;
					i++;
				}
				
				//loads
				this->ll_cpp = new int[ this->nlcase_cpp ];
				this->loads_n_cpp = new int[ this->nlcase_cpp ];
				this->loads_node_cpp = new int*[ this->nlcase_cpp ];
				this->loads_direction_cpp = new int*[ this->nlcase_cpp ];
				this->loads_cpp = new double*[ this->nlcase_cpp ];
				
				for(int load=0; load < this->nlcase_cpp; load++) {
					
					if (inputFile.good()) {
						getline(inputFile, line); 
						lineC = new char[line.length() + 1];
						strcpy(lineC, line.c_str());
						this->ll_cpp[load] = atoi(strtok(lineC, " "));
						this->loads_n_cpp[load] = atoi(strtok(NULL, " "));
						delete[] lineC;
					} else {
						cout << "Some problem with input file -- information about the loads." << endl;
						exit( EXIT_FAILURE );
					}
					this->loads_node_cpp[load] = new int[this->loads_n_cpp[load]];
					this->loads_direction_cpp[load] = new int[this->loads_n_cpp[load]];
					this->loads_cpp[load] = new double[this->loads_n_cpp[load]];
					i = 0;
					while (inputFile.good() && i < this->loads_n_cpp[load]) {
						getline(inputFile, line);
						lineC = new char[line.length() + 1];
						strcpy(lineC, line.c_str());
						this->loads_node_cpp[load][i] = atoi(strtok(lineC, " "));
						this->loads_direction_cpp[load][i] = atoi(strtok(NULL, " "));
						this->loads_cpp[load][i] = atof(strtok(NULL, " "));
						delete[] lineC;
						//cout << this->loads_node_cpp[i] << "\t" << this->loads_direction_cpp[i] << "\t" << this->loads_cpp[i] << endl;
						i++;
					}
				}
				
				//elements
				if (inputFile.good()) {
					getline(inputFile, line); 
					lineC = new char[line.length() + 1];
					strcpy(lineC, line.c_str());
					this->elements_npar_1 = atoi(strtok(lineC, " "));
					this->elements_npar_2 = atoi(strtok(NULL, " "));
					this->elements_npar_3 = atoi(strtok(NULL, " "));
					delete[] lineC;
				} else {
					cout << "Some problem with input file -- information about the elements." << endl;
					exit(1);
				}
				element_id_cpp = new int[this->elements_npar_2];
				element_cpp = new double[this->elements_npar_2];
				element_extra_cpp = new double[this->elements_npar_2];
				i = 0;
				while (inputFile.good() && i < this->elements_npar_2) {
					getline(inputFile, line);
					lineC = new char[line.length() + 1];
					strcpy(lineC, line.c_str());
					this->element_id_cpp[i] = atoi(strtok(lineC, " "));
					this->element_cpp[i] = atof(strtok(NULL, " "));
					this->element_extra_cpp[i] = atof(strtok(NULL, " "));
					delete[] lineC;
					//cout << this->element_id_cpp[i] << "\t" << this->element_cpp[i] << "\t" << this->element_extra_cpp[i] << endl;
					i++;
				}
				this->m_cpp = new int[this->elements_npar_2];
				this->ii_cpp = new int[this->elements_npar_2];
				this->jj_cpp = new int[this->elements_npar_2];
				this->mtyp_cpp = new int[this->elements_npar_2];
				this->kg_cpp = new int[this->elements_npar_2];
				i = 0;
				while (inputFile.good() && i < this->elements_npar_2) {
					getline(inputFile, line);
					lineC = new char[line.length() + 1];
					strcpy(lineC, line.c_str());
					this->m_cpp[i] = atoi(strtok(lineC, " "));
					this->ii_cpp[i] = atoi(strtok(NULL, " "));
					this->jj_cpp[i] = atoi(strtok(NULL, " "));
					this->mtyp_cpp[i] = atoi(strtok(NULL, " "));
					this->kg_cpp[i] = atoi(strtok(NULL, " "));
					delete[] lineC;
					//cout << this->m_cpp[i] << "\t" << this->ii_cpp[i] << "\t" << this->jj_cpp[i] << "\t" << this->mtyp_cpp[i] << "\t" << this->kg_cpp[i] << endl;
					i++;
				}
			}
		}
		inputFile.close();

	}
	
	int const TrussBarStructureStaticProblem::getNumberOfBars() {
		return this->numberOfBars;
	}
	
	int* const TrussBarStructureStaticProblem::getGrouping() {
		return NULL;
	}
	
	double* const TrussBarStructureStaticProblem::getAreasAux() {
		return this->areasAux;
	}
	
	int const TrussBarStructureStaticProblem::getNLCase() {
		return this->nlcase_cpp;
	}
	
	double const TrussBarStructureStaticProblem::getDisplacementConstraint() {
		return this->displacementConstraint;
	}
	
	double* const TrussBarStructureStaticProblem::getStressDisplacementAux() {
		return this->stressDisplacementAux;
	}
	
	TrussBarStructureStaticSimulatorPtr const TrussBarStructureStaticProblem::getSimulator() {
		return this->simulator;
	}
	

}
