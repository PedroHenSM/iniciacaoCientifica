/* 
 * File:   F103Truss25Bar.cpp
 * Author: hedersb
 * 
 * Created on January 21, 2011, 11:15 AM
 */

#include <cmath>
#include <iostream>
#include <memory.h>
#include <stdio.h>
#include "F103Truss25Bar.h"

using namespace std;

namespace problem {

    F103Truss25Bar::F103Truss25Bar(): TrussBarStructureStaticProblem(8, NULL, 240000, 1, 29, 25,  0.1, 40000, 0.35, "input25.dat", 0.1, 3.4), grouping(NULL) {
        
    }

    F103Truss25Bar::F103Truss25Bar(const F103Truss25Bar& orig): TrussBarStructureStaticProblem(orig), grouping(NULL) {
        
    }

    F103Truss25Bar::~F103Truss25Bar() {
        delete[] this->grouping;
    }
	
	int* const F103Truss25Bar::getGrouping () {
		
		if (this->grouping==NULL) {
			this->grouping = new int[ this->getDimension() ];
			this->grouping[0] = 1;
			this->grouping[1] = 5;
			this->grouping[2] = 9;
			this->grouping[3] = 11;
			this->grouping[4] = 13;
			this->grouping[5] = 17;
			this->grouping[6] = 21;
			this->grouping[7] = 25;
		}
		
		return this->grouping;
		
	}
	
	void F103Truss25Bar::evaluation(void* vector, void* values) {
		
		TrussBarStructureStaticProblem::evaluation(vector, values);
		
		double* val = (double*) values;
		
		int i = this->getNLCase()*this->getNumberOfBars();
		val[i+1+2] = abs(this->getStressDisplacementAux()[i+3])/this->getDisplacementConstraint() - 1;
		val[i+1+3] = abs(this->getStressDisplacementAux()[i+4])/this->getDisplacementConstraint() - 1;
		
	}

}
