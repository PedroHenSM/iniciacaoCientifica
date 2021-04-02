/* 
 * File:   F107Truss72Bar.cpp
 * Author: Heder S. Bernardino
 * 
 * Created on January 22, 2011, 6:33 PM
 */

#include <cmath>
#include <memory.h>
#include "F107Truss72Bar.h"

using namespace std;

namespace problem {

    // F107Truss72Bar::F107Truss72Bar(): TrussBarStructureStaticProblem(16, NULL, 35000, 1, 168, 72, 0.1, 25000, 0.25, "input72.dat", 0.1, 2.5), grouping(NULL) { // original line
    // F107Truss72Bar::F107Truss72Bar(): TrussBarStructureStaticProblem(16, NULL, 35000, 1, 168, 72, 0.1, 25000, 0.25, "../functions/eureka/input_data/input72.dat", 0.1, 2.5), grouping(NULL) {
    F107Truss72Bar::F107Truss72Bar(): TrussBarStructureStaticProblem(16, NULL, 35000, 1, 168, 72, 0.1, 25000, 0.25, "ins/input_data/input72.dat", 0.1, 2.5), grouping(NULL) {
        
    }

    F107Truss72Bar::F107Truss72Bar(const F107Truss72Bar& orig): TrussBarStructureStaticProblem(orig), grouping(NULL) {
        
    }

    F107Truss72Bar::~F107Truss72Bar() {
        delete[] this->grouping;
    }
	
	int* const F107Truss72Bar::getGrouping () {
		if (this->grouping==NULL) {
			this->grouping = new int[ this->getNumberOfBars() ];
			this->grouping[0] = 4;
			this->grouping[1] = 12;
			this->grouping[2] = 16;
			this->grouping[3] = 18;
			this->grouping[4] = 22;
			this->grouping[5] = 30;
			this->grouping[6] = 34;
			this->grouping[7] = 36;
			this->grouping[8] = 40;
			this->grouping[9] = 48;
			this->grouping[10] = 52;
			this->grouping[11] = 54;
			this->grouping[12] = 58;
			this->grouping[13] = 66;
			this->grouping[14] = 70;
			this->grouping[15] = 72;
		}
		return this->grouping;
	}
	
	void F107Truss72Bar::evaluation(void* vector, void* values) {
        
		TrussBarStructureStaticProblem::evaluation(vector, values);
		
        double* val = (double*) values;
        
		int i=this->getNLCase()*this->getNumberOfBars();
		int max = i+12;
        
		for(; i< max; i++) {
            val[i+1+12] = abs( ( this->getStressDisplacementAux()[ i+this->getNumberOfDisplacements() ] ) ) / this->getDisplacementConstraint() - 1;
        }
		
    }

}
