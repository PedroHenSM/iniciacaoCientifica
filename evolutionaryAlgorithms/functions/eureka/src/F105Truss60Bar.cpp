/* 
 * File:   F105Truss60Bar.cpp
 * Author: Heder S. Bernardino
 * 
 * Created on January 22, 2011, 6:32 PM
 */

#include <cmath>
#include <iostream>
#include <memory.h>
#include "F105Truss60Bar.h"


using namespace std;

namespace problem {

    F105Truss60Bar::F105Truss60Bar(): TrussBarStructureStaticProblem(25, NULL, 150000, 1, 198, 60, 0.1, 10000, 1.75, "input60.dat", 0.5, 5), grouping(NULL) {
		this->displacementConstraint2 = 2.25;
		this->displacementConstraint3 = 2.75;
    }

    F105Truss60Bar::F105Truss60Bar(const F105Truss60Bar& orig): TrussBarStructureStaticProblem(orig), grouping(NULL) {
        
    }

    F105Truss60Bar::~F105Truss60Bar() {
        delete[] this->grouping;
    }

    int* const F105Truss60Bar::getGrouping () {
		if (this->grouping==NULL) {
			this->grouping = new int[ this->getNumberOfBars() ];
			this->grouping[0] = 1;
			this->grouping[1] = 2;
			this->grouping[2] = 3;
			this->grouping[3] = 4;
			this->grouping[4] = 5;
			this->grouping[5] = 6;
			this->grouping[6] = 7;
			this->grouping[7] = 8;
			this->grouping[8] = 9;
			this->grouping[9] = 10;
			this->grouping[10] = 11;
			this->grouping[11] = 12;
			this->grouping[12] = 1;
			this->grouping[13] = 2;
			this->grouping[14] = 3;
			this->grouping[15] = 4;
			this->grouping[16] = 5;
			this->grouping[17] = 6;
			this->grouping[18] = 7;
			this->grouping[19] = 8;
			this->grouping[20] = 9;
			this->grouping[21] = 10;
			this->grouping[22] = 11;
			this->grouping[23] = 12;
			this->grouping[24] = 13;
			this->grouping[25] = 14;
			this->grouping[26] = 15;
			this->grouping[27] = 16;
			this->grouping[28] = 17;
			this->grouping[29] = 18;
			this->grouping[30] = 19;
			this->grouping[31] = 20;
			this->grouping[32] = 21;
			this->grouping[33] = 22;
			this->grouping[34] = 23;
			this->grouping[35] = 24;
			this->grouping[36] = 13;
			this->grouping[37] = 14;
			this->grouping[38] = 15;
			this->grouping[39] = 16;
			this->grouping[40] = 17;
			this->grouping[41] = 18;
			this->grouping[42] = 19;
			this->grouping[43] = 20;
			this->grouping[44] = 21;
			this->grouping[45] = 22;
			this->grouping[46] = 23;
			this->grouping[47] = 24;
			this->grouping[48] = 0;
			this->grouping[49] = 0;
			this->grouping[50] = 0;
			this->grouping[51] = 0;
			this->grouping[52] = 0;
			this->grouping[53] = 0;
			this->grouping[54] = 0;
			this->grouping[55] = 0;
			this->grouping[56] = 0;
			this->grouping[57] = 0;
			this->grouping[58] = 0;
			this->grouping[59] = 0;
		}
		return this->grouping;
	}
	
	void F105Truss60Bar::fillAreasAux(double* x) {
		
		double* areas = this->getAreasAux();
		int* grouping = this->getGrouping();
		
		for(int i = 0; i < this->getNumberOfBars(); i++) {
			areas[i] = x[ grouping[i] ];
		}
		
	}
	
	void F105Truss60Bar::evaluation(void* vector, void* values) {
		
		TrussBarStructureStaticProblem::evaluation(vector, values);
        double* val = (double*) values;
        
		int i = this->getNLCase()*this->getNumberOfBars();
		
		/*for(; i< this->numberConstraints; i++) {
            val[i+1] = abs(this->getStressDisplacementAux()[i])/this->getDisplacementConstraint() - 1;
        }*/
		
		
		/*int max = (this->getNLCase()*(this->getNumberOfBars() + this->getNumberOfDisplacements()));
		int base = this->getNLCase()*this->getNumberOfDisplacements();
		int base2 = this->getNLCase()*this->getNumberOfDisplacements();
		for(i = base; i< base+base2; i++) {
            cout << i << "\t" << (abs(this->getStressDisplacementAux()[i])/this->displacementConstraint3 - 1) << endl;
        }
		max = 0;*/
		
		
		
		
		int numberOfDisplacements = this->getNumberOfDisplacements();
		val[i+1] = abs(this->getStressDisplacementAux()[i+6])/this->getDisplacementConstraint() - 1; //only some displacements are constraints
		val[i+2] = abs(this->getStressDisplacementAux()[i+7])/this->getDisplacementConstraint() - 1;
		val[i+3] = abs(this->getStressDisplacementAux()[i+22])/this->displacementConstraint2 - 1;
		val[i+4] = abs(this->getStressDisplacementAux()[i+23])/this->displacementConstraint2 - 1;
		val[i+5] = abs(this->getStressDisplacementAux()[i+33])/this->displacementConstraint3 - 1;
		val[i+6] = abs(this->getStressDisplacementAux()[i+34])/this->displacementConstraint3 - 1;
		val[i+7] = abs(this->getStressDisplacementAux()[i+6+numberOfDisplacements])/this->getDisplacementConstraint() - 1; //only some displacements are constraints
		val[i+8] = abs(this->getStressDisplacementAux()[i+7+numberOfDisplacements])/this->getDisplacementConstraint() - 1;
		val[i+9] = abs(this->getStressDisplacementAux()[i+22+numberOfDisplacements])/this->displacementConstraint2 - 1;
		val[i+10] = abs(this->getStressDisplacementAux()[i+23+numberOfDisplacements])/this->displacementConstraint2 - 1;
		val[i+11] = abs(this->getStressDisplacementAux()[i+33+numberOfDisplacements])/this->displacementConstraint3 - 1;
		val[i+12] = abs(this->getStressDisplacementAux()[i+34+numberOfDisplacements])/this->displacementConstraint3 - 1;
		val[i+13] = abs(this->getStressDisplacementAux()[i+6+numberOfDisplacements+numberOfDisplacements])/this->getDisplacementConstraint() - 1; //only some displacements are constraints
		val[i+14] = abs(this->getStressDisplacementAux()[i+7+numberOfDisplacements+numberOfDisplacements])/this->getDisplacementConstraint() - 1;
		val[i+15] = abs(this->getStressDisplacementAux()[i+22+numberOfDisplacements+numberOfDisplacements])/this->displacementConstraint2 - 1;
		val[i+16] = abs(this->getStressDisplacementAux()[i+23+numberOfDisplacements+numberOfDisplacements])/this->displacementConstraint2 - 1;
		val[i+17] = abs(this->getStressDisplacementAux()[i+33+numberOfDisplacements+numberOfDisplacements])/this->displacementConstraint3 - 1;
		val[i+18] = abs(this->getStressDisplacementAux()[i+34+numberOfDisplacements+numberOfDisplacements])/this->displacementConstraint3 - 1;
		
		/*v[181] = (d__1 = d1[7] / ubarra1, abs(d__1)) - 1.;
		v[182] = (d__1 = d1[8] / ubarra1, abs(d__1)) - 1.;
		v[183] = (d__1 = d1[23] / ubarra2, abs(d__1)) - 1.;
		v[184] = (d__1 = d1[24] / ubarra2, abs(d__1)) - 1.;
		v[185] = (d__1 = d1[34] / ubarra3, abs(d__1)) - 1.;
		v[186] = (d__1 = d1[35] / ubarra3, abs(d__1)) - 1.;
		v[187] = (d__1 = d2[7] / ubarra1, abs(d__1)) - 1.; //45 is the number of displacements for each case
		v[188] = (d__1 = d2[8] / ubarra1, abs(d__1)) - 1.;
		v[189] = (d__1 = d2[23] / ubarra2, abs(d__1)) - 1.;
		v[190] = (d__1 = d2[24] / ubarra2, abs(d__1)) - 1.;
		v[191] = (d__1 = d2[34] / ubarra3, abs(d__1)) - 1.;
		v[192] = (d__1 = d2[35] / ubarra3, abs(d__1)) - 1.;
		v[193] = (d__1 = d3[7] / ubarra1, abs(d__1)) - 1.;
		v[194] = (d__1 = d3[8] / ubarra1, abs(d__1)) - 1.;
		v[195] = (d__1 = d3[23] / ubarra2, abs(d__1)) - 1.;
		v[196] = (d__1 = d3[24] / ubarra2, abs(d__1)) - 1.;
		v[197] = (d__1 = d3[34] / ubarra3, abs(d__1)) - 1.;
		v[198] = (d__1 = d3[35] / ubarra3, abs(d__1)) - 1.;*/
		
    }

}
