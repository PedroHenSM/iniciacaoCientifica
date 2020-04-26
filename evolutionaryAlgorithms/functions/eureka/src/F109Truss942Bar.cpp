/* 
 * File:   F109Truss942Bar.cpp
 * Author: hedersb
 * 
 * Created on January 22, 2011, 6:34 PM
 */

#include "F109Truss942Bar.h"
#include <memory.h>

using namespace std;

namespace problem {

    // F109Truss942Bar::F109Truss942Bar(): TrussBarStructureStaticProblem(59, NULL, 150000, 1, 954/*696+942?*/, 942, 0.1, 25000, 15, "input942.dat", 1, 200), grouping(NULL) { // original line
    F109Truss942Bar::F109Truss942Bar(): TrussBarStructureStaticProblem(59, NULL, 150000, 1, 954/*696+942?*/, 942, 0.1, 25000, 15, "../functions/eureka/input_data/input942.dat", 1, 200), grouping(NULL) {
        
    }

    F109Truss942Bar::F109Truss942Bar(const F109Truss942Bar& orig): TrussBarStructureStaticProblem(orig), grouping(NULL) {
        
    }

    F109Truss942Bar::~F109Truss942Bar() {
        delete[] this->grouping;
    }
	
	int* const F109Truss942Bar::getGrouping () {
		
		if (this->grouping==NULL) {
			//information of the 942-bar truss problem
			this->grouping = new int[ this->getDimension() ];
			this->grouping[0] = 2;
			this->grouping[1] = 10;
			this->grouping[2] = 18;
			this->grouping[3] = 34;
			this->grouping[4] = 46;
			this->grouping[5] = 58;
			this->grouping[6] = 82;
			this->grouping[7] = 86;
			this->grouping[8] = 90;
			this->grouping[9] = 98;
			this->grouping[10] = 106;
			this->grouping[11] = 122;
			this->grouping[12] = 130;
			this->grouping[13] = 162;
			this->grouping[14] = 170;
			this->grouping[15] = 186;
			this->grouping[16] = 194;
			this->grouping[17] = 226;
			this->grouping[18] = 234;
			this->grouping[19] = 258;
			this->grouping[20] = 270;
			this->grouping[21] = 318;
			this->grouping[22] = 330;
			this->grouping[23] = 338;
			this->grouping[24] = 342;
			this->grouping[25] = 350;
			this->grouping[26] = 358;
			this->grouping[27] = 366;
			this->grouping[28] = 382;
			this->grouping[29] = 390;
			this->grouping[30] = 398;
			this->grouping[31] = 430;
			this->grouping[32] = 446;
			this->grouping[33] = 462;
			this->grouping[34] = 486;
			this->grouping[35] = 498;
			this->grouping[36] = 510;
			this->grouping[37] = 558;
			this->grouping[38] = 582;
			this->grouping[39] = 606;
			this->grouping[40] = 630;
			this->grouping[41] = 642;
			this->grouping[42] = 654;
			this->grouping[43] = 702;
			this->grouping[44] = 726;
			this->grouping[45] = 750;
			this->grouping[46] = 774;
			this->grouping[47] = 786;
			this->grouping[48] = 798;
			this->grouping[49] = 846;
			this->grouping[50] = 870;
			this->grouping[51] = 894;
			this->grouping[52] = 902;
			this->grouping[53] = 906;
			this->grouping[54] = 910;
			this->grouping[55] = 918;
			this->grouping[56] = 926;
			this->grouping[57] = 934;
			this->grouping[58] = 942;
		}
		
		return this->grouping;
		
	}
	
	void F109Truss942Bar::fillAreasAux(double* x) {

		double* areas = this->getAreasAux();
		int integer;
		
		int j = 0;
		for(int i = 0; i < this->getDimension(); i++) {
			integer = (int) x[i];
			integer = integer<200? integer: 200;
			for(; j< this->getGrouping()[i]; j++) {
				areas[j] = integer;
			}
		}
		
	}

        //solution of Hasancebi2002
//        x[0] = 1;
//        x[1] = 1;
//        x[2] = 3;
//        x[3] = 1;
//        x[4] = 1;
//        x[5] = 17;
//        x[6] = 3;
//        x[7] = 7;
//        x[8] = 20;
//        x[9] = 1;
//        x[10] = 8;
//        x[11] = 7;
//        x[12] = 19;
//        x[13] = 2;
//        x[14] = 5;
//        x[15] = 1;
//        x[16] = 22;
//        x[17] = 3;
//        x[18] = 9;
//        x[19] = 1;
//        x[20] = 34;
//        x[21] = 3;
//        x[22] = 19;
//        x[23] = 27;
//        x[24] = 42;
//        x[25] = 1;
//        x[26] = 12;
//        x[27] = 16;
//        x[28] = 19;
//        x[29] = 14;
//        x[30] = 42;
//        x[31] = 4;
//        x[32] = 4;
//        x[33] = 4;
//        x[34] = 1;
//        x[35] = 1;
//        x[36] = 62;
//        x[37] = 3;
//        x[38] = 2;
//        x[39] = 4;
//        x[40] = 1;
//        x[41] = 2;
//        x[42] = 77;
//        x[43] = 3;
//        x[44] = 2;
//        x[45] = 3;
//        x[46] = 2;
//        x[47] = 3;
//        x[48] = 100;
//        x[49] = 4;
//        x[50] = 1;
//        x[51] = 4;
//        x[52] = 6;
//        x[53] = 3;
//        x[54] = 49;
//        x[55] = 1;
//        x[56] = 62;
//        x[57] = 1;
//        x[58] = 3;

}
