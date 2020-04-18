/* 
 * File:   F101Truss10Bar.cpp
 * Author: Heder S. Bernardino
 * 
 * Created on January 7, 2011, 11:08 AM
 */

#include "F101Truss10Bar.h"
#include <stdio.h>
#include <memory.h>
#include <iostream>

using namespace std;

namespace problem {

    F101Truss10Bar::F101Truss10Bar(): TrussBarStructureStaticProblem(10, NULL, 280000, 1, 18, 10, 0.1, 25000, 2, "input10.dat", 0.1, 40) {
        
    }

    F101Truss10Bar::F101Truss10Bar(const F101Truss10Bar& orig): TrussBarStructureStaticProblem(orig) {
        
    }

    F101Truss10Bar::~F101Truss10Bar() {
        
    }

}
