#include "foo.h"  /* Include the header (not strictly necessary here) */

namespace testNamespace {

    int testClass::bar(int x){
        return x + 10;
    }

    testClass::~testClass(){
        
    }
    // testClass(){};

}

int foo(int x)    /* Function definition */
{
    return x + 5;
}