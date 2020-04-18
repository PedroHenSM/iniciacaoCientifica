/* 
 * File:   F105Truss60Bar.h
 * Author: Heder S. Bernardino
 *
 * Created on January 22, 2011, 6:32 PM
 */

#ifndef _F105TRUSS60BAR_H
#define	_F105TRUSS60BAR_H

#include "TrussBarStructureStaticProblem.h"

namespace problem {

    class F105Truss60Bar: public TrussBarStructureStaticProblem {
    public:
        F105Truss60Bar();
        F105Truss60Bar(const F105Truss60Bar& orig);
        virtual ~F105Truss60Bar();
		virtual void evaluation(void* vector, void* values);
	protected:
        virtual int* const getGrouping();
		virtual void fillAreasAux(double* x);
    private:
        int* grouping;
		double displacementConstraint2;
		double displacementConstraint3;
    };

}

#endif	/* _F105TRUSS60BAR_H */

