/* 
 * File:   F109Truss942Bar.h
 * Author: hedersb
 *
 * Created on January 22, 2011, 6:34 PM
 */

#ifndef _F109TRUSS942BAR_H
#define	_F109TRUSS942BAR_H

#include "TrussBarStructureStaticProblem.h"

namespace problem {

    class F109Truss942Bar: public TrussBarStructureStaticProblem {
    public:
        F109Truss942Bar();
        F109Truss942Bar(const F109Truss942Bar& orig);
        virtual ~F109Truss942Bar();
    protected:
        virtual int* const getGrouping();
		virtual void fillAreasAux(double* x);
	private:
		int* grouping;
    };

}

#endif	/* _F109TRUSS942BAR_H */

