/* 
 * File:   F103Truss25Bar.h
 * Author: Heder S. Bernardino
 *
 * Created on January 21, 2011, 11:15 AM
 */

#ifndef _F103TRUSS25BAR_H
#define	_F103TRUSS25BAR_H

#include "TrussBarStructureStaticProblem.h"

namespace problem {

    class F103Truss25Bar: public TrussBarStructureStaticProblem {
    public:
        F103Truss25Bar();
        F103Truss25Bar(const F103Truss25Bar& orig);
        virtual ~F103Truss25Bar();
		virtual void evaluation(void* vector, void* values);
	protected:
		virtual int* const getGrouping();
    private:
        int* grouping;
    };

}

#endif	/* _F103TRUSS25BAR_H */

