/* 
 * File:   F101Truss10Bar.h
 * Author: Heder S. Bernardino
 *
 * Created on January 7, 2011, 11:08 AM
 */

#ifndef _F101TRUSS10BAR_H
#define	_F101TRUSS10BAR_H

#include "TrussBarStructureStaticProblem.h"

namespace problem {

    class F101Truss10Bar: public TrussBarStructureStaticProblem {
    public:
        F101Truss10Bar();
        F101Truss10Bar(const F101Truss10Bar& orig);
        virtual ~F101Truss10Bar();
    };

}

#endif	/* _F101TRUSS10BAR_H */

