/* 
 * File:   F107Truss72Bar.h
 * Author: Heder S. Bernardino
 *
 * Created on January 22, 2011, 6:33 PM
 */

#ifndef _F107TRUSS72BAR_H
#define	_F107TRUSS72BAR_H

#include "TrussBarStructureStaticProblem.h"

namespace problem {

    class F107Truss72Bar: public TrussBarStructureStaticProblem {
    public:
        F107Truss72Bar();
        F107Truss72Bar(const F107Truss72Bar& orig);
        virtual ~F107Truss72Bar();
		virtual void evaluation(void* vector, void* values);
	protected:
        virtual int* const getGrouping();
    private:
        int* grouping;
    };

}

#endif	/* _F107TRUSS72BAR_H */

