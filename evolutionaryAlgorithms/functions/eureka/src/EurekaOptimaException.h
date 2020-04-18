/* 
 * File:   EurekaOptimaException.h
 * Author: hedersb
 *
 * Created on November 7, 2010, 1:53 PM
 */

#ifndef _EUREKAOPTIMAEXCEPTION_H
#define	_EUREKAOPTIMAEXCEPTION_H

#include <exception>

namespace exception {

    class EurekaOptimaException : public std::exception {
    public:
        EurekaOptimaException(const char* message);
        EurekaOptimaException(const EurekaOptimaException& orig);
        virtual ~EurekaOptimaException() throw ();
        EurekaOptimaException& operator= (const EurekaOptimaException& orig) throw();
        virtual const char* what() const throw();
    private:
        char* message;
    };

}

#endif	/* _EUREKAOPTIMAEXCEPTION_H */

