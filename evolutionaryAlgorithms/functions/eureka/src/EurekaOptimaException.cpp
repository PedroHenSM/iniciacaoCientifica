/* 
 * File:   EurekaOptimaException.cpp
 * Author: hedersb
 * 
 * Created on November 7, 2010, 1:53 PM
 */

#include <string.h>

#include "EurekaOptimaException.h"

namespace exception {


    EurekaOptimaException::EurekaOptimaException(const char* message) {
        this->message = new char[strlen(message)+1];
        strcpy(this->message, message);
    }

    EurekaOptimaException::EurekaOptimaException(const EurekaOptimaException& orig): std::exception(orig) {
        this->message = new char[strlen(orig.message)+1];
        strcpy(this->message, orig.message);
    }

    EurekaOptimaException::~EurekaOptimaException() throw () {
        delete[] this->message;
    }

    EurekaOptimaException& EurekaOptimaException::operator= (const EurekaOptimaException& orig) throw() {
        EurekaOptimaException* newException = new EurekaOptimaException(orig);
        return *newException;
    }

    const char* EurekaOptimaException::what() const throw() {
        return this->message;
    }

}