/* 
 * File:   Problem.cpp
 * Author: hedersb
 * 
 * Created on November 19, 2010, 10:43 AM
 */

#include <iostream>
#include <sstream>
#include "Problem.h"

using namespace std;

namespace problem {

    Problem::Problem(int dimension, void* bounds, int maxNumberObjectiveFunctionEvaluations, int numberObjectives, int numberConstraints) {
        this->dimension = dimension;
        this->numberObjectiveFunctionEvaluations = 0;
        this->maxNumberObjectiveFunctionEvaluations = maxNumberObjectiveFunctionEvaluations;
        this->bounds = bounds;
        this->numberObjectives = numberObjectives;
        this->numberConstraints = numberConstraints;
    }

    Problem::Problem(const Problem& orig) {
        this->dimension = orig.dimension;
        this->numberObjectiveFunctionEvaluations = orig.numberObjectiveFunctionEvaluations;
        this->maxNumberObjectiveFunctionEvaluations = orig.maxNumberObjectiveFunctionEvaluations;
        this->bounds = orig.bounds;
        this->numberObjectives = orig.numberObjectives;
        this->numberConstraints = orig.numberConstraints;
    }

    Problem::~Problem() {
        
    }

    void Problem::evaluate(void* vector, void* values) {
        numberObjectiveFunctionEvaluations++;
        this->evaluation(vector, values);
    }
	
	//TODO - Should I include this here or in the ParallelEurekaOptima project ?
	/*void Problem::evaluate(void** vectors, unsigned int vectorsLength, void* values, unsigned int valuesLength) {
        numberObjectiveFunctionEvaluations += vectorsLength;
        this->evaluation(vectors, vectorsLength, values, valuesLength);
    }*/

    int Problem::getDimension() const {
        return this->dimension;
    }

    int Problem::getNumberObjectiveFunctionEvaluations() const {
        return this->numberObjectiveFunctionEvaluations;
    }

    int Problem::getMaxNumberObjectiveFunctionEvaluations() const {
        return this->maxNumberObjectiveFunctionEvaluations;
    }

    void* Problem::getBounds() const {
        return this->bounds;
    }

    int Problem::getNumberObjectives() const {
        return this->numberObjectives;
    }

    int Problem::getNumberConstraints() const {
        return this->numberConstraints;
    }
	
	string Problem::toString() const {
		ostringstream out;
        out << "----- Problem -----\ndimension: " << dimension << "\nnumber of objectives: " << this->numberObjectives << "\nnumber of constraints: " << this->numberConstraints << "\nmaximum number of objective function evaluations: " << maxNumberObjectiveFunctionEvaluations << endl;
        return out.str();
	}


}
