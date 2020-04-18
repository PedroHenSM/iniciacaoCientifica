/*
 * File:   Problem.h
 * Author: hedersb
 *
 * Created on November 19, 2010, 10:43 AM
 */

#ifndef _PROBLEM_H
#define	_PROBLEM_H

#include <string>
#include <tr1/memory>

using namespace std;
using namespace std::tr1;


namespace problem {

    class Problem {
    public:
        Problem(int dimension, void* bounds, int maxNumberObjectiveFunctionEvaluations, int numberObjectives, int numberConstraints);
        Problem(const Problem& orig);
        virtual ~Problem();
        void evaluate(void* vector, void* values);
		// void evaluate(void** vectors, unsigned int vectorsLength, void* values, unsigned int valuesLength);
        int getNumberObjectiveFunctionEvaluations() const;
        int getMaxNumberObjectiveFunctionEvaluations() const;
        void* getBounds() const;
        int getNumberObjectives() const;
        int getNumberConstraints() const;
        int getDimension() const;
		virtual string toString() const;
    protected:
        int dimension;
        virtual void evaluation(void* vector, void* values) = 0;
		//TODO - Should I include this here or in the ParallelEurekaOptima project ?
		//virtual void evaluation(void** vectors, unsigned int vectorsLength, void** values, unsigned int valuesLength) = 0;
        int numberObjectiveFunctionEvaluations;
        int maxNumberObjectiveFunctionEvaluations;
        void* bounds;
        int numberObjectives;
        int numberConstraints;
    };

    // typedef shared_ptr<Problem> ProblemPtr; // TODO: MODIFIED. 'shared_ptr' is ambiguous
    typedef std::shared_ptr<Problem> ProblemPtr;
}

#endif	/* _PROBLEM_H */

