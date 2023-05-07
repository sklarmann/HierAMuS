// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause



#ifndef FEMEXCEPTIONS_H_
#define FEMEXCEPTIONS_H_

#include <exception>


namespace HierAMuS {

	class FEMExceptions : public std::exception {
	public:
		virtual const char* what() const throw(){
			return "Basic FEM Project Exception happened!";
		}
	};
	
	
	class dofHandlerGenericException : public FEMExceptions {
	public:
		virtual const char* what() const throw(){
			return "Basic Dof Handler Exception happend!";
		}
	};
	
	class dofHandlernewNodeSetException : public dofHandlerGenericException {
	public:
		virtual const char* what() const throw(){
			return "An Error happend when trying to add a new node set!";
		}
	};
	
	
	class nodeSetGenericException : public FEMExceptions {
	public:
		virtual const char* what() const throw(){
			return "Basic Node Set Exception happend!";
		}
	};
	
	class nodeSetSetTypeException : public nodeSetGenericException {
	public:
		virtual const char* what() const throw(){
			return "Error: Tried to change already set node type!";
		}
	};
	
	class nodeSetNumberOfNodesException : public nodeSetGenericException {
	public:
		virtual const char* what() const throw(){
			return "Error: Tried to change the number of nodes in the set!";
		}
	};

} /* namespace HierAMuS */

#endif /* FEMEXCEPTIONS_H_ */
