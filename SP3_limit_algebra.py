### Michael Jarrett
### 3D SP3 Limit
### To calculate SP3 limit of 2D/1D method with anisotropic TL
#########################################
# This code performs algebra to determine the theoretical limit of 
# accuracy for a 2D/1D method with linear and quadratic axial and 
# radial TL moments.
# There are 4 linear and 6 quadratic moments. Equations for each 
# in terms of the others are determined by hand, and hardcoded into 
# this script. The script then solves the 10x10 system algebraically 
# to determine each moment in terms of derivatives of the isotropic moment.
#########################################

import numpy as np
import math
import sys
import time

###################
### Code design ###
###################

# The system of equations is defined hierarchically.
# Each TL component has an expression. e.g.,
# dEzz/dz = [ 1/3 Lz + 1/5 Lz^2 + 1/7 Lz^3 + ...] [ F - dJ/dz ]
# Each expression is a sum of several terms.
# Each term consists of several subterms, and a variable. 
# Each subterm consists of a coefficient, and an order for the radial and axial operator (0-3), 
# and a component. (e.g., 1/3 Lz)
# If the sum of the order of the operators is greater than (3), the subterm is dropped.

# Two terms can be multiplied. Coefficients and operator orders are multiplied accordingly.

# An expression can be simplified by combining like terms (i.e., terms with the same variable).
# An expression can be solved if there are terms with the same variable as the self-variable.
# A component will have a method to substitute another component in for the variable in one of the terms. 

# The component definitions will be hardcorded. The solution order will also be hardcorded. 

class ComponentList:
    def __init__(self, component_list):
        self.component_list = component_list

    def __del__(self):
        self.component_list = []

class Component:
    def __init__(self, term_list, myTL_ID):
        self.term_list = term_list
        self.CID = myTL_ID
        self.LHS = Term([Subterm(1.0,0,0)], myTL_ID)

    def __del__(self):
        self.term_list = []
        self.CID = 0

    @property
    def CID(self):
        return self.CID
    @CID.setter
    def nterms(self,ID):
        self.CID = ID

    def nterms(self):
         return len(self.term_list)

    def substitute_component(self,thatComponent):
        # substitute the expression for that Component 
        # where it appears in this Component
        print "that Component ID = %i" % thatComponent.CID
        itx = 0 
        newTermList = [] 
        while(itx < self.nterms()):
            tmpTerm = self.term_list[itx]
            print "tmpTerm.TID = %i" % tmpTerm.TID
            if(tmpTerm.TID == thatComponent.CID):
                print "Match at ID = %i" % tmpTerm.TID
                # calculate new terms by distributing
                for thatTerm in thatComponent.term_list:
                    newTerm = term_mult(thatTerm,tmpTerm)
                   
                    newTermList.append(newTerm)

                #delete the old term
                del self.term_list[itx]

            else:
                itx += 1

        #add new term list to Component expression
        self.term_list.extend(newTermList)

    def combine_terms(self):
        itx = 0
        while (itx < self.nterms()):
            tmpTerm1 = self.term_list[itx]
            x = itx+1
            while (x < self.nterms()):
                tmpTerm2 = self.term_list[x]
                if(tmpTerm1.match_ID(tmpTerm2)):
                    tmpTerm1.term_add(tmpTerm2)
                    del self.term_list[x]
                else:
                    x += 1 
            itx += 1

    def remove_high_order(self):
        for myterm in self.term_list:
            myterm.remove_high_order()

    # solve component if the term appears on left and right sides of the equation
    def solve(self):
        # first, ensure the component is "clean" for solving
        self.remove_high_order()
        self.combine_terms()

        ### move self terms to the LHS
        itx = 0
        while(itx < self.nterms()):
            myTerm = self.term_list[itx]
            newSubTermList = []
            if(myTerm.TID == self.CID):
                for mySubTerm in myTerm.subterm_list:
                    newSubTermList.append(Subterm(-1.0*mySubTerm.coeff,mySubTerm.r_order,mySubTerm.z_order))
                self.LHS.subterm_list.extend(newSubTermList)
                del self.term_list[itx]
            itx += 1

        ### invert the operator on the LHS and apply to the RHS
        ### (1 - x)^-1 = 1 + x + x^2 + x^3 + ...

        ### calc extra term x, then calculate x^2, x^3, etc.
         
        extraTerm = Term([],0)

        if(self.LHS.subterm_list[0].isIdentity):
            for mySubTerm in self.LHS.subterm_list[1:]:
                extraTerm.subterm_list.append(mySubTerm)
        else:
            ### TODO: if the first subterm is not I, then adjust the whole term
            print "We've run into a non-identity! Have to code for this case."
        
        ### negate each subterm
        for mySubTerm in extraTerm.subterm_list:
            mySubTerm.coeff = -1.0*mySubTerm.coeff

        ### calculate x^2, x^3
        extraTerm_squared = term_mult(extraTerm,extraTerm)
        #extraTerm_squared.remove_high_order()
        #extraTerm_squared.combine_subterms()

        extraTerm_cubed = term_mult(extraTerm,extraTerm_squared)
        #extraTerm_cubed.remove_high_order()
        #extraTerm_cubed.combine_subterms()
         
        extraTerm.term_add(extraTerm_squared)
        extraTerm.term_add(extraTerm_cubed)

        extraTerm.subterm_list.insert(0,Subterm(1.0,0,0))

        ### now, extra term = I + x + x^2 + x^3
        print "Extra term is:"
        extraTerm.print_info()

        ### next, multiply each term in the component by the new extraTerm
        newTermList = []
        for myTerm in self.term_list:
            tmpTerm = term_mult(myTerm,extraTerm)
            #tmpTerm.remove_high_order() 
            #tmpTerm.combine_subterms() 
            newTermList.append(tmpTerm)

        self.term_list = newTermList

        self.reset_LHS()

    def reset_LHS(self):
        self.LHS.subterm_list = []
        self.LHS = Term([Subterm(1.0,0,0)],0)


    def print_info(self):
        print "My component terms: ID = %i (%s)" % (self.CID,variable_names[self.CID])
        for myterm in self.term_list:
            myterm.print_info()


class Term:
    def __init__(self, subterm_list, variable_ID):
        self.subterm_list = subterm_list
        self.TID = variable_ID
        self.max_order = 3 

    def __del__(self):
        self.subterm_list = []
        self.TID = 0

    
    @property
    def TID(self):
        return self.TID
    @TID.setter
    def TID(self,ID):
        self.TID = ID

    @property
    def max_order(self):
        return self.max_order
    @max_order.setter
    def max_order(self,number):
        self.max_order = number

    def nsubterms(self):
        return len(self.subterm_list)

    def term_add(self, otherTerm):
        if(self.TID == otherTerm.TID):
            self.subterm_list.extend(otherTerm.subterm_list)
        else:
            print "These terms don't add! myTermID = %i, otherTermID = %i" % (self.TID,otherTerm.TID)
        del otherTerm

    # combine subterms with the same r_order and z_order
    def combine_subterms(self):
        ist = 0
        while (ist < self.nsubterms()):
            #print "ist = %i, self.nsubterms = %i" % (ist,self.nsubterms())
            tmpSubTerm1 = self.subterm_list[ist]
            x = ist+1
            while (x < self.nsubterms()):
                #print "x = %i, self.nsubterms = %i" % (x,self.nsubterms())
                tmpSubTerm2 = self.subterm_list[x]
                if(tmpSubTerm1.match_order(tmpSubTerm2)):
                    tmpSubTerm1.subterm_add(tmpSubTerm2)
                    del self.subterm_list[x]
                else:
                    x += 1 
            ist += 1

    # remove terms with total order higher than max_order (3)
    def remove_high_order(self):
        ist = 0
        while (ist < self.nsubterms()):
            tmpSubTerm1 = self.subterm_list[ist]
            if(tmpSubTerm1.total_order() > self.max_order):
                del self.subterm_list[ist]
            else:
                ist += 1 

    # check that this term's component TL ID matches another term's
    def match_ID(self,thatTerm):
        if(self.TID == thatTerm.TID):
            return True
        else:
            return False

    def print_info(self):
        print "Term variable ID: %i (%s)" % (self.TID,variable_names[self.TID])
        for mySubTerm in self.subterm_list:
            mySubTerm.print_info()

class Subterm:
    def __init__(self, coeff, r_order, z_order):
        self.coeff = coeff
        self.r_order = r_order
        self.z_order = z_order

    def __del__(self):
        self.coeff = 0.0
        self.r_order = 0
        self.z_order = 0

    @property
    def coeff(self):
        return self.coeff
    @coeff.setter
    def coeff(self,number):
        self.coeff = number

    @property
    def r_order(self):
        return self.r_order
    @r_order.setter
    def r_order(self,number):
        self.r_order = number

    @property
    def z_order(self):
        return self.z_order
    @z_order.setter
    def z_order(self,number):
        self.z_order = number

    @property
    def coeff(self):
        return self.coeff
    @coeff.setter
    def coeff(self,number):
        self.coeff = number

    def add_r_order(self,number):
        self.r_order += number

    def add_z_order(self,z_order):
        self.z_order += number

    def total_order(self):
        return self.r_order + self.z_order

    # add two subterms
    def subterm_add(self,thatSubTerm):
        if(self.match_order(thatSubTerm)):
            self.coeff = self.coeff + thatSubTerm.coeff 
        else:
            print "These terms don't add! my (r,z) = (%i,%i), their (r,z) = (%i,%i) " % (
                  (self.r_order,self.z_order,thatSubTerm.r_order,thatSubTerm.z_order)   )

    # determine if r_order and z_order matches that of another subterm
    def match_order(self,thatSubTerm):
        if(self.r_order == thatSubTerm.r_order and 
           self.z_order == thatSubTerm.z_order):
            return True
        else:
            return False

    def isIdentity(self):
        if(self.coeff == 1.0 and self.r_order == 0 and self.z_order == 0):
            return True
        else:
            return False
        
    def print_info(self):
        print "Subterm: coeff = %.5f, Lr = %i, Lz = %i" % (self.coeff,self.r_order,self.z_order)

# multiply two terms
# retain TL component ID from first term
def term_mult(myTerm, thatTerm):
    ### this function multiplies two terms
    ### each term has an associated TL component ID.
    ### when multiplying, the ID of the self is retained. Only the 
    ### subterms of otherTerm are retained.

    newSubTermList = []
    for mySubTerm in myTerm.subterm_list:
        for thatSubTerm in thatTerm.subterm_list:
            newSubTerm = subterm_mult(mySubTerm,thatSubTerm)
            newSubTermList.append(newSubTerm)
    
    newTerm = Term(newSubTermList,myTerm.TID)
    newTerm.remove_high_order()
    newTerm.combine_subterms()

    return newTerm


# multiply two subterms
def subterm_mult(mySubTerm,thatSubTerm):
    coeff = mySubTerm.coeff * thatSubTerm.coeff
    r_order = mySubTerm.r_order + thatSubTerm.r_order
    z_order = mySubTerm.z_order + thatSubTerm.z_order

    return Subterm(coeff,r_order,z_order)


### indexing for each TL component 

LIN_AX_POL   =  1
LIN_AX_AZI   =  2
LIN_RAD_POL  =  3
LIN_RAD_AZI  =  4
QUAD_AX_POL  =  5
QUAD_AX_AZI  =  6
QUAD_AX_XXX  =  7
QUAD_RAD_POL =  8
QUAD_RAD_AZI =  9
QUAD_RAD_XXX = 10
SOURCE_2D    = 11
SOURCE_1D    = 12


# names for each TL component variable
variable_names = ['null',
                  'Linear polar axial TL',
                  'Linear azimuthal axial TL',
                  'Linear polar radial TL',
                  'Linear azimuthal radial TL',
                  'Quadratic polar axial TL',
                  'Quadratic azimuthal axial TL',
                  'Quadratic cross moment axial TL',
                  'Quadratic polar radial TL',
                  'Quadratic azimuthal radial TL',
                  'Quadratic cross moment radial TL',
                  'Isotropic 2D MOC source',
                  'Isotropic 1D axial source']

if __name__ == '__main__':

    #########################
    #### TESTING METHODS ####
    #########################

    test_methods = True
    if(test_methods):

        print variable_names

        tmpTermList = []

        tmpSubTermList = []
        tmpSubTermList.append(Subterm(1.0/3.0,1,0))
        tmpSubTermList.append(Subterm(1.0/5.0,2,0))
        tmpSubTermList.append(Subterm(1.0/7.0,3,0))

        tmpTermList.append(Term(tmpSubTermList,SOURCE_2D))

        tmpSubTermList = []
        tmpSubTermList.append(Subterm(5.0/21.0,1,1))
        tmpSubTermList.append(Subterm(5.0/63.0,1,2))
        tmpSubTermList.append(Subterm(5.0/99.0,2,1))

        tmpTermList.append(Term(tmpSubTermList,SOURCE_1D))
        
        myComponent = Component(tmpTermList,LIN_AX_POL)

        tmpSubTermList = []
        tmpSubTermList.append(Subterm(1.0/7.0,0,1))
        tmpSubTermList.append(Subterm(2.0/35.0,0,2))
        tmpSubTermList.append(Subterm(2.0/105.0,1,1))

        myTerm = Term(tmpSubTermList,QUAD_AX_XXX)

        myComponent.term_list.append(myTerm)

        print "My component terms: "
        for myterm in myComponent.term_list:
            thisID = myterm.TID
            print "Term variable ID: %i (%s)" % (thisID,variable_names[thisID])
            for mySubTerm in myterm.subterm_list:
                print "Subterm: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm.coeff,mySubTerm.r_order,mySubTerm.z_order)


        ### test mult and add subterms
        myTerm1 = myComponent.term_list[1]
        mySubTerm1 = myTerm1.subterm_list[0]

        myTerm2 = myComponent.term_list[2]
        mySubTerm2 = myTerm2.subterm_list[2]

        ### test add_term and combine_subterms
        myTerm1.term_add(myTerm2) ### should not go, terms have different TIDs

        myTerm1.TID = myTerm2.TID
        myTerm1.term_add(myTerm2) 
        print "myTerm1 + myTerm2 = "
        for mySubTerm in myTerm1.subterm_list:
            print "Subterm: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm.coeff,mySubTerm.r_order,mySubTerm.z_order)

        myTerm1.combine_subterms()
        print "myTerm1 + myTerm2 = "
        for mySubTerm in myTerm1.subterm_list:
            print "Subterm: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm.coeff,mySubTerm.r_order,mySubTerm.z_order)

        print "Subterm 1: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm1.coeff,mySubTerm1.r_order,mySubTerm1.z_order)
        print "Subterm 2: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm2.coeff,mySubTerm2.r_order,mySubTerm2.z_order)
        mySubTerm1.subterm_add(mySubTerm2)
        print "Subterm sum: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm1.coeff,mySubTerm1.r_order,mySubTerm1.z_order)

        mySubTerm1 = myTerm1.subterm_list[0]
        print "Term->Subterm sum: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm1.coeff,mySubTerm1.r_order,mySubTerm1.z_order)

        print "Subterm 1: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm1.coeff,mySubTerm1.r_order,mySubTerm1.z_order)
        print "Subterm 2: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm2.coeff,mySubTerm2.r_order,mySubTerm2.z_order)
        newSubTerm = subterm_mult(mySubTerm1,mySubTerm2)
        mySubTerm1 = newSubTerm
        print "Subterm product: coeff = %.5f, Lr = %i, Lz = %i" % (mySubTerm1.coeff,mySubTerm1.r_order,mySubTerm1.z_order)

        ### test component subsitution

        print " \n \n Testing component substitution!!!! \n "

        tmpTermList = []

        tmpSubTermList = []
        tmpSubTermList.append(Subterm(1.0/63.0,1,0))
        tmpSubTermList.append(Subterm(8.0/105.0,0,1))
        tmpSubTermList.append(Subterm(3.0/147.0,1,1))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_AZI))

        tmpSubTermList = []
        tmpSubTermList.append(Subterm(1.0/7.0,1,0))
        tmpSubTermList.append(Subterm(3.0/175.0,1,1))
        tmpSubTermList.append(Subterm(5.0/33.0,2,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_AZI))

        tmpSubTermList = []
        tmpSubTermList.append(Subterm(5.0/21.0,0,1))
        tmpSubTermList.append(Subterm(5.0/63.0,0,2))
        tmpSubTermList.append(Subterm(5.0/99.0,2,1))

        tmpTermList.append(Term(tmpSubTermList,LIN_AX_POL))

        newComponent = Component(tmpTermList,QUAD_AX_XXX)
        
        myComponent.term_list[2].TID = 8
        myComponent.remove_high_order()
        print "Component #1"
        myComponent.print_info()
        print "Component #2"
        newComponent.print_info()

        ### substitute new component into myComponent
        myComponent.substitute_component(newComponent)

        print "Combination:"
        myComponent.print_info()

        print "Solving component..."
        myComponent.solve()
        print "New expression:"
        myComponent.print_info()
       

    logname = "logfile"
    
    logfile = open(logname,"w")

    #logstr = "iter: %i, keff = %.5f, kdiff = %.4e, fissdiff = %.4e \n" % (t+1,keff,kdiff,fissdiff)
    t = 0
    keff = 1.0
    kdiff = 0.0
    fissdiff = 0.0
    logstr = "iter: %i, keff = %.5f, kdiff = %.4e, fissdiff = %.4e \n" % (t+1,keff,kdiff,fissdiff)
    logfile.write(logstr)
    #print "iter: %i, keff = %.5f, kdiff = %.4e, fissnorm = %.4e, fissdiff = %.4e " % (t+1,keff,kdiff,fiss_norm,fissdiff)
    t=t+1

    logfile.close()
