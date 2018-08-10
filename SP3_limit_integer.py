### Michael Jarrett
### 3D SP3 Limit
### To calculate SP3 limit of 2D/1D method with anisotropic TL
#########################################
# This code performs algebra to determine the theoretical limit of 
# accuracy for a 2D/1D method with linear axial and radial TL moments.
# There are 4 linear moments. Equations for each 
# in terms of the others are determined by hand, and hardcoded into 
# this script. The script then solves the 10x10 system algebraically 
# to determine each moment in terms of derivatives of the isotropic moment.
#########################################

import numpy as np
import math
import sys
import time
from subprocess import call

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

    def ncomponents(self):
        return len(self.component_list)

    def add_component(self,thatComponent):
        ic = 0
        while(ic < self.ncomponents()):
            if(self.component_list[ic].CID == thatComponent.CID):
                del self.component_list[ic]
            else:
                ic += 1

        self.component_list.append(thatComponent)

    def find_component(self,thatID):
        for myComponent in self.component_list:
            if(myComponent.CID == thatID):
                return myComponent

class Component:
    def __init__(self, term_list, myTL_ID):
        self.term_list = term_list
        self.CID = myTL_ID
        self.LHS = Term([Subterm(1.0,0,0)], myTL_ID)
        #self.LHS = Term([Subterm([],[],0,0)], myTL_ID)

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
        #print "that Component ID = %i" % thatComponent.CID
        itx = 0 
        newTermList = [] 
        #print "thatComponenet.CID = %i" % thatComponent.CID
        while(itx < self.nterms()):
            tmpTerm = self.term_list[itx]
            #print "tmpTerm.TID = %i" % tmpTerm.TID
            if(tmpTerm.TID == thatComponent.CID):
                #print "Match at ID = %i" % tmpTerm.TID
                # calculate new terms by distributing
                for thatTerm in thatComponent.term_list:
                    newTerm = term_mult(thatTerm,tmpTerm)
                   
                    newTermList.append(newTerm)

                #delete the old term
                del self.term_list[itx]

            else:
                itx += 1

        # add new term list to Component expression
        self.term_list.extend(newTermList)

        # combine any like terms that result from substitution
        self.combine_terms()

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
                    newSubTerm = copy_subterm(mySubTerm)
                    newSubTerm.sign = -1*mySubTerm.sign
                    newSubTermList.append(newSubTerm)
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
            print "*********************************************************."
            print "We've run into a non-identity! Have to code for this case."
            print "*********************************************************."
        
        ### negate each subterm
        for mySubTerm in extraTerm.subterm_list:
            mySubTerm.sign = -1*mySubTerm.sign

        ### calculate x^2, x^3
        extraTerm_squared = term_mult(extraTerm,extraTerm)
        extraTerm_cubed = term_mult(extraTerm,extraTerm_squared)
         
        extraTerm.term_add(extraTerm_squared)
        extraTerm.term_add(extraTerm_cubed)

        extraTerm.subterm_list.insert(0,Subterm(1.0,0,0))
        #extraTerm.subterm_list.insert(0,Subterm([],[],0,0))

        ### now, extra term = I + x + x^2 + x^3
        print "Extra term is:"
        extraTerm.print_info()

        ### next, multiply each term in the component by the new extraTerm
        newTermList = []
        for myTerm in self.term_list:
            tmpTerm = term_mult(myTerm,extraTerm)
            newTermList.append(tmpTerm)

        self.term_list = newTermList

        self.reset_LHS()

    # perform final solution for expression F = [ ... ] \phi
    # requires that the component be an expression in the form \phi = [ ... ] F
    def solve_F(self):
        self.remove_high_order()
        self.combine_terms()

        if(self.nterms() == 1 and 
           self.term_list[0].TID == SOURCE_3D and
           self.CID == SCALAR_FLUX):

            extraTerm = Term([],SCALAR_FLUX)
            FTerm = self.term_list[0]
            if(FTerm.subterm_list[0].isIdentity):
                for mySubTerm in FTerm.subterm_list[1:]:
                    newSubTerm = copy_subterm(mySubTerm)
                    #newSubTerm.sign = -1*mySubTerm.sign
                    extraTerm.subterm_list.append(newSubTerm)

            else:
                ### TODO: if the first subterm is not I, then adjust the whole term
                print "*********************************************************."
                print "We've run into a non-identity! Have to code for this case."
                print "*********************************************************."
        else:
            print "called solve_F but the component is not ready!"
            self.print_component()

        ### negate each subterm
        for mySubTerm in extraTerm.subterm_list:
            #mySubTerm.coeff = -1.0*mySubTerm.coeff
            #mySubTerm.numerator.insert(0,[-1])
            mySubTerm.sign = -1*mySubTerm.sign

        ### calculate x^2, x^3
        print "extraTerm subterms:"
        extraTerm.print_term()
        print "calculating extraTerm_squared. subterms:"
        extraTerm_squared = term_mult(extraTerm,extraTerm)
        extraTerm_squared.print_term()
        print "calculating extraTerm_cubed. subterms:"
        extraTerm_cubed = term_mult(extraTerm,extraTerm_squared)
        extraTerm_cubed.print_term()

        extraTerm.term_add(extraTerm_squared)
        extraTerm.term_add(extraTerm_cubed)

        extraTerm.subterm_list.insert(0,Subterm(1.0,0,0))
        #extraTerm.subterm_list.insert(0,Subterm([],[],0,0))

        ### now, extra term = I + x + x^2 + x^3
        ### combine subterms
        extraTerm.combine_subterms()
        ### sort subterms
        extraTerm.sort_subterms()
        print "Extra term is:"
        extraTerm.print_info()

        newComponent = Component([extraTerm],SOURCE_3D)
        print "The solution is:"
        newComponent.print_component()
        newComponent.write_latex_equation()

    ### multiplies entire component by halfLZ
    def mult_halfLZ(self):
        halfLZterm = Term([Subterm(1.0,0,0.5)],0)
        #halfLZterm = Term([Subterm([],[],0,0.5)],0)
        newTermList = []
        for myTerm in self.term_list:
            tmpTerm = term_mult(myTerm,halfLZterm)
            newTermList.append(tmpTerm)

        self.term_list = newTermList

    def reset_LHS(self):
        self.LHS.subterm_list = []
        self.LHS = Term([Subterm(1.0,0,0)],0)
        #self.LHS = Term([Subterm([],[],0,0)],0)


    def print_info(self):
        print "My component terms: ID = %i (%s)" % (self.CID,variable_names[self.CID])
        for myterm in self.term_list:
            myterm.print_info()

    def print_component(self):
        print "My component terms: ID = %i (%s)" % (self.CID,variable_names[self.CID])
        for myterm in self.term_list:
            myterm.print_term()

    def write_latex_equation(self):
        texfile = open(texname,"a")

        texstr = "\\begin{alignat}{2} \n"
        texfile.write(texstr)
        texstr = "\\label{eq:%s_%i} \n" % (equation_labels[self.CID],equation_counter.get_index(self.CID))
        texfile.write(texstr)
        #texstr = "%s = & \\nonumber \\\\ \n" % (latex_expressions[self.CID])
        texstr = "& %s = \\\\ \n" % (latex_expressions[self.CID])
        texfile.write(texstr)

        for myTerm in self.term_list:
            print "nsubterms = %i" % myTerm.nsubterms()
            print "Printing term info:"
            myTerm.print_info()
            print "done printing term info:"
            texstr = myTerm.write_term_latex()
            texfile.write(texstr)

        texstr = "\\end{alignat} \n"
        texfile.write(texstr)
        texstr = "\n"
        texfile.write(texstr)
        texfile.write(texstr)

class Term:
    def __init__(self, subterm_list, variable_ID):
        self.subterm_list = subterm_list
        self.TID = variable_ID
        #self.max_order = 3 

    def __del__(self):
        self.subterm_list = []
        self.TID = 0

    
    @property
    def TID(self):
        return self.TID
    @TID.setter
    def TID(self,ID):
        self.TID = ID

    #@property
    #def max_order(self):
    #    return self.max_order
    #@max_order.setter
    #def max_order(self,number):
    #    self.max_order = number

    def nsubterms(self):
        return len(self.subterm_list)

    def term_add(self, otherTerm):
        if(self.TID == otherTerm.TID):
            self.subterm_list.extend(otherTerm.subterm_list)
            self.combine_subterms()
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

        self.remove_zero()

    def remove_zero(self):
        ist = 0
        while (ist < self.nsubterms()):
             tmpSubTerm = self.subterm_list[ist]
             if(abs(tmpSubTerm.coeff_float()) < THRESHOLD):
                 del self.subterm_list[ist]
             else:
                 ist += 1

    # remove terms with total order higher than max_order (3)
    def remove_high_order(self):
        ist = 0
        while (ist < self.nsubterms()):
            tmpSubTerm1 = self.subterm_list[ist]
            if(tmpSubTerm1.total_order() > MAX_ORDER):
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

    def print_term(self):
        print "Term variable ID: %i (%s)" % (self.TID,variable_names[self.TID])
        print "( ",
        for mySubTerm in self.subterm_list:
            #if(mySubTerm.coeff_float() > 0.0):
            if(mySubTerm.sign == 1):
                print "+ ",
            else:
                print "- ",
            #print "%7.5f" % (mySubTerm.coeff_float()),
            print "%i/%i" % (mySubTerm.get_numerator(),mySubTerm.get_denominator()),
            if(mySubTerm.r_order > 0):
                if(mySubTerm.r_order == 1):
                    print "Lr",
                else:
                    print "Lr^%i" % mySubTerm.r_order,
            if(mySubTerm.z_order > 0):
                if(mySubTerm.z_order == 1):
                    print "Lz",
                else:
                    print "Lz^%i" % mySubTerm.z_order,
        print ") "

    def write_term_latex(self):
        tmpstr = []

        ### determine if whole term is negative
        allneg = True
        for mySubTerm in self.subterm_list:
            if(mySubTerm.sign == 1):
                allneg = False

        if(allneg):
            tmpstr.append('& -')
        else:
            tmpstr.append('& +')

        tmpstr.append(' \\l[ ')
        firstTerm = True
        for mySubTerm in self.subterm_list:

            if(mySubTerm.coeff_float() == 1.0):
                mystr = ""
            elif(mySubTerm.coeff_float() == -1.0):
                if(allneg == False):
                    mystr = "- "
                else:
                    mystr = ""
            else:
                if(mySubTerm.sign == 1 or allneg):
                    if(not firstTerm): 
                        tmpstr.append('+ ')
                else:
                    if(not firstTerm):
                        tmpstr.append('- ')

                mystr = '\\frac{%i}{%i} ' % (mySubTerm.get_numerator(),mySubTerm.get_denominator())
            tmpstr.append(mystr)

            if(mySubTerm.r_order > 0):
                if(mySubTerm.r_order == 1):
                    tmpstr.append('\\opL_r ')
                else:
                    mystr = '\\opL_r^{%i} ' % mySubTerm.r_order
                    tmpstr.append(mystr)

            intzorder = int(mySubTerm.z_order)
            if(intzorder > 0):
                if(intzorder == 1):
                    tmpstr.append('\\opL_z ')
                else:
                    mystr = '\\opL_z^{%i} ' % intzorder
                    tmpstr.append(mystr)

            if(mySubTerm.r_order == 0 and intzorder == 0):
                mystr = "I "
                tmpstr.append(mystr)
            
            firstTerm = False
                
        tmpstr.append(' \\r] ')

        if((2*self.subterm_list[0].z_order)%2 == 1): # half LZ is present
            tmpstr.append('\halfLZ ') 

        tmpstr.append(latex_expressions[self.TID])
        tmpstr.append(' \\nonumber \\\\ \n')

        #print tmpstr
        nospace = ''
        #str_return = nospace.join(['\\l[','\\r]'])
        return nospace.join(tmpstr)

    def sort_subterms(self):
        # first, ensure higher order terms are removed
        self.remove_high_order() 
        # and like terms are combined
        self.combine_subterms() 

        ist = 0
        st_order = []
        while(ist < self.nsubterms()):
            thisSubTerm = self.subterm_list[ist]
            st_order.append(thisSubTerm.r_order + thisSubTerm.z_order)
            ist += 1

        # st_order is a list of the order of each subterm
        # sorted_order is a list of the old subterm indices 
        # in the order they should be placed
        sorted_order = [0] * self.nsubterms()
        sorted_operator_order = [0] * self.nsubterms()
        
        x = 0
        # go up by 1/2, orders can be non-zero
        tot_order = 0.0
        while(tot_order <= MAX_ORDER):
            for ist in range(0,self.nsubterms()):
                if(st_order[ist] == tot_order):
                    sorted_order[x] = ist
                    sorted_operator_order[x] = st_order[ist] 
                    x += 1
            tot_order += 0.5

        ### pass 10 times to get it right
        for i in range(0,10):
            for x in range(0,self.nsubterms()-1):
                ### sort within a group of similar order
                if(sorted_operator_order[x] == sorted_operator_order[x+1]):
                    ist = sorted_order[x]
                    istp = sorted_order[x+1]
                    ### put non-cross terms first
                    if((self.subterm_list[ist].r_order > 0.0 and self.subterm_list[ist].z_order > 0.0) and
                       (self.subterm_list[istp].r_order == 0.0 or self.subterm_list[istp].z_order == 0.0)):
                        tmp = sorted_order[x]
                        sorted_order[x] = sorted_order[x+1]
                        sorted_order[x+1] = tmp
                    ### put higher r_order first
                    elif((self.subterm_list[ist].z_order > self.subterm_list[ist].r_order) and
                         (self.subterm_list[istp].r_order > self.subterm_list[istp].z_order)):
                        tmp = sorted_order[x]
                        sorted_order[x] = sorted_order[x+1]
                        sorted_order[x+1] = tmp

        ### generate new subterm list in the correct order
        newSubTermList = []
        for x in range(0,self.nsubterms()):
            tmpSubTerm = copy_subterm(self.subterm_list[sorted_order[x]])
            newSubTermList.append(tmpSubTerm)

        self.subterm_list = newSubTermList
            

class Subterm:
    #def __init__(self, numerator, denominator, r_order, z_order):
    def __init__(self, coeff, r_order, z_order):
        if(coeff < 0.0): 
            self.sign = -1
        else:
            self.sign = 1

        i = compare_fraction(coeff)
        if(i >= 0):
            myPair = common_fractions[i]
            self.numerator = compute_prime_factors(myPair[0])
            self.denominator = compute_prime_factors(myPair[1])
            self.r_order = r_order
            self.z_order = z_order
            #print "Found the fraction! value = %8.6f (%i/%i)" % (self.coeff_float(),myPair[0],myPair[1])
        else:
            print "Did not find the fraction! value = %8.6f" % coeff
            self.numerator = []
            self.denominator = []
            self.r_order = 0
            self.z_order = 0

    def __del__(self):
        #self.coeff = 0.0
        self.numerator = []
        self.denominator = []
        self.r_order = 0
        self.z_order = 0
        self.sign = 1

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
    def sign(self):
        return self.sign
    @z_order.setter
    def sign(self,sign):
        self.sign = sign

    def nnumer(self):
        return len(self.numerator)
    def ndenom(self):
        return len(self.denominator)

    def coeff_float(self):
        numval = self.get_numerator()
        denomval = self.get_denominator()
        return float(self.sign)*float(numval)/float(denomval)

    def get_numerator(self):
        numval = 1
        for val in self.numerator:
            numval *= val
        return numval

    def get_denominator(self):
        denomval = 1
        for val in self.denominator:
            denomval *= val
        return denomval

    def reduce_fraction(self):
        self.sort_num_denom()
        ival = 0
        while(ival < self.nnumer()):
            jval = 0
            while(jval < self.ndenom() and ival < self.nnumer()):
                #print "ival = %i, nnum = %i, jval = %i, ndenom = %i" % (ival,self.nnumer(),jval,self.ndenom())
                if(self.numerator[ival] == self.denominator[jval]):
                ### common multiple, eliminate this term
                    del self.numerator[ival]
                    del self.denominator[jval]
                else:
                    jval += 1
            ival += 1

    def sort_num_denom(self):
        self.numerator.sort()
        self.denominator.sort()

    def add_r_order(self,number):
        self.r_order += number

    def add_z_order(self,z_order):
        self.z_order += number

    def total_order(self):
        return self.r_order + self.z_order

    # add two subterms
    def subterm_add(self,thatSubTerm):
        if(self.match_order(thatSubTerm)):

            tmpStore    = Subterm(1.0,0,0)
            tmpStore.numerator = list(self.denominator)
            tmpStore.denominator = list(thatSubTerm.denominator)

            ### reduce fractions
            common_factors = tmpStore.find_common_factors()
            #tmpStore.reduce_fraction()

            ### compute sum of the two fractions
            numeratorA = self.get_numerator() * tmpStore.get_denominator()
            numeratorB = tmpStore.get_numerator() * thatSubTerm.get_numerator()
            tmpNum = self.sign*numeratorA + thatSubTerm.sign*numeratorB 
            if(tmpNum >= 0.0):
                self.sign = 1
            else:
                self.sign = -1
            #tmpSubTerm1.denominator.extend(tmpSubTerm2.denominator)
            tmpStore.denominator.extend(tmpStore.numerator)
            tmpStore.denominator.extend(common_factors)
            self.denominator = tmpStore.denominator

            ### old way
            #print "Adding two terms: term1 = %i/%i, term2 = %i/%i" % (self.get_numerator(),
            #   self.get_denominator(),thatSubTerm.get_numerator(),thatSubTerm.get_denominator())
            ### calculate numerator
            #numeratorA = self.get_numerator() * thatSubTerm.get_denominator()
            #numeratorB = self.get_denominator() * thatSubTerm.get_numerator()
            #print "numeratorA = %6.1f, numeratorB = %6.1f" % (numeratorA,numeratorB)
            #tmpNum = self.sign*numeratorA + thatSubTerm.sign*numeratorB

            #tmpDenom = self.get_denominator() * thatSubTerm.get_denominator()
            #self.denominator.extend(thatSubTerm.denominator)

            #print "tmpnum = %i, tmpDenom = %i" % (tmpNum, tmpDenom)
            self.numerator = compute_prime_factors(int(abs(tmpNum)))
            #self.denominator = compute_prime_factors(tmpDenom)
            self.reduce_fraction()

        else:
            print "These terms don't add! my (r,z) = (%i,%i), their (r,z) = (%i,%i) " % (
                  (self.r_order,self.z_order,thatSubTerm.r_order,thatSubTerm.z_order)   )

    def find_common_factors(self):
        # find common factors between numerator and denominator
        # this is used for adding two fractions
        common_factors = []
        ifac = 0
        while(ifac < self.nnumer()):
            jfac = 0
            while(jfac < self.ndenom() and ifac < self.nnumer()):
                if(self.numerator[ifac] == self.denominator[jfac]):
                    common_factors.append(int(self.numerator[ifac]))
                    del self.numerator[ifac]
                    del self.denominator[jfac]
                else:
                    jfac += 1
            ifac += 1

        return common_factors

    # determine if r_order and z_order matches that of another subterm
    def match_order(self,thatSubTerm):
        if(self.r_order == thatSubTerm.r_order and 
           self.z_order == thatSubTerm.z_order):
            return True
        else:
            return False

    def isIdentity(self):
        if(self.coeff_float() == 1.0 and self.r_order == 0 and self.z_order == 0):
            return True
        else:
            return False
        
    def print_info(self):
        #print "Subterm: coeff = %.5f, Lr = %i, Lz = %i" % (self.coeff,self.r_order,self.z_order)
        print "Subterm: coeff = %.5f, (%i/%i), Lr = %i, Lz = %i" % (self.coeff_float(),
                  self.get_numerator(),self.get_denominator(),self.r_order,self.z_order)

class Counter:
    def __init__(self,nvalues):
        self.nvalues = nvalues
        self.count_list = [0] * (nvalues+1)

    def __del__(self):
        self.nvalues = 0
        self.count_list = []

    def get_index(self,ID):
        index = self.count_list[ID]
        self.count_list[ID] += 1
        return index


# multiply two terms
# retain TL component ID from first term
def term_mult(myTerm, thatTerm):
    ### this function multiplies two terms
    ### each term has an associated TL component ID.
    ### when multiplying, the ID of the self is retained. Only the 
    ### subterms of otherTerm are retained.

    #print "Multiplying myTerm: "
    #myTerm.print_term()
    #print "by thatTerm: "
    #thatTerm.print_term()
    newSubTermList = []
    for mySubTerm in myTerm.subterm_list:
        for thatSubTerm in thatTerm.subterm_list:
            newSubTerm = subterm_mult(mySubTerm,thatSubTerm)
            newSubTermList.append(newSubTerm)
    
    newTerm = Term(newSubTermList,myTerm.TID)
    newTerm.remove_high_order()
    ### combining subterms
    newTerm.print_info()
    newTerm.sort_subterms()

    return newTerm


# multiply two subterms
def subterm_mult(mySubTerm,thatSubTerm):
    newSubTerm = copy_subterm(mySubTerm)
    thatCopy = copy_subterm(thatSubTerm)

    newSubTerm.numerator.extend(thatCopy.numerator)
    newSubTerm.denominator.extend(thatCopy.denominator)
    newSubTerm.reduce_fraction()

    newSubTerm.sign = mySubTerm.sign*thatSubTerm.sign

    newSubTerm.r_order = mySubTerm.r_order + thatSubTerm.r_order
    newSubTerm.z_order = mySubTerm.z_order + thatSubTerm.z_order

    #print "Subterm mult result:"
    #print newSubTerm.numerator,newSubTerm.denominator

    return newSubTerm

def copy_subterm_list(thatSubTermList):
    newSubTermList = []
    for mySubTerm in thatSubTermList:
        newSubTermList.append(copy_subterm(mySubTerm))
    return newSubTermList

def copy_subterm(thatSubTerm):
    newSubTermCopy = Subterm(1.0,thatSubTerm.r_order,thatSubTerm.z_order)
    newSubTermCopy.numerator = list(thatSubTerm.numerator)
    newSubTermCopy.denominator = list(thatSubTerm.denominator)
    newSubTermCopy.sign = thatSubTerm.sign
    return newSubTermCopy

def compare_fraction(thatFloat):
    for i in range(0,nfrac):
        myPair = common_fractions[i]
        myFloat = float(myPair[0])/float(myPair[1])
        if(abs(myFloat - abs(thatFloat)) < THRESHOLD*abs(thatFloat)): # this is a match
            return i           
    return -1

def compare_fraction_sum(thatFloat):
    for i in range(0,nfrac):
        for j in range(0,nfrac):
            myPair1 = common_fractions[i]
            myPair2 = common_fractions[j]
            myFloat1 = float(myPair1[0])/float(myPair1[1])
            myFloat2 = float(myPair2[0])/float(myPair2[1])
            if(abs(myFloat1 + myFloat2 - abs(thatFloat)) < THRESHOLD*abs(thatFloat)): # this is a match
                return [i,j]           
    return [-1,-1]

def compute_prime_factors(myInt):
    factor_list = []
    i = 2
    while(myInt > 1):
        if(myInt%i == 0):
            ### i is a factor of myInt. is it prime? it should be.
            myInt = myInt/i 
            factor_list.append(i)
        else:
            i += 1
            #if(is_prime(i)):

    return factor_list

#def is_prime(n):
#  if(n == 2 or n == 3): return True
#  if(n < 2 or n%2 == 0): return False
#  if(n < 9): return True
#  if(n%3 == 0): return False
#  r = int(n**0.5)
#  f = 5
#  while(f <= r):
#    if(n%f == 0): return False
#    if(n%(f+2) == 0): return False
#    f +=6
#  return True 

### indexing for each TL component 

### linear TL terms
LIN_AX_POL   =  1
LIN_AX_AZI   =  2
LIN_RAD_POL  =  3
LIN_RAD_AZI  =  4
### quadratic TL terms
QUAD_AX_POL  =  5
QUAD_AX_AZI  =  6
QUAD_AX_XXX  =  7
QUAD_RAD_POL =  8
QUAD_RAD_AZI =  9
QUAD_RAD_XXX = 10
### isotropic TL terms
ISO_AX_TL    = 11 # dJ/dz
ISO_RAD_TL   = 12 # (dJ/dx + dJ/dy)
### isotropic source terms
SOURCE_2D    = 13 # F - dJ/dz
SOURCE_1D    = 14 # F - (dJ/dx + dJ/dy)
### flux 
SCALAR_FLUX  = 15
### F: 3D source (this is the component we ultimately need to solve for)
SOURCE_3D  = 16

nID = 16

MAX_ORDER = 3
THRESHOLD = 1.0E-6

if(int(sys.argv[1]) == 1):
    texname = "equations_output_integer_linear.tex"
if(int(sys.argv[1]) == 2):
    texname = "equations_output_integer.tex"

equation_counter = Counter(nID)

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
                  'Isotropic axial TL',
                  'Isotropic radial TL',
                  'Isotropic 2D MOC source',
                  'Isotropic 1D axial source',
                  'Scalar Flux',
                  'Source (F)']

latex_expressions = ['',
                     '\\axpolTL',
                     '\\axaziTL',
                     '\\radpolTL',
                     '\\radaziTL',
                     '\\quadaxpolTL',
                     '\\quadaxaziTL',
                     '\\quadaxcrossTL',
                     '\\quadradpolTL',
                     '\\quadradaziTL',
                     '\\quadradcrossTL',
                     '\\ATL',
                     '\\RTL',
                     '\\l[ F - \ATL \\r]',
                     '\\l[ F - \RTL \\r]',
                     '\phi',
                     'F']

equation_labels   = ['',
                     'linaxpolTL',
                     'linaxaziTL',
                     'linradpolTL',
                     'linradaziTL',
                     'quadaxpolTL',
                     'quadaxaziTL',
                     'quadaxcrossTL',
                     'quadradpolTL',
                     'quadradaziTL',
                     'quadradcrossTL',
                     'isoATL',
                     'isoRTL',
                     'transport2D',
                     'transport1D',
                     'scalarFlux',
                     'transport3D']

##### commonly encountered fractions for LaTeX output
common_fractions = [(1,  1),
                    (1,  3),
                    (1,  5),
                    (1,  7),
                    (1,  9),
                    (1, 11),
                    (1, 15),
                    (1, 21),
                    (1, 25),
                    (1, 35),
                    (1, 63),
                    (1, 99),
                    (1,105),
                    (1,126),
                    (3,  3),
                    (3,  5),
                    (3,  7),
                    (3,  9),
                    (3, 11),
                    (3, 15),
                    (3, 21),
                    (3, 35),
                    (3, 63),
                    (3, 99),
                    (3,105),
                    (3,126),
                    (5,  3),
                    (5,  5),
                    (5,  7),
                    (5,  9),
                    (5, 11),
                    (5, 15),
                    (5, 21),
                    (5, 33),
                    (5, 35),
                    (5, 63),
                    (5, 99),
                    (5,105),
                    (5,126),
                    (7,  3),
                    (7,  5),
                    (7,  7),
                    (7,  9),
                    (7, 11),
                    (7, 15),
                    (7, 21),
                    (7, 35),
                    (7, 63),
                    (7, 99),
                    (7,105),
                    (7,126)]

nfrac = len(common_fractions)

if __name__ == '__main__':


    ### clear the equations_output.tex file
    call(["rm",texname])
    call(["touch",texname])

    ### linear
    order = int(sys.argv[1])

    if(order == 1): # linear
        # set up component list
        tmpComponentList = [] 

        ############## SCALAR FLUX (\phi) #############
        tmpTermList      = [] 
        ### SOURCE_2D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0,0,0))
        tmpSubTermList.append(Subterm(1.0/3.0,1,0))
        tmpSubTermList.append(Subterm(1.0/5.0,2,0))
        tmpSubTermList.append(Subterm(1.0/7.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_2D))

        ### LIN_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0,0,0))
        tmpSubTermList.append(Subterm(3.0/5.0,1,0))
        tmpSubTermList.append(Subterm(3.0/7.0,2,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_AZI))

        tmpComponentList.append(Component(tmpTermList,SCALAR_FLUX))
        tmpComponentList[-1].write_latex_equation()


        ########### ISO_AX_TL component (dJ/dz) ##########
        tmpTermList      = [] 
        ### SOURCE_1D term
        tmpSubTermList      = [] 
        tmpSubTermList.append(Subterm(-1.0/3.0,0,1))
        tmpSubTermList.append(Subterm(-1.0/5.0,0,2))
        tmpSubTermList.append(Subterm(-1.0/7.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_1D))

        ### LIN_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0,0,0.5))
        tmpSubTermList.append(Subterm(-3.0/5.0,0,1.5))
        tmpSubTermList.append(Subterm(-3.0/7.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_POL))

        tmpComponentList.append(Component(tmpTermList,ISO_AX_TL))
        tmpComponentList[-1].write_latex_equation()

        ######### LIN_RAD_POL component (dJ/dz) #########
        tmpTermList      = [] 
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/5.0,1,0))
        tmpSubTermList.append(Subterm(3.0/35.0,2,0))
        tmpSubTermList.append(Subterm(1.0/21.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_POL))

        tmpComponentList.append(Component(tmpTermList,LIN_RAD_POL))
        tmpComponentList[-1].write_latex_equation()

        ######### LIN_AX_POL component (dJ/dz) #########
        tmpTermList      = [] 

        ### SOURCE_1D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/3.0,0,0.5))
        tmpSubTermList.append(Subterm(1.0/5.0,0,1.5))
        tmpSubTermList.append(Subterm(1.0/7.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_1D))

        ### LIN_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(3.0/5.0,0,1))
        tmpSubTermList.append(Subterm(3.0/7.0,0,2))
        tmpSubTermList.append(Subterm(1.0/3.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_POL))
        tmpComponentList.append(Component(tmpTermList,LIN_AX_POL))
        tmpComponentList[-1].write_latex_equation()

        ### initialize the ComponentList object from the tmpComponentList 
        myComponentList = ComponentList(tmpComponentList)

        ### Eliminate the LIN_AX_POL component by substitution into the LIN_RAD_POL component 
        ### (in the document, this is Eq. (33) into Eq. (34)
        linradpolComp = myComponentList.find_component(LIN_RAD_POL)
        linradpolComp.write_latex_equation()
        linaxpolComp = myComponentList.find_component(LIN_AX_POL)
        linaxpolComp.write_latex_equation()
        linradpolComp.substitute_component(linaxpolComp)
        linradpolComp.write_latex_equation()

        print "Printing the new LIN_RAD_POL component:"
        linradpolComp.print_component()
        print "Solving the new LIN_RAD_POL component:"
        linradpolComp.solve()
        linradpolComp.print_component()
        linradpolComp.write_latex_equation()

        ### we actually need to know d\dz 1/Sigma_t * LIN_RAD_POL = Lz^(1/2) * LIN_RAD_POL
        #linradpolComp.mult_halfLZ()
        #print "After multiplying through by halfLZ:"
        #linradpolComp.print_component()

        ### Eliminate the LIN_RAD_POL component in the ISO_AX_TL term
        ### in the document, this is Eq. (39) into Eq. (31)
        isoaxtlComp = myComponentList.find_component(ISO_AX_TL)
        isoaxtlComp.substitute_component(linradpolComp)
        print "Printing the ISO_AX_TL (dJ/dz) component:"
        isoaxtlComp.print_component()
        isoaxtlComp.write_latex_equation()

        ### SOURCE_1D (F-(dJ/dx + dJ/dy)) is equivalent to SCALAR_FLUX + ISO_AX_TL ( \phi + dJ/dz )
        newTermList = []
        tmpSubTermList = isoaxtlComp.term_list[0].subterm_list
        newTermList.append(Term(tmpSubTermList,SCALAR_FLUX))
        newTermList.append(Term(tmpSubTermList,ISO_AX_TL))
        newisoaxTLComp = Component(newTermList,ISO_AX_TL)
        ### delete old component, overwrite with new one
        myComponentList.add_component(newisoaxTLComp)
        newisoaxTLComp.write_latex_equation()

        ##E solve for ISO_AX_TL
        newisoaxTLComp.solve()
        print "Printing the solved ISO_AX_TL (dJ/dz) component:"
        newisoaxTLComp.print_component()

        ########### LIN_AX_AZI comonent #############
        tmpTermList      = [] 
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/5.0,0,1))
        tmpSubTermList.append(Subterm(3.0/35.0,0,2))
        tmpSubTermList.append(Subterm(1.0/21.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_AZI))

        myComponentList.add_component(Component(tmpTermList,LIN_AX_AZI))
        myComponentList.component_list[-1].write_latex_equation()
        
        ########### LIN_RAD_AZI comonent #############
        tmpTermList      = [] 

        ### SOURCE_2D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/3.0,1,0))
        tmpSubTermList.append(Subterm(1.0/5.0,2,0))
        tmpSubTermList.append(Subterm(1.0/7.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_2D))

        ### LIN_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(3.0/5.0,1,0))
        tmpSubTermList.append(Subterm(3.0/7.0,2,0))
        #tmpSubTermList.append(Subterm(1.0/3.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_AZI))

        myComponentList.add_component(Component(tmpTermList,LIN_RAD_AZI))
        myComponentList.component_list[-1].write_latex_equation()

        ### solve for LIN_AX_AZI
        ### in document, substitute Eq. (45) into Eq. (44) and obtain (47)
        linaxaziComp = myComponentList.find_component(LIN_AX_AZI)
        linradaziComp = myComponentList.find_component(LIN_RAD_AZI)
        linaxaziComp.substitute_component(linradaziComp)
        linaxaziComp.solve()
        print "Printing the solved LIN_AX_AZI component:"
        linaxaziComp.print_component()
        linaxaziComp.write_latex_equation()

        ### solve for phi
        scalarflxComp = myComponentList.find_component(SCALAR_FLUX)
        scalarflxComp.substitute_component(linaxaziComp)
        print "Printing the scalar flux component:"
        scalarflxComp.print_component()
        scalarflxComp.write_latex_equation()

        ### SOURCE_2D = F - dJ/dz
        newTermList = []
        tmpSubTermList1 = copy_subterm_list(scalarflxComp.term_list[0].subterm_list)
        tmpSubTermList2 = copy_subterm_list(scalarflxComp.term_list[0].subterm_list)
        newTermList.append(Term(tmpSubTermList1,SOURCE_3D))
        axTLterm = Term(tmpSubTermList2,ISO_AX_TL)
        ### negate terms for axTLterm
        for mySubTerm in axTLterm.subterm_list:
            #mySubTerm.coeff = -1.0*mySubTerm.coeff
            mySubTerm.sign = -1*mySubTerm.sign

        print "tmpSubTermList 1:"
        print newTermList[0].print_term()
        print "tmpSubTermList 2:"
        print axTLterm.print_term()

        newTermList.append(axTLterm)
        newphiComp = Component(newTermList,SCALAR_FLUX)
        newphiComp.write_latex_equation()

        newisoaxtlComp = myComponentList.find_component(ISO_AX_TL)
        print "Iso AX TL component:"
        newisoaxtlComp.print_component()
        newphiComp.substitute_component(newisoaxtlComp)
        newisoaxtlComp.write_latex_equation()

        ### Expression for PHI and F
        newphiComp.print_component()
        newphiComp.solve()
        newphiComp.write_latex_equation()

        newphiComp.solve_F()


    elif(order == 2): # quadratic
        # set up component list
        tmpComponentList = [] 

        ############## SCALAR FLUX (\phi) #############
        tmpTermList      = [] 
        ### SOURCE_2D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0,0,0))
        tmpSubTermList.append(Subterm(1.0/3.0,1,0))
        tmpSubTermList.append(Subterm(1.0/5.0,2,0))
        tmpSubTermList.append(Subterm(1.0/7.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_2D))

        ### LIN_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0    ,0,0))
        tmpSubTermList.append(Subterm(3.0/5.0,1,0))
        tmpSubTermList.append(Subterm(3.0/7.0,2,0))
        tmpSubTermList.append(Subterm(1.0/3.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_AZI))

        ### QUAD_AX_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-5.0/3.0 ,0,0))
        tmpSubTermList.append(Subterm(-1.0/3.0 ,1,0))
        tmpSubTermList.append(Subterm(-1.0/7.0 ,2,0))
        #tmpSubTermList.append(Subterm(-5.0/63.0,3,0))
        tmpSubTermList.append(Subterm(-5.0/63.0,2,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_POL))

        #print "Test subterm_add:"
        #tmpSubTermList[3].subterm_add(tmpSubTermList[2])

        #print "Test subterm_mult:"
        #temporary = subterm_mult(tmpSubTermList[2],tmpSubTermList[3])

        ### QUAD_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-5.0/3.0,0,0))
        tmpSubTermList.append(Subterm(-1.0    ,1,0))
        tmpSubTermList.append(Subterm(-5.0/7.0,2,0))
        tmpSubTermList.append(Subterm(-5.0/9.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_AZI))

        tmpComponentList.append(Component(tmpTermList,SCALAR_FLUX))
        tmpComponentList[-1].write_latex_equation()


        ############## dJ/dz (iso AXIAL TL) ############
        tmpTermList      = [] 
        ### SOURCE_1D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/3.0,0,1))
        tmpSubTermList.append(Subterm(-1.0/5.0,0,2))
        tmpSubTermList.append(Subterm(-1.0/7.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_1D))

        ### LIN_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0    ,0,0.5))
        tmpSubTermList.append(Subterm(-3.0/5.0,0,1.5))
        tmpSubTermList.append(Subterm(-3.0/7.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_POL))

        ### QUAD_RAD_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/3.0, 0,1))
        tmpSubTermList.append(Subterm(1.0/7.0, 0,2))
        tmpSubTermList.append(Subterm(5.0/63.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_AZI))

        ### QUAD_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0    ,0,1))
        tmpSubTermList.append(Subterm(5.0/7.0,0,2))
        tmpSubTermList.append(Subterm(5.0/9.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_POL))

        tmpComponentList.append(Component(tmpTermList,ISO_AX_TL))
        tmpComponentList[-1].write_latex_equation()


        ############## LIN_RAD_POL component ############
        tmpTermList      = [] 
        ### LIN_AX_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/5.0 ,1,0))
        tmpSubTermList.append(Subterm(3.0/35.0,2,0))
        tmpSubTermList.append(Subterm(1.0/21.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_POL))

        ### QUAD_AX_XXX term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0     ,0,0.5))
        tmpSubTermList.append(Subterm(-3.0/7.0 ,0,1.5))
        tmpSubTermList.append(Subterm(-5.0/21.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_XXX))

        tmpComponentList.append(Component(tmpTermList,LIN_RAD_POL))
        tmpComponentList[-1].write_latex_equation()

        ######## QUAD_AX_XXX ############
        tmpTermList      = [] 
        ### LIN_RAD_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/5.0 ,0,0))
        tmpSubTermList.append(Subterm(-3.0/35.0,0,1))
        tmpSubTermList.append(Subterm(-1.0/21.0,0,2))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_AZI))

        ### QUAD_RAD_XXX term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(3.0/7.0 ,0,0.5))
        tmpSubTermList.append(Subterm(5.0/21.0,0,1.5))
        tmpSubTermList.append(Subterm(5.0/33.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_XXX))

        tmpComponentList.append(Component(tmpTermList,QUAD_AX_XXX))
        tmpComponentList[-1].write_latex_equation()


        ############ LIN_RAD_AZI component #############
        tmpTermList      = [] 
        ### SOURCE_2D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/3.0,1,0))
        tmpSubTermList.append(Subterm(1.0/5.0,2,0))
        tmpSubTermList.append(Subterm(1.0/7.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_2D))

        ### LIN_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(3.0/5.0,1,0))
        tmpSubTermList.append(Subterm(3.0/7.0,2,0))
        tmpSubTermList.append(Subterm(1.0/3.0,2,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_AZI))

        ### QUAD_AX_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/3.0 ,1,0))
        tmpSubTermList.append(Subterm(-1.0/7.0 ,2,0))
        tmpSubTermList.append(Subterm(-5.0/63.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_POL))

        ### QUAD_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0    ,1,0))
        tmpSubTermList.append(Subterm(-5.0/7.0,2,0))
        tmpSubTermList.append(Subterm(-5.0/9.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_AZI))

        tmpComponentList.append(Component(tmpTermList,LIN_RAD_AZI))
        tmpComponentList[-1].write_latex_equation()

        ############ QUAD_RAD_XXX component ############
        tmpTermList      = [] 
        ### LIN_AX_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/5.0 ,1,0))
        tmpSubTermList.append(Subterm(-3.0/35.0,2,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_POL))

        ### QUAD_AX_XXX term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(3.0/7.0 ,1,0.5))
        tmpSubTermList.append(Subterm(5.0/21.0,2,0.5))
        tmpSubTermList.append(Subterm(5.0/33.0,3,0.5))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_XXX))

        tmpComponentList.append(Component(tmpTermList,QUAD_RAD_XXX))
        tmpComponentList[-1].write_latex_equation()
        
        ############ QUAD_AX_POL component #############
        tmpTermList      = [] 
        ### SOURCE_1D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/5.0,0,1))
        tmpSubTermList.append(Subterm(-1.0/7.0,0,2))
        tmpSubTermList.append(Subterm(-1.0/9.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_1D))

        ### LIN_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-3.0/5.0,0,0.5))
        tmpSubTermList.append(Subterm(-3.0/7.0,0,1.5))
        tmpSubTermList.append(Subterm(-1.0/3.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_POL))

        ### QUAD_RAD_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/7.0 ,0,1))
        tmpSubTermList.append(Subterm(5.0/63.0,0,2))
        tmpSubTermList.append(Subterm(5.0/99.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_AZI))

        ### QUAD_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(5.0/7.0 ,0,1))
        tmpSubTermList.append(Subterm(5.0/9.0 ,0,2))
        tmpSubTermList.append(Subterm(5.0/11.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_POL))

        tmpComponentList.append(Component(tmpTermList,QUAD_AX_POL))
        tmpComponentList[-1].write_latex_equation()

        ############ LIN_AX_POL component #############
        tmpTermList      = [] 
        ### SOURCE_1D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/3.0,0,0.5))
        tmpSubTermList.append(Subterm(1.0/5.0,0,1.5))
        tmpSubTermList.append(Subterm(1.0/7.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_1D))

        ### LIN_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(3.0/5.0,0,1))
        tmpSubTermList.append(Subterm(3.0/7.0,0,2))
        tmpSubTermList.append(Subterm(1.0/3.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_POL))

        ### QUAD_RAD_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/3.0 ,0,0.5))
        tmpSubTermList.append(Subterm(-1.0/7.0 ,0,1.5))
        tmpSubTermList.append(Subterm(-5.0/63.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_AZI))

        ### QUAD_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0    ,0,0.5))
        tmpSubTermList.append(Subterm(-5.0/7.0,0,1.5))
        tmpSubTermList.append(Subterm(-5.0/9.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_POL))

        tmpComponentList.append(Component(tmpTermList,LIN_AX_POL))
        tmpComponentList[-1].write_latex_equation()

        ############ QUAD_RAD_AZI component #############
        tmpTermList      = [] 
        ### SOURCE_2D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/5.0,1,0))
        tmpSubTermList.append(Subterm(-1.0/7.0,2,0))
        tmpSubTermList.append(Subterm(-1.0/9.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_2D))

        ### LIN_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-3.0/5.0 ,0,0))
        tmpSubTermList.append(Subterm(-3.0/7.0 ,1,0))
        tmpSubTermList.append(Subterm(-1.0/3.0 ,2,0))
        tmpSubTermList.append(Subterm(-3.0/11.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_AZI))

        ### QUAD_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(5.0/7.0 ,1,0))
        tmpSubTermList.append(Subterm(5.0/9.0 ,2,0))
        tmpSubTermList.append(Subterm(5.0/11.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_AZI))

        ### QUAD_AX_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(5.0/35.0,1,0))
        tmpSubTermList.append(Subterm(5.0/63.0,2,0))
        tmpSubTermList.append(Subterm(5.0/99.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_POL))

        tmpComponentList.append(Component(tmpTermList,QUAD_RAD_AZI))
        tmpComponentList[-1].write_latex_equation()

        ############ QUAD_RAD_POL component #############
        tmpTermList      = [] 
        ### SOURCE_2D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/15.0,1,0))
        tmpSubTermList.append(Subterm(-1.0/35.0,2,0))
        tmpSubTermList.append(Subterm(-1.0/63.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_2D))

        ### LIN_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/5.0 ,0,0))
        tmpSubTermList.append(Subterm(-3.0/35.0,1,0))
        tmpSubTermList.append(Subterm(-1.0/21.0,2,0))
        tmpSubTermList.append(Subterm(-1.0/33.0,2,0))
        tmpTermList.append(Term(tmpSubTermList,LIN_AX_AZI))

        ### QUAD_AX_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(5.0/35.0,1,0))
        tmpSubTermList.append(Subterm(5.0/63.0,2,0))
        tmpSubTermList.append(Subterm(5.0/99.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_AZI))

        ### QUAD_AX_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(5.0/35.0 ,1,0))
        tmpSubTermList.append(Subterm(5.0/105.0,2,0))
        tmpSubTermList.append(Subterm(5.0/126.0,3,0))
        tmpTermList.append(Term(tmpSubTermList,QUAD_AX_POL))

        tmpComponentList.append(Component(tmpTermList,QUAD_RAD_POL))
        tmpComponentList[-1].write_latex_equation()

        ############ QUAD_AX_AZI component #############
        tmpTermList      = [] 
        ### SOURCE_1D term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/15.0,0,1))
        tmpSubTermList.append(Subterm(-1.0/35.0,0,2))
        tmpSubTermList.append(Subterm(-1.0/63.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_1D))

        ### LIN_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0/5.0 ,0,0.5))
        tmpSubTermList.append(Subterm(-3.0/35.0,0,1.5))
        tmpSubTermList.append(Subterm(-1.0/21.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_POL))

        ### QUAD_RAD_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(5.0/35.0 ,0,1))
        tmpSubTermList.append(Subterm(1.0/21.0 ,0,2))
        tmpSubTermList.append(Subterm(5.0/126.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_AZI))

        ### QUAD_RAD_POL term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(5.0/35.0,0,1))
        tmpSubTermList.append(Subterm(5.0/63.0,0,2))
        tmpSubTermList.append(Subterm(5.0/99.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_POL))

        tmpComponentList.append(Component(tmpTermList,QUAD_AX_AZI))
        tmpComponentList[-1].write_latex_equation()

        ############ LIN_AX_AZI component ############
        ############ with halfLR included!!!! ########
        tmpTermList      = [] 
        ### LIN_RAD_AZI term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(1.0/5.0 ,0,1))
        tmpSubTermList.append(Subterm(3.0/35.0,0,2))
        tmpSubTermList.append(Subterm(1.0/21.0,0,3))
        tmpTermList.append(Term(tmpSubTermList,LIN_RAD_AZI))

        ### QUAD_RAD_XXX term
        tmpSubTermList   = [] 
        tmpSubTermList.append(Subterm(-1.0     ,0,0.5))
        tmpSubTermList.append(Subterm(-3.0/7.0 ,0,1.5))
        tmpSubTermList.append(Subterm(-5.0/21.0,0,2.5))
        tmpTermList.append(Term(tmpSubTermList,QUAD_RAD_XXX))

        tmpComponentList.append(Component(tmpTermList,LIN_AX_AZI))
        tmpComponentList[-1].write_latex_equation()

        myComponentList = ComponentList(tmpComponentList)

        ########## print all of the components #########
        #print "PRINTING ALL OF THE COMPONENTS NOW!"
        #for myComponent in myComponentList.component_list:
        #    print " "
        #    myComponent.print_component()
        #    print " "

        #############################################################
        ######### Substitution and solution is performed  ###########
        ######### in this section of the code. The method ###########
        ######### for solving the equations is outlined.  ###########
        #############################################################

        ### solution process / algorithm / steps
        ### In indexing, we have the following dependencies:
        ### Components (1,5,6) all depend on (4,8,9)
        ### Components (2,7) both depend on (3,10)

        ### Components (3,8,9) all depend on (2,5,6)
        ### Components (4,10) both depend on (1,7)

        ### Start by taking equations for (1,5,6). Subtitute in the 
        ### expressions for (4,8,9) into each of these equations.
        ### Now, (1,5,6) depend on (1,2,5,6,7).
        ### Take the equations for (2,7). Substitute in expressions for (3,10).
        ### Now, (2,7) depend on (1,2,5,6,7)

        ### Now, we have a 5x5 system of (1,2,5,6,7) in terms of (1,2,5,6,7)
        ### Solve this system of equations so that none of the higher moments 
        ### (1,2,5,6,7) depend on one another. Only dependnece will be on 
        ### SOURCE_1D,SOURCE_2D, ISO_AX_TL.

        ### Substitute in these new expressions for (2,5,6) into (3,8,9), and 
        ### substitute (1,7) into equations for (4,10).

        ### Now, all components 1-10 are known in terms of SOURCE_1D, SOURCE_2D, ISO_AX_TL
        ### Substitute the necessary components into equation for SCALAR_FLUX
        ### solve this equation for F

        ### 1) solve for LIN_AX_AZI
        linaxaziTL_comp     = myComponentList.find_component(  LIN_AX_AZI)
        linradaziTL_comp    = myComponentList.find_component( LIN_RAD_AZI)
        quadradcrossTL_comp = myComponentList.find_component(QUAD_RAD_XXX)

        #print "LIN_AX_AZI:"
        #linaxaziTL_comp.print_component()
        #print "LIN_RAD_AZI:"
        #linradaziTL_comp.print_component()
        #print "QUAD_RAD_CROSS:"
        #quadradcrossTL_comp.print_component()
        #print "This part is over"

        linaxaziTL_comp.substitute_component(linradaziTL_comp)
        linaxaziTL_comp.substitute_component(quadradcrossTL_comp)
        ### 
        linaxaziTL_comp.solve()
        linaxaziTL_comp.write_latex_equation()

        ### 2) solve for QUAD_AX_XXX
        quadaxcrossTL_comp  = myComponentList.find_component(QUAD_AX_XXX)

        quadaxcrossTL_comp.substitute_component(linradaziTL_comp)
        quadaxcrossTL_comp.substitute_component(quadradcrossTL_comp)
        quadaxcrossTL_comp.solve()
        quadaxcrossTL_comp.write_latex_equation()

        ### 3) solve for LIN_AX_POL
        linaxpolTL_comp     = myComponentList.find_component(  LIN_AX_POL)
        linradpolTL_comp    = myComponentList.find_component( LIN_RAD_POL)
        quadradaziTL_comp   = myComponentList.find_component(QUAD_RAD_AZI)
        quadradpolTL_comp   = myComponentList.find_component(QUAD_RAD_POL)

        linaxpolTL_comp.substitute_component(linradpolTL_comp)
        linaxpolTL_comp.substitute_component(quadradaziTL_comp)
        linaxpolTL_comp.substitute_component(quadradpolTL_comp)
        linaxpolTL_comp.solve()
        linaxpolTL_comp.write_latex_equation()

        ### 4) solve for QUAD_AX_POL
        quadaxpolTL_comp     = myComponentList.find_component( QUAD_AX_POL)

        quadaxpolTL_comp.substitute_component(linradpolTL_comp)
        quadaxpolTL_comp.substitute_component(quadradaziTL_comp)
        quadaxpolTL_comp.substitute_component(quadradpolTL_comp)
        quadaxpolTL_comp.solve()
        quadaxpolTL_comp.write_latex_equation()

        ### 5) solve for QUAD_AX_AZI
        quadaxaziTL_comp     = myComponentList.find_component( QUAD_AX_AZI)

        quadaxaziTL_comp.substitute_component(linradpolTL_comp)
        quadaxaziTL_comp.substitute_component(quadradaziTL_comp)
        quadaxaziTL_comp.substitute_component(quadradpolTL_comp)
        quadaxaziTL_comp.solve()
        quadaxaziTL_comp.write_latex_equation()

        ### solve the whole system
        
        ### 6) substitute LIN_AX_POL into LIN_AX_AZI, QUAD_AX_POL, QUAD_AX_AZI, QUAD_AX_XXX
        linaxaziTL_comp.substitute_component(   linaxpolTL_comp)
        quadaxpolTL_comp.substitute_component(  linaxpolTL_comp)
        quadaxaziTL_comp.substitute_component(  linaxpolTL_comp)
        quadaxcrossTL_comp.substitute_component(linaxpolTL_comp)

        linaxaziTL_comp.solve()
        ### 7) substitute LIN_AX_AZI into QUAD_AX_POL, QUAD_AX_AZI, QUAD_AX_XXX
        quadaxpolTL_comp.substitute_component(  linaxaziTL_comp)
        quadaxaziTL_comp.substitute_component(  linaxaziTL_comp)
        quadaxcrossTL_comp.substitute_component(linaxaziTL_comp)

        quadaxpolTL_comp.solve()
        ### 8) substitute QUAD_AX_POL into QUAD_AX_AZI, QUAD_AX_XXX
        quadaxaziTL_comp.substitute_component(  quadaxpolTL_comp)
        quadaxcrossTL_comp.substitute_component(quadaxpolTL_comp)

        quadaxaziTL_comp.solve()
        ### 9) substitute QUAD_AX_AZI into QUAD_AX_XXX
        quadaxcrossTL_comp.substitute_component(quadaxaziTL_comp)

        quadaxcrossTL_comp.solve()
        quadaxcrossTL_comp.write_latex_equation()
        ### 10) solve for QUAD_AX_XXX and substitute into QUAD_AX_AZI
        quadaxaziTL_comp.substitute_component(quadaxcrossTL_comp) 
        quadaxaziTL_comp.solve()
        quadaxcrossTL_comp.write_latex_equation()
        
        ### 11) substitute QUAD_AX_AZI and QUAD_AX_XXX into QUAD_AX_POL
        quadaxpolTL_comp.substitute_component(quadaxcrossTL_comp) 
        quadaxpolTL_comp.substitute_component(  quadaxaziTL_comp) 
        quadaxpolTL_comp.solve()
        quadaxpolTL_comp.write_latex_equation()
        
        ### 12) substitute QUAD_AX_POL, QUAD_AX_AZI and QUAD_AX_XXX into LIN_AX_AZI
        linaxaziTL_comp.substitute_component(quadaxcrossTL_comp) 
        linaxaziTL_comp.substitute_component(  quadaxaziTL_comp) 
        linaxaziTL_comp.substitute_component(  quadaxpolTL_comp) 
        linaxaziTL_comp.solve()
        linaxaziTL_comp.write_latex_equation()
        
        ### 13) substitute LIN_AX_AZI, QUAD_AX_POL, QUAD_AX_AZI and QUAD_AX_XXX into LIN_AX_POL
        linaxpolTL_comp.substitute_component(quadaxcrossTL_comp) 
        linaxpolTL_comp.substitute_component(  quadaxaziTL_comp) 
        linaxpolTL_comp.substitute_component(  quadaxpolTL_comp) 
        linaxpolTL_comp.substitute_component(   linaxaziTL_comp) 
        linaxpolTL_comp.solve()
        linaxpolTL_comp.write_latex_equation()

        ### 14) Substitute LIN_AX_AZI, QUAD_AX_POL, and QUAD_AX_AZI into LIN_RAD_AZI
        linradaziTL_comp.substitute_component(linaxaziTL_comp)
        linradaziTL_comp.substitute_component(quadaxpolTL_comp)
        linradaziTL_comp.substitute_component(quadaxaziTL_comp)
        linradaziTL_comp.write_latex_equation()

        ### 15) Substitute LIN_AX_AZI, QUAD_AX_POL, and QUAD_AX_AZI into QUAD_RAD_AZI
        quadradaziTL_comp.substitute_component(linaxaziTL_comp)
        quadradaziTL_comp.substitute_component(quadaxpolTL_comp)
        quadradaziTL_comp.substitute_component(quadaxaziTL_comp)
        quadradaziTL_comp.write_latex_equation()

        ### 16) Substitute LIN_AX_AZI, QUAD_AX_POL, and QUAD_AX_AZI into QUAD_RAD_POL
        quadradpolTL_comp.substitute_component(linaxaziTL_comp)
        quadradpolTL_comp.substitute_component(quadaxpolTL_comp)
        quadradpolTL_comp.substitute_component(quadaxaziTL_comp)
        quadradpolTL_comp.write_latex_equation()

        ### 17) Substitute LIN_AX_POL, QUAD_AX_XXX into LIN_RAD_POL
        linradpolTL_comp.substitute_component(linaxpolTL_comp)
        linradpolTL_comp.substitute_component(quadaxcrossTL_comp)
        linradpolTL_comp.write_latex_equation()

        ### 18) Substitute LIN_AX_POL, QUAD_AX_XXX into QUAD_RAD_XXX
        quadradcrossTL_comp.substitute_component(linaxpolTL_comp)
        quadradcrossTL_comp.substitute_component(quadaxcrossTL_comp)
        quadradcrossTL_comp.write_latex_equation()

        ### 19) Substitute LIN_AX_AZI, QUAD_AX_POL, QUAD_AX_AZI into SCALAR_FLUX
        scalarFlux_comp = myComponentList.find_component(SCALAR_FLUX)
        scalarFlux_comp.substitute_component( linaxaziTL_comp)
        scalarFlux_comp.substitute_component(quadaxpolTL_comp)
        scalarFlux_comp.substitute_component(quadaxaziTL_comp)
        scalarFlux_comp.write_latex_equation()

        ### 20) SOURCE_1D term (F - \RTL) = (Sigma_t \phi + dJ/dz)
        ########### ISO_RAD_TL component (dJ/dx + dJ/dy) ##########
        tmpTermList      = [] 
        tmpSubTermList      = [] 
        tmpSubTermList.append(Subterm(1.0,0,0))
        tmpTermList.append(Term(tmpSubTermList,SCALAR_FLUX))
        tmpTermList.append(Term(tmpSubTermList,ISO_AX_TL))
        myComponentList.add_component(Component(tmpTermList,SOURCE_1D))

        ### 21) Substitute ISO_RAD_TL into SCALAR_FLUX
        #isoradTL_comp = myComponentList.find_component(ISO_RAD_TL)
        #"print isoradTL_comp:"
        #isoradTL_comp.print_component()
        #scalarFlux_comp.substitute_component(isoradTL_comp)
        #scalarFlux_comp.solve()
        #scalarFlux_comp.write_latex_equation()
        source1D_comp = myComponentList.find_component(SOURCE_1D)
        "print isoradTL_comp:"
        source1D_comp.print_component()
        scalarFlux_comp.substitute_component(source1D_comp)
        scalarFlux_comp.solve()
        scalarFlux_comp.write_latex_equation()

        ### 22) SOURCE_2D = SOURCE_3D - ISO_AX_TL
        tmpTermList      = [] 
        tmpSubTermList      = [] 
        tmpSubTermList.append(Subterm(1.0,0,0))
        tmpTermList.append(Term(tmpSubTermList,SOURCE_3D))
        tmpSubTermList      = [] 
        tmpSubTermList.append(Subterm(-1.0,0,0))
        tmpTermList.append(Term(tmpSubTermList,ISO_AX_TL))
        myComponentList.add_component(Component(tmpTermList,SOURCE_2D))
        source2D_comp = myComponentList.find_component(SOURCE_2D)

        ### 23) Substitute SOURCE_2D into SCALAR_FLUX
        scalarFlux_comp.substitute_component(source2D_comp)
        scalarFlux_comp.write_latex_equation()

        ### 24) Solve ISO_AX_TL (substitute LIN_RAD_AZI, QUAD_RAD_AZI, QUAD_RAD_POL)
        ###     substitute SOURCE_2D = SOURCE_3D - ISO_AX_TL to solve
        isoaxTL_comp = myComponentList.find_component(ISO_AX_TL)
        isoaxTL_comp.substitute_component( linradpolTL_comp)
        isoaxTL_comp.substitute_component(quadradaziTL_comp)
        isoaxTL_comp.substitute_component(quadradpolTL_comp)

        isoaxTL_comp.substitute_component(source1D_comp)
        isoaxTL_comp.substitute_component(source2D_comp)
        isoaxTL_comp.solve()
        isoaxTL_comp.write_latex_equation()
  
        ### 25) Substitute ISO_AX_TL into SCALAR_FLUX
        scalarFlux_comp.substitute_component(isoaxTL_comp)
        scalarFlux_comp.solve()
        scalarFlux_comp.write_latex_equation()
        scalarFlux_comp.solve_F()
        scalarFlux_comp.write_latex_equation()

    #logname = "logfile"
    #
    #logfile = open(logname,"w")

    ##logstr = "iter: %i, keff = %.5f, kdiff = %.4e, fissdiff = %.4e \n" % (t+1,keff,kdiff,fissdiff)
    #kdiff = 0.0
    #fissdiff = 0.0
    #logstr = "iter: %i, keff = %.5f, kdiff = %.4e, fissdiff = %.4e \n" % (t+1,keff,kdiff,fissdiff)
    #logfile.write(logstr)
    ##print "iter: %i, keff = %.5f, kdiff = %.4e, fissnorm = %.4e, fissdiff = %.4e " % (t+1,keff,kdiff,fiss_norm,fissdiff)
    #t=t+1

    #logfile.close()
