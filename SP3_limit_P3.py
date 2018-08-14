### Michael Jarrett ### 3D SP3 Limit
### To calculate SP3 limit of 2D/1D method with anisotropic TL
#########################################
# This code performs algebra to determine the theoretical limit of 
# accuracy for a 2D/1D method with linear axial and radial TL moments.
# There are 4 linear moments. 
# Equations for each  in terms of the others are determined 
# integration by quadrature rule. There will be 
# We are actually analyzing a 1D/1D system here, with 1D radial dependence and 
# 1D axial dependence. Thus, it is a x-z or r-z system, not an (xy)-z system.
# Effectively, we have assumed that the problem is infinite in y and all y 
# derivatives are zero. This simplifies the analysis.
# With this assumption, there are still a number of coefficients. For both the 
# axial and radial TL, there will be an azimuthally isotropic moment, a cosine 
# moment. The sine moments will be zero. 
# Because cos(2x) = cos(x)^2 - sin(x)^2, the cosine moments contain both cosine 
# and sine components. cos(4x) = cos(2x)^2 - sin(2x)^2 = cos(x)^4 + sin(x)^4 - 4 * sin(x)^2 * cos(x)^2
#########################################

import numpy as np
import math
import sys
import time
import matplotlib.pyplot as plt
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
# The component definitions are hardcoded as infinite expansions and angular integrals. The value of the 
# integrals is determined by quadrature rule, using 128 polar angles and 64 azimuthal angles.

# the solution order is hardcoded based on examination of the form of the system of equations.
# the solution order is general, meaning that it works for any order of TL expansion Pn, Fn

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

    def remove_component(self,thatID):
        ic = 0
        while(ic < self.ncomponents()):
            if(self.component_list[ic].CID == thatID):
                del self.component_list[ic]
            else:
                ic += 1

    def find_component(self,thatID):
        for myComponent in self.component_list:
            if(myComponent.CID == thatID):
                return myComponent
        return -1

    def print_dependence_matrix(self,suffix):
        ### print out a matrix that shows the dependence of all of the terms
        ### the matrix should be a square matrix of size NxN, where N = total number of moments
        suffix = str(suffix)

        #total_ID = nID + 2*nmom*nfourier
        total_ID = self.ncomponents()
        dependence_matrix = np.zeros([total_ID+1,total_ID+1])
        for idx in range(0,total_ID+1):
            print "idx = %i" % idx
            thisComponent = self.find_component(idx)
            if(thisComponent == -1): 
                # cycle
                print "no component for this term. idx = %i" % idx
            else:
                for term in thisComponent.term_list: 
                    dependence_matrix[idx,term.TID] = 1
    
        fig = plt.figure(1)
        plt.spy(dependence_matrix)
        #plt.show()
        filetype = 'png'
        if(filetype == 'eps'):
            filename = "matrix_dependence_P%i_F%i_%s.eps" % (nmom-1,nfourier-1,suffix)
            savename = "./%s" % (filename)
            fig.savefig(savename, format='eps', dpi=300)
            plt.close(fig)
        elif(filetype == 'png'):
            filename = "matrix_dependence_P%i_F%i_%s.png" % (nmom-1,nfourier-1,suffix)
            savename = "./%s" % (filename)
            fig.savefig(savename, format='png')
            plt.close(fig)
            

class Component:
    def __init__(self, term_list, myTL_ID):
        self.term_list = term_list
        self.CID = myTL_ID
        thisMu = muTerm(0,0,[])
        thisOmega = omegaTerm([],[])
        self.LHS = Term([Subterm(1.0,0,0,thisMu,thisOmega)], myTL_ID)
        self.angles_are_integrated = False

    def __del__(self):
        self.term_list = []
        self.CID = 0

    @property
    def CID(self):
        return self.CID
    @CID.setter
    def CID(self,ID):
        self.CID = ID

    def nterms(self):
         return len(self.term_list)

    def substitute_component(self,thatComponent):
        # substitute the expression for that Component 
        # where it appears in this Component
        itx = 0 
        newTermList = [] 
        while(itx < self.nterms()):
            tmpTerm = self.term_list[itx]
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
            tmpTerm1.combine_subterms()
            #tmpTerm1.remove_halves()

    def remove_high_order(self):
        for myterm in self.term_list:
            myterm.remove_high_order()
            #myterm.remove_halves()

    # solve component if the term appears on left and right sides of the equation
    def solve(self):
        # first, ensure the component is "clean" for solving
        self.remove_high_order()
        self.combine_terms()

        ### move self terms to the LHS
        itx = 0
        null_muterm = muTerm(0,0,[])
        null_omegaterm = omegaTerm([],[])
        while(itx < self.nterms()):
            myTerm = self.term_list[itx]
            newSubTermList = []
            if(myTerm.TID == self.CID):
                for mySubTerm in myTerm.subterm_list:
                    newSubTermList.append(Subterm(-1.0*mySubTerm.coeff,mySubTerm.r_order,mySubTerm.z_order,copy_mu(null_muterm),copy_omega(null_omegaterm)))
                self.LHS.subterm_list.extend(newSubTermList)
                del self.term_list[itx]
            itx += 1
        self.LHS.subterm_list.extend(newSubTermList)

        tmpMu = muTerm(0,0,[])
        tmpOmega = omegaTerm([],[])

        ### move self terms to the LHS
        itx = 0
        while(itx < self.nterms()):
            myTerm = self.term_list[itx]
            newSubTermList = []
            if(myTerm.TID == self.CID):
                for mySubTerm in myTerm.subterm_list:
                    newSubTermList.append(Subterm(-1.0*mySubTerm.coeff,mySubTerm.r_order,mySubTerm.z_order,copy_mu(tmpMu),copy_omega(tmpOmega)))
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
            mySubTerm.coeff = -1.0*mySubTerm.coeff

        ### calculate x^2, x^3
        extraTerm_squared = term_mult(extraTerm,extraTerm)
        extraTerm_cubed = term_mult(extraTerm,extraTerm_squared)
         
        extraTerm.term_add(extraTerm_squared)
        extraTerm.term_add(extraTerm_cubed)

        tmpMu = muTerm(0,0,[])
        tmpOmega = omegaTerm([],[])
        extraTerm.subterm_list.insert(0,Subterm(1.0,0,0,tmpMu,tmpOmega))

        extraTerm.combine_subterms()

        ### now, extra term = I + x + x^2 + x^3
        #print "Extra term is:"
        #extraTerm.print_info()

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

        if(self.nterms() == 1 and 
           self.term_list[0].TID == SOURCE_3D and
           self.CID == SCALAR_FLUX):

            extraTerm = Term([],SCALAR_FLUX)
            FTerm = self.term_list[0]
            if(FTerm.subterm_list[0].isIdentity):
                for mySubTerm in FTerm.subterm_list[1:]:
                    extraTerm.subterm_list.append(mySubTerm)
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
            mySubTerm.coeff = -1.0*mySubTerm.coeff

        ### calculate x^2, x^3
        extraTerm_squared = term_mult(extraTerm,extraTerm)

        extraTerm_cubed = term_mult(extraTerm,extraTerm_squared)
         
        extraTerm.term_add(extraTerm_squared)
        extraTerm.term_add(extraTerm_cubed)

        tmpMu = muTerm(0,0,[])
        tmpOmega = omegaTerm([],[])
        extraTerm.subterm_list.insert(0,Subterm(1.0,0,0,copy_mu(tmpMu),copy_omega(tmpOmega)))

        ### now, extra term = I + x + x^2 + x^3
        ### combine subterms
        extraTerm.combine_subterms()
        ### sort subterms
        extraTerm.sort_subterms()
        #print "Extra term is:"
        #extraTerm.print_info()

        newComponent = Component([extraTerm],SOURCE_3D)
        print "The solution is:"
        newComponent.print_component()
        newComponent.write_latex_equation()

    def integrate_angles(self):
        if(self.angles_are_integrated):
            print "Angles should only be integrated once!"
        else:
            ### iterate through terms in the component, subterms in the term
            for term in self.term_list:
                for subterm in term.subterm_list:
                    subterm.integrate_mu()
                    subterm.integrate_omega()

            itx = 0
            while(itx < self.nterms()):
                ### clear out zeroes
                ist = 0
                while(ist < self.term_list[itx].nsubterms()):
                    mycoeff = self.term_list[itx].subterm_list[ist].coeff
                    if(abs(mycoeff) < THRESHOLD):
                        del self.term_list[itx].subterm_list[ist]
                    else:
                        ist += 1
                
                ### delete term if no subterms remain
                if(self.term_list[itx].nsubterms() == 0):
                    del self.term_list[itx]
                else:
                    itx += 1

            self.angles_are_integrated = True

    def reset_LHS(self):
        self.LHS.subterm_list = []
        tmpMu = muTerm(0,0,[])
        tmpOmega = omegaTerm([],[])
        self.LHS = Term([Subterm(1.0,0,0,copy_mu(tmpMu),copy_omega(tmpOmega))],0)


    def print_info(self):
        print "My component terms: ID = %i (%s)" % (self.CID,variable_names[self.CID])
        for myterm in self.term_list:
            myterm.print_info()

    def print_component(self):
        if(self.CID <= nID):
            print "My component terms: ID = %i (%s)" % (self.CID,variable_names[self.CID])
        else: # this is a TL moment
            [dimension,ileg,ifourier] = getMomOrder(self.CID)
            if(dimension == AXIAL_TL):
                print "AXIAL TL COMPONENT: l = %i, p = %i" % (ileg,ifourier)
            elif(dimension == RADIAL_TL):
                print "RADIAL TL COMPONENT: l = %i, p = %i" % (ileg,ifourier)
        for myterm in self.term_list:
            myterm.print_term()

    def write_latex_equation(self):
        texfile = open(texname,"a")

        texstr = "\\begin{alignat}{2} \n"
        texfile.write(texstr)
        if(self.CID <= nID):
            texstr = "\\label{eq:%s_%i} \n" % (equation_labels[self.CID],equation_counter.get_index(self.CID))
        else:
            [dimension,ileg,ifourier] = getMomOrder(self.CID)
            if(dimension == AXIAL_TL):
                texstr = "\\label{eq:f_%i%i_%i} \n" % (ileg,ifourier,equation_counter.get_index(self.CID))
            elif(dimension == RADIAL_TL):
                texstr = "\\label{eq:g_%i%i_%i} \n" % (ileg,ifourier,equation_counter.get_index(self.CID))

        texfile.write(texstr)
        texstr = "& %s = \\\\ \n" % latex_expr(self.CID)
        texfile.write(texstr)

        for myTerm in self.term_list:
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

        self.remove_zero()

    def remove_zero(self):
        ist = 0
        while (ist < self.nsubterms()):
             tmpSubTerm = self.subterm_list[ist]
             if(abs(tmpSubTerm.coeff) < THRESHOLD):
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
        if(self.TID <= nID):
            print "Term variable ID: %i (%s)" % (self.TID,variable_names[self.TID])
        else: # this is a TL moment
            [dimension,ileg,ifourier] = getMomOrder(self.TID)
            if(dimension == AXIAL_TL):
                print "Axial TL: l = %i, p = %i" % (ileg,ifourier)
            elif(dimension == RADIAL_TL):
                print "Radial TL: l = %i, p = %i" % (ileg,ifourier)
        print "( ",
        for mySubTerm in self.subterm_list:
            if(mySubTerm.coeff > 0.0):
                print "+ ",
            else:
                print "",
            print "%7.5f" % (mySubTerm.coeff),
            if(mySubTerm.r_order > 0):
                if(mySubTerm.r_order == 1):
                    print "Lr",
                else:
                    print "Lr^%3.1f" % mySubTerm.r_order,
            if(mySubTerm.z_order > 0):
                if(mySubTerm.z_order == 1):
                    print "Lz",
                else:
                    print "Lz^%3.1f" % mySubTerm.z_order,

            ### mu terms
            if(mySubTerm.mu_term.cos_order > 0):
                print "u^%i" % mySubTerm.mu_term.cos_order,
            if(mySubTerm.mu_term.sin_order > 0):
                print "sqrt(1-u^2)^%i" % mySubTerm.mu_term.sin_order,
            if(mySubTerm.mu_term.nlegterms > 0):
                for ileg in mySubTerm.mu_term.leg_list:
                    print "P_%i(u)" % ileg,

            ### omega terms
            myWterm = mySubTerm.omega_term
            for icos in myWterm.cos_list:
                print "cos(%iw)^%i" % (icos[0],icos[1])
            for isin in myWterm.sin_list:
                print "sin(%iw)^%i" % (isin[0],isin[1])

        print ") "

    def write_term_latex(self):
        tmpstr = []

        ### determine if whole term is negative
        allneg = True
        for mySubTerm in self.subterm_list:
            if(mySubTerm.coeff > 0.0):
                allneg = False

        if(allneg):
            tmpstr.append('& -')
        else:
            tmpstr.append('& +')

        tmpstr.append(' \\l[ ')
        firstTerm = True
        for mySubTerm in self.subterm_list:

            if(mySubTerm.coeff == 1.0):
                mystr = ""
            elif(mySubTerm.coeff == -1.0):
                if(allneg == False):
                    mystr = "- "
                else:
                    mystr = ""
            else:
                ifrac = compare_fraction(mySubTerm.coeff)
                if(ifrac >= 0):
                    if(mySubTerm.coeff > 0.0 or allneg):
                        if(not firstTerm): 
                            tmpstr.append('+ ')
                    else:
                        if(not firstTerm):
                            tmpstr.append('- ')

                    thisPair = common_fractions[ifrac] 
                    mystr = '\\frac{%i}{%i} ' % (thisPair[0],thisPair[1])
                else: # no match, write a floating-point coefficient
                    if(mySubTerm.coeff > 0.0 or allneg):
                        if(not firstTerm):
                            tmpstr.append('+ ')
                    else:
                        if(not firstTerm):
                            tmpstr.append(' ')
                    if(allneg):
                        mystr = ' %7.5f ' % (abs(mySubTerm.coeff))
                    else:
                        mystr = ' %7.5f ' % (mySubTerm.coeff)
            tmpstr.append(mystr)

            if(mySubTerm.r_order > 0):
                if(mySubTerm.r_order == 1):
                    tmpstr.append('\\opL_r ')
                else:
                    mystr = '\\opL_r^{%3.1f} ' % mySubTerm.r_order
                    tmpstr.append(mystr)

            #intzorder = int(mySubTerm.z_order)
            if(mySubTerm.z_order > 0):
                if(mySubTerm.z_order == 1):
                    tmpstr.append('\\opL_z ')
                else:
                    mystr = '\\opL_z^{%3.1f} ' % mySubTerm.z_order 
                    tmpstr.append(mystr)

            if(mySubTerm.r_order == 0 and mySubTerm.z_order == 0):
                mystr = "I "
                tmpstr.append(mystr)
            
            firstTerm = False
                
        tmpstr.append(' \\r] ')

        #if((2*mySubTerm.z_order)%2 == 1): # half LZ is present
        #    tmpstr.append('\halfLZ ') 

        #if(self.TID < len(latex_expressions)):
        #    tmpstr.append(latex_expressions[self.TID])
        #else:
        #    if(self.TID <= nID + (nmom*nfourier)):
        #        [dimension,ileg,ifourier] = getMomOrder(self.TID)
        #        if(dimension == AXIAL_TL):
        #            if(ifourier == 0):
        #                momstr = "f_{%i,0}" % (ileg)
        #            else:
        #                momstr = "f_{c,%i,%i}" % (ileg,ifourier)
        #        elif(dimension == RADIAL_TL):
        #            if(ifourier == 0):
        #                momstr = "g_{%i,0}" % (ileg)
        #            else:
        #                momstr = "g_{c,%i,%i}" % (ileg,ifourier)
        #    tmpstr.append(momstr)
        #tmpstr.append(' \\nonumber \\\\ \n')

        momstr = latex_expr(self.TID)
        tmpstr.append(momstr)

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
    def __init__(self, coeff, r_order, z_order, muTerm, omegaTerm):
        self.coeff = coeff
        self.r_order = r_order
        self.z_order = z_order
        self.mu_term = muTerm
        self.omega_term = omegaTerm

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
            oldcoeff = float(self.coeff)
            self.coeff = self.coeff + thatSubTerm.coeff 

            have_frac = compare_fraction(self.coeff)
            if((self.r_order + self.z_order) <= MAX_ORDER): # we need to find a fraction
                if(have_frac == -1): # we didn't find the fraction! 
                    if(abs(self.coeff) >= THRESHOLD):
                        #print "Didn't find the fraction in our list! The float %8.6e is the sum of " % abs(self.coeff)
                        #print "%16.14e and %16.14e" % (abs(oldcoeff),abs(thatSubTerm.coeff))
                        i = compare_fraction(oldcoeff)
                        j = compare_fraction(thatSubTerm.coeff)
                        if(i >= 0):
                            myPair1 = common_fractions[i]
                        #else:
                        #    print "did not find fraction for coeff 1"
                        if(j >= 0):
                            myPair2 = common_fractions[j]
                        #else:
                        #    print "did not find fraction for coeff 2"
                        #if(i >= 0 and j >= 0): 
                        #    print "%i/%i + %i/%i = %i/%i " %  (myPair1[0],myPair1[1],myPair2[0],
                        #                 myPair2[1],myPair1[0]*myPair2[1]+myPair2[0]*myPair1[1],myPair1[1]*myPair2[1])

                        ### search for sum of two common fractions
                        #[frac1,frac2] = compare_fraction_sum(coeff)
                        #if(frac1 >= 0): # we found a match
                        #    print "The float is the sum of two fractions:"
                        #    
                        #    myPair1 = common_fractions[frac1]
                        #    myPair2 = common_fractions[frac2]
                        #    print "%i/%i + %i/%i = %i/%i " %  (myPair1[0],myPair1[1],myPair2[0],
                        #                 myPair2[1],myPair1[0]*myPair2[1]+myPair2[0]*myPair1[1],myPair1[1]*myPair2[1])
 
                


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

    def integrate_mu(self):
        muint = self.mu_term.eval_quad_mu()
        self.coeff = self.coeff*muint
        null_muterm = muTerm(0,0,[])
        self.mu_term = null_muterm

    def integrate_omega(self):
        omegaint = self.omega_term.eval_quad_omega()
        self.coeff = self.coeff*omegaint
        null_omegaterm = omegaTerm([],[])
        self.omega_term = null_omegaterm
        
    def print_info(self):
        print "Subterm: coeff = %.5f,Lr = %3.1f, Lz = %3.1f" % (self.coeff,self.r_order,self.z_order),
        print "u^{%i}, q(1-u^2)^{%i}," % (self.mu_term.cos_order,self.mu_term.sin_order),
        for ileg in self.mu_term.leg_list:
            print "P{%i}(u)" % ileg,
        for icos in self.omega_term.cos_list:
            print "cos(%iw)^{%i}" % (icos[0],icos[1]),
        for isin in self.omega_term.sin_list:
            print "sin(%iw)%{%i}" % (isin[0],isin[1]),
        print ""

class muTerm:
    def __init__(self,cos_order,sin_order,legendre_list):
        self.cos_order = cos_order
        self.sin_order = sin_order
        self.leg_list = legendre_list

    def __del__(self):
        self.cos_order = 0
        self.sin_order = 0
        self.leg_list = []

    @property
    def cos_order(self):
        return self.cos_order
    @cos_order.setter
    def cos_order(self,val):
        self.cos_order = val

    @property
    def sin_order(self):
        return self.cos_order
    @sin_order.setter
    def sin_order(self,val):
        self.cos_order = val

    def nlegterms(self):
        return len(self.leg_list)

    def eval_quad_mu(self):
        f = 0.0
        for ipol in range(0,npol): 
            temp = 1.0
            temp = temp*sinmu[ipol]**self.sin_order
            temp = temp*mu[ipol]**self.cos_order
            for ileg in self.leg_list:
                temp = temp*legendre[ileg,ipol]
            f += temp*wtpol[ipol]
        return f

class omegaTerm:
    def __init__(self,cos_list,sin_list):
        self.cos_list = cos_list
        self.sin_list = sin_list

    def __del__(self):
        self.cos_list = []
        self.sin_list = []

    def eval_quad_omega(self):
        f = 0.0
        for iazi in range(0,nazi): 
            tmpazi = angazi[iazi]
            temp = 1.0
            for icos in self.cos_list:
                frequency = icos[0]
                order = icos[1]
                temp = temp*math.cos(frequency*tmpazi)**order

            for isin in self.sin_list:
                frequency = isin[0]
                order = isin[1]
                temp = temp*math.sin(frequency*tmpazi)**order

            f += temp*wtazi
        return f

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

    newSubTermList = []
    for mySubTerm in myTerm.subterm_list:
        for thatSubTerm in thatTerm.subterm_list:
            newSubTerm = subterm_mult(mySubTerm,thatSubTerm)
            newSubTermList.append(newSubTerm)
    
    newTerm = Term(newSubTermList,myTerm.TID)
    newTerm.remove_high_order()
    newTerm.combine_subterms()
    newTerm.sort_subterms()

    return newTerm

# multiply two subterms
def subterm_mult(mySubTerm,thatSubTerm):
    coeff = mySubTerm.coeff * thatSubTerm.coeff
    r_order = mySubTerm.r_order + thatSubTerm.r_order
    z_order = mySubTerm.z_order + thatSubTerm.z_order
    have_frac = compare_fraction(coeff)
    if((r_order + z_order) <= MAX_ORDER): # we need to find a fraction
        if(have_frac == -1): # we didn't find the fraction! 
            if(abs(coeff) >= THRESHOLD):
                #print "Didn't find the fraction in our list! The float %8.6e is the product of " % abs(coeff)
                #print "%16.14e and %16.14e" % (abs(mySubTerm.coeff),abs(thatSubTerm.coeff))
                i = compare_fraction(mySubTerm.coeff)
                j = compare_fraction(thatSubTerm.coeff)
                if(i >= 0):
                    myPair1 = common_fractions[i]
                #else:
                #    print "did not find fraction for coeff 1"
                if(j >= 0):
                    myPair2 = common_fractions[j]
                #else:
                #    print "did not find fraction for coeff 2"
                #if(i >= 0 and j >= 0): 
                #    print "%i/%i x %i/%i = %i/%i " %  (myPair1[0],myPair1[1],myPair2[0],
                #                 myPair2[1],myPair1[0]*myPair2[0],myPair1[1]*myPair2[1])

                ### search for sum of two common fractions
                #[frac1,frac2] = compare_fraction_sum(coeff)
                #if(frac1 >= 0): # we found a match
                #    #print "The float is the sum of two fractions:"
                #    
                #    myPair1 = common_fractions[frac1]
                #    myPair2 = common_fractions[frac2]
                #    #print "%i/%i + %i/%i = %i/%i " %  (myPair1[0],myPair1[1],myPair2[0],
                #    #             myPair2[1],myPair1[0]*myPair2[1]+myPair2[0]*myPair1[1],myPair1[1]*myPair2[1])
 
                
    cos_order = mySubTerm.mu_term.cos_order + thatSubTerm.mu_term.cos_order
    sin_order = mySubTerm.mu_term.sin_order + thatSubTerm.mu_term.sin_order
    leg_list = []
    leg_list.extend(mySubTerm.mu_term.leg_list)
    leg_list.extend(thatSubTerm.mu_term.leg_list)

    cos_list = []
    cos_list.extend(mySubTerm.omega_term.cos_list)
    cos_list.extend(thatSubTerm.omega_term.cos_list)

    sin_list = []
    sin_list.extend(mySubTerm.omega_term.sin_list)
    sin_list.extend(thatSubTerm.omega_term.sin_list)

    tmpMuTerm = muTerm(cos_order,sin_order,leg_list)
    tmpOmegaTerm = omegaTerm(cos_list,sin_list)

    return Subterm(coeff,r_order,z_order,tmpMuTerm,tmpOmegaTerm)

def copy_subterm_list(thatSubTermList):
    newSubTermList = []
    for mySubTerm in thatSubTermList:
        newSubTermList.append(Subterm(mySubTerm.coeff,mySubTerm.r_order,mySubTerm.z_order) )
    return newSubTermList

def copy_subterm(thatSubTerm):
    return Subterm(thatSubTerm.coeff,thatSubTerm.r_order,thatSubTerm.z_order,thatSubTerm.mu_term,thatSubTerm.omega_term)

def copy_mu(thatMuTerm):
    return muTerm(thatMuTerm.cos_order,thatMuTerm.sin_order,thatMuTerm.leg_list)

def copy_omega(thatOmegaTerm):
    return omegaTerm(thatOmegaTerm.cos_list,thatOmegaTerm.sin_list)

def compare_fraction(thatFloat):
    for i in range(0,nfrac):
        myPair = common_fractions[i]
        myFloat = float(myPair[0])/float(myPair[1])
        if(abs(myFloat - abs(thatFloat)) < THRESHOLD*abs(thatFloat)): # this is a match
            return i           
    return -1

#def compare_fraction_sum(thatFloat):
#    for i in range(0,nfrac):
#        for j in range(0,nfrac):
#            myPair1 = common_fractions[i]
#            myPair2 = common_fractions[j]
#            myFloat1 = float(myPair1[0])/float(myPair1[1])
#            myFloat2 = float(myPair2[0])/float(myPair2[1])
#            if(abs(myFloat1 + myFloat2 - abs(thatFloat)) < THRESHOLD*abs(thatFloat)): # this is a match
#                return [i,j]           
#    return [-1,-1]

def getMomID(ileg,ifourier,dimension):
    if(dimension == AXIAL_TL):
        return nID + (nmom*ifourier) + ileg + 1
    elif(dimension == RADIAL_TL):
        return nID + (nmom*nfourier) + (nmom*ifourier) + ileg + 1
    else:
        print "Moment ID for dimension that is not AXIAL or RADIAL!!!"

def getMomOrder(IDX):
     
    if(IDX <= (nID + nmom*nfourier)):
        dimension = AXIAL_TL
        momID = IDX - nID - 1
        ileg = momID%nmom
        ifourier = ( momID - ileg ) / nmom
    elif(IDX <= (nID + 2*nmom*nfourier)):
        dimension = RADIAL_TL
        momID = IDX - nID - 1
        momID = momID - nmom*nfourier
        ileg = momID%nmom
        ifourier = ( momID - ileg ) / nmom
    else:
        print "This isn't a TL moment!"

    return [dimension,ileg,ifourier]

def latex_expr(ID):
    tmpstr = []
    if(ID < len(latex_expressions)):
        tmpstr.append(latex_expressions[ID])
    else:
        [dimension,ileg,ifourier] = getMomOrder(ID)
        if(dimension == AXIAL_TL):
            if(ifourier == 0):
                momstr = "f_{%i,0}" % (ileg)
            else:
                momstr = "f_{c,%i,%i}" % (ileg,ifourier)
        elif(dimension == RADIAL_TL):
            if(ifourier == 0):
                momstr = "g_{%i,0}" % (ileg)
            else:
                momstr = "g_{c,%i,%i}" % (ileg,ifourier)
        tmpstr.append(momstr)
    tmpstr.append(' \\nonumber \\\\ \n')

    return ''.join(tmpstr)

### indexing for each TL component 

#quadfile = "polquad_n8.txt"
#weightfile = "weightquad_n8.txt"
quadfile = "polquad_n128.txt"
weightfile = "weightquad_n128.txt"

### flux 
SCALAR_FLUX  = 1
### F: 3D source (this is the component we ultimately need to solve for)
SOURCE_3D  = 2

AXIAL_TL = 0
RADIAL_TL = 1
nID = 2

MAX_ORDER = 3
LEG_ORDER = 3
FOURIER_ORDER = 3
THRESHOLD = 1.0E-5

#if(int(sys.argv[1]) == 1):
#    texname = "equations_output_linear.tex"
#if(int(sys.argv[1]) == 2):
#    texname = "equations_output_quadratic.tex"
texname = "equations_output_p%i.tex" % (int(sys.argv[1])-1)

countmom = int(sys.argv[1])
countfourier = int(sys.argv[2])

equation_counter = Counter(nID + 2*countmom*countfourier)

# names for each TL component variable
variable_names = ['null',
                  'Scalar Flux',
                  'Source (F)']

latex_expressions = ['',
                     '\phi',
                     'F']

equation_labels   = ['',
                     'scalarFlux',
                     'transport3D']

##### commonly encountered fractions for LaTeX output
common_fractions = [(1,  3),
                    (1,  5),
                    (1,  7),
                    (1,  9),
                    (4, 45),
                    (6, 45),
                    (8, 45),
                    (9, 45),
                   (44,945),
                   (44,315)]
#                    (1, 11),
#                    (1, 15),
#                    (1, 21),
#                    (1, 25),
#                    (1, 27),
#                    (1, 33),
#                    (1, 35),
#                    (1, 42),
#                    (1, 49),
#                    (1, 63),
#                    (1, 75),
#                    (1, 99),
#                    (1,105),
#                    (1,125),
#                    (1,147),
#                    (1,175),
#                    (1,189),
#                    (1,245),
#                    (1,375),
#                    (2,  3),
#                    (2,  5),
#                    (2, 25),
#                    (2, 35),
#                    (2, 45),
#                    (2, 49),
#                    (2, 75),
#                    (2,105),
#                    (2,125),
#                    (2,135),
#                    (2,147),
#                    (2,175),
#                    (2,245),
#                    (3,  5),
#                    (3,  7),
#                    (3, 25),
#                    (3, 35),
#                    (3, 49),
#                    (3,125),
#                    (3,175),
#                    (4, 35)]
#common_fractions = [(1,  3),
#                    (1,  5),
#                    (1,  7),
#                    (1,  9),
#                    (1, 11),
#                    (1, 15),
#                    (1, 21),
#                    (1, 25),
#                    (1, 27),
#                    (1, 33),
#                    (1, 35),
#                    (1, 42),
#                    (1, 49),
#                    (1, 63),
#                    (1, 75),
#                    (1, 99),
#                    (1,105),
#                    (1,125),
#                    (1,147),
#                    (1,175),
#                    (1,189),
#                    (1,231),
#                    (1,245),
#                    (1,343),
#                    (1,375),
#                    (1,735),
#                    (1,1715),
#                    (2,  3),
#                    (2,  5),
#                    (2, 25),
#                    (2, 35),
#                    (2, 45),
#                    (2, 49),
#                    (2, 75),
#                    (2,105),
#                    (2,125),
#                    (2,135),
#                    (2,147),
#                    (2,175),
#                    (2,245),
#                    (3,  5),
#                    (3,  7),
#                    (3, 25),
#                    (3, 35),
#                    (3, 49),
#                    (3,125),
#                    (3,175),
#                    (3,245),
#                    (3,343),
#                    (3,1715),
#                    (4, 35),
#                    (4, 45),
#                    (4,105),
#                    (4,175),
#                    (4,245),
#                    (4,525),
#                    (5,  3),
#                    (5,  7),
#                    (5,  9),
#                    (5, 11),
#                    (5, 21),
#                    (5, 27),
#                    (5, 33),
#                    (5, 49),
#                    (5, 63),
#                    (5, 99),
#                    (5,126),
#                    (5,147),
#                    (5,189),
#                    (5,231),
#                    (5,343),
#                    (5,441),
#                    (5,1715),
#                    (6, 35),
#                    (6, 45),
#                    (6, 77),
#                    (6,125),
#                    (6,175),
#                    (6,245),
#                    (6,539),
#                    (6,875),
#                    (6,1125),
#                    (8, 45),
#                    (8,105),
#                    (8,175),
#                    (8,225),
#                    (8,245),
#                    (8,315),
#                    (8,525),
#                    (8,875),
#                    (8,1575),
#                    (9, 35),
#                    (9, 49),
#                    (9,125),
#                    (9,175),
#                    (9,245),
#                   (10, 49),
#                   (10,189),
#                   (11,175),
#                   (12,175),
#                   (12,875),
#                   (12,245),
#                   (15,147),
#                   (16,175),
#                   (18,245),
#                   (18,875),
#                   (19,105),
#                   (22,2625),
#                   (24, 35),
#                   (24,245),
#                   (24,875),
#                   (25,441),
#                   (29,875),
#                   (30,539),
#                   (32,2625),
#                   (37,700),
#                   (39,225),
#                   (44,945),
#                   (46,675),
#                   (46,875),
#                   (62,147),
#                   (97,147),
#                   (97,735),
#                   (18,1225),
#                   (18,4375),
#                   (27,1225),
#                   (96,4900),
#                   (52,1575),
#                   (87,7875),
#                  (161,735),
#                  (101,1575),
#                  (101,3675),
#                  (106,4725),
#                  (112,1715),
#                  (162,7875),
#                  (252,8557),
#                  (252,8575),
#                  (372,1225),
#                  (348,6125),
#                  (592,3885),
#                  (696,6125),
#                  (548,18375),
#                  (938,42875),
#                  (3304,42875),
#                  (980,7203)]

nfrac = len(common_fractions)

if __name__ == '__main__':

    global mu
    global sinmu
    global wtpol
    global legendre
    global npol
    global nfourier
    global nmom
    global nazi
    global wtazi
    global angazi

    ### load quadrature information
    #npol = 8
    npol = 128
    nmom = int(sys.argv[1])
    nfourier = int(sys.argv[2])
    mu = np.empty([npol])
    sinmu = np.empty([npol])
    wtpol = np.empty([npol])
    legendre = np.empty([nmom,npol])
    nazi = 64
    wtazi = 2.0*math.pi/float(nazi)
    angazi = np.empty([nazi])

    quad = open(weightfile)
    ipol=0
    for line in quad: 
        tempweight = float(line)
        wtpol[ipol]=tempweight
        ipol+=1
    quad.close()

    quad = open(quadfile)
    ipol=0
    for line in quad: 
        theta = float(line)
        mutemp = math.cos(theta)
        if(theta < 0.0):
            mutemp = -1.0*abs(mutemp)
        mu[ipol]=mutemp
        legendre[0,ipol]=1.0
        if(nmom > 1):
            legendre[1,ipol]=mutemp
        if(nmom > 2):
            legendre[2,ipol]=(3.0*mutemp*mutemp-1.0)/2.0
        if(nmom > 3):
            legendre[3,ipol]=(5.0*mutemp**3-3.0*mutemp)/2.0
        if(nmom > 4):
            legendre[4,ipol]=(35.0*mutemp**4-30.0*mutemp**2+3)/8.0
        if(nmom > 5):
            legendre[5,ipol]=(63.0*mutemp**5-70.0*mutemp**3+15.0*mutemp)/8.0
        if(nmom > 6):
            legendre[6,ipol]=(231.0*mutemp**6-315.0*mutemp**4+105.0*mutemp**2-5.0)/16.0
        sinmu[ipol] = math.sqrt(1.0-mutemp**2)
        ipol+=1
    quad.close()
   
    # setup azimuthal quadrature
    halfwtazi = wtazi/2.0
    angazi = np.linspace(halfwtazi,2.0*math.pi+halfwtazi,nazi+1)
    nazi_half = nazi/2

##########################################################
########## nothing should go before this in main #########
##########################################################

    ### clear the equations_output.tex file
    call(["rm",texname])
    call(["touch",texname])


    ### moment identification
    ### 4 legendre moments, 4 fourier moments
    ### define the terms, with their angle-dependence
    ### nMomID = nleg*nfourier = 16
    ### polar_momentID = nleg*ifourier+ileg
    ### azi_momentID = nMomID + nleg*ifourier+ileg
    ### 32 total moments

    null_muterm = muTerm(0,0,[])
    null_omegaterm = omegaTerm([],[])
    null_cos_list = []
    null_sin_list = []
    null_leg_list = []
    null_r_order = 0
    null_z_order = 0

    # set up component list
    tmpComponentList = [] 

    ############## SCALAR FLUX (\phi) #############
    #tmpTermList      = [] 
    ### SOURCE_2D term

    tmpSubTermList   = [] 
    ### infinite series expansion term for omega_x
    z_order = 0 
    for n in range(0,(2*MAX_ORDER+1)):
        #tmpMuTerm = copy_mu(null_muterm)
        tmpMuTerm = muTerm(0,n,[])
        cos_list =  [(1,n)]
        r_order = n/2.0
        coeff = (-1.0)**n
        thisOmegaTerm = omegaTerm(cos_list,null_sin_list)
        newSubTerm = Subterm(coeff,r_order,z_order,tmpMuTerm,thisOmegaTerm)
        tmpSubTermList.append(newSubTerm)

    ### this term does not have an associated moment ID
    infexpTerm_omegax = Term(tmpSubTermList,0) 

    tmpSubTermList   = [] 
    ### infinite series expansion term for omega_z = mu
    r_order = 0 
    for n in range(0,(2*MAX_ORDER+1)):
        tmpMuTerm = muTerm(n,0,[])
        z_order = n/2.0
        coeff = (-1.0)**n
        thisOmegaTerm = omegaTerm(null_cos_list,null_sin_list)
        newSubTerm = Subterm(coeff,r_order,z_order,tmpMuTerm,thisOmegaTerm)
        tmpSubTermList.append(newSubTerm)

    ### this term does not have an associated moment ID
    infexpTerm_omegaz = Term(tmpSubTermList,0) 

    ### the rest of the terms have only one subterm, and are associated with a moment ID 
    tmpTermList = []
    ### SOURCE_3D
    tmpTerm = Term([Subterm(1.0/(4.0*math.pi),0,0,copy_mu(null_muterm),copy_omega(null_omegaterm))],SOURCE_3D)
    tmp_dist = term_mult(tmpTerm,infexpTerm_omegax)
    tmpTermList.append(tmp_dist)

    ### axial TL term
    ### iterate through Legendre and Fourier moments, creating a term for each
    for ileg in range(0,nmom):
        tmpMuTerm = muTerm(1,0,[ileg])
        for ifourier in range(0,nfourier):
            if(ifourier == 0):
                coeff = - 1.0/(2.0*math.pi) * (2.0*ileg + 1.0) / 2.0
                tmpOmegaTerm = null_omegaterm
            else:
                coeff = - 1.0/(math.pi) * (2.0*ileg + 1.0) / 2.0
                cos_list = [(ifourier,1)]
                tmpOmegaTerm = omegaTerm(cos_list,[])

            momID = getMomID(ileg,ifourier,AXIAL_TL)
            tmpTerm = Term([Subterm(coeff,0,0.5,tmpMuTerm,tmpOmegaTerm)],momID)
            tmp_dist = term_mult(tmpTerm,infexpTerm_omegax)
            tmpTermList.append(tmp_dist)

    ### set up SCALAR_FLUX component
    scalarFlux_comp = Component(tmpTermList,SCALAR_FLUX)

    ### solve the scalar flux component by first integrating over mu and omega to determine coefficients
    scalarFlux_comp.integrate_angles() 

    #scalarFlux_comp.print_component()
    scalarFlux_comp.write_latex_equation()

    ### generate a null component for SOURCE_3D
    tmpTermList = []
    ### SOURCE_3D
    tmpTerm = Term([Subterm(1.0,0,0,copy_mu(null_muterm),copy_omega(null_omegaterm))],SOURCE_3D)
    tmpTermList.append(tmpTerm)
    source3D_comp = Component(tmpTermList,SOURCE_3D)

    myComponentList = ComponentList([scalarFlux_comp,source3D_comp])

    ### generate a component for each axial TL moment and each radial TL moment

    ### 1D angular flux moments
    ### iterate through the Legendre and Fourier moments, creating a component for each
    cos_order = 0
    sin_order = 0
    for jleg in range(0,nmom):
        # for the jth moment, integrate over expression for 1D angular flux,
        # operated on by (\mu d/dz) and P_j(\mu), and cos(pw)
        #thisLegList = [jleg]
        for jfourier in range(0,nfourier):
            print "Generating axial TL component for l = %i, p = %i" % (jleg,jfourier)
            #if(jfourier > 0):
            #    this_cos_list = [(jfourier,1)]
            #else:
            #    this_cos_list = []

            ### iterate through Legendre and Fourier moments, creating a term for each
            tmpTermList = []

            if(jfourier == 0):
                ### SOURCE_3D
                tmpMuTerm = muTerm(cos_order,sin_order,[jleg])
                tmpOmegaTerm = omegaTerm([],[])
                tmpTerm = Term([Subterm(1.0/(4.0*math.pi),0,0,tmpMuTerm,tmpOmegaTerm)],SOURCE_3D)
                tmp_dist = term_mult(tmpTerm,infexpTerm_omegaz)
                tmpTermList.append(tmp_dist)
            #else:
            #    tmpMuTerm = muTerm(0,0,[])
            #    tmpOmegaTerm = omegaTerm(this_cos_list,[])
            #    tmpTerm = Term([Subterm(1.0/(4.0*math.pi),0,0,tmpMuTerm,tmpOmegaTerm)],SOURCE_3D)
            #    tmp_dist = term_mult(tmpTerm,infexpTerm_omegaz)
            #    tmpTermList.append(tmp_dist)
            #    tmpTermList = []

            for ileg in range(0,nmom):
                tmpMuTerm = muTerm(cos_order,sin_order,[jleg,ileg])
                #for ifourier in range(0,nfourier):
                if(jfourier == 0):
                    coeff = - 1.0/(2.0*math.pi) * (2.0*ileg + 1.0) / 2.0
                    #coeff = - (2.0*ileg + 1.0) / 2.0
                    tmpOmegaTerm = null_omegaterm
                else:
                    #coeff = - 1.0/(math.pi) * (2.0*ileg + 1.0) / 2.0
                    coeff = - 1.0/(2.0*math.pi) * (2.0*ileg + 1.0) / 2.0
                    #cos_list = [(ifourier,1)]
                    #cos_list.extend(this_cos_list) 
                    #tmpOmegaTerm = omegaTerm(cos_list,[])
                    #coeff = - (2.0*ileg + 1.0) / 2.0
                    tmpOmegaTerm = null_omegaterm

                momID = getMomID(ileg,jfourier,RADIAL_TL)
                tmpTerm = Term([Subterm(coeff,0,0,tmpMuTerm,tmpOmegaTerm)],momID)
                tmp_dist = term_mult(tmpTerm,infexpTerm_omegaz)
                tmpTermList.append(tmp_dist)

            ### set up component for this TL moment
            thisMomID = getMomID(jleg,jfourier,AXIAL_TL)
            newComponent = Component(tmpTermList,thisMomID)
            newComponent.integrate_angles()
            myComponentList.add_component(newComponent)


    ##############################################
    ### Set up equations for radial TL moments ###
    ##############################################
    ### iterate through the Legendre and Fourier moments, creating a component for each
    cos_order = 0
    sin_order = 1
    for jleg in range(0,nmom):
        # for the jth moment, integrate over expression for 1D angular flux,
        # operated on by (\Ox d/dx) and P_j(\mu), and cos(pw)
        for jfourier in range(0,nfourier):
            print "Generating radial TL component for l = %i, p = %i" % (jleg,jfourier)

            if(jfourier > 0):
                this_cos_list = [(jfourier,1)]
            else:
                this_cos_list = []
            this_cos_list.append((1,1))

            ### iterate through Legendre and Fourier moments, creating a term for each
            tmpTermList = []

            ### SOURCE_3D
            tmpMuTerm = muTerm(0,sin_order,[jleg])
            tmpOmegaTerm = omegaTerm(this_cos_list,[])
            tmpTerm = Term([Subterm(1.0/(4.0*math.pi),0.5,0,tmpMuTerm,tmpOmegaTerm)],SOURCE_3D)
            tmp_dist = term_mult(tmpTerm,infexpTerm_omegax)
            tmpTermList.append(tmp_dist)

            ### axial TL moments have an additional mu, halfLZ
            for ileg in range(0,nmom):
                tmpMuTerm = muTerm(1,sin_order,[jleg,ileg])
                for ifourier in range(0,nfourier):
                    if(ifourier == 0):
                        coeff = - 1.0/(2.0*math.pi) * (2.0*ileg + 1.0) / 2.0
                        tmpOmegaTerm = omegaTerm(this_cos_list,[])
                    else:
                        coeff = - 1.0/(math.pi) * (2.0*ileg + 1.0) / 2.0
                        cos_list = [(ifourier,1)]
                        cos_list.extend(this_cos_list) 
                        tmpOmegaTerm = omegaTerm(cos_list,[])

                    momID = getMomID(ileg,ifourier,AXIAL_TL)
                    tmpTerm = Term([Subterm(coeff,0.5,0.5,tmpMuTerm,tmpOmegaTerm)],momID)
                    tmp_dist = term_mult(tmpTerm,infexpTerm_omegax)
                    tmpTermList.append(tmp_dist)

            ### set up component for this TL moment
            thisMomID = getMomID(jleg,jfourier,RADIAL_TL)
            newComponent = Component(tmpTermList,thisMomID)
            newComponent.integrate_angles()
            myComponentList.add_component(newComponent)

    #print "Print all of the components!"
    for myComponent in myComponentList.component_list:
        myComponent.combine_terms()
        #myComponent.print_component()

    ### write the latex equations for 1D angular flux moments
    for ileg in range(0,nmom):
        for ifourier in range(0,nfourier):
            thisMomID = getMomID(ileg,ifourier,AXIAL_TL)
            thisAxTLComp = myComponentList.find_component(thisMomID)
            thisAxTLComp.write_latex_equation()

    ### write the latex equations for radial TL moments
    for ileg in range(0,nmom):
        for ifourier in range(0,nfourier):
            thisMomID = getMomID(ileg,ifourier,RADIAL_TL)
            thisRadTLComp = myComponentList.find_component(thisMomID)
            thisRadTLComp.write_latex_equation()

    ### first, substitute the radial TL moments into the 1D angular flux moments
    ###

    for ileg in range(0,nmom):
        for ifourier in range(0,nfourier):
            ### find the 1D angular flux component 
            thisMomID = getMomID(ileg,ifourier,AXIAL_TL)
            #print "Substituting into ID = %i" % thisMomID
            thisAngFluxComp = myComponentList.find_component(thisMomID)

            ### figure out what it depends on, subsitute those components
            dependence_list = [] 
            for term in thisAngFluxComp.term_list:
                if(term.TID > SOURCE_3D): 
                    dependence_list.append(term.TID)

            for tid in dependence_list:
                # don't sub scalar flux or SOURCE_3D yet
                tmpComponent = myComponentList.find_component(tid) 
                thisAngFluxComp.substitute_component(tmpComponent)

            thisAngFluxComp.combine_terms()
            #thisAngFluxComp.print_component()
            thisAngFluxComp.write_latex_equation()

    ### at this point, the radial TL moments are eliminated, and we can safely remove them

    for ileg in range(0,nmom):
        for ifourier in range(0,nfourier):
            ### find the 1D angular flux component 
            thisMomID = getMomID(ileg,ifourier,RADIAL_TL)
            myComponentList.remove_component(thisMomID)

    #myComponentList.print_dependence_matrix("step3")

    ### solve each component to eliminate the diagonal term
    ### after solving, substitute this into the rest of the terms
    ### 
    count = 4
    for ileg in range(0,nmom):
        for ifourier in range(0,nfourier):
            ### find the 1D angular flux component 
            thisMomID = getMomID(ileg,ifourier,AXIAL_TL)
            thisAngFluxComp = myComponentList.find_component(thisMomID)

            #print "Solving component for ID = %i" % thisMomID
            thisAngFluxComp.solve()
            ### substitute the solved expression into other terms
            for jleg in range(0,nmom):
                for jfourier in range(0,nfourier):
                    otherMomID = getMomID(jleg,jfourier,AXIAL_TL)
                    if(not (otherMomID == thisMomID)):
                        otherAngFluxComp = myComponentList.find_component(otherMomID)
                        otherAngFluxComp.substitute_component(thisAngFluxComp)

            #count += 1
            #suff = "step%i" % count
            #myComponentList.print_dependence_matrix(suff)

            #thisAngFluxComp.print_component()
            thisAngFluxComp.write_latex_equation()
    #count += 1
    #suff = "step%i" % count
    #myComponentList.print_dependence_matrix(suff)

    ### solve for the polar isotropic moments
    ileg = 0
    for ifourier in range(0,nfourier):
        thisMomID = getMomID(ileg,ifourier,AXIAL_TL)
        thisAngFluxComp = myComponentList.find_component(thisMomID)
        for term in thisAngFluxComp.term_list:
            if(term.TID == thisAngFluxComp.CID):
                print "Solving component for ID = %i" % thisMomID
                thisAngFluxComp.solve()
                ### substitute the solved expression into other relevant terms
                for jleg in range(0,nmom):
                    for jfourier in range(0,nfourier):
                        otherMomID = getMomID(jleg,jfourier,AXIAL_TL)
                        if(not (otherMomID == thisMomID)):
                            otherAngFluxComp = myComponentList.find_component(otherMomID)
                            otherAngFluxComp.substitute_component(thisAngFluxComp)

    ### substitute terms into the scalar flux
    scalarFluxComp = myComponentList.find_component(SCALAR_FLUX) 

    for ileg in range(0,nmom):
        for ifourier in range(0,nfourier):
            thisMomID = getMomID(ileg,ifourier,AXIAL_TL)
            thisAngFluxComp = myComponentList.find_component(thisMomID)
            scalarFluxComp.substitute_component(thisAngFluxComp)

    ## we have phi = ( ... ) F
    scalarFluxComp.print_component()
    scalarFluxComp.write_latex_equation()
    ## solve for F = ( ... ) phi
    scalarFluxComp.solve_F()
