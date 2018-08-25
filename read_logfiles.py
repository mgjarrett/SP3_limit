### Michael Jarrett ### 3D SP3 Limit
### To calculate SP3 limit of 2D/1D method with anisotropic TL
#########################################
# This is a short script to read output from SP3_limit_P3.py
# and turn it into a table.
#########################################

import numpy as np
import math
import sys
import subprocess
import matplotlib.pyplot as plt

def isFloat(element):
    try:
        float(element)
        return True
    except ValueError:
        return False

if __name__ == '__main__':

    ### load quadrature information
    #nmom = int(sys.argv[1])
    #nfourier = int(sys.argv[2])

    casename = sys.argv[1]
    #casename = "gauss8_chebyshev16"
    #casename = sys.argv[1]

    max_nmom = 6
    max_nfourier = 6
    #max_nmom = 3
    #max_nfourier = 3

    ### correct SP3 limit coefficients
    lr1 = 1.0/3.0
    lr2 = 4.0/45.0
    lr3 = 44.0/945.0

    lrlz = 8.0/45.0
    lr2lz = 44.0/315.0
    lrlz2 = 44.0/315.0

    sp3_coeffs = [1.000,lr1,lr1,lr2,lr2,lrlz,lr3,lr3,lr2lz,lrlz2]
    ncoeffs = len(sp3_coeffs)

    #coeff_locs = [2,3,5,7,9,11,14,16,18,21]

    coeff_table = np.zeros([max_nmom,ncoeffs])
    rel_table = np.zeros([max_nmom,ncoeffs])

    for nmom in range(0,max_nmom):
        #for nfourier in range(0,max_nfourier):
        nfourier = nmom
         
        logname = "./%s/P%i_F%i_limit.log" % (casename,nmom+1,nfourier+1)
        #logname = "./%s/P%i_F%i_limit.log" % (casename,3,nfourier+1)
        #print subprocess.check_output("sed","-i",":a","-e","'/^\n*$/{$d;N;};/\n$/ba'",logname)

        x = subprocess.check_output(["tail", "-2", logname])

        #terms = x.split("-")
        terms = x.split()

        coeffs = []

        #for ix in coeff_locs:
        for ix in range(1,len(terms)):
            floatbool = isFloat(terms[ix])
            if(floatbool):
                coeffs.append(-float(terms[ix]))

        coeff_table[nmom,:] = coeffs
        for ix in range(0,ncoeffs):
            rel_table[nmom,ix] = coeffs[ix]/sp3_coeffs[ix]

    #print coeff_table
    #print rel_table

    ### print table with non-cross moments if yamamoto quadrature
    if(casename.find("yamamoto") == -1):

        ### write a latex table
        textable = "./tables/%s_sp3_limit_abs.tex" % casename
        texfile = open(textable,"w")  

        texfile.write("\\begin{tabular}{|cc|ccc|} \n")
        hline = "\\hline \n"
        texfile.write(hline)
        #texfile.write(" N & F & $\opL_x \opL_z$ & $\opL_x^2 \opL_z$ & $\opL_x \opL_z^2 \\\\ \n")
        texfile.write(" N & F & A & B & C \\\\ \n")
        texfile.write(hline)
        for imom in range(0,max_nmom):
            texline = " %i & %i & %7.5f & %7.5f & %7.5f \\\\ \n" % (imom+1,imom+1,coeff_table[imom,5],coeff_table[imom,8],coeff_table[imom,9])
            #texline = " %i & %i & %7.5f & %7.5f & %7.5f \\\\ \n" % (3,imom+1,coeff_table[imom,5],coeff_table[imom,8],coeff_table[imom,9])
            texfile.write(texline)
        texfile.write(hline)
        texline = " SP3 & SP3 & %7.5f & %7.5f & %7.5f \\\\ \n" % (sp3_coeffs[5],sp3_coeffs[8],sp3_coeffs[9])
        texfile.write(texline)
        texfile.write(hline)
        texfile.write("\end{tabular}")
        texfile.close()


        ### write a latex table
        textable = "./tables/%s_sp3_limit_rel.tex" % casename
        texfile = open(textable,"w")  

        texfile.write("\\begin{tabular}{|cc|ccc|} \n")
        hline = "\\hline \n"
        texfile.write(hline)
        #texfile.write(" N & F & $\opL_x \opL_z$ & $\opL_x^2 \opL_z$ & $\opL_x \opL_z^2 \\\\ \n")
        texfile.write(" N & F & A & B & C \\\\ \n")
        texfile.write(hline)
        for imom in range(0,max_nmom):
            #texline = " %i & %i & %5.3f & %5.3f & %5.3f \\\\ \n" % (imom+1,imom+1,rel_table[imom,5],rel_table[imom,8],rel_table[imom,9])
            texline = " %i & %i & %8.6f & %8.6f & %8.6f \\\\ \n" % (imom+1,imom+1,rel_table[imom,5],rel_table[imom,8],rel_table[imom,9])
            #texline = " %i & %i & %5.3f & %5.3f & %5.3f \\\\ \n" % (3,imom+1,rel_table[imom,5],rel_table[imom,8],rel_table[imom,9])
            texfile.write(texline)
        texfile.write(hline)
        texline = " SP$_3$ & SP$_3$ & %5.3f & %5.3f & %5.3f \\\\ \n" % (1.0,1.0,1.0)
        texfile.write(texline)
        texfile.write(hline)
        texfile.write("\end{tabular}")
        texfile.close()

    else:

        ### write a latex table
        textable = "./tables/%s_sp3_limit_abs.tex" % casename
        texfile = open(textable,"w")  

        texfile.write("\\begin{tabular}{|cc|ccccccccc|} \n")
        hline = "\\hline \n"
        texfile.write(hline)
        texfile.write(" N & F & $\opL_x$ & $\opL_z$ & $\opL_x^2$ & $\opL_z^2$ & $\opL_x \opL_z$ & $\opL_x^3$ & $\opL_z^3$ & $\opL_x^2 \opL_z$ & $\opL_x \opL_z^2$ \\\\ \n")
        texfile.write(hline)
        for imom in range(0,max_nmom):
            tmpstr = []
            texstr = " %i & %i" % (imom+1,imom+1)
            #texstr = " %i & %i" % (3,imom+1)
            tmpstr.append(texstr)
            for icoeff in range(1,ncoeffs):
                texstr = " & %7.5f" % coeff_table[imom,icoeff]
                tmpstr.append(texstr)
            tmpstr.append("\\\\ \n")

            nospace = ''
            texline = nospace.join(tmpstr)
            texfile.write(texline)

        texfile.write(hline)
        texfile.write("\end{tabular}")
        texfile.close()


        ### write a latex table
        textable = "./tables/%s_sp3_limit_rel.tex" % casename
        texfile = open(textable,"w")  

        #texfile.write("\\begin{tabular}{|cc|ccc|} \n")
        texfile.write("\\begin{tabular}{|cc|ccccccccc|} \n")
        hline = "\\hline \n"
        texfile.write(hline)
        #texfile.write(" N & F & $\opL_x \opL_z$ & $\opL_x^2 \opL_z$ & $\opL_x \opL_z^2$ \\\\ \n")
        texfile.write(" N & F & $\opL_x$ & $\opL_z$ & $\opL_x^2$ & $\opL_z^2$ & $\opL_x \opL_z$ & $\opL_x^3$ & $\opL_z^3$ & $\opL_x^2 \opL_z$ & $\opL_x \opL_z^2 \\\\ \n")
        texfile.write(hline)
        for imom in range(0,max_nmom):
            tmpstr = []
            texstr = " %i & %i" % (imom+1,imom+1)
            #texstr = " %i & %i" % (3,imom+1)
            tmpstr.append(texstr)
            for icoeff in range(1,ncoeffs):
                texstr = " & %7.5f" % rel_table[imom,icoeff]
                tmpstr.append(texstr)
            tmpstr.append("  \\\\ \n")

            nospace = ''
            texline = nospace.join(tmpstr)
            texfile.write(texline)

        texfile.write(hline)
        texfile.write("\end{tabular}")
        texfile.close()


    # plot the table
    table_lookup = [5,8,9]
    nplotcoeff = len(table_lookup)
    plotdata = np.zeros([max_nmom+1,nplotcoeff+2])
    for imom in range(0,max_nmom):
        plotdata[imom+1,0] = imom+1
        for ic in range(0,nplotcoeff):
            icoeff = table_lookup[ic]
            plotdata[imom+1,ic+1] = rel_table[imom,icoeff]
        plotdata[imom,-1] = 1.0
    plotdata[max_nmom,-1] = 1.0

    plt.rc('text', usetex=True)
    font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 16}
    plt.rc('font', **font)
    fig = plt.figure(1)
    lxlz  = plt.plot(plotdata[:,0],plotdata[:,1], 'bo', label='$L_x L_z$')
    lx2lz = plt.plot(plotdata[:,0],plotdata[:,2], 'rs', label='$L_x^2 L_z$')
    lxlz2 = plt.plot(plotdata[:,0],plotdata[:,3], 'gD', label='$L_x L_z^2$')
    unity = plt.plot(plotdata[:,0],plotdata[:,4], 'k:', label='Unity')
    
    titlestr = "Convergence of 2D/1D to SP3 Limit"
    plt.title(titlestr)
    plt.axis([0,6,0,1.5])
    plt.xlabel('Legendre / Fourier Expansion Order')
    plt.ylabel('Fraction of Correct SP$_3$ Limit')
    plt.legend(['$L_x L_z$','$L_x^2 L_z$','$L_x L_z^2$'])
    #plt.legend(numpoints=1)
    #plt.legend(handles=[lxlz, lx2lz, lxlz2])


    #plt.show()
    filetype = 'png'
    if(filetype == 'eps'):
        filename = "sp3_limit_convergence_%s.eps" % (casename)
        dirname = "figures"
        savename = "%s/%s" % (dirname,filename)
        fig.savefig(savename, format='eps', dpi=600)
        plt.close(fig)
    elif(filetype == 'png'):
        filename = "sp3_limit_convergence_%s.png" % (casename)
        dirname = "figures"
        savename = "%s/%s" % (dirname,filename)
        fig.savefig(savename, format='png')
        plt.close(fig)
