from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import random
import sys
import math


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-w','--workspace',dest="workspace",default="wMultiJet",type="string",
                  help="input workspace name")
    parser.add_option('-f','--fit-result',dest="fitResult", default="fitresult_extRazorPdf_data_obs",type="string",
                  help="input fit result name")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input ROOT file containing workspace containing fit result")

    
    (options,args) = parser.parse_args()
     
    fitFile = rt.TFile.Open(options.inputFitFile)
    w = fitFile.Get(options.workspace)
    fr = w.obj(options.fitResult)
    
    #get the fit parmeter names 
    parList = fr.floatParsFinal()
    paramNames = []
    for p in rootTools.RootIterator.RootIterator(parList):
        paramNames.append(p.GetName())

    print "\nINFO: %s fit result!\n"%options.fitResult
    fr.Print("v")
    
    print "\nINFO: fit parameters are", paramNames
    
    print "\nINFO: retreiving %s covariance matrix\n"%options.fitResult
    covMatrix = fr.covarianceMatrix()
    covMatrix.Print("")
    paramList = rt.RooArgList()
    for paramName in paramNames:
        paramList.add(w.var(paramName))
        
    covMatrixE = rt.TMatrixDSymEigen(covMatrix)
    eigenVal = covMatrixE.GetEigenValues()
    eigenVect = covMatrixE.GetEigenVectors()
        
    eigenVectT = eigenVect.Clone()
    eigenVectT.Transpose(eigenVect)

    print "\nINFO: diagonalizing covariance matrix\n"
    diag = eigenVectT * (covMatrix *  eigenVect)
    eigenVal.Print("")
    
    eigenVal.Sqrt()

    rotEigenVal =  eigenVal.Clone()
    rotEigenVal *=  eigenVect
        
    variation = []
    for j in range(0,len(paramNames)):
        variation.append([eigenVal[j]*eigenVect[i][j] for i in range(0,len(paramNames))])
        
    cen = [w.var(paramName).getVal() for paramName in paramNames]
    err = [w.var(paramName).getError() for paramName in paramNames]

    sign = {}
    sigma = 1.0
    for p in range(0,len(paramNames)):
        print "\nINFO: Now varying fit parameters\n"
        print "eigenvector #%02d"%p
        variationName = paramNames[p]
        sign["Up",variationName] = sigma
        sign["Down",variationName] = -sigma
        if True:
            for syst in ["Up","Down"]:
                for q in range(0,len(paramNames)):
                    paramName = paramNames[q]
                    relErr = sign[syst,variationName]*variation[p][q]/(err[q])
                    print paramName, syst, " = ", cen[q]+sign[syst,variationName]*variation[p][q], " -> ", "%.2fsigma"%(relErr)
                    w.var(paramName).setVal(cen[q]+sign[syst,variationName]*variation[p][q])
