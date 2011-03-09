#lsalib.py -- Library of Local Similarity Analysis(LSA) Package

#LICENSE: BSD

#Copyright (c) 2008 Li Charles Xia
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions
#are met:
#1. Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#3. The name of the author may not be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
#IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
#OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
#IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
#NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""lsalib.py -- Library of Local Similarity Analysis(LSA) Package

    NOTE: numpy and scipy is required to use thie module
    NOTE: accepts input sequence table as delimited text file 
        with first row is the "Date" and other factor labels and
        first column is the time spot labels
"""
#import public resources
from numpy import *
import numpy
from scipy import *
import scipy
from scipy.stats import norm
from scipy.stats import stats
from scipy.stats.distributions import t
from numpy.random import shuffle
import pdb
import csv

#import lower level resource
from lsa import compcore

#################################
# normalTransform performs rank normal score transform (Ker-chau Li, PNAS 2002)
# rawMatrix is an original data matrix, with each column as factor and each row as a time point
# rankScoreMatrix(i,:)=inverse normal cumulative distribution of rank(rawMatrix(i,:))/(rowNum+1)
#################################

def normalTransform( rawMatrix ):
    """normal transformation of the raw data

        rawMatrix - matrix of the raw data
    """
    rankScoreMatrix = -ones(rawMatrix.shape, dtype=float64)
    for i in range(rawMatrix.shape[1]):
        rankScoreMatrix[:,i] = norm.ppf( stats.rankdata(rawMatrix[:,i])/(rawMatrix.shape[0]+1) )
    return rankScoreMatrix

def rowNormalTransform( rawMatrix ):
    rankScoreMatrix = -ones(rawMatrix.shape, dtype=float64)
    for i in range(rawMatrix.shape[0]):
      rankScoreMatrix[i,:] = norm.ppf( stats.rankdata(rawMatrix[i,:])/(rawMatrix.shape[1]+1) )
    return rankScoreMatrix

###############################
# temporalChangeSeq 
#
###############################

def temporalChangeSeq( dataX ):
    """calculate the temporal change between pairs of sequences

        dataX - matrix of sequence data, raw data
    """
    data1 = zeros((dataX.shape[0], dataX.shape[1]-1), dtype=float64)
    data2 = zeros(data1.shape, dtype=float64)
    data1 = dataX[:,1:dataX.shape[1]] - dataX[:,0:dataX.shape[1]-1]
    data2 = dataX[:,1:dataX.shape[1]] + dataX[:,0:dataX.shape[1]-1]
    data2 = where(data2==0, 108, data2) # not ought to be 108, any value except 0 will work
    return data1/data2
            
###############################
# colSum20ne
#
############################### 

def colSum20ne( dataX ):
    """normalize data by column sum
        
        dataX - matrix of sequence data, either raw data or transformed
    """
    dataP = zeros(dataX.shape)
    csum1 = dataX.sum(0) 
    dataP = dataX / csum1
    return dataP

###############################
# normalization
#
############################### 

def dataNormalize( dataX ): # assume missing value has already been dealt
    """normal(Z) transformation of raw data
        
        dataX - matrix of sequence data, raw data
    """
    if rank(dataX) == 1:
        dataX = (dataX - mean(dataX))/std(dataX)
    if rank(dataX) == 2:
        for i in range(dataX.shape[1]):
            dataX[:,i] = (dataX[:,i] - (dataX[:,i]).mean())/(dataX[:,i]).std()
    return dataX

def rowDataNormalize( dataX ):
    if rank(dataX) == 1:
        dataX = (dataX - mean(dataX))/std(dataX)
    if rank(dataX) == 2:
        for i in range(dataX.shape[0]):
            dataX[i,:] = (dataX[i,:] - (dataX[i,:]).mean())/(dataX[i,:]).std()
    return dataX
    
###############################
# localSimilarity
# return a list = [ L-score, start-X, start-Y, length, Ahead/Lag ]
############################### 

def localSimilarity( lsTS1, lsTS2, maxDelay=3, scale=False ):
    """do local simularity alignment and return a list in following format:
        [ Local Similarity(LS) Score, Seq X's Align Start Position, 
            Seq Y's Align Start Position, Length of Aligned Potion, 
            Positive/Negtive Correlated ]

    lsTS1 - sequence data of Seq X
    lsTS2 - sequence data of Seq Y
    maxDelay - maximum time unit of delayed response allowed
    scale - scale the data before alignment, default: False
    """
    #if scale == True:
    #   lsTS1 = (lsTS1-mean(lsTS1))/std(lsTS1)
    #   lsTS2 = (lsTS2-mean(lsTS2))/std(lsTS2)
    posScoreMatrix = zeros( (lsTS1.shape[0]+1,lsTS2.shape[0]+1), dtype=float64 )
    negScoreMatrix = zeros( (lsTS1.shape[0]+1,lsTS2.shape[0]+1), dtype=float64 )
    scoreMax = 0.
    posOrNeg = 0 # 1 for positive, 0 for negative, -1 for unknown 
    startX = 0
    startY = 0
    for i in range(1, lsTS1.shape[0]+1):
        #for j in range(1, lsTS2.shape[0]+1):
        for j in range(maximum(1, i-maxDelay), minimum(lsTS2.shape[0]+1, i+maxDelay+1)):
            #if abs(i-j) > maxDelay:
            #   continue
            #else:
            posScoreMatrix[i,j]=max(0, posScoreMatrix[i-1,j-1]+lsTS1[i-1]*lsTS2[j-1])
            negScoreMatrix[i,j]=max(0, negScoreMatrix[i-1,j-1]-lsTS1[i-1]*lsTS2[j-1])
            if posScoreMatrix[i,j] > scoreMax:
                scoreMax = posScoreMatrix[i,j]
                startX = i
                startY = j
                posOrNeg = 1
            if negScoreMatrix[i,j] > scoreMax:
                scoreMax = negScoreMatrix[i,j]
                startX = i
                startY = j
                posOrNeg = -1
    #thresh = .00001
    length = 0
    if posOrNeg == 1:
        while posScoreMatrix[startX - length, startY - length] != 0:
            length = length + 1
    elif posOrNeg == -1:
        while negScoreMatrix[startX - length, startY - length] != 0:
            length = length + 1
    #pdb.set_trace()
    return [scoreMax, startX - length + 1, startY - length + 1, length, posOrNeg]

def newLocalSimilarity( lsTS1, lsTS2, maxDelay=3, scale=False ):
    """do local simularity alignment and return a list in following format:
        [ Local Similarity(LS) Score, Seq X's Align Start Position, 
            Seq Y's Align Start Position, Length of Aligned Potion, 
            Positive/Negtive Correlated ]
        new implemented similarity alignment using external C++ routine

    lsTS1 - sequence data of Seq X
    lsTS2 - sequence data of Seq Y
    maxDelay - maximum time unit of delayed response allowed
    scale - scale the data before alignment, default: False
    """
    lsad=compcore.LSA_Data(maxDelay,lsTS1,lsTS2)
    lsar=compcore.DP_lsa(lsad)
    len=lsar.trace.size()
    if lsar.score > 0:
        posOrNeg = 1
    elif lsar.score < 0:
        posOrNeg = -1
    else:
        posOrNeg = 0
    #for i in range(0, len):
    #  print lsar.trace[i][0], "\t", lsar.trace[i][1]
    #quit()
    return [lsar.score, lsar.trace[len-1][0], lsar.trace[len-1][1], len, posOrNeg]

###############################
# sigTest
# return a float = significance(P-value) based on a permutation test
############################### 

def sigTest( lsTS1, lsTS2, refScore=0, maxDelay=3, permu=1000, scale=False ):
    """significance test by permutation, return a P-value for the LS score.
    lsTS1 - sequence data of Seq X
    lsTS2 - sequence data of Seq Y
    refScore - score from local similarity alignment
    maxDelay - maximum time unit of delayed response allowed
    permu - number of permutation time
    scale - scale the data before alignment, default: False
    """
    rTS1=zeros_like(lsTS1)+lsTS1
    rTS2=zeros_like(lsTS2)+lsTS2
    highScoreNum = 0
    for i in range(permu):
        shuffle(rTS1) 
        shuffle(rTS2) 
        scoreRandom = localSimilarity( rTS1, rTS2, maxDelay, scale )[0]
        if scoreRandom >= refScore:
            highScoreNum = highScoreNum + 1
    pValue = highScoreNum / (float64)(permu)
    return pValue

def newSigTest( lsTS1, lsTS2, refScore=0, maxDelay=3, permu=1000, scale=False ):
    """significance test by permutation, return a P-value for the LS score.
    lsTS1 - sequence data of Seq X
    lsTS2 - sequence data of Seq Y
    refScore - score from local similarity alignment
    maxDelay - maximum time unit of delayed response allowed
    permu - number of permutation time
    scale - scale the data before alignment, default: False

    using new C++ routine 
    """
    lsad=compcore.LSA_Data(maxDelay,lsTS1,lsTS2)
    lsar=compcore.DP_lsa(lsad)
    lsa_dpr=compcore.LSA_PT(lsad, lsar, permu, compcore.LSA_test)
    return lsa_dpr.pvalue

###############################
# sigTestMaxDelay
# return a list-of-list, where i-th list = [ s1, s2, L-score, X, Y, length, P/N, P-value, corr, corPval
############################### 

def sigTestMaxDelay( dataX, maxDelay=3, permu=1000 ):   
    """calculate pair to pair LS and P-value and return a lsa table.
    each row in the table is a list in following format:
    [ Seq X's Idx, Seq Y's Idx, LS Score, X's Start Position, 
        Y's Start Position, Alignment Length, Delay in Time Unit,
        P-value, Pearson' Correlation, Correlation P-value ]

    dataX - transformed data matrix
    maxDelay - maximum time unit of delayed response allowed
    permu - number of permutation time
    """
    
    lsMatrix = []
    counter = 0
    for i in range( dataX.shape[1]-1 ):
        for j in range( i+1, dataX.shape[1] ):
            pairLs = localSimilarity( dataX[:,i], dataX[:,j], maxDelay, False )
            pValue = sigTest( dataX[:,i], dataX[:,j], pairLs[0], maxDelay, permu, False )
            corr = corrcoef( dataX[:,i], dataX[:,j] )[0,1] 
            # or use sd.correlation(X,Y), see discussion online by search "correlation coefficient numpy"
            #corpVal = .5 + sign(corr)*(.5 - t.cdf( corr*sqrt((dataX.shape[0]-1)/(float)(1-corr**2)), dataX .shape[0]-1))  
            tcdf = t.cdf( corr*sqrt((dataX.shape[0]-1)/(float64)(1.000000001-corr**2)), (dataX.shape[0]-1))                 
            #adhot approach for dealing with perfect correlation
            corpVal = .5 + sign(corr)*(.5 - tcdf )
            #corpVal = .5 + sign(corr)*(.5 - t.cdf( corr*sqrt((dataX.shape[0]-1)/(float)(1-corr**2)), \
            #  (float)(dataX.shape[0]-1)))  
            pairLs[0] = "%.4f" % (pairLs[0]/(float)(dataX.shape[0])*int(pairLs[4])) # normalize and already signed
            pValue = "%.5f" % pValue
            corr = "%.4f" % corr
            corpVal = "%.5f" % corpVal
            lsMatrix.append( [i+1, j+1] + pairLs[0:4] + [pairLs[1]-pairLs[2]] + [pValue, corr, corpVal] )
            counter = counter + 1
    return lsMatrix

def newSigTestMaxDelay( dataX, maxDelay=3, permu=1000 ):   
    """calculate pair to pair LS and P-value and return a lsa table.
    each row in the table is a list in following format:
    [ Seq X's Idx, Seq Y's Idx, LS Score, X's Start Position, 
        Y's Start Position, Alignment Length, Delay in Time Unit,
        P-value, Pearson' Correlation, Correlation P-value ]

    dataX - transformed data matrix
    maxDelay - maximum time unit of delayed response allowed
    permu - number of permutation time

    using newSigTest and newLocalSimilarity
    """
    
    lsMatrix = []
    counter = 0
    for i in range( dataX.shape[1]-1 ):
        for j in range( i+1, dataX.shape[1] ):
            pairLs = newLocalSimilarity( dataX[:,i], dataX[:,j], maxDelay, False )
            pValue = newSigTest( dataX[:,i], dataX[:,j], pairLs[0], maxDelay, permu, False )
            corr = corrcoef( dataX[:,i], dataX[:,j] )[0,1] 
            tcdf = t.cdf( corr*sqrt((dataX.shape[0]-1)/(float64)(1.000000001-corr**2)), (dataX.shape[0]-1)) 
            #adhot approach for dealing with perfect correlation
            corpVal = .5 + sign(corr)*(.5 - tcdf )  
            pairLs[0] = "%.4f" % ((pairLs[0]/(float)(dataX.shape[0]))) # normalize and add sign
            pValue = "%.5f" % pValue
            corr = "%.4f" % corr
            corpVal = "%.5f" % corpVal
            lsMatrix.append( [i+1, j+1] + pairLs[0:4] + [pairLs[1]-pairLs[2]] + [pValue, corr, corpVal] )
            counter = counter + 1
    return lsMatrix
###############################
#linearFill
###############################
def linearFillRow( csvReader, csvWriter, placeHolder="na", mode="zero", skiprows=1):
    """linearly fill the missing value by in/extrapolation, serie data in column.

    csvReader - csv file descriptor of the raw data file
    csvWriter - csv file descriptor of a file for temporary storage of data
    placeHolder - place holder in data file for missing value, default: na
    mode - in/extrapolation to be used, default: intropolation only
    skiprows - rows to skip before data
    """
    idx=0
    for row in csvReader:
        if idx < skiprows:
            csvWriter.writerow(row)
            idx = idx + 1
        else:
            for i in range(1, len(row)):
                j = i
                if row[i] == "na" or row[i] == '':
                    while j < len(row):
                        if row[j] != "na" and row[j] != '':
                            break
                        j=j+1
                if i == 1: # na from start
                    if mode == "in/ex" and j <= len(row)-2 and (row[j+1] != "na" and row[j+1] != ''): # inex mode and extropolate possible
                        for k in range(i,j):
                            row[k]="%4.f" % (float(row[j])+(float(row[j])-float(row[j+1]))*(j-k))
                    else:
                        for k in range(i,j):
                            row[k]="0"
                elif j == len(row): # na at end
                    if mode == "in/ex" and i >= 2: # inex mode and extropolate possible
                        for k in range(i,j):
                            row[k]="%.4f" % (float(row[i-1])+(float(row[i-1])-float(row[i-2]))*(k-(i-1)))
                    else:
                        for k in range(i,j):
                            row[k]="0"
                else: # intropolate or zero
                    for k in range(i,j):
                        if mode != "zero":
                            row[k]="%.4f" % (float(row[i-1])+(float(row[j])-float(row[i-1]))/(j-i+1)*(k-i+1))
                        else:
                            row[k]="0"
                i = j
            csvWriter.writerow(row)  

def linearFillColumn( csvReader, csvWriter, placeHolder="na", mode="zero", skiprows=1):
    """linearly fill the missing value by in/extrapolation, serie data in column.

    csvReader - csv file descriptor of the raw data file
    csvWriter - csv file descriptor of a file for temporary storage of data
    placeHolder - place holder in data file for missing value, default: na
    mode - in/extrapolation to be used, default: intropolation only
    skiprows - rows to skip before data
    """
    table=[]
    idx=0
    for row in csvReader:
        if idx < skiprows:
            csvWriter.writerow(row)
            idx = idx + 1
            continue
        else:
            table.append(row)
    # exception handler: table empty?
    for l in range(1, len(table[0])):
        for i in range(0,len(table)):
            j = i
            if table[i][l] == "na" or table[i][1] == '':
                while j < len(table): 
                    if table[j][l] != "na" and table[j][1] != '': 
                        break
                    j = j + 1
            if i == 0: # na from start
                if mode == "in/ex" and j <= len(table)-2 and (table[j+1][l] != "na" and table[j+1][1] != ''): # inex mode and extropolate possible
                    for k in range(i,j):
                        table[k][l]="%.4f" % (float(table[j][l])+(float(table[j][l])-float(table[j+1][l]))*(j-k))
                else:
                    for k in range(i,j):
                        table[k][l]="0"
            elif j == len(table): # na from end
                if mode == "in/ex" and i >= 2: # inex mode and extropolate possible
                    for k in range(i,j):
                        table[k][l]="%.4f" % (float(table[i-1][l])+(float(table[i-1][l])-float(table[i-2][l]))*(k-(i-1)))
                else:
                    for k in range(i,j):
                        table[k][l]="0"
            else: # intropolate
                for k in range(i,j):
                    if mode != "zero":
                        table[k][l]="%.4f" % (float(table[i-1][l])+(float(table[j][l])-float(table[i-1][l]))/(j-i+1)*(k-i+1))
                    else:
                        table[k][l]="0"
            i = j
    csvWriter.writerows(table)  
