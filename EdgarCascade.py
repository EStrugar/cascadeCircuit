##########################################################################################################################################################################################

#File name:                             EdgarCascade.py
#Summary:                               Program to analyse cascade circuits and export the circuit characteristics, with both real and imaginary impedances
#
#
#Date:                                  04.05.2021
#Author:                                EES42
#

##########################################################################################################################################################################################
#python EdgarCascade.py a_Test_Circuit_1.net testing.csv
##########################################################################################################################################################################################

#Import required libraries:
import os
import sys                                                                                  ###sys reads the command line arguments
import csv                                                                                  ###importer .csv library to write to csv files more easily
import re                                                                                   ###support for regular expression; used to concatanate strings
import numpy as np                                                                          ###numpy allows the use of multidimensional arrays
from math import pi                                                                         ###math supports the use of pi and j

##########################################################################################################################################################################################

#input the command line arguments:

if (len(sys.argv) != 3):                                                                      ###initial test to check if the correct number of arguments are provided, reporting the error and exiting
    print("\n\tError, program needs two arguments to run\n")                                 
    sys.exit(1)

print('\tScript name: ', sys.argv[0])                                                        ###display file name in terminal for easy checking
print('\tThe arguments are: ', str(sys.argv),'\n')
inputfilename = sys.argv[1]                                                                  ###takes the first argument as the input -- expected to be a text file
outputfilename = sys.argv[2]   
                                                              ###takes the second argument as the output -- expected to be a .csv file                      
try:
        fin = open(inputfilename, 'rt')   
except FileNotFoundError:
        print('File <%s> not found'%(inputfilename))
        current_location=os.getcwd()
        # gets the directory where the program is executing
        print("executing program in directory: "+current_location)
        sys.exit(2)  # exits the program, returning 2 to the operating system
 
sourcefile = fin.read()                           
fout = open(outputfilename, 'wt') 
               
              
    
    
##########################################################################################################################################################################################

#Define class for components

class Components:                                                         ###class is used to instantiate seperate objects for each component
    def __init__(self, n1, n2, R, G, C, L):                               ###Wwhich gives easy access to the information
        self.n1 = n1
        self.n2 = n2
        self.R = R
        self.G = G
        self.C = C
        self.L = L

##########################################################################################################################################################################################

#define the functions required
                                                                                                    
def blockBreaker(sourcefile):                                                                       
    """
    exports three arrays containing the stripped CIRCUIT,TERMS and OUTPUTS variables"                            
    :param sourcefile : the input file used
    :ptype string
    :return componentVariables, sourceVariables, outputVariables
    :rtype list[strings]
    """                                                                                               ###breaks the input file into the three sub blocks seperated by the xtml type
    holder = re.compile('<CIRCUIT>(.*?)</CIRCUIT>', re.DOTALL).findall(sourcefile)                  ###finds all information between the describe markers, and temporarily stores
    componentVariables = " "                                                                        ###it. It is then stripped of comments and stored in a variable, allowing the 
    componentVariables = componentVariables.join(holder)                                            ###code to replace the data stored with that of the next block
    stripped = (' ')
    for line in componentVariables.splitlines():
        if not line.startswith('#'):
            stripped = stripped +(" ")+ line
    componentVariables = stripped.split()
    componentVariables.pop(0)                                                                       ###pop(0) used to remove the initial 'and' in the blocks

    holder = re.compile('<TERMS>(.*?)</TERMS>', re.DOTALL).findall(sourcefile)     
    sourceVariables = ""
    sourceVariables = sourceVariables.join(holder)
    stripped = (' ')
    for line in sourceVariables.splitlines():
        if not line.startswith('#'):
            stripped = stripped +(" ")+ line
    sourceVariables = stripped.split()
    sourceVariables.pop(0)

    holder = re.compile('<OUTPUT>(.*?)</OUTPUT.', re.DOTALL).findall(sourcefile)     
    outputVariables = ""
    outputVariables = outputVariables.join(holder)
    for line in outputVariables.splitlines():
        if not line.startswith('#'):
            stripped = stripped +(" ")+ line
    outputVariables = stripped.split()
    outputVariables.pop(0)
    outputVariables.pop(-1)
    return componentVariables, sourceVariables, outputVariables

def componentExtraction(componentVariables):
    """
    extracts a list of objects using a class for the CIRCUIT components
    :param componentVariables
    :ptype list[strings]
    :return objectList
    :rtype list
    """ 
    
    inputArray = []                                                                                 ###uses a blank array to iterate through the CIRCUITS block from the previous
    for i in range(len(componentVariables)):                                                        ###function and split the information, creating a dictionary
        iVariableValue = componentVariables[i]
        iVariableValue = iVariableValue.split('=')
        inputArray.append(iVariableValue)

    objectList = []
    length = int(len(inputArray)/3)                                                                 ###since any line of the input block can nly have three seperate details;
                                                                                                    ### n1, n2, and a component, the number of lines is therefore the total number
    S=0                                                                                             ### of keys divided by 3
    for i in range(0,length,1):
        ###The component will always take the third slot, so by indexing using a value for S it allows all the information to be handled in the same way
        ### all the data is converted into float for ease later. S is iterated by 3 to keep the relative index position

        if 'R' in inputArray[S+2][0]:                                                               
            objectList.append(Components(float(inputArray[S][1]), float(inputArray[S+1][1]), float(inputArray[S+2][1]), 1/float(inputArray[S+2][1]), 0, 0))
            S = S+3
        elif 'G' in inputArray[S+2][0]:
            objectList.append(Components(float(inputArray[S][1]), float(inputArray[S+1][1]), 1/float(inputArray[S+2][1]), float(inputArray[S+2][1]), 0, 0))   
            S = S+3
        elif 'C' in inputArray[S+2][0]:
            objectList.append(Components(float(inputArray[S][1]), float(inputArray[S+1][1]), 0, 0, float(inputArray[S+2][1]), 0))
            S = S+3
        elif 'L' in inputArray[S+2][0]:
            objectList.append(Components(float(inputArray[S][1]), float(inputArray[S+1][1]), 0, 0, 0, float(inputArray[S+2][1])))
            S = S+3
        else:
            print('\n\tNo component value provided')                                            ###if no component information is provided the used is alerted
        
    return objectList

def sourcesExtraction(sourceVariables):
    """
    function to extract the TERMS information
    :param sourceVariables
    :ptype list[strings]
    :return:sourceData, sourcesDict, blankSource : list of dictionaries, an array and a duplicate for modification
    :rtype list[strings], dict[any, float], list[strings]
    """
    sourcesArray = []                                                                                 ###use the same blank array technique to seperate the information
    for i in range(len(sourceVariables)):
        sVariableValue = sourceVariables[i]
        sVariableValue = sVariableValue.split('=')
        sourcesArray.append(sVariableValue)

    sourcesDict = {variable:(float(value)) for variable, value in sourcesArray}
    
    blankSource = ['VT', 'IN', 'RS', 'GS', 'RL', 'Fstart', 'Fend', 'Nfreq']                           ###using a blank array with placeholder variables, populate it using              
    blankSource[0] = sourcesDict.get('VT','no')                                                       ###the dictionary method .get() to retrive the information represented by that key
    blankSource[1] = sourcesDict.get('IN','no')
    blankSource[2] = sourcesDict.get('RS','no')
    blankSource[3] = sourcesDict.get('GS','no')
    blankSource[4] = sourcesDict.get('RL','no')
    blankSource[5] = float(sourcesDict.get('Fstart','no'))  
    blankSource[6] = float(sourcesDict.get('Fend','no'))
    blankSource[7] = float(sourcesDict.get('Nfreqs','no'))
    
    if blankSource[2] == 'no':                                                                      ###resistance to calculate admittance and vice versa, allowing voltage and 
        blankSource[2] = 1/blankSource[3]                                                           ###current inputs to be equated from each otehr
    else:
        blankSource[3] = 1/blankSource[2]

    if blankSource[0] == 'no':
        blankSource[0] = blankSource[1]*blankSource[2]      
    else:
        blankSource[1] = blankSource[0]/blankSource[2]

    sourceData = blankSource

    return sourceData, sourcesDict, blankSource

def frequencyCalcs(sourcesDict, sourceData):
    """
    calculate the frequencies required for the complex impedences
    :param : sourcesDict, sourcesData : outputs from sourcesExtraction
    :ptype dict[any, float], tuple[list[strings]]
    :return : testfrequencies : tuple containing required frequencies
    :rtype : tuple
    """
    if ('LFstart' and 'LFend') in sourcesDict:                                                      ###for logarithmic sweep variable name begins with an L. Then use geomspace from
        testFrequnecies = np.geomspace(int(sourceData[5]), int(sourceData[6]), int(sourceData[7]))  ### the numpy library, taking the three terms as its arguments
    elif ('Fstart' and 'Fend') in sourcesDict:
        testFrequnecies = np.linspace(int(sourceData[5]), int(sourceData[6]), int(sourceData[7]))   ###if no L found, the code takes the lin space method with the same arguments
    elif (int(sourceData[5]) == int(sourceData[6])):
        print('\n\t start and end frequency are the same. Only one frequency used to sample')
    else:
        print('\n\tFrequency data insufficiently defined\n')                                        ###in case the frequencies don't match, an error is generated and system exits
        sys.exit(2)
    return testFrequnecies

def orderNode1(unorganisedComponentsList):
    """
    order the components with respect to their n1 values
    :param : unorganisedComponentsList : generated in the main program
    :ptype : list
    :return : sortedComponents
    :rtype : list
    """
    sortedComponents = unorganisedComponentsList.copy()                                             ###python automaticaly orders in an ascending manner so only n1 needs to be altered
    P1=0                                                                                            ###by iterating through the length of the component list, compare adjacent n1 values and swap accordingly
    while(P1<len(sortedComponents)-1):
        P2 = P1 + 1
        while P2 < len(sortedComponents):
            if sortedComponents[P1]['n1'] > sortedComponents[P2]['n1']:
                sortedComponents[P1], sortedComponents[P2] = sortedComponents[P2], sortedComponents[P1]
                P2 += 1
            else:
                P1 += 1
                break

            return sortedComponents

def calculateImpedance(sortedComponents, testFrequencies):
    """
    Impedance caluclated here by applying various equations for resistance, capacitance and inductance
    :param : sortedComponents : orderNode1 function return
    :ptype : list
    :return : impedance
    :rtype : list
    """
    impedance=[]
    for f in testFrequencies:
        for i in range(len(sortedComponents)):                                                           ###differentiatie between shunt and series components by filtering n2 values
            if sortedComponents[i].get('n2') == 0:                                                       ###The impedance equations are then applied, 
                sortedComponents[i]['type'] = 'shunt'
            else:
                sortedComponents[i]['type'] = 'series'
            if sortedComponents[i]['n2'] == 0:
                if sortedComponents[i]['R'] != 0:
                    componentZ = 1/sortedComponents[i]['R']
                elif sortedComponents[i]['C'] != 0:
                    componentZ = (sortedComponents[i]['C']*2*f*pi)*1j
                elif sortedComponents[i]['L'] != 0:
                    componentZ = (1/(sortedComponents[i]['L']*2*f*pi))*1j
            else: sortedComponents[i]['n2'] != 0                                                        ###alternative equations used for series compoennts
            if sortedComponents[i]['R'] != 0:
                    componentZ = sortedComponents[i]['R']
            elif sortedComponents[i]['C'] != 0:
                    componentZ = (1/(sortedComponents[1]['C']*2*f*pi))*1j
            elif sortedComponents[i]['L'] != 0:
                    componentZ = (sortedComponents[i]['L']*2*f*pi)*1j
            impedance.append(componentZ)
        return impedance

def ABCDMul (impedance, sortedComponents):
    """
    Function to create the overall ABCD matrix of the circuit
    :param : impedance, sortedComponents
    :ptype : list, list
    :return : ABCD, A, B, C, D : broken up for easier manipulation when calculating equations
    :rtype : tuples
    """
    ABCD = np.array([[1,0], [0,1]], dtype = complex)                                    ###using an identity matrix and establish it as complex
    for i in range(len(sortedComponents)):                                              ###iterate through the sorted component list to guarantee the correct order
        if sortedComponents[i]['n2'] != 0:                                              ###populate the temporary matrix differently for series and shunt components
            holdABCD = np.array([[1, impedance[i]], [0, 1]])
        else:
            holdABCD = np.array([[1, 0], [(1/impedance[i]), 1]])
        ABCD = np.matmul(ABCD, holdABCD)                                                ###matmul used from the numpy library to easily multiply and store back in ABCD
    
    A = ABCD[0,0] 
    B = ABCD[0,1]
    C = ABCD[1,0]
    D = ABCD[1,1]

    return ABCD, A, B, C, D

##########################################################################################################################################################################################

#Define output equations using functions to improve versatility

##########################################################################################################################################################################################

def vIn(VT, Zin, RS):
    Vin = VT * (Zin/(Zin + RS))
    return Vin

# def iIn(VT, Zin):
#     Iin = (VT)/(Zin(A,B,C,D,RL)+Zin)
#     print(Iin)
#     return Iin 

def pIn(Vin, Iin):
    Pin = Vin*Iin
    return Pin

def pOut(Vout, Iout):
    Pout = Vout*Iout
    return Pout

def vOut(Vin, Vgain):
    Vout = Vin*Vgain
    return Vout

def iOut(Vout, Zout):
    Iout = Vout/Zout
    return Iout

def zIn(A,B,C,D, RL):
    Zin = ((A*RL)+B)/((C*RL)+D)
    return Zin

def zOut(A,B,C,D, RS):
    Zout = ((D*RS)+B)/((C*RS)+A)
    return Zout

def vGain(A,B,RL):
     Vgain = 1/(A+(B*(1/RL)))
     return Vgain

def iGain(C,D, RL):
    Igain = 1/((C*RL)+D)
    return Igain

##########################################################################################################################################################################################

#MAIN PROGRAM

##########################################################################################################################################################################################

componentVariables, sourceVariables, outputVariables = blockBreaker(sourcefile)                             ###order of functions called crucial: first extract list of components
objectList = componentExtraction(componentVariables)  

unorganisedComponentsList = []
for i in range(len(objectList)):
    unorganisedComponentsList.append(objectList[i].__dict__)                                                ###convert this into an array of dictionaries for easier indexing and manipulation using keys
    
sortedComponents = orderNode1(unorganisedComponentsList)
print('\n\tThe components organised chronologically: ',sortedComponents)                                    ###next sort the array by n1 value and display in the terminal using indenting and spacing

sourceData, sourcesDict, blankSource = sourcesExtraction(sourceVariables)                                   ###extract TERMS components and display the array in the terminal using indenting and spacing
print('\n\tThe source data is: ', sourceData)

testFrequencies = frequencyCalcs (sourcesDict, sourceData)                                                  ###calculate the requried testing frequencies for both logarithmic and linear sweeps
print('\n\t\The frequency range is :', testFrequencies)

VT = float(blankSource[0])    
IN = float(blankSource[1])
RS = float(blankSource[2])
GS = float(blankSource[3])                                                                                  ###TERMS data seperated into components for easier use when executing individual equations
RL = float(blankSource[4])
Fstart = int(blankSource[5]) 
Fend = int(blankSource[6])
Nfreqs = int(blankSource[7])

impedance = calculateImpedance(sortedComponents, testFrequencies)                                                            ###calcuate the impedances for the components inputted
print('\n\tThe component impedances are :',impedance)

ABCD, A, B, C, D = ABCDMul(impedance, sortedComponents)                                                     ###the circuit ABCD matrix is printed in the termianl
print('\n\tCircuit matrix :', ABCD)

#Iin = 0.00277                                                                                               ##the current value provided difficulties so a placeholder was used whilst trying to troubleshoot other equations 

Zin = zIn(A,B,C,D,RL)

Zout = zOut(A,B,C,D,RS)

Vgain = vGain(A,B,RL)

Igain = iGain(C,D,RL)

Vin = vIn(VT,Zout,RS)

Vout = vOut(Vin,Vgain)
""" 
Iin = iIn(VT,Zin)
 """
Iout = iOut(Vout,Zout)

#Pin = pIn(Vin,Iin)

Pout = pOut(Vout,Iout)

outputfinal = (Zout,Zin,Vgain)
outputfilename = open('testing.csv', 'wt')                                                      ###output to the argv[2] file described in the command line

try:
    outputfilename.write(outputfinal)

except TypeError:
    print('\n\tIncorrect data type for output')

outputfilename.close()                                                                          ###close the file, ending the program




