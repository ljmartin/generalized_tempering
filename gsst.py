from __future__ import print_function

"""
Generalized serial simulated tempering. The skeleton of this script 
comes from Peter Eastman's simulatedtempering.py, available 
through OpenMM www.github.com/openmm/openmm

A portion of the original license:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

import simtk.unit as unit
from simtk.openmm import CustomIntegrator
import math
import random
from sys import stdout
from scipy.special import logsumexp
import numpy as np

try:
    import bz2
    have_bz2 = True
except: have_bz2 = False

try:
    import gzip
    have_gzip = True
except: have_gzip = False

class GSST(object):
    """This script implements generalized serial simulated tempering.
    
    Some relevant references:
    - Estimates weights using Wang-Landau https://doi.org/10.1103/PhysRevLett.86.2050
    - Proposes state changes using metropolized independence sampling:
    https://doi.org/10.1063/1.3660669
    - generalized REST - this describes the 'generalized' formulation used here,
    but was for replica exchange, not serial tempering: https://doi.org/10.1063/1.5016222
    
    Usage goes like this:
      - set up an OpenMM simulation that has some custom force applied. Tempering
      will be applied to the parameter / collective variable described by this
      custom force. Supply the simulation to GSST, and tempering will be added as 
      a reporter. 
      - create a function yourself that takes a simulation object and a value for
      the parameter to be tempered, and then sets the tempering state of the simulation.
      Most commonly this will boil down to:
        def scaling_function(simulation, parameterValue):
          simulation.context.setParameter('parameterName', parameterValue)
      - define some levels across which to temper, for example lambda values between
      [0,1] (like alchemical perturbation), or r0 length values across which you 
      want to calculate a PMF (like umbrella sampling). 
      - set a cutoff. This is for estimating weights using the Wang-Landau algorithm. 
      1e-8 is safe for almost anything, go exponentially bigger (i.e. 1e-5) to 
      save time during equilibration, but at the cost of slightly incorrect weights,
      meaning the simulation might not maintain detailed balance. 
      - baseTemp is just the simulation temperature (i.e. 298*kelvin), which is 
      used to calculate kT.
      - weights - you may supply pre-equilibrated weights, otherwise they will be
      estimated using Wang-Landau algorithm. 
      - changeInterval - how often, in number of timesteps, to swap from dynamics 
      to Monte Carlo sampling through tempering space. Ideally this is long enough
      such that adjacently sampled potential energies are uncorrelated, i.e. 
      this should be longer than the autocorrelation time of the PE for the 
      process being tempered. 
      
    """

    def __init__(self, simulation, scaling_function, levels, cutoff, baseTemp,
                    weights=None, changeInterval=1000, reportInterval=1000, reportFile=stdout,):
        """Create a new SimulatedTempering.

        Parameters
        ----------
        simulation: Simulation
            The Simulation defining the System, Context, and Integrator to use
        scaling_function: function taking a simulation object, a parameterName (str),
            and a parameterValue (float).
        levels: list or ndarray
            the possible values for the parameter defined in `scaling_function`
        cutoff: float
            when to stop adjusting weights. When _weightUpdateFactor reduces below this,
            the _updateWeights flag turns off and weights no longer update (i.e. equilibration is over).
        baseTemp: temperature (kelvins)
            the temperature of the integrator used in the simulation. 
        weights: list
            The weight factor for each temperature.  If none, weights are selected automatically.
        changeInterval: int
            The interval (in time steps) at which to attempt transitions between levels
        reportInterval: int
            The interval (in time steps) at which to write information to the report file
        reportFile: string or file
            The file to write reporting information to, specified as a file name or file object
        """
        self.cutoff = cutoff
        self.scaling_function=scaling_function
        self.levels = levels
        self.baseTemp = baseTemp
        self.inverseBaseTemperature = 1.0/(unit.MOLAR_GAS_CONSTANT_R*self.baseTemp)
        self.simulation = simulation
        self.changeInterval = changeInterval
        self.reportInterval = reportInterval
        self.levels = levels

        # If necessary, open the file we will write reports to.

        self._openedFile = isinstance(reportFile, str)
        if self._openedFile:
            # Detect the desired compression scheme from the filename extension
            # and open all files unbuffered
            if reportFile.endswith('.gz'):
                if not have_gzip:
                    raise RuntimeError("Cannot write .gz file because Python could not import gzip library")
                self._out = gzip.GzipFile(fileobj=open(reportFile, 'wb', 0))
            elif reportFile.endswith('.bz2'):
                if not have_bz2:
                    raise RuntimeError("Cannot write .bz2 file because Python could not import bz2 library")
                self._out = bz2.BZ2File(reportFile, 'w', 0)
            else:
                self._out = open(reportFile, 'w', 1)
        else:
            self._out = reportFile

        # Initialize the weights.

        if weights is None:
            self._weights = [0.0]*len(self.levels)
            self._updateWeights = True
            self._weightUpdateFactor = 1.0
            self._histogram = [0]*len(self.levels)
            self._hasMadeTransition = False
        else:
            self._weights = weights
            self._updateWeights = False

        # Select the starting level of tempering.

        self.currentLevel = 0
        self.scaling_function(self.simulation, self.levels[self.currentLevel])
        #self.simulation.context.setParameter('scalingFactor', self.scalingFactors[self.currentTemperature])

        # Add a reporter to the simulation which will handle the updates and reports.

        class STReporter(object):
            def __init__(self, st, fg):
                self.st = st
                self.fg = fg

            def describeNextReport(self, simulation):
                st = self.st
                steps1 = st.changeInterval - simulation.currentStep%st.changeInterval
                steps2 = st.reportInterval - simulation.currentStep%st.reportInterval
                steps = min(steps1, steps2)
                isUpdateAttempt = (steps1 == steps)
                return (steps, False, False, False, False)

            def report(self, simulation, state):
                #calculate energies from each level:
                energies = list()
                for level in self.st.levels:
                    self.st.scaling_function(simulation, level)
                    pe = simulation.context.getState(getEnergy=True).getPotentialEnergy() #this is in kj/mole now
                    energies.append(pe)
                #set the tempered parameter back to what it used to be:
                self.st.scaling_function(simulation, self.st.levels[self.st.currentLevel])
                st = self.st
                if st._weightUpdateFactor<st.cutoff:
                    st._updateWeights=False
                if simulation.currentStep%st.changeInterval == 0:
                    st._attemptLevelChange(energies)
                if simulation.currentStep%st.reportInterval == 0:
                    st._writeReport(energies[self.st.currentLevel])

        simulation.reporters.append(STReporter(self, self.forceGroup))

        # Write out the header line.

        headers = ['Steps', 'weightUpdate', 'Level', 'PotentialEnergy']
        for t in self.levels:
            headers.append('%g Weight' % t)
        print('"%s"' % ('"\t"').join(headers), file=self._out)
        #print('"'+('"\t"').join(headers)+'"', file=self._out)


    def __del__(self):
        if self._openedFile:
            self._out.close()

    @property
    def weights(self):
        return [x-self._weights[0] for x in self._weights]

    def step(self, steps):
        """Advance the simulation by integrating a specified number of time steps."""
        self.simulation.step(steps)

    def _attemptStateChange(self, probability):
        #this is the metropolized independence sampling
        r = random.random()
        for j in range(len(probability)):
            if r < probability[j]:
                return j
            r -= probability[j]
        return self.currentTemperature

    def _attemptLevelChange(self, energies):
        """Attempt to move to a different temperature."""

        #this turns the PE into a 'reduced potential', then calculates the normalized probability that the
        #present configuration would have come from each of the levels. It then uses metropolized
        #independence sampling to pick a level. 
      
        logProbability = [(self._weights[i]-self.inverseBaseTemperature*energies[i]) for i in range(len(self._weights))]
        logProbability -= logsumexp(logProbability)
        probability = np.exp(logProbability)

        #peter eastmans no-numpy way:
        #maxLogProb = max(logProbability)
        #offset = maxLogProb + math.log(sum(math.exp(x-maxLogProb) for x in logProbability))
        #probability = [math.exp(x-offset) for x in logProbability]


        proposal = self._attemptStateChange(probability)
        if proposal != self.currentLevel:
            self._hasMadeTransition = True
            self.currentLevel = proposal
            self.scaling_function(self.simulation, self.levels[self.currentLevel])

        if self._updateWeights:
            self._weights[self.currentLevel] -= self._weightUpdateFactor
            self._histogram[self.currentLevel] += 1
            minCounts = min(self._histogram)
            if minCounts > 20 and minCounts >= 0.2*sum(self._histogram)/len(self._histogram):
                # Reduce the weight update factor and reset the histogram.
                self._weightUpdateFactor *= 0.5
                self._histogram = [0]*len(self.levels)
                self._weights = [x-self._weights[0] for x in self._weights]
            elif not self._hasMadeTransition and probability[self.currentLevel] > 0.99:
                # Rapidly increase the weight update factor at the start of the simulation to find
                # a reasonable starting value.
                self._weightUpdateFactor *= 2.0
                self._histogram = [0]*len(self.levels)
        return

    def _writeReport(self, nrg):
        """Write out a line to the report."""
        
        potEnergy = nrg.value_in_unit(unit.kilojoule/unit.mole)
        values = [self._weightUpdateFactor]+[self.levels[self.currentLevel]]+[potEnergy]+self.weights
        print(('%d\t' % self.simulation.currentStep) + '\t'.join('%g' % v for v in values), file=self._out)
