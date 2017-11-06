# -*- coding: utf-8 -*-
"""
Created on Oct 2017

@author: ben
"""
# Using encoding
# -*- coding: utf-8 -*-
__author__ = "Benjamin Santos"
__copyright__ = "Copyright 2017, Benjamin Santos"
__license__ = "Apache v2.0"
__version__ = "0.1.0"
__email__ = "caos21@gmail.com"
__status__ = "Development"

import sys
from collections import OrderedDict
from collections import defaultdict
from antlr4 import *
from antlr4.InputStream import InputStream
from reactLexer import reactLexer
from reactParser import reactParser
from reactVisitor import reactVisitor
from scipy.constants import constants
from scipy.constants import physical_constants
from scipy.integrate import ode
import matplotlib.pyplot as plt
import seaborn as sns

kbev = physical_constants['Boltzmann constant in eV/K'][0]

import numpy as np

def updateReplace(d, key, value):
  """ Update dictionary with only new keys
  """
  if key not in d:
    d[key] = value
  else:
    #print("[ii] Reassigning value ", d[key], " by ", value ," for key ", key)
    d[key] = value

class ReactVisitor(reactVisitor):
  """ Reactions and SODE

      Attributes
      ----------
      constants : stores the constants for the system
      diffusions : stores the diffusion constants for the species
      species : stores all the species of the system
  """
  def __init__(self):
    self.constants = OrderedDict({})
    self.diffusions = OrderedDict({})
    self.species = OrderedDict({})
    self.variablespecies = OrderedDict({})
    self.reactants = OrderedDict({})
    self.products = OrderedDict({})
    self.reactantslist = []
    self.uniquereactantslist = []
    self.productslist = []
    self.uniqueproductslist = []
    self.rates =  OrderedDict({})

    self.nreactions = []

    self.relements = []
    self.qelements = []
    self.aelements = []
    self.pvector = []
    self.jelements = []
    self.jacelements = []

    self.ratevalues = []

    self.nsode = []
    self.jacobian = []

    self.odesystem = []
    self.odeindices = OrderedDict({})

  def uniqueTupleListSum(self, mylist):
    """  Returns a list of tuple, where each tuple is unique by its first
         element, and the value of the second element is the sum of
         repeated tuples. Ex [('a', 1), ('b', 2), ('a', 3)]
                           = [('a', 4), ('b', 2)]
    """
    ddict = defaultdict(float)
    for firsttupleval, secondtupleval in mylist:
        ddict[firsttupleval] += secondtupleval

    return list(ddict.items())

  def uniqueTripleListSub(self, qlist, rlist):
    """  Returns a list of tuples (triple), where each tuple is unique by its
         first and second element, and the value of the third element is the
         substraction of repeated tuples. [(0, 1, 2), (1, 1, 2)] - [(0, 1, 4)]
                                          = [(0, 1, -2), (1, 1, 2)]
    """
    ddict = defaultdict(float)
    for qfirsttupleval, qsecondtupleval, qthirdtupleval in qlist:
      ddict[(qfirsttupleval, qsecondtupleval)] += qthirdtupleval
    for rfirsttupleval, rsecondtupleval, rthirdtupleval in rlist:
      ddict[(rfirsttupleval, rsecondtupleval)] -= rthirdtupleval

    return list(ddict.items())

  def genUniqueList(self, reactantslist):
    """ Returns a list of unique reactants/products
    """
    uniquelist = []
    for reactants in reactantslist:
      uniques = self.uniqueTupleListSum(reactants)
      uniquelist.append(uniques)
    return uniquelist

  def genListOfReactantTuples(self, uniquelist):
    """ Generates a list of reactant tuples
        (#reaction, #specie, #quantity of specie)
    """
    listofreactanttuples = []
    for ireaction, reactants in enumerate(uniquelist):
      for reactant in reactants:
        rtuple = (ireaction,
                  list(self.species.keys()).index(reactant[0]),
                  reactant[1])
        listofreactanttuples.append(rtuple)
    return listofreactanttuples

  def genAElements(self):
    """ Generates the A matrix elements
    """
    aelem = self.uniqueTripleListSub(self.qelements, self.relements)
    # transpose and eliminate zeroes
    return [((second, first), third) for ((first, second), third) in aelem  if third != 0.0]


  def filterSpecies(self):
    """ Filter variable species from species list
        no constants
    """
    for specname, specval in self.species.items():
      if specname not in self.constants:
        self.variablespecies.update({specname: specval})
        print(specname)

  def arrheniusRate(self, energy, rateval):
    """ Update rates
    """
    a = rateval['aconstant']

    print(len(rateval))
    if len(rateval) == 1:
      return a

    if len(rateval) == 2:
      ea = rateval['eactivation']
      efactor = np.exp(-3.0*ea/(2.0*energy))
      return a * efactor

    if len(rateval) == 3:
      ea = rateval['eactivation']
      beta = rateval['beta']
      tfactor = (2.0*energy/(3.0*kbev))**beta
      efactor = np.exp(-3.0*ea/(2.0*energy))
      return a*tfactor*efactor

  def updateRates(self, energy):
    """ Update rates
    """
    for irate, (ratename, rateval) in enumerate(self.rates.items()):
      self.ratevalues.append(self.arrheniusRate(energy, rateval))

  def setDensity(self, species, density):
    """ Set density for named species
    """
    self.species[species] = density

  def printSODE(self):
    """ Prints the system of differential equations
    """
    for ispec, species in enumerate(self.species):
      if species not in self.constants:
        sode = str()
        strdiff = (str(' - D[' + species + ']' + species) if species in self.diffusions else '')
        #sode.append("d[" + species + "]/dt  = ")
        sode += "d[" + species + "]/dt  ="
        for ia, aelem in enumerate(self.aelements):
          # check row
          if aelem[0][0] == ispec:
            #for ip, pcomp in enumerate(self.pvector):
              #[(0, 1.0), (1, 1.0), (2, 1.0)]
            strprefac = self.strsign(int(aelem[1]))

            sode +=  strprefac + str(list(self.rates.items())[aelem[0][1]][0])
            # select component of p
            for p in self.pvector[aelem[0][1]]:
              sode += (' ' +  str(list(self.species.items())[p[0]][0])
                       + (('^' + str(int(p[1]))) if p[1]>1 else ''))
        print(sode+strdiff)

  def checkValues(self):
    """ Checks if constants, diffusions and rates have values
    """
    for cname, cvalue in self.constants.items():
      if cvalue is None:
        for rname, rvalue in self.reactants.items():
          if cname == rname:
            print("[ee] Value for constant ", cname, "is needed")
            return -1

    for dname, dvalue in self.diffusions.items():
      if dvalue is None:
        print("[ee] Value for diffusion ", dname, "is needed")
        return -2

    for rname, rodict in self.rates.items():
      if not rodict:
        print("[ee] Aconstant for rate ", rname, "is needed")
        return -3

      if rodict['aconstant'] == 0:
        print("[ee] Aconstant value for rate ", rname, "is ", rodict['aconstant'])
        return -4

      for sname, svalue in self.species.items():
        if sname in self.reactants:
          if svalue is None:
            print("[ee] Density ", sname, "is None")
            return -5
    # all checks passed
    return 0

  def genSODE(self):
    """ Generates the system of differential equations (SODE)

        Stores the SODE in a nsode. Each entry is a tuple with the following
        structure:
        (specie_index, monomial, monomial, ..., diffusion_index)
        monomial consists in a tuple
        (prefactor, ratevalue_index, species, species, ...)
        species is the tuple
        ('name', value)
    """

    # assigns the constant value to tuple species
    for ispec, species in enumerate(self.species):
      if species in self.constants:
        self.species[species] = self.constants[species]

    # checks if we have all values
    if self.checkValues() < 0:
      print("[ii] Nothing to do!")
      return -100

    self.nsode = []
    for ispec, species in enumerate(self.species):
      if species not in self.constants:
        monomial = []
        #strdiff = (str(' - D[' + species + ']' + species) if species in self.diffusions else '')
        #diffterm = self.diffusions[species] if species in self.diffusions else None
        diffterm = list(self.diffusions.keys()).index(species) if species in self.diffusions else None
        for ia, aelem in enumerate(self.aelements):
          # check row
          eqdiffelem = []
          if aelem[0][0] == ispec:
            #eqdiffelem.append(list(self.rates.items())[aelem[0][1]]) <- rate dict
            #eqdiffelem.append(self.ratevalues[aelem[0][1])
            eqdiffelem.append(aelem[1])#<- prefactor
            eqdiffelem.append(aelem[0][1])#<- rate index
            for p in self.pvector[aelem[0][1]]:
              eqdiffelem.append(p)
          if eqdiffelem:
            #eqdiffelem.append(diffterm)
            monomial.append(tuple(eqdiffelem))
            #monomial.append(diffterm)
        tmptuple = (list(self.species.keys()).index(species), *(tuple(monomial)), diffterm)
        #self.nsode.append(tuple(monomial))
        self.nsode.append(tmptuple)

    #print(self.nsode)
    #for n in self.nsode:
      #print(n)

  def sysODE(self, t, ndensity):
    """ Generates the system of differential equations
    """
    #self.odesystem = []
    #vars = np.zeros(len(self.nsode))

    neq = len(self.nsode)
    odesystem = np.zeros(neq)
    # iterate in equation r.h.s line
    # ns = (eq-specie monome monome monome ... diffusion)
    # tns = monome/diffusion = (rate specie specie ...)
    # specie = (nspecie density)

    for ins, nsline in enumerate(self.nsode):
      # sum of monomials
      polynomial = 0.0
      #pstr = ''
      # nsline[0] is the index of species, iterate in monomials
      diffusion = 0.0
      for monos in nsline[1:]:
        # monos[0] is a prefactor
        # monos[1] is the index of rate
        # if tns is a monomial (tuple) multiply it
        monomial = 0
        #mstr = ''
        if isinstance(monos, tuple):
          monomial = monos[0] * self.ratevalues[monos[1]]
          #mstr = str(monos[0]) + ' x '+ str(self.ratevalues[monos[1]]) + ' '
          #print(self.ratevalues[monos[0]])
          #print(self.odeindices)
          for spec in monos[2:]:
            if spec[0] in self.odeindices:
              monomial *= np.power(ndensity[self.odeindices[spec[0]]], spec[1])
              #mstr += ' x ' + str(np.power(ndensity[self.odeindices[spec[0]]], spec[1]))
            else:
              monomial *= np.power(list(self.species.items())[spec[0]][1], spec[1])
              #print('MONO')
              ##mstr += ' x ' + str(np.power(list(self.species.items())[spec[0]][1], spec[1]))
        else:# maybe is diffusion
          if isinstance(monos, float):
            diffusion = -self.ratevalues[0] * ndensity[ins]
        polynomial += monomial
        #pstr += ' + ' + mstr
      odesystem[ins] = polynomial + diffusion
      #print(pstr + ' + ' + str(diffusion))
    #res = np.copy(self.s)
    #print(vars)
    #print(odesystem)
    return odesystem

  def sysJacobian(self):
    """ Returns the evaluated Jacobian
    """
    # iterate in species, rows
    redrow = 0
    self.jacobian = []
    for ispec, ispecies in enumerate(self.species):
      # check if species is constant
      jacdiffusion = None
      if ispecies not in self.constants:
        jacode = str()
        # iterate in species, cols
        redcol = 0
        #strdiff = (str(' - D[' + ispecies + ']') if ispecies in self.diffusions else '')
        #if strdiff:
          #print('J', (redrow, redrow), '=', strdiff)
        if ispecies in self.diffusions:
          jacdiffusion = list(self.diffusions.keys()).index(ispecies)
          self.jacobian.append(((redrow, redrow), jacdiffusion))
        for jspec, jspecies in enumerate(self.species):
          # check if species is constant
          if jspecies not in self.constants:
            # Jacobian indices
            pair = (ispec, jspec)
            # iterate in Jacobian matrix elements
            for ijelem, jelem in enumerate(self.jacelements):
              # if the indices == to element indices, we have an element
              if pair == jelem[0]:
                redpair = (redrow, redcol)
                # computes the prefactor
                #strfactor = self.strsign(int(jelem[1]*jelem[2]))
                jacfactor = jelem[1]*jelem[2]
                # if jelem[3] is empty, we have a constant times the rate
                if not jelem[3]:
                  #print('J', redpair, '=', strfactor, str(list(self.rates.items())[jelem[4]][0]))
                  self.jacobian.append((redpair, (jacfactor, jelem[4], None)))
                else:
                  # in this case we have a list of species
                  strksp = ''
                  jacterms = []
                  for kspec, kspecies in enumerate(jelem[3]):
                    #print(kspecies[0])
                    #strksp += (str(list(self.species.items())[kspecies[0]][0])
                            #+ (('^' + str(int(kspecies[1]))) if kspecies[1]>1 else ''))
                    jacterms.append(kspecies)
                    #print(kspecies)
                  #print('J', redpair, '=', strfactor, str(list(self.rates.items())[jelem[4]][0]), strksp)
                  self.jacobian.append((redpair, (jacfactor, jelem[4], *tuple(jacterms))))
            redcol += 1
        redrow += 1
    for j in self.nsode:
      print(j)

    for j in self.jacobian:
      print(j)


  #def stepSODE(self):
    #""" Solves step of the system of differential equations
    #"""

    ## WARNING FIXME
    ## resolves specie : #eq specie (without constants)
    #for ins, ns in enumerate(self.nsode):
      #self.odeindices.update({ns[0]: ins})

    #neq = len(self.nsode)
    #solver = ode(self.sysODE)
    ##solver.set_integrator('vode', method='bdf')
    ##solver.set_integrator('dop853')
    #solver.set_integrator('lsoda')
    ##solver.set_integrator('dopri5')
    ##[  7.15827069e-01   9.18553452e-06   2.84163746e-01]
    ##parameters to solver
    ##solver.set_f_params()

    #indensity = [list(self.species.values())[ns[0]] for ns in self.nsode]
    #solver.set_initial_value(indensity, t0)

    #t = np.linspace(t0, tf, N)
    #print(t)
    #outdensity = np.zeros((N, neq))
    #outdensity[0] = indensity
    #print(indensity)
    #k = 1
    #while solver.successful() and solver.t < tf:
      #solver.integrate(t[k])
      #outdensity[k] = solver.y
      #print('sum = ', np.sum(outdensity[k]))
      #k += 1

    ## Plot

    #for n, ns in enumerate(self.nsode):
      #fig = plt.figure()
      #label = list(self.species.keys())[ns[0]]
      #plt.plot(t, outdensity[:, n], label = label)
      #plt.ylabel('Density')
      #plt.xlabel('t')
      #plt.grid(True)
      #plt.legend()
      #plt.show()

    #print(outdensity[-1])
    ##return res


  def solveSODE(self, tf, N = 1000, t0 = 0.0):
    """ Solves the system of differential equations
    """

    # resolves specie : #eq specie (without constants)
    for ins, ns in enumerate(self.nsode):
      self.odeindices.update({ns[0]: ins})

    neq = len(self.nsode)
    solver = ode(self.sysODE)
    #solver.set_integrator('vode', method='bdf')
    #solver.set_integrator('dop853')
    solver.set_integrator('lsoda')
    #solver.set_integrator('dopri5')
    #[  7.15827069e-01   9.18553452e-06   2.84163746e-01]
    #parameters to solver
    #solver.set_f_params()

    indensity = [list(self.species.values())[ns[0]] for ns in self.nsode]
    solver.set_initial_value(indensity, t0)

    t = np.linspace(t0, tf, N)
    print(t)
    outdensity = np.zeros((N, neq))
    outdensity[0] = indensity
    print(indensity)
    k = 1
    while solver.successful() and solver.t < tf:
      solver.integrate(t[k])
      outdensity[k] = solver.y
      print('sum = ', np.sum(outdensity[k]))
      k += 1

    # Plot

    fig = plt.figure()
    for n, ns in enumerate(self.nsode):

      label = list(self.species.keys())[ns[0]]
      plt.plot(t, outdensity[:, n], label = label)
    plt.ylabel('Density')
    #plt.yscale('log')
    plt.xlabel('t')
    plt.grid(True)
    plt.legend()
    plt.show()

    #for n, ns in enumerate(self.nsode):
      #fig = plt.figure()
      #label = list(self.species.keys())[ns[0]]
      #plt.plot(t, outdensity[:, n], label = label)
      #plt.ylabel('Density')
      #plt.xlabel('t')
      #plt.grid(True)
      #plt.legend()
      #plt.show()

    print(outdensity[-1])
    #return res


  def printJacobian(self):
    # iterate in species, rows
    redrow = 0
    for ispec, ispecies in enumerate(self.species):
      # check if species is constant
      if ispecies not in self.constants:
        jacode = str()
        # iterate in species, cols
        redcol = 0
        strdiff = (str(' - D[' + ispecies + ']') if ispecies in self.diffusions else '')
        if strdiff:
          print('J', (redrow, redrow), '=', strdiff)
        for jspec, jspecies in enumerate(self.species):
          # check if species is constant
          if jspecies not in self.constants:
            # Jacobian indices
            pair = (ispec, jspec)
            # iterate in Jacobian matrix elements
            for ijelem, jelem in enumerate(self.jacelements):
              # if the indices == to element indices, we have an element
              if pair == jelem[0]:
                redpair = (redrow, redcol)
                # computes the prefactor
                strfactor = self.strsign(int(jelem[1]*jelem[2]))
                # if jelem[3] is empty, we have a constant times the rate
                if not jelem[3]:
                  print('J', redpair, '=', strfactor, str(list(self.rates.items())[jelem[4]][0]))
                else:
                  # in this case we have a list of species
                  strksp = ''
                  for kspec, kspecies in enumerate(jelem[3]):
                    #print(kspecies[0])
                    strksp += (str(list(self.species.items())[kspecies[0]][0])
                            + (('^' + str(int(kspecies[1]))) if kspecies[1]>1 else ''))
                    #print(kspecies)
                  print('J', redpair, '=', strfactor, str(list(self.rates.items())[jelem[4]][0]), strksp)
            redcol += 1
        redrow += 1

  def genJacobian(self):
    """ Populate the list with Jacobian terms
    """
    # iterate in species, rows
    redrow = 0
    self.jacobian = []
    for ispec, ispecies in enumerate(self.species):
      # check if species is constant
      jacdiffusion = None
      if ispecies not in self.constants:
        jacode = str()
        # iterate in species, cols
        redcol = 0
        #strdiff = (str(' - D[' + ispecies + ']') if ispecies in self.diffusions else '')
        #if strdiff:
          #print('J', (redrow, redrow), '=', strdiff)
        if ispecies in self.diffusions:
          jacdiffusion = list(self.diffusions.keys()).index(ispecies)
          self.jacobian.append(((redrow, redrow), jacdiffusion))
        for jspec, jspecies in enumerate(self.species):
          # check if species is constant
          if jspecies not in self.constants:
            # Jacobian indices
            pair = (ispec, jspec)
            # iterate in Jacobian matrix elements
            for ijelem, jelem in enumerate(self.jacelements):
              # if the indices == to element indices, we have an element
              if pair == jelem[0]:
                redpair = (redrow, redcol)
                # computes the prefactor
                #strfactor = self.strsign(int(jelem[1]*jelem[2]))
                jacfactor = jelem[1]*jelem[2]
                # if jelem[3] is empty, we have a constant times the rate
                if not jelem[3]:
                  #print('J', redpair, '=', strfactor, str(list(self.rates.items())[jelem[4]][0]))
                  self.jacobian.append((redpair, (jacfactor, jelem[4], None)))
                else:
                  # in this case we have a list of species
                  strksp = ''
                  jacterms = []
                  for kspec, kspecies in enumerate(jelem[3]):
                    #print(kspecies[0])
                    #strksp += (str(list(self.species.items())[kspecies[0]][0])
                            #+ (('^' + str(int(kspecies[1]))) if kspecies[1]>1 else ''))
                    jacterms.append(kspecies)
                    #print(kspecies)
                  #print('J', redpair, '=', strfactor, str(list(self.rates.items())[jelem[4]][0]), strksp)
                  self.jacobian.append((redpair, (jacfactor, jelem[4], *tuple(jacterms))))
            redcol += 1
        redrow += 1
    for j in self.nsode:
      print(j)

    for j in self.jacobian:
      print(j)

  def strsign(self, number):
    """ Return a string with sign of number and the number if abs(number) > 1
    """
    absnumber = abs(number)
    if number == 1:
      return ' + '
    if number == -1:
      return ' - '
    if number > 1:
      return ' + '+ str(absnumber)
    if number < 1:
      return ' - '+ str(absnumber)

  def genElements(self):
    """ Generates the elements of R, Q and A matrices
    """
    # hack FIXME set the number of reactions
    self.nreactions = len(self.productslist)
    # Filter reactants list
    self.uniquereactantslist = self.genUniqueList(self.reactantslist)
    #print("unique reactants: ", self.uniquereactantslist)
    self.relements = self.genListOfReactantTuples(self.uniquereactantslist)
    #print("reactant tuples : ", self.relements)
    #
    self.uniqueproductslist = self.genUniqueList(self.productslist)
    #print("unique products : ", self.uniqueproductslist)
    self.qelements = self.genListOfReactantTuples(self.uniqueproductslist)
    #print("product tuples : ", self.qelements)
    self.aelements = self.genAElements()
    #print("A tuples : ", self.aelements)

    self.pvector = []
    for ireaction in np.arange(self.nreactions):
      pcomponents = [(ritem[1], ritem[2]) for ritem in self.relements if ritem[0] == ireaction]
      self.pvector.append(pcomponents)

    # list of tuple, float, list
    # (row, col), factor, optional [species pair]
    self.jelements = []
    # iterate in components of pvector (row number)
    for ireaction, pcomponents in enumerate(self.pvector):
      # for each component
      for pci, pcomponent in enumerate(pcomponents):
        # and for each species (column number)
        for ispec, species in enumerate(self.species):
          elements = []
          # if the component is dependent of species
          if pcomponent[0] == ispec:
            # get constants with respect to species
            elements.append((ireaction, ispec,))
            constants = [p for p in pcomponents if p[0] != ispec]
            elementpairs = []
            if constants:
              # store constants if any
              elementpairs.append(constants)
            # derivative
            diff = pcomponent[1]-1
            if diff > 0:
              # store specie and exponent
              elementpairs.append([(pcomponent[0], diff)])
            # store multiplicative factor (former exponent of species)
            elements.append((pcomponent[1]))

            # flatten elements is list of lists
            if any(isinstance(el, list) for el in elementpairs):
              elementpairs = [elem for subrow in elementpairs for elem in subrow]
            # add elementpairs to pair row column
            #if elementpairs:
            elements.append(elementpairs)
            # store elements
            self.jelements.append(elements)

    self.jacelements = []
    for iaelem, aelem in enumerate(self.aelements):
      for ijelem, jelem in enumerate(self.jelements):
        if aelem[0][1] == jelem[0][0]:
          #print(aelem[0], jelem[0], '=', aelem, jelem)
          self.jacelements.append([(aelem[0][0], jelem[0][1]), aelem[1], jelem[1], jelem[2], jelem[0][0]])

  #def visitEntry(self, ctx):
    #reaction = self.visit(ctx.reaction())
    #return 0

  def visitReaction(self, ctx):
    r = self.visit(ctx.reactants())
    p = self.visit(ctx.products())
    self.reactantslist.append(r)
    self.productslist.append(p)

    # check if rate name was provided
    name = None
    if ctx.rate():
      name = self.visit(ctx.rate())
      # append string to name if rate name exists
      while name in self.rates:
        name += 'l' + str(len(self.reactantslist)-1)
      self.rates[name] = None
    else :# gives a default dummy name for rate
      # default name is k + nreaction
      name = 'k' + str(len(self.reactantslist)-1)
      # append string to name if rate name exists
      while name in self.rates:
        name += 'k' + str(len(self.reactantslist)-1)
      self.rates[name] = None

    # check if aconstant was provided
    if ctx.aconstant():
      acvalue = self.visit(ctx.aconstant())
      self.rates[name] = {"aconstant" : acvalue}

    # check if eactivation was provided
    if ctx.eactivation():
      eavalue = self.visit(ctx.eactivation())
      self.rates[name]["eactivation"] = eavalue

    # check if beta was provided
    if ctx.beta():
      beta = self.visit(ctx.beta())
      self.rates[name]["beta"] = beta

    # check if more constants were provided
    if ctx.scientific():
      arrayct = []
      for ctxsci in ctx.scientific():
        arrayct.append(self.visit(ctxsci))
      self.rates[name]["misc"] = arrayct

    return 0

  def addSpecies(self, ctx):
    symbol = self.visit(ctx.symbol())
    factor = 1.0
    if ctx.scientific():
      factor = self.visit(ctx.scientific())
      #print("[dd] have prefactor ", factor)
    updateReplace(self.species, symbol, None)
    return symbol, factor

  def retReactants(self, ctx):
    inreactants = []
    for rctx in ctx:
      symbol, value = self.visit(rctx)
      inreactants.append((symbol, value))
    return inreactants

  def visitProducts(self, ctx):
    return self.retReactants(ctx.product())

  def visitReactants(self, ctx):
    return self.retReactants(ctx.reactant())

  #def visitCombineOperator(self, ctx):
    #return (ctx.getText())

  def visitReactant(self, ctx):
    symbol, value = self.addSpecies(ctx)
    updateReplace(self.reactants, symbol, value)
    return symbol, value

  def visitBeta(self, ctx):
    value = self.visit(ctx.scientific())
    return value

  def visitEactivation(self, ctx):
    value = self.visit(ctx.scientific())
    return value

  def visitAconstant(self, ctx):
    value = self.visit(ctx.scientific())
    return value

  def visitRate(self, ctx):
    symbol = self.visit(ctx.symbol())
    return symbol

  def visitProduct(self, ctx):
    symbol, value = self.addSpecies(ctx)
    updateReplace(self.products, symbol, value)
    return symbol, value

  def visitConstant(self, ctx):
    symbol = self.visit(ctx.symbol())
    if ctx.cvalue():
      value = self.visit(ctx.cvalue())
    else :
      value = None
      print("[ii] Warning value for constant ", symbol, " is None")
    updateReplace(self.constants, symbol, value)
    return 0

  def visitCvalue(self, ctx):
    value = self.visit(ctx.scientific())
    return value

  def visitDvalue(self, ctx):
    value = self.visit(ctx.scientific())
    return value

  def visitDiffusion(self, ctx):
    symbol = self.visit(ctx.symbol())
    if ctx.dvalue():
      value = self.visit(ctx.dvalue())
    else :
      value = None
      print("[ii] Warning diffusion constant for ", symbol, " is None")
    updateReplace(self.diffusions, symbol, value)
    return 0

  def visitSymbol(self, ctx):
    return ctx.getText()

  def visitScientific(self, ctx):
    return float(ctx.getText())

    #def visitParens(self, ctx):
        #return self.visit(ctx.expr())


def readAndPrint(str_stream):

  input_stream = InputStream(str_stream)

  lexer = reactLexer(input_stream)
  token_stream = CommonTokenStream(lexer)
  parser = reactParser(token_stream)
  #tree = parser.reaction()
  tree = parser.entries()

  #print(tree.toStringTree(recog=parser))
  visitor = ReactVisitor()
  visitor.visit(tree)

  print()
  visitor.genElements()

  print()
  visitor.printSODE()

  print()
  visitor.printJacobian()

  print()
if __name__ == '__main__':
    if len(sys.argv) > 1:
        input_stream = FileStream(sys.argv[1])
    else:
        input_stream = InputStream(sys.stdin.readline())

    lexer = reactLexer(input_stream)
    token_stream = CommonTokenStream(lexer)
    parser = reactParser(token_stream)
    #tree = parser.reaction()
    tree = parser.entries()

    #print(tree.toStringTree(recog=parser))
    react = ReactVisitor()
    react.visit(tree)

    print("Constants  : ", react.constants)
    print("Diffusions : ", react.diffusions)
    print("Species    : ", react.species)
    #print("Reactants  : ", react.reactants)
    #print("Products   : ", react.products)

    #print('reactants list :', react.reactantslist)
    #print('products list  :', react.productslist)

    #print()

    print('rates  :', react.rates)

    print()
    react.genElements()

    print()
    react.printSODE()

    print()
    #react.printJacobian()


    react.updateRates(1.36)
    print(react.ratevalues)

##
    idens = 1.0e10
    react.setDensity('e', idens)
    react.setDensity('Ar+', idens)
    react.setDensity('Ar*', idens)

    print()
    react.genSODE()

    print()
    react.solveSODE(10.0)
##
    #react.setDensity('A', 1.0)
    #react.setDensity('B', 0.0)
    #react.setDensity('C', 0.0)

    #print()
    #react.genSODE()

    #print()
    #react.solveSODE(40.0)
##
    print()
    react.genJacobian()

    #print()
    #print(react.aelements)

    #print()
    #print(react.pvector)
