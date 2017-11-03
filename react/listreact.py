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

import numpy as np

def updateReplace(d, key, value):
  """ Update dictionary with only new keys
  """
  if key not in d:
    d[key] = value
  else:
    #print("[ii] Reassigning value ", d[key], " by ", value ," for key ", key)
    d[key] = value

class MyVisitor(reactVisitor):
  def __init__(self):
    self.constants = OrderedDict({})
    self.diffusions = OrderedDict({})
    self.species = OrderedDict({})
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

    self.variablespecies = []

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


  def printSODE(self):
    #sode = []
    for ispec, species in enumerate(self.species):
      if species not in self.constants:
        sode = str()
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
              sode += ' ' +  str(list(self.species.items())[p[0]][0]) + (('^' + str(int(p[1]))) if p[1]>1 else '')
            #print(pcomp[isp])
          #for ip, pcomp in enumerate(self.pvector):
            #print
            #print(aelem, pcomp)
        print(sode)

  def printJacobian(self):
    # iterate in species, rows
    redrow = 0
    for ispec, ispecies in enumerate(self.species):
      # check if species is constant
      if ispecies not in self.constants:
        jacode = str()
        # iterate in species, cols
        redcol = 0
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

    #print("p")
    #for p in self.pvector:
      #print(p)
    #print()

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

    #print("J")
    #for j in self.jacelements:
      #print(j)


    #print(self.jelements)

    #for j in self.jelements:
      #print(j)
      #for r in self.relements:
        #print(r)

    #for ispec, species in enumerate(self.species):
    ##sode = []
    #for ispec, species in enumerate(self.species):
      #sode = str()
      ##sode.append("d[" + species + "]/dt  = ")
      #sode += "d[" + species + "]/dt  ="
      #for ia, aelem in enumerate(self.aelements):
        ## check row
        #if aelem[0][0] == ispec:
          ##for ip, pcomp in enumerate(self.pvector):
            ##[(0, 1.0), (1, 1.0), (2, 1.0)]
          #strprefac = self.strsign(int(aelem[1]))

          #sode +=  strprefac + str(list(self.rates.items())[aelem[0][1]][0])
          ## select component of p
          #for p in self.pvector[aelem[0][1]]:
            #sode += ' ' +  str(list(self.species.items())[p[0]][0]) + (('^' + str(int(p[1]))) if p[1]>1 else '')
          ##print(pcomp[isp])
        ##for ip, pcomp in enumerate(self.pvector):
          ##print
          ##print(aelem, pcomp)
      #print(sode)
    #for ireaction, reactants in enumerate(self.uniquereactantslist):
      #for p in self.pvector:
        #for t in p:
          #if t[0] == 
            #print(t)
  #def visitEntry(self, ctx):
    #reaction = self.visit(ctx.reaction())
    #return 0

  def genSysDE(self):
    """ Generates the system of differential equations
    """
    return 0

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
  visitor = MyVisitor()
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
    visitor = MyVisitor()
    visitor.visit(tree)

    #print("Constants  : ", visitor.constants)
    #print("Diffusions : ", visitor.diffusions)
    #print("Species    : ", visitor.species)
    #print("Reactants  : ", visitor.reactants)
    #print("Products   : ", visitor.products)

    #print('reactants list :', visitor.reactantslist)
    #print('products list  :', visitor.productslist)

    #print()

    #print('rates  :', visitor.rates)

    print()
    visitor.genElements()

    print()
    visitor.printSODE()

    print()
    visitor.printJacobian()

    print()
