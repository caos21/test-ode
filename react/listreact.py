import sys
from collections import OrderedDict
from collections import defaultdict
from antlr4 import *
from antlr4.InputStream import InputStream
from reactLexer import reactLexer
from reactParser import reactParser
from reactVisitor import reactVisitor

def updateReplace(d, key, value):
  """ Update dictionary with only new keys
  """
  if key not in d:
    d[key] = value
  else:
    print("[ii] Reassigning value ", d[key], " by ", value ," for key ", key)
    d[key] = value

class MyVisitor(reactVisitor):
  def __init__(self):
    self.memory = {}
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
    return [((second, first), third) for ((first, second), third) in aelem  ]


  def genElements(self):
    """ Generates the elements of R, Q and A matrices
    """
    # hack FIXME set the number of reactions
    self.nreactions = len(self.productslist)
    # Filter reactants list
    self.uniquereactantslist = self.genUniqueList(self.reactantslist)
    print("unique reactants: ", self.uniquereactantslist)
    self.relements = self.genListOfReactantTuples(self.uniquereactantslist)
    print("reactant tuples : ", self.relements)
    #
    self.uniqueproductslist = self.genUniqueList(self.productslist)
    print("unique products : ", self.uniqueproductslist)
    self.qelements = self.genListOfReactantTuples(self.uniqueproductslist)
    print("product tuples : ", self.qelements)
    self.aelements = self.genAElements()
    print("A tuples : ", self.aelements)

    pv = []
    for ireaction, reactants in enumerate(self.uniquereactantslist):
      row = list(self.rates.keys())[0]
      for ispec, spec in enumerate(self.species):
        for reactant in reactants:
          if ispec == reactant[1]:
            row += spec

      pv.append(row)
    print(pv)
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
      print("[dd] have prefactor ", factor)
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

    print("Constants  : ", visitor.constants)
    print("Diffusions : ", visitor.diffusions)
    print("Species    : ", visitor.species)
    print("Reactants  : ", visitor.reactants)
    print("Products   : ", visitor.products)

    print('reactants list :', visitor.reactantslist)
    print('products list  :', visitor.productslist)

    print()

    print('rates  :', visitor.rates)

    print()
    visitor.genElements()

