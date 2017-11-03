# Generated from react.g4 by ANTLR 4.6
from antlr4 import *
if __name__ is not None and "." in __name__:
    from .reactParser import reactParser
else:
    from reactParser import reactParser

# This class defines a complete generic visitor for a parse tree produced by reactParser.

class reactVisitor(ParseTreeVisitor):

    # Visit a parse tree produced by reactParser#entries.
    def visitEntries(self, ctx:reactParser.EntriesContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#entry.
    def visitEntry(self, ctx:reactParser.EntryContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#reaction.
    def visitReaction(self, ctx:reactParser.ReactionContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#diffusion.
    def visitDiffusion(self, ctx:reactParser.DiffusionContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#dvalue.
    def visitDvalue(self, ctx:reactParser.DvalueContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#beta.
    def visitBeta(self, ctx:reactParser.BetaContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#eactivation.
    def visitEactivation(self, ctx:reactParser.EactivationContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#aconstant.
    def visitAconstant(self, ctx:reactParser.AconstantContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#rate.
    def visitRate(self, ctx:reactParser.RateContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#products.
    def visitProducts(self, ctx:reactParser.ProductsContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#product.
    def visitProduct(self, ctx:reactParser.ProductContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#reactants.
    def visitReactants(self, ctx:reactParser.ReactantsContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#reactant.
    def visitReactant(self, ctx:reactParser.ReactantContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#cvalue.
    def visitCvalue(self, ctx:reactParser.CvalueContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#constant.
    def visitConstant(self, ctx:reactParser.ConstantContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#symbol.
    def visitSymbol(self, ctx:reactParser.SymbolContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#combineOperator.
    def visitCombineOperator(self, ctx:reactParser.CombineOperatorContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by reactParser#scientific.
    def visitScientific(self, ctx:reactParser.ScientificContext):
        return self.visitChildren(ctx)



del reactParser