# Generated from Selection.g4 by ANTLR 4.7.1
from antlr4 import *
if __name__ is not None and "." in __name__:
    from .SelectionParser import SelectionParser
else:
    from SelectionParser import SelectionParser

# This class defines a complete generic visitor for a parse tree produced by SelectionParser.

class SelectionVisitor(ParseTreeVisitor):

    # Visit a parse tree produced by SelectionParser#sele.
    def visitSele(self, ctx:SelectionParser.SeleContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#ExprAndOrExpr.
    def visitExprAndOrExpr(self, ctx:SelectionParser.ExprAndOrExprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#ParenthesizedExpr.
    def visitParenthesizedExpr(self, ctx:SelectionParser.ParenthesizedExprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#NotExpr.
    def visitNotExpr(self, ctx:SelectionParser.NotExprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#SelectionExpr.
    def visitSelectionExpr(self, ctx:SelectionParser.SelectionExprContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#SelectionByDistance.
    def visitSelectionByDistance(self, ctx:SelectionParser.SelectionByDistanceContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#SelectionByID.
    def visitSelectionByID(self, ctx:SelectionParser.SelectionByIDContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#SelectionByName.
    def visitSelectionByName(self, ctx:SelectionParser.SelectionByNameContext):
        return self.visitChildren(ctx)


    # Visit a parse tree produced by SelectionParser#cutoff.
    def visitCutoff(self, ctx:SelectionParser.CutoffContext):
        return self.visitChildren(ctx)



del SelectionParser