# Generated from Selection.g4 by ANTLR 4.7.1
# encoding: utf-8
from antlr4 import *
from io import StringIO
from typing.io import TextIO
import sys

def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u608b\ua72a\u8133\ub9ed\u417c\u3be7\u7786\u5964\3\20")
        buf.write("9\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\3\2\7\2\f\n\2\f\2\16")
        buf.write("\2\17\13\2\3\2\3\2\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\5\3")
        buf.write("\33\n\3\3\3\3\3\3\3\7\3 \n\3\f\3\16\3#\13\3\3\4\3\4\3")
        buf.write("\4\3\4\3\4\5\4*\n\4\3\4\3\4\3\4\3\4\3\4\3\4\7\4\62\n\4")
        buf.write("\f\4\16\4\65\13\4\3\5\3\5\3\5\2\4\4\6\6\2\4\6\b\2\4\3")
        buf.write("\2\3\4\4\2\b\t\f\f\2:\2\r\3\2\2\2\4\32\3\2\2\2\6)\3\2")
        buf.write("\2\2\b\66\3\2\2\2\n\f\5\4\3\2\13\n\3\2\2\2\f\17\3\2\2")
        buf.write("\2\r\13\3\2\2\2\r\16\3\2\2\2\16\20\3\2\2\2\17\r\3\2\2")
        buf.write("\2\20\21\7\2\2\3\21\3\3\2\2\2\22\23\b\3\1\2\23\24\7\5")
        buf.write("\2\2\24\25\5\4\3\2\25\26\7\6\2\2\26\33\3\2\2\2\27\30\7")
        buf.write("\7\2\2\30\33\5\4\3\4\31\33\5\6\4\2\32\22\3\2\2\2\32\27")
        buf.write("\3\2\2\2\32\31\3\2\2\2\33!\3\2\2\2\34\35\f\6\2\2\35\36")
        buf.write("\t\2\2\2\36 \5\4\3\7\37\34\3\2\2\2 #\3\2\2\2!\37\3\2\2")
        buf.write("\2!\"\3\2\2\2\"\5\3\2\2\2#!\3\2\2\2$%\b\4\1\2%&\t\3\2")
        buf.write("\2&*\7\r\2\2\'(\t\3\2\2(*\7\16\2\2)$\3\2\2\2)\'\3\2\2")
        buf.write("\2*\63\3\2\2\2+,\f\3\2\2,-\7\n\2\2-.\5\b\5\2./\7\13\2")
        buf.write("\2/\60\5\6\4\4\60\62\3\2\2\2\61+\3\2\2\2\62\65\3\2\2\2")
        buf.write("\63\61\3\2\2\2\63\64\3\2\2\2\64\7\3\2\2\2\65\63\3\2\2")
        buf.write("\2\66\67\7\17\2\2\67\t\3\2\2\2\7\r\32!)\63")
        return buf.getvalue()


class SelectionParser ( Parser ):

    grammarFileName = "Selection.g4"

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]

    sharedContextCache = PredictionContextCache()

    literalNames = [ "<INVALID>", "'and'", "'or'", "'('", "')'", "'not'", 
                     "'residue'", "'atom'", "'within'", "'of'", "'chain'" ]

    symbolicNames = [ "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                      "<INVALID>", "<INVALID>", "<INVALID>", "<INVALID>", 
                      "<INVALID>", "<INVALID>", "CHAIN", "IDENTIFIER", "INDEX", 
                      "NUMBER", "WS" ]

    RULE_sele = 0
    RULE_expr = 1
    RULE_selection = 2
    RULE_cutoff = 3

    ruleNames =  [ "sele", "expr", "selection", "cutoff" ]

    EOF = Token.EOF
    T__0=1
    T__1=2
    T__2=3
    T__3=4
    T__4=5
    T__5=6
    T__6=7
    T__7=8
    T__8=9
    CHAIN=10
    IDENTIFIER=11
    INDEX=12
    NUMBER=13
    WS=14

    def __init__(self, input:TokenStream, output:TextIO = sys.stdout):
        super().__init__(input, output)
        self.checkVersion("4.7.1")
        self._interp = ParserATNSimulator(self, self.atn, self.decisionsToDFA, self.sharedContextCache)
        self._predicates = None



    class SeleContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def EOF(self):
            return self.getToken(SelectionParser.EOF, 0)

        def expr(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(SelectionParser.ExprContext)
            else:
                return self.getTypedRuleContext(SelectionParser.ExprContext,i)


        def getRuleIndex(self):
            return SelectionParser.RULE_sele

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSele" ):
                listener.enterSele(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSele" ):
                listener.exitSele(self)




    def sele(self):

        localctx = SelectionParser.SeleContext(self, self._ctx, self.state)
        self.enterRule(localctx, 0, self.RULE_sele)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 11
            self._errHandler.sync(self)
            _la = self._input.LA(1)
            while (((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << SelectionParser.T__2) | (1 << SelectionParser.T__4) | (1 << SelectionParser.T__5) | (1 << SelectionParser.T__6) | (1 << SelectionParser.CHAIN))) != 0):
                self.state = 8
                self.expr(0)
                self.state = 13
                self._errHandler.sync(self)
                _la = self._input.LA(1)

            self.state = 14
            self.match(SelectionParser.EOF)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx

    class ExprContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser


        def getRuleIndex(self):
            return SelectionParser.RULE_expr

     
        def copyFrom(self, ctx:ParserRuleContext):
            super().copyFrom(ctx)


    class ExprAndOrExprContext(ExprContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SelectionParser.ExprContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def expr(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(SelectionParser.ExprContext)
            else:
                return self.getTypedRuleContext(SelectionParser.ExprContext,i)


        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterExprAndOrExpr" ):
                listener.enterExprAndOrExpr(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitExprAndOrExpr" ):
                listener.exitExprAndOrExpr(self)


    class ParenthesizedExprContext(ExprContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SelectionParser.ExprContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def expr(self):
            return self.getTypedRuleContext(SelectionParser.ExprContext,0)


        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterParenthesizedExpr" ):
                listener.enterParenthesizedExpr(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitParenthesizedExpr" ):
                listener.exitParenthesizedExpr(self)


    class NotExprContext(ExprContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SelectionParser.ExprContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def expr(self):
            return self.getTypedRuleContext(SelectionParser.ExprContext,0)


        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterNotExpr" ):
                listener.enterNotExpr(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitNotExpr" ):
                listener.exitNotExpr(self)


    class SelectionExprContext(ExprContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SelectionParser.ExprContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def selection(self):
            return self.getTypedRuleContext(SelectionParser.SelectionContext,0)


        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSelectionExpr" ):
                listener.enterSelectionExpr(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSelectionExpr" ):
                listener.exitSelectionExpr(self)



    def expr(self, _p:int=0):
        _parentctx = self._ctx
        _parentState = self.state
        localctx = SelectionParser.ExprContext(self, self._ctx, _parentState)
        _prevctx = localctx
        _startState = 2
        self.enterRecursionRule(localctx, 2, self.RULE_expr, _p)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 24
            self._errHandler.sync(self)
            token = self._input.LA(1)
            if token in [SelectionParser.T__2]:
                localctx = SelectionParser.ParenthesizedExprContext(self, localctx)
                self._ctx = localctx
                _prevctx = localctx

                self.state = 17
                self.match(SelectionParser.T__2)
                self.state = 18
                self.expr(0)
                self.state = 19
                self.match(SelectionParser.T__3)
                pass
            elif token in [SelectionParser.T__4]:
                localctx = SelectionParser.NotExprContext(self, localctx)
                self._ctx = localctx
                _prevctx = localctx
                self.state = 21
                self.match(SelectionParser.T__4)
                self.state = 22
                self.expr(2)
                pass
            elif token in [SelectionParser.T__5, SelectionParser.T__6, SelectionParser.CHAIN]:
                localctx = SelectionParser.SelectionExprContext(self, localctx)
                self._ctx = localctx
                _prevctx = localctx
                self.state = 23
                self.selection(0)
                pass
            else:
                raise NoViableAltException(self)

            self._ctx.stop = self._input.LT(-1)
            self.state = 31
            self._errHandler.sync(self)
            _alt = self._interp.adaptivePredict(self._input,2,self._ctx)
            while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                if _alt==1:
                    if self._parseListeners is not None:
                        self.triggerExitRuleEvent()
                    _prevctx = localctx
                    localctx = SelectionParser.ExprAndOrExprContext(self, SelectionParser.ExprContext(self, _parentctx, _parentState))
                    self.pushNewRecursionContext(localctx, _startState, self.RULE_expr)
                    self.state = 26
                    if not self.precpred(self._ctx, 4):
                        from antlr4.error.Errors import FailedPredicateException
                        raise FailedPredicateException(self, "self.precpred(self._ctx, 4)")
                    self.state = 27
                    _la = self._input.LA(1)
                    if not(_la==SelectionParser.T__0 or _la==SelectionParser.T__1):
                        self._errHandler.recoverInline(self)
                    else:
                        self._errHandler.reportMatch(self)
                        self.consume()
                    self.state = 28
                    self.expr(5) 
                self.state = 33
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,2,self._ctx)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.unrollRecursionContexts(_parentctx)
        return localctx

    class SelectionContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser


        def getRuleIndex(self):
            return SelectionParser.RULE_selection

     
        def copyFrom(self, ctx:ParserRuleContext):
            super().copyFrom(ctx)


    class SelectionByDistanceContext(SelectionContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SelectionParser.SelectionContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def selection(self, i:int=None):
            if i is None:
                return self.getTypedRuleContexts(SelectionParser.SelectionContext)
            else:
                return self.getTypedRuleContext(SelectionParser.SelectionContext,i)

        def cutoff(self):
            return self.getTypedRuleContext(SelectionParser.CutoffContext,0)


        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSelectionByDistance" ):
                listener.enterSelectionByDistance(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSelectionByDistance" ):
                listener.exitSelectionByDistance(self)


    class SelectionByIDContext(SelectionContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SelectionParser.SelectionContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def INDEX(self):
            return self.getToken(SelectionParser.INDEX, 0)

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSelectionByID" ):
                listener.enterSelectionByID(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSelectionByID" ):
                listener.exitSelectionByID(self)


    class SelectionByNameContext(SelectionContext):

        def __init__(self, parser, ctx:ParserRuleContext): # actually a SelectionParser.SelectionContext
            super().__init__(parser)
            self.copyFrom(ctx)

        def IDENTIFIER(self):
            return self.getToken(SelectionParser.IDENTIFIER, 0)

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterSelectionByName" ):
                listener.enterSelectionByName(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitSelectionByName" ):
                listener.exitSelectionByName(self)



    def selection(self, _p:int=0):
        _parentctx = self._ctx
        _parentState = self.state
        localctx = SelectionParser.SelectionContext(self, self._ctx, _parentState)
        _prevctx = localctx
        _startState = 4
        self.enterRecursionRule(localctx, 4, self.RULE_selection, _p)
        self._la = 0 # Token type
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 39
            self._errHandler.sync(self)
            la_ = self._interp.adaptivePredict(self._input,3,self._ctx)
            if la_ == 1:
                localctx = SelectionParser.SelectionByNameContext(self, localctx)
                self._ctx = localctx
                _prevctx = localctx

                self.state = 35
                _la = self._input.LA(1)
                if not((((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << SelectionParser.T__5) | (1 << SelectionParser.T__6) | (1 << SelectionParser.CHAIN))) != 0)):
                    self._errHandler.recoverInline(self)
                else:
                    self._errHandler.reportMatch(self)
                    self.consume()
                self.state = 36
                self.match(SelectionParser.IDENTIFIER)
                pass

            elif la_ == 2:
                localctx = SelectionParser.SelectionByIDContext(self, localctx)
                self._ctx = localctx
                _prevctx = localctx
                self.state = 37
                _la = self._input.LA(1)
                if not((((_la) & ~0x3f) == 0 and ((1 << _la) & ((1 << SelectionParser.T__5) | (1 << SelectionParser.T__6) | (1 << SelectionParser.CHAIN))) != 0)):
                    self._errHandler.recoverInline(self)
                else:
                    self._errHandler.reportMatch(self)
                    self.consume()
                self.state = 38
                self.match(SelectionParser.INDEX)
                pass


            self._ctx.stop = self._input.LT(-1)
            self.state = 49
            self._errHandler.sync(self)
            _alt = self._interp.adaptivePredict(self._input,4,self._ctx)
            while _alt!=2 and _alt!=ATN.INVALID_ALT_NUMBER:
                if _alt==1:
                    if self._parseListeners is not None:
                        self.triggerExitRuleEvent()
                    _prevctx = localctx
                    localctx = SelectionParser.SelectionByDistanceContext(self, SelectionParser.SelectionContext(self, _parentctx, _parentState))
                    self.pushNewRecursionContext(localctx, _startState, self.RULE_selection)
                    self.state = 41
                    if not self.precpred(self._ctx, 1):
                        from antlr4.error.Errors import FailedPredicateException
                        raise FailedPredicateException(self, "self.precpred(self._ctx, 1)")
                    self.state = 42
                    self.match(SelectionParser.T__7)
                    self.state = 43
                    self.cutoff()
                    self.state = 44
                    self.match(SelectionParser.T__8)
                    self.state = 45
                    self.selection(2) 
                self.state = 51
                self._errHandler.sync(self)
                _alt = self._interp.adaptivePredict(self._input,4,self._ctx)

        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.unrollRecursionContexts(_parentctx)
        return localctx

    class CutoffContext(ParserRuleContext):

        def __init__(self, parser, parent:ParserRuleContext=None, invokingState:int=-1):
            super().__init__(parent, invokingState)
            self.parser = parser

        def NUMBER(self):
            return self.getToken(SelectionParser.NUMBER, 0)

        def getRuleIndex(self):
            return SelectionParser.RULE_cutoff

        def enterRule(self, listener:ParseTreeListener):
            if hasattr( listener, "enterCutoff" ):
                listener.enterCutoff(self)

        def exitRule(self, listener:ParseTreeListener):
            if hasattr( listener, "exitCutoff" ):
                listener.exitCutoff(self)




    def cutoff(self):

        localctx = SelectionParser.CutoffContext(self, self._ctx, self.state)
        self.enterRule(localctx, 6, self.RULE_cutoff)
        try:
            self.enterOuterAlt(localctx, 1)
            self.state = 52
            self.match(SelectionParser.NUMBER)
        except RecognitionException as re:
            localctx.exception = re
            self._errHandler.reportError(self, re)
            self._errHandler.recover(self, re)
        finally:
            self.exitRule()
        return localctx



    def sempred(self, localctx:RuleContext, ruleIndex:int, predIndex:int):
        if self._predicates == None:
            self._predicates = dict()
        self._predicates[1] = self.expr_sempred
        self._predicates[2] = self.selection_sempred
        pred = self._predicates.get(ruleIndex, None)
        if pred is None:
            raise Exception("No predicate with index:" + str(ruleIndex))
        else:
            return pred(localctx, predIndex)

    def expr_sempred(self, localctx:ExprContext, predIndex:int):
            if predIndex == 0:
                return self.precpred(self._ctx, 4)
         

    def selection_sempred(self, localctx:SelectionContext, predIndex:int):
            if predIndex == 1:
                return self.precpred(self._ctx, 1)
         




