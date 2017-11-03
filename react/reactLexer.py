# Generated from react.g4 by ANTLR 4.6
from antlr4 import *
from io import StringIO


def serializedATN():
    with StringIO() as buf:
        buf.write("\3\u0430\ud6d1\u8206\uad2d\u4417\uaef1\u8d80\uaadd\2\16")
        buf.write("k\b\1\4\2\t\2\4\3\t\3\4\4\t\4\4\5\t\5\4\6\t\6\4\7\t\7")
        buf.write("\4\b\t\b\4\t\t\t\4\n\t\n\4\13\t\13\4\f\t\f\4\r\t\r\4\16")
        buf.write("\t\16\4\17\t\17\4\20\t\20\4\21\t\21\4\22\t\22\3\2\3\2")
        buf.write("\7\2(\n\2\f\2\16\2+\13\2\3\3\3\3\3\4\3\4\7\4\61\n\4\f")
        buf.write("\4\16\4\64\13\4\3\5\3\5\3\5\5\59\n\5\3\5\3\5\5\5=\n\5")
        buf.write("\3\6\3\6\3\6\5\6B\n\6\5\6D\n\6\3\6\3\6\3\6\5\6I\n\6\3")
        buf.write("\7\6\7L\n\7\r\7\16\7M\3\b\3\b\3\t\3\t\5\tT\n\t\3\n\3\n")
        buf.write("\3\13\3\13\3\f\3\f\3\r\3\r\3\r\3\16\3\16\3\17\3\17\3\17")
        buf.write("\3\20\3\20\3\21\3\21\3\22\3\22\3\22\3\22\2\2\23\3\3\5")
        buf.write("\4\7\5\t\6\13\2\r\2\17\2\21\2\23\2\25\7\27\b\31\t\33\n")
        buf.write("\35\13\37\f!\r#\16\3\2\7\4\2\f\f\17\17\4\2C\\c|\t\2%%")
        buf.write("((,-/;B\\aac|\4\2GGgg\5\2\13\f\17\17\"\"n\2\3\3\2\2\2")
        buf.write("\2\5\3\2\2\2\2\7\3\2\2\2\2\t\3\2\2\2\2\25\3\2\2\2\2\27")
        buf.write("\3\2\2\2\2\31\3\2\2\2\2\33\3\2\2\2\2\35\3\2\2\2\2\37\3")
        buf.write("\2\2\2\2!\3\2\2\2\2#\3\2\2\2\3%\3\2\2\2\5,\3\2\2\2\7.")
        buf.write("\3\2\2\2\t\65\3\2\2\2\13H\3\2\2\2\rK\3\2\2\2\17O\3\2\2")
        buf.write("\2\21S\3\2\2\2\23U\3\2\2\2\25W\3\2\2\2\27Y\3\2\2\2\31")
        buf.write("[\3\2\2\2\33^\3\2\2\2\35`\3\2\2\2\37c\3\2\2\2!e\3\2\2")
        buf.write("\2#g\3\2\2\2%)\7%\2\2&(\n\2\2\2\'&\3\2\2\2(+\3\2\2\2)")
        buf.write("\'\3\2\2\2)*\3\2\2\2*\4\3\2\2\2+)\3\2\2\2,-\7=\2\2-\6")
        buf.write("\3\2\2\2.\62\t\3\2\2/\61\t\4\2\2\60/\3\2\2\2\61\64\3\2")
        buf.write("\2\2\62\60\3\2\2\2\62\63\3\2\2\2\63\b\3\2\2\2\64\62\3")
        buf.write("\2\2\2\65<\5\13\6\2\668\5\17\b\2\679\5\21\t\28\67\3\2")
        buf.write("\2\289\3\2\2\29:\3\2\2\2:;\5\r\7\2;=\3\2\2\2<\66\3\2\2")
        buf.write("\2<=\3\2\2\2=\n\3\2\2\2>C\5\r\7\2?A\5!\21\2@B\5\r\7\2")
        buf.write("A@\3\2\2\2AB\3\2\2\2BD\3\2\2\2C?\3\2\2\2CD\3\2\2\2DI\3")
        buf.write("\2\2\2EF\5!\21\2FG\5\r\7\2GI\3\2\2\2H>\3\2\2\2HE\3\2\2")
        buf.write("\2I\f\3\2\2\2JL\5\23\n\2KJ\3\2\2\2LM\3\2\2\2MK\3\2\2\2")
        buf.write("MN\3\2\2\2N\16\3\2\2\2OP\t\5\2\2P\20\3\2\2\2QT\5\25\13")
        buf.write("\2RT\5\27\f\2SQ\3\2\2\2SR\3\2\2\2T\22\3\2\2\2UV\4\62;")
        buf.write("\2V\24\3\2\2\2WX\7-\2\2X\26\3\2\2\2YZ\7/\2\2Z\30\3\2\2")
        buf.write("\2[\\\7F\2\2\\]\7]\2\2]\32\3\2\2\2^_\7_\2\2_\34\3\2\2")
        buf.write("\2`a\7/\2\2ab\7@\2\2b\36\3\2\2\2cd\7?\2\2d \3\2\2\2ef")
        buf.write("\7\60\2\2f\"\3\2\2\2gh\t\6\2\2hi\3\2\2\2ij\b\22\2\2j$")
        buf.write("\3\2\2\2\f\2)\628<ACHMS\3\b\2\2")
        return buf.getvalue()


class reactLexer(Lexer):

    atn = ATNDeserializer().deserialize(serializedATN())

    decisionsToDFA = [ DFA(ds, i) for i, ds in enumerate(atn.decisionToState) ]


    COMMENT = 1
    SEMICOLON = 2
    SYMBOL = 3
    SCIENTIFIC_NUMBER = 4
    PLUS = 5
    MINUS = 6
    LDIFFBRACKET = 7
    RBRACKET = 8
    PRODUCE = 9
    EQUAL = 10
    POINT = 11
    WHITESPACE = 12

    modeNames = [ "DEFAULT_MODE" ]

    literalNames = [ "<INVALID>",
            "';'", "'+'", "'-'", "'D['", "']'", "'->'", "'='", "'.'" ]

    symbolicNames = [ "<INVALID>",
            "COMMENT", "SEMICOLON", "SYMBOL", "SCIENTIFIC_NUMBER", "PLUS", 
            "MINUS", "LDIFFBRACKET", "RBRACKET", "PRODUCE", "EQUAL", "POINT", 
            "WHITESPACE" ]

    ruleNames = [ "COMMENT", "SEMICOLON", "SYMBOL", "SCIENTIFIC_NUMBER", 
                  "NUMBER", "INTEGER", "E", "SIGN", "DIGIT", "PLUS", "MINUS", 
                  "LDIFFBRACKET", "RBRACKET", "PRODUCE", "EQUAL", "POINT", 
                  "WHITESPACE" ]

    grammarFileName = "react.g4"

    def __init__(self, input=None):
        super().__init__(input)
        self.checkVersion("4.6")
        self._interp = LexerATNSimulator(self, self.atn, self.decisionsToDFA, PredictionContextCache())
        self._actions = None
        self._predicates = None


