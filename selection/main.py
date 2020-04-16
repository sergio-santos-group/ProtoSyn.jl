
from src.SelectionLexer import SelectionLexer
from src.SelectionParser import SelectionParser
from CustomListener import CustomSelectionListener

from antlr4 import FileStream, CommonTokenStream, ParseTreeWalker

def main(src_name):
    finput = FileStream(src_name)
    lexer = SelectionLexer(finput)
    stream = CommonTokenStream(lexer)
    parser = SelectionParser(stream)
    tree = parser.sele()

    listener = CustomSelectionListener()
    walker = ParseTreeWalker()
    walker.walk(listener, tree)


if __name__ == '__main__':
    import sys
    
    main(sys.argv[1])