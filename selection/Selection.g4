grammar Selection;


sele
    : (expr)* EOF
    ;

expr
    : expr ('and'|'or') expr        #ExprAndOrExpr
    | '(' expr ')'                  #ParenthesizedExpr
    | 'not' expr                    #NotExpr
    | selection                     #SelectionExpr
    ;

selection
    : ('chain'|'residue'|'atom') IDENTIFIER     #SelectionByName
    | ('chain'|'residue'|'atom') INDEX          #SelectionByID
    | selection 'within' cutoff 'of' selection  #SelectionByDistance
    ;

cutoff
    : NUMBER
    ;

/* --------------------------------------------------------
    LEXER RULES
-------------------------------------------------------- */

// RESERVED
CHAIN: 'chain' ;




fragment INT: [0-9];

fragment CHAR: [a-zA-Z_];

IDENTIFIER
    : CHAR (CHAR | INT)*
    ;

INDEX
    : INT+
    ;

// xxx (. xxx)? (eE (+-)? xxxx)?
NUMBER
    : INT+ ('.' INT+)? (('e'|'E') ('+'|'-')? INT+)?
    ;

WS : [ \t\r\n\f]+ -> channel(HIDDEN) ;
