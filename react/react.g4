/*
* Grammar for reactions
*
*/

grammar react;

// @header {
// }
// 
// @parser::members {
// def eval(self, left, op, right):
//   if reactParser.combineOperator == op.type:
//     return left, right
//   else:
//     return 0
// }


// A collection of entry
entries
  : entry+
  ;

// Each entry can be:
// - comment, starts with # ( # comment )
// - diffusion ( D[Ar+] = 3.5; )
// - constant ( Ar+ = 1.5e-21; ) or
// - reaction ( Ar* + Ar* -> Ar + Ar+ + 2e; )
entry
  : comment | ((constant | (diffusion | reaction)) SEMICOLON)
  ;

// A reaction is one of more reactants which produce
// products, and a rate name with constants aconstant eactivation and beta
// Ar* + Ar* -> Ar + Ar+ + 2e 2.0e25 1.0  15.60
reaction
  : reactants PRODUCE products (rate)? (aconstant)? (eactivation)? (beta scientific*)?
  ;

// diffusion coefficient
diffusion
  : LDIFFBRACKET symbol RBRACKET (EQUAL dvalue)?
  ;

// diffusion value
dvalue
  : scientific
  ;

// beta exponent Arrhenius: A eps^\beta exp(-eactivation/eps)
beta
  : scientific
  ;

// Activation energy Arrhenius: A eps^\beta exp(-eactivation/eps)
eactivation
  : scientific
  ;

// A constant Arrhenius: A eps^\beta exp(-eactivation/eps)
aconstant
  : scientific
  ;
// rate, name for a rate constant
// ki
rate
  : symbol
  ;

// products, a combination of several product
// Ar + Ar+ + 2e
products
  : product (combineOperator product)*
  ;

// product follows the same logic as reactants
// 2e
product
  : (scientific)?symbol
  ;

// reactants, reactant (+ reactant ...)
// Ar* + Ar*
reactants
  : reactant (combineOperator reactant)*
  ;

// reactant a factor times a symbol
// Ar*
reactant
  : (scientific)?symbol
  ;

// constant value
cvalue
  : scientific
  ;

// A constant may be followed by a value
// ng = 5.0e21; (or)
// Ar
constant
  : symbol (EQUAL cvalue)?
  ;

symbol
  : SYMBOL
  ;

// Operator for combination of reactants
combineOperator
  : PLUS
  ;

// scientific notation for numbers
scientific
  : (PLUS | MINUS)?SCIENTIFIC_NUMBER
  ;

comment
  : COMMENT
  ;

// comment line
COMMENT
  : '#' ~( '\r' | '\n' )*
  ;

// end of entry
SEMICOLON
  : ';'
  ;

// // a string
// STRING
//   : '"' (~ '"')* '"'
//   ;

// chemical symbol: starts with a letter and followed by a letter or numbers
SYMBOL
  : [a-zA-Z][a-zA-Z0-9_*@./#&+-]*
  ;


SCIENTIFIC_NUMBER
  : NUMBER (E SIGN? INTEGER)?
  ;

// number 10.10 or .10
fragment NUMBER
  : INTEGER (POINT (INTEGER)?)? | POINT INTEGER
  ;

fragment INTEGER
  : DIGIT+
  ;

fragment E
  : 'E' | 'e'
  ;


fragment SIGN
  : (PLUS | MINUS)
  ;

// fragment only see by lexer
fragment DIGIT
  : ('0' .. '9')
  ;

PLUS
  : '+'
  ;

MINUS
  : '-'
  ;

LDIFFBRACKET
  : 'D['
  ;


RBRACKET
  : ']'
  ;
// STAR
//   : '*'
//   ;

// reaction produces
PRODUCE
  : '->'
  ;

// assigns constant
EQUAL
  : '='
  ;

// floats
POINT
  : '.'
  ;

// ignore any whitespace
WHITESPACE
  : [ \t\r\n] -> skip
  ;
