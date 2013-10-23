plast
=====

parsing AST analyizer

## Usage

Plast code has been copied into the example directory and the module name plast has been stripped from cmd.js for testing.
Go ahead and try it.

```bash
example>> plast -g -d 1 -f index.js
```
This tells plast to parse the AST of `myscript.js` and print all the functions at a functional depth of 1 in the AST.

-  `-f` is mandatory
-  `-d` defaults to `0`
-  `-g` tells plast to print function names


```bash
example>> echo "wrapfuncs getfuncs" | plast -w "plast" -f cmd.js
```
Wrap `wrapfuncs` and `getfuncs` functions with wrapper `plast` so that we get `plast.wrapfuncs`, `plast.getfuncs`.
Output is a console log of the transformed source file.


### Bonus points
```bash
example>> plast -g -d 1 -f index.js | plast -w "plast" -f cmd.js
```
**yeh!**


## Installation
'''bash
>> npm install -g plast
```
