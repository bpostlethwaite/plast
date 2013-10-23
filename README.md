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
example>> echo "wrapfuncs getfuncs | plast -w "plast" -f cmd.js
```
Wrap `func1 func2` and `func3` with `wrapper` so that we get `wrapper.func1`, `wrapper.func2` and `wrapper.func3`.
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
