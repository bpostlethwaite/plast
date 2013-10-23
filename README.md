plast
=====

parsing AST analyizer

## Usage


```bash
> plast -g -d 1 -f myscript.js
```
This tells plast to parse the AST of `myscript.js` and print all the functions at a functional depth of 1 in the AST.

`-f` is mandatory
`-d` defaults to `0`
`-g` tells plast to print function names


```bash
> echo "func1 func2 func3" | plast -w "wrapper" -f myscript.js
```
Wrap `func1 func2` and `func3` with `wrapper` so that we get `wrapper.func1`, `wrapper.func2` and `wrapper.func3`.
Output is a console log of the transformed source file.


### Bonus points
```bash
> plast -g -d 0 -f sourcescript.js | pasta -w "wrapper" -f externalscript.js
```
**yeh!**


## Installation
'''bash
>> npm install -g plast
```
