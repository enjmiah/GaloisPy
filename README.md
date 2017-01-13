GaloisPy
===================

<a href='https://travis-ci.org/enjmiah/GaloisPy'>
  <img src='https://travis-ci.org/enjmiah/GaloisPy.svg?branch=master' />
</a>

A Python library for computations involving finite Galois fields.

Written in Python 3.  Experimental support for Python 2.

[Try it online.](https://repl.it/COuf/3)



### Jump to:

+ [Usage](#usage)
  + [Creating a field](#creating-a-field)
  + [Basic operations](#basic-operations)
  + [Encoding and decoding](#encoding-and-decoding)
  + [Step-by-step solutions](#step-by-step-solutions)
+ [Contributing](#contributing)



<br>

Usage
-------------------

### Creating a field

Use the constructor for the `GF` class to create a finite Galois field, with the size of the field as its first argument.

```python
>>> from Galois import *
>>> GF4 = GF(4)
```

The `GF` constructor also accepts an optional second argument, `verbose`, which determines whether or not you wish to see step-by-step solutions.  For more information, see the [Step-by-step solutions](#step-by-step-solutions) section.



### Basic operations

In **GaloisPy**, vectors are represented by Python lists: `[0, 1, 0, 1]`

Matrices are represented by lists of row vectors:

```
[[1, 0, 0],
 [0, 1, 0],
 [0, 0, 1]]
```

Note that **GaloisPy** uses Python primitives as field elements.  The elements *a*
and *b* in *GF(4)* are represented with the strings `'a'` and `'b'`
respectively.  Integers are represented by Python numbers.

To add or multiply two scalars in a field:
```python
>>> GF11 = GF(11)
>>> GF11.add(1, 1)
2
>>> GF11.add(5, 6)
0
>>> GF11.mult_scalar(2, 10)
9
>>> GF4 = GF(4)
>>> GF4.add(1, "b")
'a'
>>> GF4.mult_scalar("a", "b")
1
```

Many of the functions in **GaloisPy** work with vectors as well:

```python
>>> GF4.add([0, 1, "a", "b"], ["a", 1, "b", 0])
['a', 0, 1, 'b']
>>> GF11.add_inverse([2, 10, 8, 7, 0])
[9, 1, 3, 4, 11]
```

Notable functions:
+ `add_inverse(x)` returns the additive inverse of vector or scalar `x`
+ `negative(x)` is another name for `add_inverse(x)`
+ `mult_inverse(a)` returns the multiplicative inverse of scalar `a`
+ `exp_scalar(a, n)` returns a<sup>n</sup> over the field
+ `scale_vec(a, v)` returns scalar `a` times vector `v`
+ `is_lin_indep(S)` returns `True` if and only if `S` is a linearly independent set of vectors
+ `dot_vec(u, v)` returns the dot product of `u` and `v`
+ `rref(M)` returns the RREF of the matrix `M`
+ `rank(M)` returns the [rank](http://en.wikipedia.org/wiki/Rank_%28linear_algebra%29) of matrix `M`

All methods are pure functions (they do not have side effects).



### Encoding and decoding

+ `encode(G, w)` returns the codeword from encoding word `w` with generator matrix `G`
+ `is_generator_matrix(M)` returns `True` if and only if `M` is a valid generator matrix
+ `is_standard_form(M, 'g')` returns `True` if and only if `M` is a valid generator matrix in standard form
+ `is_standard_form(M, 'p')` returns `True` if and only if `M` is a valid parity-check matrix in standard form
+ `create_pc_matrix(G)` creates a parity-check matrix from generator matrix `G`
+ `is_pc_matrix(A, B)` returns `True` if and only if `A` and `B` are parity-check matrices of each other

All methods are pure functions (they do not have side effects).



### Step-by-step solutions

An important feature in **GaloisPy** is the ability to see step-by-step solutions.  Whether or not you see step-by-step solutions is determined by the `verbose` member of the created `GF` instance.

You can specify this when you create the field.  Alternately, you can change the `verbose` member at any time.

```python
>>> GF4 = GF(4, verbose=True)
>>> GF7 = GF(7)
>>> GF7.verbose = True
```

You can also pass an additional `verbose` argument to methods that support it.
This takes precedence over the value of the `verbose` member.

```python
>>> GF3 = GF(3)
>>> M = [[1, 1, 2, 1, 2],
         [1, 0, 1, 1, 0],
         [1, 2, 0, 1, 1],
         [1, 1, 2, 0, 2],
         [2, 2, 1, 2, 1]]
>>> M = GF3.rref(M, verbose=True)
# You will see the solution here...
```

The following functions feature step-by-step solutions:

+ `rref()`
+ `mult_inverse()`

The following functions will print some extra information if `verbose` is
`True`:

+ `is_pc_matrix()`

Step-by-step solutions use the 1-based indexing convention of mathematics.

Here's some examples:

```python
>>> GF4 = GF(4)
>>> a = "a"; b = "b"
>>> M_2 = [[0, 0, b, 0],
           [0, 0, 0, 0],
           [a, 0, b, 1],
           [1, 0, a, b]]
>>> M_2rref = GF4.rref(M_2, verbose=True)

| 0 0 b 0 |   Original matrix.
| 0 0 0 0 |
| a 0 b 1 |
| 1 0 a b |

| a 0 b 1 |   Exchange rows 1 and 3.
| 0 0 0 0 |
| 0 0 b 0 |
| 1 0 a b |

| a 0 b 1 |   PLAN: Pivot down from position (1, 1)
| 0 0 0 0 |
| 0 0 b 0 |
| 1 0 a b |

| a 0 b 1 |   Added b times row 1 to row 4.
| 0 0 0 0 |
| 0 0 b 0 |
| 0 0 0 0 |

| 1 0 a b |   Scale row 1 by b.
| 0 0 0 0 |
| 0 0 b 0 |
| 0 0 0 0 |

| 1 0 a b |   Exchange rows 2 and 3.
| 0 0 b 0 |
| 0 0 0 0 |
| 0 0 0 0 |

| 1 0 a b |   PLAN: Pivot down from position (2, 3)
| 0 0 b 0 |
| 0 0 0 0 |
| 0 0 0 0 |

| 1 0 a b |   Scale row 2 by a.
| 0 0 1 0 |
| 0 0 0 0 |
| 0 0 0 0 |

| 1 0 a b |   PLAN: Pivot up from position (2, 3)
| 0 0 1 0 |
| 0 0 0 0 |
| 0 0 0 0 |

| 1 0 0 b |   Added a times row 2 to row 1.
| 0 0 1 0 |
| 0 0 0 0 |
| 0 0 0 0 |
>>> GF983 = GF(983)
>>> GF983.mult_inverse(444, verbose=True)

983 = 2 * 444 + 95   ==>   95 = 983 - 2 * 444   ==>   t = 0 - 2 * 1 = -2
444 = 4 * 95 + 64   ==>   64 = 444 - 4 * 95   ==>   t = 1 - 4 * -2 = 9
95 = 1 * 64 + 31   ==>   31 = 95 - 1 * 64   ==>   t = -2 - 1 * 9 = -11
64 = 2 * 31 + 2   ==>   2 = 64 - 2 * 31   ==>   t = 9 - 2 * -11 = 31
31 = 15 * 2 + 1   ==>   1 = 31 - 15 * 2   ==>   t = -11 - 15 * 31 = -476
2 = 2 * 1 + 0   ==>   0 = 2 - 2 * 1
-476 mod 983 = 507
Multiplicative inverse is 507
507
```



<br>

Contributing
-------------------

See [Contributing section](CONTRIBUTING.md).
