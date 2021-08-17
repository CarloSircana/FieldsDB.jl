# FieldsDB.jl

## About

FieldsDB is a software package providing an interface to a database of number fields maintained by Tommy Hofmann and Carlo Sircana.
It is written in [julia](https://www.julialang.org) and is based on the computer algebra system [OSCAR](https://oscar.computeralgebra.de/). 
The development is supported by the Deutsche Forschungsgemeinschaft DFG within the Collaborative Research Center TRR 195.
  
  ## Quick start

Here is a quick example of using FieldsDB:

```julia
julia> using FieldsDB
...
julia> db = fields_database();
julia> lf = load_fields(db, degree = 2, discriminant_range = (fmpz(8), fmpz(8)))
1-element Vector{FieldsDB.DBField}:
 Record of a number field
julia> K = number_field(lf[1])
Number field over Rational Field with defining polynomial x^2 - 2
julia> class_group(lf[1])
GrpAb: Z/1
julia> regulator(lf[1])
[0.8813735870195430252 +/- 6.76e-20]

julia> Qx, x = FlintQQ["x"];
julia> K, a = number_field(x^2+5, cached = false)
(Number field over Rational Field with defining polynomial x^2 + 5, _a)
julia> r = find_DBfield(db, K)
Record of a number field defined by x^2 + 5
julia> class_group(r)
GrpAb: Z/2
julia> h = hilbert_class_field(r)
Record of a number field
julia> number_field(h)
Number field over Rational Field with defining polynomial x^4 - 2*x^3 + x^2 + 5
julia> factor(discriminant(h))
1 * 5^2 * 2^4




```
