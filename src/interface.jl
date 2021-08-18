export fields_database, isregulator_known, isclass_group_known, isgalois_group_known, aresubfields_known,
        ishilbert_class_field_known, isnormal_closure_known, ramified_primes

Hecke.add_verbose_scope(:FieldsDB)

@doc Markdown.doc"""
    fields_database(password::String)-> Database

Returns the database.
"""
function fields_database(password::String = "")
  if isempty(password)
    return LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag")
  else
    return LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag password=" * password)
  end
end

mutable struct DBField
  id::Int
  connection::LibPQ.Connection
  degree::Int
  polynomial::fmpq_poly
  number_field::AnticNumberField
  class_group::GrpAbFinGen
  discriminant::fmpz
  ramified_primes::Vector{fmpz}
  
  regulator::arb
  signature::Tuple{Int, Int}

  galois_group::PermGroup
  transitive_id::Int

  subfields::Vector{DBField}
  normal_closure::DBField
  hilbert_class_field::DBField

  torsion_size::Int
  automorphisms_order::Int 
  #Those should be Bool, but there is a problem with the initialization.
  #Therefore, -1: Not set in the database
  #           0 : Not known
  #           1 : True
  #           2 : False 
  is_canonical_poly::Int
  GRH::Int 
  is_cm::Int 

  function DBField(conn::LibPQ.Connection, id::Int)
    z = new()
    z.connection = conn
    z.id = id
    z.degree = -1
    z.signature = (-1, 0)
    z.torsion_size = -1
    z.transitive_id = -1
    z.automorphisms_order = -1
    z.is_canonical_poly = 0
    z.GRH = 0
    z.is_cm = 0
    return z
  end
end

function LibPQ.pqparse(::Type{Vector{BigInt}}, x::String)
  if x == "{}"
    return BigInt[]
  end
  s = split(x[2:end-1], ",")
  return BigInt[parse(BigInt, split(ss, ".")[1]) for ss in s]
end

################################################################################
#
#  Hash
#
################################################################################

function Base.hash(x::DBField, h::UInt)
  return hash(x.id, h)
end

################################################################################
#
#  Show
#
################################################################################

function Base.show(io::IO, x::DBField)
  print(io, "Record of a number field")
  if isdefined(x, :polynomial)
    print(io, " defined by ")
    print(io, x.polynomial)
  end
  print(io, "\n")
  return nothing
end

################################################################################
#
#  IsSet function
#
################################################################################

function areramified_primes_known(x::DBField)
  query = "SELECT ramified_primes FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id], column_types = Dict(:ramified_primes => Vector{BigInt}))
  return result[1, 1] !== missing
end

function isregulator_known(x::DBField)
  query = "SELECT regulator FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id], column_types = Dict(:regulator => BigFloat))
  return result[1, 1] !== missing
end

function isclass_group_known(x::DBField)
  query = "SELECT class_group_id FROM field WHERE field_id = \$1"
  data = x.id
  result = execute(x.connection, query, [data])
  return result[1, 1] !== missing
end

function isgalois_group_known(x::DBField)
  query = "SELECT group_id FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing
end

function iscm_property_known(x::DBField)
  query = "SELECT CM FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing
end

function isautomorphism_order_known(x::DBField)
  query = "SELECT automorphisms_order FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing
end

function istorsion_unit_size_known(x::DBField)
  query = "SELECT torsion_size FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing
end

function aresubfields_known(x::DBField)
  query = "SELECT subfields FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing
end

function iscanonical_polynomial_known(x::DBField)
  query = "SELECT is_canonical_poly FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing 
end

function ishilbert_class_field_known(x::DBField)
  query = "SELECT hilbert_class_field FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing 
end

function isnormal_closure_known(x::DBField)
  query = "SELECT normal_closure FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  return result[1, 1] !== missing 
end


################################################################################
#
#  Importing infos from the database
#
################################################################################

@doc Markdown.doc"""
    defining_polynomial(x::DBField) -> fmpq_poly

Returns the polynomial defining the field corresponding to the record.
"""
function Oscar.defining_polynomial(x::DBField; cached::Bool = true)
  if cached && isdefined(x, :polynomial)
    return x.polynomial
  end
  query = "SELECT polynomial FROM field WHERE field_id = \$1"
  data = x.id
  result = execute(x.connection, query, [data], column_types = Dict(:polynomial => Vector{BigInt}))
  coeffs = map(fmpz, result[1, 1])
  Qx = PolynomialRing(FlintQQ, "x", cached = false)[1]
  x.polynomial = Qx(coeffs)
  return x.polynomial
end

@doc Markdown.doc"""
    degree(x::DBField) -> Int
Returns the degree of the field corresponding to the record.
"""
function Oscar.degree(x::DBField)
  if x.degree != -1
    return x.degree
  end
  query = "SELECT degree FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  x.degree = result[1, 1]
  return x.degree
end

@doc Markdown.doc"""
    number_field(x::DBField) -> AnticNumberField

Returns the corresponding field.
"""
function Oscar.number_field(x::DBField)
  if isdefined(x, :number_field)
    return x.number_field
  end
  f = defining_polynomial(x)
  x.number_field = number_field(f, cached = false, check = false)[1]
  return x.number_field
end

@doc Markdown.doc"""
    signature(x::DBField) -> Tuple{Int, Int}

Returns the signature of the corresponding field.
"""
function Oscar.signature(x::DBField)
  if x.signature[1] != -1
    return x.signature
  end
  n = degree(x)
  query = "SELECT real_embeddings FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  real_embs = result[1, 1]
  x.signature = (real_embs, divexact(n-real_embs, 2))
  return x.signature
end

@doc Markdown.doc"""
    discriminant(x::DBField) -> fmpz

Returns the discriminant of the maximal order of the corresponding field.
"""
function Oscar.discriminant(x::DBField)
  if isdefined(x, :discriminant)
    return x.discriminant
  end
  query = "SELECT discriminant FROM field WHERE field_id = \$1"
  data = x.id
  result = execute(x.connection, query, [data], column_types = Dict(:discriminant => BigInt))
  x.discriminant = fmpz(result[1, 1])
  return x.discriminant
end

@doc Markdown.doc"""
    ramified_primes(x::DBField) -> Vector{fmpz}

Returns the prime numbers dividing the discriminant of the maximal order of the corresponding field.
"""
function ramified_primes(x::DBField)
  if isdefined(x, :ramified_primes)
    return x.ramified_primes
  end
  query = "SELECT ramified_primes FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id], column_types = Dict(:ramified_primes => Vector{BigInt}))
  if result[1, 1] === missing
    return missing
  end
  res = map(fmpz, result[1, 1]::Vector{BigInt})
  x.ramified_primes = res
  return res
end

@doc Markdown.doc"""
    class_group(x::DBField) -> GrpAbFinGen

Returns, if known, an abstract group isomorphic to the class group of the corresponding field.
"""
function Oscar.class_group(x::DBField; cached::Bool = true)
  if cached  && isdefined(x, :class_group)
    return x.class_group
  end
  query = "SELECT class_group_id FROM field WHERE field_id = $(x.id)"
  result = execute(x.connection, query)
  if result[1, 1] === missing
    return missing
  end
  query1 = "SELECT structure FROM class_group WHERE class_group_id = $(result[1, 1])"
  result1 = execute(x.connection, query1, column_types = Dict(:structure => Vector{BigInt}))
  r = result1[1, 1]::Vector{BigInt}
  x.class_group = GrpAbFinGen(r)
  return x.class_group
end

@doc Markdown.doc"""
    regulator(x::DBField) -> arb

Returns, if known, the regulator of the corresponding field.
"""
function Oscar.regulator(x::DBField)
  if isdefined(x, :regulator)
    return x.regulator
  end
  query = "SELECT regulator FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id], column_types = Dict(:regulator => Decimal))
  if result[1, 1] === missing
    return missing
  end
  s = string(result[1, 1])
  R = ArbField(64, cached = false)
  return R(s)
end

@doc Markdown.doc"""
    galois_group(x::DBField) -> PermGroup

Returns, if known, a permutation group which is isomorphic (as a permutation group) to the Galois group of the corresponding field.
"""
function Oscar.galois_group(x::DBField)
  if isdefined(x, :galois_group)
    return x.galois_group
  end
  #I need to retrieve it from the database.
  query = "SELECT group_id FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id], column_types = Dict(:group_id => Int))
  if result[1, 1] === missing
    return missing
  end
  G = find_group(x.connection, result[1, 1])
  x.galois_group = G
  return G
end

@doc Markdown.doc"""
    is_cm(x::DBField) -> Bool

Returns true whether the corresponding field is a CM field, false otherwise.
"""
function is_cm(x::DBField)
  if x.is_cm != 0
    if x.is_cm == -1
      return missing
    elseif x.is_cm == 1
      return true
    else
      return false
    end
  end
  query = "SELECT CM FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = result[1, 1]
  if data === missing
    return missing
  elseif data
    x.is_cm = 1
    return true
  else
    x.is_cm = 2
    return false
  end
end

@doc Markdown.doc"""
    automorphisms_order(x::DBField) -> Int

Returns the order of the automorphism group of the corresponding field.
"""
function automorphisms_order(x::DBField)
  if x.automorphisms_order != -1
    return x.automorphisms_order
  end
  query = "SELECT automorphisms_order FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = result[1, 1]
  if data === missing
    return missing
  end
  x.automorphisms_order = data
  return data
end

@doc Markdown.doc"""
    torsion_units_size(x::DBField) -> Int

Returns the order of the torsion units of the corresponding field.
"""
function torsion_units_size(x::DBField)
  if x.torsion_size!= -1
    return x.torsion_size
  end
  query = "SELECT torsion_size FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = result[1, 1]
  if data === missing
    return missing
  end
  x.torsion_size = data
  return data
end

function assumes_GRH(x::DBField)
  if x.GRH != 0
    if x.GRH == -1
      return missing
    elseif x.GRH == 1
      return true
    else
      return false
    end
  end
  query = "SELECT GRH FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = result[1, 1]
  if data === missing
    x.GRH = -1
    return missing
  elseif data
    x.GRH = 1
    return true
  else
    x.GRH = 2
    return false
  end
end

function has_canonical_defining_polynomial(x::DBField)
  if x.is_canonical_poly != 0
    if x.is_canonical_poly == -1
      return missing
    elseif x.is_canonical_poly == 1
      return true
    else
      return false
    end
  end
  query = "SELECT is_canonical_poly FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = result[1, 1]
  if data === missing
    x.is_canonical_poly = -1
    return missing
  elseif data
    x.is_canonical_poly = 1
    return true
  else
    x.is_canonical_poly = 2
    return false
  end
end

@doc Markdown.doc"""
    subfields(x::DBField) -> Vector{DBField}

Returns, if known, the subfields of the corresponding field.
"""
function Oscar.subfields(x::DBField)
  if isdefined(x, :subfields)
    return x.subfields
  end
  query = "SELECT subfields FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  sub = result[1, 1]
  if sub === missing
    return missing
  end
  res = Vector{DBField}(undef, length(sub))
  for i = 1:length(res)
    res[i] = DBField(x.connection, sub[i])
  end
  return res
end

@doc Markdown.doc"""
    normal_closure(x::DBField) -> DBField

Returns, if known, the normal closure of the corresponding field.
"""
function Oscar.normal_closure(x::DBField)
  if isdefined(x, :normal_closure)
    return x.normal_closure
  end
  query = "SELECT normal_closure FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  sub = result[1, 1]
  if sub === missing
    return missing
  end
  res = DBField(x.connection, sub)
  return res
end

@doc Markdown.doc"""
    hilbert_class_field(x::DBField) -> DBField

Returns, if known, the hilbert_class_field of the corresponding field.
"""
function Oscar.hilbert_class_field(x::DBField)
  if isdefined(x, :hilbert_class_field)
    return x.hilbert_class_field
  end
  query = "SELECT hilbert_class_field FROM field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  sub = result[1, 1]
  if sub === missing
    return missing
  end
  res = DBField(x.connection, sub)
  return res
end

################################################################################
#
#  Set functions - Setting a property in the database 
#
################################################################################

function set_polynomial(x::DBField, f::fmpq_poly; is_canonical::Bool = false)
  pol = BigInt[BigInt(numerator(coeff(f, i))) for i = 0:degree(f)]
  query = "UPDATE field SET polynomial = \$1, is_canonical_poly = \$2 WHERE field_id = \$3"
  execute(x.connection, query, (pol, is_canonical, x.id))
  x.polynomial = f
  return nothing
end

function set_ramified_primes(x::DBField, lp::Vector{fmpz})
  rp = BigInt[BigInt(x) for x in lp]
  query = "UPDATE field  SET ramified_primes = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (rp, x.id))
  return nothing
end

function set_galois_group(x::DBField, G::PermGroup)
  id = _find_group_id(x.connection, G)
  if id === missing
    insert_group(x.connection, G)
    id = _find_group_id(x.connection, G)::Int
  end
  query = "UPDATE field  SET group_id = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (id, x.id))
  return nothing
end

function set_class_group(x::DBField, C::GrpAbFinGen; GRH::Bool = true)
  id = _find_class_group_id(x.connection, C)
  if id === missing
    insert_class_group(x.connection, C)
    id = _find_class_group_id(x.connection, C)
  end 
  if GRH
    x.GRH = 1
  else
    x.GRH = 2
  end
  query = "UPDATE field SET class_group_id = \$1, GRH = \$2  WHERE field_id = \$3"
  execute(x.connection, query, (id, GRH, x.id))
  return nothing
end

function set_CM_property(x::DBField, iscm::Bool)
  query = "UPDATE field  SET CM = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (iscm, x.id))
  return nothing
end

function set_torsion_size(x::DBField, n::Int)
  query = "UPDATE field  SET torsion_size = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (n, x.id))
  return nothing
end

function set_automorphisms_order(x::DBField, n::Int)
  query = "UPDATE field SET automorphisms_order = \$1 WHERE field_id = \$2"
  execute(x.connection, query, (n, x.id))
  return nothing
end

################################################################################
#
#  "Fill" functions
#
################################################################################

function set_canonical_defining_polynomial(x::DBField)
  K = number_field(x)
  K1 = simplify(K, canonical = true, cached = false, save_LLL_basis = false)[1]
  set_polynomial(x, defining_polynomial(K1), is_canonical = true)
  x.is_canonical_poly = 1
  return nothing
end

function set_ramified_primes(x::DBField)
  d = discriminant(x)
  rp = collect(keys(factor(d).fac))
  set_ramified_primes(x, rp)
end

function set_automorphisms_order(x::DBField)
  K = number_field(x)
  set_automorphisms_order(x, length(automorphisms(K)))
end

function set_torsion_unit_size(x::DBField)
  K = number_field(x)
  set_torsion_size(x, Hecke.torsion_units_order(K))
end

function set_iscm(x::DBField)
  K = number_field(x)
  set_CM_property(x, Hecke.iscm_field(K)[1])
end

function set_class_group(x::DBField, GRH::Bool = true)
  K = number_field(x)
  OK = maximal_order(K)
  C = class_group(OK, GRH = GRH)[1]
  set_class_group(x, C, GRH = GRH)
end

function set_regulator(x::DBField)
  if degree(x) == 2 && signature(x) == (0, 1)
    r = Decimal(1)
  else
    r = _regulator_as_decimal(number_field(x))
  end
  query = "UPDATE field SET regulator = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (r, x.id))
  return nothing
end

function set_class_group_and_regulator(x::DBField, GRH::Bool = true)
  K = number_field(x)
  OK = maximal_order(K)
  C = class_group(OK, GRH = GRH)[1]
  set_class_group(x, C, GRH = GRH)
  set_regulator(x)
end

function set_galois_group(x::DBField)
  G = galois_group(number_field(x))[1]
  set_galois_group(x, G)
end

function set_subfields(x::DBField)
  K = number_field(x)
  ao = automorphisms_order(x)
  if ao !== missing && ao == degree(K)
    lll(maximal_order(K))
    subfs = Hecke.subfields_normal(K, true)
  else
    subfs = subfields(K)
  end
  subs = AnticNumberField[y[1] for y in subfs]
  for i = 1:length(subs)
    subs[i] = simplify(subs[i], cached = false, save_LLL_basis = false)[1]
    @assert isdefining_polynomial_nice(subs[i])
  end
  lS = isomorphism_class_representatives(subs)
  ids = Vector{BigInt}(undef, length(lS))
  for i = 1:length(lS)
    K = lS[i]
    y = find_DBfield(x.connection, K)
    if y === missing
      insert_field(x.connection, K)
      y = find_DBfield(x.connection, K)::DBField
    end
    ids[i] = y.id
  end
  query = "UPDATE field SET subfields = \$1 WHERE field_id = \$2"
  execute(x.connection, query, (ids, x.id))
  return nothing
end

function set_normal_closure(x::DBField)
  f = defining_polynomial(x)
  L = splitting_field(f)
  insert_fields(x.connection, AnticNumberField[L])
  id = find_DBfield(x.connection, L)
  x.normal_closure = id
  query = "UPDATE field SET normal_closure = $(id.id) WHERE field_id = $(x.id)"
  execute(x.connection, query)
  return nothing
end

function set_hilbert_class_field(x::DBField)
  K = number_field(x)
  L = absolute_simple_field(number_field(hilbert_class_field(K)))[1]
  insert_fields(x.connection, AnticNumberField[L])
  id = find_DBfield(x.connection, L)
  x.normal_closure = id
  query = "UPDATE field SET hilbert_class_field = $(id.id) WHERE field_id = $(x.id)"
  execute(x.connection, query)
  return nothing
end

################################################################################
#
#  Auxiliary functions
#
################################################################################
@doc Markdown.doc"""
    isomorphic_transitive_perm_group(G::Oscar.GAPGroup, degree::Int) -> PermGroup

Given a group $G$, returns a permutation group $GP$ isomorphic to $G$ that is a
transitive subgroup of the symmetric group over 'degree' points. Throws an error if
such a group does not exist.
"""
function isomorphic_transitive_perm_group(G::Oscar.GAPGroup, deg::Int)
  oG = order(G)
  @assert iszero(mod(oG, deg))
  target_order = divexact(oG, deg)
  lC = conjugacy_classes_subgroups(G)
  for x in lC
    r = representative(x)
    if order(r) != target_order
      continue
    end
    C, mC = core(G, r)
    if order(C) != 1
      continue
    end
    mp = GAP.Globals.FactorCosetAction(G.X, r.X)
    H = PermGroup(GAP.Globals.Image(mp))
    return H
  end
  error("Required representation does not exist")
end

@doc Markdown.doc"""
    find_automorphisms_order(G::PermGroup) -> Int

Returns the order of the automorphism group of any field with Galois group $G$.
"""
function find_automorphisms_order(G::PermGroup)
  oG = order(G)
  deg = degree(G)
  @assert iszero(mod(oG, deg))
  target_order = divexact(oG, deg)
  lC = conjugacy_classes_subgroups(G)
  for x in lC
    r = representative(x)
    if order(r) != target_order
      continue
    end
    C, mC = core(G, r)
    if order(C) != 1
      continue
    end
    #I have the right group.
    return divexact(order(normalizer(G, r)[1]), target_order)
  end
  error("Something went wrong")
end

@doc Markdown.doc"""
    possible_signatures(G::PermGroup) -> Vector{Tuple{Int, Int}}

Returns the signatures of any field with Galois group $G$.
"""
function possible_signatures(G::PermGroup)
  oG = order(G)
  deg = degree(G)
  @assert iszero(mod(oG, deg))
  target_order = divexact(oG, deg)
  lC = conjugacy_classes_subgroups(G)
  ind_right = 0
  involutions = typeof(lC)()
  for i = 1:length(lC)
    r = representative(lC[i])
    if order(r) == 2
      push!(involutions, lC[i])
    end
    if ind_right != 0
      continue
    end
    if order(r) != target_order
      continue
    end
    C, mC = core(G, r)
    if order(C) == 1
      ind_right = i
    end
  end

  #Now, I have the conjugacy classes of involutions and the subgroup H.
  #Given a conjugacy class, the number of real places of the fixed field is equal
  #to the number of conjugates contained in H
  possible_real_plcs = Set{Int}()
  H = representative(lC[ind_right])
  corr = order(H)
  for i = 1:length(involutions)
    elts = collect(involutions[i])
    rp = 0
    for j = 1:length(elts)
      I = intersect(H, elts[j])[1]
      if order(I) == 2
        rp += 1
      end
    end 
    push!(possible_real_plcs, Int(divexact(rp*order(G), order(H)*length(elts))))
  end
  d = degree(G)
  possible_sigs = Vector{Tuple{Int, Int}}()
  push!(possible_sigs, (d, 0))
  for r1 in possible_real_plcs
    r2 = div(d - r1, 2)
    push!(possible_sigs, (r1, r2))
  end
  return possible_sigs
end

function _pdtype_shape(OK::NfOrd, p::Int)
  pd = prime_decomposition_type(OK, p)
  res = Dict{Tuple{Int, Int}, Int}()
  for x in pd
    if haskey(res, x)
      res[x] += 1
    else
      res[x] = 1
    end
  end
  return res
end

@doc Markdown.doc"""
    isomorphism_class_representatives(fields::Vector{AnticNumberField}) -> Vector{AnticNumberField}

Returns a vector of pairwise non isomorphic fields that are representatives for the isomorphism class
of the fields contained in 'fields'.
"""
function isomorphism_class_representatives(v::Vector{AnticNumberField})
  #First, I use signature and discriminant as sieving parameters 
  first_sieve = Dict{Tuple{Int, Int, fmpz}, Vector{Int}}()
  for i = 1:length(v)
    x = v[i]
    d = degree(x)
    sign = signature(x)[1]
    disc = discriminant(maximal_order(x))
    if haskey(first_sieve, (d, sign, disc))
      push!(first_sieve[(d, sign, disc)], i)
    else
      first_sieve[(d, sign, disc)] = Int[i]
    end
  end
  #Now, we sieve by prime splitting
  clusters = Vector{Vector{Int}}()
  for x in values(first_sieve)
    if length(x) == 1
      push!(clusters, x)
      continue
    end
    S = PrimesSet(100, -1)
    ind = 0
    cx = Vector{Int}[x]
    for p in S
      cx_new = Vector{Int}[]
      for y in cx
        shapes = [_pdtype_shape(maximal_order(v[t]), p) for t in y]
        sshapes = unique(shapes)
        if length(sshapes) == 1
          push!(cx_new, y)
          continue
        end
        #We have length(sshapes) iso classes (at least)
        new_clusters = Vector{Vector{Int}}(undef, length(sshapes))
        for i = 1:length(new_clusters)
          new_clusters[i] = Vector{Int}()
        end
        for j = 1:length(y)
          for k = 1:length(sshapes)
            if shapes[j] == sshapes[k]
              push!(new_clusters[k], y[j])
              break
            end
          end
        end   
        append!(cx_new, new_clusters)     
      end
      cx = cx_new
      ind += 1
      if ind == 10
        append!(clusters, cx)
        break
      end
    end
  end
  #Now, check isomorphisms.
  res = Vector{AnticNumberField}()
  for x in clusters
    if length(x) == 1
      push!(res, v[x[1]])
      continue
    end
    push!(res, v[x[1]])
    for j = 2:length(x)
      found = false
      for k = 1:j-1
        if isisomorphic(v[x[j]], v[x[k]])[1]
          found = true
          break
        end
      end
      if !found
        push!(res, v[x[j]])
      end
    end
  end
  return res
end

################################################################################
#
#  To populate the database
#
################################################################################

function _get_fields_for_class_group_computation(connection::LibPQ.Connection)
  query = "SELECT field_id FROM field WHERE class_group_id IS NULL LIMIT 20"
  result = Tables.rows(execute(connection, query))
  res = Vector{DBField}(undef, 20)
  ind = 1
  for x in result
    res[ind] = DBField(connection, x[1])
    ind += 1
  end
  return res  
end

function _get_fields_for_class_group_computation(connection::LibPQ.Connection, degree::Int)
  query = "SELECT field_id FROM field WHERE degree = $degree AND class_group_id IS NULL LIMIT 20"
  result = Tables.rows(execute(connection, query))
  res = Vector{DBField}(undef, 20)
  ind = 1
  for x in result
    res[ind] = DBField(connection, x[1])
    ind += 1
  end
  return res  
end


function _get_fields_for_subfields_computation(connection::LibPQ.Connection)
  query = "SELECT field_id FROM field WHERE subfields IS NULL LIMIT 20"
  result = Tables.rows(execute(connection, query))
  res = Vector{DBField}(undef, 20)
  ind = 1
  for x in result
    res[ind] = DBField(connection, x[1])
    ind += 1
  end
  return res  
end

function _get_fields_for_subfields_computation(connection::LibPQ.Connection, degree::Int)
  query = "SELECT field_id FROM field WHERE degree = $degree AND subfields IS NULL LIMIT 20"
  result = Tables.rows(execute(connection, query))
  res = Vector{DBField}(undef, 20)
  ind = 1
  for x in result
    res[ind] = DBField(connection, x[1])
    ind += 1
  end
  return res  
end

function _get_fields_for_ramified_primes_computation(connection::LibPQ.Connection)
  query = "SELECT field_id FROM field WHERE ramified_primes IS NULL"
  result = Tables.rows(execute(connection, query))
  res = Vector{DBField}()
  for x in result
    push!(res, DBField(connection, x[1]))
  end
  return res  
end

################################################################################
#
#  Helper functions for regulator
#
################################################################################

function _string(x::arb, digits::Int)
   cstr = ccall((:arb_get_str, Hecke.Nemo.libarb), Ptr{UInt8},
                (Ref{arb}, Int, UInt),
                x, Int(digits), UInt(2))
   res = unsafe_string(cstr)
   ccall((:flint_free, Hecke.Nemo.libflint), Nothing,
         (Ptr{UInt8},),
         cstr)
   return res
end

function _regulator_as_string(K::AnticNumberField, scale::Int = 20)
  p = Int(ceil(3.3 * (scale + 2)))
  # make radius less than 10^(-(scale + 2))
  U, mU = unit_group_fac_elem(maximal_order(K))
  r = regulator([mU(U[i]) for i in 2:ngens(U)], p)
  return _string(r, scale)
end

function _regulator_as_decimal(K::AnticNumberField, scale::Int = 20)
  return decimal(_regulator_as_string(K, scale))
end

