export fields_database, load_fields, insert_field, insert_fields, insert_complete_table

function fields_database()
  return LibPQ.Connection("host=localhost dbname=fields port=5432 user=postgres")
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

  subfields::Vector

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
    print(io, "defined by ")
    print(io, x.polynomial)
  end
  print(io, "\n")
  return nothing
end

################################################################################
#
#  Importing infos from the database
#
################################################################################

function LibPQ.pqparse(::Type{Vector{BigInt}}, x::String)
  s = split(x[2:end-1], ",")
  return [parse(BigInt, split(ss, ".")[1]) for ss in s]
end

function Oscar.defining_polynomial(x::DBField)
  if isdefined(x, :polynomial)
    return x.poly
  end
  query = "SELECT polynomial FROM fields.field WHERE field_id = \$1"
  data = x.id
  result = execute(x.connection, query, [data], column_types = Dict(:polynomial => Vector{BigInt}))
  data = columntable(result)[1][1]
  @show data
  coeffs = Vector{fmpz}(undef, length(data))
  for i = 1:length(coeffs)
    coeffs[i] = fmpz(BigInt(data[i]))
  end
  Qx = PolynomialRing(FlintQQ, "x", cached = false)[1]
  x.polynomial = Qx(coeffs)
  return x.polynomial
end

function Oscar.degree(x::DBField)
  if x.degree != -1
    return x.degree
  end
  query = "SELECT degree FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  x.degree = columntable(result)[1][1]
  return x.degree
end

function Oscar.number_field(x::DBField)
  if isdefined(x, :number_field)
    return x.number_field
  end
  f = defining_polynomial(x)
  x.number_field = number_field(f, cached = false, check = false)[1]
  return x.number_field
end

function Oscar.signature(x::DBField)
  if x.signature[1] != -1
    return x.signature
  end
  n = degree(x)
  query = "SELECT real_embeddings FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.data])
  real_embs = columntable(result)[1][1]
  x.signature = (real_embs, divexact(n-real_embs, 2))
  return x.signature
end

function Oscar.discriminant(x::DBField)
  if isdefined(x, :discriminant)
    return x.discriminant
  end
  query = "SELECT discriminant FROM fields.field WHERE field_id = \$1"
  data = x.id
  result = execute(x.connection, query, [data], column_types = Dict(:discriminant => BigInt))
  x.discriminant = fmpz(BigInt(columntable(result)[1][1]))
  return x.discriminant
end

function ramified_primes(x::DBField)
  if isdefined(x, :ramified_primes)
    return x.ramified_primes
  end
  query = "SELECT ramified_primes FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id], column_types = Dict(:ramified_primes => Vector{BigInt}))
  tb = columntable(result)[1][1]
  if tb === missing
    return missing
  end
  res = Vector{fmpz}(undef, length(tb))
  for i = 1:length(tb)
    res[i] = fmpz(tb[i])
  end
  x.ramified_primes = res
  return res
end

function Oscar.class_group(x::DBField)
  if isdefined(x, :class_group)
    return x.class_group
  end
  query = "SELECT class_group_id FROM fields.field WHERE field_id = \$1"
  data = x.id
  result = execute(x.connection, query, [data])
  tb = columntable(result)[1][1]
  if tb === missing
    return missing
  end
  query1 = "SELECT structure FROM fields.class_group WHERE class_group_id = \$1"
  result1 = execute(x.connection, query1, [tb], column_types = Dict(:structure => Vector{BigInt}))
  str = columntable(result)[1][1]
  invs = Vector{fmpz}(undef, length(str))
  for i = 1:length(invs)
    invs[i] = fmpz(str[i])
  end
  x.class_group = abelian_group(invs)
  return x.class_group
end

function Oscar.regulator(x::DBField)
  if isdefined(x, :regulator)
    return x.regulator
  end
  query = "SELECT regulator " * "FROM fields.field" * " WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id], column_types = Dict(:regulator => BigFloat))
  data = columntable(result)[1][1]
  if data === missing
    return missing
  end
  R = ArbField(64, cached = false)
  return R(data)
end

function Oscar.galois_group(x::DBField)
  if isdefined(x, :galois_group)
    return x.galois_group
  end
  #I need to retrieve it from the database.
  query = "SELECT group_id FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = columntable(result)[1][1]
  if data === missing
    return missing
  end
  query1 = "SELECT transitive_group_id FROM fields.group WHERE group_id = \$1"
  result1 = execute(x.connection, query1, [data])
  data1 = columntable(result1)[1][1]
  if data1 != missing
    x.galois_group = transitive_group(degree(x), data1)
    return x.galois_group
  end
  query2 = "SELECT group_order, small_group_id FROM fields.group WHERE group_id = \$1"
  result2 = execute(x.connection, query2, [data])
  data2 = columntable(result2)
  if data2[2][1] !== missing
    #I have a polycyclic group and I want to get a perm group of the right degree.
    PC = small_group(data[1][1], data[2][1])
    H = isomorphic_transitive_perm_group(PC, degree(x))
    return H
  end
  query3 = "SELECT generators FROM fields.group WHERE group_id = \$1"
  result3 = execute(x.connection, query3, [data])
  data3 = columntable(result2)[1][1]
  S = symmetric_group(degree(x))
  g = eval(Meta.parse(data3))
  perms = Vector{PermGroupElem}(undef, length(g))
  for i = 1:length(g)
    perms[i] = S(g[i])
  end
  H, mH = sub(S, perms)
  x.galois_group = H
  return H
end


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
  query = "SELECT CM FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = columntable(result)[1][1]
  if data === missing
    return missing
  end
  if data
    x.is_cm = 1
    return true
  else
    x.is_cm = 2
    return false
  end
end

function automorphisms_order(x::DBField)
  if x.automorphisms_order != -1
    return x.automorphisms_order
  end
  query = "SELECT automorphisms_order FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = columntable(result)[1][1]
  if data === missing
    return missing
  end
  x.automorphisms_order = data
  return data
end

function torsion_units_size(x::DBField)
  if x.torsion_size!= -1
    return x.torsion_size
  end
  query = "SELECT torsion_size FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = columntable(result)[1][1]
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
  query = "SELECT GRH FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = columntable(result)[1][1]
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
  query = "SELECT is_canonical_poly FROM fields.field WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  data = columntable(result)[1][1]
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

function Oscar.subfields(x::DBField)
  if isdefined(x, :subfields)
    return x.subfields::Vector{DBFields}
  end
  query = "SELECT subfields " * "FROM fields.field" * " WHERE field_id = \$1"
  result = execute(x.connection, query, [x.id])
  sub = columntable(result)[1]
  if data[1] === missing
    return missing
  end
  v = data[1]
  res = Vector{DBFields}(undef, length(v))
  for i = 1:length(v)
    res[i] = DBField(x.connection, v[i])
  end
  return res
end

################################################################################
#
#  Load function
#
################################################################################

function load_fields(conn::LibPQ.Connection; degree_range::Tuple{Int, Int} = (0, -1), discriminant_range::Tuple{fmpz, fmpz} = (fmpz(0), fmpz(-1)), 
                         signature::Tuple{Int, Int} = (-1, 0), unramified_outside::Vector{fmpz} = fmpz[-1], ramified_at::Vector{fmpz} = fmpz[-1],
                         galois_group::PermGroup = symmetric_group(1), class_number::Int = -1, 
                         class_group_structure::Vector{fmpz} = fmpz[-1], 
                         class_group_ranks_range::Dict{fmpz, Tuple{Int, Int}} = Dict{fmpz, Tuple{Int, Int}}(), only_count::Type{Val{T}} = Val{false}) where T


  parameters = String[]
  values = []
  ind = 1
  if degree_range[2] != -1
    if degree_range[1] == degree_range[2]
      push!(parameters, "degree = \$$(ind)")
      push!(values, degree_range[1])
      ind += 1
    else
      push!(parameters, "degree BETWEEN \$$(ind) AND \$$(ind+1) ")
      push!(values, degree_range[1])
      push!(values, degree_range[2])
      ind += 2
    end
  end
  if discriminant_range[2] != -1
    if discriminant_range[1] == discriminant_range[2]
      push!(parameters, "discriminant = \$$(ind)")
      push!(values, BigInt(discriminant_range[1]))
      ind += 1
    else
      push!(parameters, "discriminant BETWEEN \$$(ind) AND \$$(ind+1)")
      push!(values, BigInt(discriminant_range[1]))
      push!(values, BigInt(discriminant_range[2]))
      ind += 2
    end
  end
  if signature[1] != -1
    push!(parameters, "real_embeddings = \$$(ind)")
    push!(values, signature[1])
    ind += 1
  end
  if unramified_outside[1] != -1
    push!(parameters, "ramified_primes <@ \$$(ind)")
    push!(values, [BigInt(x) for x in unramified_outside])
    ind += 1
  end
  if ramified_at[1] != -1
    push!(parameters, "ramified_primes @> \$$(ind)")
    push!(values, [BigInt(x) for x in ramified_at])
    ind += 1
  end
  if !isone(order(galois_group))
    id_group = _find_group_id(conn, galois_group)
    if id_group === missing
      return DBField[]
    end
    push!(parameters, "group_id = \$$(ind)")
    push!(values, id_group)
    ind += 1
  end
  if class_group_structure[1] != -1
    id_class_group = _find_class_group_id(conn, abelian_group(class_group_structure))
    if id_class_group !== missing
      push!(parameters, "class_group_id = \$$(ind)")
      push!(values, id_class_group)
      ind += 1
    end
  end
  if class_number != -1
    abgroups = abelian_groups(class_number)
    class_group_ids = Vector{Int}()
    for C in abgroups
      idC = _find_class_group_id(conn, C)
      if idC !== missing
        push!(class_group_ids, idC)
      end
    end
    if isempty(class_group_ids)
      return Vector{DBField}()
    end
    push!(parameters, "class_group_id = ANY( \$$(ind) ) ")
    push!(values, class_group_ids)
    ind += 1
  end
  if !isempty(class_group_ranks_range)
    ids_class_group = _find_class_group_ids(conn, class_group_ranks_range)
    push!(parameters, "class_group_id = ANY( \$$(ind)) ")
    push!(values, ids_class_group)
    ind += 1
  end
  if discriminant_range[2] != -1
    #I want to print the completeness data.
    if order(galois_group) != 1
      if signature[1] != -1
        res = find_completeness_data(conn, galois_group, signature)
        if res === missing
          println("The data might not be complete, even assuming GRH")
        end
        if res[2] > max(abs(discriminant_range[1]), abs(discriminant_range[2]))
          println("The data is complete even without assuming GRH")
        elseif res[1] > max(abs(discriminant_range[1]), abs(discriminant_range[2]))
          println("The data is complete assuming GRH")
        end
      else
        ps = possible_signatures(galois_group)
        res = find_completeness_data(conn, galois_group, ps[1])
        if res === missing
          println("The data might not be complete, even assuming GRH")
        else
          ismissing = false
          for i = 2:length(ps)
            resi = find_completeness_data(conn, galois_group, ps[i])
            if resi === missing
              ismissing = true
              break
            end
            res[1] = min(res[1], resi[1])
            res[2] = min(res[2], resi[2])
          end
          if ismissing
            println("The data might not be complete, even assuming GRH")
          elseif res[2] > max(abs(discriminant_range[1]), abs(discriminant_range[2]))
            println("The data is complete even without assuming GRH")
          elseif res[1] > max(abs(discriminant_range[1]), abs(discriminant_range[2]))
            println("The data is complete assuming GRH")
          end
        end
      end
    end
  end

  #Now, I can do the query
  if only_count == Val{true}
    query = "SELECT COUNT(*) FROM fields.field"
  else
    query = "SELECT field_id FROM fields.field"
  end
  if !isempty(values)
    query *= " WHERE "
    for i = 1:length(parameters)-1
      query = query * parameters[i] * " AND "
    end
    query =  query * parameters[end]
  end
  result = execute(conn, query, values)
  data = columntable(result)[1]
  if only_count == Val{true}
    return data[1]
  end
  fields = Vector{DBField}(undef, length(data))
  for j = 1:length(fields)
    fields[j] = DBField(conn, data[j])
  end
  return fields
end

function _find_group_id(conn::LibPQ.Connection, G::PermGroup)
  d = degree(G)
  id = transitive_group_identification(G)
  if id != -1
    #GREAT! Unique identification.
    query = "SELECT group_id FROM fields.group WHERE degree = \$1 AND transitive_group_id = \$2"
    values = Int[d, id]
    result = execute(conn, query, values)
    data = columntable(result)[1][1]
    return data
  end
  
  o = order(G)
  if o < 2000 && o != 1024
    #Now, I try the small group id
    id = small_group_identification(G)
    query = "SELECT group_id FROM fields.group WHERE degree = \$1 AND group_order = \$2 AND small_group_id = \$3"
    values = [d, order, id[2]]
    result = execute(conn, query, values)
    data = columntable(result)[1][1]
    return data
  end
  #Sad. That's quite hard now. I have to search for all the groups
  #with a given order and check isomorphism.
  query = "SELECT group_id, generators FROM fields.group WHERE degree = \$1 AND group_order = \$2"
  result = execute(conn, query, [d, o])
  data = Tables.rows(result)
  ind = 0
  S = symmetric_group(deg)
  for r in data
    vects = eval(Meta.parse(r[2]))
    perms = PermGroupElem[S(x) for x in vects]
    H, mH = sub(S, perms)
    if isisomorphic(H, G)
      return r[1]
    end
  end
  return missing
end

function _find_class_group_id(conn::LibPQ.Connection, C::GrpAbFinGen)
  if isone(order(C))
    query = "SELECT class_group_id FROM fields.class_group WHERE group_order = \$1"
    result = execute(conn, query, [1])
    return columntable(result)[1][1]
  end
  invs = snf(C)[1].snf
  query = "SELECT class_group_id FROM fields.class_group WHERE structure = \$1"
  result = execute(conn, query, [invs])
  return columntable(result)[1][1]
end

function _find_class_group_ids(conn::LibPQ.Connection, ranks::Dict{fmpz, Tuple{Int, Int}})
  #I want to find the class group ids of the class group having the rank as required.
  divs = [BigInt(x) for x in keys(ranks)]
  sort!(divs)
  query = "SELECT class_group_id, prime_divisors, ranks FROM fields.class_group WHERE prime_divisors @> \$1"
  result = execute(conn, query, [divs])
  data = Tables.rows(result)
  #I assume that the divisors are ordered.
  ids = Vector{Int}()
  for x in data
    divsi = [BigInt(y) for y in x[2]]
    ranksi = x[3]
    acceptable = true
    ind = 1
    for i = 1:length(divs)
      while divsi[ind] < divs[i]
        ind += 1
      end
      @assert divsi[ind] == divs[i]
      vs = ranks[fmpz(divs[i])]
      if vs[1] > ranksi[ind] || vs[2] < ranksi[ind] 
        acceptable = false
        break
      end
    end
    if acceptable
      push!(ids, x[1])
    end
  end
  return ids
end

function find_completeness_data(conn::LibPQ.Connection, G::PermGroup, signature::Tuple{Int, Int})
  idG = _find_group_id(conn, G)
  @assert idG !== missing
  query = "SELECT GRH, discriminant_bound FROM fields.completeness WHERE group_id = \$1 AND real_embeddings = \$2"
  result = execute(conn, query, [idG, signature[1]])
  data = Tables.rows(result)
  if isempty(data)
    return missing
  end
  #First result under GRH, second without.
  res = Vector{fmpz}(undef, 2)
  for x in data
    if x[1]
      res[1] = fmpz(BigInt(x[2]))
    else
      res[2] = fmpz(BigInt(x[2]))
    end
  end
  return res
end

function find_DBfield(conn::LibPQ.Connection, K::AnticNumberField)
  d = degree(K)
  disc = discriminant(maximal_order(K))
  lf = load_fields(conn, degree_range = (d, d), discriminant_range = (disc, disc), signature = signature(K))
  if isempty(lf)
    return missing
  end
  for x in lf
    if isisomorphic(K, number_field(x))[1]
      return x
    end
  end
  return missing
end

################################################################################
#
#  Insertion
#
################################################################################

function insert_complete_table(connection::LibPQ.Connection, fields::Vector{AnticNumberField}, galois_group::PermGroup, discriminant_bound::fmpz, GRH::Bool = true, signature::Tuple{Int, Int} = (-1, 0); check::Bool = true)
  g_id = _find_group_id(connection, galois_group)
  if g_id === missing
    insert_group(connection, galois_group)
    g_id = _find_group_id(connection, galois_group)
  end
  #I compute the number of automorphisms of the fields
  #If H is the subgroup fixing K, the number of automorphisms is equal to
  #the order normalizer of H divided by the order of H
  aut_order = Int(find_automorphisms_order(galois_group))
  execute(conn, "BEGIN;")
  real_embs = signature[1]
  have_signature = (real_embs != -1)
  for K in fields
    d = discriminant(maximal_order(K))
    if !have_signature
      real_embs = Hecke.signature(K)[1]
    end
    pol = BigInt[BigInt(numerator(coeff(K.pol, i))) for i = 0:degree(K)]
    deg = degree(K)
    if check
      lf = load_fields(connection, degree_range = (degree(K), degree(K)), discriminant_range = (d, d), signature = Hecke.signature(K))
      if !isempty(lf)
        found = false
        for x in lf
          if isisomorphic(K, number_field(x))[1]
            println("field already in database!")
            found = true
            break
          end
        end
        if found
          continue
        end
      end
    end
    LibPQ.load!(
      (real_embeddings = [real_embs], 
      polynomial = [pol], 
      discriminant = [BigInt(d)], 
      degree = [deg],
      group_id = [g_id],
      automorphisms_order = [aut_order]
      ),
      connection,
      "INSERT INTO fields.field (
        real_embeddings, 
        polynomial,
        discriminant, 
        degree,
        group_id,
        automorphisms_order
      ) VALUES (\$1, \$2, \$3, \$4, \$5, \$6);",
    )
  end
  execute(conn, "COMMIT;")
  #Finally, we save the completeness data
  if have_signature
    insert_completeness_data(connection, galois_group, signature, discriminant_bound, GRH)
  else
    for sg in possible_signatures(galois_group)
      insert_completeness_data(connection, galois_group, sg, discriminant_bound, GRH)
    end
  end
  return nothing
end


function insert_fields(fields::Vector{AnticNumberField}, conn; check::Bool = true, galois_group = symmetric_group(1))
  if order(galois_group) > 1
    g_id = _find_group_id(connection, galois_group)
    if g_id === missing
      insert_group(connection, galois_group)
      g_id = _find_group_id(connection, galois_group)
    end
  else
    g_id = 0
  end
  
  execute(conn, "BEGIN;")
  for K in fields
    @assert isdefining_polynomial_nice(K)
    d = discriminant(maximal_order(K))
    real_embs = Hecke.signature(K)[1]
    pol = BigInt[BigInt(numerator(coeff(K.pol, i))) for i = 0:degree(K)]
    deg = degree(K)
    if check
      lf = load_fields(conn, degree_range = (degree(K), degree(K)), discriminant_range = (d, d), signature = signature(K))
      if !isempty(lf)
        found = false
        for x in lf
          if isisomorphic(K, number_field(x))[1]
            println("field already in database!")
            found = true
            break
          end
        end
        if found
          continue
        end
      end
    end
    if !iszero(g_id)
      LibPQ.load!(
        (real_embeddings = [real_embs], 
        polynomial = [pol], 
        discriminant = [BigInt(d)], 
        degree = [deg],
        group_id = [g_id]
        ),
        conn,
        "INSERT INTO fields.field (
          real_embeddings, 
          polynomial,
          discriminant, 
          degree,
          group_id
        ) VALUES (\$1, \$2, \$3, \$4, \$5);",
      )
    else
      LibPQ.load!(
        (real_embeddings = [real_embs], 
        polynomial = [pol], 
        discriminant = [BigInt(d)], 
        degree = [deg]
        ),
        conn,
        "INSERT INTO fields.field (
          real_embeddings, 
          polynomial,
          discriminant, 
          degree
        ) VALUES (\$1, \$2, \$3, \$4);",
      )
    end
  end
  execute(conn, "COMMIT;")
  return nothing
end

function insert_field(K::AnticNumberField, connection::LibPQ.Connection; check::Bool = true)
  return insert_fields(AnticNumberField[K], connection, check = check)
end

function insert_class_group(conn::LibPQ.Connection, C::GrpAbFinGen)
  o = BigInt(order(C))
  str = map(BigInt, snf(C)[1].snf)
  lf = factor(str[end])
  divs = BigInt[BigInt(x) for x in keys(lf.fac)]
  sort!(divs)
  ranks = Vector{Int}(undef, length(divs))
  for i = 1:length(divs)
    ind = 0
    for j = 1:length(str)
      if iszero(mod(str[j], divs[i]))
        ind = j
        break
      end
    end
    ranks[i] = length(str)-ind+1
  end
  LibPQ.load!(
    (group_order = [o], 
    structure = [str],
    prime_divisors = [divs],
    ranks = [ranks], 
    ),
    conn,
    "INSERT INTO fields.class_group (
      group_order, 
      structure,
      prime_divisors,
      ranks
    ) VALUES (\$1, \$2, \$3, \$4);"
  )
  return nothing
end

function insert_group(conn::LibPQ.Connection, G::PermGroup)
  o = order(G)
  isab = isabelian(G)
  issolv = issolvable(G)
  issimp = issimple(G)
  isnilp = isnilpotent(G) 
  isperf = isperfect(G)
  isprim = isprimitive(G)
  d = degree(G)
  if d <= 33
    id = transitive_group_identification(G)
    LibPQ.load!(
      (group_order = [o], 
      degree = [d],
      transitive_group_id = [id],
      abelian = [isab], 
      nilpotent = [isnilp], 
      solvable = [issolv],
      issimple = [issimp],
      perfect = [isperf],
      primitive = [isprim]
      ),
      conn,
      "INSERT INTO fields.group (
        group_order, 
        degree,
        transitive_group_id,
        abelian, 
        nilpotent, 
        solvable,
        issimple,
        perfect,
        primitive
      ) VALUES (\$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9);"
    )
    return nothing
  end
  if order < 2000 && order != 1024
    id = small_group_identification(G)
    LibPQ.load!(
      (group_order = [o], 
      degree = [d],
      small_group_id = [id[2]],
      abelian = [isab], 
      nilpotent = [isnilp], 
      solvable = [issolv],
      issimple = [issimp],
      perfect = [isperf],
      primitive = [isprim]
      ),
      conn,
      "INSERT INTO fields.group (
        group_order, 
        degree,
        small_group_id,
        abelian, 
        nilpotent, 
        solvable,
        issimple,
        perfect,
        primitive
      ) VALUES (\$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9);"
    )
    return nothing
  end
  #Bad case. Cache a good set of generators.
  g = gens(G)
  s = "["
  for i = 1:length(gens)-1
    vi = Int[g[i][j] for j = 1:degree(G)]
    s = s * "$vi , "
  end
  vi = Int[g[length(gens)][j] for j = 1:degree(G)]
  s = s * "$vi]"
  LibPQ.load!(
    (group_order = [o], 
    degree = [d],
    generators = [s],
    abelian = [isab], 
    nilpotent = [isnilp], 
    solvable = [issolv],
    issimple = [issimp],
    perfect = [isperf],
    primitive = [isprim]
    ),
    conn,
    "INSERT INTO fields.group (
      group_order, 
      degree,
      transitive_group_id,
      abelian, 
      nilpotent, 
      solvable,
      issimple,
      perfect, 
      primitive
    ) VALUES (\$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9);"
  )
  return nothing
end

function insert_completeness_data(conn::LibPQ.Connection, group::PermGroup, signature::Tuple{Int, Int}, discriminant_bound::fmpz, GRH::Bool)
  gid = _find_group_id(conn, group)
  @assert gid !== missing
  LibPQ.load!(
    (group_id = [gid], 
    GRH = [GRH],
    real_embeddings = [s],
    discriminant_bound = BigInt(discriminant_bound)
    ),
    conn,
    "INSERT INTO fields.completeness (
      group_id, 
      GRH,
      real_embeddings,
      discriminant_bound
    ) VALUES (\$1, \$2, \$3, \$4);"
  )
  return nothing
end

################################################################################
#
#  Set functions - Setting a property in the database 
#
################################################################################

function set_polynomial(x::DBField, f::fmpq_poly; is_canonical::Bool = false)
  query = "UPDATE fields.field SET polynomial = \$1, is_canonical_poly = \$2 WHERE field_id = \$3"
  pol = BigInt[BigInt(coeff(f, i)) for i = 0:degree(f)]
  execute(x.connection, query, (pol, is_canonical, x.id))
  return nothing
end

function set_ramified_primes(x::DBField, lp::Vector{fmpz})
  rp = BigInt[BigInt(x) for x in lp]
  query = "UPDATE fields.field  SET ramified_primes = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (rp, x.id))
  return nothing
end

function set_galois_group(x::DBField, G::PermGroup)
  id = _find_group_id(x.connection, G)
  if id === missing
    insert_group(x.connection, G)
    id = _find_group_id(x.connection, G)
  end
  query = "UPDATE fields.field  SET group_id = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (id, x.id))
  return nothing
end

function set_regulator(x::DBField, r::arb)
  query = "UPDATE fields.field SET regulator = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (BigFloat(r), x.id))
  return nothing
end

function set_class_group(x::DBField, C::GrpAbFinGen; GRH::Bool = true)
  id = _find_class_group_id(x.connection, C)
  if id === missing
    insert_class_group(x.connection, C)
    id = _find_class_group_id(x.connection, C)
  end 
  query = "UPDATE fields.field  SET class_group_id = \$1, GRH = \$2  WHERE field_id = \$3"
  execute(x.connection, query, (id, GRH, x.id))
  return nothing
end

function set_CM_property(x::DBField, iscm::Bool)
  query = "UPDATE fields.field  SET CM = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (iscm, x.id))
  return nothing
end

function set_torsion_size(x::DBField, n::Int)
  query = "UPDATE fields.field  SET torsion_size = \$1  WHERE field_id = \$2"
  execute(x.connection, query, (n, x.id))
  return nothing
end

function set_automorphisms_order(x::DBField, n::Int)
  query = "UPDATE fields.field SET automorphisms_order = \$1 WHERE field_id = \$2"
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
  set_torsion_size(x, torsion_units_order(K))
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
  set_regulator(x, regulator(number_field(x)))
end

function set_galois_group(x::DBField)
  G = galois_group(number_field(x))
  set_galois_group(x, G)
end

function set_subfields(x::DBField)
  lS = isomorphism_class_representatives(AnticNumberField[y[1] for y in subfields(number_field(x))])
  ids = Vector{BigInt}(undef, length(lS))
  for i = 1:length(lS)
    insert_field(lS[i], x.connection, check = true)
    y = find_DBfield(lS[i])
    ids[i] = y.id
  end
  query = "UPDATE fields.field SET subfields = \$1 WHERE field_id = \$2"
  execute(x.connection, query, (ids, x.id))
  return nothing
end


################################################################################
#
#  Auxiliary functions
#
################################################################################

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
    elts = elements(involutions[i])
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
      res += 1
    else
      res = 1
    end
  end
  return res
end

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

  query = "SELECT field_id FROM fields.field WHERE class_group_id = \$1 LIMIT 20"
  result = rows(execute(connection, query, [missing]))
  res = Vector{DBField}(undef, 20)
  ind = 1
  for x in result
    res[i] = DBField(connection, x[1])
  end
  return res  
end


