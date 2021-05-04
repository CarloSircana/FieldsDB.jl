################################################################################
#
#  Load function
#
################################################################################

function load_fields_with_discriminant(db::LibPQ.Connection, discriminants::Vector{fmpz}, galois_group::PermGroup)
  discs = BigInt[BigInt(x) for x in discriminants]
  params = Vector{BigInt}[discs]
  if degree(galois_group) == 1
    query = "SELECT field_id, polynomial, discriminant FROM field WHERE discriminant = ANY(\$1)"
  else
    id = _find_group_id(galois_group)
    @assert id !== missing
    query = "SELECT field_id, polynomial, discriminant FROM field WHERE discriminant = ANY(\$1) AND degree = $(degree(galois_group)) AND (group_id IS NULL OR group_id = $id)"
  end
  result = execute(db, query, params, column_types = Dict(:polynomial => Vector{BigInt}, :discriminant => BigInt))
  data = Tables.rows(result)
  res = DBField[]
  Qx = PolynomialRing(FlintQQ, "x")[1]
  for r in data
    f = DBField(db, r[1])
    f.degree = length(r[2])-1
    f.polynomial = Qx(map(fmpz, r[2]))
    f.discriminant = fmpz(r[3])
    push!(res, f)
  end
  return res
end

function load_fields_with_discriminant(db::LibPQ.Connection, discriminants::Vector{fmpz}, deg::Int = -1)
  discs = BigInt[BigInt(x) for x in discriminants]
  params = Vector{BigInt}[discs]
  if deg == -1
    query = "SELECT field_id, polynomial, discriminant FROM field WHERE discriminant = ANY(\$1)"
  else
    query = "SELECT field_id, polynomial, discriminant FROM field WHERE discriminant = ANY(\$1) AND degree = $(deg)"
  end
  result = execute(db, query, params, column_types = Dict(:polynomial => Vector{BigInt}, :discriminant => BigInt))
  data = Tables.rows(result)
  res = DBField[]
  Qx = PolynomialRing(FlintQQ, "x")[1]
  for r in data
    f = DBField(db, r[1])
    f.degree = length(r[2])-1
    f.polynomial = Qx(map(fmpz, r[2]))
    f.discriminant = fmpz(r[3])
    push!(res, f)
  end
  return res
end


function load_fields(connection::LibPQ.Connection; degree::Int = -1, degree_range::UnitRange{Int} = -1:-1, discriminant_range::Tuple{fmpz, fmpz} = (fmpz(0), fmpz(-1)), 
                         signature::Tuple{Int, Int} = (-1, 0), unramified_outside::Vector{fmpz} = fmpz[], ramified_at::Vector{fmpz} = fmpz[],
                         galois_group::PermGroup = symmetric_group(1), class_number::Int = -1, 
                         class_group_structure::Vector{fmpz} = fmpz[-1], 
                         class_group_ranks_range::Dict{fmpz, Tuple{Int, Int}} = Dict{fmpz, Tuple{Int, Int}}(), only_count::Type{Val{T}} = Val{false}) where T


  parameters = String[]
  values = []
  ind = 1

  if degree != -1
    push!(parameters, "degree = $(degree)")
  end
  if first(degree_range) != -1
    if first(degree_range) == last(degree_range)
      push!(parameters, "degree = $(first(degree_range))")
    else
      push!(parameters, "degree BETWEEN $(first(degree_range)) AND $(last(degree_range)) ")
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
    push!(parameters, "real_embeddings = $(signature[1])")
  end
  if !isempty(unramified_outside)
    push!(parameters, "ramified_primes <@ \$$(ind)")
    push!(values, BigInt[BigInt(x) for x in unramified_outside])
    ind += 1
  end
  if !isempty(ramified_at)
    push!(parameters, "ramified_primes @> \$$(ind)")
    push!(values, BigInt[BigInt(x) for x in ramified_at])
    ind += 1
  end
  if !isone(order(galois_group))
    id_group = _find_group_id(connection, galois_group)
    if id_group === missing
      return DBField[]
    end
    push!(parameters, "group_id = $(id_group)")
  end
  if class_group_structure[1] != -1
    id_class_group = _find_class_group_id(connection, abelian_group(class_group_structure))
    if id_class_group !== missing
      push!(parameters, "class_group_id = $(id_class_group)")
    end
  end
  if class_number != -1
    abgroups = abelian_groups(class_number)
    class_group_ids = Vector{Int}()
    for C in abgroups
      idC = _find_class_group_id(connection, C)
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
    ids_class_group = _find_class_group_ids(connection, class_group_ranks_range)
    push!(parameters, "class_group_id = ANY( \$$(ind)) ")
    push!(values, ids_class_group)
    ind += 1
  end
  #Now, I can do the query
  if only_count == Val{true}
    query = "SELECT COUNT(*) FROM field"
  else
    query = "SELECT field_id FROM field"
  end
  if !isempty(parameters)
    query *= " WHERE "
    for i = 1:length(parameters)-1
      query = query * parameters[i] * " AND "
    end
    query =  query * parameters[end]
  end
  result = execute(connection, query, values)
  data1 = columntable(result)
  if only_count == Val{true}
    return data1[1][1]::Int
  end
  data = data1[1]
  l = length(data)::Int
  fields = Vector{DBField}(undef, l)
  for j = 1:length(fields)
    fields[j] = DBField(connection, data[j])
  end
  return fields
end

function _find_group_id(connection::LibPQ.Connection, G::PermGroup)
  d = degree(G)
  id = transitive_group_identification(G)
  if id != -1
    #GREAT! Unique identification.
    query = "SELECT group_id FROM galois_group WHERE degree = \$1 AND transitive_group_id = \$2"
    values = Int[d, id]
    result = execute(connection, query, values, column_types = Dict(:group_id => Int64))
    data = columntable(result)[1][1]
    return data::Union{Missing, Int}
  end
  o = order(G)
  if o < 2000 && o != 1024
    #Now, I try the small group id
    id = small_group_identification(G)
    query = "SELECT group_id FROM galois_group WHERE degree = \$1 AND group_order = \$2 AND small_group_id = \$3"
    values = [d, o, id[2]]
    result = execute(connection, query, values, column_types = Dict(:group_id => Int64))
    data = columntable(result)[1][1]
    return data::Union{Missing, Int}
  end
  #Sad. That's quite hard now. I have to search for all the groups
  #with a given order and check isomorphism.
  query = "SELECT group_id, generators FROM galois_group WHERE degree = \$1 AND group_order = \$2"
  result = execute(connection, query, [d, o], column_types = Dict(:group_id => Int64, :generators => String))
  data = Tables.rows(result)
  S = symmetric_group(d)
  for r in data
    vects = eval(Meta.parse(r[2]))
    perms = PermGroupElem[S(x) for x in vects]
    H, mH = sub(S, perms)
    if isisomorphic(G, H)[1]
      return r[1]::Int
    end
  end
  return missing
end

function _find_class_group_id(connection::LibPQ.Connection, C::GrpAbFinGen)
  if isone(order(C))
    query = "SELECT class_group_id FROM class_group WHERE group_order = 1"
    result = execute(connection, query, column_types = Dict(:class_group_id => Int64))
    return columntable(result)[1][1]::Union{Missing, Int}
  end
  Csnf = snf(C)[1]
  invs = Vector{BigInt}(undef, ngens(Csnf))
  for i = 1:ngens(Csnf)
    invs[i] = BigInt(Csnf.snf[i])
  end
  query = "SELECT class_group_id FROM class_group WHERE structure = \$1"
  result = execute(connection, query, [invs], column_types = Dict(:class_group_id => Int64))
  res = columntable(result)
  return res[1][1]::Union{Missing, Int}
end

function _find_class_group_ids(connection::LibPQ.Connection, ranks::Dict{fmpz, Tuple{Int, Int}})
  #I want to find the class group ids of the class group having the rank as required.
  divs = [BigInt(x) for x in keys(ranks)]
  sort!(divs)
  query = "SELECT class_group_id, prime_divisors, ranks FROM class_group WHERE prime_divisors @> \$1"
  result = execute(connection, query, [divs])
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

function find_group(connection::LibPQ.Connection, id::Int)
  #First, I try with the transitive_group_id
  query1 = "SELECT degree, transitive_group_id FROM galois_group WHERE group_id = \$1"
  result = execute(connection, query1, [id], column_types = Dict(:degree => Int64, :transitive_group_id => Int64))
  data = columntable(result)
  deg = data[1][1]
  if data[2][1] !== missing
    return transitive_group(deg, data[2][1])
  end
  query2 = "SELECT group_order, small_group_id FROM galois_group WHERE group_id = \$1"
  result2 = execute(connection, query2, [id], column_types = Dict(:group_order => Int64, :small_group_id => Int64))
  data2 = columntable(result2)
  if data2[2][1] !== missing
    #I have a polycyclic group and I want to get a perm group of the right degree.
    PC = small_group(data2[1][1], data2[2][1])
    H = isomorphic_transitive_perm_group(PC, deg)
    return H
  end
  query3 = "SELECT generators FROM galois_group WHERE group_id = \$1"
  result3 = execute(connection, query3, [id])
  data3 = columntable(result3)[1][1]
  S = symmetric_group(deg)
  g = eval(Meta.parse(data3))
  perms = Vector{PermGroupElem}(undef, length(g))
  for i = 1:length(g)
    perms[i] = S(g[i])
  end
  H, mH = sub(S, perms)
  return H
end

function find_DBfield(connection::LibPQ.Connection, K::AnticNumberField; already_in_DB::Bool = false)
  #First, I check whether the polynomial defining K is in the database.
  if !isdefining_polynomial_nice(K)
    K = simplify(K, cached = false, save_LLL_basis = false)[1]
  end
  pol = K.pol
  coeffs = BigInt[BigInt(numerator(coeff(pol, i))) for i = 0:degree(K)]
  sig = signature(K)
  d = degree(K)
  disc = discriminant(maximal_order(K))
  #=
  query = "SELECT field_id FROM field WHERE degree = $d AND discriminant = $disc AND real_embeddings = $(sig[1]) AND polynomial = \$1"
  @time result = columntable(execute(connection, query, Vector{BigInt}[coeffs]))
  if result[1][1] !== missing
    return DBField(connection, result[1][1])
  end
  #Bad luck. Now we need to check isomorphism
  =#
  #TODO: Think about it. Would it make sense to compute the canonical defining equation?
  lf = load_fields(connection, degree = d, discriminant_range = (disc, disc), signature = sig)
  if isempty(lf)
    return missing
  end
  if length(lf) == 1 && already_in_DB
    return lf[1]
  end
  ck = coefficients(defining_polynomial(K))
  for x in lf
    if coefficients(defining_polynomial(x)) == ck
      return x
    end
  end
  for x in lf
    fl = isisomorphic(K, number_field(x))[1]
    if fl
      return x
    end
  end
  return missing
end