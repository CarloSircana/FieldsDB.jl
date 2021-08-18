export completeness_data

################################################################################
#
#  Completeness data
#
################################################################################

function _print_group(db, G::PermGroup)
  d = degree(G)
  try
    id = transitive_group_identification(G)
    if id != -1
      return "$(d)T$(id)"
    end
  catch e 

  end
  try 
    id1 = small_group_identification(G)
    return "SmallGrp$(id1)"
  catch e
    
  end
  return "Degree $(d) DBid: $(_find_group_id(db, G))"
end

function _sort_groups!(db, Gs::Vector{PermGroup})
  ltt = function(G)
    d = degree(G)
    try id = transitive_group_identification(G)
      if id != -1
        return d * 100000 + id * 10000
      end
    catch e 
      rethrow(e)
    end

    try 
      id1 = small_group_identification(G)
      return d * 10000000 + id1[2] * 1000000
    catch e
      rethrow(e)
    end
    return _find_group_id(db, G)
  end
  return sort!(Gs, by = ltt)
end

function _construct_matrix(tb, G, db)
  #Now, I need to organize the table for the pretty printing.
  rows_to_consider = Int[]
  dict = Dict{Tuple{Int, Int}, Vector{Tuple{Bool, fmpz}}}()
  for i = 1:length(tb)
    sign = (tb[i, 3], div(degree(G)-tb[i, 3], 2))
    fl = tb[i, 1]
    bound = fmpz(tb[i, 2])
    if haskey(dict, sign)
      if dict[sign][1][1]
        if dict[sign][1][2] <= bound
          dict[sign][1] = (fl, bound)
        else
          push!(dict[sign], (fl, bound))
        end
      else
        if dict[sign][1][2] <= bound
          push!(dict[sign], (fl, bound))
        end
      end
    else
      dict[sign] = Tuple{Bool, fmpz}[(fl, bound)]
    end
  end
  to_print = [(sign, GRH, bound) for (sign, y) in dict for (GRH, bound) in y]
  M = Array{String, 2}(undef, length(to_print), 5)
  M[1, 1] = _print_group(db, G)
  for i = 1:length(to_print)
    if i > 1
      M[i, 1] = ""
    end
    M[i, 2] = string(to_print[i][2])
    M[i, 3] = string(to_print[i][1])
    disc = to_print[i][3]
    e, n = ispower(disc)
    s = string(n)
    if e > 1
      s *= "^"*string(e)
    else
      f = factor(disc)
      s = string(f)
    end
    M[i, 4] = s
    M[i, 5] = string(count_fields(db, galois_group = G, signature = to_print[i][1]))
  end
  return M
end

function _construct_matrix(db::LibPQ.Connection, G::PermGroup)
  id = _find_group_id(db, G)
  query = "SELECT GRH, discriminant_bound, real_embeddings FROM completeness WHERE group_id = \$1"
  result = execute(db, query, [id], column_types = Dict(:GRH => Bool, :discriminant_bound => BigInt, :real_embeddings => Int))
  if isempty(result)
    return Array{String, 2}(undef, 0, 5)
  end
  return _construct_matrix(result, G, db)
end

@doc Markdown.doc"""
    completeness_data(db::Database, G::PermGroup)

Prints a table giving information on the completeness of the table of fields with Galois group G.
For each possible signature, the discriminant bound column shows the value for which it has been proven that
the table is complete. The GRH column consists of a boolean value, stating whether GRH was assumed or not when
computing the data (GRH = true means that GRH was assumed).
"""
function completeness_data(db::LibPQ.Connection, G::PermGroup)
  id = _find_group_id(db, G)
  if id === missing
    println("Data not available in the database")
    return nothing
  end
  M = _construct_matrix(db, G)
  if iszero(size(M)[1])
    println("Data not available in the database")
    return nothing
  end
  t = Tables.table(M)
  pretty_table(t, ["Group", "GRH", "Signature", "Discriminant bound", "Number of fields"], compact_printing = true, crop = :none)
  return nothing
end

@doc Markdown.doc"""
    completeness_data(db::Database, degree::Int)

Prints a table giving information on the completeness of the table of fields of the degree given as input.
For each possible Galois group and possible signature, the discriminant bound column shows the value for which it has been proven that
the table is complete. The GRH column consists of a boolean value, stating whether GRH was assumed or not when
computing the data (GRH = true means that GRH was assumed).
"""
function completeness_data(db::LibPQ.Connection, degree::Int)
  #First, we find the groups.
  query = "SELECT group_id FROM galois_group WHERE degree = \$1"
  result = execute(db, query, [degree], column_types = Dict(:group_id => Int))
  if isempty(result) 
    println("Data not available in the database")
    return nothing
  end
  groups = PermGroup[]
  for i in 1:length(result)
    r = find_group(db, result[i, 1]) 
    if r !== missing
      push!(groups, r)
    end
  end
  groups = _sort_groups!(db, groups)
  M = Array{String, 2}(undef, 0, 5)
  for G in groups
    M = vcat(M, _construct_matrix(db, G))
  end
  t = Tables.table(M)
  pretty_table(t, ["Group", "GRH", "Signature", "Discriminant bound", "Number of fields"], compact_printing = true, crop = :none)
  return nothing
end

function find_completeness_data(connection::LibPQ.Connection, G::PermGroup, signature::Tuple{Int, Int})
  idG = _find_group_id(connection, G)
  if idG === missing
    return missing
  end
  query = "SELECT GRH, discriminant_bound FROM completeness WHERE group_id = $(idG) AND real_embeddings = $(signature[1])"
  result = execute(connection, query, column_types = Dict(:GRH => Bool, :discriminant_bound => BigInt))
  if isempty(result)
    return missing
  end
  #First result under GRH, second without.
  res = Vector{fmpz}(undef, 2)
  for i = 1:length(result)
    if result[i, 1]
      res[1] = fmpz(result[i, 2])
    else
      res[2] = fmpz(result[i, 2])
    end
  end
  if !isassigned(res, 1)
    res[1] = fmpz()
  end
  if !isassigned(res, 2)
    res[2] = fmpz()
  end
  return res
end


function insert_completeness_data(connection::LibPQ.Connection, group::PermGroup, signature::Tuple{Int, Int}, discriminant_bound::fmpz, GRH::Bool)
  cd = find_completeness_data(connection, group, signature)
  if cd !== missing
    @show "here"
    if GRH
      if cd[1] > discriminant_bound
        return nothing
      elseif !iszero(cd[1]) 
        gid = _find_group_id(connection, group)
        query = "UPDATE completeness SET discriminant_bound = \$1 WHERE GRH = \$2 AND group_id = \$3 AND real_embeddings = \$4"
        execute(connection, query, [BigInt(discriminant_bound), GRH, gid, signature[1]])
        return nothing
      end
    else
      if cd[2] > discriminant_bound
        return nothing
      elseif !iszero(cd[2])
        gid = _find_group_id(connection, group)
        query = "UPDATE completeness SET discriminant_bound = \$1 WHERE GRH = \$2 AND group_id = \$3 AND real_embeddings = \$4"
        execute(connection, query, [BigInt(discriminant_bound), GRH, gid, signature[1]])
        return nothing
      end
    end
  end
  gid = _find_group_id(connection, group)
  @assert gid !== missing
  LibPQ.load!(
    (group_id = [gid], 
    GRH = [GRH],
    real_embeddings = [signature[1]],
    discriminant_bound = [BigInt(discriminant_bound)]
    ),
    connection,
    "INSERT INTO completeness (
      group_id, 
      GRH,
      real_embeddings,
      discriminant_bound
    ) VALUES (\$1, \$2, \$3, \$4);"
  )
  return nothing
end

function _remove_completeness_data(db::LibPQ.Connection, G::PermGroup, signature::Tuple{Int, Int})
  gid = _find_group_id(db, G)
  @assert gid !== missing
  real_embeddings = signature[1]
  query = "DELETE FROM completeness WHERE group_id = $(gid)::int AND real_embeddings = $(real_embeddings)::int"
  execute(db, query)
end
