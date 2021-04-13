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
    return "Degree $(d) SGid: $(id1)"
  catch e
    
  end
  return "Degree $(d) DBid: $(_find_group_id(db, G))"
end

function _construct_matrix(tb, G, db)
  #Now, I need to organize the table for the pretty printing.
  M = Array{String, 2}(undef, length(tb[1]), 5)
  M[1, 1] = _print_group(db, G)
  for i = 1:length(tb[1])
    if i > 1
      M[i, 1] = ""
    end
    M[i, 2] = string(tb[1][i])
    signatr = (tb[3][i], div(degree(G)-tb[3][i], 2))
    M[i, 3] = string(signatr)
    e = Int(round(log(10, tb[2][i])))
    M[i, 4] = string("~10^$e")
    M[i, 5] = string(load_fields(db, galois_group = G, signature = signatr, only_count = Val{true}))
  end
  return M
end

function _construct_matrix(db::LibPQ.Connection, G::PermGroup)
  id = _find_group_id(db, G)
  query = "SELECT GRH, discriminant_bound, real_embeddings FROM completeness WHERE group_id = \$1"
  result = execute(db, query, [id], column_types = Dict(:GRH => Bool, :discriminant_bound => BigInt, :real_embeddings => Int))
  tb = columntable(result)
  if tb[1][1] === missing
    return Array{String, 2}(undef, 0, 5)
  end
  return _construct_matrix(tb, G, db)
end

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
  pretty_table(t, ["Group", "GRH", "Signature", "Discriminant bound", "Number of fields"], compact_printing = true)
  return nothing
end

function completeness_data(db::LibPQ.Connection, degree::Int)
  #First, we find the groups.
  query = "SELECT group_id FROM galois_group WHERE degree = \$1"
  result = execute(db, query, [degree], column_types = Dict(:group_id => Int))
  tg = columntable(result)
  if tg[1][1] === missing
    println("Data not available in the database")
    return nothing
  end
  groups = PermGroup[]
  for x in tg[1]
    r = find_group(db, x) 
    if r !== missing
      push!(groups, r)
    end
  end
  M = Array{String, 2}(undef, 0, 5)
  for G in groups
    M = vcat(M, _construct_matrix(db, G))
  end
  t = Tables.table(M)
  pretty_table(t, ["Group", "GRH", "Signature", "Discriminant bound", "Number of fields"], compact_printing = true)
  return nothing
end

function find_completeness_data(connection::LibPQ.Connection, G::PermGroup, signature::Tuple{Int, Int})
  idG = _find_group_id(connection, G)
  if idG === missing
    return missing
  end
  query = "SELECT GRH, discriminant_bound FROM completeness WHERE group_id = \$1 AND real_embeddings = \$2"
  result = execute(connection, query, [idG, signature[1]])
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