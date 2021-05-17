
function fields_by_degree(db::LibPQ.Connection)
  query = "SELECT degree, real_embeddings, COUNT(field_id) "
  query *= "FROM field GROUP BY degree, real_embeddings ORDER BY degree, real_embeddings"
  result = Tables.columntable(execute(db, query, column_types = Dict(:degree => Int, :real_embeddings => Int)))
  M = Array{String, 2}(undef, length(result[1]), 3)
  for i = 1:length(result[1])
    if i == 1 || result[1][i-1] != result[1][i]
      M[i, 1] = string(result[1][i])
    else
      M[i, 1] = " "
    end
    s = result[2][i]
    M[i, 2] = "($s, $(div(result[1][i]-s, 2)))"
    M[i, 3] = string(result[3][i])
  end
  t = Tables.table(M)
  pretty_table(t, ["Degree", "Signature", "Number of fields"])
  return nothing
end

function class_numbers_by_degree(db::LibPQ.Connection)
  M = Array{String, 2}(undef, 48, 12)
  for i = 1:48
    for j = 1:12
      M[i, j] = ""
    end
  end
  query = "SELECT degree, COUNT(field_id) FROM field WHERE class_group_id = (SELECT class_group_id FROM class_group WHERE group_order = 1) GROUP BY degree"
  result = Tables.columntable(execute(db, query))
  for j = 1:length(result[1])
    M[result[1][j], 1] = string(result[2][j])
  end
  for i = 2:11
    query = "SELECT degree, COUNT(field_id) FROM field WHERE class_group_id = ANY(SELECT class_group_id FROM class_group WHERE group_order BETWEEN $(10^(i-2)+1) AND $(10^(i-1))) GROUP BY degree"
    result = Tables.columntable(execute(db, query))
    for j = 1:length(result[1])
      M[result[1][j], i] = string(result[2][j])
    end
  end
  t = Tables.table(M)
  header1 = Vector{String}(undef, 12)
  header2 = Vector{String}(undef, 12)
  header3 = Vector{String}(undef, 12)
  header1[1] = " "
  header2[1] = "h = 1"
  header3[1] = " "
  for i = 2:12
    header1[i] = "10^$(i-2)"
    header2[i] = "< h  <="
    header3[i] = "10^$(i-1)"
  end
  pretty_table(t, header = (header1, header2, header3), row_names = ["$i" for i in 1:48])
  return nothing
end

function _missing_class_groups_by_degree(db::LibPQ.Connection)
  query = "SELECT degree, COUNT(*) FROM field WHERE class_group_id IS NULL GROUP BY degree"
  result = Tables.columntable(execute(db, query))
  pretty_table(result)
  return nothing
end

function _missing_subfields_by_degree(db::LibPQ.Connection)
  query = "SELECT degree, COUNT(*) FROM field WHERE subfields IS NULL GROUP BY degree"
  result = Tables.columntable(execute(db, query))
  pretty_table(result)
  return nothing
end

function minimal_discriminant(db::LibPQ.Connection, GP::PermGroup, signature::Tuple{Int, Int})
  id = _find_group_id(db, GP)
  if id === missing
    error("Group not in database")
  end
  query1 = "SELECT MIN(ABS(discriminant)) FROM field where group_id = $id"
  min_disc = Tables.columtable(execute(db, query))[1][1]
  ps = possible_signatures(GP)
  b = fmpz(0)
  for s in ps
    cd = FieldsDB.find_completeness_data(db, GP, s)
    b = min(b, max(cd))
    if b < min_disc
      break
    end
  end
  print("The minimal discriminant is $min_disc")
  if b < min_disc
    println(" (not proven)")
  else
    println(" (proven)")
  end
end