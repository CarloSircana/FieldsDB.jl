using FieldsDB, Printf, ArgParse

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--degree"
      help = "Degree of the fields"
      arg_type = Int
      default = -1
    "--batch_size"
      help = "Number of fields per process"
      arg_type = Int
      default = 1
    "--id"
      help = "Transitive group id"
      arg_type = Int
      default = -1
  end
  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()

  deg = 1
  transitive_id = 0
  batch_size = 1

  for (arg, val) in parsed_args
    println("$arg => $val")
    if arg == "degree"
      deg = val
    elseif arg == "id"
      transitive_id = val
    elseif arg == "batch_size"
      batch_size = val
    end
  end

  file = open("./password.log", "r")
  if !isfile(file)
    throw(error("Password not found!"))
  end
  s = readline(file)
  db = fields_database(s)
  cnt = 1
  println("Batch number $cnt")
  GP = transitive_group(deg, transitive_id)
  grp_id = FieldsDB._find_group_id(db, GP)
  target_group = FieldsDB.isomorphic_transitive_perm_group(GP, Int(order(GP)))
  batch = get_batch(db, grp_id, batch_size)
  while !isempty(batch)
    main_loop(db, target_group, batch_size, batch)
    cnt += 1
    println("Batch number $cnt")
    batch = get_batch(db, grp_id, batch_size)
  end
  close(db)
end

function get_batch(db::FieldsDB.LibPQ.Connection, grp_id::Int, batch_size::Int)
  println("Retrieving fields")
  query = "SELECT field_id, polynomial
           FROM field 
           WHERE group_id = $grp_id AND normal_closure IS NULL AND random() < 0.2
           LIMIT $batch_size"
  @time result = FieldsDB.LibPQ.execute(db, query, column_types = Dict(:polynomial => Vector{BigInt}))
  Qx, x = PolynomialRing(FlintQQ, "x", cached = false)
  res = Vector{FieldsDB.DBField}(undef, length(result))
  for i = 1:length(result)
    res[i] = FieldsDB.DBField(db, result[i, 1])
    res[i].polynomial = Qx(map(fmpz, result[i, 2]))
  end
  return res
end

function main_loop(db::FieldsDB.LibPQ.Connection, target_group::PermGroup, batch_size::Int, res::Vector{FieldsDB.DBField})
  flds = AnticNumberField[number_field(x) for x in res]
  closures = AnticNumberField[splitting_field(defining_polynomial(x)) for x in flds]
  for i = 1:length(closures)
    @assert degree(closures[i]) == degree(target_group)
  end
  println("Inserting fields")
  @time FieldsDB.insert_fields(db, closures, galois_group = target_group)
  #TODO: Check if I can get the ids immediately
  ids = Int[FieldsDB.find_DBfield(db, x, already_in_DB = true).id for x in closures]
  println("Set closures")
  @time set_closures(res, ids)
  return nothing
end

function set_closures(flds::Vector{FieldsDB.DBField}, ids::Vector{Int})
  @assert length(flds) == length(ids)
  values_string = ""
  for i = 1:length(flds)-1
    values_string *= "($(flds[i].id),  $(ids[i])::int), "
  end
  values_string *= "($(flds[end].id), $(ids[end])::int) "
  query = " UPDATE field SET normal_closure = c.normal_closure FROM (VALUES " 
  query *= values_string 
  query *= " ) as c(field_id, normal_closure) WHERE field.field_id = c.field_id"
  FieldsDB.LibPQ.execute(flds[1].connection, query)
end

main()