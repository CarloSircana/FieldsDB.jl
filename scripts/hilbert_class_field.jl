using FieldsDB, Printf, ArgParse

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--degree"
      help = "Degree of the base fields"
      arg_type = Int
      default = -1
    "--batch_size"
      help = "Number of fields per process"
      arg_type = Int
      default = 1
    "--degree_bound"
      help = "Bound for the degree of the Hilbert class field"
      arg_type = Int
      default = -1
  end
  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()

  deg = 1
  degree_bound = 0
  batch_size = 1

  for (arg, val) in parsed_args
    println("$arg => $val")
    if arg == "degree"
      deg = val
    elseif arg == "degree_bound"
      degree_bound = val
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
  class_group_ids = Int[]
  for i = 2:div(degree_bound, deg)
    for C in abelian_groups(i)
      id = FieldsDB._find_class_group_id(db, C)
      if id !== missing
        push!(class_group_ids, id)
      end
    end
  end
  println("Batch number $cnt")
  batch = get_batch(db, deg, batch_size, class_group_ids)
  while !isempty(batch)
    main_loop(db, batch_size, batch)
    cnt += 1
    println("Batch number $cnt")
    batch = get_batch(db, deg, batch_size, class_group_ids)
  end
  close(db)
end

function get_batch(db::FieldsDB.LibPQ.Connection, degree::Int, batch_size::Int, class_group_ids::Vector{Int})
  println("Retrieving fields")
  query = "SELECT field_id, polynomial
           FROM field 
           WHERE degree = $(degree) AND hilbert_class_field IS NULL AND class_group_id = ANY(\$1) AND random() < 0.2
           LIMIT $batch_size"
  @time result = FieldsDB.LibPQ.execute(db, query, [class_group_ids], column_types = Dict(:polynomial => Vector{BigInt}))
  Qx, x = PolynomialRing(FlintQQ, "x", cached = false)
  res = Vector{FieldsDB.DBField}(undef, length(result))
  for i = 1:length(result)
    res[i] = FieldsDB.DBField(db, result[i, 1])
    res[i].polynomial = Qx(map(fmpz, result[i, 2]))
  end
  return res
end

function main_loop(db::FieldsDB.LibPQ.Connection, batch_size::Int, res::Vector{FieldsDB.DBField})
  flds = AnticNumberField[number_field(x) for x in res]
  class_fields = AnticNumberField[absolute_simple_field(number_field(hilbert_class_field(x), using_norm_relation = true, over_subfield = true), simplify = true)[1] for x in flds]
  println("Inserting fields")
  @time FieldsDB.insert_fields(db, class_fields)
  ids = Int[FieldsDB.find_DBfield(db, x, already_in_DB = true).id for x in class_fields]
  println("Set Hilbert class fields")
  @time set_hcf(res, ids)
  return nothing
end

function set_hcf(flds::Vector{FieldsDB.DBField}, ids::Vector{Int})
  @assert length(flds) == length(ids)
  values_string = ""
  for i = 1:length(flds)-1
    values_string *= "($(flds[i].id),  $(ids[i])::int), "
  end
  values_string *= "($(flds[end].id), $(ids[end])::int) "
  query = " UPDATE field SET hilbert_class_field = c.hilbert_class_field FROM (VALUES " 
  query *= values_string 
  query *= " ) as c(field_id, hilbert_class_field) WHERE field.field_id = c.field_id"
  FieldsDB.LibPQ.execute(flds[1].connection, query)
end

main()