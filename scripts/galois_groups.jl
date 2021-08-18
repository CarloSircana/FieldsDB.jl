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
  end
  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()

  deg = 1
  batch_size = 1

  for (arg, val) in parsed_args
    println("$arg => $val")
    if arg == "degree"
      deg = val
    elseif arg == "batch_size"
      batch_size = val
    end
  end

  file = open("./password.log", "r")
  if !isfile(file)
    throw(error("Password not found!"))
  end
  s = readline(file)
  db = FieldsDB.LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag password =" * s)
  cnt = 1
  println("Batch number $cnt")
  batch = get_batch(db, deg, batch_size)
  while !isempty(batch)
    main_loop(db, deg, batch_size, batch)
    cnt += 1
    println("Batch number $cnt")
    batch = get_batch(db, deg, batch_size)
  end
  close(db)
end

function get_batch(db::FieldsDB.LibPQ.Connection, degree::Int, batch_size::Int)
  println("Retrieving fields")
  query = "SELECT field_id, polynomial
           FROM field 
           WHERE degree = $degree AND group_id IS NULL AND random() < 0.3
           LIMIT $batch_size"
  @time result = FieldsDB.LibPQ.execute(db, query, column_types = Dict(:polynomial => Vector{BigInt}))
  Qx, x = PolynomialRing(FlintQQ, "x", cached = false)
  res = Vector{FieldsDB.DBField}(undef, length(result))
  for i = 0:length(res)
    res[i] = FieldsDB.DBField(db, result[i, 1])
    res[i].polynomial = Qx(map(fmpz, result[i, 2]))
  end
  return res[1:ind-1]
end


function main_loop(db::FieldsDB.LibPQ.Connection, deg::Int, batch_size::Int, res::Vector{FieldsDB.DBField})
  flds = AnticNumberField[number_field(x) for x in res]
  for (i, K) in enumerate(flds)
    println(K.pol)
    automorphisms(K)
    @time FieldsDB.set_galois_group(res[i])
  end
  return nothing
end

main()