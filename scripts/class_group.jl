using FieldsDB, Decimals, Tables, Printf, ArgParse

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
    "--simplify"
      help = "Simplify the fields or not"
      action = :store_true
  end
  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()

  deg = 1
  batch_size = 1
  simplify = false

  for (arg, val) in parsed_args
    println("$arg => $val")
    if arg == "degree"
      deg = val
    elseif arg == "batch_size"
      batch_size = val
    elseif arg == "simplify"
      simplify = true
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
    main_loop(db, deg, batch_size, batch, simplify)
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
           WHERE degree = $degree AND class_group_id IS NULL AND random() < 0.1
           LIMIT $batch_size"
  @time result = Tables.rows(FieldsDB.LibPQ.execute(db, query, column_types = Dict(:polynomial => Vector{BigInt})))
  if length(result) < batch_size
    query = "SELECT field_id, polynomial
           FROM field 
           WHERE degree = $degree AND class_group_id IS NULL
           LIMIT $batch_size"
    @time result = Tables.rows(FieldsDB.LibPQ.execute(db, query, column_types = Dict(:polynomial => Vector{BigInt})))
  end
  Qx, x = PolynomialRing(FlintQQ, "x", cached = false)
  res = Vector{FieldsDB.DBField}(undef, batch_size)
  ind = 1
  for x in result
    res[ind] = FieldsDB.DBField(db, x[1])
    res[ind].polynomial = Qx(map(fmpz, x[2]))
    ind += 1
  end
  return res[1:ind-1]
end


function main_loop(db::FieldsDB.LibPQ.Connection, deg::Int, batch_size::Int, res::Vector{FieldsDB.DBField}, simplify::Bool)
  flds = AnticNumberField[number_field(x) for x in res]
  clgps = GrpAbFinGen[]
  for i in 1:length(flds)
    K = flds[i]
    println(K.pol)
    automorphisms(K)
    if simplify
      lll(maximal_order(K))
      K1 = Hecke.simplify(K, cached = false)[1]
      nbK = sum(nbits(numerator(x)) for x in coefficients(defining_polynomial(K)))
      nbK1 = sum(nbits(numerator(x)) for x in coefficients(defining_polynomial(K1)))
      if nbK1 < nbK
        println("Found a better defining polynomial")
        println(defining_polynomial(K1))
        FieldsDB.set_polynomial(res[i], defining_polynomial(K1))
        K = K1
      end
    end

    c = nothing

    if deg > 30
      OK = lll(maximal_order(K))
      fl, N = Hecke.norm_relation(K, small_degree = false)
      if fl
        @time c, = Hecke.NormRel.class_group_via_brauer(OK, N)
      else
        @time c, = class_group(K)
      end
    else
      @time c, = class_group(K)
    end
    push!(clgps, c)
  end

  @assert length(flds) == length(clgps)

  #First, I set the regulators.
  println("Setting regulators")
  regs = Decimal[]
  for x in flds
    if deg == 2 && signature(x) == (0, 1)
      r = Decimal(1)
    else
      r = FieldsDB._regulator_as_decimal(x)
    end
    push!(regs, r)
  end
  @time set_regulators(res, regs)
  #TODO: Put more thought in it.
  println("Setting class groups")
  @time for i = 1:length(res)
    FieldsDB.set_class_group(res[i], clgps[i])
  end
  return nothing
end

function set_regulators(flds::Vector{FieldsDB.DBField}, regs::Vector{Decimal})
  @assert length(flds) == length(regs)
  values_string = ""
  for i = 1:length(flds)-1
    values_string *= "($(flds[i].id), $(regs[i])), "
  end
  values_string *= "($(flds[end].id), $(regs[end]))"
  query = " UPDATE field
            SET regulator = c.regulator
            FROM (VALUES " * values_string * 
            "            
            ) as c(field_id, regulator)
          WHERE field.field_id = c.field_id
          "
  FieldsDB.LibPQ.execute(flds[1].connection, query)
end

main()
