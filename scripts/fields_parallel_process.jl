using FieldsDB, Printf, ArgParse

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--n", "-o"
      help = "Order of the group"
      arg_type = Int
      default = -1
    "--id", "-i"
      help = "Id of the group"
      arg_type = Int
      default = -1
    "--batch"
      help = "Batch number"
      arg_type = Int
      default = 1
    "--rt"
      help = "Root discriminant"
      arg_type = Int
      default = 1
    "--only_real"
      help = "Only totally real fields flag"
      action = :store_true
  end
  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()

  n = 1 
  i = 1
  n_batch = 1
  root_disc = 1
  only_real = false

  for (arg, val) in parsed_args
    println("$arg => $val")
    @show arg
    if arg == "n"
      n = val
    elseif arg == "id"
      i = val
    elseif arg == "batch"
      n_batch = val
    elseif arg == "rt"
      root_disc = val
    elseif arg == "only_real"
      only_real = val
    end
  end
  
  grp_id = (n, i)
  discriminant_bound = fmpz(root_disc)^n

  file = open("./password.log", "r")
  if !isfile(file)
    throw(error("Password not found!"))
  end
  s = readline(file)
  db = FieldsDB.LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag password =" * s)
  Hecke.set_verbose_level(:Fields, 3)
  f_batch = open("./batch_$(n)_$(i)_$(n_batch).log", "r")
  ids = eval(Meta.parse(readline(f_batch)))
  close(f_batch)
  flds = FieldsDB.DBField[FieldsDB.DBField(db, s) for s in ids]
  ctxs = Hecke.FieldsTower[Hecke.field_context(number_field(x)) for x in flds]
  close(db)
  l = Hecke.fields(n, i, ctxs, discriminant_bound, only_real = only_real)
  flds_to_insert = AnticNumberField[number_field(x) for x in l]
  f = open("./check_$(n)_$(i)_$(n_batch).log", "w")
  println(f, length(flds_to_insert))
  for K in flds_to_insert
    println(f, collect(coefficients(K.pol)))
  end
  close(f)
  G = small_group(n, i)
  GP = FieldsDB.isomorphic_transitive_perm_group(G, n)
  db = FieldsDB.LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag password =" * s)
  FieldsDB.insert_fields(db, flds_to_insert, galois_group = GP)
  close(db)
end

main()