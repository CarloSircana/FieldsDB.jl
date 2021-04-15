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
  end
  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()

  n = 1 
  i = 1
  n_batch = 1
  root_disc = 1

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
    end
  end
  
  grp_id = (n, i)
  discriminant_bound = fmpz(root_disc)^n

  
  Hecke.set_verbose_level(:Fields, 3)
  f_batch = open("./batch_$(n)_$(i)_$(n_batch).log", "r")
  conds = eval(Meta.parse(readline(f_batch)))
  close(f_batch)
  G = GAP.Globals.SmallGroup(n, i)
  li = GAP.gap_to_julia(Vector{Int}, GAP.Globals.AbelianInvariants(G))
  li = map(Int, snf(abelian_group(li))[1].snf)
  lf = abelian_fields(li, conds, discriminant_bound)
  flds = AnticNumberField[]
  for x in lf
    L = number_field(x)
    Lns = number_field(NfAbsNS, L)[1]
    Ls = Hecke.simplified_simple_extension(Lns)[1]
    push!(flds, Ls)
  end
  f_check = open("check_$(n)_$(i)_$(n_batch).log", "w")
  for x in flds
    println(f_check, coefficients(defining_polynomial(x)))
  end
  close(f_check)
  file = open("./password.log", "r")
  if !isfile(file)
    throw(error("Password not found!"))
  end
  s = readline(file)
  db = FieldsDB.LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag password =" * s)
  Goscar = small_group(n, i)
  GP = FieldsDB.isomorphic_transitive_perm_group(Goscar, n)
  FieldsDB.insert_fields(db, flds, galois_group = GP)
  close(db)
end

main()