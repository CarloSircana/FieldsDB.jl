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
    "--t"
      help = "Number of processes"
      arg_type = Int
      default = 1
    "--batch_size"
      help = "Number of fields per process"
      arg_type = Int
      default = 1
    "--only_real"
      help = "Only real fields flag"
      action = :store_true
    "--rt"
      help = "Root discriminant"
      arg_type = Int
      default = 1
  end
  return parse_args(s)
end

function main()
  @show parsed_args = parse_commandline()

  n = 1 
  i = 1
  n_proc = 1
  root_disc = 1
  batch_size = 1
  only_real = false

  for (arg, val) in parsed_args
    println("$arg => $val")
    @show arg
    if arg == "n"
      n = val
    elseif arg == "id"
      i = val
    elseif arg == "t"
      n_proc = val
    elseif arg == "rt"
      root_disc = val
    elseif arg == "batch_size"
      batch_size = val
    elseif arg == "only_real"
      only_real = val
    end
  end

  G = small_group(n, i)
  julia_exe = Base.julia_cmd()
  if isabelian(G)
    fields_abelian_control(n, i, root_disc, batch_size, n_proc, only_real)
  else
    fields_nonabelian_control(n, i, root_disc, batch_size, n_proc, only_real)
  end
  return nothing
end

function fields_nonabelian_control(n::Int, i::Int, root_disc::Int, batch_size::Int, n_proc::Int, only_real::Bool)
  discriminant_bound = fmpz(root_disc)^n
  G = small_group(n, i)
  L = derived_series(G)
  G1 = quo(G, L[end-1])[1]
  discriminant_bound_subfield = fmpz(root_disc)^order(G1)
  id_G1 = small_group_identification(G1)
  oG1 = Int(order(G1))
  G1P = FieldsDB.isomorphic_transitive_perm_group(G1, oG1)
  file = open("./password.log", "r")
  if !isfile(file)
    throw(error("Password not found!"))
  end
  s = readline(file)
  db = FieldsDB.LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag password =" * s)
  only_real_base_field = only_real || isodd(oG1) || (Hecke._real_level(GAP.Globals.DerivedSeries(G.X)) == length(L)-1)
  if only_real_base_field
    cd = FieldsDB.find_completeness_data(db, G1P, (oG1, 0))
    if maximum(cd) < discriminant_bound_subfield
      error("The database does not contain enough fields!")
    end
  else
    cd1 = FieldsDB.find_completeness_data(db, G1P, (oG1, 0))
    cd2 = FieldsDB.find_completeness_data(db, G1P, (0, divexact(oG1, 2)))
    if maximum(cd1) < discriminant_bound_subfield || maximum(cd2) < discriminant_bound_subfield
      error("The database does not contain enough fields!")
    end 
  end
  if only_real_base_field
    lf = load_fields(db, galois_group = G1P, discriminant_range = (-discriminant_bound_subfield, discriminant_bound_subfield), signature = (oG1, 0))
  else
    lf = load_fields(db, galois_group = G1P, discriminant_range = (-discriminant_bound_subfield, discriminant_bound_subfield))
  end
  ids = Int[x.id for x in lf]
  @show "Number of fields: $(length(ids))"
  #Now, we create the batch of fields
  julia_exe = Base.julia_cmd()
  errored_fields = Int[]
  total_number = div(length(ids), batch_size) +1
  procs = Cmd[]
  path_to_file = joinpath(@__DIR__, "fields_parallel_process.jl")
  for s = 1:total_number
    idsx_start = (s-1)*batch_size+1
    idsx_end = min(length(ids), s*batch_size)
    idsx = ids[idsx_start:idsx_end]
    f = open("./batch_$(n)_$(i)_$(s).log", "w")
    print(f, idsx)
    close(f)
    if only_real
      push!(procs, `$(julia_exe) $(path_to_file) --n=$n --id=$i --batch=$s --rt=$root_disc --only_real`)
    else
      push!(procs, `$(julia_exe) $(path_to_file) --n=$n --id=$i --batch=$s --rt=$root_disc`)
    end
  end
  ind = 1
  started_procs = []
  while ind <= length(procs) || count(process_exited, started_procs) != length(procs)
    @show number_running_procs = count(process_running, started_procs)
    if ind <= length(procs) && number_running_procs < n_proc
      push!(started_procs, run(procs[ind], wait = false))
      ind += 1
    else
      sleep(30)
    end
  end
  if all(success, started_procs)
    GP = FieldsDB.isomorphic_transitive_perm_group(G, n)
    FieldsDB.insert_completeness_data(db, GP, (n, 0), discriminant_bound, true)
    if iseven(n) && !only_real 
      FieldsDB.insert_completeness_data(db, GP, (0, div(n, 2)), discriminant_bound, true)
    end
  else
    f_err = open("./errored_$(n)_$(i).log", "w")
    for s = 1:length(procs)
      if !success(started_procs[s])
        println(f_err, s)
      else
        rm("./batch_$(n)_$(i).log")
        rm("./check_$(n)_$(i).log")
      end
    end
    close(f_err)
  end
  close(db)

end

function fields_abelian_control(n::Int, i::Int, root_disc::Int, batch_size::Int, n_proc::Int, only_real::Bool)
  G = GAP.Globals.SmallGroup(n, i)
  li = GAP.gap_to_julia(Vector{Int}, GAP.Globals.AbelianInvariants(G))
  li = map(Int, snf(abelian_group(li))[1].snf)
  discriminant_bound = fmpz(root_disc)^n
  KQ = Hecke.rationals_as_number_field()[1]
  OKQ = maximal_order(KQ)
  conds = Hecke.conductorsQQ(OKQ, li, discriminant_bound)
  @show "Number of conductors: $(length(conds))"
  batch = open("./batch", "w")
  julia_exe = Base.julia_cmd()
  errored_fields = Int[]
  total_number = div(length(conds), batch_size)+1
  procs = Cmd[]
  path_to_file = joinpath(@__DIR__, "fields_abelian_parallel_process.jl")
  for s = 1:total_number
    idsx_start = (s-1)*batch_size+1
    idsx_end = min(length(conds), s*batch_size)
    idsx = conds[idsx_start:idsx_end]
    f = open("./batch_$(n)_$(i)_$(s).log", "w")
    print(f, idsx)
    close(f)
    if only_real
      push!(procs, `$(julia_exe) $(path_to_file) --n=$n --id=$i --batch=$s --rt=$root_disc --only_real`)
    else
      push!(procs, `$(julia_exe) $(path_to_file) --n=$n --id=$i --batch=$s --rt=$root_disc`)
    end
  end
  ind = 1
  started_procs = []
  while ind <= length(procs) || count(process_exited, started_procs) != length(procs)
    @show number_running_procs = count(process_running, started_procs)
    if ind <= length(procs) && number_running_procs < n_proc
      push!(started_procs, run(procs[ind], wait = false))
      ind += 1
    else
      sleep(30)
    end
  end
  if all(success, started_procs)
    file = open("./password.log", "r")
    if !isfile(file)
      throw(error("Password not found!"))
    end
    s = readline(file)
    db = FieldsDB.LibPQ.Connection("host=tabularix dbname=fields port=5432 user=agag password =" * s)
    Goscar = small_group(n, i)
    GP = FieldsDB.isomorphic_transitive_perm_group(Goscar, n)
    FieldsDB.insert_completeness_data(db, GP, (n, 0), discriminant_bound, false)
    if !only_real && iseven(n)
      FieldsDB.insert_completeness_data(db, GP, (0, div(n, 2)), discriminant_bound, false)
    end
    close(db)
  else
    f_err = open("./errored_$(n)_$(i).log", "w")
    for s = 1:length(procs)
      if !success(started_procs[s])
        println(f_err, s)
      end
    end
    close(f_err)
  end

end

main()
