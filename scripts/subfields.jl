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
           WHERE degree = $degree AND subfields IS NULL
           LIMIT $batch_size"
  @time result = Tables.rows(FieldsDB.LibPQ.execute(db, query, column_types = Dict(:polynomial => Vector{BigInt})))
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

function _get_subfields(flds::Vector{AnticNumberField}, res::Vector{FieldsDB.DBField})
  subs = Vector{AnticNumberField}[]
  deg = degree(flds[1])
  for i = 1:length(flds)
    K1 = flds[i]
    OK1 = lll(maximal_order(K1))
    K = simplify(K1, cached = false)[1]
    nbK = sum(nbits(numerator(x)) for x in coefficients(defining_polynomial(K)))
    nbK1 = sum(nbits(numerator(x)) for x in coefficients(defining_polynomial(K1)))
    if nbK < nbK1
      FieldsDB.set_polynomial(res[i], defining_polynomial(K))
    end
    println(K.pol)
    auts = automorphisms(K)
    if length(auts) == deg
      @time subf = Hecke.subfields_normal(K, true)
    else
      @time subf = subfields(K)
    end
    to_insert = AnticNumberField[]
    @time for (x, mx) in subf
      if degree(x) == 1
        push!(to_insert, Hecke.rationals_as_number_field()[1])
      elseif degree(x) == degree(K)
        push!(to_insert, K)
      else
        push!(to_insert, simplify(x, cached = false, save_LLL_basis = false)[1])
      end
    end
    push!(subs, FieldsDB.isomorphism_class_representatives(to_insert))
  end
  return subs
end

function main_loop(db::FieldsDB.LibPQ.Connection, deg::Int, batch_size::Int, res::Vector{FieldsDB.DBField})
  flds = AnticNumberField[number_field(x) for x in res]
  subs = _get_subfields(flds, res)
  degs = Set{Int}([degree(x) for y in subs for x in y])
  ids = Tuple{AnticNumberField, Int}[]
  fields_to_insert = AnticNumberField[]
  for d in degs
    if d == 1
      lf = load_fields(db, degree = 1)
      push!(ids, (Hecke.rationals_as_number_field()[1], lf[1].id))
      continue
    elseif d == deg
      for j = 1:length(flds)
          for i = 1:length(subs[j])
            if degree(subs[j][i]) == d
              push!(ids, (subs[j][i], res[j].id))
              break
            end
          end
      end
      continue
    end
    flds_d = AnticNumberField[]
    for x in subs
      for s = 1:length(x)
        y = x[s]
        if degree(y) == d
          found = false
          for j = 1:length(flds_d)
            if isisomorphic(flds_d[j], y)[1]
              found = true
              x[s] = flds_d[j]
              break
            end
          end
          if !found
            push!(flds_d, y)
          end
        end
      end
    end
    discs = Set(fmpz[discriminant(maximal_order(x)) for x in flds_d])
    println("Check if the fields of degree $d are already in the database")
    @time lf = FieldsDB.load_fields_with_discriminant(db, collect(discs), d)
    for x in flds_d
      polx = coefficients(defining_polynomial(x))
      found = false
      for y in lf
        if coefficients(defining_polynomial(y)) == polx
          found = true
          push!(ids, (x, y.id))
          break
        end
      end
      if found 
        continue
      end
      dx = discriminant(maximal_order(x))
      same_disc = FieldsDB.DBField[]
      for y in lf
        if discriminant(y) == dx
          push!(same_disc, y)
        end
      end
      if isempty(same_disc)
        push!(fields_to_insert, x)
        continue
      end
      for y in same_disc
        if isisomorphic(x, number_field(y))[1]
          found = true
          push!(ids, (x, y.id))
          break
        end
      end
      if !found
        push!(fields_to_insert, x)
      end
    end
  end
  println("Inserting fields")
  @time FieldsDB._insert_fields(fields_to_insert, db)
  pols_to_insert = Vector{BigInt}[]
  #COULD BE OPTIMIZED
  println("Retrieving id from database")
  @time for x in fields_to_insert
    push!(ids, (x, FieldsDB.find_DBfield(db, x, already_in_DB = true).id))
  end
  #Now, I have all the ids I need.
  subfs_ids = Vector{Vector{Int}}()
  for i = 1:length(subs)
    ids_i = Vector{Int}()
    for j = 1:length(subs[i])
      if degree(subs[i][j]) == 1
        ind = findfirst(x -> degree(x[1]) == 1, ids)
      else
        ind = findfirst(x -> x[1] == subs[i][j], ids)
      end
      push!(ids_i, ids[ind][2])
    end
    push!(subfs_ids, ids_i)
  end
  @time set_subfields(res, subfs_ids)
  return nothing
end

function set_subfields(flds::Vector{FieldsDB.DBField}, subfields::Vector{Vector{Int}})
  @assert length(flds) == length(subfields)
  values_string = ""
  for i = 1:length(flds)-1
    values_string *= "($(flds[i].id),  \$$(i)::bigint[]), "
  end
  values_string *= "($(flds[end].id), \$$(length(flds))::bigint[]) "
  query = " UPDATE field SET subfields = c.subfields FROM (VALUES " 
  query *= values_string 
  query *= " ) as c(field_id, subfields) WHERE field.field_id = c.field_id"
  FieldsDB.LibPQ.execute(flds[1].connection, query, subfields)
end

main()
