################################################################################
#
#  Insertion
#
################################################################################

function insert_complete_table(db::LibPQ.Connection, fields::Vector{AnticNumberField}, galois_group::PermGroup, discriminant_bound::fmpz, GRH::Bool = true, signature::Tuple{Int, Int} = (-1, 0))
  g_id = _find_group_id(db, galois_group)
  if g_id === missing
    insert_group(db, galois_group)
    g_id = _find_group_id(db, galois_group)::Int
  end
  insert_fields(db, fields, galois_group = galois_group)
  #Finally, we save the completeness data
  have_signature = signature[1] != -1
  if have_signature
    insert_completeness_data(connection, galois_group, signature, discriminant_bound, GRH)
  else
    for sg in possible_signatures(galois_group)
      insert_completeness_data(connection, galois_group, sg, discriminant_bound, GRH)
    end
  end
  return nothing
end


function insert_fields(db::LibPQ.Connection, fields::Vector{AnticNumberField}; galois_group::PermGroup = symmetric_group(1))
  l = div(length(fields), 1000)+1
  for i = 1:l
    println("Inserting batch $i / $l")
    istart = (i-1)*1000+1
    iend = min(i*1000, length(fields))
    fieldsi = fields[istart:iend]
    @time insert_fields_split(db, fieldsi, galois_group)
  end
  return nothing
end

function insert_fields_split(db::LibPQ.Connection, fields::Vector{AnticNumberField}, galois_group::PermGroup = symmetric_group(1))
  if !isone(degree(galois_group))
    flds = _sieve_fields(db, fields, degree(galois_group))
  else
    flds = _sieve_fields(db, fields)
  end
  return _insert_fields(flds, db, galois_group = galois_group)
end

function _sieve_fields(db::LibPQ.Connection, fields::Vector{AnticNumberField}, deg::Int = -1)
  discs = collect(Set(fmpz[discriminant(maximal_order(x)) for x in fields]))
  if deg == -1
    lf = load_fields_with_discriminant(db, discs)
  else
    lf = load_fields_with_discriminant(db, discs, deg)
  end
  #First, we sieve the fields so that we remove the duplicates
  fields_to_insert = AnticNumberField[]
  for x in fields
    dx = degree(x)
    polx = defining_polynomial(x)
    found = false
    for y in lf
      if degree(y) == degree(x) && defining_polynomial(y) == polx
        found = true
        break
      end
    end
    if found 
      continue
    end
    dx = discriminant(maximal_order(x))
    same_disc = DBField[]
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
        break
      end
    end
    if !found
      push!(fields_to_insert, x)
    end
  end
  return fields_to_insert
end

function _insert_fields(fields::Vector{AnticNumberField}, connection::LibPQ.Connection; galois_group = symmetric_group(1))
  if order(galois_group) > 1
    g_id = _find_group_id(connection, galois_group)
    if g_id === missing
      insert_group(connection, galois_group)
      g_id = _find_group_id(connection, galois_group)::Int
    end
    aut_order = Int(find_automorphisms_order(galois_group))
  else
    g_id = 0
  end
  #First, I need to sieve the fields. I don't want to have duplicates in the database.
  
  
  execute(connection, "BEGIN;")
  for K1 in fields
    if isdefining_polynomial_nice(K1)
      K = K1
    else
      K = simplify(K1)[1]
      @assert isdefining_polynomial_nice(K)
    end
    d = discriminant(maximal_order(K))
    real_embs = Hecke.signature(K)[1]
    pol = BigInt[BigInt(numerator(coeff(K.pol, i))) for i = 0:degree(K)]
    deg = degree(K)
    if !iszero(g_id)
      LibPQ.load!(
        (real_embeddings = [real_embs], 
        polynomial = [pol], 
        discriminant = [BigInt(d)], 
        degree = [deg],
        group_id = [g_id],
        automorphisms_order = [aut_order]
        ),
        connection,
        "INSERT INTO field (
          real_embeddings, 
          polynomial,
          discriminant, 
          degree,
          group_id,
          automorphisms_order
        ) VALUES (\$1, \$2, \$3, \$4, \$5);",
      )
    else
      LibPQ.load!(
        (real_embeddings = [real_embs], 
        polynomial = [pol], 
        discriminant = [BigInt(d)], 
        degree = [deg]
        ),
        connection,
        "INSERT INTO field (
          real_embeddings, 
          polynomial,
          discriminant, 
          degree
        ) VALUES (\$1, \$2, \$3, \$4);",
      )
    end
  end
  execute(connection, "COMMIT;")
  return nothing
end

function insert_field(connection::LibPQ.Connection, K::AnticNumberField)
  return insert_fields(connection, AnticNumberField[K])
end

function insert_class_group(connection::LibPQ.Connection, C::GrpAbFinGen)
  o = BigInt(order(C))
  str = map(BigInt, snf(C)[1].snf)
  lf = factor(str[end])
  divs = BigInt[BigInt(x) for x in keys(lf.fac)]
  sort!(divs)
  ranks = Vector{Int}(undef, length(divs))
  for i = 1:length(divs)
    ind = 0
    for j = 1:length(str)
      if iszero(mod(str[j], divs[i]))
        ind = j
        break
      end
    end
    ranks[i] = length(str)-ind+1
  end
  LibPQ.load!(
    (group_order = [o], 
    structure = [str],
    prime_divisors = [divs],
    ranks = [ranks], 
    ),
    connection,
    "INSERT INTO class_group (
      group_order, 
      structure,
      prime_divisors,
      ranks
    ) VALUES (\$1, \$2, \$3, \$4);"
  )
  return nothing
end

function insert_group(connection::LibPQ.Connection, G::PermGroup)
  o = order(G)
  isab = isabelian(G)
  issolv = issolvable(G)
  issimp = issimple(G)
  isnilp = isnilpotent(G) 
  isperf = isperfect(G)
  isprim = isprimitive(G)
  d = degree(G)
  id = transitive_group_identification(G)
  if id != -1
    LibPQ.load!(
      (group_order = [o], 
      degree = [d],
      transitive_group_id = [id],
      abelian = [isab], 
      nilpotent = [isnilp], 
      solvable = [issolv],
      issimple = [issimp],
      perfect = [isperf],
      primitive = [isprim]
      ),
      connection,
      "INSERT INTO galois_group (
        group_order, 
        degree,
        transitive_group_id,
        abelian, 
        nilpotent, 
        solvable,
        issimple,
        perfect,
        primitive
      ) VALUES (\$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9);"
    )
    return nothing
  end
  try 
    id = small_group_identification(G)
    LibPQ.load!(
      (group_order = [o], 
      degree = [d],
      small_group_id = [id[2]],
      abelian = [isab], 
      nilpotent = [isnilp], 
      solvable = [issolv],
      issimple = [issimp],
      perfect = [isperf],
      primitive = [isprim]
      ),
      connection,
      "INSERT INTO galois_group (
        group_order, 
        degree,
        small_group_id,
        abelian, 
        nilpotent, 
        solvable,
        issimple,
        perfect,
        primitive
      ) VALUES (\$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9);"
    )
    return nothing
  catch e 

  end
  #Bad case. Cache a good set of generators.
  g = gens(G)
  s = "["
  for i = 1:length(g)-1
    vi = Int[g[i](j) for j = 1:degree(G)]
    s = s * "$vi , "
  end
  vi = Int[g[length(g)](j) for j = 1:degree(G)]
  s = s * "$vi]"
  LibPQ.load!(
    (group_order = [o], 
    degree = [d],
    generators = [s],
    abelian = [isab], 
    nilpotent = [isnilp], 
    solvable = [issolv],
    issimple = [issimp],
    perfect = [isperf],
    primitive = [isprim]
    ),
    connection,
    "INSERT INTO galois_group (
      group_order, 
      degree,
      generators,
      abelian, 
      nilpotent, 
      solvable,
      issimple,
      perfect, 
      primitive
    ) VALUES (\$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9);"
  )
  return nothing
end
