export fields

function field_context(x::DBField)
  return Hecke.field_context(number_field(x))
end

function fields(a::Int, b::Int, db::LibPQ.Connection, absolute_bound::fmpz; only_real::Bool = false, unramified_outside::Vector{fmpz} = fmpz[])
  G = GAP.Globals.SmallGroup(a, b)
  GP = isomorphic_transitive_perm_group(PcGroup(G), a)
  ps = Tuple{Int, Int}[(a, 0)]
  if !only_real && iseven(a)
    push!(ps, (0, divexact(a, 2)))
  end
  is_complete = true
  bound_db = fmpz(10)^10000
  for sgn in ps
    dbound_completeness = find_completeness_data(db, GP, sgn)
    if dbound_completeness === missing
      bound_db = fmpz()
      is_complete = false
      break
    end
    bound_db = min(bound_db, dbound_completeness[1])
    if bound_db < absolute_bound
      is_complete = false
    end
  end
  if is_complete
    if only_real
      lf = load_fields(db, discriminant_range = (-absolute_bound, absolute_bound), galois_group = GP, signature = (a, 0), unramified_outside = unramified_outside)
    else
      lf = load_fields(db, discriminant_range = (-absolute_bound, absolute_bound), galois_group = GP, unramified_outside = unramified_outside)
    end
    return map(field_context, lf)
  end
  L = GAP.Globals.DerivedSeries(G)
  invariants = GAP.gap_to_julia(Vector{Int}, GAP.Globals.AbelianInvariants(L[end-1]))
  lG = snf(abelian_group(invariants))[1]
  invariants = map(Int, lG.snf)
  if GAP.Globals.IsAbelian(G)
    @vprint :Fields 1 "computing abelian extension of Q with invariants $(invariants) and bound ~10^$(clog(absolute_bound, 10))\n"
    @vprint :FieldsNonFancy 1 "Doing Group ($a, $b) with bound $absolute_bound\n"
    fcs = Hecke.abelian_extensionsQQ(invariants, absolute_bound, only_real, unramified_outside = unramified_outside)
    #Now, we insert the fields in the database.
    if isempty(unramified_outside)
      if only_real
        insert_complete_table(db, map(number_field, fcs), GP, absolute_bound, true, (a, 0), check = true)
      else
        insert_complete_table(db, map(number_field, fcs), GP, absolute_bound, true, check = true)
      end
    end
    return fcs
  end
  must_be_ram_surely, must_be_ram_maybe = Oscar.Hecke.must_be_ramified(L, length(L)-1)
  lvl = Hecke._real_level(L)
  G1 = GAP.Globals.FactorGroup(L[1], L[end-1])
  IdGroupGAP = GAP.Globals.IdGroup(G1)
  IdGroup = GAP.gap_to_julia(Vector{Int}, IdGroupGAP)
  pinvariants = prod(invariants)
  if must_be_ram_surely
    #The extension must be ramified. Find a constant...
    cd = 1
    if iszero(mod(invariants[end], 2))
      #2 must be wildly ramified
      #The conductor must have at least valuation 2 at every prime over 2...
      cd = 2^pinvariants
    else
      #2 is not wildly ramified. Then we only have the boring bound...
      d = Int(minimum(keys(factor(invariants[end]).fac)))
      cd = 2^((d-1)*div(pinvariants, d))
    end 
    #But I want the minimum. So I have to look at the other primes..
    SP = PrimesSet(3, -1)
    for p in SP
      if p >= cd
        break
      end
      if iszero(mod(invariants[end], p))
        #p must be wildly ramified
        #The conductor must have at least valuation 2 at every prime over p...
        s = valuation(invariants[end], p)
        cd1 = p^(2*(p^s-p^(s-1))*divexact(pinvariants, p^s))
        if cd > cd1
          cd = cd1
        end
      else
        #p is not wildly ramified. Then we only have the boring bound...
        d = Int(minimum(keys(factor(invariants[end]).fac)))
        cd1 = p^((d-1)*div(pinvariants, d))
        if cd > cd1
          cd = cd1
        end
      end 
    end
    bound = root(div(absolute_bound, cd), prod(invariants))
  else
    bound = root(absolute_bound, prod(invariants))
  end
  list = fields(IdGroup[1], IdGroup[2], db, bound; only_real = (only_real || lvl == length(L)-1), unramified_outside = unramified_outside)
  if isempty(list)
    return list
  end
  @vprint :Fields 1 "computing extensions with Galois group ($a, $b) and bound ~10^$(clog(absolute_bound, 10))\n"
  @vprint :Fields 1 "Abelian invariants of the relative extension: $(invariants)\n"
  @vprint :Fields 1 "Number of fields at this step: $(length(list)) \n"
  @vprint :FieldsNonFancy 1 "Number of fields at this step: $(length(list)) \n"
  
  @vprint :Fields 1 "Computing obstructions\n"
  @vprint :FieldsNonFancy 1 "Computing obstructions\n"
  #@vtime :Fields 1 
  list = Hecke.check_obstruction(list, L, length(L)-1, invariants)
  @vprint :Fields 1 "Fields to check: $(length(list))\n\n"
  @vprint :FieldsNonFancy 1 "Fields to check: $(length(list))\n\n"
  if isempty(list)
    return FieldsTower[]
  end
  Id = GAP.Globals.IdGroup(G)
  fld_res = Hecke.field_extensions(list, absolute_bound, Id, invariants, only_real, unramified_outside = unramified_outside)
  #We insert the fields in the database
  if isempty(unramified_outside)
    if only_real
      insert_complete_table(db, map(number_field, fld_res), GP, absolute_bound, true, (a, 0), check = true)
    else
      insert_complete_table(db, map(number_field, fld_res), GP, absolute_bound, true, check = true)
    end
  else
    insert_fields(map(number_field, fld_res), db, check = true, galois_group = GP)
  end
  return fld_res
end