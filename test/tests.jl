@testset "Interface" begin

  #Now, we can test everything.
  FieldsDB.fields(8, 4, db, fmpz(10)^12)
  flds = load_fields(db)
  @test length(flds) == 138
  for i = 1:length(flds)
    FieldsDB.set_ramified_primes(flds[i])
    FieldsDB.set_class_group(flds[i])
    FieldsDB.set_regulator(flds[i])
    FieldsDB.set_subfields(flds[i])
    FieldsDB.set_galois_group(flds[i])
    FieldsDB.set_canonical_defining_polynomial(flds[i])
    FieldsDB.set_automorphisms_order(flds[i])
    FieldsDB.set_iscm(flds[i])
    FieldsDB.set_torsion_unit_size(flds[i])
  end
  @test load_fields(db, signature = (0, 4), only_count = Val{true}) == 1
  @test load_fields(db, signature = (8, 0), only_count = Val{true}) == 1
  @test load_fields(db, signature = (4, 0), only_count = Val{true}) == 136
  @test load_fields(db, discriminant_range = (fmpz(12230590464), fmpz(12230590464)), only_count = Val{true}) == 2
  @test load_fields(db, signature = (8, 0), class_group_ranks_range = Dict(fmpz(2) => (2, 2)), only_count = Val{true}) == 0
  @test load_fields(db, signature = (0, 4), class_group_ranks_range = Dict(fmpz(2) => (1, 2)), only_count = Val{true}) == 1
  GP = FieldsDB.isomorphic_transitive_perm_group(small_group(4, 2), 4)
  @test load_fields(db, galois_group = GP, only_count = Val{true}) == 136
end
  

