cnt = 0; while true
  @show global cnt +=1
  @time GC.gc(true)
  l = FieldsDB._get_fields_for_class_group_computation(db)
  for x in l
    @show K = number_field(x)
    automorphisms(K)
    if degree(K) > 30
      OK = lll(maximal_order(K))
      fl, N = norm_relation(K, small_degree = false)
      if fl
        @time Hecke.NormRel.class_group_via_brauer(OK, N)
      end
    end
    @time FieldsDB.set_class_group_and_regulator(x, true)
  end
end
