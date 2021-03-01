module FieldsDB


  export Hecke, Oscar
  using LibPQ, Tables, Oscar
  
  exclude = Symbol[:DirectProductOfElem, :coset_decomposition, :disc_log,
  :generator, :isdihedral_group, :visual, :weight  ]

  for i in names(Oscar)
    if i in exclude
      continue
    end
    eval(Meta.parse("import Oscar." * string(i)))
    eval(Expr(:export, i))
  end

  include("./interface.jl")
  include("./fields.jl")
end # module
