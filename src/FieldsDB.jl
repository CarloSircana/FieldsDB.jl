module FieldsDB

  import Oscar, Pkg
  export Oscar
  using LibPQ, Tables, PrettyTables, Decimals, Markdown
  
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
  include("./queries.jl")
  include("./insertion.jl")
  include("./completeness.jl")
  include("./fields.jl")
  include("./stats.jl")
end # module
