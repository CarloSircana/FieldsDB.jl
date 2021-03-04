using FieldsDB
using Test
using LibPQ, Oscar.Hecke

DATABASE_USER = get(ENV, "LIBPQJL_DATABASE_USER", "postgres")
if get(ENV, "CI", "false") == "true"
  db = LibPQ.Connection("dbname=postgres user=$DATABASE_USER")
else
  db = LibPQ.Connection("host=localhost dbname=postgres user=$DATABASE_USER")
end

  #We create temporary tables.
  execute(db, """
    CREATE TEMP TABLE class_group(
      class_group_id SERIAL,
      group_order NUMERIC NOT NULL,
      structure NUMERIC[] NOT NULL,
      prime_divisors NUMERIC[] NOT NULL,
      ranks SMALLINT[] NOT NULL, 
      PRIMARY KEY(class_group_id)
    );

    CREATE TEMP TABLE galois_group(
      group_id SERIAL PRIMARY KEY,
      group_order NUMERIC NOT NULL, 
      degree INT NOT NULL, 
      transitive_group_id INT,
      small_group_id INT,
      generators VARCHAR UNIQUE,
      abelian BOOLEAN NOT NULL,
      nilpotent BOOLEAN NOT NULL, 
      solvable BOOLEAN NOT NULL,
      primitive BOOLEAN NOT NULL,
      perfect BOOLEAN NOT NULL,
      issimple BOOLEAN NOT NULL
    );

    CREATE TEMP TABLE completeness(
      GRH BOOLEAN NOT NULL,
      group_id INT REFERENCES galois_group(group_id) NOT NULL,
      real_embeddings SMALLINT NOT NULL,
      discriminant_bound NUMERIC NOT NULL, 
      PRIMARY KEY (GRH, group_id, real_embeddings)
    );


    CREATE TEMP TABLE field(
      field_id BIGSERIAL PRIMARY KEY,
      polynomial NUMERIC[] UNIQUE NOT NULL,
      degree SMALLINT NOT NULL,
      real_embeddings SMALLINT NOT NULL,
      class_group_id INT REFERENCES class_group(class_group_id),
      ramified_primes NUMERIC[],
      regulator NUMERIC,
      discriminant NUMERIC NOT NULL,
      GRH BOOLEAN,
      group_id INT REFERENCES galois_group(group_id),
      CM BOOLEAN, 
      torsion_size INT, 
      automorphisms_order SMALLINT,
      is_canonical_poly BOOLEAN,
      subfields BIGINT[]
    );
  """)

@time include("tests.jl")