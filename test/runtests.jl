using FieldsDB
using Test
using LibPQ, Oscar.Hecke

db = LibPQ.Connection("host=localhost dbname=fields_test port=5432 user=postgres")

  #We create temporary tables.
  execute(db, """
    CREATE TABLE fields.class_group(
      class_group_id SERIAL,
      group_order NUMERIC NOT NULL,
      structure NUMERIC[] NOT NULL,
      prime_divisors NUMERIC[] NOT NULL,
      ranks SMALLINT[] NOT NULL, 
      PRIMARY KEY(class_group_id)
    );

    CREATE TABLE fields.group(
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

    CREATE TABLE fields.completeness(
      GRH BOOLEAN NOT NULL,
      group_id INT REFERENCES fields.group(group_id) NOT NULL,
      real_embeddings SMALLINT NOT NULL,
      discriminant_bound NUMERIC NOT NULL, 
      PRIMARY KEY (GRH, group_id, real_embeddings)
    );


    CREATE TABLE fields.field(
      field_id BIGSERIAL PRIMARY KEY,
      polynomial NUMERIC[] UNIQUE NOT NULL,
      degree SMALLINT NOT NULL,
      real_embeddings SMALLINT NOT NULL,
      class_group_id INT REFERENCES fields.class_group(class_group_id),
      ramified_primes NUMERIC[],
      regulator NUMERIC,
      discriminant NUMERIC NOT NULL,
      GRH BOOLEAN,
      group_id INT REFERENCES fields.group(group_id),
      CM BOOLEAN, 
      torsion_size INT, 
      automorphisms_order SMALLINT,
      is_canonical_poly BOOLEAN,
      subfields BIGINT[]
    );
  """)
  
try
  @time include("tests.jl")
finally 
  execute(db, """
    DROP TABLE fields.field;
    DROP TABLE fields.completeness;
    DROP TABLE fields.class_group;
    DROP TABLE fields.group;
  """)
  close(db)
end