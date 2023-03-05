{-
Welcome to a Spago project!
You can edit this file as you like.
-}
{ name = "my-project"
, dependencies =
  [ "arrays"
  , "assert"
  , "cartesian"
  , "effect"
  , "foldable-traversable"
  , "integers"
  , "maybe"
  , "numbers"
  , "ordered-collections"
  , "prelude"
  , "rationals"
  , "sparse-polynomials"
  , "tuples"
  ]
, packages = ./packages.dhall
, sources = [ "src/**/*.purs", "test/**/*.purs" ]
}
