module Test.Main where

import Prelude
import Effect (Effect)
import Data.Number (abs)
import Data.Sparse.Matrix 
  ( Matrix(..)
  , (??)
  , (~)
  , (~+)
  , (~*)
  , determinant
  , ringDeterminant
  , eye
  , transpose
  , trace
  , diagonalize
  , plu
  )
import Data.Sparse.Polynomial (Polynomial, (^))
import Data.Ratio (Ratio, (%))
import Test.Assert (assert')

mat1 :: Matrix Int
mat1 = 
  Matrix 
    { height: 3
    , width: 2
    , coefficients: 1^0^0 + 2^0^1 
                  + 3^1^0 + 4^1^1
                  + 5^2^0  +6^2^1
    }
mat2 :: Matrix Int
mat2 = 
  Matrix 
    { height: 2
    , width: 1
    , coefficients: 7^0^0 + 9^1^0
    }

mat3 :: Matrix (Ratio Int)
mat3 = 
  Matrix 
    { height: 3
    , width: 3
    , coefficients: (2%1)^0^0 + (1%1)^0^1 
                              + (3%1)^1^1 
                                          + (5%1)^2^2
    }

mat4 :: Matrix (Ratio Int)
mat4 = 
    Matrix 
      { height: 3
      , width: 3
      , coefficients: (1%1)^0^0 + (4%1)^1^1 + (9%1)^2^2
      }

sym :: Polynomial (Polynomial Number)
sym = 1.0^0^0 + 2.0^0^1 + 3.0^0^2
    + 2.0^1^0 + 5.0^1^1 + 6.0^1^2
    + 3.0^2^0 + 6.0^2^1 + 9.0^2^2
      
mat5 :: Matrix Number
mat5 = Matrix { height:3, width:3, coefficients: sym }

type P3 = Polynomial (Polynomial (Polynomial (Ratio Int)))

main :: Effect Unit
main = do
  assert' "matrix transpose" $ (transpose mat1)??[1,0] == 2
  assert' "matrix replacement" $ 
    mat2 ?? [0,0] == 7 && ( 8^0^0 ~ mat2) ?? [0,0] == 8 
  assert' "matrix element increment" $ 
    mat2 ?? [0,0] == 7 && ( 8^0^0 ~+ mat2) ?? [0,0] == 15 
  assert' "matrix element scale" $ 
    mat3 ?? [0,0] == 2%1 && ( (1%4)^0^0 ~* mat3) ?? [0,0] == 1%2 
  assert' "matrix trace" $ trace mat3 == 10%1
  assert' "matrix sum" $ (mat3 + mat4) ?? [2,2] == 14%1
  assert' "matrix difference" $ (mat3 - mat4) ?? [2,2] == -4%1
  let m@(Matrix mat12) = mat1 * mat2
  assert' "matrix product" $ mat12.width == 1 
                         && mat12.height == 3 
                         && m ?? [2,0] == 89
  assert' "square matrix ring determinant" $ ringDeterminant (transpose mat1 * mat1) == 24
  assert' "square matrix division-ring determinant" $ determinant mat4 == 36%1 
  let a = transpose mat3 - ((_ * (2%1)) <$> mat4)
      { p, parity: _, l, u } = plu a
  assert' "LU decomposition" $ l * u == p * a 
  assert' "square matrix inversion" $ mat3 * recip mat3 == eye 3
  let {val, vec} = diagonalize mat5
  assert' "real diagonalization" $ 
    (abs $ determinant $ mat5 - vec * val * recip vec) < 1e-12
  
  let x :: P3
      x = (1%1) ^ 0 ^ 0 ^ 1
      y :: P3
      y = (1%1) ^ 0 ^ 1 ^ 0
      z :: P3
      z = (1%1) ^ 1 ^ 0 ^ 0
      
  assert' "symbolic computation" $ 
    ( ringDeterminant $ 
      Matrix 
        { height: 2
        , width: 2
        , coefficients:   x^0^0 + z^0^1
                        + z^1^0 + y^1^1 
        }
    ) == x*y-z*z

