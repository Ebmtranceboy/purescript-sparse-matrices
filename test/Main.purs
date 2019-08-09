module Test.Main where

import Prelude
import Effect (Effect)
import Data.Sparse.Matrix(Matrix(..),(??),(~),(~+),(~*)
                          ,determinant,eye,transpose)
import Data.Sparse.Polynomial((^))
-- import Data.Complex(Cartesian, i, magnitudeSquared)
import Data.Ratio(Ratio, (%))

foreign import assert :: String -> Boolean -> Effect Unit

mat1 :: Matrix Int
mat1 = Matrix { height: 3, width: 2
              , coefficients: 1^0^0+2^0^1+3^1^0+4^1^1+5^2^0+6^2^1}
mat2 :: Matrix Int
mat2 = Matrix {height: 2, width: 1, coefficients: 7^0^0+9^1^0}

mat3 :: Matrix (Ratio Int)
mat3 = Matrix {height: 3, width: 3,
  coefficients: (2%1)^0^0 + (1%1)^0^1 + (3%1)^1^1 + (5%1)^2^2}

mat4 :: Matrix (Ratio Int)
mat4 = Matrix {height: 3, width: 3,
  coefficients: (1%1)^0^0 + (4%1)^1^1 + (9%1)^2^2}

main :: Effect Unit
main = do
  assert "matrix transpose" $ (transpose mat1)??[1,0] == 2
  assert "matrix replacement" $ 
    mat2 ?? [0,0] == 7 && ( 8^0^0 ~ mat2) ?? [0,0] == 8 
  assert "matrix element increment" $ 
    mat2 ?? [0,0] == 7 && ( 8^0^0 ~+ mat2) ?? [0,0] == 15 
  assert "matrix element scale" $ 
    mat3 ?? [0,0] == 2%1 && ( (1%4)^0^0 ~* mat3) ?? [0,0] == 1%2 
  assert "matrix sum" $ (mat3 + mat4) ?? [2,2] == 14%1
  assert "matrix difference" $ (mat3 - mat4) ?? [2,2] == -4%1
  let m@(Matrix mat12) = mat1 * mat2
  assert "matrix product" $ mat12.width == 1 
                         && mat12.height == 3 
                         && m ?? [2,0] == 89
  assert "square matrix determinant" $ determinant mat4 == 36%1 
  assert "square matrix inversion" $ mat3 * recip mat3 == eye 3
