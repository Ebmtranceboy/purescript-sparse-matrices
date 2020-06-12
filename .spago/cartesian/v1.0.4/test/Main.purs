module Test.Main where

import Prelude (Unit, discard, gcd, mod, negate, one, show, ($), (&&), (*), (+), (-), (/), (<), (<$>), (==))
import Effect (Effect)
import Data.Complex

foreign import assert :: String -> Boolean -> Effect Unit

z1 :: Cartesian Int
z1 = Cartesian 1 2

z2 :: Cartesian Number
z2 = Cartesian 2.0 (-3.0)

z3 :: Cartesian Number
z3 = Cartesian 8.0 1.0

z4 :: Cartesian Int
z4 = Cartesian 5 (-4)

main :: Effect Unit
main = do
  assert "Cartesian Int is showable" $ show z1 == "1+2i"
  assert "Cartesian Number is showable" $ show z2 == "2.0-3.0i"
  assert "Cartesian Number is equatable" $ z1 == z1 && z2 == z2
  assert "Cartesian components" $ real z1 == 1 && imag z1 == 2
  assert "Complex basis" $ one + ((_*2)<$>i) == z1
  assert "Complex subtraction" $ z3-z2 == Cartesian 6.0 4.0
  assert "Complex conjugaison" $ conjugate z1 == Cartesian 1 (-2)
  assert "Complex division" $ z3/z2 == Cartesian 1.0 2.0
  let n = magnitudeSquared $ normalize z3
  assert "Complex normalization" $ 1.0 - n < 1e-6  
  let g = gcd z1 z4
  let z1' = z1/g
  let z4' = z4/g
  let m = mod z1 z4
  assert "Gauss integers" $ z1' * g == z1 
                         && z4' * g == z4
                         && magnitudeSquared m < magnitudeSquared z4
  assert "Complex power" $ 
    magnitudeSquared (pow (Cartesian 1.0 1.0) 2.0 - ((_*2.0) <$>i)) < 1e-6 
