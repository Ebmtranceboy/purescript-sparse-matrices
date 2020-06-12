module Data.Complex (module Data.Complex) where

import Prelude
import Data.Int(toNumber, round,ceil) as Int
import Math(sqrt, cos, sin, atan2)
import Math(pow) as Math

-- | Complex number defined by its real and imaginary parts
-- | ```
-- | z1 :: Cartesian Int
-- | z1 = Cartesian 1 2
-- | 
-- | z2 :: Cartesian Number
-- | z2 = Cartesian 2.0 (-3.0)
-- | 
-- | z3 :: Cartesian Number
-- | z3 = Cartesian 8.0 1.0
-- | 
-- | z4 :: Cartesian Int
-- | z4 = Cartesian 5 (-4)
-- | 
-- | main :: Effect Unit
-- | main = do
-- |   assert "Cartesian Int is showable" $ show z1 == "1+2i"
-- |   assert "Cartesian Number is showable" $ show z2 == "2.0-3.0i"
-- |   assert "Cartesian Number is equatable" $ z1 == z1 && z2 == z2
-- |   assert "Cartesian components" $ real z1 == 1 && imag z1 == 2
-- |   assert "Complex basis" $ one + ((_*2)<$>i) == z1
-- |   assert "Complex subtraction" $ z3-z2 == Cartesian 6.0 4.0
-- |   assert "Complex conjugaison" $ conjugate z1 == Cartesian 1 (-2)
-- |   assert "Complex division" $ z3/z2 == Cartesian 1.0 2.0
-- |   let n = magnitudeSquared $ normalize z3
-- |   assert "Complex normalization" $ 1.0 - n < 1e-6  
-- |   let g = gcd z1 z4
-- |   let z1' = z1/g
-- |   let z4' = z4/g
-- |   let m = mod z1 z4
-- |   assert "Gauss integers" $ z1' * g == z1 
-- |                          && z4' * g == z4
-- |                          && magnitudeSquared m < magnitudeSquared z4
-- |   assert "Complex power" $ 
-- |     magnitudeSquared (pow (Cartesian 1.0 1.0) 2.0 - ((_*2.0) <$>i)) < 1e-6 
-- | ```

data Cartesian a = Cartesian a a

-- | Real part
real :: forall a. Cartesian a -> a
real (Cartesian x _) = x

-- | Imaginary part
imag :: forall a. Cartesian a -> a
imag (Cartesian _ y) = y

-- | Imaginary unit
i :: forall a. Semiring a => Cartesian a
i = Cartesian zero one

instance showCartesian :: (Show a, Ord a, Semiring a, Ring a) => Show (Cartesian a) where
  show (Cartesian a b) = show a <> (if b<zero then "-" <> show (negate b) else "+" <> show b) <> "i"  

instance eqCartesian :: Eq a => Eq (Cartesian a) where
  eq (Cartesian a b) (Cartesian c d) = a==c && b==d

instance semiringCartesian :: Ring a => Semiring (Cartesian a) where
  add (Cartesian a b) (Cartesian c d) = Cartesian (a+c) (b+d)
  zero = Cartesian zero zero
  mul (Cartesian a b) (Cartesian c d) = Cartesian (a*c-b*d) (a*d+b*c)
  one = Cartesian one zero

instance functorCartesian :: Functor Cartesian where
  map f (Cartesian x y) = Cartesian (f x) (f y)

instance ringCartesian :: Ring a => Ring (Cartesian a) where
  sub z1 z2 = add z1 (map (_ * (sub zero one)) z2)

-- | Conjugate
conjugate :: forall a. Ring a => Cartesian a -> Cartesian a
conjugate (Cartesian x y) = Cartesian x (negate y)

-- | Magnitude Squared
magnitudeSquared :: forall a. Ring a => Cartesian a -> a
magnitudeSquared z = real $ z * (conjugate z)

-- | Normalize to norm 1
normalize :: Cartesian Number -> Cartesian Number
normalize z = map (_ / (sqrt $ magnitudeSquared z)) z

instance commutativeRingCartesian :: CommutativeRing a => CommutativeRing (Cartesian a)

instance divisionRingCartesian :: DivisionRing a => DivisionRing (Cartesian a) where
  recip z = map (_ * (recip $ magnitudeSquared z)) (conjugate z)

instance euclideanRingCartesianNumber :: EuclideanRing (Cartesian Number) where
  degree = Int.ceil <<< magnitudeSquared
  div z z' = z * (recip z')
  mod z z' = zero

instance euclideanRingCartesianInt :: EuclideanRing (Cartesian Int) where
  degree = magnitudeSquared
  div z z' = map Int.round $ div (map Int.toNumber z) (map Int.toNumber z')
  mod z z' = sub z $ mul z' $ div z z' 

-- | From radius and angle in radians
fromPolar :: Number -> Number -> Cartesian Number
fromPolar r theta = Cartesian (r * cos theta) (r * sin theta)

-- | Angle in radians
angle :: Cartesian Number -> Number
angle (Cartesian x y) = atan2 y x

-- | Real power of a complex
pow :: Cartesian Number -> Number -> Cartesian Number
pow z n = fromPolar (Math.pow (sqrt $ magnitudeSquared z) n) (angle z * n)

