module Data.Sparse.Matrix(module Data.Sparse.Matrix) where

import Prelude
import Data.Sparse.Polynomial(Polynomial(..),(?),roots)
import Data.Sparse.Polynomial(monoPol) as Poly
import Data.Maybe(Maybe(..),fromMaybe,maybe)
import Data.Map(Map,toUnfoldable, fromFoldable, keys)
import Data.Foldable(foldr,sum,product)
import Data.FoldableWithIndex (foldrWithIndex)
import Data.Array((..),(:))
import Data.Tuple(Tuple(..))
import Data.Set(findMin)
import Data.Int (toNumber)
import Data.Complex (Cartesian, real)

-- | Imported from Data.Sparse.Polynomial
monoPol :: forall a. a -> Int -> Polynomial a
monoPol = Poly.monoPol
infixl 8 monoPol as ^

type Poly2 a = Polynomial (Polynomial a)

-- | Representation of a matrix without storage of the zero coefficients.
-- | 
-- | 
-- | ``` 
-- | > -- Working with Sparse Matrices
-- | > import Data.Sparse.Matrix
-- | > -- Define a matrix with sum of terms (see Data.Sparse.Polynomial)
-- | > a = Matrix {height:3, width:2, coefficients: 5.0^2^0+7.0^0^1}
-- | > a
-- |  0.0 7.0
-- |  0.0 0.0
-- |  5.0 0.0
-- |
-- | > -- Modify all the non-zero elements with the functor
-- | > b = (_ + 2.0) <$> a
-- | > b
-- |  0.0 9.0
-- |  0.0 0.0
-- |  7.0 0.0
-- |
-- | > -- Set a specific element
-- | > c = 2.0^1^0 ~ b
-- | > c
-- |  0.0 9.0
-- |  2.0 0.0
-- |  7.0 0.0
-- |
-- | > -- Inspect a specific location
-- | > c??[1,0]
-- | 2.0
-- |
-- | > -- transpose, multiply
-- | > d = a * transpose c
-- | > d
-- |  63.0 0.0 0.0
-- |  0.0 0.0 0.0
-- |  0.0 10.0 35.0
-- |
-- | > -- Weight a specific element
-- | > e = 0.5^0^0 ~* d
-- | > e
-- |  31.5 0.0 0.0
-- |  0.0 0.0 0.0
-- |  0.0 10.0 35.0
-- |
-- |
-- | > -- Increment a specific element
-- | > f = 4.0^1^1 ~+ e
-- | > f
-- |  31.5 0.0 0.0
-- |  0.0 4.0 0.0
-- |  0.0 10.0 35.0
-- |
-- | > -- Identity, subtraction
-- | > g = f - eye 3
-- | > g
-- |  30.5 0.0 0.0
-- |  0.0 3.0 0.0
-- |  0.0 10.0 34.0
-- |
-- | > -- Vectors are just column matrices
-- | > v = Matrix {height:3, width:1, coefficients: 1.0^0^0+2.0^1^0+3.0^2^0}
-- | > v
-- |  1.0
-- |  2.0
-- |  3.0
-- |
-- | > -- Solve the linear system
-- | > luSolve g v
-- |  0.03278688524590164
-- |  0.6666666666666666
-- |  -0.1078431372549019
-- |
-- | ``` 
newtype Matrix a = 
  Matrix { height :: Int
         , width :: Int
         , coefficients :: Poly2 a}



-- | Coefficient extraction
extract :: forall a. Eq a => Semiring a => Matrix a -> Array Int -> a
extract (Matrix m) [i, j] = m.coefficients ? j ? i
extract _ _ = zero

-- | Coefficient extraction infix notation
infixl 8 extract as ??

instance showMatrix :: (Show a, Semiring a, Eq a) => 
  Show (Matrix a) where
  show mt@(Matrix m) = 
    let isEndOfLine j = j == m.width - 1
        isDefaultI i = m.height > 6 && i == 3
        isDefaultJ j = m.width > 6 && j == 3
        isHidI i = m.height > 6 && i > 3 && i < m.height - 3
        isHidJ j = m.width > 6 && j > 3 && j < m.width - 3
        showElem e i j = next
          where
          next
            | isEndOfLine j && not (isHidI i) 
                               && not (isDefaultI i) =
                                 " " <> show e <> "\n"
            | isEndOfLine j && isDefaultI i = "\n"
            | isDefaultI i && isDefaultJ j = ""
            | isDefaultI i = " ."
            | isHidI i = ""
            | isDefaultJ j = " . . ."
            | isHidJ j = ""
            | otherwise = " " <> show e
      in foldr (<>) "" $ map (\i -> 
           foldr (<>) "" $ map (\j -> 
             showElem (mt ?? [i, j]) i j) 
               (0..(m.width-1))) (0..(m.height-1))

instance eqMatrix :: Eq a => Eq (Matrix a) where
  eq (Matrix a) (Matrix b) = a==b

-- | Returns the polynomial internal structure
internalMap :: forall a. Polynomial a -> Map Int a
internalMap (Poly m) = m

-- | Matrix transposition 
transpose :: forall a. Eq a => Semiring a => Ring a =>  Matrix a -> Matrix a
transpose a@(Matrix m) = 
  let f i j b
        | i == m.height && j == m.width = b
        | i == m.height = f 0 (j+1) b
        | otherwise = f (i+1) j ((a??[i,j])^j^i ~ b)
   in f 0 0 $ Matrix {height: m.width, width: m.height, coefficients: zero}

trace :: forall a. Eq a => Semiring a => Matrix a -> a
trace a@(Matrix m) = sum $ (\i -> a??[i,i]) <$> 0..(m.width-1)
  

instance semiringMatrix :: (Eq a, Semiring a) => 
  Semiring (Matrix a) where
  add (Matrix a) (Matrix b) = 
    Matrix { height: max a.height b.height
           , width: max a.width b.width
           , coefficients: a.coefficients + b.coefficients}
  zero = Matrix {height: 1, width: 1, coefficients: zero}
  one = Matrix {height:1, width:1, coefficients: one}
  mul (Matrix a)(Matrix b) =
    Matrix { height: a.height, width: b.width , coefficients}
      where Poly columns = b.coefficients
            coefficients :: Poly2 a
            coefficients = sum $ map (\(Tuple j (Poly w)) ->
               (_ ^ j) $ sum $ map (\ (Tuple i v) ->
                 Poly $ (fromFoldable :: Array (Tuple Int _) -> 
                    Map Int a) $ 
                   map (\(Tuple k u) -> Tuple k (u*v)) $ 
                     toUnfoldable (internalMap $ a.coefficients ? i) 
                     ) $ (toUnfoldable w :: Array (Tuple Int _)) 
                       ) $ (toUnfoldable columns :: 
                              Array (Tuple Int _))

instance functorMatrix :: Functor Matrix where
  map f (Matrix m) = 
    Matrix { height: m.height
           , width: m.width
           , coefficients: map ((<$>) f) m.coefficients}

instance ringMatrix :: (Eq a, Ring a) => Ring (Matrix a) where
  sub a b = a + ((_ * (-one)) <$> b)

parseMonom :: forall a. Eq a => Semiring a => 
  Poly2 a -> Maybe {i :: Int, j :: Int, v :: a}
parseMonom r@(Poly p) =  
  let j = fromMaybe (-1) $ findMin $ keys p 
      i = fromMaybe (-1) $ findMin $ keys $ internalMap $ r ? j
   in if i<0 || j<0
      then Nothing
      else Just {i,j,v: r ? j ? i}
 
-- | Element replacement
replace :: forall a. Eq a => Semiring a => Ring a =>
  Poly2 a -> Matrix a -> Matrix a
replace r a@(Matrix m) = 
  maybe a (\{i,j,v} -> a + Matrix { height: m.height
                 , width: m.width
                 , coefficients: (v - a ?? [i,j])^i^j}) $
                   parseMonom r
-- | Element replacement infix notation
infix 7 replace as ~

-- | Element increment
increment :: forall a. Eq a => Semiring a => Ring a =>
  Poly2 a -> Matrix a -> Matrix a
increment r a@(Matrix m) = 
  maybe a (\{i,j,v} -> a + Matrix { height: m.height
                 , width: m.width
                 , coefficients: v^i^j}) $ parseMonom r

-- | Element increment infix notation
infix 7 increment as ~+

-- | Element weighting
scale :: forall a. Eq a => Semiring a => Ring a => EuclideanRing a =>
  Poly2 a -> Matrix a -> Matrix a
scale r a@(Matrix m) = 
  maybe a (\{i,j,v} -> a + Matrix { height: m.height
                 , width: m.width
                 , coefficients: (a??[i,j]*(v-one))^i^j}) $ parseMonom r

-- | Element weighting infix notation
infix 7 scale as ~*

type Square a = Matrix a
type Vector a = Matrix a

-- | Identity matrix
eye :: forall a. Eq a => Semiring a => Int -> Square a
eye n = 
  let cs 0 = []
      cs k = one^(k-1)^(k-1) : cs (k-1)
   in Matrix {height: n, width: n, coefficients: sum $ cs n}  

-- | Zero matrix
zeros :: forall a. Eq a => Semiring a => Int -> Int -> Matrix a
zeros m n = Matrix {height: m, width: n, coefficients: zero}  

height :: forall a. Matrix a -> Int
height (Matrix m) = m.height

width :: forall a. Matrix a -> Int
width (Matrix m) = m.width

-- | Returns the 2 square matrices L and U such that
-- | * LU = A where A is a square matrix
-- | * L is triangular inferior with unitary diagonal elements
-- | * U is triangular superior
lu :: forall a. Eq a => Semiring a => Ring a => EuclideanRing a =>
  Square a -> {l :: Square a, u :: Square a}
lu a = 
  let n = height a
      go l_ u_ i 
        | i==n      = {l:l_, u:u_}
        | otherwise = 
            let row u' j 
                  | j==n      = u'
                  | otherwise =
                      let row' u'' q 
                            | q==i      = row u'' (j+1)
                            | otherwise = 
                                let x = -l_??[i,q] * u''??[q,j]
                                in row' (x^i^j ~+ u'') (q+1)
                        in row' ((a??[i,j])^i^j ~ u') 0
                column l' u' j
                  | j == n    = l'
                  | otherwise =
                      let column' l'' q
                            | q==i      = column ((one/u'??[i,i])^j^i ~* l'') 
                                                 u' (j+1)
                            | otherwise = 
                                let x = -l''??[j,q] * u'??[q,i]
                                in column' (x^j^i ~+ l'') (q+1)
                      in column' ((a??[j,i])^j^i ~ l') 0
                u = row u_ i
                l = column l_ u (i+1)
            in go l u (i+1)
  in go (eye n) (eye n) 0

-- | Determinant of a square matrix
determinant :: forall a. Eq a => Semiring a => Ring a => EuclideanRing a =>
  Square a -> a
determinant a@(Matrix m) =
  let { l: _l, u } = lu a
   in product $ map (\i -> u??[i,i]) $ 0..(m.width-1)

-- | Characteristic polynomial of a real square matrix
faddeev :: Square Number -> Polynomial Number
faddeev a =
  let n = height a
      go 0 _ _ p = p
      go i m c p =
        let k = n - i + 1
            am = a * m + ((_ * c) <$> eye n)
            coef = - trace (a * am) / toNumber k
        in go (i-1) am coef (p + coef^(i-1))
  in go n (Matrix {height: n, width: n, coefficients: zero}) 1.0 (1.0^n)

-- | Eigen vlaues of a real square matrix
eigenValues :: Square Number -> Array (Cartesian Number)
eigenValues = roots <<< faddeev

-- | Integer power of a square matrix
pow :: forall a. Eq a => Semiring a => Square a -> Int -> Square a
pow m 0 = eye (height m)
pow m i = m * pow m (i-1)

-- | Polynomial application
applyPoly :: forall a. Eq a => Semiring a => Polynomial (Square a) -> Square a -> Square a
applyPoly (Poly coeffs) m =
   foldrWithIndex (\i v acc -> acc + v * pow m i) zero coeffs

type Symmetric = Square Number

-- | Square real matrix diagonalization such that m = vec * val * recip vec
diagonalize :: Symmetric -> { val :: Symmetric, vec :: Square Number}
diagonalize m =
  let n = height m
      vs = real <$> eigenValues m
      toMatrix p = Matrix { height: n
                  , width: n
                  , coefficients: p
                  }
      val = toMatrix 
            $ foldrWithIndex (\i v acc -> acc + v^i^i) 
                            zero 
                            vs
      fromCst c = (_ * c) <$> eye n
      pol = faddeev m
      f v j =
        let p = fromCst <$> pol / ((-v)^0+1.0^1)
            b = applyPoly p m 
          in b * (toMatrix $ foldr (\i acc -> acc + (toNumber i)^i^j) zero $ 0..(n-1))
      vec = foldrWithIndex (\j v acc -> acc + f v j) (toMatrix zero) vs
  in {val, vec}
  

-- | Solves the system Ax=B for x, for A square inversible and any B
luSolve :: forall a. Eq a => Semiring a => Ring a => EuclideanRing a =>
  Square a -> Matrix a -> Matrix a
luSolve a b = 
  let m = height b
      n = width b
      {l,u} = lu a
      lSolve i j k acc ml
        | i == m = ml
        | k == i = lSolve (i+1) j 0 zero ((b??[i,j]-acc)^i^j ~ ml)
        | otherwise = lSolve i j (k+1) (acc + l??[i,k] * ml??[k,j]) ml
      uSolve i j k acc ml mu
        | i < 0 = {ml, mu}
        | k == i = uSolve (i-1) j (m-1) zero ml 
                          (((ml??[i,j]-acc)/u??[i,i])^i^j ~ mu)
        | otherwise = uSolve i j (k-1) (acc + u??[i,k] * mu??[k,j]) ml mu
      solve j ml mu  
        | j==n = mu
        | otherwise = 
           let {ml, mu} = uSolve (m-1) j (m-1) zero (lSolve 0 j 0 zero ml) mu
            in solve (j+1) ml mu
  in solve 0 (zeros m n) (zeros m n)

instance divisionRingMatrix :: (Eq a, DivisionRing a, EuclideanRing a) => 
  DivisionRing (Matrix a) where
  recip m = luSolve m (eye $ width m)


