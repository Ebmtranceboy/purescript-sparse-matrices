module Data.Sparse.Matrix (module Data.Sparse.Matrix) where

import Prelude

import Data.Array ((..), (:))
import Data.Complex (Cartesian, real)
import Data.Foldable (foldr, sum, product)
import Data.FoldableWithIndex (foldrWithIndex)
import Data.Int (toNumber)
import Data.Map (Map, fromFoldable, toUnfoldable, keys, delete, update, isEmpty, alter, singleton)
import Data.Maybe (Maybe (..), maybe, fromMaybe)
import Data.Ord (abs)
import Data.Set (findMin)
import Data.Sparse.Polynomial 
  ( class Divisible
  , class Leadable
  , Polynomial (..)
  , (!)
  , (^)
  , roots
  , xchng
  )
import Data.Sparse.Polynomial (pow) as Poly
import Data.Tuple (Tuple (..))

-- | Representation of a matrix without storage of the zero coefficients.
-- | 
-- | 
-- | ``` 
-- | > -- Working with Sparse Matrices
-- | > import Data.Sparse.Matrix
-- | > import Data.Sparse.Polynomial ((^))
-- | > -- Define a matrix with sum of terms (see Data.Sparse.Polynomial)
-- | > a = Matrix { height:3, width:2, coefficients: 5.0^2^0+7.0^0^1 }
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
-- | > -- Coefficient extraction
-- | > c !! [1,0]
-- | 2.0
-- |
-- | > -- Column(s) extraction
-- | > c * Matrix { height: 3, width: 1, coefficients: 1.0^1^0 }
-- |  9.0
-- |  0.0
-- |  0.0
-- |
-- | > -- Row(s) extraction
-- | > Matrix { height: 2, width: 3, coefficients: 1.0^0^0 + 1.0^1^2 } * c
-- |  0.0 9.0
-- |  7.0 0.0
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
-- | > v = Matrix { height:3, width:1, coefficients: 1.0^0^0+2.0^1^0+3.0^2^0 }
-- | > v
-- |  1.0
-- |  2.0
-- |  3.0
-- |
-- | > -- Solve the linear system
-- | > pluSolve g v
-- |  0.03278688524590164
-- |  0.6666666666666666
-- |  -0.10784313725490192
-- |
-- | ``` 

type Poly2 a = Polynomial (Polynomial a)

newtype Matrix a = 
  Matrix 
    { height :: Int
    , width :: Int
    , coefficients :: Poly2 a
    }

-- | Coefficient extraction
extract :: forall a. Eq a => Semiring a => Matrix a -> Array Int -> a
extract (Matrix m) [i, j] = m.coefficients ! j ! i
extract _ _ = zero

-- | Coefficient extraction infix notation
infixl 8 extract as !!

instance showMatrix :: 
  ( Show a
  , Semiring a
  , Eq a
  ) => Show (Matrix a) where
  show mt@(Matrix m) = 
    let sizeLimit = 8
        isEndOfLine j = j == m.width - 1
        isDefaultI i = m.height > sizeLimit && i == sizeLimit / 2
        isDefaultJ j = m.width > sizeLimit && j == sizeLimit / 2
        isHidI i = 
          m.height > sizeLimit && 
          i > sizeLimit / 2 && 
          i < m.height - sizeLimit / 2
        isHidJ j = 
          m.width > sizeLimit && 
          j > sizeLimit / 2 && 
          j < m.width - sizeLimit / 2
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
             showElem (mt !! [i, j]) i j) 
               (0..(m.width-1))) (0..(m.height-1))
    
instance eqMatrix :: Eq a => Eq (Matrix a) where
  eq (Matrix a) (Matrix b) = a == b

-- | Returns the polynomial internal structure
internalMap :: forall a. Polynomial a -> Map Int a
internalMap (Poly m) = m

-- | Matrix transposition 
transpose :: forall a. 
  Eq a => 
  Semiring a => 
  Matrix a -> Matrix a
transpose (Matrix m) = 
  Matrix 
    { height: m.width
    , width: m.height
    , coefficients: xchng m.coefficients 
    }

trace :: forall a. Eq a => Semiring a => Square a -> a
trace a@(Matrix m) = sum $ (\i -> a !! [i,i]) <$> 0..(m.width-1)
  
instance semiringMatrix :: 
  ( Eq a
  , Semiring a
  ) => Semiring (Matrix a) where
  add (Matrix a) (Matrix b) = 
    Matrix 
      { height: max a.height b.height
      , width: max a.width b.width
      , coefficients: a.coefficients + b.coefficients
      }
  zero = Matrix { height: 1, width: 1, coefficients: zero }
  one = Matrix { height: 1, width: 1, coefficients: one }
  mul (Matrix a) (Matrix b) =
    Matrix { height: a.height, width: b.width, coefficients: coeffs }
      where Poly columns = b.coefficients
            coeffs :: Poly2 a
            coeffs = sum $ map (\(Tuple j (Poly w)) ->
               (_ ^ j) $ sum $ map (\ (Tuple i v) ->
                 Poly $ (fromFoldable :: Array (Tuple Int _) -> 
                    Map Int a) $ 
                   map (\(Tuple k u) -> Tuple k (u*v)) $ 
                     toUnfoldable (internalMap $ a.coefficients ! i) 
                     ) $ (toUnfoldable w :: Array (Tuple Int _)) 
                       ) $ (toUnfoldable columns :: 
                              Array (Tuple Int _))

instance functorMatrix :: Functor Matrix where
  map f (Matrix m) = 
    Matrix 
      { height: m.height
      , width: m.width
      , coefficients: map ((<$>) f) m.coefficients
      }

instance ringMatrix :: (Eq a, Ring a) => Ring (Matrix a) where
  sub a b = a + ((_ * (-one)) <$> b)

parseMonom :: forall a. 
  Eq a => Semiring a => 
  Poly2 a -> Maybe { i :: Int, j :: Int, v :: a }
parseMonom r@(Poly p) =  
  let j = fromMaybe (-1) $ findMin $ keys p 
      i = fromMaybe (-1) $ findMin $ keys $ internalMap $ r ! j
  in if i < 0 || j < 0
      then Nothing
      else Just { i, j, v: r ! j ! i }
 
-- | Element replacement
replace :: forall a. 
  Eq a => Semiring a =>
  Poly2 a -> Matrix a -> Matrix a
replace r a@(Matrix m) = 
    maybe a 
      ( \ { i, j, v } ->
        Matrix $ m { coefficients = replace' i j v m.coefficients } 
      ) $ parseMonom r
    
replace' :: forall a. 
  Eq a => Semiring a =>
  Int -> Int -> a -> Poly2 a -> Poly2 a
replace' i j v (Poly p) =  -- p + (- p ! j ! i + v)^i^j 
  if v == zero
    then 
      Poly $
        update (\ (Poly m) -> 
          let m' = delete i m
          in if isEmpty m'
                then Nothing
            else Just $ Poly m') j p
     else 
      Poly $ alter (\ arg -> 
        case arg of
             Nothing -> Just $ Poly $ singleton i v
             Just (Poly m) -> Just $ Poly $ alter (\ _ -> Just v) i m) j p

 
-- | Element replacement infix notation
infix 7 replace as ~

-- | Element increment
increment :: forall a. 
  Eq a => Semiring a => Ring a =>
  Poly2 a -> Matrix a -> Matrix a
increment r a@(Matrix m) = 
  maybe a 
    ( \ { i, j, v } -> 
      a + Matrix 
        { height: m.height
        , width: m.width
        , coefficients: v^i^j
        }
    ) $ parseMonom r

-- | Element increment infix notation
infix 7 increment as ~+

-- | Element weighting
scale :: forall a. Eq a => Semiring a => Ring a =>
  Poly2 a -> Matrix a -> Matrix a
scale r a@(Matrix m) = 
  maybe a 
    ( \ { i, j, v } -> 
      a + Matrix 
        { height: m.height
        , width: m.width
        , coefficients: (a!![i,j]*(v-one))^i^j
        }
    ) $ parseMonom r

-- | Element weighting infix notation
infix 7 scale as ~*

type Square a = Matrix a
type Vector a = Matrix a

-- | Identity matrix
eye :: forall a. 
  Eq a => Semiring a => 
  Int -> Square a
eye n = 
  let cs 0 = []
      cs k = one^(k-1)^(k-1) : cs (k-1)
   in 
     Matrix 
      { height: n
      , width: n
      , coefficients: sum $ cs n
      }  

-- | Zero matrix
zeros :: forall a. 
  Eq a => Semiring a => 
  Int -> Int -> Matrix a
zeros m n = 
  Matrix 
    { height: m
    , width: n
    , coefficients: zero 
    }  

height :: forall a. Matrix a -> Int
height (Matrix m) = m.height

width :: forall a. Matrix a -> Int
width (Matrix m) = m.width

coefficients :: forall a. Matrix a -> Poly2 a
coefficients (Matrix m) = m.coefficients

-- | As the name suggests. Note that the two rows are provided as 0-indices
swapRows :: forall a. 
  Eq a => 
  Divisible a => 
  Ring a => 
  Leadable a =>
  Int -> Int -> Matrix a -> Matrix a
swapRows ir is m =
  if ir > is 
     then swapRows is ir m
     else let r = ( Matrix 
                    { height: 1
                    , width: width m
                    , coefficients: one^ir^ir 
                    }) * m
              s = ( Matrix 
                    { height: 1
                    , width: width m
                    , coefficients: one^is^is 
                    }) * m
              delta = Poly.pow (one^1^0) (is-ir)
          in Matrix 
              { height: height m
              , width: width m
              , coefficients: 
                coefficients m - coefficients r - coefficients s 
                + coefficients r * delta
                + coefficients s / delta
              }

-- | `rowCombination orig ko dest kd` 
-- | performs : dest <- ko * orig + kd * dest
-- | which only alters dest
rowCombination :: forall a.
  Eq a => 
  Divisible a => 
  Ring a => 
  Leadable a =>
  Int -> a -> Int -> a -> Matrix a -> Matrix a
rowCombination orig ko dest kd m =
  if orig > dest 
     then let m' = swapRows dest orig m
          in swapRows dest orig $ rowCombination dest ko orig kd m' 
     else let o = ( Matrix 
                    { height: 1
                    , width: width m
                    , coefficients: one^orig^orig
                    }) * m
              d = ( Matrix
                    { height: 1
                    , width: width m
                    , coefficients: one^dest^dest 
                    }) * m
              delta = Poly.pow (one^1^0) (dest-orig)
          in Matrix 
              { height: height m
              , width: width m
              , coefficients: 
                coefficients m - coefficients d
                + ko^0^0 * coefficients o * delta
                + kd^0^0 * coefficients d
              }
   
-- | Returns the row index of the maximum element (in magnitude) 
-- | among all the elements below (possibly as high as) the diagonal element
-- | of the column whose index is provided.
pivot :: forall a.
  Eq a =>
  Ord a =>
  Ring a =>
  Int -> Matrix a -> Int
pivot col m@(Matrix r) =
  let column = r.coefficients ! col
  in foldr 
    ( \i acc -> 
        if abs(column ! i) > abs(column ! acc) 
          then i 
          else acc
    ) col $ col..(height m - 1)

-- | Returns the 3 square matrices P, L and U 
-- | and a Boolean such that
-- | * LU = PA where A is the given square matrix with not null determinant
-- | * P is a permutation matrix with parity given by the Boolean (true: 1, false: -1)
-- | * L is triangular inferior with unitary diagonal elements
-- | * U is triangular superior
plu :: forall a. 
  Eq a => 
  Ring a =>
  Ord a =>
  DivisionRing a =>
  Divisible a =>
  Leadable a =>
  Square a -> { p :: Square a, parity :: Boolean, l :: Square a, u :: Square a }
plu m =
  let n = height m
      go rec i
        | i == n-1  = rec { l = rec.p * rec.l }
        | otherwise = 
            let piv = pivot i rec.u
                swp = swapRows i piv (eye n)
                u' = swapRows i piv rec.u
                p' = swapRows i piv rec.p
                goLU = 
                  foldr 
                    ( \k { l, u } ->
                      let coef = u !! [k,i] * (recip $ u !! [i,i])
                      in 
                        { l: coef^k^i ~ l 
                        , u: rowCombination i (-coef) k one u 
                        }
                    ) { u: u', l: eye n } $ (i+1)..(n-1)
              in 
                go 
                  { p: p'
                  , parity: (i == piv) == rec.parity
                  , l: rec.l * swp * goLU.l
                  , u: goLU.u 
                  } 
                  (i+1)
    in go { p: eye n, parity: true, l: eye n, u: m } 0

pivot' :: forall a.
  Eq a =>
  Ord a =>
  Ring a =>
  Int -> Int -> Polynomial a -> Maybe Int
pivot' i h p =
  if i < h
    then 
      if p ! i == zero 
         then pivot' (i+1) h p 
         else Just i
    else Nothing

  -- | Determinant of a square matrix in a ring
ringDeterminant :: forall a. 
  Ring a => 
  Eq a => 
  Ord a => 
  Divisible a => 
  Leadable a => 
  EuclideanRing a => 
  Matrix a -> a

ringDeterminant m =
  let n = height m
      go rec i
        | i == n  = 
            let d = rec.piv in if rec.perm then d else -d
        | otherwise = 
            let mpiv = pivot' i n (coefficients rec.m ! i)
            in 
              maybe 
                zero 
                (\piv ->
                  let m' = coefficients (if i == piv then rec.m else swapRows i piv rec.m)
                      p' = if i == piv then rec.perm else not rec.perm
                      goM = 
                        foldr 
                          ( \k m'' ->
                            foldr
                              (\l m''' ->
                                replace' k l ((m''' ! i ! i * m''' ! l ! k - m''' ! i ! k * m''' ! l ! i) / rec.piv)
                                m'''
                              ) m'' $ (i+1)..(n-1)
                          ) m' $ (i+1)..(n-1)
                  in 
                    go 
                      { perm: p'
                      , m: Matrix { height: n, width: n, coefficients: goM }
                      , piv: m' ! i ! i 
                      } 
                      (i+1)
                  ) 
                  mpiv
    in go { perm: true, m, piv: one } 0

-- | Determinant of a square matrix in a field
determinant :: forall a. 
  Eq a => 
  Ord a => 
  Ring a => 
  DivisionRing a =>
  Divisible a =>
  Leadable a =>
  Square a -> a
determinant a@(Matrix m) =
  let { p: _p, parity, l: _l, u } = plu a
      d = product $ map (\i -> u!![i,i]) $ 0..(m.width-1)
  in if parity then d else -d

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
  in go n (Matrix { height: n, width: n, coefficients: zero }) 1.0 (1.0^n)

-- | Eigen values of a real square matrix
eigenValues :: Square Number -> Array (Cartesian Number)
eigenValues = roots <<< faddeev

-- | Integer power of a square matrix
pow :: forall a. Eq a => Semiring a => Square a -> Int -> Square a
pow m 0 = eye (height m)
pow m i 
  | i `mod` 2 == 0 = 
    let q = m `pow` (i `div` 2)
    in q * q
  | otherwise = 
    let q = m `pow` ((i-1) `div` 2)
    in m * q * q

-- | Polynomial application
applyPoly :: forall a. 
  Eq a => 
  Semiring a => 
  Polynomial (Square a) -> Square a -> Square a
applyPoly (Poly coeffs) m =
   foldrWithIndex (\i v acc -> acc + v * pow m i) zero coeffs

type Symmetric = Square Number

-- | Square real matrix diagonalization such that m = vec * val * recip vec
diagonalize :: Symmetric -> { val :: Symmetric, vec :: Square Number }
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
        in b * ( toMatrix $ 
                  foldr 
                    (\i acc -> acc + (toNumber i)^i^j) 
                    zero $ 
                    0..(n-1))
      vec = foldrWithIndex (\j v acc -> acc + f v j) (toMatrix zero) vs
  in { val, vec }

-- | Solves the system Ax=B for x, for A square inversible and any B
pluSolve :: forall a. 
  Eq a => 
  Ring a => 
  Ord a =>
  Divisible a =>
  Leadable a =>
  DivisionRing a =>
  Square a -> Matrix a -> Matrix a
pluSolve a b =
  let { p, parity: _, l, u } = plu a
  in luSolve l u (p*b)

-- | Solves the system LUx=B for x, for L,U square inversible and any B
luSolve :: forall a. 
  Eq a => 
  Ring a => 
  DivisionRing a =>
  Square a -> Square a -> Matrix a -> Matrix a
luSolve l u b = 
  let m = height b
      n = width b
      lSolve i j k acc ml
        | i == m = ml
        | k == i = lSolve (i+1) j 0 zero ((b!![i,j]-acc)^i^j ~ ml)
        | otherwise = lSolve i j (k+1) (acc + l!![i,k] * ml!![k,j]) ml
      uSolve i j k acc ml mu
        | i < 0 = { ml, mu }
        | k == i = uSolve (i-1) j (m-1) zero ml 
                          (((ml!![i,j]-acc) * (recip $ u!![i,i]))^i^j ~ mu)
        | otherwise = uSolve i j (k-1) (acc + u!![i,k] * mu!![k,j]) ml mu
      solve j ml mu  
        | j == n = mu
        | otherwise = 
           let { ml, mu } = uSolve (m-1) j (m-1) zero (lSolve 0 j 0 zero ml) mu
            in solve (j+1) ml mu
  in solve 0 (zeros m n) (zeros m n)

instance divisionRingMatrix :: 
  ( Eq a
  , DivisionRing a
  , Ord a
  , Divisible a
  , Leadable a
  ) => DivisionRing (Square a) where
  recip m = pluSolve m (eye $ width m)
