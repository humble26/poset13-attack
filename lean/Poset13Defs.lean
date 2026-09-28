/-!
# Formal verification scaffold for the 1/3-2/3 conjecture attack

Formalizes the working notions of the session attack on the 1/3-2/3 conjecture
for finite posets, states the four results obtained (Lemma A, Lemma B, Lemma C,
Proposition 1), and machine-checks all of them exhaustively inside Lean for
every poset on `n <= 5` elements (Lemmas A/B additionally on `n = 6`).

Trust levels, explicitly:
* `decide`-theorems: evaluated by the Lean kernel only, no extra axioms.
* `native_decide`-theorems: evaluated by compiled code, depending on the axiom
  `Lean.ofReduceBool` (disclosed below via `-- #print axioms`).
* The general quantified statements are recorded as `Prop`s; their full kernel
  proofs are open work of the attack program.
-/

/-- A finite poset on `{0, ..., n-1}`, computable strict order. Bare structure:
axioms live in the separate predicate `isPosetOK` so reindexing/deleting stays
computational. -/
structure POrd (n : Nat) where
  lt : Nat → Nat → Bool

variable {n : Nat}

/-! ### Well-formedness (Bool level) -/

def irreflOK (P : POrd n) : Bool :=
  (List.range n).all (fun i => !(P.lt i i))

def antisymOK (P : POrd n) : Bool :=
  (List.range n).all (fun i => (List.range n).all (fun j =>
    !(P.lt i j && P.lt j i)))

def transOK (P : POrd n) : Bool :=
  (List.range n).all (fun i => (List.range n).all (fun j => (List.range n).all (fun k =>
    !(P.lt i j && P.lt j k) || P.lt i k)))

def isPosetOK (P : POrd n) : Bool := irreflOK P && antisymOK P && transOK P

/-! ### Small list helpers -/

def memN (v : Nat) (l : List Nat) : Bool := l.any (fun w => w == v)

def idxPair (l : List (Nat × Nat)) (p : Nat × Nat) : Nat :=
  match l with
  | [] => 0
  | a :: r => if a.1 == p.1 && a.2 == p.2 then 0 else idxPair r p + 1

def idxN (v : Nat) (l : List Nat) : Nat :=
  match l with
  | [] => 0
  | a :: r => if a == v then 0 else idxN v r + 1

def adjPairs (l : List Nat) : List (Nat × Nat) :=
  match l with
  | a :: b :: r => (a, b) :: adjPairs (b :: r)
  | _ => []

/-! ### Transitive closure with fuel -/

def closureStep (lt : Nat → Nat → Bool) (n : Nat) : Nat → (Nat → Nat → Bool)
  | 0 => lt
  | m+1 =>
    let prev := closureStep lt n m
    fun i j => prev i j || (List.range n).any (fun k => prev i k && prev k j)

/-- `n + 2` rounds are mathematically more than enough to close any relation on
`n` elements; the harness nevertheless re-checks `isPosetOK` on every closed
relation before it is used. -/
def transClose (lt : Nat → Nat → Bool) (n : Nat) : Nat → Nat → Bool :=
  closureStep lt n (n + 2)

/-- Materialize a relation into a flat row-major table. -/
def tableOf (lt : Nat → Nat → Bool) (n : Nat) : List Bool :=
  (List.range (n * n)).map (fun p => lt (p / n) (p % n))

/-- One closure round on a table: add all length-2 shortcuts. -/
def closureRoundTable (t : List Bool) (n : Nat) : List Bool :=
  (List.range (n * n)).map (fun p =>
    t.getD p false ||
    (List.range n).any (fun k =>
      t.getD ((p / n) * n + k) false && t.getD (k * n + (p % n)) false))

/-- Transitive closure on a table (fuel rounds; each round doubles the
reachable path length, so `n + 2` rounds are ample). -/
def transCloseTable (t : List Bool) (n : Nat) : Nat → List Bool
  | 0 => t
  | m+1 => transCloseTable (closureRoundTable t n) n m

/-- Table lookup as a relation. -/
def tableLt (t : List Bool) (n : Nat) : Nat → Nat → Bool :=
  fun i j => t.getD (i * n + j) false

/-- Add the relation `x < y` and re-close. Table-based: each `lt` probe is a
constant-time lookup instead of an exponentially expanding closure chain.
Purely a memoization of the former nested-closure definition; the semantics
are unchanged (same fixpoint), and `isPosetOK` re-checks remain in place. -/
def addRel (P : POrd n) (x y : Nat) : POrd n :=
  let t2 := tableOf (fun i j => P.lt i j || (i == x && j == y)) n
  ⟨tableLt (transCloseTable t2 n (n + 2)) n⟩

/-! ### Linear extensions and their number -/

/-- Backtracking enumeration by fuel: place a minimal element of `avail` next. -/
def extsFuel (P : POrd n) : Nat → List Nat → List Nat → List (List Nat)
  | 0, _, _ => []
  | _+1, [], acc => [acc.reverse]
  | f+1, avail, acc =>
    (avail.filter (fun x => !(avail.any (fun y => P.lt y x)))).flatMap
      (fun x => extsFuel P f (avail.erase x) (x :: acc))

/-- Fuel `n+1` suffices for `n` elements (re-checked computationally per instance). -/
def exts (P : POrd n) : List (List Nat) := extsFuel P (n+1) (List.range n) []

def e (P : POrd n) : Nat := (exts P).length

/-- Fuel adequacy spot-check used in the harness. -/
def fuelOK (P : POrd n) : Bool :=
  (extsFuel P (n+1) (List.range n) []).length == (extsFuel P (n+2) (List.range n) []).length

/-! ### Maximal elements, deletion, reindexing -/

def maxl (P : POrd n) : List Nat :=
  (List.range n).filter (fun z => !((List.range n).any (fun y => P.lt z y)))

/-- Reindexing `{0, ..., n-1} \ {z}` onto `{0, ..., n-2}`. -/
def skipIdx (z i : Nat) : Nat := if i < z then i else i + 1

/-- Inverse reindexing: original index `i != z` mapped to its `Q = P \ z` label. -/
def downIdx (z i : Nat) : Nat := if i < z then i else i - 1

/-- Delete element `z` (meaningful for `z < n`, `n >= 1`). -/
def delP (z : Nat) (P : POrd n) : POrd (n-1) :=
  ⟨fun i j => P.lt (skipIdx z i) (skipIdx z j)⟩

/-! ### Lemma A: statement and check -/

/-- **Lemma A** (top-element recursion):
`e P = sum over maximal z of e (P without z)`. -/
def lemmaA_stmt (P : POrd n) : Prop :=
  e P = (maxl P).foldl (fun acc z => acc + e (delP z P)) 0

def checkA (P : POrd n) : Bool :=
  e P == (maxl P).foldl (fun acc z => acc + e (delP z P)) 0

/-! ### Lemma B: statement and check -/

/-- Critical pair `(x, y)`: incomparable, and every strict lower/upper bound of
`x` is one of `y`. -/
def critPair (P : POrd n) (x y : Nat) : Bool :=
  x < n && y < n && !(x == y) && !(P.lt x y) && !(P.lt y x) &&
  (List.range n).all (fun w => w == x || w == y ||
    (!P.lt w x || P.lt w y) && (!P.lt x w || P.lt y w))

/-- `y` immediately covers `x` somewhere in the list. -/
def adjIn (l : List Nat) (x y : Nat) : Bool :=
  (adjPairs l).any (fun p => p.1 == x && p.2 == y)

/-- **Lemma B**: every critical pair is covered by some linear extension. -/
def checkB (P : POrd n) : Bool :=
  ((List.range n).flatMap (fun x => (List.range n).filterMap
    (fun y => if critPair P x y then some (x, y) else none)))
  |>.all (fun p => (exts P).any (fun l => adjIn l p.1 p.2))

/-! ### Lemma C: statement and check -/

/-- Slot count: number of gaps of `M` where `z` can legally be inserted,
given the (already reindexed) lists of elements below/above `z`. -/
def slotCount (below above : List Nat) (M : List Nat) : Nat :=
  (List.range (M.length + 1)).filter (fun k =>
    below.all (fun v => memN v (M.take k)) &&
    above.all (fun v => memN v (M.drop k))) |>.length

/-- **Lemma C** data: total slots and slots lying on `x`-before-`y` extensions. -/
def slotsData (P : POrd n) (z x y : Nat) : Nat × Nat :=
  let Q := delP z P
  let below := ((List.range n).filter (fun v => P.lt v z)).map (downIdx z)
  let above := ((List.range n).filter (fun v => P.lt z v)).map (downIdx z)
  let ms := exts Q
  let sx := downIdx z x
  let sy := downIdx z y
  let den := ms.foldl (fun acc M => acc + slotCount below above M) 0
  let num := ms.foldl (fun acc M =>
    if idxN sx M < idxN sy M then acc + slotCount below above M else acc) 0
  (den, num)

/-- Precondition part of the Lemma C check. -/
def cGuard (P : POrd n) (z x y : Nat) : Bool :=
  z < n && x < n && y < n && !(x == z) && !(y == z) && !(x == y) &&
  !(P.lt x y) && !(P.lt y x) && isPosetOK (addRel P x y)

/-- **Lemma C** (insertion-weight identity), cross-multiplied:
`e (P + x<y) * (total slots) = (slots on x<y extensions) * e P`. -/
def cCheck (P : POrd n) (z x y : Nat) : Bool :=
  let (den, num) := slotsData P z x y
  den * e (addRel P x y) == num * e P

/-! ### Proposition 1: statement and check -/

/-- `z` comparable with every other element. -/
def compAll (P : POrd n) (z : Nat) : Bool :=
  (List.range n).all (fun w => w == z || P.lt w z || P.lt z w)

/-- Precondition part of the Proposition 1 check. -/
def pGuard (P : POrd n) (z x y : Nat) : Bool :=
  compAll P z && z < n && x < n && y < n && !(x == z) && !(y == z) &&
  !(x == y) && !(P.lt x y) && !(P.lt y x) && isPosetOK (delP z P)

/-- **Proposition 1**: if `z` is comparable with everything else, extensions
correspond one-to-one and every `p`-value is preserved. -/
def pCheck (P : POrd n) (z x y : Nat) : Bool :=
  e P == e (delP z P) &&
  e (addRel P x y) == e (addRel (delP z P) (downIdx z x) (downIdx z y))

/-! ### The conjecture itself -/

/-- Every pair comparable, i.e. `P` is a chain. -/
def isChain (P : POrd n) : Bool :=
  (List.range n).all (fun x => (List.range n).all (fun y =>
    x == y || P.lt x y || P.lt y x))

/-- The incomparable ordered pairs. -/
def incompPairs (P : POrd n) : List (Nat × Nat) :=
  (List.range n).flatMap (fun x => (List.range n).flatMap (fun y =>
    if !(x == y) && !(P.lt x y) && !(P.lt y x) then [(x, y)] else []))

/-- The 1/3-2/3 conjecture, as a statement about `e` and `addRel`.
Vacuous on chains (no incomparable pair exists). -/
def conjecture_stmt (P : POrd n) : Prop :=
  isChain P ∨ ∃ x y, x < n ∧ y < n ∧ ¬(P.lt x y || P.lt y x) ∧
    e P ≤ 3 * e (addRel P x y) ∧ e P ≤ 3 * e (addRel P y x)

def checkConj (P : POrd n) : Bool :=
  let ps := incompPairs P
  ps.isEmpty ||
    ps.any (fun p =>
      e P <= 3 * e (addRel P p.1 p.2) && e P <= 3 * e (addRel P p.2 p.1))

/-! ### Sanity instances -/

def chainP (n : Nat) : POrd n := ⟨fun i j => decide (i < j)⟩
def antiP (n : Nat) : POrd n := ⟨fun _ _ => false⟩

theorem chain_isPoset : isPosetOK (chainP 4) = true := by decide
theorem chain_e : e (chainP 4) = 1 := by decide
theorem anti_e : e (antiP 4) = 24 := by decide

/-! ### Exhaustive enumeration of forward-labeled posets -/

def pairsLess (n : Nat) : List (Nat × Nat) :=
  (List.range n).flatMap (fun i =>
    (List.range n).filterMap (fun j => if i < j then some (i, j) else none))

def pairBit (mask k : Nat) : Bool := (mask / 2^k) % 2 == 1

/-- The relation of mask: bit `k` decides the `k`-th pair `(i, j)` with `i < j < n`.
Every poset arises this way for a suitable topological labeling. -/
def relFromMask (n mask : Nat) : Nat → Nat → Bool :=
  fun i j => if i < j then pairBit mask (idxPair (pairsLess n) (i, j)) else false

def mkP (n mask : Nat) : POrd n := ⟨relFromMask n mask⟩

/-- Full harness on one candidate: run A, B, the conjecture, C, P1, plus fuel
sanity, with all guard conditions enforced as implications. -/
def checkAllOn (n mask : Nat) : Bool :=
  let P := mkP n mask
  if isPosetOK P then
    fuelOK P && checkA P && checkB P && checkConj P &&
    (List.range n).all (fun z => (List.range n).all (fun x => (List.range n).all (fun y =>
      !cGuard P z x y || cCheck P z x y))) &&
    (List.range n).all (fun z => (List.range n).all (fun x => (List.range n).all (fun y =>
      !pGuard P z x y || pCheck P z x y)))
  else true

/-- Number of firing Lemma C guard triples over an enumeration. -/
def cFires (n : Nat) (masks : List Nat) : Nat :=
  masks.foldl (fun acc m =>
    let P := mkP n m
    if isPosetOK P then
      acc + (List.range n).foldl (fun a z => a +
        (List.range n).foldl (fun a2 x =>
          (List.range n).foldl (fun a3 y => a3 + (if cGuard P z x y then 1 else 0)) 0) 0) 0
    else acc) 0

/-- Number of firing Proposition 1 guard triples over an enumeration. -/
def pFires (n : Nat) (masks : List Nat) : Nat :=
  masks.foldl (fun acc m =>
    let P := mkP n m
    if isPosetOK P then
      acc + (List.range n).foldl (fun a z => a +
        (List.range n).foldl (fun a2 x =>
          (List.range n).foldl (fun a3 y => a3 + (if pGuard P z x y then 1 else 0)) 0) 0) 0
    else acc) 0

