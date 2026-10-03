import Poset13Defs

/-!
# Session-3 verification module: theta-structure, heavy layer, refined conservation

Machine-checks the session-3 structural results on every target-class candidate
(`P = Q + z`, `width(Q) <= 2` via an explicit chain partition, `P` not a chain):

* `fStructCheck` (Lemma F structure): with `theta(a_k) = M_A(k)` the cumulative
  box-column mass, `M_A` is monotone; on each chain `D ∩ A` is a prefix and
  `U ∩ A` a suffix of the chain order, so the z-window is the contiguous index
  range `(iD, iU]`.
* `fCutCheck` (Lemma F cut, under fail): window theta-values avoid `[1/3, 2/3]`
  and everything above the last `< 1/3` index exceeds `2/3`.
* `gCheck` (Lemma G, under fail): if the A-window is nonempty some box column
  has mass `> 1/3` (dually rows), and every jump column is squeezed:
  `2/3 - colmass < M_A(c) < 1/3`.
* `s1Check` (refined conservation, under fail): some incomparable cross pair
  `(a_k, b_l)` is balanced and lies in a same-side block
  `(k ≤ i* ∧ l ≤ j*) ∨ (k ≥ i*+1 ∧ l ≥ j*+1)`.
* `s2Check` (mixed corner orientation, under fail): `p(a_{i*}, b_{j*+1}) > 2/3`
  and `p(a_{i*+1}, b_{j*}) < 1/3` whenever those cells exist and are incomparable.
* `sandwichOK` (Lemma H, all posets): for every incomparable pair,
  `e(P∖y) ≤ e(P+x<y)` and `e(P+x<y) + e(P∖x) ≤ e(P)`.

All comparisons are cross-multiplied `Nat` inequalities against `eP = e P`
(mass numerators), so everything is `native_decide`-friendly. Sharding:
`step3All n lo hi` checks masks in `[lo, hi)`.
-/

variable {n : Nat}

/-! ### Bit helpers and subset-DP extension counting (from t_lemmad) -/

def bitOf (i S : Nat) : Bool := S / 2^i % 2 == 1

def clearBit (S i : Nat) : Nat := if bitOf i S then S - 2^i else S

def popcnt (m S : Nat) : Nat := (List.range m).filter (fun i => bitOf i S) |>.length

def countSub (P : POrd n) : Nat → Nat → Nat
  | 0, avail => if avail == 0 then 1 else 0
  | f+1, avail =>
    if avail == 0 then 1 else
      (List.range n).filter (fun x => bitOf x avail && isMinIn P avail x)
        |>.foldl (fun acc x => acc + countSub P f (clearBit avail x)) 0
where isMinIn (P : POrd n) (avail : Nat) (x : Nat) : Bool :=
  (List.range n).all (fun y => y == x || !(P.lt y x) || !bitOf y avail)

/-- Extension count via subset DP (equals `e P`; sanity-checked in harness). -/
def eSub (P : POrd n) : Nat := countSub P n (2^n - 1)

def downOK (Q : POrd m) (S : Nat) : Bool :=
  (List.range m).all (fun i => !bitOf i S ||
    (List.range m).all (fun j => !(Q.lt j i) || bitOf j S))

def boxOK (P : POrd n) (z : Nat) (S : Nat) : Bool :=
  let Dmask := (List.range n).foldl (fun acc v =>
    if P.lt v z then acc + 2^(downIdx z v) else acc) 0
  let Umask := (List.range n).foldl (fun acc v =>
    if P.lt z v then acc + 2^(downIdx z v) else acc) 0
  (S &&& Dmask) == Dmask && (S &&& Umask) == 0

def toOrig (z m : Nat) (S : Nat) : Nat :=
  (List.range m).foldl (fun acc i =>
    if bitOf i S then acc + 2^(skipIdx z i) else acc) 0

/-! ### Chain partitions of Q (bitmasks over Q-labels 0..m-1) -/

def chainMaskOK (Q : POrd m) (S : Nat) : Bool :=
  (List.range m).all (fun i => !bitOf i S ||
    (List.range m).all (fun j => !bitOf j S || i == j ||
      Q.lt i j || Q.lt j i))

def partOK (Q : POrd m) (Amask : Nat) : Bool :=
  let full := 2^m - 1
  Amask != 0 && Amask != full &&
    chainMaskOK Q Amask && chainMaskOK Q (full ^^^ Amask)

/-! ### Per-candidate data -/

/-- The box: ideals of Q containing D(z), avoiding U(z), as Q-label bitmasks. -/
def boxList (P : POrd n) (z : Nat) : List Nat :=
  let m := n - 1
  let Q := delP z P
  (List.range (2^m)).filter (fun S => downOK Q S && boxOK P z S)

/-- Lemma D weight of a box ideal: `e(I) * e(Q∖I)` on original labels. -/
def boxW (P : POrd n) (z m : Nat) (S : Nat) : Nat :=
  let full := 2^m - 1
  countSub P (popcnt m S) (toOrig z m S) *
    countSub P (popcnt m (full ^^^ S)) (toOrig z m (full ^^^ S))

/-- Rank of `x` inside its chain mask = number of chain elements below it. -/
def rankIn (Q : POrd m) (Amask : Nat) (x : Nat) : Nat :=
  (List.range m).filter (fun u => bitOf u Amask && Q.lt u x) |>.length

/-- Cumulative column mass `M_A(k) = μ{S in box : |S∩A| < k}` (numerator over eP). -/
def colCum (P : POrd n) (z : Nat) (Amask : Nat) (box : List Nat) (k : Nat) : Nat :=
  let m := n - 1
  box.foldl (fun acc S =>
    if popcnt m (S &&& Amask) < k then acc + boxW P z m S else acc) 0

/-- Total box mass; equals `e P` by Lemma D (sanity-checked in harness). -/
def boxTotal (P : POrd n) (z m : Nat) (box : List Nat) : Nat :=
  box.foldl (fun acc S => acc + boxW P z m S) 0

/-! ### The fail guard: every z-pair unbalanced (p ∉ [1/3, 2/3]) -/

def failQ (P : POrd n) (z : Nat) (nP : Nat) : Bool :=
  (List.range n).all (fun v =>
    v == z || P.lt v z || P.lt z v ||
    3 * eSub (addRel P z v) < nP || 3 * eSub (addRel P z v) > 2 * nP)

/-! ### Lemma H sandwich (all posets, all incomparable pairs)

`P(y last) ≤ p(x,y) ≤ 1 − P(x last)`; for a MAXIMAL endpoint `w`,
`P(w last) = e(P∖w)/e(P)`, which is the only case the counts encode directly.
-/

def isMaxW (P : POrd n) (w : Nat) : Bool :=
  !((List.range n).any (fun u => !(u == w) && P.lt w u))

def sandwichOK (P : POrd n) : Bool :=
  (List.range n).all (fun x => (List.range n).all (fun y =>
    x == y || P.lt x y || P.lt y x ||
    (let nP := eSub P
     (!isMaxW P y || eSub (delP y P) <= eSub (addRel P x y)) &&
     (!isMaxW P x || eSub (addRel P x y) + eSub (delP x P) <= nP))))

/-! ### The per-(P, z, partition) checks -/

def step3Part (P : POrd n) (z : Nat) (nP : Nat) (box : List Nat)
    (Q : POrd (n-1)) (Amask : Nat) : Bool :=
  let m := n - 1
  let full := 2^m - 1
  let Bmask := full ^^^ Amask
  let s := popcnt m Amask
  let t := popcnt m Bmask
  let Dmask := (List.range n).foldl (fun acc v =>
    if P.lt v z then acc + 2^(downIdx z v) else acc) 0
  let Umask := (List.range n).foldl (fun acc v =>
    if P.lt z v then acc + 2^(downIdx z v) else acc) 0
  let iD := popcnt m (Dmask &&& Amask)
  let jD := popcnt m (Dmask &&& Bmask)
  let iU := s - popcnt m (Umask &&& Amask)
  let jU := t - popcnt m (Umask &&& Bmask)
  let MA := fun k => colCum P z Amask box k
  let MB := fun l => colCum P z Bmask box l
  -- Lemma F structure: chain order respects D-prefix / U-suffix / window
  let fStruct :=
    (List.range m).all (fun x => !bitOf x Amask ||
      (if P.lt (skipIdx z x) z then rankIn Q Amask x < iD
       else if P.lt z (skipIdx z x) then rankIn Q Amask x >= iU
       else rankIn Q Amask x >= iD && rankIn Q Amask x < iU)) &&
    (List.range m).all (fun x => !bitOf x Bmask ||
      (if P.lt (skipIdx z x) z then rankIn Q Bmask x < jD
       else if P.lt z (skipIdx z x) then rankIn Q Bmask x >= jU
       else rankIn Q Bmask x >= jD && rankIn Q Bmask x < jU)) &&
    (List.range (s+1)).all (fun k => MA k <= MA (k+1)) &&
    (List.range (t+1)).all (fun l => MB l <= MB (l+1)) &&
    boxTotal P z m box == nP
  -- Lemma F cut + Lemma G + S1 + S2, all under the fail hypothesis
  let iStar := (List.range (s+1)).foldl (fun acc k =>
    if 3 * MA k < nP then k else acc) 0
  let jStar := (List.range (t+1)).foldl (fun acc l =>
    if 3 * MB l < nP then l else acc) 0
  let fCut :=
    (List.range (s+1)).all (fun k =>
      !(iD < k && k <= iU) || 3 * MA k < nP || 3 * MA k > 2 * nP) &&
    (List.range (t+1)).all (fun l =>
      !(jD < l && l <= jU) || 3 * MB l < nP || 3 * MB l > 2 * nP) &&
    (List.range (s+1)).all (fun k =>
      !(k > iStar) || 3 * MA k > 2 * nP) &&
    (List.range (t+1)).all (fun l =>
      !(l > jStar) || 3 * MB l > 2 * nP)
  -- Lemma G: heavy column/row + squeeze at jump columns/rows
  let colmass := fun c => box.foldl (fun acc S =>
    if popcnt m (S &&& Amask) == c then acc + boxW P z m S else acc) 0
  let rowmass := fun r => box.foldl (fun acc S =>
    if popcnt m (S &&& Bmask) == r then acc + boxW P z m S else acc) 0
  let gCheck :=
    (!(iD < iU) || (List.range (s+1)).any (fun c => 3 * colmass c > nP)) &&
    (!(jD < jU) || (List.range (t+1)).any (fun r => 3 * rowmass r > nP)) &&
    (List.range (s+1)).all (fun c =>
      !(iD < c + 1 && c + 1 <= iU && 3 * colmass c > nP && 3 * MA c < nP) ||
        3 * (MA c + colmass c) > 2 * nP) &&
    (List.range (t+1)).all (fun r =>
      !(jD < r + 1 && r + 1 <= jU && 3 * rowmass r > nP && 3 * MB r < nP) ||
        3 * (MB r + rowmass r) > 2 * nP)
  -- S1: balanced incomparable cross pair inside a same-side block
  let balAt (x y : Nat) : Bool :=
    let ox := skipIdx z x; let oy := skipIdx z y
    !(P.lt ox oy) && !(P.lt oy ox) &&
      let pn := eSub (addRel P ox oy)
      3 * pn >= nP && 3 * pn <= 2 * nP
  let s1 :=
    (List.range m).any (fun x => bitOf x Amask &&
      (List.range m).any (fun y => bitOf y Bmask &&
        let k := rankIn Q Amask x + 1
        let l := rankIn Q Bmask y + 1
        balAt x y &&
          ((k <= iStar && l <= jStar) || (k >= iStar + 1 && l >= jStar + 1))))
  -- S2: mixed corner orientation
  let s2 :=
    ((List.range m).all (fun x => !(bitOf x Amask && rankIn Q Amask x + 1 == iStar) ||
        (List.range m).all (fun y => !(bitOf y Bmask && rankIn Q Bmask y + 1 == jStar + 1) ||
          let ox := skipIdx z x; let oy := skipIdx z y
          P.lt ox oy || P.lt oy ox || 3 * eSub (addRel P ox oy) > 2 * nP))) &&
    ((List.range m).all (fun x => !(bitOf x Amask && rankIn Q Amask x + 1 == iStar + 1) ||
        (List.range m).all (fun y => !(bitOf y Bmask && rankIn Q Bmask y + 1 == jStar) ||
          let ox := skipIdx z x; let oy := skipIdx z y
          P.lt ox oy || P.lt oy ox || 3 * eSub (addRel P ox oy) < nP)))
  fStruct && (!failQ P z nP || (fCut && gCheck && s1 && s2))

/-- Full check on one candidate poset (given by mask): sandwich everywhere,
then for each z admitting a chain partition of Q = P∖z, all partition checks. -/
def step3One (n mask : Nat) : Bool :=
  let P := mkP n mask
  if isPosetOK P then
    sandwichOK P &&
    (isChain P ||
      (List.range n).all (fun z =>
        let m := n - 1
        let Q := delP z P
        let box := boxList P z
        (List.range (2^m)).all (fun Amask =>
          !partOK Q Amask || step3Part P z (eSub P) box Q Amask)))
  else true

/-- Shard over forward-labeled candidates with masks in `[lo, hi)`. -/
def step3All (n lo hi : Nat) : Bool :=
  (List.range (hi - lo)).all (fun d => step3One n (lo + d))
