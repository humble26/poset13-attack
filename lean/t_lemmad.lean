import Poset13Defs

/-!
## Lemma D (ideal decomposition) and the strip formula, machine-checked

`e P = Σ_{I ideal of Q=P\z, D(z) ⊆ I, I ∩ U(z) = ∅} e(I) · e(Q \ I)`

`Σ_{I in box, u ∉ I} e(I) · e(Q \ I) = e (P + z<u)`  for u ∥ z  (strip formula)

Subset-DP (`countSub`) computes extension counts of induced posets directly on
bitmasks of the ORIGINAL labels, avoiding all reindexing.
-/

variable {n : Nat}

def bitOf (i S : Nat) : Bool := S / 2^i % 2 == 1

def clearBit (S i : Nat) : Nat := if bitOf i S then S - 2^i else S

def popcnt (m S : Nat) : Nat := (List.range m).filter (fun i => bitOf i S) |>.length

/-- Extensions of the induced poset on `avail` (subset bitmask of original labels):
count topological arrangements placing P-minimal elements of `avail` next. -/
def countSub (P : POrd n) : Nat → Nat → Nat
  | 0, avail => if avail == 0 then 1 else 0
  | f+1, avail =>
    if avail == 0 then 1 else
      (List.range n).filter (fun x => bitOf x avail && isMinIn P avail x)
        |>.foldl (fun acc x => acc + countSub P f (clearBit avail x)) 0
where isMinIn (P : POrd n) (avail : Nat) (x : Nat) : Bool :=
  (List.range n).all (fun y => y == x || !(P.lt y x) || !bitOf y avail)

/-- Down-closed subsets (ideals) of `Q`, as bitmasks over `0..m-1`. -/
def downOK (Q : POrd m) (S : Nat) : Bool :=
  (List.range m).all (fun i => !bitOf i S ||
    (List.range m).all (fun j => !(Q.lt j i) || bitOf j S))

/-- The "box" of ideals containing D(z) and avoiding U(z). -/
def boxOK (P : POrd n) (z : Nat) (S : Nat) : Bool :=
  let Dmask := (List.range n).foldl (fun acc v =>
    if P.lt v z then acc + 2^(downIdx z v) else acc) 0
  let Umask := (List.range n).foldl (fun acc v =>
    if P.lt z v then acc + 2^(downIdx z v) else acc) 0
  (S &&& Dmask) == Dmask && (S &&& Umask) == 0

/-- Convert a Q-label bitmask to an original-label bitmask. -/
def toOrig (z m : Nat) (S : Nat) : Nat :=
  (List.range m).foldl (fun acc i =>
    if bitOf i S then acc + 2^(skipIdx z i) else acc) 0

/-- Lemma D check (cross-multiplied form is an exact equality of Nats). -/
def lemmaDCheck (P : POrd n) (z : Nat) : Bool :=
  let m := n - 1
  let Q := delP z P
  let full := 2^m - 1
  let total := (List.range (2^m)).filter (fun S => downOK Q S && boxOK P z S)
    |>.foldl (fun acc S =>
      acc + countSub P (popcnt m S) (toOrig z m S) *
            countSub P (popcnt m (full ^^^ S)) (toOrig z m (full ^^^ S))) 0
  total == e P

/-- Strip formula check: mass of box ideals avoiding u equals e(P + z<u). -/
def stripCheck (P : POrd n) (z u : Nat) : Bool :=
  !(u == z) && !(P.lt z u) && !(P.lt u z) &&
  (let m := n - 1
   let Q := delP z P
   let full := 2^m - 1
   let du := downIdx z u
   let total := (List.range (2^m)).filter (fun S =>
      downOK Q S && boxOK P z S && !bitOf du S)
     |>.foldl (fun acc S =>
      acc + countSub P (popcnt m S) (toOrig z m S) *
            countSub P (popcnt m (full ^^^ S)) (toOrig z m (full ^^^ S))) 0
   total == e (addRel P z u))

def lemmaDAll (n : Nat) : Bool := (List.range (2^(n*(n-1)/2))).all (fun mask =>
  let P := mkP n mask
  if isPosetOK P then
    (List.range n).all (fun z =>
      lemmaDCheck P z &&
      (List.range n).all (fun u => (u == z) || (P.lt z u) || (P.lt u z) || stripCheck P z u))
  else true)


def ld4 : Bool := lemmaDAll 4
theorem ld4_true : ld4 = true := by native_decide
#print axioms ld4_true

def ld5 : Bool := lemmaDAll 5
theorem ld5_true : ld5 = true := by native_decide
#print axioms ld5_true
