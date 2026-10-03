import Poset13Step3

/-! n=4 smoke: subset-DP `eSub` anchored against list enumeration `e`,
then the full session-3 harness over all 64 candidates. -/

def eSubAll (n : Nat) : Bool := (List.range (2^(n*(n-1)/2))).all (fun mask =>
  let P := mkP n mask
  if isPosetOK P then eSub P == e P else true)

def s4a : Bool := eSubAll 4
theorem s4a_true : s4a = true := by native_decide
#print axioms s4a_true

def s4 : Bool := step3All 4 0 64
theorem s4_true : s4 = true := by native_decide
#print axioms s4_true
