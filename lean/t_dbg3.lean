import Poset13Step3

/-! Locate the first failing (mask, z, Amask) and the failing sub-check. -/

def subChecks (n mask : Nat) : List (String × Bool) :=
  let P := mkP n mask
  if isPosetOK P then
    let nP := eSub P
    let sandwich := ("sandwich", sandwichOK P)
    if isChain P then [sandwich] else
      (List.range n).flatMap (fun z =>
        let m := n - 1
        let Q := delP z P
        let box := boxList P z
        (List.range (2^m)).filterMap (fun Amask =>
          if partOK Q Amask then
            some (s!"z={z} A={Amask}", step3Part P z nP box Q Amask)
          else none)
        |>.filterMap (fun p => if p.2 then none else some p))
      |>.map (fun p => p) |>.take 1 |>.map id
      |>.append [sandwich]
  else []

#eval (List.range 64).filterMap (fun mask =>
  let res := subChecks 4 mask
  if res.any (fun p => !p.2) then some (mask, res) else none)
