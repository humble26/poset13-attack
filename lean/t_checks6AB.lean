import Poset13Defs

def checks6AB : Bool := (List.range (2^15)).all (fun m =>
  let P := mkP 6 m
  if isPosetOK P then fuelOK P && checkA P && checkB P else true)
theorem checks6AB_true : checks6AB = true := by native_decide
#print axioms checks6AB_true
